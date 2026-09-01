"""design() -- one protein into one backbone, verified, end to end.

This is PLAN.md's walking skeleton: a single design, proven and annotated, with
every lane wired together for the first time. It is deliberately NARROW. It ships
one candidate, no gallery. It scores nothing -- every objective is reported
`unavailable` with a reason, never omitted, so the scorecard cannot look complete
when it is not. It carries no baseline: `native_baseline` stays None and the
"the native sequence is often the right answer" sentence never renders, because
that claim is about a real wild-type CDS and the skeleton has none. What it DOES
do is real: it refuses to emit a construct the validator cannot prove, it never
touches the user's backbone (I9 armed), and it exports the GenBank the user takes
away.

The one subtlety worth reading twice is the flank orientation. Tier A seeds its
automaton with the immutable backbone on either side of the insert so a forbidden
site formed half by the vector and half by the first codon is excluded by
construction. Those flanks are in CODING orientation -- for a reverse-strand
insertion site the coding-5' neighbour is the reverse complement of the
backbone's DOWNSTREAM bases -- and this is the one place a sign error comes back
silently clean, so it has its own test.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from bt5.codon.tables import FileTableProvider
from bt5.core.context import (
    BiosecurityVerdict,
    ContextSlot,
    DesignContext,
    HostId,
    Modality,
)
from bt5.core.result import Candidate, DesignResult, ObjectiveScore, ScoreCard
from bt5.core.services import FoldEngine
from bt5.core.types import Construct, Interval, reverse_complement
from bt5.design.catalog import partition_forbidden, scored_objectives
from bt5.design.errors import DesignError
from bt5.design.provenance import (
    build_provenance,
    constraint_set_hash,
    design_hash_context,
)
from bt5.design.sites import choose_site
from bt5.rules.vendors import DEFAULT_SELECTION, VendorSelection
from bt5.score import build_report, design_hash, render
from bt5.score.report import QcReport
from bt5.solver.catalog import RuleSet, build_rule_set, default_services
from bt5.solver.pipeline import OptimizeResult, optimize
from bt5.structure.vienna import degradation_reason
from bt5.vector import Assembly, annotate, assemble, construct_to_record, write_genbank
from bt5.vector.backbone import InsertionSite, VectorBackbone

_GC_WINDOW = 50


@dataclass(frozen=True)
class DesignInputs:
    """The resolved inputs of one design, recorded so a result can be read
    without the call site."""

    protein: str
    table_id: int
    modality: Modality
    hosts: tuple[HostId, ...]
    seed: int
    vendors: VendorSelection
    preset_id: str | None
    site_label: str


@dataclass(frozen=True)
class SkeletonResult:
    """One design, end to end. Everything a caller or a test needs to see."""

    inputs: DesignInputs
    result: DesignResult
    report: QcReport
    rendered: str
    assembly: Assembly
    optimize_result: OptimizeResult
    genbank: str
    notes: tuple[str, ...]


def _bases_before(seq: str, pos: int, k: int) -> str:
    """`k` backbone bases ending just before `pos`, wrapping a circular origin."""
    n = len(seq)
    if k <= 0 or n == 0:
        return ""
    return "".join(seq[(pos - k + i) % n] for i in range(k))


def _bases_after(seq: str, pos: int, k: int) -> str:
    """`k` backbone bases starting at `pos`, wrapping a circular origin."""
    n = len(seq)
    if k <= 0 or n == 0:
        return ""
    return "".join(seq[(pos + i) % n] for i in range(k))


def _coding_flanks(
    backbone: VectorBackbone, site: InsertionSite, usable_forbidden: Sequence[str]
) -> tuple[str, str]:
    """The immutable neighbours of the insert, in CODING orientation.

    Exactly `longest_pattern - 1` bases, which is all a junction-spanning motif
    can reach -- never the whole backbone, because Tier A consumes the flank on
    every iteration and `optimal_back_translate` refuses a flank that itself
    contains a forbidden motif. For a reverse-strand site the coding order runs
    against the plus strand, so each flank is the reverse complement of the
    backbone on the OPPOSITE side from the forward case.
    """
    longest = max((len(m) for m in usable_forbidden), default=1)
    k = longest - 1
    if k <= 0:
        return "", ""
    seq = backbone.sequence
    start = site.interval.start
    end = site.interval.end
    if site.strand == 1:
        return _bases_before(seq, start, k), _bases_after(seq, end, k)
    # Reverse strand: coding 5' neighbour is the revcomp of the plus-strand bases
    # 3' of the insert; coding 3' neighbour the revcomp of those 5' of it.
    left = reverse_complement(_bases_after(seq, end, k))
    right = reverse_complement(_bases_before(seq, start, k))
    return left, right


def _context(
    *,
    modality: Modality,
    hosts: Sequence[HostId],
    table_id: int,
    cassette_orientation: int,
    seed: int = 0,
) -> DesignContext:
    """A single producer slot for the walking skeleton.

    The compound three-slot case (propagate in E. coli, produce in HEK293,
    transduce a target) is what the data model exists for, but one producer slot
    is enough to drive one design end to end and to make a modality-dependent rule
    like d4 resolve to HARD_REPAIR. Multi-slot context is the ranking increment's.
    """
    slot = ContextSlot(role="producer", host=hosts[0], modality=modality, table_id=table_id)
    return DesignContext(
        slots=(slot,),
        cassette_orientation=cassette_orientation,  # type: ignore[arg-type]
        seed=seed,
        screen=BiosecurityVerdict("not_run", None, "protein-level screening did not run"),
    )


def design(
    *,
    backbone: VectorBackbone,
    protein: str,
    table_id: int,
    modality: Modality,
    hosts: Sequence[HostId],
    site: InsertionSite | None = None,
    site_interval: Interval | None = None,
    site_label: str | None = None,
    seed: int = 0,
    vendors: VendorSelection = DEFAULT_SELECTION,
    preset_id: str | None = None,
    fold: FoldEngine | None = None,
    max_candidates: int = 256,
) -> SkeletonResult:
    """Design one protein into one backbone and return a proven, annotated result.

    `table_id` is required and never defaulted (CLAUDE.md 3.1). `protein` must
    start with the initiator, because the assembler builds the cassette
    `starts_at_initiator=True`.
    """
    if not hosts:
        raise DesignError("at least one host is required to build a design context")
    if not protein.startswith("M"):
        raise DesignError(
            f"protein must start with the initiator M; got {protein[:1]!r}. The "
            f"assembler builds the cassette with starts_at_initiator=True."
        )

    services = default_services(seed=seed, fold=fold)
    # From FileTableProvider directly so the type is NcbiGeneticCode (what the
    # solver's DP requires), not the wider GeneticCode the Services protocol names.
    code = FileTableProvider().genetic_code(table_id)

    chosen = choose_site(
        backbone,
        table_id=table_id,
        site=site,
        site_interval=site_interval,
        site_label=site_label,
    )

    ctx = _context(
        modality=modality,
        hosts=hosts,
        table_id=table_id,
        cassette_orientation=chosen.strand,
        seed=seed,
    )

    # The solver lane runs the catalog: one rule set, from which the forbidden
    # motifs, the breach finder, the per-rule policies and the oracle's GC band
    # all derive so they cannot disagree (#68). The design lane adds only what the
    # solver has no backbone to know: which forbidden motifs the vector already
    # carries and cannot recode away.
    rules = build_rule_set(ctx, services, vendors=vendors)
    usable_forbidden, carried_forbidden = partition_forbidden(rules.forbidden(), backbone, chosen)

    gc_bounds = rules.oracle_bounds().gc_bounds
    left_flank, right_flank = _coding_flanks(backbone, chosen, usable_forbidden)

    def assembler(cds: str) -> Construct:
        return assemble(backbone, cds, protein=protein, table_id=table_id, site=chosen).construct

    # The I9 reference is length-only (its CDS span is filler), so it can be built
    # from any valid CDS of the emitted length before the real one exists.
    placeholder = "ATG" + "AAA" * (len(protein) - 1) + "TAA"
    reference = assemble(
        backbone, placeholder, protein=protein, table_id=table_id, site=chosen
    ).reference

    # optimize() directly, not solver.optimize_with, because the forbidden set it
    # is given is the PARTITIONED one -- the carried motifs excluded from both the
    # automaton and the validator.
    result = optimize(
        protein,
        code,
        assemble=assembler,
        find_breaches=rules.breach_finder(),
        forbidden=usable_forbidden,
        score=None,
        gc_bounds=gc_bounds,
        gc_window=_GC_WINDOW,
        left_flank=left_flank,
        right_flank=right_flank,
        seed=seed,
        table_id=table_id,
        policies=rules.policies(_GC_WINDOW),
        max_candidates=max_candidates,
        original_backbone=reference,
    )

    final_cds = result.cds
    final_assembly = assemble(backbone, final_cds, protein=protein, table_id=table_id, site=chosen)

    # --- hashes -------------------------------------------------------------
    dhash = design_hash(
        final_cds,
        context=design_hash_context(
            backbone_sequence=backbone.sequence,
            table_id=table_id,
            forbidden=usable_forbidden,
            gc_bounds=gc_bounds,
            vendors=vendors,
        ),
    )
    constraint_hash = constraint_set_hash(
        rules.specs, ctx, forbidden=usable_forbidden, gc_bounds=gc_bounds, vendors=vendors
    )

    # --- annotate + GenBank round-trip -------------------------------------
    annotated, _report = annotate(final_assembly, source_name=backbone.name, design_hash=dhash)
    genbank = write_genbank(construct_to_record(annotated, name=backbone.name))

    # --- scorecard: every objective unavailable, honestly ------------------
    scored = tuple(
        ObjectiveScore.unavailable(
            spec.id, spec.unit, "ranking not computed in the walking skeleton"
        )
        for spec in scored_objectives(rules)
    )
    # HARD_CHECK findings (an over-length fragment, an ITR palindrome): real,
    # reported, never chased by the solver.
    hard_check_breaches = rules.advise()(annotated)
    scorecard = ScoreCard(scores=scored, hard_checks=hard_check_breaches, total=0.0)

    # --- advisories + degradations -----------------------------------------
    advisories = tuple(b.message for b in result.repair_outcome.advisory) + tuple(
        f"backbone carries the forbidden motif {m}; it cannot be recoded away and "
        f"is excluded from enforcement"
        for m in carried_forbidden
    )
    degradations = _degradations(rules, carried_forbidden)

    provenance = build_provenance(
        seed=seed,
        table_id=table_id,
        fold=services.fold,
        constraint_hash=constraint_hash,
        degradations=degradations,
    )

    candidate = Candidate(
        label="design",
        construct=annotated,
        cds=final_cds,
        scorecard=scorecard,
        design_hash=dhash,
    )
    design_result = DesignResult(
        candidates=(candidate,),
        native_baseline=None,  # deferred; the wild-type sentence never renders
        conflicts=(),
        provenance=provenance,
    )

    report = build_report(
        design_result,
        candidate,
        translation_table_id=table_id,
        preset_id=preset_id or "",
        vendor=vendors.keys[0],
        advisories=advisories,
    )

    inputs = DesignInputs(
        protein=protein,
        table_id=table_id,
        modality=modality,
        hosts=tuple(hosts),
        seed=seed,
        vendors=vendors,
        preset_id=preset_id,
        site_label=chosen.label,
    )
    notes = advisories + tuple(
        n.render(final_assembly.construct.length) for n in final_assembly.notes
    )

    return SkeletonResult(
        inputs=inputs,
        result=design_result,
        report=report,
        rendered=render(report),
        assembly=final_assembly,
        optimize_result=result,
        genbank=genbank,
        notes=notes,
    )


def _degradations(rules: RuleSet, carried: Sequence[str]) -> tuple[str, ...]:
    """Every reason this skeleton run is not a complete evaluation, stated.

    A set the report renders and a test pins by equality: a new silent
    degradation must fail that test rather than slip in unremarked.
    """
    out: list[str] = []
    fold = degradation_reason()
    if fold is not None:
        out.append(fold)
    out.append("ranking not computed: no null distribution and no percentiles")
    out.append("protein-level biosecurity screening did not run")
    out.append("single candidate only: no gallery")
    for spec_id in rules.unrunnable:
        out.append(
            f"rule {spec_id} not run: its thresholds are calibrated against a folding "
            f"engine that is not available"
        )
    for motif in carried:
        out.append(f"forbidden motif {motif} carried by the backbone, excluded from enforcement")
    return tuple(out)
