"""design() -- one protein into one backbone, ranked, verified, end to end.

PLAN.md's v1 bar: protein -> validated -> screened -> CDS planned -> spliced into
the circular backbone -> Tier-A DP -> Tier-B repair -> independent
verify_construct -> normalized scorecard -> 5-candidate gallery -> annotated
GenBank + order CSV, under 10 s. This is that path. The walking skeleton it
replaces shipped one candidate and scored nothing, and said so in its own
docstring; what it did NOT do was the ranking, and ranks are the product.

What the ranking is allowed to claim is narrower than a reader expects, and the
narrowness is deliberate. Every number that leaves here is a PERCENTILE against
a distribution of random synonymous variants of the same protein in the same
construct -- "this design's CAI sits above 94% of random synonymous versions of
this gene in this vector". It is not a predicted expression level, titer, yield
or fold-improvement; all computable design features together explain 5-31% of
protein-level variance, and nine benchmarked commercial optimizers were a coin
flip against native sequence. `native_baseline` -- "do not optimize" -- is a
first-class candidate for exactly that reason: populated when the caller
supplies a real wild-type CDS, and left None with the reason stated when they do
not. It is never manufactured by back-translating.

Two seams are worth reading twice.

The flank orientation, inherited from the skeleton. Tier A seeds its automaton
with the immutable backbone on either side of the insert so a forbidden site
formed half by the vector and half by the first codon is excluded by
construction. Those flanks are in CODING orientation -- for a reverse-strand
insertion site the coding-5' neighbour is the reverse complement of the
backbone's DOWNSTREAM bases -- and this is the one place a sign error comes back
silently clean, so it has its own test.

And the biosecurity verdict, which this lane RENDERS and does not produce.
`design()` takes whatever `BiosecurityVerdict` it is handed and defaults to
`not_run`; the report never prints "clear" for a screen that did not run, and a
non-clear verdict is a degradation, which keeps `QcReport.is_complete` False
until a real screen has actually run. Making the screen real is M8's lane.
"""

from __future__ import annotations

import functools
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from bt5.codon.tables import CodonUsage, FileTableProvider, NcbiGeneticCode
from bt5.core.context import (
    BiosecurityVerdict,
    ContextSlot,
    DesignContext,
    HostId,
    Modality,
)
from bt5.core.registry import all_specs
from bt5.core.result import Candidate, DesignResult, ScoreCard, VerificationError
from bt5.core.services import FoldEngine
from bt5.core.spec import Enforcement
from bt5.core.types import Construct, Interval, reverse_complement
from bt5.design.catalog import partition_forbidden, scored_objectives
from bt5.design.errors import DesignError
from bt5.design.gallery import DEFAULT_SWEEP_STEPS, SolveSpace, sweep_designs
from bt5.design.order import entries_for, order_csv
from bt5.design.provenance import (
    build_provenance,
    constraint_set_hash,
    design_hash_context,
)
from bt5.design.ranking import Nulls, build_nulls, score_candidate
from bt5.design.sites import choose_site
from bt5.rules.vendors import DEFAULT_SELECTION, VendorSelection
from bt5.score import build_report, design_hash, render
from bt5.score.distance import distance_matrix
from bt5.score.gallery import Gallery
from bt5.score.order import OrderEntry
from bt5.score.report import QcReport
from bt5.solver.catalog import RuleSet, build_rule_set, default_services
from bt5.solver.pipeline import OptimizeResult
from bt5.structure.vienna import degradation_reason
from bt5.vector import Assembly, annotate, assemble, construct_to_record, write_genbank
from bt5.vector.backbone import InsertionSite, VectorBackbone
from bt5.verify import verify_construct

_GC_WINDOW = 50

#: The panel size PLAN asks for. `build_gallery` accepts 3-8.
DEFAULT_GALLERY_SIZE = 5

#: The verdict a caller who has not screened gets. Never "clear": a report that
#: implies a clean screen nobody ran is the failure the whole biosecurity
#: posture exists to prevent.
UNSCREENED = BiosecurityVerdict("not_run", None, "protein-level screening did not run")


@dataclass(frozen=True)
class DesignInputs:
    """The resolved inputs of one design, recorded so a result can be read
    without the call site.

    Fields added since the walking skeleton are defaulted, so every existing
    construction of this type still works.
    """

    protein: str
    table_id: int
    modality: Modality
    hosts: tuple[HostId, ...]
    seed: int
    vendors: VendorSelection
    preset_id: str | None
    site_label: str
    gallery_size: int = DEFAULT_GALLERY_SIZE
    sweep_steps: int = DEFAULT_SWEEP_STEPS
    native_supplied: bool = False


@dataclass(frozen=True)
class SkeletonResult:
    """One design run, end to end. Everything a caller or a test needs to see.

    Fields are only ever ADDED here. S5's CLI is built against this type, so a
    removal or a rename breaks a lane that cannot see this file.

    `result.candidates` is ranked best-first by `ScoreCard.total`, and
    `assembly`, `optimize_result` and `genbank` all describe `candidates[0]` --
    the top-ranked design, which is the one a user takes away.
    """

    inputs: DesignInputs
    result: DesignResult
    report: QcReport
    rendered: str
    assembly: Assembly
    optimize_result: OptimizeResult
    genbank: str
    notes: tuple[str, ...]
    #: The sweep's own evidence: how many weight vectors solved, how many
    #: distinct designs they produced, and the minimum pairwise codon distance
    #: across the panel -- the number gate G4 reads. None when the sweep produced
    #: nothing and the design fell back to a single unsteered solve.
    gallery: Gallery | None = None
    #: One order line per orderable candidate, and the CSV of them. PLAN locks
    #: v1's output as "annotated construct + vendor order CSV"; this is the
    #: second half of it.
    orders: tuple[OrderEntry, ...] = ()
    order_csv: str = ""
    #: The nulls every percentile on the report was measured against, kept so a
    #: caller can show the distribution behind a number and not only the number.
    nulls: Nulls | None = None

    @property
    def meets_g4(self) -> bool:
        """Did the panel clear gate G4's 15% minimum pairwise codon distance?

        A False here invalidates a PRODUCT decision, not a technical one: if the
        sweep cannot produce genuinely distinct designs then the gallery is not
        a feature and a UI built on it is a lie. It is reported, never lowered.
        """
        return self.gallery is not None and self.gallery.meets_g4


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
    screen: BiosecurityVerdict = UNSCREENED,
) -> DesignContext:
    """A single producer slot, and whatever screen verdict the caller holds.

    The compound three-slot case (propagate in E. coli, produce in HEK293,
    transduce a target) is what the data model exists for, but one producer slot
    is enough to drive one design end to end and to make a modality-dependent
    rule like d4 resolve to HARD_REPAIR. Multi-slot context is a later
    increment's.

    `screen` is passed through untouched. This lane consumes a verdict; it does
    not produce one, and it may not upgrade a `not_run` into anything.
    """
    slot = ContextSlot(role="producer", host=hosts[0], modality=modality, table_id=table_id)
    return DesignContext(
        slots=(slot,),
        cassette_orientation=cassette_orientation,  # type: ignore[arg-type]
        seed=seed,
        screen=screen,
    )


def _host_usage(host: HostId) -> tuple[CodonUsage | None, str | None]:
    """The host's codon-usage table, or None plus the degradation to report.

    Absence is a gap in BT5's data, not an error in the request, so it degrades
    the null from host-frequency sampling to uniform-synonymous sampling and
    SAYS so. Inventing a usage table, or silently borrowing a related host's,
    would make every percentile in the report a measurement against a
    distribution the user was never told about.
    """
    try:
        return FileTableProvider().usage(str(host)), None
    except (FileNotFoundError, KeyError):
        return None, (
            f"no codon usage table on file for host {host}; the null was sampled "
            f"uniformly over synonymous codons rather than at host frequency, and "
            f"the percentiles are against that null"
        )


def _unrunnable_objectives(rules: RuleSet) -> dict[str, tuple[str, str]]:
    """SOFT rules `check_engine_calibration` dropped, so they can still be named.

    `build_rule_set` instantiates only the runnable rules, so an objective whose
    thresholds need a folding engine that is not installed is absent from
    `rules.specs` entirely -- and an objective absent from the scorecard looks
    exactly like one that was never configured. Recovering the class from the
    registry costs one lookup and turns a silent omission into a stated one.
    """
    by_id = {cls.id: cls for cls in all_specs()}
    out: dict[str, tuple[str, str]] = {}
    for spec_id in rules.unrunnable:
        cls = by_id.get(spec_id)
        if cls is None or cls.enforcement is not Enforcement.SOFT:
            continue
        out[spec_id] = (
            cls.unit,
            "not evaluated: this objective's thresholds are calibrated against a "
            "folding engine that is not available",
        )
    return out


def _check_native_encodes(
    native_cds: str, *, protein: str, code: NcbiGeneticCode, table_id: int
) -> None:
    """Refuse a native CDS that encodes something else, BEFORE anything is solved.

    A caller error that costs a 20-vector sweep and a 200-variant null before it
    is reported is a caller error reported badly, so this runs at the top of
    `design()` rather than where the baseline is assembled.
    """
    if len(native_cds) % 3:
        raise DesignError(
            f"native_cds is {len(native_cds)} nt, not a whole number of codons. A "
            f"baseline that is not in frame is not a baseline."
        )
    translated = code.translate(native_cds)
    if not translated.endswith("*") or translated[:-1] != protein:
        raise DesignError(
            f"native_cds does not encode the given protein under NCBI table "
            f"{table_id}: it translates to {translated[:12]!r}... with the protein "
            f"starting {protein[:12]!r}. A baseline encoding a different protein is "
            f"not a baseline."
        )


def _native_assembly(
    native_cds: str,
    *,
    protein: str,
    backbone: VectorBackbone,
    table_id: int,
    site: InsertionSite,
    forbidden: Sequence[str],
    gc_bounds: tuple[float, float] | None,
    reference: Construct,
) -> tuple[Assembly | None, str | None]:
    """Assemble and PROVE the caller's wild-type CDS, or say why it was not used.

    "Do not optimize" is a real answer and frequently the right one for
    homologous mammalian expression, so the native sequence is a candidate a
    user can actually order -- which means it goes through the same independent
    validator as every design. A native CDS the validator refuses is reported as
    refused, with the invariant named; it is never quietly dropped, and it is
    never replaced by a back-translation. A manufactured "native" sequence is not
    a baseline, it is a design wearing the word.
    """
    native = assemble(backbone, native_cds, protein=protein, table_id=table_id, site=site)
    try:
        verify_construct(
            native.construct,
            protein=protein,
            table_id=table_id,
            forbidden=forbidden,
            gc_bounds=gc_bounds,
            original_backbone=reference,
        )
    except VerificationError as exc:
        return None, (
            f"the native CDS was supplied but is not offered as a candidate: the "
            f"independent validator refused it ({exc.invariant}: {exc.detail}). BT5 "
            f"does not offer a sequence it cannot prove."
        )
    return native, None


def _distances(sequences: Sequence[str], labels: Sequence[str]) -> dict[str, dict[str, float]]:
    """label -> {other label -> codon distance}, for `Candidate.codon_distance_to`.

    Codon distance rather than nucleotide distance, and that choice is the whole
    point: two sequences can differ at 30% of their BASES while encoding every
    difference as third-position wobble, and a panel selected on base distance
    would present five of those as five options.
    """
    matrix = distance_matrix(list(sequences))
    return {
        label: {other: matrix[i][j] for j, other in enumerate(labels) if j != i}
        for i, label in enumerate(labels)
    }


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
    native_cds: str | None = None,
    screen: BiosecurityVerdict = UNSCREENED,
    gallery_size: int = DEFAULT_GALLERY_SIZE,
    sweep_steps: int = DEFAULT_SWEEP_STEPS,
    null_sizes: Mapping[str, int] | None = None,
) -> SkeletonResult:
    """Design one protein into one backbone and return a ranked, proven result.

    `table_id` is required and never defaulted (CLAUDE.md 3.1). `protein` must
    start with the initiator, because the assembler builds the cassette
    `starts_at_initiator=True`.

    Every parameter added since the walking skeleton is keyword-only with a
    default, so the frozen signature this lane publishes to the CLI still accepts
    every call that worked before.

    `native_cds` is the caller's real wild-type CDS. There is deliberately no way
    to ask BT5 to invent one.
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
    if native_cds is not None:
        _check_native_encodes(native_cds, protein=protein, code=code, table_id=table_id)

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
        screen=screen,
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

    usage, usage_degradation = _host_usage(hosts[0])

    # optimize() through SolveSpace, not solver.optimize_with, because the
    # forbidden set it is given is the PARTITIONED one -- the carried motifs
    # excluded from both the automaton and the validator.
    space = SolveSpace(
        protein=protein,
        code=code,
        assemble=assembler,
        forbidden=usable_forbidden,
        seed=seed,
        table_id=table_id,
        usage=dict(usage.w) if usage is not None else {},
        find_breaches=rules.breach_finder(),
        gc_bounds=gc_bounds,
        gc_window=_GC_WINDOW,
        left_flank=left_flank,
        right_flank=right_flank,
        policies=rules.policies(_GC_WINDOW),
        max_candidates=max_candidates,
        reference=reference,
    )

    gallery, picks, solved, gallery_degradation = _panel(
        space, steps=sweep_steps, k=gallery_size
    )

    assemblies = {
        cds: assemble(backbone, cds, protein=protein, table_id=table_id, site=chosen)
        for cds in picks
    }
    # ONE catalog pass per candidate: `findings()` yields both the HARD_CHECK
    # findings and every rule's raw score, and evaluating the catalog twice per
    # candidate is the difference between a design inside PLAN's 10 s bar and one
    # outside it.
    findings = {cds: rules.findings(assemblies[cds].construct) for cds in picks}
    evaluations = {
        cds: {ev.spec_id: ev for ev in found.evaluations} for cds, found in findings.items()
    }

    # --- the null, built once and shared by every candidate ------------------
    # Memoised because each objective's null re-draws the SAME seeded variant
    # stream: a moderate objective's 40 variants are a prefix of a cheap one's
    # 200, so without this the assembler runs once per objective per variant.
    build_variant = functools.lru_cache(maxsize=None)(assembler)
    objectives = scored_objectives(rules)
    nulls = build_nulls(
        picks[0],
        objectives,
        anchor=evaluations[picks[0]],
        code=code,
        ctx=ctx,
        svc=services,
        build=build_variant,
        seed=seed,
        usage=usage,
        sizes=null_sizes,
    )
    extra_unavailable = _unrunnable_objectives(rules)

    # --- score every candidate, then rank ------------------------------------
    hash_context = design_hash_context(
        backbone_sequence=backbone.sequence,
        table_id=table_id,
        forbidden=usable_forbidden,
        gc_bounds=gc_bounds,
        vendors=vendors,
    )

    cards: dict[str, ScoreCard] = {
        cds: score_candidate(
            evaluations[cds],
            objectives,
            nulls=nulls,
            hard_checks=findings[cds].hard_check,
            extra_unavailable=extra_unavailable,
        )
        for cds in picks
    }
    # Best total first; the sweep's own order is the tie-break, so two candidates
    # that score identically rank deterministically rather than by sort accident.
    ordered_cds = sorted(picks, key=lambda cds: (-cards[cds].total, picks.index(cds)))
    labels = {cds: f"design_{rank + 1}" for rank, cds in enumerate(ordered_cds)}

    native, native_degradation = (
        _native_assembly(
            native_cds,
            protein=protein,
            backbone=backbone,
            table_id=table_id,
            site=chosen,
            forbidden=usable_forbidden,
            gc_bounds=gc_bounds,
            reference=reference,
        )
        if native_cds is not None
        else (None, None)
    )

    distance_sequences = [*ordered_cds]
    distance_labels = [labels[cds] for cds in ordered_cds]
    if native is not None and native_cds is not None:
        distance_sequences.append(native_cds)
        distance_labels.append("native_baseline")
    distances = _distances(distance_sequences, distance_labels)

    candidates = tuple(
        Candidate(
            label=labels[cds],
            construct=assemblies[cds].construct,
            cds=cds,
            scorecard=cards[cds],
            design_hash=design_hash(cds, context=hash_context),
            codon_distance_to=distances[labels[cds]],
        )
        for cds in ordered_cds
    )

    native_baseline: Candidate | None = None
    if native is not None and native_cds is not None:
        native_findings = rules.findings(native.construct)
        native_baseline = Candidate(
            label="native_baseline",
            construct=native.construct,
            cds=native_cds,
            scorecard=score_candidate(
                {ev.spec_id: ev for ev in native_findings.evaluations},
                objectives,
                nulls=nulls,
                hard_checks=native_findings.hard_check,
                extra_unavailable=extra_unavailable,
            ),
            design_hash=design_hash(native_cds, context=hash_context),
            codon_distance_to=distances["native_baseline"],
        )

    # --- annotate the winner, and export ------------------------------------
    winner = candidates[0]
    final_assembly = assemblies[winner.cds]
    optimize_result = solved[winner.cds]
    annotated, _annotation_report = annotate(
        final_assembly, source_name=backbone.name, design_hash=winner.design_hash
    )
    genbank = write_genbank(construct_to_record(annotated, name=backbone.name))

    advisories = tuple(b.message for b in optimize_result.repair_outcome.advisory) + tuple(
        f"backbone carries the forbidden motif {m}; it cannot be recoded away and "
        f"is excluded from enforcement"
        for m in carried_forbidden
    )
    degradations = _degradations(
        rules,
        carried_forbidden,
        screen=screen,
        gallery=gallery_degradation,
        usage=usage_degradation,
        native=native_degradation,
        native_supplied=native_cds is not None,
        nulls=nulls,
    )
    provenance = build_provenance(
        seed=seed,
        table_id=table_id,
        fold=services.fold,
        constraint_hash=constraint_set_hash(
            rules.specs, ctx, forbidden=usable_forbidden, gc_bounds=gc_bounds, vendors=vendors
        ),
        degradations=degradations,
    )
    design_result = DesignResult(
        candidates=candidates,
        native_baseline=native_baseline,
        conflicts=(),
        provenance=provenance,
    )
    report = build_report(
        design_result,
        winner,
        translation_table_id=table_id,
        preset_id=preset_id or "",
        vendor=vendors.keys[0],
        advisories=advisories,
    )

    # The native baseline is on the order file too. "Do not optimize" is only a
    # real option if the user can actually order the tube.
    orderable = [*candidates, *([native_baseline] if native_baseline is not None else [])]
    inputs = DesignInputs(
        protein=protein,
        table_id=table_id,
        modality=modality,
        hosts=tuple(hosts),
        seed=seed,
        vendors=vendors,
        preset_id=preset_id,
        site_label=chosen.label,
        gallery_size=gallery_size,
        sweep_steps=sweep_steps,
        native_supplied=native_cds is not None,
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
        optimize_result=optimize_result,
        genbank=genbank,
        notes=notes,
        gallery=gallery,
        orders=entries_for(orderable, construct_name=backbone.name),
        order_csv=order_csv(orderable, construct_name=backbone.name),
        nulls=nulls,
    )


def _panel(
    space: SolveSpace, *, steps: int, k: int
) -> tuple[Gallery | None, list[str], dict[str, OptimizeResult], str | None]:
    """The candidate panel, or the honest single candidate that replaces it.

    Three outcomes, and collapsing any two of them would hide a real finding:

    - a panel meeting G4 -- no degradation;
    - a panel that solved but is not diverse enough -- reported with the
      distance it actually reached, NEVER by lowering `G4_MIN_PAIRWISE_DISTANCE`,
      because a G4 failure invalidates a product decision rather than a technical
      one;
    - nothing at all, when every weight vector was infeasible or refused. Then
      the unsteered solve stands alone and the report says there is no gallery.
    """
    gallery, solved = sweep_designs(space, steps=steps, k=k)
    picks = [point.cds for point in gallery.picks]
    if picks:
        if gallery.meets_g4:
            return gallery, picks, solved, None
        return (
            gallery,
            picks,
            solved,
            f"the {len(picks)}-candidate panel does not meet gate G4: its minimum "
            f"pairwise codon distance is {gallery.min_pairwise_distance:.1%}, below "
            f"the 15% at which candidates are genuinely different designs rather "
            f"than one design shown several times",
        )
    fallback = space.solve(None)
    if fallback is None:
        raise DesignError(
            "no candidate survived: every weight vector in the sweep and the "
            "unsteered solve were either infeasible or refused by the independent "
            "validator"
        )
    solved[fallback.cds] = fallback
    return (
        None,
        [fallback.cds],
        solved,
        "single candidate only: no weight vector in the sweep produced a design, "
        "so there is no gallery to choose from",
    )


def _degradations(
    rules: RuleSet,
    carried: Sequence[str],
    *,
    screen: BiosecurityVerdict,
    gallery: str | None,
    usage: str | None,
    native: str | None,
    native_supplied: bool,
    nulls: Nulls,
) -> tuple[str, ...]:
    """Every reason this run is not a complete evaluation, stated.

    A set the report renders and a test pins by equality: a new silent
    degradation must fail that test rather than slip in unremarked.
    `QcReport.is_complete` reads this tuple, so an entry here is what stops a
    partial run from looking finished -- and an empty tuple is what finally lets
    a genuinely complete one say so. The skeleton's three unconditional entries
    ("ranking not computed", "no gallery", "screening did not run") are now each
    conditional on the thing actually being absent, which is the whole point of
    this increment.
    """
    out: list[str] = []
    fold = degradation_reason()
    if fold is not None:
        out.append(fold)
    if screen.status != "clear":
        detail = f" ({screen.detail})" if screen.detail else ""
        out.append(f"protein-level biosecurity screening: {screen.status}{detail}")
    for optional in (gallery, usage, native):
        if optional is not None:
            out.append(optional)
    if native is None and not native_supplied:
        out.append(
            "no native baseline: the caller supplied no wild-type CDS, and BT5 will "
            "not manufacture one by back-translation. 'Do not optimize' is frequently "
            "the right answer and cannot be evaluated here."
        )
    for spec_id, reason in sorted(nulls.unavailable.items()):
        out.append(f"objective {spec_id} not ranked: {reason}")
    for spec_id in rules.unrunnable:
        out.append(
            f"rule {spec_id} not run: its thresholds are calibrated against a folding "
            f"engine that is not available"
        )
    for motif in carried:
        out.append(f"forbidden motif {motif} carried by the backbone, excluded from enforcement")
    return tuple(out)
