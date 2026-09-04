"""The QC report: what BT5 did, what it could not do, and what to do next.

The hard part of this file is what it refuses to say. All computable design
features together explain 5-31% of protein-level variance (mean ~14%, Cambray
2018, 244,000 designed sequences), and nine benchmarked commercial optimizers
were a coin flip against native sequence (Ranaghan 2021). So there is no
predicted expression level here, no titer, no yield, no fold-improvement --
only ranks, percentiles against a random-synonymous null, and statements about
what was and was not evaluated. A CI grep enforces the vocabulary; this module
is where the temptation to break it actually lives.

Three things the report must state out loud, because omitting any of them turns
an honest ranking into a dishonest one:

- The genetic code table, printed. NCBI table 12 reassigns CTG to Ser rather
  than Leu; a wrong table is a silently wrong protein no assay catches for
  months, and the report is the last place a human can catch it.
- Every objective that could NOT be evaluated. A scorecard missing its
  highest-weight term looks exactly like one where that term was never
  configured, and the difference is whether the ranking means anything.
- The native baseline. For homologous mammalian expression "use the native
  sequence" is frequently the right answer and no vendor tool will say so.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

from bt5.core.result import Candidate, Conflict, DesignResult, ObjectiveScore
from bt5.core.types import Provenance
from bt5.rules.vendors import DEFAULT_VENDOR, orderable

#: Mean error-free synthesis length, in bp, from each vendor's published
#: per-base fidelity. P(perfect clone) = exp(-L/E). These are VENDOR-ASSERTED
#: and drift, so they carry a verification date like any other vendor number.
#:
#: **Keyed on `vendors.PROFILES` keys, and that is load-bearing.** This table
#: used to be keyed on "twist" and "idt_eblocks" -- a THIRD vendor namespace,
#: overlapping the real registry by one key and by coincidence. It validated its
#: own keys, so it raised for `idt_gblocks`, which is BT5's own DEFAULT_VENDOR,
#: and accepted "twist", which names no orderable product. PR #53 merged two such
#: namespaces into one; this was the third. `test_no_fourth_vendor_namespace`
#: keeps it from reopening.
#:
#: A profile ABSENT from this table is not an error -- it is a vendor whose
#: fidelity BT5 has not got on file. See `screening_burden`.
ERROR_FREE_BP: dict[str, int] = {
    # https://www.twistbioscience.com/faq/gene-synthesis
    "twist_gene_fragment": 7500,
    # Inherited from the plain Gene Fragment: adapter-on is a checkout option on
    # the same product and the same synthesis process, not a different one. Same
    # within-a-vendor inheritance policy vendors.py states for GC bands and run
    # limits, and here it is inheritance within a single PRODUCT, which is safer
    # still.
    "twist_gene_fragment_adapter_on": 7500,
    # https://www.idtdna.com/pages/products/genes-and-gene-fragments/double-stranded-dna-fragments/eblocks-gene-fragments
    "idt_eblocks": 5000,
    # https://www.idtdna.com/pages/products/genes-and-gene-fragments/double-stranded-dna-fragments/gblocks-gene-fragments
    # "gBlocks Gene Fragments are double-stranded DNA fragments 125-3000 bp in
    # length with a median error rate of less than 1:5000." Retrieved 2026-09-04.
    #
    # INDEPENDENTLY PUBLISHED FOR gBLOCKS, NOT INHERITED FROM eBLOCKS -- and that
    # distinction is the entire reason #56 was an issue rather than a
    # two-character edit. The two products genuinely carry the same figure, and
    # IDT states both side by side in one comparison table on this page, so the
    # coincidence is theirs and not ours. The citation above is the gBlocks page;
    # it must never be re-pointed at the eBlocks URL, because a fidelity number
    # comes out of a manufacturing process and these are two processes at one
    # company. `test_the_vendor_changes_the_answer`: "reporting one vendor's
    # number under another's name is the kind of quiet error a lab pays for in
    # plates."
    #
    # NOT gBlocks HiFi, which is a SEPARATE product (1000-3000 bp, NGS
    # sequence-verified, "median error rate of less than 1:12,000" on the same
    # page). `vendors.py` gives idt_gblocks length_bp=(125, 3000), which is the
    # standard product's range and is how the two are told apart here. A HiFi
    # profile would need its own key and its own 12000.
    #
    # MEDIAN, AND A BOUND. IDT publishes "less than 1:5000" as a median; the
    # model below wants a mean rate. Reading one as the other is an
    # interpretation -- the same one the eBlocks entry already makes -- and the
    # "less than" puts E = 5000 at the CONSERVATIVE end, which over-estimates
    # colonies to pick. That is the safe direction for a bench instruction.
    "idt_gblocks": 5000,
}
ERROR_FREE_LAST_VERIFIED = "2026-09-04"

#: How sure the user wants to be of having at least one perfect clone.
DEFAULT_CONFIDENCE = 0.95


@dataclass(frozen=True, slots=True)
class ScreeningBurden:
    """How many colonies to pick, which is the decision the number is for.

    P(perfect clone) alone is a statistic. "Pick 14 colonies" is an instruction,
    and it is the only form of this number anyone acts on.
    """

    vendor: str
    length_bp: int
    error_free_bp: int
    p_perfect: float
    confidence: float
    colonies_to_pick: int
    last_verified: str = ERROR_FREE_LAST_VERIFIED


def screening_burden(
    length_bp: int,
    *,
    vendor: str = DEFAULT_VENDOR,
    confidence: float = DEFAULT_CONFIDENCE,
) -> ScreeningBurden | None:
    """Colonies to pick for `confidence` of at least one error-free clone.

    P(perfect) = exp(-L/E); n = ceil(ln(1-confidence) / ln(1-P)).

    Returns None -- not a raise, and not a guess -- when `vendor` is a real
    orderable configuration whose error-free length BT5 has not got on file.
    Two failures live here and they are NOT the same failure:

    - "acme" is not a vendor anyone can order from, and `none` is "no vendor
      chosen". Both are caller bugs and both raise.
    - `idt_gblocks` IS orderable, and is BT5's own default, and no published
      fidelity figure exists for it in this repo. That is a gap in BT5's data,
      not an error in the request, and the honest answer is "I do not know how
      many colonies to pick" rather than eBlocks' number wearing gBlocks' name.

    Collapsing the two is how the report ends up quietly authoritative about a
    number nobody measured. `build_report` turns the None into a degradation, so
    the absence is stated rather than merely missing.
    """
    if length_bp <= 0:
        raise ValueError(f"length must be positive, got {length_bp}")
    if not 0.0 < confidence < 1.0:
        raise ValueError(f"confidence must be in (0, 1), got {confidence}")
    # THREE outcomes, and the registry draws two of the lines for us. `orderable`
    # raises for a key nobody sells ("acme") and for `none`, which is "no vendor
    # chosen" -- and how many colonies to pick is a question only a vendor can
    # answer, so shrugging at it would be the guessed-vendor bug wearing a shrug.
    orderable(vendor)
    e = ERROR_FREE_BP.get(vendor)
    if e is None:
        return None

    p = math.exp(-length_bp / e)
    # p == 1.0 only for a construct short enough that the log below divides by
    # zero. One colony is the honest answer there, not an infinity.
    colonies = 1 if p >= 1.0 else math.ceil(math.log(1.0 - confidence) / math.log(1.0 - p))
    return ScreeningBurden(
        vendor=vendor,
        length_bp=length_bp,
        error_free_bp=e,
        p_perfect=p,
        confidence=confidence,
        colonies_to_pick=max(1, colonies),
    )


@dataclass(frozen=True, slots=True)
class QcReport:
    """Everything the run can honestly say, as data rather than prose.

    Structured because the server and the UI both consume it, and because a
    report assembled as a string is a report nothing can test.
    """

    design_hash: str
    translation_table_id: int
    construct_length_bp: int
    candidates: int
    #: Objectives that WERE evaluated, and where this design sits against the
    #: null. Never an absolute number.
    scored: tuple[ObjectiveScore, ...] = ()
    #: Objectives that could not be evaluated, each with its reason.
    unavailable: tuple[ObjectiveScore, ...] = ()
    conflicts: tuple[Conflict, ...] = ()
    #: Findings no codon can fix -- a uAUG in the user's own 5'UTR, an AAV size
    #: overflow. Separated because the action is different.
    advisories: tuple[str, ...] = ()
    degradations: tuple[str, ...] = ()
    strain_protocol: tuple[str, ...] = ()
    burden: ScreeningBurden | None = None
    native_baseline_available: bool = False
    preset_id: str = ""
    engine_versions: tuple[tuple[str, str], ...] = ()

    @property
    def is_complete(self) -> bool:
        """Did every configured objective actually get evaluated?"""
        return not self.unavailable and not self.degradations


def build_report(
    result: DesignResult,
    candidate: Candidate,
    *,
    translation_table_id: int,
    preset_id: str = "",
    vendor: str = DEFAULT_VENDOR,
    advisories: Sequence[str] = (),
    strain_protocol: Sequence[str] = (),
    extra_degradations: Sequence[str] = (),
) -> QcReport:
    """Assemble the report for one candidate of a design run.

    `translation_table_id` is passed explicitly and has no default. Defaulting
    it to 1 here would be the same silent-wrong-protein bug the contract refuses
    everywhere else, one layer further out, and the report is the last place a
    human sees it before ordering DNA.
    """
    provenance: Provenance | None = result.provenance
    degradations = tuple(provenance.degradations) if provenance else ()
    engines = tuple(sorted(provenance.engine_versions.items())) if provenance else ()

    # A vendor with no error-free length on file is REPORTED as unknown, never
    # answered with a neighbour's number. `is_complete` reads `degradations`, so
    # this also stops a report missing its screening line from looking complete.
    burden = screening_burden(len(candidate.cds), vendor=vendor)
    if burden is None:
        degradations += (
            f"screening burden unavailable: no published error-free length on file "
            f"for {vendor}, so BT5 cannot say how many colonies to pick",
        )

    return QcReport(
        design_hash=candidate.design_hash,
        translation_table_id=translation_table_id,
        construct_length_bp=candidate.construct.length,
        candidates=len(result.candidates),
        scored=candidate.scorecard.available,
        unavailable=candidate.scorecard.unavailable,
        conflicts=result.conflicts,
        advisories=tuple(advisories),
        degradations=degradations + tuple(extra_degradations),
        strain_protocol=tuple(strain_protocol),
        burden=burden,
        native_baseline_available=result.native_baseline is not None,
        preset_id=preset_id,
        engine_versions=engines,
    )


def render(report: QcReport) -> str:
    """A plain-text rendering, for the CLI and the GenBank note.

    Deliberately boring. Every line is either a fact about this run or a
    statement about what was not done; there is nowhere for a predicted number
    to hide.
    """
    lines: list[str] = [
        f"BT5 design {report.design_hash}",
        f"  genetic code       NCBI translation table {report.translation_table_id}",
        f"  construct          {report.construct_length_bp} bp",
        f"  candidates         {report.candidates}",
    ]
    if report.preset_id:
        lines.append(f"  preset             {report.preset_id}")
    for name, version in report.engine_versions:
        lines.append(f"  engine             {name} {version}")

    lines.append("")
    lines.append("Objectives (percentile against a random-synonymous null; not a prediction)")
    if report.scored:
        for score in report.scored:
            lines.append(
                f"  {score.spec_id:<28} {score.percentile:>6.1%}  "
                f"(raw {score.raw:g} {score.unit}, n={score.null_n})"
            )
    else:
        lines.append("  none evaluated")

    if report.unavailable:
        lines.append("")
        lines.append("NOT evaluated -- the ranking above does not account for these")
        for score in report.unavailable:
            lines.append(f"  {score.spec_id:<28} {score.unavailable_reason}")

    if report.conflicts:
        lines.append("")
        lines.append("Conflicts")
        for conflict in report.conflicts:
            lines.append(
                f"  [{conflict.interval.start}, {conflict.interval.end}) "
                f"{conflict.kind}: {', '.join(conflict.spec_ids)} "
                f"(binding: {conflict.binding_spec_id})"
            )

    if report.advisories:
        lines.append("")
        lines.append("Not fixable by codon choice")
        for advisory in report.advisories:
            lines.append(f"  {advisory}")

    if report.degradations:
        lines.append("")
        lines.append("Degradations")
        for degradation in report.degradations:
            lines.append(f"  {degradation}")

    if report.burden is not None:
        b = report.burden
        lines.append("")
        lines.append(
            f"Synthesis screening ({b.vendor}, error-free length {b.error_free_bp} bp, "
            f"verified {b.last_verified})"
        )
        lines.append(
            f"  {b.length_bp} bp: P(perfect clone) {b.p_perfect:.1%}; "
            f"pick {b.colonies_to_pick} colonies for {b.confidence:.0%} confidence"
        )

    if report.strain_protocol:
        lines.append("")
        lines.append("Propagation")
        for step in report.strain_protocol:
            lines.append(f"  {step}")

    lines.append("")
    lines.append(
        "BT5 reports ranks and percentiles, never a predicted expression level, "
        "titer or yield. All computable design features together explain 5-31% "
        "of protein-level variance."
    )
    if report.native_baseline_available:
        lines.append(
            "The native sequence is included as a candidate. For homologous "
            "mammalian expression it is frequently the right answer."
        )
    return "\n".join(lines)
