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

#: Mean error-free synthesis length, in bp, from each vendor's published
#: per-base fidelity. P(perfect clone) = exp(-L/E). These are VENDOR-ASSERTED
#: and drift, so they carry a verification date like any other vendor number.
ERROR_FREE_BP: dict[str, int] = {
    "twist": 7500,
    "idt_eblocks": 5000,
}
ERROR_FREE_LAST_VERIFIED = "2026-08-28"

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
    vendor: str = "idt_eblocks",
    confidence: float = DEFAULT_CONFIDENCE,
) -> ScreeningBurden:
    """Colonies to pick for `confidence` of at least one error-free clone.

    P(perfect) = exp(-L/E); n = ceil(ln(1-confidence) / ln(1-P)).
    """
    if length_bp <= 0:
        raise ValueError(f"length must be positive, got {length_bp}")
    if not 0.0 < confidence < 1.0:
        raise ValueError(f"confidence must be in (0, 1), got {confidence}")
    try:
        e = ERROR_FREE_BP[vendor]
    except KeyError:
        raise ValueError(
            f"no error-free length on file for {vendor!r}; have {sorted(ERROR_FREE_BP)}"
        ) from None

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
    vendor: str = "idt_eblocks",
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
        burden=screening_burden(len(candidate.cds), vendor=vendor),
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
