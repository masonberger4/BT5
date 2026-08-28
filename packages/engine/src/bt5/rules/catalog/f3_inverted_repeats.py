"""F3 -- inverted repeats: cruciform extrusion, fork stalling, SbcCD cleavage.

An inverted repeat is NOT a direct repeat with a minus sign on it, and treating
the two as one rule is the mistake this file exists to avoid. They are different
mechanisms with different answers:

  direct repeat     lost by deletion (slipped-strand / SSA); codon choice fixes it
  inverted repeat   extrudes a cruciform, stalls the fork, is cleaved by SbcCD;
                    the answer is usually a strain and a temperature

So F3 does not reuse F1's (length, spacer) bands. Those are calibrated on the
deletion literature and applying them here would be a category error. This rule
bands on stem length and loop size, which is what the palindrome literature
actually measures -- `bt5/vector/kmers.py` says in as many words that banding
inverted repeats "needs its own citations. Until the rules lane does that pass,
this reports geometry and nothing else." This is that pass.

The scale that matters:

  stem >= 30 bp, or total palindrome >= 60 bp   hard
  perfect palindrome >= 150 bp                  do not build in a standard host
  stem 15-29 bp with loop <= 100 nt             warn
  IR < 100 bp                                   usually stable
  IR 100-1200 bp                                transiently unstable

Loop size matters because extrusion needs the arms to find each other: a tight
loop is a hairpin that forms readily, and a perfect palindrome -- arms abutting,
loop 0 -- is the most extrudable configuration there is.

AAV is the strict case and gets its own threshold. AAV-GPseq mapped truncation
hotspots EXACTLY to inverted repeats in the CMV enhancer, the CB promoter and
the EGFP ORF, so a stem that a plasmid tolerates can still truncate a packaged
genome.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import ClassVar

from bt5.core.context import ContextSlot, DesignContext, Modality
from bt5.core.registry import register
from bt5.core.services import Services
from bt5.core.spec import (
    Breach,
    Citation,
    Direction,
    Enforcement,
    Evaluation,
    Evidence,
    LatticeTerms,
    LocalizationPolicy,
    RepairPolicy,
)
from bt5.core.types import Construct, Interval
from bt5.rules.exempt import both_arms_exempt

#: Shortest stem worth reporting, paired with a tight loop.
WARN_STEM_BP = 15
#: Hard: a perfect stem at least this long.
HARD_STEM_BP = 30
#: Hard: arms plus loop at least this long.
HARD_TOTAL_BP = 60
#: "Do not build in a standard host" -- SbcCD destroys the replicon.
UNBUILDABLE_BP = 150
#: AAV is stricter: truncation hotspots map exactly to inverted repeats.
AAV_STEM_BP = 20
#: Above this loop the arms are too far apart for the warn band to mean much.
MAX_LOOP_BP = 100

#: A palindromic region must not flood the report.
MAX_FINDINGS = 200


def _loop(first: Interval, second: Interval, length: int, circular: bool) -> int:
    """Bases between the two arms.

    Wrap-aware: on a circular construct the 3' arm can sit past the origin, and
    a negative or whole-plasmid loop would put the most extrudable hairpins in
    the least alarming band.
    """
    gap = second.start - first.end
    if gap >= 0:
        return gap
    return gap + length if circular else 0


@register
class InvertedRepeats:
    id: ClassVar[str] = "f3_inverted_repeats"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "Inverted repeats and palindromes (cruciform liability)"
    enforcement: ClassVar[Enforcement] = Enforcement.HARD_CHECK
    evidence: ClassVar[Evidence] = Evidence.EVIDENCE_BACKED
    direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
    unit: ClassVar[str] = "weighted inverted repeats"
    band: ClassVar[tuple[float, float] | None] = None
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "Long palindromes are cleaved by SbcCD and destroy the replicon; IRs under "
            "100 bp are usually stable, 100-1200 bp transiently unstable. Rank by the "
            "nearest-neighbour dG of the basal 20 bp of the stem (Sinden's determinant)",
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5426353/",
            2017,
            sign="supports",
        ),
        Citation(
            "AAV-GPseq mapped truncation hotspots EXACTLY to inverted repeats in the CMV "
            "enhancer, CB promoter and EGFP ORF -- a stem a plasmid tolerates can still "
            "truncate a packaged genome",
            "https://www.cell.com/molecular-therapy-family/advances/fulltext/S2329-0501(20)30156-X",
            2020,
            sign="supports",
        ),
        Citation(
            "Stbl2 to Stbl3 alone rescued an HIV vector lost entirely in 0.5 L Stbl2 "
            "cultures -- the strain and temperature protocol is the answer for an "
            "unavoidable IR, which is why this rule is HARD_CHECK and not HARD_REPAIR",
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3563744/",
            2013,
            sign="supports",
        ),
    )
    last_verified: ClassVar[str] = "2026-08-28"
    weight_provenance: ClassVar[str] = ""  # hard rule; never weighted
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.0
    #: HARD_CHECK is reported, never chased by the solver, so steering is a
    #: category error here: the dominant cases (ITRs, shRNA hairpins, an IR in
    #: the user's own backbone) cannot be steered away by codon choice at all.
    steering_weight: ClassVar[float] = 0.0
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.PAIRED_SEGMENTS
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "moderate"
    conflicts_with: ClassVar[tuple[str, ...]] = ()
    brief_ref: ClassVar[str] = "2.F3"
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "min_stem": {"type": "integer", "default": WARN_STEM_BP, "minimum": 8},
            "max_loop": {"type": "integer", "default": MAX_LOOP_BP, "minimum": 0},
        },
    }

    def __init__(self, min_stem: int = WARN_STEM_BP, max_loop: int = MAX_LOOP_BP) -> None:
        if min_stem < 8:
            raise ValueError(
                f"min_stem {min_stem} is short enough to occur by chance everywhere; "
                f"an 8 bp stem is noise, not a finding"
            )
        if max_loop < 0:
            raise ValueError(f"max_loop must be >= 0, got {max_loop}")
        self.min_stem = min_stem
        self.max_loop = max_loop

    def gate(self, slot: ContextSlot) -> bool:
        # Every construct passes through a cloning host, and an IVT template is
        # still a plasmid before it is a transcript.
        return True

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """None, and not for the usual reason. An inverted repeat is not a fixed
        motif set, but more importantly this rule is HARD_CHECK: the finding is
        real and mostly NOT fixable by codon choice, so there is nothing for the
        automaton to make unreachable."""
        return None

    def _stem_threshold(self, ctx: DesignContext) -> int:
        """AAV truncates on stems a plasmid tolerates."""
        if any(s.modality is Modality.AAV for s in ctx.active_slots):
            return AAV_STEM_BP
        return HARD_STEM_BP

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        index = svc.kmer.of(c, self.min_stem)
        hard_stem = self._stem_threshold(ctx)
        breaches: list[Breach] = []
        worst_stem = 0

        for first, second in index.revcomp_pairs(self.min_stem, self.max_loop):
            if len(breaches) >= MAX_FINDINGS:
                break
            if both_arms_exempt(c, first, second):
                continue
            stem = first.length
            loop = _loop(first, second, c.length, c.is_circular)
            total = 2 * stem + loop
            worst_stem = max(worst_stem, stem)

            if stem >= UNBUILDABLE_BP:
                severity, magnitude = "unbuildable", 4.0
                advice = (
                    "Do not build this in a standard host: SbcCD cleaves long "
                    "palindromes and destroys the replicon. Use a sbcC-deficient "
                    "strain at 42 C."
                )
            elif stem >= hard_stem or total >= HARD_TOTAL_BP:
                severity, magnitude = "hard", 2.0
                advice = (
                    "Propagate in a recA- strain at 30 C; for a palindrome this long "
                    "prefer a sbcC-deficient host at 42 C."
                )
            else:
                severity, magnitude = "warn", 0.5
                advice = "Usually stable below 100 bp, but worth knowing about."

            span = Interval(first.start, max(first.end, second.end))
            breaches.append(
                Breach(
                    spec_id=self.id,
                    interval=span,
                    magnitude=magnitude,
                    message=(
                        f"inverted repeat: {stem} bp stem with a {loop} bp loop "
                        f"({total} bp total) at {first.start}"
                        + (" -- PERFECT PALINDROME, the most extrudable form" if loop == 0 else "")
                        + f"; {severity}. {advice}"
                    ),
                    # An IR whose arms lie in the CDS can be broken by recoding.
                    # One in the user's backbone cannot, and that is the common
                    # case -- which is why the class is HARD_CHECK.
                    fixable_by_codon_choice=c.overlaps_editable(first)
                    and c.overlaps_editable(second),
                    detail={
                        "stem": float(stem),
                        "loop": float(loop),
                        "total": float(total),
                        "severity": severity,
                        "perfect_palindrome": "yes" if loop == 0 else "no",
                    },
                )
            )

        return Evaluation(
            spec_id=self.id,
            passes=worst_stem < hard_stem,
            raw_score=float(sum(b.magnitude for b in breaches)),
            breaches=tuple(breaches),
            n_evaluated=c.length,
        )
