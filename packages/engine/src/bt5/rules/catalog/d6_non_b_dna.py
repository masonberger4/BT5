"""D6 -- G-quadruplexes and telomere repeats, weighted by whether they actually fold.

The reason this rule has a per-modality enforcement rather than one level is
that the evidence genuinely points in opposite directions depending on the host,
and both directions are well measured:

  E. coli      G4s fold, stall replication forks and impair growth. Mutation
               rates across G4 variants span 5.5e-5 to 2.7e-10 per cell per
               generation, with up to 8-fold dependence on orientation relative
               to the fork.
  mammalian    rG4s are globally UNFOLDED in cells -- median folding score 0.06.
               The sequence is there; the structure mostly is not.

One frozen level would be wrong in one of those two jobs, so `enforcement_for`
escalates to HARD_REPAIR on the bacterial modality and leaves it SOFT elsewhere.
The two citations carry opposite `sign` values for the same reason: a single URL
would make the badge dishonest on precisely the rule where the disagreement is
the point.

WHAT THIS RULE DELIBERATELY DOES NOT COVER. Brief row 2.D6 also lists Z-DNA and
triplex/H-DNA. They are not here, because a Spec carries ONE evidence badge and
theirs is not the same: Z-DNA's ranking as the most destabilising non-B
structure rests on a single 2004 study (evidence C). Folding them into a rule
badged EVIDENCE_BACKED would launder a weak claim through a strong one, which is
the failure the badge exists to prevent. They belong in their own rule, badged
`contested`.

Honesty about the detector: the regex misses ~37% of experimentally detected
rG4s. It is a screen, not a census, and `n_evaluated` reports what was scanned
so the report can say so.
"""

from __future__ import annotations

import re
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
from bt5.core.types import Construct, Interval, Strand, reverse_complement

#: The canonical four-tract G4 motif. Loops 1-7 nt, four runs of 3+ G.
G4_PATTERN = re.compile(r"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}")

#: Telomere-like tandem repeats, two units or more.
TELOMERE_PATTERNS: Mapping[str, re.Pattern[str]] = {
    "vertebrate_TTAGGG": re.compile(r"(?:TTAGGG){2,}"),
    "plant_TTTAGGG": re.compile(r"(?:TTTAGGG){2,}"),
    "ciliate_TTGGGG": re.compile(r"(?:TTGGGG){2,}"),
    "ciliate_TTTTGGGG": re.compile(r"(?:TTTTGGGG){2,}"),
}

#: G4Hunter: window 25, |score| >= 1.2 flags, >= 1.5 is severe.
G4HUNTER_WINDOW = 25
G4HUNTER_FLAG = 1.2
G4HUNTER_SEVERE = 1.5

#: Longest motif this rule can match, for the circular scan overlap. Four runs
#: capped at 12 G plus three 7 nt loops is comfortably inside this.
MAX_MOTIF = 100


def g4hunter_score(seq: str) -> float:
    """Mean G4Hunter score: each base scores by the length of its own G or C run.

    G-runs count positive and C-runs negative, so a strongly negative score means
    a G4 on the OTHER strand. Callers take the absolute value; the sign says
    which strand.
    """
    if not seq:
        return 0.0
    scores: list[int] = []
    i = 0
    while i < len(seq):
        base = seq[i]
        run = 1
        while i + run < len(seq) and seq[i + run] == base:
            run += 1
        value = min(run, 4)
        if base == "G":
            scores.extend([value] * run)
        elif base == "C":
            scores.extend([-value] * run)
        else:
            scores.extend([0] * run)
        i += run
    return sum(scores) / len(scores)


def peak_g4hunter(seq: str, window: int = G4HUNTER_WINDOW) -> float:
    """The strongest |score| over any window. 0.0 when the sequence is shorter."""
    if len(seq) < window:
        return 0.0
    return max(abs(g4hunter_score(seq[i : i + window])) for i in range(len(seq) - window + 1))


@register
class NonBDna:
    id: ClassVar[str] = "d6_non_b_dna"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "G-quadruplexes and telomere repeats"
    #: The FLOOR. `enforcement_for` escalates on the bacterial modality, where
    #: these structures are measured to fold.
    enforcement: ClassVar[Enforcement] = Enforcement.SOFT
    evidence: ClassVar[Evidence] = Evidence.EVIDENCE_BACKED
    direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
    unit: ClassVar[str] = "weighted motifs"
    band: ClassVar[tuple[float, float] | None] = None
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "G4s fold in E. coli and impair growth: mutation rates span 5.5e-5 to "
            "2.7e-10 per cell per generation across G4 variants, up to 8-fold dependent "
            "on orientation relative to the replication fork",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC10530614/",
            2023,
            sign="supports",
        ),
        Citation(
            "rG4s are globally UNFOLDED in mammalian cells, median folding score 0.06 -- "
            "the sequence is present and the structure mostly is not, which is why this "
            "rule is soft rather than hard in a mammalian context",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC5367264/",
            2017,
            sign="refutes",
        ),
        Citation(
            "G4Hunter: window 25, |score| >= 1.2 flags and >= 1.5 is severe; the "
            "four-tract regex misses ~37% of experimentally detected rG4s",
            "https://academic.oup.com/nar/article/44/4/1746/1852481",
            2016,
            sign="qualifies",
        ),
    )
    last_verified: ClassVar[str] = "2026-08-28"
    weight_provenance: ClassVar[str] = (
        "Moderate, and deliberately NOT high in the default (mammalian) case. The "
        "sequence motif is common and the structure is mostly absent in mammalian "
        "cells, so weighting it heavily would spend the sequence's freedom removing "
        "motifs that do not fold -- and every base spent there is a base not spent on "
        "repeats, which are the measured top predictor of synthesis failure. Where the "
        "structures do fold, enforcement_for escalates out of the weighted sum instead "
        "of relying on weight."
    )
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.35
    steering_weight: ClassVar[float] = 0.2
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.MOTIF_LEN_MINUS_1
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "cheap"
    conflicts_with: ClassVar[tuple[str, ...]] = ()
    brief_ref: ClassVar[str] = "2.D6"
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "g4hunter_flag": {"type": "number", "default": G4HUNTER_FLAG, "minimum": 0.0},
            "telomeres": {
                "type": "boolean",
                "default": True,
                "description": "Also report telomere-like tandem repeats.",
            },
        },
    }

    def __init__(self, g4hunter_flag: float = G4HUNTER_FLAG, telomeres: bool = True) -> None:
        if g4hunter_flag <= 0:
            raise ValueError(f"g4hunter_flag must be positive, got {g4hunter_flag}")
        self.g4hunter_flag = g4hunter_flag
        self.telomeres = telomeres

    def gate(self, slot: ContextSlot) -> bool:
        # An IVT mRNA is not replicated, so the fork-stalling argument does not
        # apply; the template plasmid is covered by that construct's own
        # propagation slot.
        return slot.modality is not Modality.IVT_MRNA

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        """Hard where the structures are measured to fold, soft where they are
        measured not to. One frozen level would be wrong in one of the two."""
        if slot.modality is Modality.BACTERIAL_EXPRESSION:
            return Enforcement.HARD_REPAIR
        return Enforcement.SOFT

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """None: a G4 is a variable-length pattern with 1-7 nt loops and
        unbounded G-runs, so it is not a finite motif set the automaton can
        make unreachable. Enumerating it would be thousands of patterns for a
        guarantee Tier-B repair gives directly."""
        return None

    def _hits(self, c: Construct) -> list[tuple[str, Interval, float, str]]:
        """(label, interval, g4hunter peak, class) on both strands, wrap-aware."""
        n = c.length
        found: list[tuple[str, Interval, float, str]] = []
        strands: tuple[Strand, ...] = (1, -1)
        for strand in strands:
            seq = c.sequence if strand == 1 else reverse_complement(c.sequence)
            scan = seq + seq[:MAX_MOTIF] if c.is_circular else seq

            patterns: list[tuple[str, re.Pattern[str], str]] = [("G4", G4_PATTERN, "quadruplex")]
            if self.telomeres:
                patterns += [(k, v, "telomere") for k, v in TELOMERE_PATTERNS.items()]

            for label, pattern, kind in patterns:
                for match in pattern.finditer(scan):
                    if match.start() >= n:
                        continue  # the wrap overlap, already reported at its real position
                    lo = match.start() if strand == 1 else n - match.end()
                    iv = Interval(max(0, lo), max(0, lo) + (match.end() - match.start()), strand)
                    found.append((label, iv, peak_g4hunter(match.group()), kind))
        return found

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        breaches: list[Breach] = []
        for label, iv, peak, kind in self._hits(c):
            severe = peak >= G4HUNTER_SEVERE
            flagged = peak >= self.g4hunter_flag
            breaches.append(
                Breach(
                    spec_id=self.id,
                    interval=iv,
                    # A motif the sequence-level regex finds but G4Hunter scores
                    # low is reported at low magnitude rather than dropped: the
                    # two detectors disagree often and neither is a census.
                    magnitude=2.0 if severe else (1.0 if flagged else 0.3),
                    message=(
                        f"{kind} motif {label} at {iv.start} on the "
                        f"{'+' if iv.strand == 1 else '-'} strand, "
                        f"G4Hunter peak {peak:.2f}"
                        + (" (severe)" if severe else "" if flagged else " (below the flag)")
                    ),
                    fixable_by_codon_choice=c.overlaps_editable(iv),
                    detail={"motif": label, "class": kind, "g4hunter": round(peak, 3)},
                )
            )
        return Evaluation(
            spec_id=self.id,
            passes=not any(b.magnitude >= 1.0 for b in breaches),
            raw_score=float(sum(b.magnitude for b in breaches)),
            breaches=tuple(breaches),
            # Both strands scanned; the report should say the regex is a screen
            # that misses ~37% of experimentally detected rG4s, not a census.
            n_evaluated=2 * c.length,
        )


__all__ = ["NonBDna", "g4hunter_score", "peak_g4hunter"]
