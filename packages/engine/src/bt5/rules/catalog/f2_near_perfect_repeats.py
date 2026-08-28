"""F2 -- near-perfect repeats: the ones F1 cannot see.

F1 finds EXACT repeats. The repeats that actually appear in the constructs BT5
serves are usually not exact: two 2A peptides that differ at a few wobble
positions, two copies of a linker with one substitution, a duplicated domain
that drifted. They recombine anyway -- single-strand annealing does not require
perfect identity -- and an exact-match scan reports nothing at all.

So this rule seeds on exact matches and EXTENDS THROUGH mismatches, reporting
pairs at least 40 bp long and at least 90% identical. Pairs that come out
perfectly identical are dropped: those are F1's, and reporting them twice would
double every repeat in the panel.

WHAT THIS MISSES, MEASURED. Seed-and-extend can only find a near-perfect repeat
whose copies share SOME exact stretch. Pigeonhole gives the bound: 40 bp at 90%
identity is at most 4 mismatches, so some exact run of at least 40 // (4 + 1) = 8
must exist, and a seed of 8 would be complete at the stated threshold.

It is not the default, because complete-in-principle is not complete in fact.
`duplicates()` returns at most 200 pairs, longest first, and chance 8-mer matches
are all length 8 -- so they tie with the real short seeds and the cap drops them
arbitrarily. Measured on random sequence:

    5 kb    seed 8: 132 chance pairs      seed 10:   9
    10 kb   seed 8: 200 (cap saturated)   seed 10:  30
    20 kb   seed 8: 200 (cap saturated)   seed 10: 128

At the sizes this tool is for -- an 8-12 kb lentiviral transfer plasmid -- a seed
of 8 spends the entire budget on noise and crowds out the findings it was
supposed to guarantee. Seed 10 catches copies whose mismatches are 10 or more
apart, which is roughly 3 mismatches in 40 bp (92.5% identity), and stays well
inside the cap at 20 kb. The gap between 90% and 92.5% is the honest cost, the
seed is a parameter, and this docstring is the disclosure.

Alignment is UNGAPPED. An indel between two copies shifts the register and this
rule will report the two halves separately rather than one longer repeat. Gapped
alignment is a different piece of machinery and claiming it here would overstate
what the code does.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import ClassVar

from bt5.core.context import ContextSlot, DesignContext
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
from bt5.rules.exempt import both_arms_exempt, overlap_length

#: The brief's threshold: pairs at least this identical over at least this long.
MIN_IDENTITY = 0.90
MIN_LENGTH_BP = 40

#: Exact seed length. See the module docstring for why this is 10 and not the
#: pigeonhole-complete 8.
DEFAULT_SEED_BP = 10

#: Risk decays exponentially with separation; the brief puts "high" under 3 kb,
#: so that is the decay length -- at 3 kb the factor is 1/e.
SPACER_DECAY_BP = 3000

#: How far to try extending from a seed. Beyond this the pair is already a
#: finding at any threshold and the extra precision buys nothing.
MAX_EXTENT_BP = 600

#: Stop extending after this many consecutive mismatches: past it, identity
#: cannot recover within the extent that matters.
MAX_RUN_OF_MISMATCHES = 6

#: Cost of a mismatch during extension, in matches. Chosen so the score is
#: positive inside a 90%-identical repeat (+0.5 per base) and sharply negative
#: in random flanking sequence (-2.75 per base), which puts the maximum at the
#: repeat's boundary instead of at the identity floor.
MISMATCH_PENALTY = 4

MAX_FINDINGS = 200


def _at(seq: str, pos: int, n: int, circular: bool) -> str | None:
    """Base at `pos`, wrapping on a circular construct. None when off the end."""
    if circular:
        return seq[pos % n]
    return seq[pos] if 0 <= pos < n else None


def _extend(
    seq: str,
    a: int,
    b: int,
    seed: int,
    *,
    circular: bool,
    min_identity: float,
) -> tuple[int, int, int, int] | None:
    """Grow a seeded exact match through mismatches, ungapped.

    Returns (a_start, b_start, length, mismatches) for the maximum-scoring
    extent, or None when what it grows into does not meet `min_identity`.

    The extent is chosen by score, NOT by "longest that still clears the
    identity floor" -- see `grow` for why those differ and why it matters. The
    final identity is recomputed and checked rather than assumed: reporting a
    pair as 90% identical when it is not is the one thing this rule must not do.
    """
    n = len(seq)

    def grow(rightward: bool, length: int, mismatches: int) -> tuple[int, int, int]:
        """Walk one direction to the MAXIMUM-SCORING extent, not the longest.

        Score is +1 per match and -MISMATCH_PENALTY per mismatch, and the best
        score seen is what is kept. That is the whole difference between this
        and "extend while identity stays above the floor", which sounds
        equivalent and is not: extending while merely above the floor keeps
        bleeding into flanking sequence until enough mismatches accumulate to
        hit it, so every finding comes back at exactly the threshold identity
        and several bases too long. The reported identity would then be an
        artefact of the stopping rule rather than a property of the repeat.

        The penalty is what puts the maximum at the repeat's real boundary.
        Inside a 90%-identical repeat the expected score per base is
        0.9(+1) + 0.1(-4) = +0.5, so extension continues; in random flanking
        sequence it is 0.25(+1) + 0.75(-4) = -2.75, so it stops at once.
        """
        best = (0, length, mismatches)
        score = best_score = 0
        run = 0
        for step in range(1, MAX_EXTENT_BP):
            if rightward:
                pa, pb = a + seed + step - 1, b + seed + step - 1
            else:
                pa, pb = a - step, b - step
            x, y = _at(seq, pa, n, circular), _at(seq, pb, n, circular)
            if x is None or y is None:
                break
            length += 1
            if x == y:
                score += 1
                run = 0
            else:
                mismatches += 1
                score -= MISMATCH_PENALTY
                run += 1
                if run > MAX_RUN_OF_MISMATCHES:
                    break
            if score > best_score:
                best_score, best = score, (step, length, mismatches)
        return best

    # Right first, then left from the extent the right pass settled on, so the
    # two directions share one running count rather than each being locally
    # optimal and jointly out of band.
    _, length, mismatches = grow(True, seed, 0)
    left_steps, length, mismatches = grow(False, length, mismatches)

    if length <= 0 or (length - mismatches) / length < min_identity:
        return None
    a_start = (a - left_steps) % n if circular else a - left_steps
    b_start = (b - left_steps) % n if circular else b - left_steps
    return a_start, b_start, length, mismatches


def risk_score(length: int, spacer: int) -> float:
    """Monotone in length, exponentially decaying in spacer.

    A relative score for ordering findings, in arbitrary units. It is not a
    recombination rate and must never be printed as one -- the literature
    supports the SHAPE of this relationship and not a calibrated frequency.
    """
    return (length / MIN_LENGTH_BP) * math.exp(-spacer / SPACER_DECAY_BP)


@register
class NearPerfectRepeats:
    id: ClassVar[str] = "f2_near_perfect_repeats"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "Near-perfect direct repeats (>=90% identity over >=40 bp)"
    enforcement: ClassVar[Enforcement] = Enforcement.SOFT
    evidence: ClassVar[Evidence] = Evidence.EVIDENCE_BACKED
    direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
    unit: ClassVar[str] = "risk score (relative, not a rate)"
    band: ClassVar[tuple[float, float] | None] = None
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "Single-strand annealing and slipped-strand mispairing do not require "
            "perfect identity, and below ~200 bp the route is RecA-independent, so a "
            "recA- strain does not suppress it",
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5426353/",
            2017,
            sign="supports",
        ),
        Citation(
            "Recombination between imperfect repeats is strongly proximity-sensitive: "
            "inserting sequence between the copies suppresses it, which is what the "
            "exponential spacer term encodes",
            "https://www.pnas.org/doi/10.1073/pnas.111008398",
            2001,
            sign="supports",
        ),
        Citation(
            "Two copies of the same 2A peptide in one ORF is the canonical near-perfect "
            "repeat in this workflow; the literature requires DIFFERENT 2A peptides and "
            "<=85% nucleotide identity between them",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC5438344/",
            2017,
            sign="supports",
        ),
        Citation(
            "The specific 90%/40 bp cut and the 3 kb decay length are conventions from "
            "the design brief, not measured constants: the evidence fixes the ordering "
            "and the proximity dependence, not a threshold",
            "https://link.springer.com/article/10.1007/BF00290109",
            1983,
            sign="qualifies",
        ),
    )
    last_verified: ClassVar[str] = "2026-08-28"
    weight_provenance: ClassVar[str] = (
        "Substantial, and below F1's steering weight rather than above it. The "
        "mechanism is the same one F1 covers and the evidence for it is the same, but "
        "the THRESHOLD here is a convention (90% over 40 bp) rather than a measured "
        "constant, and the detector is a seed-and-extend approximation with a stated "
        "blind spot. Weighting it level with an exact-repeat finding would give a "
        "conventional cut the authority of a measured one. It is not lower still "
        "because the canonical case -- two similar 2A peptides in one ORF -- is "
        "invisible to every exact-match tool and is a documented failure."
    )
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.6
    steering_weight: ClassVar[float] = 0.5
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.PAIRED_SEGMENTS
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "moderate"
    conflicts_with: ClassVar[tuple[str, ...]] = ()
    brief_ref: ClassVar[str] = "2.F2"
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "min_identity": {
                "type": "number",
                "default": MIN_IDENTITY,
                "minimum": 0.5,
                "maximum": 1.0,
            },
            "min_length": {"type": "integer", "default": MIN_LENGTH_BP, "minimum": 12},
            "seed": {
                "type": "integer",
                "default": DEFAULT_SEED_BP,
                "minimum": 6,
                "description": (
                    "Exact seed length. 8 is pigeonhole-complete at 90%/40 bp but "
                    "saturates the pair budget with chance matches above ~10 kb."
                ),
            },
        },
    }

    def __init__(
        self,
        min_identity: float = MIN_IDENTITY,
        min_length: int = MIN_LENGTH_BP,
        seed: int = DEFAULT_SEED_BP,
    ) -> None:
        if not 0.5 <= min_identity <= 1.0:
            raise ValueError(
                f"min_identity {min_identity} outside [0.5, 1.0]; below half identity "
                f"two sequences are not repeats of each other in any useful sense"
            )
        if min_length < 12:
            raise ValueError(f"min_length {min_length} is short enough to occur by chance")
        if seed < 6:
            raise ValueError(f"seed {seed} would match everywhere; 4^6 is already only 4096")
        if seed > min_length:
            raise ValueError(
                f"seed {seed} exceeds min_length {min_length}: no seed could ever grow "
                f"into a reportable pair"
            )
        self.min_identity = min_identity
        self.min_length = min_length
        self.seed = seed

    def gate(self, slot: ContextSlot) -> bool:
        return True

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """None. Approximate identity between two distant spans is not decidable
        from a bounded suffix at all, let alone expressible as a motif set."""
        return None

    def _spacer(self, a_end: int, b_start: int, n: int, circular: bool) -> int:
        gap = b_start - a_end
        if not circular:
            return max(0, gap)
        return max(0, min(gap % n, (-gap) % n))

    def _already_reported(self, arm: Interval, kept: list[Interval], c: Construct) -> bool:
        """Has a previous seed already grown into this same physical repeat?

        Adjacent seeds inside one repeat each extend to nearly the same extent,
        differing by a base or two at the ends, so a key on exact coordinates
        deduplicates nothing -- one 60 bp repeat reported three times, at 197,
        198 and 200. Overlap is the right test, and generously: any substantial
        shared span means the same finding.
        """
        return any(
            overlap_length(arm, other, c.length, c.is_circular) >= 0.5 * arm.length
            for other in kept
        )

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        seq, n = c.sequence, c.length
        index = svc.kmer.of(c, self.seed)
        breaches: list[Breach] = []
        kept: list[Interval] = []

        for first, second in index.duplicates(self.seed):
            if len(breaches) >= MAX_FINDINGS:
                break
            grown = _extend(
                seq,
                first.start,
                second.start,
                first.length,
                circular=c.is_circular,
                min_identity=self.min_identity,
            )
            if grown is None:
                continue
            a_start, b_start, length, mismatches = grown
            if length < self.min_length:
                continue
            if mismatches == 0:
                continue  # perfectly identical: F1's finding, not this rule's

            arm_a = Interval(a_start, a_start + length)
            arm_b = Interval(b_start, b_start + length)
            if both_arms_exempt(c, arm_a, arm_b):
                continue
            if self._already_reported(arm_a, kept, c):
                continue
            kept.append(arm_a)

            identity = (length - mismatches) / length
            spacer = self._spacer(a_start + length, b_start, n, c.is_circular)
            score = risk_score(length, spacer)

            breaches.append(
                Breach(
                    spec_id=self.id,
                    interval=Interval(a_start, max(a_start + length, b_start + length)),
                    magnitude=score,
                    message=(
                        f"{length} bp repeat at {identity:.0%} identity "
                        f"({mismatches} mismatch{'es' if mismatches != 1 else ''}) "
                        f"between {a_start} and {b_start}, {spacer} bp apart. "
                        f"Exact-match scans do not see this. A recA- strain does not "
                        f"cover it either: below ~200 bp the route is RecA-independent."
                    ),
                    fixable_by_codon_choice=c.overlaps_editable(arm_a)
                    or c.overlaps_editable(arm_b),
                    detail={
                        "length": float(length),
                        "identity": round(identity, 4),
                        "mismatches": float(mismatches),
                        "spacer": float(spacer),
                        "risk_score": round(score, 4),
                    },
                )
            )

        return Evaluation(
            spec_id=self.id,
            passes=not breaches,
            raw_score=float(sum(b.magnitude for b in breaches)),
            breaches=tuple(breaches),
            n_evaluated=c.length,
        )
