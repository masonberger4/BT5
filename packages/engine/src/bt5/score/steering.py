"""Codon scorers the gallery sweep can actually sweep over.

`build_gallery` needs a `solve(weights)` whose output MOVES as the weights move.
The solver ships two scorers -- `cai_scorer` and `repeat_breaking_scorer` -- and
sweeping a mixture of just those two gives a clustered gallery, because both are
dominated by the same relative-adaptiveness table `w`: at every weight the
argmin is very nearly the same codon. A gallery of near-duplicates is the exact
failure gate G4 exists to catch, and the honest fix is a genuinely orthogonal
axis rather than a lower threshold.

The axis added here is GC lean. Third positions are where synonymous choice
lives, and GC-rich against GC-poor synonyms differ at most wobble positions, so
leaning the two ways produces designs that are far apart in CODON space -- which
is the space `distance.codon_distance` measures and G4 gates on. Both leans stay
inside E2's band, because `optimize()` wraps whatever base score it is given in
`solve_with_gc_steering`; the lean chooses WHERE in the band a design sits, and
the band itself is still enforced by Tier B repair plus the independent
validator.

**Nothing here enforces anything.** These are steering weights in CLAUDE.md
3.5's sense, not enforcement: every candidate the sweep produces goes through
the same Tier-B repair and the same `verify_construct` refusal as any other
design. That is also why `REPEAT_STEERING_PENALTY` is 4.0 rather than
`repeat_breaking_scorer`'s 100.0. At 100 the repeat term swamps every other
term at any non-zero weight, so most of the simplex collapses onto one design
and the sweep stops sweeping. 4.0 still outranks the whole [0, 1] range of a
codon-preference difference, so a repeat-extending codon is never preferred to a
unique one at equal preference -- and the hard guarantee against repeats was
never this scorer's job in the first place. E5 and F1 are HARD_REPAIR; they hold
regardless of what weight the sweep picked.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from bt5.solver.reference import CodonScorer

#: The sweep's axes. Names, not rule ids: these are directions to pull the DP
#: in, and no one of them is any single rule's objective. `gc_lean_at` and
#: `gc_lean_gc` are two ends of ONE axis rather than a signed one, because
#: `simplex_weights` produces non-negative weights summing to 1 and cannot
#: express a sign.
SWEEP_AXES: tuple[str, ...] = (
    "codon_adaptation",
    "repeat_avoidance",
    "gc_lean_at",
    "gc_lean_gc",
)

#: See the module docstring: deliberately NOT `repeat_breaking_scorer`'s 100.0.
REPEAT_STEERING_PENALTY = 4.0

#: The k-mer whose repetition the steering term penalises. Same 9 as
#: `repeat_breaking_scorer`, and for the same reason: repetitive 9-mers per
#: 100 bp is one of the two highest-importance features in the published
#: synthesis-success model.
REPEAT_KMER = 9


def gc_fraction(codon: str) -> float:
    """Fraction of G/C bases in one codon, in [0, 1]."""
    if not codon:
        return 0.0
    return sum(1 for base in codon.upper() if base in "GC") / len(codon)


def blended_scorer(
    weights: Mapping[str, float],
    *,
    usage: Mapping[str, float],
    kmer: int = REPEAT_KMER,
    repeat_penalty: float = REPEAT_STEERING_PENALTY,
) -> CodonScorer:
    """One `CodonScorer` mixing the sweep axes at the given weights.

    A COST, like every `CodonScorer` in BT5: the DP minimises it, so a term is
    written as a penalty and a preference is written as a negative one. An axis
    absent from `weights` contributes nothing, so a caller sweeping two axes
    pays for two.
    """
    adaptation = float(weights.get("codon_adaptation", 0.0))
    repeats = float(weights.get("repeat_avoidance", 0.0))
    lean_at = float(weights.get("gc_lean_at", 0.0))
    lean_gc = float(weights.get("gc_lean_gc", 0.0))

    def score(_i: int, codon: str, prefix: str) -> float:
        cost = -adaptation * usage.get(codon, 0.0)
        if repeats:
            candidate = prefix + codon
            if len(candidate) >= kmer and candidate[-kmer:] in candidate[:-1]:
                cost += repeats * repeat_penalty
        if lean_at or lean_gc:
            gc = gc_fraction(codon)
            cost += lean_at * gc + lean_gc * (1.0 - gc)
        return cost

    return score


def live_axes(usage: Mapping[str, float], axes: Sequence[str] = SWEEP_AXES) -> tuple[str, ...]:
    """The axes that can actually move the answer, given this usage table.

    `codon_adaptation` reads `usage` and nothing else. With no host codon-usage
    table on file its term is identically zero, so a weight vector pushing only
    that axis solves to exactly the design an unsteered solve produces -- and
    `sweep` pays a full `optimize()`, Tier B repair included, to rediscover it.
    Measured at 500 aa on the reference backbone that is ~2.4 s bought for
    nothing, against a 10 s end-to-end budget.

    Dropping it is not a reduction in coverage. An axis whose contribution to the
    objective is the zero function cannot reach a point on the front that the
    remaining axes cannot; `Gallery.distinct` confirms it empirically, reporting
    the same count with and without.
    """
    if any(weight for weight in usage.values()):
        return tuple(axes)
    return tuple(axis for axis in axes if axis != "codon_adaptation")
