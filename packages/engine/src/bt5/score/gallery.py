"""Choosing 3-8 candidates a user would actually see as different options.

This is gate G4, and its failure invalidates a PRODUCT decision rather than a
technical one: if the gallery cannot produce genuinely distinct designs, the
gallery is not a feature and the UI built on it is a lie.

The selection rule is not obvious, and the obvious one is wrong. Das & Dennis
(1997) proved two things about weighted sums that together rule out the natural
approach: non-convex regions of the Pareto front are UNREACHABLE by any weight
vector, and -- worse for a gallery -- an evenly spaced set of weights produces an
UNEVENLY spaced set of solutions, because the weight tracks the local slope of
the front rather than position along it. So sweeping weights evenly and showing
the results gives a clustered gallery that looks diverse only in the weights.

  https://link.springer.com/article/10.1007/BF01197559

What works instead, and what this implements: sweep densely, then select on
distance in SEQUENCE space by greedy max-min. The sweep's job is to generate
candidates; the selection's job is to spread them. Neither can do the other's.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass

from bt5.score.distance import codon_distance, pairwise_minimum

#: docs/PLAN.md: a panel of 3-8 genuinely different candidates, not one "optimal".
MIN_GALLERY = 3
MAX_GALLERY = 8

#: Gate G4's threshold: pairwise codon distance at or above this, or the sweep
#: is not producing a gallery and the plan says to re-plan around
#: epsilon-constraint enumeration BEFORE a UI is built on it.
G4_MIN_PAIRWISE_DISTANCE = 0.15


@dataclass(frozen=True, slots=True)
class SweepPoint:
    """One solved design and the weights that produced it."""

    weights: Mapping[str, float]
    cds: str


@dataclass(frozen=True, slots=True)
class Gallery:
    """The selected panel, and the evidence that it is one."""

    picks: tuple[SweepPoint, ...]
    swept: int
    distinct: int
    min_pairwise_distance: float

    @property
    def meets_g4(self) -> bool:
        return (
            len(self.picks) >= MIN_GALLERY
            and self.min_pairwise_distance >= G4_MIN_PAIRWISE_DISTANCE
        )


def simplex_weights(objectives: Sequence[str], steps: int) -> tuple[dict[str, float], ...]:
    """Dense lattice over the weight simplex, summing to 1.

    Dense on purpose. These weights are the SAMPLING of the front, never the
    selection from it -- picking the gallery by taking evenly spaced weights is
    exactly the mistake Das & Dennis rules out, so the only requirement here is
    coverage.
    """
    if not objectives:
        raise ValueError("cannot sweep with no objectives")
    if steps < 1:
        raise ValueError(f"steps must be positive, got {steps}")

    def compositions(k: int, total: int) -> list[list[int]]:
        if k == 1:
            return [[total]]
        return [[i, *rest] for i in range(total + 1) for rest in compositions(k - 1, total - i)]

    return tuple(
        dict(zip(objectives, [c / steps for c in combo], strict=True))
        for combo in compositions(len(objectives), steps)
    )


def sweep(
    objectives: Sequence[str],
    solve: Callable[[Mapping[str, float]], str | None],
    *,
    steps: int = 8,
) -> tuple[SweepPoint, ...]:
    """Solve at every lattice point. A weight vector that cannot be solved is
    dropped rather than substituted -- an infeasible corner of the simplex is a
    real answer about the constraint set, not a gap to fill in."""
    out: list[SweepPoint] = []
    for w in simplex_weights(objectives, steps):
        cds = solve(w)
        if cds is not None:
            out.append(SweepPoint(weights=w, cds=cds))
    return tuple(out)


def greedy_max_min(
    points: Sequence[SweepPoint], k: int, *, distance: Callable[[str, str], float] = codon_distance
) -> tuple[SweepPoint, ...]:
    """Pick `k` points maximising the minimum pairwise distance, greedily.

    Seeded with the farthest-apart PAIR rather than an arbitrary first point.
    Greedy max-min is only as good as where it starts, and starting from the
    extremes of the set is what keeps the two ends of the trade-off in a small
    gallery -- exactly the designs a user most needs to see side by side.
    """
    if k < 1:
        raise ValueError(f"a gallery of {k} is not a gallery")
    if len(points) <= k:
        return tuple(points)

    best_pair = max(
        (
            (distance(a.cds, b.cds), i, j)
            for i, a in enumerate(points)
            for j, b in enumerate(points)
            if i < j
        ),
        default=(0.0, 0, 1),
    )
    chosen = [points[best_pair[1]], points[best_pair[2]]]
    remaining = [p for idx, p in enumerate(points) if idx not in (best_pair[1], best_pair[2])]

    while len(chosen) < k and remaining:
        far = max(remaining, key=lambda p: min(distance(p.cds, c.cds) for c in chosen))
        chosen.append(far)
        remaining.remove(far)
    return tuple(chosen)


def build_gallery(
    objectives: Sequence[str],
    solve: Callable[[Mapping[str, float]], str | None],
    *,
    steps: int = 8,
    k: int = 5,
) -> Gallery:
    """Sweep densely, deduplicate, then spread. The three steps are separate
    because conflating any two of them is a documented way to get this wrong."""
    if not MIN_GALLERY <= k <= MAX_GALLERY:
        raise ValueError(f"gallery size must be {MIN_GALLERY}-{MAX_GALLERY}, got {k}")
    swept = sweep(objectives, solve, steps=steps)
    # Deduplicate BEFORE selecting: many weight vectors land on the same design,
    # and a greedy pick over duplicates spends slots proving they are identical.
    seen: dict[str, SweepPoint] = {}
    for point in swept:
        seen.setdefault(point.cds, point)
    unique = tuple(seen.values())
    picks = greedy_max_min(unique, k)
    return Gallery(
        picks=picks,
        swept=len(swept),
        distinct=len(unique),
        min_pairwise_distance=pairwise_minimum([p.cds for p in picks]),
    )
