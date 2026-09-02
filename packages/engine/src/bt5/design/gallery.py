"""Solving the same protein several ways, and keeping the ones that differ.

`bt5.score.gallery` already knows how to CHOOSE a gallery -- sweep densely,
deduplicate, then spread by greedy max-min on codon distance, because Das &
Dennis (1997) rules out picking evenly spaced weights and calling the results
evenly spread. What it does not know is how to solve anything. This module is
the join: it turns a weight vector into a real, repaired, VERIFIED design and
hands `build_gallery` a `solve` it can sweep.

Every candidate is a full `optimize()` -- Tier A, Tier B, and the independent
validator's refusal. A gallery of five designs where four were never verified
would be worse than no gallery, because the user picks one of them to order.
That is also why an infeasible weight vector is DROPPED rather than repaired
into something else: `sweep`'s docstring calls an infeasible corner of the
simplex "a real answer about the constraint set, not a gap to fill in", and the
count of what survived travels on `Gallery.swept` / `Gallery.distinct`.

The sweep is deliberately coarse, and the number is measured rather than
guessed. `simplex_weights` over four axes at the shipped `steps=8` is 165
lattice points, so 165 full solves; at 500 aa on the reference backbone that is
far outside PLAN's 10 s G7 bar. But the reason `DEFAULT_SWEEP_STEPS` is 1 is not
only cost -- it is that density bought NOTHING. Measured on that fixture, 3, 4, 6
and 20 weight vectors each yielded exactly **3 distinct designs**, at 5.2 s,
7.5 s, 10.4 s and 33.6 s. The extra 17 solves rediscovered designs the first
three had already found.

At `steps=1` the lattice is the one-hot corners: each axis pushed alone. Those
are the EXTREMES of the trade-off, which is what `greedy_max_min` seeds itself
from anyway -- "starting from the extremes of the set is what keeps the two ends
of the trade-off in a small gallery". Mixtures interpolate between designs the
corners already reach.

Cost is dominated by Tier B, not by the DP. A solve whose CDS is already clean
costs ~0.09 s; one that needs repair costs ~2.5 s, because `repair()` evaluates
up to `max_candidates` candidates per iteration and each is a full catalog pass.
So the lever that matters is the NUMBER OF SOLVES, and `sweep` calling `solve`
once per lattice point is why the lattice is small.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field

from bt5.codon.tables import NcbiGeneticCode
from bt5.core.result import InfeasibleConstraints, VerificationError
from bt5.core.types import Construct
from bt5.score.gallery import Gallery, build_gallery
from bt5.score.steering import blended_scorer, live_axes
from bt5.solver.pipeline import OptimizeResult, optimize
from bt5.solver.repair import BreachFinder, RulePolicy

#: Lattice density per axis. 1 is the one-hot corners -- each axis pushed alone.
#: See the module docstring: on the reference fixture, raising this to 2 or 3
#: multiplied the solve count by 2x and 6.7x and produced the same number of
#: distinct designs. Raise it for a protein whose front is richer; it is a
#: parameter of `design()` for exactly that reason.
DEFAULT_SWEEP_STEPS = 1


@dataclass(frozen=True)
class SolveSpace:
    """Everything `optimize()` needs except the objective weights.

    Frozen and passed whole so a sweep cannot drift: every candidate in a
    gallery must be solved against the SAME forbidden set, the same GC band, the
    same policies and the same I9 reference, or the designs are not alternatives
    to each other and the codon distances between them mean nothing.
    """

    protein: str
    code: NcbiGeneticCode
    assemble: Callable[[str], Construct]
    forbidden: tuple[str, ...]
    seed: int
    table_id: int
    usage: Mapping[str, float]
    find_breaches: BreachFinder | None = None
    gc_bounds: tuple[float, float] | None = None
    gc_window: int = 50
    left_flank: str = ""
    right_flank: str = ""
    policies: Mapping[str, RulePolicy] = field(default_factory=dict)
    max_candidates: int = 256
    reference: Construct | None = None

    def solve(self, weights: Mapping[str, float] | None) -> OptimizeResult | None:
        """One verified design at these weights, or None if there is none.

        `None` for `weights` is the unsteered solve -- `score=None`, which is
        what the walking skeleton did and what the design falls back to when the
        sweep yields nothing.

        Both failures caught here are answers rather than errors.
        `InfeasibleConstraints` carries a proof that this corner of the simplex
        has no solution; `VerificationError` means the validator refused to emit
        what the search produced, and refusing is the guarantee working. Neither
        may be allowed to take down a design that has 19 other weight vectors to
        try -- and neither may be swallowed into a candidate, which is why this
        returns None instead of a partial result.
        """
        score = (
            blended_scorer(weights, usage=self.usage, gc_bounds=self.gc_bounds)
            if weights is not None
            else None
        )
        try:
            return optimize(
                self.protein,
                self.code,
                assemble=self.assemble,
                find_breaches=self.find_breaches,
                forbidden=self.forbidden,
                score=score,
                gc_bounds=self.gc_bounds,
                gc_window=self.gc_window,
                left_flank=self.left_flank,
                right_flank=self.right_flank,
                seed=self.seed,
                table_id=self.table_id,
                policies=dict(self.policies),
                max_candidates=self.max_candidates,
                original_backbone=self.reference,
            )
        except (InfeasibleConstraints, VerificationError):
            return None


def sweep_designs(
    space: SolveSpace,
    *,
    axes: Sequence[str] | None = None,
    steps: int = DEFAULT_SWEEP_STEPS,
    k: int = 5,
) -> tuple[Gallery, dict[str, OptimizeResult]]:
    """Sweep, select, and keep each pick's proven `OptimizeResult`.

    `axes` defaults to the axes that can actually move THIS design's answer.
    With no host codon-usage table on file, `codon_adaptation`'s term is
    identically zero and its corner of the simplex solves to the unsteered
    design -- a full `optimize()` spent rediscovering a design another vector
    already found. `live_axes` drops it; see its docstring.

    `build_gallery`'s `solve` returns only a CDS, but the repair outcome and the
    assembled construct behind each pick are what the report and the GenBank are
    built from -- and re-solving to recover them would be a second search that
    could return something else. So they are recorded on the way past, keyed by
    the CDS that `build_gallery` will hand back.
    """
    solved: dict[str, OptimizeResult] = {}

    def solve(weights: Mapping[str, float]) -> str | None:
        result = space.solve(weights)
        if result is None:
            return None
        solved.setdefault(result.cds, result)
        return result.cds

    gallery = build_gallery(axes or live_axes(space.usage), solve, steps=steps, k=k)
    return gallery, solved
