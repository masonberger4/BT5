"""The three-tier optimize() entry point.

This is where the HARD_REPAIR guarantee is actually assembled. Each tier alone is
insufficient and the composition is the point:

    Tier A   exact DP          forbidden motifs unreachable BY CONSTRUCTION
             + GC steering     composition pulled toward the band
    Tier B   localized repair   windowed statistics brought into the band
    Validator                   PROVES it, and REFUSES TO EMIT otherwise

The refusal is not a formality. Without it, "hard" would mean "we tried", which
is exactly the failure mode the enforcement enum exists to prevent -- and the
independent validator lives in bt5.verify, on a different code path from the
scorer that guided the search, so it cannot rubber-stamp the optimizer's own
mistake.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from bt5.codon.tables import NcbiGeneticCode
from bt5.core.types import Construct
from bt5.solver.lattice import optimal_back_translate
from bt5.solver.reference import CodonScorer
from bt5.solver.repair import Assembler, BreachFinder, RepairOutcome, repair
from bt5.verify import verify_construct


@dataclass(frozen=True)
class OptimizeResult:
    cds: str
    construct: Construct
    repair_outcome: RepairOutcome
    verified: bool = True


def optimize(
    protein: str,
    code: NcbiGeneticCode,
    *,
    assemble: Assembler,
    find_breaches: BreachFinder | None = None,
    forbidden: Sequence[str] = (),
    score: CodonScorer | None = None,
    gc_bounds: tuple[float, float] | None = None,
    gc_window: int = 50,
    left_flank: str = "",
    right_flank: str = "",
    seed: int = 0,
    table_id: int | None = None,
    _verify: bool = True,
) -> OptimizeResult:
    """Design a CDS, then prove it before returning it.

    `_verify` exists only so the verification path itself can be tested. It
    defaults to True and a CI grep asserts that default, because an optimizer
    that can be asked politely not to check itself is not an optimizer with a
    hard-constraint guarantee.
    """
    tid = table_id if table_id is not None else code.table_id

    # --- Tier A: exact, motif-free by construction -------------------------
    cds = optimal_back_translate(
        protein,
        code,
        forbidden=forbidden,
        score=score,
        left_flank=left_flank,
        right_flank=right_flank,
    )

    # --- Tier B: localized repair of what the lattice cannot express -------
    outcome = RepairOutcome(cds, 0, True)
    if find_breaches is not None:
        outcome = repair(
            cds,
            protein,
            code,
            assemble=assemble,
            find_breaches=find_breaches,
            forbidden=forbidden,
            window=gc_window,
            left_flank=left_flank,
            right_flank=right_flank,
            seed=seed,
        )
        cds = outcome.cds

    construct = assemble(cds)

    # --- The validator: independent, and it refuses ------------------------
    if _verify:
        verify_construct(
            construct,
            protein=protein,
            table_id=tid,
            forbidden=forbidden,
            gc_bounds=gc_bounds,
            gc_window=gc_window,
        )

    return OptimizeResult(cds=cds, construct=construct, repair_outcome=outcome)
