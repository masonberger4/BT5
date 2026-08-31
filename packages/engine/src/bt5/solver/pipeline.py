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
mistake. That is also why Tier B is asked NOT to raise here: the validator, not
an exception from the search, is what decides whether a construct is emitted.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from bt5.codon.tables import NcbiGeneticCode
from bt5.core.types import Construct
from bt5.solver.lattice import optimal_back_translate, solve_with_gc_steering
from bt5.solver.reference import CodonScorer
from bt5.solver.repair import Assembler, BreachFinder, RepairOutcome, RulePolicy, repair
from bt5.verify import verify_construct


@dataclass(frozen=True)
class OptimizeResult:
    cds: str
    construct: Construct
    repair_outcome: RepairOutcome
    verified: bool = True
    #: Whether Tier B actually ran. False means the caller passed no
    #: `find_breaches`, so `repair_outcome` is a "not run" placeholder rather than
    #: a clean repair -- otherwise the two are indistinguishable, and a report
    #: could imply the HARD_REPAIR rules were checked and passed when they were
    #: never evaluated.
    tier_b_ran: bool = True


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
    priority: Mapping[str, int] | None = None,
    policies: Mapping[str, RulePolicy] | None = None,
    max_candidates: int = 256,
    original_backbone: Construct | None = None,
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
    # `gc_bounds` is the only band signal available here, so when it is present
    # steer toward it: the docstring above has always promised Tier-A steering,
    # and until now `optimize()` did not do it, so a GC-extreme protein could
    # fail on a design that steering would have reached. Steering early-returns
    # at multiplier 0.0 when the unsteered solve is already in band, so the
    # common case costs one extra DP solve and nothing more.
    if gc_bounds is not None:
        cds, _multiplier = solve_with_gc_steering(
            protein,
            code,
            gc_bounds=gc_bounds,
            base_score=score,
            forbidden=forbidden,
            left_flank=left_flank,
            right_flank=right_flank,
        )
    else:
        cds = optimal_back_translate(
            protein,
            code,
            forbidden=forbidden,
            score=score,
            left_flank=left_flank,
            right_flank=right_flank,
        )

    # --- Tier B: localized repair of what the lattice cannot express -------
    tier_b_ran = find_breaches is not None
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
            priority=priority,
            policies=policies,
            max_candidates=max_candidates,
            # The validator below is the refusal, not this. Repair returns its
            # best attempt and the honest certificate; a construct still out of
            # band is caught by verify_construct, on an independent code path.
            raise_on_infeasible=False,
        )
        cds = outcome.cds
    else:
        # Distinguishable from "ran and found nothing": a fabricated clean
        # outcome would let a report claim the HARD_REPAIR rules passed when they
        # were never run.
        outcome = RepairOutcome(cds, 0, True, stop_reason="not_run")

    construct = assemble(cds)

    # --- The validator: independent, and it refuses ------------------------
    if _verify:
        verify_construct(
            construct,
            protein=protein,
            table_id=tid,
            forbidden=forbidden,
            gc_bounds=gc_bounds,
            # I9 -- byte-identity of the backbone complement, the invariant
            # verify.py calls "highest value". Unarmed on every design path until
            # a caller supplies the parsed input backbone (E1 passes
            # `assembly.reference`).
            original_backbone=original_backbone,
        )

    return OptimizeResult(
        cds=cds, construct=construct, repair_outcome=outcome, tier_b_ran=tier_b_ran
    )
