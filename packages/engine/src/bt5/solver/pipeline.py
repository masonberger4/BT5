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

`bt5.solver.catalog` is what supplies the arguments below from the shipped rule
catalog. This module stays primitive on purpose: it takes motifs, a breach
finder, policies and validator bounds as VALUES, so it never imports the rules,
the vector lane or the folding engine, and a caller with three hand-written
constraints can still use it.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, replace

from bt5.codon.tables import NcbiGeneticCode
from bt5.core.result import InfeasibilityCertificate, InfeasibleConstraints
from bt5.core.spec import Breach, LocalizationPolicy, RepairPolicy
from bt5.core.types import Construct, Interval
from bt5.solver.lattice import optimal_back_translate
from bt5.solver.reference import CodonScorer
from bt5.solver.repair import (
    Assembler,
    BreachFinder,
    Cost,
    RepairOutcome,
    no_rules,
    repair,
)
from bt5.verify import find_motifs, verify_construct


@dataclass(frozen=True)
class OptimizeResult:
    cds: str
    construct: Construct
    repair_outcome: RepairOutcome
    verified: bool = True
    #: Hard findings NO codon can fix -- an over-length fragment, a homopolymer
    #: the user's own backbone carries, a palindrome inside an ITR. Reported
    #: rather than routed to the solver, which would exhaust the mutation space
    #: chasing a fix that does not exist and report a fine design infeasible.
    advisories: tuple[Breach, ...] = ()


def optimize(
    protein: str,
    code: NcbiGeneticCode,
    *,
    assemble: Assembler,
    find_breaches: BreachFinder,
    forbidden: Sequence[str] = (),
    score: CodonScorer | None = None,
    gc_bounds: tuple[float, float] | None = None,
    gc_window: int = 50,
    localization: LocalizationPolicy | Callable[[str], LocalizationPolicy] = (
        LocalizationPolicy.WINDOW_MINUS_1
    ),
    repair_policy: RepairPolicy = RepairPolicy.FIXED_POINT,
    cost: Callable[[Sequence[Breach]], Cost] | None = None,
    advise: Callable[[Construct], tuple[Breach, ...]] | None = None,
    max_homopolymer: int | None = None,
    original_backbone: Construct | None = None,
    left_flank: str = "",
    right_flank: str = "",
    seed: int = 0,
    table_id: int | None = None,
    _verify: bool = True,
) -> OptimizeResult:
    """Design a CDS, then prove it before returning it.

    `find_breaches` is REQUIRED and has no default. It used to default to None,
    and the None branch fabricated `RepairOutcome(cds, 0, converged=True)` --
    zero iterations, converged, nothing remaining -- which is byte-identical to
    "Tier B ran and found nothing to fix". Forgetting the argument therefore
    produced a construct that looked fully repaired and had never been checked.
    To skip Tier B deliberately, pass `catalog.no_rules`, which says so and
    reports `RepairOutcome.ran is False`.

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

    # --- The immutable-region screen ---------------------------------------
    #
    # Tier A already refuses a motif in the left flank it was seeded with
    # (`lattice.py:161-170`), but the flank is the junction context, not the
    # vector. A motif sitting deep in the backbone survives Tier A, survives
    # repair -- no codon can reach it -- and surfaces at the very end as a bare
    # `VerificationError("I6", ...)` after a full solve, naming a position and
    # nothing else. That was tolerable while `forbidden` was three sites a caller
    # typed by hand; with a catalog-derived set it is the likely failure, so it
    # gets the certificate the frozen `proof` enum has been carrying an unused
    # value for since PR #1.
    if forbidden:
        screened = assemble(cds)
        stuck = [
            (pattern, pos)
            for pattern, pos in find_motifs(screened, list(forbidden))
            if not screened.overlaps_editable(Interval(pos, pos + len(pattern)))
        ]
        if stuck:
            pattern, pos = stuck[0]
            raise InfeasibleConstraints(
                InfeasibilityCertificate(
                    interval=Interval(pos, pos + len(pattern)),
                    protein_span=(0, len(protein)),
                    minimal_conflicting_specs=tuple(sorted({p for p, _ in stuck})),
                    proof="immutable_region",
                )
            )

    # --- Tier B: localized repair of what the lattice cannot express -------
    outcome = repair(
        cds,
        protein,
        code,
        assemble=assemble,
        find_breaches=find_breaches,
        forbidden=forbidden,
        policy=localization,
        repair_policy=repair_policy,
        cost=cost,
        window=gc_window,
        left_flank=left_flank,
        right_flank=right_flank,
        seed=seed,
    )
    if find_breaches is no_rules:
        outcome = replace(outcome, ran=False)
    cds = outcome.cds

    construct = assemble(cds)

    # --- The validator: independent, and it refuses ------------------------
    #
    # `original_backbone` is now forwarded, which is how I9 stopped being dead
    # code on the only path in src/ that calls the oracle at all. I9 -- every
    # backbone base byte-identical to the input -- is the highest-value invariant
    # in the file, and it had never run from optimize(). It arrives as a VALUE:
    # verify.py may not import the vector lane that produces the reference, and
    # that independence is the point.
    #
    # `max_homopolymer` is a parameter a caller MAY pass to arm I8's homopolymer
    # scan, but the catalog path (`catalog.optimize_with`) deliberately leaves it
    # unset: E1's run limits already reach the oracle as `forbidden` motifs, and
    # I6 proves those absent on both strands WHILE honouring `Construct.exempt`,
    # which the I8 homopolymer scan does not. See `catalog.OracleBounds`.
    #
    # `max_repeat` is not passed either, for the same reason one level worse:
    # I8's repeat scan walks every k-mer with no reference to `Construct.exempt`,
    # so arming it would refuse every lentiviral vector on its own identical LTRs.
    if _verify:
        verify_construct(
            construct,
            protein=protein,
            table_id=tid,
            forbidden=forbidden,
            gc_bounds=gc_bounds,
            max_homopolymer=max_homopolymer,
            original_backbone=original_backbone,
        )

    return OptimizeResult(
        cds=cds,
        construct=construct,
        repair_outcome=outcome,
        advisories=advise(construct) if advise is not None else (),
    )
