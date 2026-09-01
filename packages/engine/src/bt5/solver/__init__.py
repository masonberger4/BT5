"""M1 solver lane: the mutation space, the three tiers, and the wiring that runs
the rule catalog against them.

`catalog` is the entry point for real use -- it turns the registered rules into
the motifs Tier A forbids, the breaches Tier B repairs, the per-rule policies each
rule declares and the bounds the independent validator is given, all from one
selection so they cannot disagree. `pipeline.optimize` stays underneath it,
taking those as plain values, so a caller with three hand-written constraints
never has to load the catalog at all.
"""

from bt5.solver.catalog import (
    Findings,
    OracleBounds,
    RuleSet,
    build_rule_set,
    default_services,
    optimize_with,
)
from bt5.solver.lattice import (
    Automaton,
    achievable_gc_range,
    cai_lattice_scorer,
    combine_scorers,
    gc_steering_scorer,
    optimal_back_translate,
    solve_with_gc_steering,
)
from bt5.solver.pipeline import OptimizeResult, optimize
from bt5.solver.reference import (
    CodonScorer,
    back_translate,
    cai_scorer,
    expand_forbidden,
    longest_repeat,
    repeat_breaking_scorer,
)
from bt5.solver.repair import (
    Assembler,
    BreachFinder,
    RepairNotConverged,
    RepairOutcome,
    RulePolicy,
    codon_span,
    localize,
    repair,
)

__all__ = [
    "Assembler",
    "Automaton",
    "BreachFinder",
    "CodonScorer",
    "Findings",
    "OptimizeResult",
    "OracleBounds",
    "RepairNotConverged",
    "RepairOutcome",
    "RulePolicy",
    "RuleSet",
    "achievable_gc_range",
    "back_translate",
    "build_rule_set",
    "cai_lattice_scorer",
    "cai_scorer",
    "codon_span",
    "combine_scorers",
    "default_services",
    "expand_forbidden",
    "gc_steering_scorer",
    "localize",
    "longest_repeat",
    "optimal_back_translate",
    "optimize",
    "optimize_with",
    "repair",
    "repeat_breaking_scorer",
    "solve_with_gc_steering",
]
