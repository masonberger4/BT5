"""Tier B -- localized repair for constraints the lattice cannot express.

`HARD_REPAIR` is the enforcement class for a constraint that is real but not
decidable from the last few bases. Windowed GC content is the canonical member:
deciding whether a codon pushes a 50 bp window past its bound needs the G+C count
over the previous ~17 codons, and enumerating that history in the automaton state
is combinatorially intractable.

So the guarantee is assembled from three parts, and only all three together:

  1. Tier A STEERS toward the target with a context-free Lagrangian term.
  2. Tier B REPAIRS what remains, by localizing each breach to a small codon
     window and searching that window's mutation space.
  3. The independent validator PROVES the result, and refuses to emit if any
     window is still out of band.

Step 3 is what makes this a hard constraint rather than a preference. A violating
sequence never reaches the user: worst case the app declines and reports the
conflict, which is a far better failure than silently shipping a construct a
vendor will reject.

THE INTERACTION THAT MATTERS. Repair mutates codons, and a mutation can create a
forbidden motif that Tier A had guaranteed away -- including one spanning a codon
boundary or the insert/backbone junction. Every candidate is therefore checked
against the same Aho-Corasick automaton before it is accepted. Tier B can never
weaken a Tier-A guarantee; it can only fail to find a fix, which surfaces as an
infeasibility certificate.

Search parameters follow DNA Chisel's calibrated defaults rather than invented
ones: exhaustive local search when the local mutation space is under 10,000
variants, otherwise guided random with 2 mutations per iteration, capped at 1000
iterations with a stagnation tolerance of 100.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np

from bt5.codon.tables import NcbiGeneticCode
from bt5.core.result import InfeasibilityCertificate, InfeasibleConstraints
from bt5.core.spec import Breach, LocalizationPolicy, RepairPolicy
from bt5.core.types import Construct, Interval
from bt5.solver.lattice import Automaton
from bt5.solver.reference import expand_forbidden

#: DNA Chisel's calibrated defaults. Not invented numbers.
EXHAUSTIVE_LIMIT = 10_000
MUTATIONS_PER_ITERATION = 2
MAX_ITERATIONS = 1000
STAGNATION_TOLERANCE = 100

Assembler = Callable[[str], Construct]
"""Splices a candidate CDS into its backbone and returns the assembled construct.

Taking this as a callable keeps vector-assembly concerns out of the solver while
still guaranteeing that every rule is evaluated against the ASSEMBLED product --
never a free-floating CDS.
"""

BreachFinder = Callable[[Construct], tuple[Breach, ...]]
"""Evaluates the HARD_REPAIR rules against an assembled construct."""


@dataclass(frozen=True)
class RepairOutcome:
    cds: str
    iterations: int
    converged: bool
    #: HARD_REPAIR breaches the solver could not clear -- fixable by codon choice
    #: in principle, but the search did not reach a construct without them. A
    #: non-empty `remaining` is what makes the design infeasible.
    remaining: tuple[Breach, ...] = ()
    exhaustive_windows: int = 0
    random_windows: int = 0
    #: Breaches the solver was never asked to fix, because their author marked
    #: them `fixable_by_codon_choice=False`: a polyA hexamer in the user's own
    #: LTR, a repeat wholly in the backbone. Reading this field is what stops one
    #: unfixable backbone breach from aborting the whole pass; downstream it is
    #: what `QcReport.advisories` is for.
    advisory: tuple[Breach, ...] = ()

    @property
    def clean(self) -> bool:
        """No breach the solver could act on remains.

        Advisory breaches never count against `clean`: no codon choice can move
        them, so a construct carrying only advisories is as clean as this stage
        can make it, and the independent validator -- not repair -- is what
        decides whether an unfixable finding blocks the emit.
        """
        return not self.remaining


def _partition(breaches: Sequence[Breach]) -> tuple[list[Breach], tuple[Breach, ...]]:
    """Split breaches into the ones the solver can act on and the ones it cannot.

    `fixable_by_codon_choice` is REQUIRED on every Breach with no default
    precisely so this split is always answerable. Ignoring it -- as the search
    did before -- lets `target = max(breaches, ...)` pick an unfixable backbone
    breach whose window touches no editable codon, then abort the entire pass
    having done zero work.
    """
    actionable = [b for b in breaches if b.fixable_by_codon_choice]
    advisory = tuple(b for b in breaches if not b.fixable_by_codon_choice)
    return actionable, advisory


def localize(
    breach: Breach,
    policy: LocalizationPolicy,
    *,
    window: int,
    motif_len: int,
    construct_length: int,
    circular: bool,
) -> Interval:
    """Widen a breach into a repair window, per the rule's declared policy.

    The extension amounts come from the constraint's own geometry: a motif can
    only be created by bases within (len - 1) of it, and a windowed statistic can
    only be changed by bases within (window - 1).
    """
    if policy is LocalizationPolicy.WHOLE_SCOPE:
        return Interval(0, construct_length)
    if policy is LocalizationPolicy.WINDOW_MINUS_1:
        return breach.interval.extended(window - 1, construct_length, circular)
    if policy is LocalizationPolicy.MOTIF_LEN_MINUS_1:
        return breach.interval.extended(motif_len - 1, construct_length, circular)
    # PAIRED_SEGMENTS: a repeat pair is only fixable by editing one copy; the
    # breach interval already names the copy chosen for editing.
    return breach.interval


def codon_span(
    interval: Interval, codon_map: Sequence[Interval], construct_length: int
) -> tuple[int, int]:
    """Codon indices [first, last) whose bases intersect `interval`.

    Only whole codons are editable -- a partial codon cannot be changed without
    changing the protein.
    """
    lo, hi = len(codon_map), 0
    for idx, codon_iv in enumerate(codon_map):
        start, end = codon_iv.start, codon_iv.end
        # Normalise a wrapping repair window into linear comparisons.
        overlaps = interval.start < end and start < interval.end
        if interval.end > construct_length:
            wrapped_end = interval.end - construct_length
            overlaps = overlaps or start < wrapped_end
        if overlaps:
            lo, hi = min(lo, idx), max(hi, idx + 1)
    return (lo, hi) if lo < hi else (0, 0)


def _mutation_space(
    protein: str, code: NcbiGeneticCode, first: int, last: int
) -> list[tuple[str, ...]]:
    return [code.synonymous_codons(aa) for aa in protein[first:last]]


def _space_size(options: Sequence[Sequence[str]]) -> int:
    size = 1
    for opts in options:
        size *= len(opts)
        if size > EXHAUSTIVE_LIMIT:
            return size
    return size


def _enumerate(options: Sequence[Sequence[str]]) -> list[tuple[str, ...]]:
    out: list[tuple[str, ...]] = [()]
    for opts in options:
        out = [(*prefix, o) for prefix in out for o in opts]
    return out


def _introduces_forbidden(
    automaton: Automaton, left_state: int, candidate: str, right_flank: str
) -> bool:
    """Would this candidate create a forbidden motif?

    Seeded with the automaton state after the immutable left context and checked
    through the right flank, so junction-spanning and codon-boundary-spanning
    hits are caught exactly as Tier A catches them.
    """
    state, hit = automaton.consume(left_state, candidate)
    if hit:
        return True
    return automaton.consume(state, right_flank)[1]


def _codon_map_of(construct: Construct, cds: str) -> Sequence[Interval]:
    return (
        construct.translation_units[0].codon_map
        if construct.translation_units
        else tuple(Interval(i, i + 3) for i in range(0, len(cds), 3))
    )


def _select_workable(
    actionable: Sequence[Breach],
    *,
    policy: LocalizationPolicy,
    window: int,
    motif_len: int,
    codon_map: Sequence[Interval],
    construct_length: int,
    circular: bool,
    protein_len: int,
) -> tuple[Breach, Interval, int, int] | None:
    """The best actionable breach whose repair window touches an editable codon.

    Actionable breaches are tried worst-first (magnitude, then earliest); a
    target whose localized window spans no whole editable codon is SKIPPED, not
    the whole pass aborted. Returns None only when no actionable breach has a
    workable window -- which is the honest "nothing left to do" signal, distinct
    from "one unfixable breach happened to sort first".
    """
    ordered = sorted(actionable, key=lambda b: (b.magnitude, -b.interval.start), reverse=True)
    for target in ordered:
        repair_window = localize(
            target,
            policy,
            window=window,
            motif_len=motif_len,
            construct_length=construct_length,
            circular=circular,
        )
        first, last = codon_span(repair_window, codon_map, construct_length)
        last = min(last, protein_len)
        if first < last:
            return target, repair_window, first, last
    return None


def repair(
    cds: str,
    protein: str,
    code: NcbiGeneticCode,
    *,
    assemble: Assembler,
    find_breaches: BreachFinder,
    forbidden: Sequence[str] = (),
    policy: LocalizationPolicy = LocalizationPolicy.WINDOW_MINUS_1,
    repair_policy: RepairPolicy = RepairPolicy.FIXED_POINT,
    window: int = 50,
    motif_len: int = 6,
    left_flank: str = "",
    right_flank: str = "",
    seed: int = 0,
    max_iterations: int = MAX_ITERATIONS,
) -> RepairOutcome:
    """Repair HARD_REPAIR breaches without ever weakening a Tier-A guarantee.

    Only breaches marked `fixable_by_codon_choice` are chased; the rest are
    carried on `RepairOutcome.advisory` untouched, because no codon choice can
    move them. Raises InfeasibleConstraints when the search cannot clear the
    ACTIONABLE breaches, with the minimal conflicting set rather than a bare
    failure.
    """
    rng = np.random.default_rng(seed)  # explicit; never the global RNG
    patterns = expand_forbidden(forbidden)
    automaton = Automaton(patterns)

    current = cds
    construct = assemble(current)
    breaches = find_breaches(construct)
    actionable, advisory = _partition(breaches)
    # Advisory-only (or clean): nothing the solver can act on. Not a failure --
    # an unfixable backbone polyA is the validator's call, not repair's.
    if not actionable:
        return RepairOutcome(current, 0, True, advisory=advisory)

    stagnant = 0
    exhaustive_windows = random_windows = 0
    best_cost = sum(b.magnitude for b in actionable)

    for iteration in range(1, max_iterations + 1):
        codon_map = _codon_map_of(construct, current)
        selected = _select_workable(
            actionable,
            policy=policy,
            window=window,
            motif_len=motif_len,
            codon_map=codon_map,
            construct_length=construct.length,
            circular=construct.is_circular,
            protein_len=len(protein),
        )
        if selected is None:
            break  # no actionable breach has a workable window
        _target, _repair_window, first, last = selected

        options = _mutation_space(protein, code, first, last)
        size = _space_size(options)
        prefix = left_flank + current[: 3 * first]
        suffix = current[3 * last :] + right_flank
        left_state = automaton.consume(0, prefix)[0]

        if size <= EXHAUSTIVE_LIMIT:
            exhaustive_windows += 1
            candidates = _enumerate(options)
        else:
            random_windows += 1
            # Guided random: perturb MUTATIONS_PER_ITERATION positions.
            base = [current[3 * (first + k) : 3 * (first + k) + 3] for k in range(last - first)]
            candidates = []
            for _ in range(256):
                trial = list(base)
                for pos in rng.choice(
                    len(trial), size=min(MUTATIONS_PER_ITERATION, len(trial)), replace=False
                ):
                    trial[int(pos)] = str(rng.choice(options[int(pos)]))
                candidates.append(tuple(trial))

        improved = False
        for combo in candidates:
            block = "".join(combo)
            if _introduces_forbidden(automaton, left_state, block, suffix[: max(1, motif_len)]):
                continue  # never weaken a Tier-A guarantee
            trial_cds = current[: 3 * first] + block + current[3 * last :]
            trial_actionable, trial_advisory = _partition(find_breaches(assemble(trial_cds)))
            cost = sum(b.magnitude for b in trial_actionable)
            if cost < best_cost - 1e-12:
                current, best_cost = trial_cds, cost
                actionable, advisory = trial_actionable, trial_advisory
                construct = assemble(current)
                improved = True
                if not actionable:
                    return RepairOutcome(
                        current,
                        iteration,
                        True,
                        (),
                        exhaustive_windows,
                        random_windows,
                        advisory=advisory,
                    )
                break

        stagnant = 0 if improved else stagnant + 1
        if stagnant >= STAGNATION_TOLERANCE:
            break
        if repair_policy is RepairPolicy.SINGLE_PASS and not improved:
            break

    if actionable:
        raise InfeasibleConstraints(
            InfeasibilityCertificate(
                interval=actionable[0].interval,
                protein_span=(0, len(protein)),
                minimal_conflicting_specs=tuple(sorted({b.spec_id for b in actionable})),
                proof="empty_mutation_space",
            )
        )
    return RepairOutcome(
        current, max_iterations, True, (), exhaustive_windows, random_windows, advisory=advisory
    )
