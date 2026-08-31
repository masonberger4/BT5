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

from collections.abc import Callable, Mapping, Sequence
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


def _aggregate(actionable: Sequence[Breach]) -> dict[str, tuple[int, float]]:
    """Per-rule (count, magnitude-sum). The cost the search actually compares.

    Breach COUNT is the cross-rule currency -- one `Breach` is "one localized
    problem" by contract, so a GC deviation of 0.05 and a repeat of magnitude 4.0
    are each worth exactly one. Summing magnitudes across rules would put
    kcal/mol, CAI and integer motif counts on one axis and let the largest-unit
    rule silently dominate. Magnitude therefore stays INSIDE its rule, as the
    within-rule gradient the tie-break uses.
    """
    agg: dict[str, tuple[int, float]] = {}
    for b in actionable:
        count, magsum = agg.get(b.spec_id, (0, 0.0))
        agg[b.spec_id] = (count + 1, magsum + b.magnitude)
    return agg


def _accepts(
    current: Mapping[str, tuple[int, float]],
    trial: Mapping[str, tuple[int, float]],
    target_spec: str,
) -> bool:
    """Is `trial` a strict improvement the search should take?

    The rules, in order, are what make the search monotone and terminating while
    keeping magnitude inside its rule:

    - Never let ANY rule's breach count rise. Count is the currency, so an
      accepted move cannot trade one rule's breach for another's. This is also
      the invariant the termination proof rests on: the per-rule count vector is
      non-increasing, and integer, so only finitely many count-reducing moves
      exist.
    - If total count falls, take it -- fewer localized problems is unambiguously
      better across rules.
    - If total count is unchanged (and, from the first rule, no single count
      rose, so every count is identical), require within-rule progress: the
      TARGET rule's magnitude-sum strictly falls and no other rule's rises. This
      is the GC-gradient step, and it never compares magnitudes between rules.
    """
    for spec, (tc, _tm) in trial.items():
        if tc > current.get(spec, (0, 0.0))[0]:
            return False

    cur_total = sum(c for c, _ in current.values())
    trial_total = sum(c for c, _ in trial.values())
    if trial_total < cur_total:
        return True
    if trial_total > cur_total:
        return False

    # Counts are identical rule-by-rule; this is a pure magnitude move.
    cur_mag = current.get(target_spec, (0, 0.0))[1]
    trial_mag = trial.get(target_spec, (0, 0.0))[1]
    if trial_mag >= cur_mag - 1e-12:
        return False
    return all(
        tm <= current.get(spec, (0, 0.0))[1] + 1e-12
        for spec, (_c, tm) in trial.items()
        if spec != target_spec
    )


def _window_for(
    breach: Breach,
    *,
    policy: LocalizationPolicy,
    window: int,
    motif_len: int,
    codon_map: Sequence[Interval],
    construct_length: int,
    circular: bool,
    protein_len: int,
) -> tuple[Interval, int, int] | None:
    """A breach's repair window and its codon span, or None if it spans no whole
    editable codon (so the search must skip it, never abort on it)."""
    repair_window = localize(
        breach,
        policy,
        window=window,
        motif_len=motif_len,
        construct_length=construct_length,
        circular=circular,
    )
    first, last = codon_span(repair_window, codon_map, construct_length)
    last = min(last, protein_len)
    if first >= last:
        return None
    return repair_window, first, last


def _select_target(
    actionable: Sequence[Breach],
    *,
    priority: Mapping[str, int],
    turns: dict[str, int],
    policy: LocalizationPolicy,
    window: int,
    motif_len: int,
    codon_map: Sequence[Interval],
    construct_length: int,
    circular: bool,
    protein_len: int,
) -> tuple[Breach, Interval, int, int] | None:
    """Pick the next breach to work, round-robin over rules, priority first.

    Worst-first over raw magnitude STARVES a rule whose findings are small: a GC
    deviation of 0.05 never outranks a repeat of magnitude 4.0, so with both
    present the GC window is never targeted and the design is emitted with it
    still out of band. Selection is therefore two-level: among the highest
    `priority` rules that currently have a workable breach, take the one served
    least often (round-robin, ties by spec_id); within that rule, take its worst
    breach (magnitude, then earliest). Default priority 0 for every rule is pure
    round-robin, so a single-rule pass is unchanged.

    Returns None only when NO actionable breach has a workable window -- the
    honest "nothing left to do" signal.
    """
    workable: dict[str, list[tuple[Breach, Interval, int, int]]] = {}
    for breach in actionable:
        found = _window_for(
            breach,
            policy=policy,
            window=window,
            motif_len=motif_len,
            codon_map=codon_map,
            construct_length=construct_length,
            circular=circular,
            protein_len=protein_len,
        )
        if found is None:
            continue
        repair_window, first, last = found
        workable.setdefault(breach.spec_id, []).append((breach, repair_window, first, last))
    if not workable:
        return None

    top = max(priority.get(spec, 0) for spec in workable)
    eligible = [spec for spec in workable if priority.get(spec, 0) == top]
    chosen = min(eligible, key=lambda spec: (turns.get(spec, 0), spec))
    turns[chosen] = turns.get(chosen, 0) + 1
    return max(workable[chosen], key=lambda t: (t[0].magnitude, -t[0].interval.start))


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
    priority: Mapping[str, int] | None = None,
    max_iterations: int = MAX_ITERATIONS,
) -> RepairOutcome:
    """Repair HARD_REPAIR breaches without ever weakening a Tier-A guarantee.

    Only breaches marked `fixable_by_codon_choice` are chased; the rest are
    carried on `RepairOutcome.advisory` untouched, because no codon choice can
    move them. Breach count is the cross-rule cost; a rule's breaches are worked
    round-robin so a small-magnitude rule is never starved by a large-magnitude
    one; `priority` (default 0 for every rule) lets a caller work some rules
    before others. Raises InfeasibleConstraints when the search cannot clear the
    ACTIONABLE breaches, with the minimal conflicting set rather than a bare
    failure.
    """
    rng = np.random.default_rng(seed)  # explicit; never the global RNG
    patterns = expand_forbidden(forbidden)
    automaton = Automaton(patterns)
    priority = priority or {}

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
    cur_agg = _aggregate(actionable)
    turns: dict[str, int] = {}

    for iteration in range(1, max_iterations + 1):
        codon_map = _codon_map_of(construct, current)
        selected = _select_target(
            actionable,
            priority=priority,
            turns=turns,
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
        target, _repair_window, first, last = selected
        target_spec = target.spec_id

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
            trial_agg = _aggregate(trial_actionable)
            if _accepts(cur_agg, trial_agg, target_spec):
                current, cur_agg = trial_cds, trial_agg
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
