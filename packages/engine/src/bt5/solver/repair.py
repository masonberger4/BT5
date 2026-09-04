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

#: Floor on codons perturbed per random trial. Only a floor: the count scales
#: with the window, see `MUTATION_WINDOW_DIVISOR`.
MUTATIONS_PER_ITERATION = 2

#: The random branch perturbs `len(window) / this` codons, floored at
#: `MUTATIONS_PER_ITERATION`.
#:
#: A fixed 2 was calibrated for nothing, and it is unreachable for a windowed
#: COMPOSITION rule. `f5_at_window` is a two-sided GC band over 100 nt with
#: `WINDOW_MINUS_1` localisation, so its repair window is the breach plus 99 nt
#: of context either side -- 67-100 codons around a 34-codon breach. Moving that
#: window's GC fraction across the band needs a coordinated multi-codon shift;
#: two substitutions cannot deliver one however many times they are drawn, and
#: re-drawing was exactly the waste `PER_TARGET_TOLERANCE` now bounds.
#:
#: Targeting was never the problem: two uniform picks land BOTH inside a
#: 34-codon breach 11-25% of the time, so dozens of well-placed trials happen
#: every iteration. Measured on the 260 aa reference design, the fourth solve:
#:
#:     2 -> 38 iterations, stagnation, 5 breaches left, candidate discarded
#:     4 -> 22 iterations, clean          16 -> 15 iterations, clean
#:     8 -> 27 iterations, clean          32 ->  3 iterations, clean
#:
#: A quarter of the window rather than a tuned constant, because the quantity
#: that must move scales with the window: 25 codons at 100, 17 at 67, and the
#: floor at 8 or fewer. Small windows are enumerated outright and never reach
#: this branch, so this only ever governs windows too large to search exactly.
MUTATION_WINDOW_DIVISOR = 4

MAX_ITERATIONS = 1000
STAGNATION_TOLERANCE = 100

#: Consecutive failed attempts on ONE breach before it stops being re-selected
#: while the sequence stands still.
#:
#: A failed iteration leaves `current` untouched, so re-targeting the same breach
#: re-runs the same window over the same options with only a fresh RNG draw --
#: a memoryless lottery, re-rolled. `STAGNATION_TOLERANCE` alone bounds only the
#: GLOBAL run of failures, so N hopeless breaches cost N x that many `find_breaches`
#: calls before anything gives up: measured at 25,601 calls for one unclearable
#: breach (#111).
#:
#: Three, because the search is memoryless: k failures are k independent
#: `max_candidates`-draw samples, so 3 x 256 = 768 consecutive misses put a
#: ~95% upper bound of roughly 0.4% on the per-draw success probability. A
#: fourth round buys a vanishing amount of search for another 256 evaluations.
#: Any improvement anywhere resets this -- see `_abandon` handling in `repair`.
PER_TARGET_TOLERANCE = 3

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
    #: Why the pass stopped. `clean` (nothing actionable left), `not_run` (Tier B
    #: was never invoked -- set by the pipeline, not here), `exhausted_targets`
    #: (every actionable breach was attempted and none remained workable),
    #: `stagnation` / `iterations` (a resource cap was hit before the space was
    #: searched out -- the search GAVE UP, and `converged` is False).
    stop_reason: str = "clean"
    #: Present only on a non-clean outcome: the honest certificate of what was
    #: worked last and searched. `None` when `clean`.
    certificate: InfeasibilityCertificate | None = None

    @property
    def clean(self) -> bool:
        """No breach the solver could act on remains.

        Advisory breaches never count against `clean`: no codon choice can move
        them, so a construct carrying only advisories is as clean as this stage
        can make it, and the independent validator -- not repair -- is what
        decides whether an unfixable finding blocks the emit.
        """
        return not self.remaining


class RepairNotConverged(InfeasibleConstraints):
    """The search hit a resource cap (stagnation or the iteration limit) before
    it could search the mutation space out.

    A subclass of InfeasibleConstraints so every existing caller and test that
    catches the base type still catches this, but a distinct type on PURPOSE: the
    base class asserts a PROOF of infeasibility, and a search that merely gave up
    has not proven anything. Until `InfeasibilityCertificate.proof` gains a
    `search_not_converged` member (a MAJOR contract change, filed separately),
    the type name and this message are what keep the certificate from claiming
    the mutation space was empty when it was only unsearched.
    """


@dataclass(frozen=True)
class RulePolicy:
    """How Tier B localizes, budgets and orders ONE rule's breaches.

    Repair works many rules at once, and they do not share a localization
    heuristic (a windowed GC statistic widens by `window - 1`; a motif by
    `motif_len - 1`), a repair discipline (SINGLE_PASS retires a breach after one
    attempt; FIXED_POINT re-targets until the rule stops producing it, which is
    what splice removal requires), or a place in the queue (a hard rule should be
    cleared before a soft one gets a turn). A caller supplies one of these per
    `spec_id`; the scalar `policy`/`repair_policy`/`window`/`motif_len` arguments
    to `repair()` are the fallback for any rule without one.

    This is caller-supplied on PURPOSE: importing the registry to read the
    policies off the rule classes would pull the whole 15-rule catalog, Biopython
    and the vendor registry into the solver, make Tier A untestable without the
    catalog, and violate the lane boundary -- and the values live on rule
    INSTANCES (`e2`'s window is `self.window`), which the classes do not carry.
    """

    localization: LocalizationPolicy
    repair: RepairPolicy
    window: int
    motif_len: int
    priority: int = 0


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


def _breach_key(breach: Breach) -> tuple[str, int, int]:
    """A breach's identity across iterations: rule plus reported span."""
    return (breach.spec_id, breach.interval.start, breach.interval.end)


def _select_target(
    actionable: Sequence[Breach],
    *,
    resolve: Callable[[str], RulePolicy],
    turns: dict[str, int],
    retired: set[tuple[str, int, int]],
    codon_map: Sequence[Interval],
    construct_length: int,
    circular: bool,
    protein_len: int,
) -> tuple[Breach, Interval, int, int, RulePolicy] | None:
    """Pick the next breach to work, round-robin over rules, priority first.

    Worst-first over raw magnitude STARVES a rule whose findings are small: a GC
    deviation of 0.05 never outranks a repeat of magnitude 4.0, so with both
    present the GC window is never targeted and the design is emitted with it
    still out of band. Selection is therefore two-level: among the highest
    `priority` rules that currently have a workable, un-retired breach, take the
    one served least often (round-robin, ties by spec_id); within that rule, take
    its worst breach (magnitude, then earliest). Each rule's window comes from
    ITS own `RulePolicy`. Default priority 0 for every rule is pure round-robin,
    so a single-rule pass is unchanged.

    Returns None only when NO actionable breach has a workable window that is not
    already retired -- the honest "nothing left to do" signal.
    """
    workable: dict[str, list[tuple[Breach, Interval, int, int]]] = {}
    priority: dict[str, int] = {}
    for breach in actionable:
        if _breach_key(breach) in retired:
            continue
        pol = resolve(breach.spec_id)
        found = _window_for(
            breach,
            policy=pol.localization,
            window=pol.window,
            motif_len=pol.motif_len,
            codon_map=codon_map,
            construct_length=construct_length,
            circular=circular,
            protein_len=protein_len,
        )
        if found is None:
            continue
        repair_window, first, last = found
        workable.setdefault(breach.spec_id, []).append((breach, repair_window, first, last))
        priority[breach.spec_id] = pol.priority
    if not workable:
        return None

    top = max(priority.values())
    eligible = [spec for spec in workable if priority[spec] == top]
    chosen = min(eligible, key=lambda spec: (turns.get(spec, 0), spec))
    turns[chosen] = turns.get(chosen, 0) + 1
    breach, repair_window, first, last = max(
        workable[chosen], key=lambda t: (t[0].magnitude, -t[0].interval.start)
    )
    return breach, repair_window, first, last, resolve(chosen)


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
    policies: Mapping[str, RulePolicy] | None = None,
    max_candidates: int = 256,
    raise_on_infeasible: bool = True,
    max_iterations: int = MAX_ITERATIONS,
) -> RepairOutcome:
    """Repair HARD_REPAIR breaches without ever weakening a Tier-A guarantee.

    Only breaches marked `fixable_by_codon_choice` are chased; the rest are
    carried on `RepairOutcome.advisory` untouched, because no codon choice can
    move them. Breach count is the cross-rule cost; a rule's breaches are worked
    round-robin so a small-magnitude rule is never starved by a large-magnitude
    one. `policies` supplies a `RulePolicy` per `spec_id` (localization, repair
    discipline, window, motif length, priority); the scalar `policy`/
    `repair_policy`/`window`/`motif_len`/`priority` arguments are the fallback for
    any rule without one.

    On failure the behaviour is the caller's choice. `raise_on_infeasible=True`
    (the default, so every existing caller is unchanged) raises: Infeasible
    Constraints when the actionable breaches were searched out and genuinely
    cannot be cleared, RepairNotConverged when a resource cap stopped the search
    first. `raise_on_infeasible=False` returns the non-clean RepairOutcome with
    the same honest certificate on it instead -- the pipeline uses this so the
    independent validator, not an exception, is what refuses to emit.
    """
    rng = np.random.default_rng(seed)  # explicit; never the global RNG
    patterns = expand_forbidden(forbidden)
    automaton = Automaton(patterns)
    priority = priority or {}
    policies = policies or {}

    def resolve(spec_id: str) -> RulePolicy:
        return policies.get(
            spec_id,
            RulePolicy(policy, repair_policy, window, motif_len, priority.get(spec_id, 0)),
        )

    # A recoded block can create a forbidden motif that runs into the immutable
    # suffix, and a motif of length L needs up to L - 1 suffix bases to do it.
    # The guard is therefore fixed by the longest PATTERN, never by a rule's
    # motif_len: working a 6 nt polyA window while GCGGCCGC (8 nt) is forbidden
    # would otherwise scan only 6 suffix bases and let a junction NotI site slip
    # past the very guarantee Tier A made.
    guard_len = max((len(p) for p in patterns), default=1) - 1

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
    retired: set[tuple[str, int, int]] = set()
    # Abandonment is NOT retirement, and the two are kept apart on purpose.
    # Retirement is SINGLE_PASS saying "this breach has had its attempt"; running
    # out of retired targets is a real `exhausted_targets` and converges.
    # Abandonment is "the search kept missing while the sequence stood still",
    # which is giving up -- so it must never be allowed to report convergence.
    abandoned: set[tuple[str, int, int]] = set()
    misses: dict[tuple[str, int, int], int] = {}
    # The certificate is built from what was ACTUALLY worked last, not from
    # `breaches[0]` and the whole protein. These record it.
    last_target: Breach | None = None
    last_span: tuple[int, int] = (0, len(protein))
    stop_reason = "iterations"  # the for-loop completing without a break
    iteration = 0

    for iteration in range(1, max_iterations + 1):
        codon_map = _codon_map_of(construct, current)
        selected = _select_target(
            actionable,
            resolve=resolve,
            turns=turns,
            retired=retired | abandoned,
            codon_map=codon_map,
            construct_length=construct.length,
            circular=construct.is_circular,
            protein_len=len(protein),
        )
        if selected is None:
            # Honest only if every target got here by retirement. If any was
            # ABANDONED, the search stopped guessing rather than ran out of
            # things to try, and `converged` must stay False.
            stop_reason = "stagnation" if abandoned else "exhausted_targets"
            break  # no actionable breach has a workable, un-retired window
        target, _repair_window, first, last, target_policy = selected
        target_spec = target.spec_id
        last_target = target
        last_span = (first, last)

        options = _mutation_space(protein, code, first, last)
        size = _space_size(options)
        prefix = left_flank + current[: 3 * first]
        suffix = current[3 * last :] + right_flank
        left_state = automaton.consume(0, prefix)[0]

        # A window whose full enumeration exceeds the budget is searched by
        # guided random, NOT enumerated -- so `find_breaches` calls per iteration
        # are bounded by `max_candidates` on both branches, and "empty mutation
        # space" is only ever claimed for a space that was actually enumerated.
        exhaustive = size <= EXHAUSTIVE_LIMIT and size <= max_candidates
        if exhaustive:
            exhaustive_windows += 1
            candidates = _enumerate(options)
        else:
            random_windows += 1
            # Guided random: perturb a FRACTION of the window, never a constant.
            # This branch is only reached when the window was too large to
            # enumerate, and how much has to move to clear a windowed statistic
            # scales with that window -- see `MUTATION_WINDOW_DIVISOR`.
            base = [current[3 * (first + k) : 3 * (first + k) + 3] for k in range(last - first)]
            mutations = min(
                len(base),
                max(MUTATIONS_PER_ITERATION, -(-len(base) // MUTATION_WINDOW_DIVISOR)),
            )
            candidates = []
            for _ in range(max_candidates):
                trial = list(base)
                for pos in rng.choice(len(trial), size=mutations, replace=False):
                    trial[int(pos)] = str(rng.choice(options[int(pos)]))
                candidates.append(tuple(trial))

        improved = False
        for combo in candidates:
            block = "".join(combo)
            if _introduces_forbidden(automaton, left_state, block, suffix[:guard_len]):
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
                        stop_reason="clean",
                    )
                break

        # SINGLE_PASS retires the breach after one attempt: the run continues to
        # other targets but never re-attacks this one. FIXED_POINT leaves it
        # eligible, so round-robin re-targets it until the rule stops producing
        # it -- what splice-site removal requires (a point mutation can create a
        # new cryptic donor, and a single pass would ship it). Retirement is also
        # the perf fix: without it, a hopeless breach is re-attacked every
        # iteration until STAGNATION_TOLERANCE.
        retired_this_iteration = target_policy.repair is RepairPolicy.SINGLE_PASS
        if retired_this_iteration:
            retired.add(_breach_key(target))

        # An improvement moved `current`, so every previous miss was measured
        # against a sequence that no longer exists. Clear the whole ledger: a
        # breach that could not be cleared before may be reachable now, which is
        # exactly the re-targeting FIXED_POINT exists to keep doing. Abandonment
        # therefore only ever accumulates across a STILL sequence.
        if improved:
            misses.clear()
            abandoned.clear()
        elif not retired_this_iteration:
            key = _breach_key(target)
            misses[key] = misses.get(key, 0) + 1
            if misses[key] >= PER_TARGET_TOLERANCE:
                abandoned.add(key)

        # Stagnation is for a FIXED_POINT rule re-targeting a breach it cannot
        # clear. A SINGLE_PASS retirement shrinks the eligible set, so it is
        # progress even without an improvement, and `_select_target` returning
        # None is what ends the pass once every breach has had its attempt.
        stagnant = 0 if (improved or retired_this_iteration) else stagnant + 1
        if stagnant >= STAGNATION_TOLERANCE:
            stop_reason = "stagnation"
            break

    # One exit. If nothing actionable remains the pass is clean regardless of why
    # the loop ended; otherwise build the honest certificate from real state.
    if not actionable:
        return RepairOutcome(
            current,
            iteration,
            True,
            (),
            exhaustive_windows,
            random_windows,
            advisory=advisory,
            stop_reason="clean",
        )

    # `converged` distinguishes "searched the space out and it cannot be cleared"
    # from "a resource cap stopped the search first". The certificate points at
    # the breach worked LAST and the codon window actually searched -- never
    # `breaches[0]` and the whole protein, which named a span the search never
    # touched.
    converged = stop_reason == "exhausted_targets"
    anchor = last_target if last_target is not None else actionable[0]
    certificate = InfeasibilityCertificate(
        interval=anchor.interval,
        protein_span=last_span,
        minimal_conflicting_specs=tuple(sorted({b.spec_id for b in actionable})),
        # The only Tier-B-relevant proof in the frozen Literal. When the search
        # merely gave up (not `converged`), RepairNotConverged is what signals
        # the space was not actually enumerated -- see that class and the RFC it
        # names. A defaulted member cannot be added here without a MAJOR change.
        proof="empty_mutation_space",
    )
    outcome = RepairOutcome(
        current,
        iteration,
        converged,
        tuple(actionable),
        exhaustive_windows,
        random_windows,
        advisory=advisory,
        stop_reason=stop_reason,
        certificate=certificate,
    )
    if not raise_on_infeasible:
        return outcome
    if converged:
        raise InfeasibleConstraints(certificate)
    raise RepairNotConverged(certificate)
