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


def no_rules(_c: Construct) -> tuple[Breach, ...]:
    """The EXPLICIT opt-out from Tier B.

    `optimize()` used to accept `find_breaches=None` and fabricate
    `RepairOutcome(cds, 0, converged=True)`, which is byte-identical to "repair
    ran and found nothing to fix" -- so a caller who forgot the argument got a
    construct that looked proven and had never been checked. Passing this is the
    same behaviour said out loud, and `optimize()` reports `ran is False` for it.
    """
    return ()


NO_RULES: BreachFinder = no_rules


@dataclass(frozen=True)
class RepairOutcome:
    cds: str
    iterations: int
    converged: bool
    remaining: tuple[Breach, ...] = ()
    exhaustive_windows: int = 0
    random_windows: int = 0
    #: Did Tier B actually run? REQUIRED to distinguish "repaired, nothing to
    #: fix" from "never checked". `optimize()` used to fabricate
    #: `RepairOutcome(cds, 0, True)` when no breach finder was supplied, which is
    #: byte-identical to a clean repair -- so a caller who forgot the argument
    #: got a construct that looked proven and had never been looked at. The two
    #: states are now different values, not the same one.
    ran: bool = True
    #: The policy the loop actually ran under, after any escalation. Recorded
    #: because the escalation is a safety decision made somewhere else -- a test
    #: that cannot see it cannot prove the join went the right way.
    effective_repair_policy: RepairPolicy = RepairPolicy.FIXED_POINT

    @property
    def clean(self) -> bool:
        return not self.remaining


#: A repair cost. Compared lexicographically, never summed. See `BreachCost`.
Cost = tuple[float, float, float, float]


class BreachCost:
    """The repair objective: a lexicographic tuple, with a scale frozen at first use.

    This loop used to minimise `sum(b.magnitude)`, but `Breach.magnitude` is
    documented rule-native: E1's is nucleotides over a run limit, E9's is a
    dimensionless 1.0/0.5/0.25 tier, E2's is a GC fraction delta. Summing them
    makes one nucleotide worth thirty-three percentage points of GC, so the
    search chases whichever rule happens to report in the largest units -- and
    the same arithmetic picked the breach to work next, so the rule with the
    biggest units also monopolised the search.

    Four terms, compared in order:

      1. how many breaches -- already unit-free, and clearing one is
         unambiguously better than shrinking one;
      2. summed magnitude, each divided by its own rule's scale;
      3. how much sequence is implicated, preferring the more localised of two
         otherwise equal states;
      4. leftmost breach, so ties break deterministically.

    THE SCALE IS FROZEN, and that is the load-bearing part. Normalising by a
    moving maximum makes the objective non-stationary: the same two breach sets
    compare one way at iteration 3 and the other way at iteration 40, so
    `best_cost` stops being a monotone bound and the stagnation counter measures
    nothing. The scale comes from the first call -- the initial breach set,
    before any candidate -- and never moves. A rule first seen later gets 1.0 and
    is recorded in `late` rather than silently rescaling the run.
    """

    __slots__ = ("_scale", "late")

    def __init__(self) -> None:
        self._scale: dict[str, float] | None = None
        self.late: tuple[str, ...] = ()

    def normalised(self, breach: Breach) -> float:
        """`breach.magnitude` in units of its own rule's initial worst case."""
        if self._scale is None:
            return abs(breach.magnitude)
        return abs(breach.magnitude) * self._scale.get(breach.spec_id, 1.0)

    def __call__(self, breaches: Sequence[Breach]) -> Cost:
        if self._scale is None:
            self._scale = _initial_scale(breaches)
        else:
            unseen = {b.spec_id for b in breaches} - set(self._scale)
            if unseen:
                self.late = tuple(sorted(set(self.late) | unseen))
        if not breaches:
            return (0.0, 0.0, 0.0, 0.0)
        return (
            float(len(breaches)),
            # Rounded so the epsilon the scalar comparison carried survives into
            # a lexicographic one: without it two states differing by 1e-16 of
            # normalised magnitude read as an improvement forever.
            round(sum(self.normalised(b) for b in breaches), 12),
            float(sum(b.interval.length for b in breaches)),
            float(min(b.interval.start for b in breaches)),
        )


def _initial_scale(breaches: Sequence[Breach]) -> dict[str, float]:
    """1 / (this rule's worst magnitude in the initial set), per rule."""
    worst: dict[str, float] = {}
    for b in breaches:
        worst[b.spec_id] = max(worst.get(b.spec_id, 0.0), abs(b.magnitude))
    return {spec_id: (1.0 / m if m > 0.0 else 1.0) for spec_id, m in worst.items()}


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


def _normalised(measure: Callable[[Sequence[Breach]], Cost], breach: Breach) -> float:
    """A breach's magnitude in its own rule's units, if the cost object knows how.

    A caller may inject any callable as the cost, so the scale is read through a
    duck-typed hook rather than requiring `BreachCost` specifically. Falling back
    to the raw magnitude reproduces the old ordering exactly, which is the right
    default for a cost function that does not model rules at all.
    """
    normalise = getattr(measure, "normalised", None)
    if normalise is None:
        return abs(breach.magnitude)
    return float(normalise(breach))


def repair(
    cds: str,
    protein: str,
    code: NcbiGeneticCode,
    *,
    assemble: Assembler,
    find_breaches: BreachFinder,
    forbidden: Sequence[str] = (),
    policy: LocalizationPolicy | Callable[[str], LocalizationPolicy] = (
        LocalizationPolicy.WINDOW_MINUS_1
    ),
    repair_policy: RepairPolicy = RepairPolicy.FIXED_POINT,
    cost: Callable[[Sequence[Breach]], Cost] | None = None,
    window: int = 50,
    motif_len: int = 6,
    left_flank: str = "",
    right_flank: str = "",
    seed: int = 0,
    max_iterations: int = MAX_ITERATIONS,
) -> RepairOutcome:
    """Repair HARD_REPAIR breaches without ever weakening a Tier-A guarantee.

    Raises InfeasibleConstraints when the search cannot clear the breaches, with
    the minimal conflicting set rather than a bare failure.

    `policy` may be one `LocalizationPolicy` for every breach, or a callable
    resolving one PER RULE. The catalog needs the second form: the four
    HARD_REPAIR rules shipped today declare three different policies between
    them -- E2 `WINDOW_MINUS_1`, E5 and F1 `PAIRED_SEGMENTS`, E7 `WHOLE_SCOPE` --
    and a repeat pair widened by `window - 1` instead of localised to the copy
    being edited searches a region that does not contain the fix.
    """
    rng = np.random.default_rng(seed)  # explicit; never the global RNG
    patterns = expand_forbidden(forbidden)
    automaton = Automaton(patterns)
    measure = cost if cost is not None else BreachCost()
    localization = policy if callable(policy) else (lambda _spec_id: policy)

    current = cds
    construct = assemble(current)
    breaches = find_breaches(construct)
    if not breaches:
        return RepairOutcome(current, 0, True, effective_repair_policy=repair_policy)

    stagnant = 0
    exhaustive_windows = random_windows = 0
    best_cost = measure(breaches)

    for iteration in range(1, max_iterations + 1):
        # Work the worst breach first; deterministic tie-break on position.
        # NORMALISED, not raw: ranking by `magnitude` across rules compares a
        # nucleotide count against a GC fraction, so the rule reporting in the
        # largest units would be the only one ever worked on.
        target = max(
            breaches,
            key=lambda b: (_normalised(measure, b), -b.interval.start),
        )
        repair_window = localize(
            target,
            localization(target.spec_id),
            window=window,
            motif_len=motif_len,
            construct_length=construct.length,
            circular=construct.is_circular,
        )
        codon_map = (
            construct.translation_units[0].codon_map
            if construct.translation_units
            else tuple(Interval(i, i + 3) for i in range(0, len(current), 3))
        )
        first, last = codon_span(repair_window, codon_map, construct.length)
        last = min(last, len(protein))
        if first >= last:
            break

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
            trial_breaches = find_breaches(assemble(trial_cds))
            trial_cost = measure(trial_breaches)
            if trial_cost < best_cost:
                current, best_cost, breaches = trial_cds, trial_cost, trial_breaches
                construct = assemble(current)
                improved = True
                if not breaches:
                    return RepairOutcome(
                        current,
                        iteration,
                        True,
                        (),
                        exhaustive_windows,
                        random_windows,
                        effective_repair_policy=repair_policy,
                    )
                break

        stagnant = 0 if improved else stagnant + 1
        if stagnant >= STAGNATION_TOLERANCE:
            break
        if repair_policy is RepairPolicy.SINGLE_PASS and not improved:
            break

    if breaches:
        raise InfeasibleConstraints(
            InfeasibilityCertificate(
                interval=breaches[0].interval,
                protein_span=(0, len(protein)),
                minimal_conflicting_specs=tuple(sorted({b.spec_id for b in breaches})),
                proof="empty_mutation_space",
            )
        )
    return RepairOutcome(
        current,
        max_iterations,
        True,
        (),
        exhaustive_windows,
        random_windows,
        effective_repair_policy=repair_policy,
    )
