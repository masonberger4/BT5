"""Tier-A: exact Viterbi DP over the codon lattice.

This is what turns "no forbidden motif" from CHECKED into GUARANTEED BY
CONSTRUCTION. The search space is every synonymous back-translation of the
protein; the DP finds the best-scoring one that never passes through an accepting
state of an Aho-Corasick automaton built over the forbidden set.

Three consequences fall out for free, and each is a bug class that a per-codon
scan cannot see:

  - motifs spanning a CODON BOUNDARY, because the automaton state carries the
    match prefix across the boundary;
  - motifs on the REVERSE STRAND, because the pattern set is closed under reverse
    complement before the automaton is built;
  - motifs spanning the CDS/BACKBONE JUNCTION, because the automaton is SEEDED
    with the state reached after consuming the immutable left flank, and the
    final state is required to survive the right flank.

Complexity is O(L x |S| x d) where L is codon count, |S| automaton states and d
the mean synonymous degeneracy (~2.9). For a 1000-codon protein against ~50
six-base motifs that is a few million relaxations -- milliseconds, not seconds.

What this tier does NOT do: windowed GC content. Deciding whether a codon pushes
a 50 bp window past its bound needs the G+C count over the previous ~17 codons,
and enumerating that history is intractable. Windowed GC is steered by a
Lagrangian term here, repaired in Tier B, and finally PROVEN by the independent
validator, which refuses to emit if a window is still out of band.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from bt5.codon.tables import NcbiGeneticCode
from bt5.core.result import InfeasibilityCertificate, InfeasibleConstraints
from bt5.core.types import Interval
from bt5.solver.reference import CodonScorer, expand_forbidden

BASES = "ACGT"
_BASE_INDEX = {b: i for i, b in enumerate(BASES)}


class Automaton:
    """Aho-Corasick with a dense goto table over {A,C,G,T}.

    Built here rather than taken from a library because the DP needs the explicit
    transition function, not just a search interface.
    """

    __slots__ = ("accepting", "goto", "n_states", "patterns")

    def __init__(self, patterns: Sequence[str]) -> None:
        self.patterns = tuple(patterns)
        # --- trie ---
        goto: list[list[int]] = [[-1] * 4]
        accepting: list[bool] = [False]
        for pattern in patterns:
            state = 0
            for ch in pattern:
                idx = _BASE_INDEX[ch]
                if goto[state][idx] == -1:
                    goto.append([-1] * 4)
                    accepting.append(False)
                    goto[state][idx] = len(goto) - 1
                state = goto[state][idx]
            accepting[state] = True

        # --- failure links, converted into a complete DFA ---
        fail = [0] * len(goto)
        queue: list[int] = []
        for idx in range(4):
            nxt = goto[0][idx]
            if nxt == -1:
                goto[0][idx] = 0
            else:
                fail[nxt] = 0
                queue.append(nxt)
        head = 0
        while head < len(queue):
            state = queue[head]
            head += 1
            # A state is accepting if it is, or if any suffix of it is.
            accepting[state] = accepting[state] or accepting[fail[state]]
            for idx in range(4):
                nxt = goto[state][idx]
                if nxt == -1:
                    goto[state][idx] = goto[fail[state]][idx]
                else:
                    fail[nxt] = goto[fail[state]][idx]
                    queue.append(nxt)

        self.goto = goto
        self.accepting = accepting
        self.n_states = len(goto)

    def step(self, state: int, base: str) -> int:
        return self.goto[state][_BASE_INDEX[base]]

    def consume(self, state: int, text: str) -> tuple[int, bool]:
        """Advance through `text`. Returns (end state, whether a match occurred)."""
        hit = False
        for ch in text:
            state = self.goto[state][_BASE_INDEX[ch]]
            if self.accepting[state]:
                hit = True
        return state, hit


def _codon_transitions(
    automaton: Automaton, codons: Sequence[str]
) -> dict[str, list[tuple[int, bool]]]:
    """Precompute, per codon, the (next state, forbidden?) for every start state.

    Done once so the inner DP loop is a table lookup rather than three automaton
    steps per relaxation.
    """
    table: dict[str, list[tuple[int, bool]]] = {}
    for codon in codons:
        row: list[tuple[int, bool]] = []
        for state in range(automaton.n_states):
            row.append(automaton.consume(state, codon))
        table[codon] = row
    return table


def optimal_back_translate(
    protein: str,
    code: NcbiGeneticCode,
    *,
    forbidden: Sequence[str] = (),
    score: CodonScorer | None = None,
    left_flank: str = "",
    right_flank: str = "",
    add_stop: bool = True,
) -> str:
    """Exact minimum-cost back-translation avoiding every forbidden motif.

    `score` is called as (codon_index, codon, "") -- the DP evaluates codons
    independently, so a scorer that depends on the emitted prefix (the
    repeat-breaking scorer) is NOT exactly optimisable here and belongs in Tier B.
    Passing one is allowed but it sees an empty prefix, so it degrades to a
    context-free preference rather than silently reporting a wrong optimum.
    """
    patterns = expand_forbidden(forbidden)
    scorer: CodonScorer = score or (lambda _i, _c, _p: 0.0)

    alphabet: set[str] = set()
    for aa in protein:
        alphabet.update(code.synonymous_codons(aa))
    if add_stop:
        alphabet.update(code.stop_codons)

    automaton = Automaton(patterns) if patterns else Automaton(())
    trans = _codon_transitions(automaton, sorted(alphabet))

    # Seed with the immutable backbone to the left: a site formed half by the
    # vector and half by the first codon is excluded by construction.
    start_state, flank_hit = automaton.consume(0, left_flank)
    if flank_hit and patterns:
        # The user's own backbone contains a forbidden motif. That is a finding
        # about the vector, not something codon choice can fix.
        raise InfeasibleConstraints(
            InfeasibilityCertificate(
                interval=Interval(0, max(1, len(left_flank))),
                protein_span=(0, 0),
                minimal_conflicting_specs=("backbone_contains_forbidden_motif",),
                proof="immutable_region",
            )
        )

    inf = math.inf
    # dp[state] -> best cost to reach that automaton state
    dp: dict[int, float] = {start_state: 0.0}
    choices: list[dict[int, tuple[int, str]]] = []  # per position: state -> (prev state, codon)

    units = list(protein) + (["*"] if add_stop else [])
    for i, aa in enumerate(units):
        options = code.stop_codons if aa == "*" else code.synonymous_codons(aa)
        nxt: dict[int, float] = {}
        back: dict[int, tuple[int, str]] = {}
        for state, cost in dp.items():
            for codon in options:
                to_state, hit = trans[codon][state]
                if hit:
                    continue  # forbidden by construction: this path does not exist
                step_cost = cost + scorer(i, codon, "")
                if step_cost < nxt.get(to_state, inf):
                    nxt[to_state] = step_cost
                    back[to_state] = (state, codon)
        if not nxt:
            raise InfeasibleConstraints(
                InfeasibilityCertificate(
                    interval=Interval(3 * i, 3 * i + 3),
                    protein_span=(i, i + 1),
                    minimal_conflicting_specs=patterns or ("unknown",),
                    proof="automaton_dead_state",
                )
            )
        choices.append(back)
        dp = nxt

    # The right flank must not complete a motif either.
    survivors = {
        state: cost for state, cost in dp.items() if not automaton.consume(state, right_flank)[1]
    }
    if not survivors:
        raise InfeasibleConstraints(
            InfeasibilityCertificate(
                interval=Interval(3 * len(units) - 3, 3 * len(units)),
                protein_span=(len(protein), len(protein)),
                minimal_conflicting_specs=("right_flank_junction",),
                proof="automaton_dead_state",
            )
        )

    # Deterministic tie-break so two runs never produce two different tubes.
    end_state = min(survivors, key=lambda s: (survivors[s], s))
    out: list[str] = []
    state = end_state
    for back in reversed(choices):
        prev, codon = back[state]
        out.append(codon)
        state = prev
    return "".join(reversed(out))


def cai_lattice_scorer(w: Mapping[str, float]) -> CodonScorer:
    """Additive log-w cost, so the DP optimum is the exact maximum-CAI sequence.

    CAI is a geometric mean, so summing -log(w) is what makes it optimisable by a
    shortest-path DP at all. Shipped for differential testing against a known
    optimum; it is NOT the default objective, for the reasons in tables.py.
    """
    floor = math.log(1e-6)

    def score(_i: int, codon: str, _prefix: str) -> float:
        weight = w.get(codon, 0.0)
        return -(math.log(weight) if weight > 0 else floor)

    return score


def gc_steering_scorer(multiplier: float) -> CodonScorer:
    """Signed Lagrangian term on GC content, for Tier A.

    This is a genuine Lagrangian relaxation, not a distance-to-target penalty,
    and the difference is not cosmetic. Penalising ``abs(gc - 0.5)`` per codon is
    SYMMETRIC: for a GC-rich protein whose synonymous families all sit at 2/3 or
    3/3 GC, it pulls toward 2/3 rather than downward, and can raise global GC
    while appearing to "steer" it. That was measured, not theorised.

    The correct form is a signed multiplier on the codon's own GC:

        cost = lambda * gc(codon)

    A POSITIVE multiplier makes GC expensive and drives composition DOWN; a
    NEGATIVE multiplier drives it UP. The multiplier is chosen by
    `solve_with_gc_steering` rather than guessed.

    This remains STEERING, not enforcement. The term is context-free -- it sees a
    codon, never a window -- so it moves global composition and thereby reduces
    Tier B's workload. Only the repair pass plus the independent validator make
    any particular window's bound hard.
    """

    def score(_i: int, codon: str, _prefix: str) -> float:
        return multiplier * (codon.count("G") + codon.count("C")) / 3.0

    return score


def combine_scorers(*scorers: CodonScorer) -> CodonScorer:
    """Sum several codon scorers. Used to add GC steering to a usage preference."""

    def score(i: int, codon: str, prefix: str) -> float:
        return sum(s(i, codon, prefix) for s in scorers)

    return score


def achievable_gc_range(
    protein: str,
    code: NcbiGeneticCode,
    *,
    forbidden: Sequence[str] = (),
    left_flank: str = "",
    right_flank: str = "",
) -> tuple[float, float]:
    """The GC fraction range this protein can actually reach, given its sequence.

    Amino-acid composition, not codon choice, is often the binding constraint. A
    peptide of Lys/Asn/Ile/Phe has at most one GC base per codon, so ~33% is its
    hard ceiling however expensive GC is made; Ala/Gly/Pro/Arg have the opposite
    floor. Establishing the range with two pure-objective solves is exact, costs
    a couple of milliseconds, and turns "the search failed" into the reportable
    fact "your target is outside this protein's achievable range of X-Y%".
    """
    lowest = optimal_back_translate(
        protein,
        code,
        forbidden=forbidden,
        score=gc_steering_scorer(1.0),
        left_flank=left_flank,
        right_flank=right_flank,
    )
    highest = optimal_back_translate(
        protein,
        code,
        forbidden=forbidden,
        score=gc_steering_scorer(-1.0),
        left_flank=left_flank,
        right_flank=right_flank,
    )

    def gc_of(seq: str) -> float:
        return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0

    return gc_of(lowest), gc_of(highest)


def solve_with_gc_steering(
    protein: str,
    code: NcbiGeneticCode,
    *,
    gc_bounds: tuple[float, float],
    base_score: CodonScorer | None = None,
    forbidden: Sequence[str] = (),
    left_flank: str = "",
    right_flank: str = "",
    max_steps: int = 24,
) -> tuple[str, float]:
    """Choose the Lagrangian multiplier by bisection, then solve.

    The multiplier is solved for, not guessed. Each Tier-A solve is milliseconds,
    so a couple of dozen is still far inside the interactive budget, and every
    result is the exact optimum FOR ITS MULTIPLIER rather than a heuristic nudge.

    Feasibility is settled up front by `achievable_gc_range`, because a bisection
    that cannot succeed would otherwise spend its whole budget discovering that.
    When the band lies outside the protein's reachable range, the closest
    attainable sequence is returned and the shortfall is left to the validator to
    refuse -- Tier B cannot fix a limit imposed by the amino-acid sequence itself.

    Returns (cds, multiplier).
    """
    lo_gc, hi_gc = gc_bounds

    def solve(mult: float) -> str:
        steer = gc_steering_scorer(mult)
        combined = combine_scorers(base_score, steer) if base_score else steer
        return optimal_back_translate(
            protein,
            code,
            forbidden=forbidden,
            score=combined,
            left_flank=left_flank,
            right_flank=right_flank,
        )

    def gc_of(seq: str) -> float:
        return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0

    unsteered = solve(0.0)
    if lo_gc <= gc_of(unsteered) <= hi_gc:
        return unsteered, 0.0

    reachable_lo, reachable_hi = achievable_gc_range(
        protein,
        code,
        forbidden=forbidden,
        left_flank=left_flank,
        right_flank=right_flank,
    )
    if reachable_lo > hi_gc:
        # Even the most GC-poor encoding is above the ceiling.
        return solve(1e6), 1e6
    if reachable_hi < lo_gc:
        # Even the most GC-rich encoding is below the floor.
        return solve(-1e6), -1e6

    # A solution exists inside the reachable range, so bracket and bisect.
    direction = 1.0 if gc_of(unsteered) > hi_gc else -1.0
    low, high = 0.0, direction
    for _ in range(max_steps):
        candidate = solve(high)
        candidate_gc = gc_of(candidate)
        if lo_gc <= candidate_gc <= hi_gc:
            return candidate, high
        if (direction > 0 and candidate_gc < lo_gc) or (direction < 0 and candidate_gc > hi_gc):
            break  # overshot: the band lies between `low` and `high`
        low, high = high, high * 2.0
        if abs(high) > 1e6:
            return candidate, high

    best_seq, best_mult = solve(high), high
    for _ in range(max_steps):
        mid = (low + high) / 2.0
        candidate = solve(mid)
        candidate_gc = gc_of(candidate)
        best_seq, best_mult = candidate, mid
        if lo_gc <= candidate_gc <= hi_gc:
            return candidate, mid
        if (direction > 0 and candidate_gc > hi_gc) or (direction < 0 and candidate_gc < lo_gc):
            low = mid
        else:
            high = mid
    return best_seq, best_mult
