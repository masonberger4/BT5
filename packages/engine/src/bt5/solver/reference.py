"""ReferenceSolver -- deliberately dumb, obviously correct, permanently useful.

This exists for two reasons and both matter:

1. It unblocks every other lane from minute one. A rules author can design and
   test a rule against a real optimizer without waiting for the Tier-A DP.
2. It survives forever as the DIFFERENTIAL ORACLE for the fast solver. When the
   automaton DP and this disagree on a small input, exactly one of them is wrong,
   and this one is small enough to read.

It is greedy with no lookahead, so it is NOT optimal and does not claim to be. It
guarantees only what it checks: the translation is correct and no forbidden motif
appears. When it cannot place a codon without creating one it raises
InfeasibleConstraints with the offending span, rather than emitting a sequence
that violates a hard constraint.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence

from bt5.codon.tables import NcbiGeneticCode
from bt5.core.result import InfeasibilityCertificate, InfeasibleConstraints
from bt5.core.types import Interval, reverse_complement

CodonScorer = Callable[[int, str, str], float]
"""score(codon_index, codon, prefix) -> penalty. LOWER IS BETTER, matching Evaluation.

`prefix` is everything already emitted, including the immutable left flank, so a
scorer can reason about context (repeats, local GC) while staying PURE -- no
mutable state to desynchronise from the search.
"""


def expand_forbidden(patterns: Sequence[str]) -> tuple[str, ...]:
    """Close the pattern set under reverse complement.

    Rules declare FORWARD motifs only. Doing the closure once, here, is what makes
    reverse-strand hits impossible to miss without every rule author remembering
    to think about them.
    """
    out: set[str] = set()
    for p in patterns:
        u = p.upper()
        out.add(u)
        out.add(reverse_complement(u))
    return tuple(sorted(out))


def _creates_forbidden(
    prefix: str, codon: str, forbidden: Sequence[str], right_flank: str = ""
) -> str | None:
    """Return the motif this codon would create, or None.

    Only the tail of the prefix can matter, so the check is bounded by the longest
    pattern rather than by sequence length.
    """
    if not forbidden:
        return None
    longest = max(len(f) for f in forbidden)
    window = (prefix[-(longest - 1) :] if longest > 1 else "") + codon + right_flank
    for motif in forbidden:
        if motif in window:
            return motif
    return None


def back_translate(
    protein: str,
    code: NcbiGeneticCode,
    *,
    forbidden: Sequence[str] = (),
    score: CodonScorer | None = None,
    left_flank: str = "",
    right_flank: str = "",
    add_stop: bool = True,
) -> str:
    """Greedy back-translation that never emits a forbidden motif.

    `left_flank` and `right_flank` are the immutable backbone on either side of
    the insert. Passing them is what makes junction-spanning motifs visible: a
    site formed half by the vector and half by the first codon is caught here,
    not discovered later by the validator.
    """
    patterns = expand_forbidden(forbidden)
    scorer: CodonScorer = score or (lambda _i, _c, _p: 0.0)

    out: list[str] = []
    for i, aa in enumerate(protein):
        options = code.synonymous_codons(aa)
        prefix = left_flank + "".join(out)
        # Only the final emitted codon can interact with the right flank.
        flank = right_flank if (i == len(protein) - 1 and not add_stop) else ""

        viable = [c for c in options if _creates_forbidden(prefix, c, patterns, flank) is None]
        if not viable:
            blocked = {c: _creates_forbidden(prefix, c, patterns, flank) for c in options}
            raise InfeasibleConstraints(
                InfeasibilityCertificate(
                    interval=Interval(3 * i, 3 * i + 3),
                    protein_span=(i, i + 1),
                    minimal_conflicting_specs=tuple(sorted({m for m in blocked.values() if m})),
                    proof="empty_mutation_space",
                )
            )
        # Deterministic tie-break: lowest score, then lexicographic codon. Two runs
        # of one protein must never produce two different tubes.
        out.append(min(viable, key=lambda c: (scorer(i, c, prefix), c)))

    if add_stop:
        prefix = left_flank + "".join(out)
        stops = [
            s
            for s in code.stop_codons
            if _creates_forbidden(prefix, s, patterns, right_flank) is None
        ]
        if not stops:
            raise InfeasibleConstraints(
                InfeasibilityCertificate(
                    interval=Interval(3 * len(protein), 3 * len(protein) + 3),
                    protein_span=(len(protein), len(protein)),
                    minimal_conflicting_specs=("stop_codon",),
                    proof="empty_mutation_space",
                )
            )
        # TAA is the preferred terminator: least read-through-prone of the three.
        out.append("TAA" if "TAA" in stops else sorted(stops)[0])

    return "".join(out)


def cai_scorer(w: Mapping[str, float]) -> CodonScorer:
    """Prefer high relative adaptiveness.

    Deliberately available but NOT the default objective. CAI does not predict
    expression (Kudla 2009: r = 0.14, not significant), and maximising it
    collapses each amino acid to one codon -- see repeat_breaking_scorer.
    """

    def score(_i: int, codon: str, _prefix: str) -> float:
        return -w.get(codon, 0.0)

    return score


def longest_repeat(seq: str, min_len: int = 8) -> tuple[str, int, int] | None:
    """Longest exact substring occurring at least twice: (kmer, first, second).

    Copies MAY overlap, which is the whole point on the sequences this exists to
    measure. One-codon-per-amino-acid back-translation of a repetitive protein
    produces a tandem array, and in a tandem array of period p the copies sit p
    apart and overlap whenever the repeat is longer than p. Bounding the search
    at n // 2 -- as though two copies had to be disjoint -- silently reports the
    bound instead of the answer: max-CAI on a (G4S)3 linker gives a true 30 bp
    repeat in 45 bp, and a n // 2 search calls it 22.

    Searched by bisection rather than by counting down. "A repeat of length L
    exists" is monotone -- a prefix of a repeat is a repeat at the same two
    starts -- so bisection finds the largest true L in O(n log n) probes instead
    of walking every length from n down, which is what the old ceiling was
    really buying. This is the reference implementation, so it is allowed to be
    slow; it is not allowed to be wrong.
    """
    n = len(seq)

    def occurring_twice(size: int) -> tuple[str, int, int] | None:
        seen: dict[str, int] = {}
        for i in range(n - size + 1):
            k = seq[i : i + size]
            if k in seen:
                return (k, seen[k], i)
            seen[k] = i
        return None

    best: tuple[str, int, int] | None = None
    lo, hi = min_len, n - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        hit = occurring_twice(mid)
        if hit is None:
            hi = mid - 1
        else:
            best = hit
            lo = mid + 1
    return best


def repeat_breaking_scorer(
    w: Mapping[str, float],
    *,
    kmer: int = 9,
    repeat_penalty: float = 100.0,
    usage_weight: float = 1.0,
) -> CodonScorer:
    """Penalise a codon that extends an exact repeat, ABOVE codon preference.

    This is the concrete form of a rule the evidence makes non-negotiable:
    repeat-breaking OVERRIDES codon optimality, it does not compete with it.

    Maximising CAI collapses each amino acid to a single codon, so naive
    back-translation of any repeat-containing protein produces PERFECT nucleotide
    repeats -- and repetitive 9-mers per 100 bp plus longest-repeat length are the
    two highest-importance features in the best published synthesis-success model
    (random forest over 1,076 real vendor outcomes, F1 0.928,
    https://pubs.acs.org/doi/10.1021/acssynbio.9b00460). Repeats also cause Gibson
    misassembly, and below ~200 bp deletion in E. coli is RecA-INDEPENDENT, so a
    recA- strain gives no protection against exactly this size class.

    The proteins this matters for are the ones people actually express:
    antibodies and scFv/CAR constructs, Fc fusions, (GGGGS)n linkers, tandem 2A
    peptides, TALEs, zinc fingers, and His tags. Twist calls out encoding His6 as
    alternating CAC/CAT by name.

    The penalty dominates `usage_weight` (100 vs 1) so a rare-but-unique codon is
    always preferred to a preferred-but-repeating one.
    """

    def score(_i: int, codon: str, prefix: str) -> float:
        candidate = prefix + codon
        penalty = 0.0
        if len(candidate) >= kmer:
            tail = candidate[-kmer:]
            # Occurs earlier in the sequence, not counting this occurrence?
            if tail in candidate[:-1]:
                penalty += repeat_penalty
        return penalty - usage_weight * w.get(codon, 0.0)

    return score
