"""How different two designs actually are.

Gate G4 asks for a gallery of candidates at least 15% apart in CODON space, and
the choice of space is the whole point. Two sequences can differ at 30% of their
BASES while encoding their differences as third-position wobble that changes
nothing a ribosome or a synthesiser cares about, and a gallery selected on base
distance will happily present five of those as five options. Counting codons that
differ is the measure that tracks what the user is actually choosing between.

Distances here are over sequences encoding the SAME protein. That is not a
limitation to work around -- it is what makes the numbers comparable, because two
CDSs of the same protein have the same codon count and a position-by-position
comparison is meaningful.
"""

from __future__ import annotations

from collections.abc import Sequence


def codon_distance(a: str, b: str) -> float:
    """Fraction of codon positions at which two CDSs differ, in [0, 1].

    Raises on a length mismatch rather than truncating. Two sequences of
    different lengths do not encode the same protein, so a distance between them
    is not the quantity anything here wants, and silently comparing the shorter
    prefix would return a plausible number for an incomparable pair.
    """
    if len(a) != len(b):
        raise ValueError(f"cannot compare CDSs of different lengths: {len(a)} vs {len(b)}")
    if len(a) % 3:
        raise ValueError(f"not a whole number of codons: {len(a)} nt")
    if not a:
        return 0.0
    n = len(a) // 3
    differing = sum(1 for i in range(n) if a[3 * i : 3 * i + 3] != b[3 * i : 3 * i + 3])
    return differing / n


def nucleotide_distance(a: str, b: str) -> float:
    """Fraction of BASES that differ. Reported alongside, never selected on.

    Kept because it is what a user reading two sequences side by side perceives,
    and because the gap between it and the codon distance is itself informative:
    a large nucleotide distance with a small codon distance means the differences
    are wobble.
    """
    if len(a) != len(b):
        raise ValueError(f"cannot compare sequences of different lengths: {len(a)} vs {len(b)}")
    if not a:
        return 0.0
    return sum(1 for x, y in zip(a, b, strict=True) if x != y) / len(a)


def pairwise_minimum(sequences: Sequence[str]) -> float:
    """The closest any two designs come. This is the number G4 gates on.

    The MINIMUM rather than the mean, because a gallery is only as diverse as its
    most similar pair: four genuinely different designs plus one near-duplicate
    of another is a gallery of four, and a mean would hide that.
    """
    if len(sequences) < 2:
        return 1.0
    return min(
        codon_distance(sequences[i], sequences[j])
        for i in range(len(sequences))
        for j in range(i + 1, len(sequences))
    )


def distance_matrix(sequences: Sequence[str]) -> tuple[tuple[float, ...], ...]:
    """Full symmetric matrix, for `Candidate.codon_distance_to` and the report."""
    n = len(sequences)
    rows: list[tuple[float, ...]] = []
    cache: dict[tuple[int, int], float] = {}
    for i in range(n):
        row: list[float] = []
        for j in range(n):
            if i == j:
                row.append(0.0)
                continue
            key = (min(i, j), max(i, j))
            if key not in cache:
                cache[key] = codon_distance(sequences[key[0]], sequences[key[1]])
            row.append(cache[key])
        rows.append(tuple(row))
    return tuple(rows)
