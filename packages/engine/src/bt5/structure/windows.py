"""Where to fold, decided without folding anything.

Kept apart from the engine on purpose. Window geometry is pure arithmetic over
intervals -- origin wraps, strand, whether a window fits at all -- and it is
where the mistakes actually live. Folding the right sequence at the wrong offset
produces a number that looks entirely reasonable, so this half is tested on its
own, with no ViennaRNA in the room.

The 5' window is the one that matters most. Kudla measured a -4..+37 window that
explains 44-59% of expression variance in bacteria, and that window SPANS THE
UTR/CDS JUNCTION: it cannot be computed from the CDS alone, which is why
`five_prime_window` takes the whole construct's geometry rather than a coding
sequence. A detector that quietly folds the CDS by itself and reports the number
as though the leader were there is the exact failure this module exists to
prevent, and the reason it returns None rather than a truncated window when the
upstream bases are not available.
"""

from __future__ import annotations

from bt5.core.types import Interval, Strand

#: Kudla 2009's window, in codon-adjacent coordinates relative to the A of ATG.
#: 4 bases of leader and 37 of coding sequence. https://pmc.ncbi.nlm.nih.gov/articles/PMC3902468/
KUDLA_UPSTREAM = 4
KUDLA_DOWNSTREAM = 37

#: Tier C's sweep, from docs/PLAN.md. 100 nt is what the plan specifies; note
#: that MFE is O(n^3), so the window width dominates the cost of any sweep --
#: 30 nt folds in 0.14 ms and 100 nt in 4.9 ms, a factor of 34 for a factor of
#: 3.3 in width. Anything that sweeps should say which width it chose and why.
SWEEP_SIZE = 100
SWEEP_STEP = 10


def five_prime_window(
    cds: Interval,
    *,
    length: int,
    circular: bool,
    upstream: int = KUDLA_UPSTREAM,
    downstream: int = KUDLA_DOWNSTREAM,
) -> Interval | None:
    """The window straddling the start codon, following the transcript.

    Returns None when the upstream bases do not exist -- a linear construct
    whose CDS starts too near an end. None is the honest answer there: a window
    silently clamped to what fits is a DIFFERENT window, and comparing it to one
    that had its full leader is comparing two different measurements.

    Strand-aware, because for a reverse-oriented cassette the 5' side of the CDS
    sits at HIGHER construct coordinates. Getting this backwards folds the 3' end
    and reports it as the 5' term, which is wrong in a way no unit looks odd.
    """
    if upstream < 0 or downstream < 1:
        raise ValueError(f"window must extend into the CDS: {upstream=} {downstream=}")
    if cds.length < downstream:
        return None

    strand: Strand = cds.strand
    start = cds.start - upstream if strand == 1 else cds.end - downstream
    if start < 0:
        if not circular:
            return None
        start += length
    if not circular and start + upstream + downstream > length:
        return None
    return Interval(start, start + upstream + downstream, strand)


def sliding_windows(
    length: int, *, size: int = SWEEP_SIZE, step: int = SWEEP_STEP, circular: bool = False
) -> tuple[Interval, ...]:
    """Tile a construct in overlapping windows.

    On a circular construct the tiling continues past the origin and the last
    windows wrap, so a structure sitting across position 0 is folded like any
    other. On a linear one the sweep stops at the last window that fits whole,
    because a short final window is not comparable to the others.
    """
    if size < 1 or step < 1:
        raise ValueError(f"{size=} and {step=} must both be positive")
    if circular:
        return tuple(Interval(i, i + size) for i in range(0, length, step))
    return tuple(Interval(i, i + size) for i in range(0, max(0, length - size + 1), step))


def windows_touching(
    changed: Interval, plan: tuple[Interval, ...], *, length: int, circular: bool
) -> tuple[int, ...]:
    """Indices of windows a change invalidates -- the incremental-update primitive.

    A single codon substitution changes a handful of windows and leaves the rest
    valid, which is the whole reason a sweep is affordable inside a search loop.
    Recomputing everything instead is correct and roughly a hundred times slower.
    """
    out: list[int] = []
    for i, w in enumerate(plan):
        if _overlaps(w, changed, length=length, circular=circular):
            out.append(i)
    return tuple(out)


def _overlaps(a: Interval, b: Interval, *, length: int, circular: bool) -> bool:
    """Wrap-aware overlap: either interval may run past the end of the sequence."""
    shifts = (0, length, -length) if circular else (0,)
    return any(a.start < b.end + s and b.start + s < a.end for s in shifts)
