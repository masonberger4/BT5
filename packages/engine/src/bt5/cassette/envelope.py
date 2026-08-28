"""What this protein CAN and CANNOT do, computed before any codon is chosen.

BT5 currently optimises toward a target and, if the target was never reachable,
reports a conflict at the end. That is backwards. The amino acid sequence fixes
the achievable range before the solver runs, so the app can open with what the
options are instead of closing with why it failed.

**The arithmetic is small and the result is exact, per window.** A window's GC
count is a SUM of per-codon contributions, and codon choices are independent
across positions, so taking each overlapping codon's cheapest contribution and
adding them gives the true minimum for that window -- not an estimate.

**But per-window exact is not jointly achievable, and the difference is the
whole honesty story.** Two windows sharing codons can each be satisfiable alone
and not together: minimising window A may be exactly what pushes window B up. So

    a target OUTSIDE a window's range is PROVEN unreachable;
    a target INSIDE every window's range is NOT PROVEN reachable.

The envelope can say no with certainty and maybe without pretending, which is
the same discipline `native_baseline` and the percentile-not-prediction rule
exist to enforce elsewhere.

**What it shows in practice.** Gly, Pro and Ala are GGN, CCN and GCN, so their
first two bases are GC whatever you choose; Trp is TGG and Met is ATG, with no
choice at all. Composed over a real sequence:

    ordinary 300 aa protein     27% - 60% achievable GC
    His6                        33% - 67%
    (GGGGS)x3 linker            60% - 93%

**A (GGGGS)3 linker cannot be brought below 60% GC by any codon assignment.**
That is the textbook repeat motif, and its hard constraint turns out to be GC
rather than repetition -- see `forced_repeat_bp`, which is almost always
negligible because repeats ARE recodeable and GC is not.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache

from bt5.codon.tables import NcbiGeneticCode

#: Windowed GC is conventionally taken at 50 bp here, matching E2/E4 and Twist's
#: published window size, so an envelope and a rule finding describe the same
#: geometry rather than two nearly-identical ones.
DEFAULT_WINDOW_BP = 50

#: How many responsible residues to name per bound. Enough to point at a linker
#: or a tag; not so many that the message becomes a listing.
NAMED_RESIDUES = 4


@dataclass(frozen=True, slots=True)
class Culprit:
    """One residue span holding a bound where it is, in PROTEIN coordinates.

    Protein coordinates rather than nucleotide ones deliberately: the actions a
    user can take -- shorten the linker, swap the tag, pick a different fusion
    partner -- are all protein edits, and handing them a nucleotide offset makes
    them do the division themselves.
    """

    start: int
    end: int
    residues: str
    #: GC percentage points this span contributes to the bound it holds.
    contribution: float


@dataclass(frozen=True, slots=True)
class GCWindow:
    """The exact achievable GC range of one window, and what pins each end."""

    start: int  # nucleotide offset into the CDS, 0-based, half-open
    end: int
    lowest: float  # percent
    highest: float
    holds_low: tuple[Culprit, ...] = ()
    holds_high: tuple[Culprit, ...] = ()

    def admits(self, gc_min: float, gc_max: float) -> bool:
        """Whether ANY codon assignment puts this window inside [gc_min, gc_max]."""
        return self.lowest <= gc_max and self.highest >= gc_min

    @property
    def span(self) -> float:
        return self.highest - self.lowest


@dataclass(frozen=True, slots=True)
class FeasibleEnvelope:
    """Everything the protein decides before the solver is allowed an opinion."""

    protein: str
    table_id: int
    window: int
    windows: tuple[GCWindow, ...]
    #: Achievable GC over the whole CDS, exact for the same reason.
    lowest: float
    highest: float
    #: Longest DNA run forced identical between two places by the protein alone,
    #: in bp. Near-zero for almost every real protein -- see `_forced_repeat`.
    forced_repeat_bp: int

    def unreachable(self, gc_min: float, gc_max: float) -> tuple[GCWindow, ...]:
        """Windows PROVEN unreachable for this band. Empty is not a guarantee.

        The name is the contract: every window returned is impossible, and a
        result of () means nothing was proven impossible -- not that the design
        is feasible.
        """
        return tuple(w for w in self.windows if not w.admits(gc_min, gc_max))


@cache
def _bounds(codons: tuple[str, ...], lo: int, hi: int) -> tuple[int, int]:
    """(min, max) G+C count over `codons`, counting only offsets [lo, hi).

    Cached because a scan revisits the same (family, offset window) pair for
    every codon of that amino acid: only three distinct offset windows exist
    (whole codon, and the two partial ones at a window's edges).
    """
    counts = [sum(c[i] in "GC" for i in range(lo, hi)) for c in codons]
    return min(counts), max(counts)


def _forced_repeat(protein: str, families: dict[str, tuple[str, ...]]) -> int:
    """Longest DNA repeat the protein forces, in bp.

    Only single-codon residues can force one: everywhere else the two copies can
    be made to differ. So this looks for the longest run of single-codon
    residues that occurs at least twice.

    It is almost always 0, and that is the finding rather than a shortcoming.
    A (GGGGS)n linker forces NO repeat -- Gly has four codons and Ser six, so
    every copy can be encoded differently -- while forcing 60% GC. Repeats are
    recodeable; composition is not. Anything the solver must actually fight over
    repeats is a preference, not a floor, which is why `e6`/`e8` are SOFT and
    why no repeat budget belongs in this envelope.
    """
    stuck = {aa for aa, codons in families.items() if len(codons) == 1}
    if not stuck:
        return 0
    best = 0
    run_start = None
    for i, aa in enumerate(protein + "\0"):
        if aa in stuck:
            if run_start is None:
                run_start = i
            continue
        if run_start is not None:
            run = protein[run_start:i]
            # Occurs somewhere else too? Then both copies are forced identical.
            if protein.count(run) > 1:
                best = max(best, len(run))
            run_start = None
    return best * 3


def _culprits(
    protein: str, per_codon: list[tuple[int, int]], first: int, last: int, *, low: bool
) -> tuple[Culprit, ...]:
    """The residues holding a bound, worst first, merged into contiguous spans.

    For the LOW bound these are the residues with the highest floor -- the ones
    you would have to remove to go lower. For the HIGH bound, the lowest ceiling.
    """
    ranked = sorted(
        range(first, last + 1),
        key=lambda i: per_codon[i][0] if low else -per_codon[i][1],
        reverse=True,
    )
    chosen = sorted(ranked[:NAMED_RESIDUES])
    out: list[Culprit] = []
    for i in chosen:
        contribution = float(per_codon[i][0] if low else per_codon[i][1])
        if out and out[-1].end == i:
            prev = out[-1]
            out[-1] = Culprit(
                prev.start, i + 1, protein[prev.start : i + 1], prev.contribution + contribution
            )
        else:
            out.append(Culprit(i, i + 1, protein[i], contribution))
    return tuple(out)


def envelope(
    protein: str, code: NcbiGeneticCode, *, window: int = DEFAULT_WINDOW_BP
) -> FeasibleEnvelope:
    """The achievable GC envelope of `protein` under `code`.

    `code` is required and never defaulted, for the reason `TranslationUnit`
    gives: NCBI table 12 reassigns CTG to Ser, so the synonymous families -- and
    therefore the envelope -- differ between tables for the same protein.
    """
    if window <= 0:
        raise ValueError(f"window must be positive, got {window}")
    if not protein:
        raise ValueError("an empty protein has no envelope")

    families = dict(code.families())
    try:
        sets = [families[aa] for aa in protein]
    except KeyError as exc:
        raise ValueError(
            f"amino acid {exc.args[0]!r} has no non-stop codon in NCBI table "
            f"{code.table_id}, so this protein cannot be encoded under it at all"
        ) from None

    whole = [_bounds(s, 0, 3) for s in sets]
    length = 3 * len(protein)
    total_low = sum(b[0] for b in whole)
    total_high = sum(b[1] for b in whole)

    windows: list[GCWindow] = []
    for start in range(0, max(1, length - window + 1)):
        end = min(start + window, length)
        first, last = start // 3, (end - 1) // 3
        lo_count = hi_count = 0
        per_codon: list[tuple[int, int]] = [(0, 0)] * len(protein)
        for i in range(first, last + 1):
            a = max(start - 3 * i, 0)
            b = min(end - 3 * i, 3)
            pair = whole[i] if (a, b) == (0, 3) else _bounds(sets[i], a, b)
            per_codon[i] = pair
            lo_count += pair[0]
            hi_count += pair[1]
        width = end - start
        windows.append(
            GCWindow(
                start=start,
                end=end,
                lowest=100.0 * lo_count / width,
                highest=100.0 * hi_count / width,
                holds_low=_culprits(protein, per_codon, first, last, low=True),
                holds_high=_culprits(protein, per_codon, first, last, low=False),
            )
        )

    return FeasibleEnvelope(
        protein=protein,
        table_id=code.table_id,
        window=window,
        windows=tuple(windows),
        lowest=100.0 * total_low / length,
        highest=100.0 * total_high / length,
        forced_repeat_bp=_forced_repeat(protein, families),
    )
