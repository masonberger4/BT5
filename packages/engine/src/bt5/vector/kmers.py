"""Exact-repeat queries over the assembled construct.

BIOSECURITY. `ConstructKmerIndex.of()` takes a `Construct` and nothing else.
There is deliberately no constructor accepting an external sequence database:
pointing a homology-minimiser at an arbitrary target turns BT5 into a
general-purpose screening-evasion tool, and constraining the index to the
assembled construct is what keeps it from becoming one. A CI grep enforces this.

Why repeats get a two-dimensional risk surface rather than a length cutoff. Two
thresholds are routinely conflated. RecBCD MEPS (23-27 bp) is the floor for
RecA-DEPENDENT recombination, but below roughly 200 bp deletion proceeds by a
RecA-INDEPENDENT route -- slipped-strand mispairing and single-strand annealing --
which a recA- strain does not suppress at all. That route is strongly
proximity-sensitive: inserting sequence between two copies suppresses it. So risk
is a function of (length, spacer), and the most dangerous configuration is a
TANDEM repeat, where the copies touch and mispairing needs no looping at all.

  Springer   https://link.springer.com/article/10.1007/BF00290109
  PMC5426353 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5426353/
  PNAS       https://www.pnas.org/doi/10.1073/pnas.111008398
  Oliveira 2008, a 28 bp pair still recombining at 7.8e-7 to 3.1e-5 in FOUR
  different recA- strains:
  https://www.genoscope.cns.fr/MGE/pubs/Oliveira_Mol_Biotechnol_2008.pdf

This module reports. It never claims a rate, and it never tells a user that a
recA- strain covers a repeat in the short regime, because it does not.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

from bt5.core.types import Construct, Interval, reverse_complement

if TYPE_CHECKING:
    from bt5.core.services import KmerIndex

RiskBand = Literal["low", "moderate", "high"]

#: Below this, a repeat is under the shortest vendor uniqueness requirement and
#: under the RecBCD MEPS floor. Twist and Gibson protocols both work to ~20 bp.
VENDOR_UNIQUENESS_BP = 20

#: Above this, recombination becomes RecA-DEPENDENT, which is what a recA-
#: strain actually suppresses. Below it the strain gives no protection.
RECA_DEPENDENT_BP = 200

#: Edge-to-edge separations at which the RecA-independent route stays efficient.
#: Proximity is the variable, not an afterthought.
NEAR_SPACER_BP = 100
FAR_SPACER_BP = 1000

#: A repeat at least this long is never reported as low risk, however far
#: apart the copies are.
SUBSTANTIAL_BP = 100

#: A cap so a tandem array cannot flood the report with thousands of pairs.
MAX_PAIRS = 200


@dataclass(frozen=True, slots=True)
class RepeatPair:
    """Two exact copies of the same sequence, with the geometry that matters."""

    first: Interval
    second: Interval
    length: int
    spacer: int
    #: True when the copies touch or overlap: a tandem array, the worst case.
    tandem: bool

    @property
    def risk(self) -> RiskBand:
        return repeat_risk(self.length, self.spacer, tandem=self.tandem)

    @property
    def reca_strain_helps(self) -> bool:
        """Only long repeats are suppressed by a recA- strain."""
        return self.length >= RECA_DEPENDENT_BP


def repeat_risk(length: int, spacer: int, *, tandem: bool = False) -> RiskBand:
    """Classify a repeat on the (length, spacer) surface.

    Deliberately coarse. The literature supports the ORDERING of these regimes
    and the claim that proximity matters; it does not support a calibrated rate,
    and presenting three bands is the honest resolution.
    """
    if tandem:
        # Checked BEFORE the length floor. Slipped-strand mispairing needs no
        # loop at all, so a short tandem array is a genuine liability even below
        # the vendor uniqueness threshold -- repetitive 9-mers per 100 bp is one
        # of the two highest-importance features in the published synthesis
        # success model (https://pubs.acs.org/doi/10.1021/acssynbio.9b00460).
        return "high" if length >= VENDOR_UNIQUENESS_BP else "moderate"
    if length < VENDOR_UNIQUENESS_BP:
        return "low"
    if length >= RECA_DEPENDENT_BP:
        # Long: real either way. This is the regime a recA- strain actually covers.
        return "moderate" if spacer > FAR_SPACER_BP else "high"
    # The RecA-INDEPENDENT regime, where the strain gives nothing.
    if spacer <= NEAR_SPACER_BP:
        return "high"
    if spacer <= FAR_SPACER_BP:
        return "moderate"
    # Distance suppresses the RecA-independent route, but a substantial exact
    # repeat is never "low" just because the copies are far apart -- a 189 bp
    # identity between two lentiviral LTRs is a real liability at any spacing.
    return "moderate" if length >= SUBSTANTIAL_BP else "low"


@dataclass(frozen=True, slots=True)
class InvertedRepeat:
    """Two arms that base-pair with each other: a hairpin or cruciform precursor.

    Mechanically distinct from a direct repeat and not interchangeable with one.
    A direct repeat is lost by deletion -- slipped-strand mispairing or
    single-strand annealing -- and codon choice fixes it. An inverted repeat
    extrudes a cruciform, stalls replication forks and is cleaved by SbcCD; the
    answer is usually a strain and a temperature, not a redesign. That is why it
    gets its own type and does NOT feed `repeat_risk`, whose bands are calibrated
    on the deletion literature and would be a category error here.

    Banding inverted repeats needs its own citations. Until the rules lane does
    that pass, this reports geometry and nothing else.
    """

    #: The 5' arm, read on the plus strand.
    first: Interval
    #: The 3' arm, reported on the MINUS strand, where it reads identically to
    #: the first. "Opposite direction on one strand" and "same direction on
    #: opposite strands" describe one physical object; strand -1 records which
    #: of the two readings this interval is.
    second: Interval
    stem: int
    loop: int

    @property
    def perfect_palindrome(self) -> bool:
        """Arms abutting with no loop: the most extrudable configuration."""
        return self.loop == 0


class ConstructKmerIndex:
    """Exact repeats and inverted repeats within one assembled construct."""

    def __init__(self, construct: Construct, k: int) -> None:
        self._construct = construct
        self._k = k
        self._n = construct.length
        # A circular construct is scanned on the doubled sequence so a repeat
        # spanning the origin is found like any other; hits are folded back.
        self._text = construct.sequence * 2 if construct.is_circular else construct.sequence

    @classmethod
    def of(cls, c: Construct, k: int) -> ConstructKmerIndex:
        """The ONLY constructor. See the biosecurity note in the module docstring."""
        if k < 1:
            raise ValueError(f"k must be positive, got {k}")
        return cls(c, k)

    @property
    def k(self) -> int:
        return self._k

    # -- the two protocol methods ------------------------------------------
    #
    # `KmerIndex` in bt5.core is what a RULE sees, and it is deliberately
    # narrow: two intervals per finding and nothing else. The rich forms below
    # (`repeat_pairs`, `inverted_repeats`) carry this lane's own value types,
    # which the frozen contract does not know about and should not have to.
    #
    # Both adapters are thin on purpose. The geometry a rule needs is
    # recoverable from the pair -- a stem is the first arm's length, a loop is
    # the gap between the arms -- so narrowing here costs a caller nothing and
    # keeps `InvertedRepeat` out of a contract that would then have to freeze it.

    def duplicates(self, min_len: int) -> Iterator[tuple[Interval, Interval]]:
        """Direct repeat pairs at least `min_len` long, satisfying the protocol."""
        for pair in self.repeat_pairs(min_len):
            yield pair.first, pair.second

    def revcomp_pairs(self, min_stem: int, max_loop: int) -> Iterator[tuple[Interval, Interval]]:
        """Inverted repeat arm pairs, satisfying the protocol.

        This name previously belonged to the rich form and returned
        `list[InvertedRepeat]`, so `ConstructKmerIndex` did not in fact satisfy
        `KmerIndex` -- a rule reaching through `Services.kmer` would have been
        handed value objects where the contract promised tuples. Nothing caught
        it because nothing had yet consumed the protocol.
        """
        for repeat in self.inverted_repeats(min_stem, max_loop):
            yield repeat.first, repeat.second

    def repeat_pairs(self, min_len: int, *, exclude: Sequence[Interval] = ()) -> list[RepeatPair]:
        """Maximal exact direct repeats, longest first, with their geometry.

        `exclude` drops pairs whose BOTH copies sit inside a listed region --
        used for ITRs and LTRs, which are reported separately as an accepted
        design feature rather than as a finding.
        """
        seeds: dict[str, list[int]] = {}
        limit = self._n if self._construct.is_circular else max(0, self._n - min_len + 1)
        for i in range(limit):
            kmer = self._text[i : i + min_len]
            if len(kmer) == min_len:
                seeds.setdefault(kmer, []).append(i)

        seen: set[tuple[int, int]] = set()
        out: list[RepeatPair] = []
        # Every seed of one match shares the diagonal `b - a`, which neither
        # extension direction changes. Remembering the span already grown on
        # each diagonal turns the redundant seeds into a lookup instead of a
        # re-walk -- the same device `inverted_repeats` uses below, and for the
        # same reason. Without it a 6 kb tandem array walks the full text once
        # per seed: 5,974 calls and 35.7M character comparisons, 1.9 seconds, to
        # return a single pair. See `_grow`.
        grown: dict[int, tuple[int, int]] = {}
        for positions in seeds.values():
            if len(positions) < 2:
                continue
            for a, b in zip(positions, positions[1:], strict=False):
                diagonal = b - a
                span = grown.get(diagonal)
                if span is not None and span[0] <= a and a + min_len <= span[1]:
                    continue  # already grown, from a seed of the same match
                start_a, end_a = self._grow(a, b, min_len)
                grown[diagonal] = (start_a, end_a)
                pair = self._pair(start_a, end_a, diagonal, min_len)
                if pair is None:
                    continue
                key = (pair.first.start, pair.second.start)
                if key in seen:
                    continue
                seen.add(key)
                if _inside_any(pair.first, exclude) and _inside_any(pair.second, exclude):
                    continue
                out.append(pair)
        out.sort(key=lambda p: (-p.length, p.first.start))
        return _drop_contained(out)[:MAX_PAIRS]

    def _grow(self, a: int, b: int, min_len: int) -> tuple[int, int]:
        """The maximal match around a seed, as (start, end) of the FIRST copy.

        The second copy is recoverable as `start + (b - a)`: both walks move the
        two positions together, so the diagonal is invariant. That is what lets
        the caller memoise by diagonal, and it is why this is split out from
        `_pair` -- the raw extent has to be recorded even when the pair is
        rejected, or a tandem array whose period is under `min_len` re-walks the
        whole text for every one of its thousands of seeds.

        Both walks are unbounded on purpose. Bounding the rightward one by the
        period looks safe, since `_pair` clamps a tandem match back to its period
        anyway, but it is only half the cost -- measured on a 6 kb tandem array
        the two directions are symmetric at 17.8M steps each -- and it would
        shorten the recorded span and so defeat the memo that fixes the other
        half. The leftward walk cannot be bounded at all without moving the
        reported start, which is what deduplicates the array down to one finding.
        """
        text = self._text
        start_a, start_b = a, b
        end_a, end_b = a + min_len, b + min_len
        while start_a > 0 and start_b > start_a and text[start_a - 1] == text[start_b - 1]:
            start_a -= 1
            start_b -= 1
        while end_b < len(text) and text[end_a] == text[end_b]:
            end_a += 1
            end_b += 1
        return start_a, end_a

    def _pair(self, start_a: int, end_a: int, diagonal: int, min_len: int) -> RepeatPair | None:
        """Turn a grown extent into a reportable pair, or reject it.

        On a circular construct the scan runs over the doubled sequence, so every
        position trivially matches itself one full length away. Those pairs are
        artefacts of the doubling, not repeats, and are dropped -- without that
        check every k-mer in the plasmid reports as a repeat.
        """
        n = self._n
        start_b = start_a + diagonal
        length = end_a - start_a
        if length < min_len:
            return None
        if start_a >= n:
            return None  # both copies live in the doubled tail
        if self._construct.is_circular and diagonal == n:
            return None  # the periodicity artefact
        if diagonal < length:
            # Overlapping copies: a tandem array. Report the period, not the
            # smeared extension, so the geometry stays interpretable.
            length = diagonal
            if length < min_len:
                return None
        return RepeatPair(
            first=Interval(start_a, start_a + length),
            second=Interval(start_b, start_b + length),
            length=length,
            spacer=start_b - (start_a + length),
            tandem=start_b - (start_a + length) <= 0,
        )

    def inverted_repeats(
        self, min_stem: int, max_loop: int, *, exclude: Sequence[Interval] = ()
    ) -> list[InvertedRepeat]:
        """Maximal inverted repeats: a stem of at least `min_stem` within `max_loop`.

        Maximal, like the direct-repeat scan and for the same reason: seeding at
        a fixed width and reporting the seed calls a 60 bp stem a 20 bp one, once
        per offset.

        The candidate window is what keeps that affordable. Indexing every k-mer
        and pairing every occurrence is quadratic, and it degrades exactly where
        a plasmid is most likely to need the answer: 800 bp of alternating AT
        produced 464,799 raw pairs, of which 1,180 survived, and took 26 seconds.
        Instead, for each 5' arm only the `max_loop` positions that could hold
        the matching 3' arm are tested -- O(n x max_loop), and 0.02s on the same
        input. The window is sound because a stem's INNERMOST seed already has
        the final loop, so it always falls inside; growing outward from there
        recovers the rest. That is what the outward pass in `_extend_stem` buys,
        and why it is not the redundant half it looks like.
        """
        text = self._text
        n, end = self._n, len(self._text) - min_stem + 1
        seen: set[tuple[int, int]] = set()
        # Every seed of one stem shares the value `a_start + b_end`, which both
        # extension directions leave untouched. Remembering the span already
        # grown on each such diagonal turns the redundant seeds into a lookup
        # instead of a re-walk: without it, 4 kb of alternating AT re-extends
        # across the whole run for every candidate and takes 27 seconds.
        grown: dict[int, tuple[int, int]] = {}
        out: list[InvertedRepeat] = []
        probe_limit = n if self._construct.is_circular else max(0, n - min_stem + 1)
        for i in range(probe_limit):
            probe = reverse_complement(text[i : i + min_stem])
            lo = i + min_stem  # closer than this and the arms would overlap
            for j in range(lo, min(lo + max_loop + 1, end)):
                if text[j : j + min_stem] != probe:
                    continue
                span = grown.get(i + j)
                if span is not None and span[0] <= i and j + min_stem <= span[1]:
                    continue  # already grown, from a seed of the same stem
                stem = self._extend_stem(i, j, min_stem, max_loop)
                if stem is None:
                    continue
                grown[i + j] = (stem.first.start, stem.second.end)
                key = (stem.first.start, stem.second.end)
                if key in seen:
                    continue
                seen.add(key)
                if _inside_any(stem.first, exclude) and _inside_any(stem.second, exclude):
                    continue
                out.append(stem)
        out.sort(key=lambda p: (-p.stem, p.first.start))
        return _drop_contained_stems(out)[:MAX_PAIRS]

    def _extend_stem(self, i: int, j: int, min_stem: int, max_loop: int) -> InvertedRepeat | None:
        """Grow a seeded stem to its full length, outward and then inward.

        A stem grows differently from a direct repeat. Its arms pair END TO END,
        so the 5' arm extends LEFT exactly when the 3' arm extends RIGHT, and
        they can also close on each other, eating the loop two bases at a time.
        Sliding both in the same direction, as the direct-repeat extension does,
        would compare the wrong two bases.

        Both directions are load-bearing given the windowed search above. The
        seed that falls inside the window is the innermost one, whose loop is
        already final; outward is the only way from there to the full stem.
        """
        text, n = self._text, self._n
        a_start, a_end = i, i + min_stem
        b_start, b_end = j, j + min_stem
        while a_start > 0 and b_end < len(text) and _complementary(text[a_start - 1], text[b_end]):
            a_start -= 1
            b_end += 1
        while b_start - a_end >= 2 and _complementary(text[a_end], text[b_start - 1]):
            a_end += 1
            b_start -= 1

        stem, loop = a_end - a_start, b_start - a_end
        if stem < min_stem or loop > max_loop:
            return None
        if a_start >= n:
            return None  # both arms live in the doubled tail
        if self._construct.is_circular and b_start - a_start == n:
            return None  # a self-complementary arm matching its own second copy
        return InvertedRepeat(
            first=Interval(a_start, a_end),
            second=Interval(b_start, b_end, -1),
            stem=stem,
            loop=loop,
        )


_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


def _complementary(x: str, y: str) -> bool:
    return _COMPLEMENT.get(x) == y


def _drop_contained_stems(
    stems: list[InvertedRepeat], limit: int = MAX_PAIRS
) -> list[InvertedRepeat]:
    """Keep only maximal stems: one hairpin otherwise reports once per seed offset.

    Stops at `limit`, which is also all the caller returns. The comparison is
    against everything kept so far, so an unbounded pass is quadratic in the
    number of findings -- 3.3 of the 9 seconds a 10 kb pure-AT sequence used to
    cost. Input arrives sorted longest-first, so the stems that survive the cap
    are the ones that would have survived anyway.
    """
    kept: list[InvertedRepeat] = []
    for stem in stems:
        if len(kept) >= limit:
            break
        if any(
            k.first.start <= stem.first.start
            and stem.first.end <= k.first.end
            and k.second.start <= stem.second.start
            and stem.second.end <= k.second.end
            for k in kept
        ):
            continue
        kept.append(stem)
    return kept


def _inside_any(iv: Interval, regions: Sequence[Interval]) -> bool:
    return any(r.start <= iv.start and iv.end <= r.end for r in regions)


def _drop_contained(pairs: list[RepeatPair]) -> list[RepeatPair]:
    """Keep only maximal pairs: a tandem array otherwise reports every offset."""
    kept: list[RepeatPair] = []
    for pair in pairs:
        if any(
            k.first.start <= pair.first.start
            and pair.first.end <= k.first.end
            and k.second.start <= pair.second.start
            and pair.second.end <= k.second.end
            for k in kept
        ):
            continue
        kept.append(pair)
    return kept


if TYPE_CHECKING:
    # `ConstructKmerIndex` is the only implementation of the frozen `KmerIndex`
    # protocol, and for a while it did not satisfy it: `revcomp_pairs` named the
    # rich form and returned `list[InvertedRepeat]` where the contract promises
    # `Iterator[tuple[Interval, Interval]]`. Nothing caught it because nothing
    # had yet consumed the protocol.
    #
    # This assertion lives in src/ rather than in a test on purpose. mypy is
    # configured over `packages/engine/src/bt5` only, so the same line written
    # in a test file type-checks nowhere and proves nothing -- which is exactly
    # the shape of the gap it exists to close.
    _protocol_conformance: type[KmerIndex] = ConstructKmerIndex
