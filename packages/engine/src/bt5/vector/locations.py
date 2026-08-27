"""GenBank locations and BT5 intervals: one conversion, used by everything.

BT5 has exactly ONE representation of an origin-spanning region -- an `Interval`
whose `end` exceeds the construct length. GenBank has a different one,
`join(1801..2000,1..80)`. Every bug where a feature silently loses its wrap, or a
plasmid map comes back with its origin somewhere else, is a bug in the translation
between those two representations. So that translation happens here and nowhere
else, and both directions are tested against each other.

A genuinely discontiguous location -- a spliced feature in the backbone, say -- is
NOT a wrap and cannot be squeezed into one `Interval`. Those keep their parts in
`ParsedLocation.parts` so export can rebuild them exactly; the single `interval`
is then the outer bound, used only for overlap questions.

Biopython ships `py.typed` but leaves its location constructors and the `Location`
base class unannotated. Rather than scatter suppressions through the lane, the
untyped boundary is confined to `_simple` and `_compound` below and to the
`GenBankLocation` protocol, so everything else in bt5.vector is strictly typed.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Protocol, cast, runtime_checkable

from Bio.SeqFeature import CompoundLocation, ExactPosition, Location, SimpleLocation

from bt5.core.types import Interval, Strand


class LocationError(ValueError):
    """A GenBank location BT5 cannot represent faithfully."""


@runtime_checkable
class GenBankSpan(Protocol):
    """One contiguous piece of a Biopython location."""

    @property
    def start(self) -> object: ...
    @property
    def end(self) -> object: ...


@runtime_checkable
class GenBankLocation(Protocol):
    """The read side of a Biopython `SimpleLocation` / `CompoundLocation`."""

    @property
    def parts(self) -> Sequence[GenBankSpan]: ...
    @property
    def strand(self) -> int | None: ...


def _simple(start: int, end: int, strand: Strand) -> Location:
    # Biopython's constructors carry no annotations; this is the only place the
    # lane calls them.
    return cast("Location", SimpleLocation(start, end, strand))  # type: ignore[no-untyped-call]


def _compound(parts: Sequence[Location]) -> Location:
    return cast("Location", CompoundLocation(list(parts)))  # type: ignore[no-untyped-call]


@dataclass(frozen=True, slots=True)
class ParsedLocation:
    """A GenBank location in BT5 coordinates.

    `parts` is empty for the ordinary case. It is populated only for a genuinely
    multi-part location that is not an origin wrap, precisely so that exporting it
    does not quietly flatten a two-exon feature into one span.
    """

    interval: Interval
    parts: tuple[Interval, ...] = ()
    wraps_origin: bool = False

    @property
    def is_compound(self) -> bool:
        return bool(self.parts)


def _strand_of(location: GenBankLocation) -> Strand:
    """A location's strand, defaulting to forward.

    Biopython reports `None` when a location carries no explicit strand, which in
    GenBank means the forward strand.
    """
    return -1 if location.strand == -1 else 1


def _exact(value: object, *, what: str) -> int:
    """Reject fuzzy positions rather than silently coercing them.

    `<1..500` and `100.200` are real GenBank; coercing either to an exact integer
    would hand the optimizer a boundary the source file never asserted.
    """
    if not isinstance(value, ExactPosition):
        raise LocationError(f"{what} is not an exact position: {value!r}")
    return int(value)


def location_to_interval(
    location: GenBankLocation, *, length: int, circular: bool
) -> ParsedLocation:
    """Convert a Biopython location into BT5 coordinates.

    An origin wrap is recognised ONLY under conditions that cannot be confused
    with a two-exon feature: the construct is circular, there are exactly two
    parts, one starts at 0 and the other ends at `length`.
    """
    strand = _strand_of(location)
    # BIOLOGICAL order, as Biopython supplies it: 5'->3' along the feature's own
    # strand, which for a minus-strand feature is DESCENDING in genomic
    # coordinates. That order is load-bearing and must not be normalised away.
    # `complement(join(A,B))` concatenates A then B and reverse-complements the
    # result, so re-emitting the parts sorted ascending describes a different
    # sequence -- an AmpR written that way translates from its signal peptide
    # rather than from its start codon, silently, in the user's exported map.
    ordered = [
        (_exact(part.start, what="location start"), _exact(part.end, what="location end"))
        for part in location.parts
    ]
    if not ordered:
        raise LocationError("location has no parts")
    # Sorted only to recognise an origin wrap and to bound the outer interval.
    spans = sorted(ordered)

    if circular and len(spans) == 2 and spans[0][0] == 0 and spans[1][1] == length:
        # The tail [0, e) is the continuation of the head [s, length).
        return ParsedLocation(
            interval=Interval(spans[1][0], length + spans[0][1], strand),
            wraps_origin=True,
        )

    if len(ordered) == 1:
        start, end = ordered[0]
        return ParsedLocation(interval=Interval(start, end, strand))

    parts = tuple(Interval(s, e, strand) for s, e in ordered)
    return ParsedLocation(interval=Interval(spans[0][0], spans[-1][1], strand), parts=parts)


def interval_to_location(iv: Interval, *, length: int) -> Location:
    """The inverse. A wrapping interval is re-split into a two-part join.

    Part order is genomic (head first, then the wrapped tail) for BOTH strands.
    Biopython writes the minus-strand case as `complement(join(1..80,1801..2000))`
    and re-parses it back to this same part order, so the round trip is exact.
    """
    if iv.end <= length:
        return _simple(iv.start, iv.end, iv.strand)
    if iv.start >= length:
        raise LocationError(f"interval {iv} starts past the construct length {length}")
    return _compound([_simple(iv.start, length, iv.strand), _simple(0, iv.end - length, iv.strand)])


def parts_to_location(parts: Sequence[Interval]) -> Location:
    """Rebuild a genuinely multi-part location from its preserved parts."""
    if not parts:
        raise LocationError("cannot build a location from zero parts")
    if len(parts) == 1:
        return _simple(parts[0].start, parts[0].end, parts[0].strand)
    return _compound([_simple(p.start, p.end, p.strand) for p in parts])
