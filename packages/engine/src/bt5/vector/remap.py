"""Coordinate remapping across an insert of a different length.

Back-translating a protein almost never yields a CDS the same length as the one
in the user's map, and the moment the length changes every downstream feature
coordinate in the vector is wrong. Silently shifted annotations are the failure
mode that makes an exported plasmid map look right and be wrong -- the promoter
label lands 60 bp into the CDS and nothing complains.

The policy here is deliberately conservative. A feature that OVERLAPS the replaced
span cannot be remapped faithfully, because part of what it described no longer
exists, so it is dropped and reported rather than clipped to something the source
file never said.
"""

from __future__ import annotations

from dataclasses import dataclass

from bt5.core.types import Interval


@dataclass(frozen=True, slots=True)
class IntervalRemapper:
    """Maps old-vector coordinates to assembled-construct coordinates.

    `replaced` must not wrap the origin; callers rotate the vector first, which is
    what `assemble()` does. That precondition is what keeps this arithmetic simple
    enough to be obviously correct.
    """

    replaced: Interval
    new_insert_length: int
    old_length: int
    circular: bool

    def __post_init__(self) -> None:
        if self.replaced.end > self.old_length:
            raise ValueError(
                f"the replaced span {self.replaced} wraps the origin; "
                f"rotate the vector before building a remapper"
            )

    @property
    def delta(self) -> int:
        """How much every coordinate after the insert moves."""
        return self.new_insert_length - self.replaced.length

    @property
    def new_length(self) -> int:
        return self.old_length + self.delta

    @property
    def insert_interval(self) -> Interval:
        """Where the new CDS lands, in assembled coordinates."""
        return Interval(
            self.replaced.start,
            self.replaced.start + self.new_insert_length,
            self.replaced.strand,
        )

    def position(self, p: int) -> int | None:
        """Remap a single base. None means the base was inside the replaced span."""
        if p < self.replaced.start:
            return p
        if p >= self.replaced.end:
            return p + self.delta
        return None

    def interval(self, iv: Interval) -> Interval | None:
        """Remap an interval, or None if it overlaps the replaced span.

        A wrapping interval shifts by `delta` at BOTH ends. Its head runs to the
        old length and its tail is expressed as `old_length + tail_end`; since the
        length itself grows by `delta`, the tail's position in the new frame is
        `new_length + tail_end`, which is the same arithmetic. A wrapping interval
        that does not overlap the insert must start after it, because its head
        already spans everything from `start` to the end of the sequence.
        """
        if iv.end <= self.old_length:
            if iv.end <= self.replaced.start:
                return iv
            if iv.start >= self.replaced.end:
                return Interval(iv.start + self.delta, iv.end + self.delta, iv.strand)
            return None

        if not self.circular:
            raise ValueError(f"interval {iv} runs past the end of a linear vector")

        tail_end = iv.end - self.old_length
        if iv.start >= self.replaced.end and tail_end <= self.replaced.start:
            return Interval(iv.start + self.delta, iv.end + self.delta, iv.strand)
        return None
