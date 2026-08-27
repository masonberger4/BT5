"""Coordinate remapping across an insert of a different length.

The failure this guards against is silent: a downstream feature that keeps its
old coordinates after the CDS changed length still exports as a valid plasmid
map, with the promoter label sitting 60 bp inside the coding sequence.
"""

from __future__ import annotations

import pytest
from bt5.core.types import Interval
from bt5.vector.remap import IntervalRemapper


def remapper(*, old_cds: Interval, new_len: int, total: int = 1000) -> IntervalRemapper:
    return IntervalRemapper(
        replaced=old_cds, new_insert_length=new_len, old_length=total, circular=True
    )


class TestDelta:
    def test_longer_insert_shifts_forward(self) -> None:
        r = remapper(old_cds=Interval(400, 700), new_len=360)
        assert r.delta == 60
        assert r.new_length == 1060
        assert r.insert_interval == Interval(400, 760)

    def test_shorter_insert_shifts_back(self) -> None:
        r = remapper(old_cds=Interval(400, 700), new_len=240)
        assert r.delta == -60
        assert r.new_length == 940

    def test_same_length_is_identity(self) -> None:
        r = remapper(old_cds=Interval(400, 700), new_len=300)
        assert r.delta == 0
        assert r.interval(Interval(800, 900)) == Interval(800, 900)


class TestPositions:
    def test_before_is_untouched(self) -> None:
        r = remapper(old_cds=Interval(400, 700), new_len=360)
        assert r.position(0) == 0
        assert r.position(399) == 399

    def test_after_shifts(self) -> None:
        r = remapper(old_cds=Interval(400, 700), new_len=360)
        assert r.position(700) == 760
        assert r.position(999) == 1059

    def test_inside_the_replaced_span_is_gone(self) -> None:
        r = remapper(old_cds=Interval(400, 700), new_len=360)
        assert r.position(400) is None
        assert r.position(699) is None


class TestIntervals:
    def test_wholly_before(self) -> None:
        r = remapper(old_cds=Interval(400, 700), new_len=360)
        assert r.interval(Interval(100, 400)) == Interval(100, 400), "abutting is not overlapping"

    def test_wholly_after(self) -> None:
        r = remapper(old_cds=Interval(400, 700), new_len=360)
        assert r.interval(Interval(700, 800)) == Interval(760, 860)

    def test_strand_is_preserved(self) -> None:
        r = remapper(old_cds=Interval(400, 700), new_len=360)
        assert r.interval(Interval(700, 800, -1)) == Interval(760, 860, -1)

    @pytest.mark.parametrize(
        "overlapping",
        [Interval(350, 450), Interval(650, 750), Interval(450, 500), Interval(300, 900)],
    )
    def test_overlapping_is_dropped_not_clipped(self, overlapping: Interval) -> None:
        """Clipping would assert a boundary the source file never contained."""
        r = remapper(old_cds=Interval(400, 700), new_len=360)
        assert r.interval(overlapping) is None

    def test_wrapping_interval_shifts_at_both_ends(self) -> None:
        """An origin-spanning feature keeps its span; both endpoints move by delta."""
        r = remapper(old_cds=Interval(400, 700), new_len=360)
        moved = r.interval(Interval(900, 1100))
        assert moved == Interval(960, 1160)
        assert moved is not None
        assert moved.length == 200, "the feature must not change length"
        assert moved.end - r.new_length == 100, "its tail still ends 100 bases past the origin"

    def test_wrapping_interval_overlapping_the_insert_is_dropped(self) -> None:
        r = remapper(old_cds=Interval(400, 700), new_len=360)
        assert r.interval(Interval(600, 1100)) is None

    def test_wrapping_tail_reaching_the_insert_is_dropped(self) -> None:
        r = remapper(old_cds=Interval(400, 700), new_len=360)
        assert r.interval(Interval(900, 1450)) is None


class TestPreconditions:
    def test_a_wrapping_replaced_span_is_refused(self) -> None:
        """assemble() rotates first; this is the guard that says so."""
        with pytest.raises(ValueError, match="wraps the origin"):
            IntervalRemapper(
                replaced=Interval(900, 1100),
                new_insert_length=100,
                old_length=1000,
                circular=True,
            )

    def test_a_wrapping_interval_on_a_linear_vector_is_an_error(self) -> None:
        r = IntervalRemapper(
            replaced=Interval(400, 700), new_insert_length=360, old_length=1000, circular=False
        )
        with pytest.raises(ValueError, match="past the end of a linear"):
            r.interval(Interval(900, 1100))
