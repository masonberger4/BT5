"""GenBank locations <-> BT5 intervals.

Both directions are tested against each other, because a conversion that is
self-consistently wrong is exactly how an origin-spanning feature silently loses
its wrap.
"""

from __future__ import annotations

import pytest
from Bio.SeqFeature import BeforePosition, CompoundLocation, SimpleLocation
from bt5.core.types import Interval
from bt5.vector.locations import (
    LocationError,
    interval_to_location,
    location_to_interval,
    parts_to_location,
)


class TestLocationToInterval:
    def test_simple_forward(self) -> None:
        parsed = location_to_interval(SimpleLocation(10, 40, 1), length=100, circular=True)
        assert parsed.interval == Interval(10, 40, 1)
        assert not parsed.wraps_origin
        assert not parsed.is_compound

    def test_simple_reverse_keeps_strand(self) -> None:
        parsed = location_to_interval(SimpleLocation(10, 40, -1), length=100, circular=True)
        assert parsed.interval == Interval(10, 40, -1)

    def test_absent_strand_means_forward(self) -> None:
        parsed = location_to_interval(SimpleLocation(10, 40, None), length=100, circular=True)
        assert parsed.interval.strand == 1

    def test_origin_wrap_becomes_one_interval(self) -> None:
        loc = CompoundLocation([SimpleLocation(90, 100, 1), SimpleLocation(0, 20, 1)])
        parsed = location_to_interval(loc, length=100, circular=True)
        assert parsed.interval == Interval(90, 120, 1)
        assert parsed.wraps_origin
        assert parsed.interval.length == 30

    def test_origin_wrap_on_the_reverse_strand(self) -> None:
        loc = CompoundLocation([SimpleLocation(90, 100, -1), SimpleLocation(0, 20, -1)])
        parsed = location_to_interval(loc, length=100, circular=True)
        assert parsed.interval == Interval(90, 120, -1)
        assert parsed.wraps_origin

    def test_two_exon_join_is_not_a_wrap(self) -> None:
        """A join touching neither end is a spliced feature, not an origin wrap."""
        loc = CompoundLocation([SimpleLocation(10, 20, 1), SimpleLocation(30, 40, 1)])
        parsed = location_to_interval(loc, length=100, circular=True)
        assert not parsed.wraps_origin
        assert parsed.parts == (Interval(10, 20, 1), Interval(30, 40, 1))
        assert parsed.interval == Interval(10, 40, 1), "the single interval is the outer bound"

    def test_a_linear_record_never_wraps(self) -> None:
        loc = CompoundLocation([SimpleLocation(90, 100, 1), SimpleLocation(0, 20, 1)])
        parsed = location_to_interval(loc, length=100, circular=False)
        assert not parsed.wraps_origin
        assert parsed.is_compound

    def test_fuzzy_positions_are_refused(self) -> None:
        """`<1..500` asserts less than an exact boundary; coercing it would invent one."""
        loc = SimpleLocation(BeforePosition(0), 500, 1)
        with pytest.raises(LocationError, match="not an exact position"):
            location_to_interval(loc, length=1000, circular=True)


class TestIntervalToLocation:
    def test_simple(self) -> None:
        loc = interval_to_location(Interval(10, 40, 1), length=100)
        assert (int(loc.start), int(loc.end), loc.strand) == (10, 40, 1)

    def test_wrap_is_resplit_into_two_parts(self) -> None:
        loc = interval_to_location(Interval(90, 120, 1), length=100)
        assert [(int(p.start), int(p.end)) for p in loc.parts] == [(90, 100), (0, 20)]

    def test_rejects_a_start_past_the_end(self) -> None:
        with pytest.raises(LocationError, match="starts past"):
            interval_to_location(Interval(100, 130, 1), length=100)


class TestRoundTrip:
    @pytest.mark.parametrize(
        "interval",
        [
            Interval(0, 10, 1),
            Interval(10, 40, -1),
            Interval(90, 120, 1),
            Interval(90, 120, -1),
            Interval(99, 101, 1),
        ],
    )
    def test_interval_survives_both_directions(self, interval: Interval) -> None:
        loc = interval_to_location(interval, length=100)
        parsed = location_to_interval(loc, length=100, circular=True)
        assert parsed.interval == interval

    def test_compound_parts_survive(self) -> None:
        parts = (Interval(10, 20, 1), Interval(30, 40, 1))
        parsed = location_to_interval(parts_to_location(parts), length=100, circular=True)
        assert parsed.parts == parts
