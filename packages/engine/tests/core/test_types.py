"""Unit tests for the construct geometry.

These are lane-owned (they live under packages/engine/tests) as distinct from the
owner-gated invariants in /tests, which are the specification and cannot be
modified without a label.
"""

from __future__ import annotations

import pytest
from bt5.core.types import (
    Construct,
    Interval,
    Segment,
    SegmentKind,
    Topology,
    TranslationUnit,
    reverse_complement,
)


class TestInterval:
    def test_length(self) -> None:
        assert Interval(4, 10).length == 6

    def test_rejects_empty_and_negative(self) -> None:
        with pytest.raises(ValueError, match="non-empty"):
            Interval(5, 5)
        with pytest.raises(ValueError, match="non-empty"):
            Interval(9, 4)
        with pytest.raises(ValueError, match="start must be"):
            Interval(-1, 4)

    def test_wraps_only_past_the_end(self) -> None:
        assert not Interval(0, 10).wraps(10)
        assert Interval(8, 14).wraps(10)

    def test_overlaps_is_half_open(self) -> None:
        assert Interval(0, 5).overlaps(Interval(4, 9), 100, False)
        assert not Interval(0, 5).overlaps(Interval(5, 9), 100, False), (
            "half-open ranges must not touch"
        )

    def test_overlaps_finds_a_hit_across_the_origin(self) -> None:
        """The case the linear predicate got wrong, in both argument orders.

        A feature stored as [4900, 5100) on a 5000 bp plasmid and a hit at
        [10, 40) occupy the same bases. Both are already in BT5's one canonical
        representation, so nothing a caller could do before the comparison makes
        plain [start, end) arithmetic say so.
        """
        wrapping = Interval(4900, 5100)
        plain = Interval(10, 40)
        assert wrapping.overlaps(plain, 5000, True)
        assert plain.overlaps(wrapping, 5000, True), "the wrap may be on either side"
        assert not wrapping.overlaps(plain, 5000, False), (
            "a linear construct has no origin to span, so the same pair must not overlap"
        )

    def test_overlaps_does_not_invent_a_hit_a_full_turn_away(self) -> None:
        wrapping = Interval(4900, 5100)
        assert not wrapping.overlaps(Interval(200, 300), 5000, True)

    def test_contains_is_wrap_aware_and_requires_the_whole_interval(self) -> None:
        outer = Interval(4900, 5100)
        assert outer.contains(Interval(10, 40), 5000, True)
        assert outer.contains(Interval(4950, 4980), 5000, True)
        assert not outer.contains(Interval(90, 110), 5000, True), (
            "the tail past 5100 is outside, so this is not wholly contained"
        )
        assert not outer.contains(Interval(10, 40), 5000, False)

    def test_the_predicates_ignore_strand(self) -> None:
        """They answer 'which bases'. A rule that also cares compares .strand
        itself -- folding it in here would hide a reverse-strand hit sitting
        inside a forward-strand feature."""
        fwd = Interval(10, 40, 1)
        rev = Interval(20, 30, -1)
        assert fwd.overlaps(rev, 100, False)
        assert fwd.contains(rev, 100, False)

    def test_extended_clamps_at_the_ends_when_linear(self) -> None:
        assert Interval(2, 6).extended(5, 20, circular=False) == Interval(0, 11)
        assert Interval(14, 18).extended(5, 20, circular=False) == Interval(9, 20)

    def test_extended_stays_inside_the_sequence_when_circular(self) -> None:
        assert Interval(10, 14).extended(5, 20, circular=True) == Interval(5, 19)

    def test_extended_wraps_rather_than_going_negative_at_the_origin(self) -> None:
        """A breach near the origin of a circular plasmid must still produce a
        usable repair window; a negative start would crash Interval."""
        original = Interval(2, 6)  # length 4
        wide = original.extended(5, 20, circular=True)
        assert wide.start == 17, "start wraps to the far end rather than going negative"
        assert wide.end == 31
        assert wide.wraps(20)
        assert wide.length == original.length + 2 * 5 == 14
        # The window covers positions 17,18,19 then 0..10, so it still contains
        # every base of the original breach.
        covered = {i % 20 for i in range(wide.start, wide.end)}
        assert set(range(original.start, original.end)) <= covered


class TestSegment:
    def test_only_designable_cds_is_editable(self) -> None:
        iv = Interval(0, 9)
        assert Segment(iv, SegmentKind.DESIGNABLE_CDS).is_editable
        assert not Segment(iv, SegmentKind.BACKBONE).is_editable
        assert not Segment(iv, SegmentKind.ANNOTATED_INTRON).is_editable

    def test_introns_and_whitelisted_repeats_are_scan_exempt(self) -> None:
        """A deliberately placed intron must survive the splice remover, and LTRs
        and ITRs violate the repeat rules by construction."""
        iv = Interval(0, 9)
        assert Segment(iv, SegmentKind.ANNOTATED_INTRON).exempt_from_scanning
        assert Segment(iv, SegmentKind.WHITELISTED_REPEAT).exempt_from_scanning
        assert not Segment(iv, SegmentKind.BACKBONE).exempt_from_scanning
        assert not Segment(iv, SegmentKind.DESIGNABLE_CDS).exempt_from_scanning


class TestConstruct:
    def build(self, seq: str = "ATGAAACCCTAAGGGCCC", circular: bool = True) -> Construct:
        return Construct(
            sequence=seq,
            topology=Topology.CIRCULAR if circular else Topology.LINEAR,
            segments=(
                Segment(Interval(0, 12), SegmentKind.DESIGNABLE_CDS, "cds"),
                Segment(Interval(12, len(seq)), SegmentKind.BACKBONE, "vector"),
            ),
        )

    def test_rejects_non_acgt(self) -> None:
        with pytest.raises(ValueError, match="non-ACGT"):
            Construct("ATGN", Topology.LINEAR, ())

    def test_slice_forward(self) -> None:
        assert self.build().slice(Interval(0, 3)) == "ATG"

    def test_slice_reverse_strand(self) -> None:
        assert self.build().slice(Interval(0, 3, -1)) == reverse_complement("ATG") == "CAT"

    def test_slice_wraps_the_origin_when_circular(self) -> None:
        c = self.build()
        # last three bases + first three bases
        assert c.slice(Interval(15, 21)) == "CCC" + "ATG"

    def test_slice_past_the_end_is_an_error_when_linear(self) -> None:
        with pytest.raises(ValueError, match="runs past the end"):
            self.build(circular=False).slice(Interval(15, 21))

    def test_tripled_offset(self) -> None:
        c = self.build()
        tripled, offset = c.tripled()
        assert offset == c.length
        assert len(tripled) == 3 * c.length
        assert tripled[offset : offset + c.length] == c.sequence

    def test_tripled_is_a_noop_when_linear(self) -> None:
        c = self.build(circular=False)
        assert c.tripled() == (c.sequence, 0)

    def test_editable_is_the_cds_and_its_complement_is_backbone(self) -> None:
        c = self.build()
        assert c.editable == (Interval(0, 12),)
        assert c.is_editable(Interval(3, 6))
        assert not c.is_editable(Interval(12, 15)), "backbone must never be editable"
        assert not c.is_editable(Interval(9, 15)), (
            "an interval straddling the junction is not editable"
        )

    def test_a_cds_spanning_the_origin_is_still_editable(self) -> None:
        """The insert is cloned across the origin, so the CDS segment wraps.

        Plain [start, end) comparison calls every codon past the origin backbone
        and hands the solver a mutation space with a hole in it -- which reads
        downstream as a design that cannot satisfy its own constraints rather
        than as a coordinate bug.
        """
        seq = "ATGAAACCCTAAGGGCCC"  # length 18
        c = Construct(
            sequence=seq,
            topology=Topology.CIRCULAR,
            segments=(
                # 15..18 then 0..6 -- six bases of CDS, straddling the origin.
                Segment(Interval(15, 24), SegmentKind.DESIGNABLE_CDS, "cds"),
                Segment(Interval(6, 15), SegmentKind.BACKBONE, "vector"),
            ),
        )
        assert c.is_editable(Interval(16, 18)), "before the origin"
        assert c.is_editable(Interval(0, 3)), "after the origin, the case that regressed"
        assert c.is_editable(Interval(15, 24)), "the segment contains itself"
        assert not c.is_editable(Interval(7, 10)), "the backbone is still off limits"

    def test_overlaps_editable_is_the_fixable_by_codon_choice_question(self) -> None:
        c = self.build()
        assert c.overlaps_editable(Interval(3, 6))
        assert not c.overlaps_editable(Interval(12, 15)), "wholly in the backbone"
        assert c.overlaps_editable(Interval(9, 15)), (
            "a hit straddling the junction IS partly recodable, which is exactly "
            "where is_editable gives the wrong answer for fixability"
        )

    def test_exempt_collects_scan_exempt_segments(self) -> None:
        seq = "ATGAAACCCTAAGGGCCC"
        c = Construct(
            seq,
            Topology.CIRCULAR,
            (
                Segment(Interval(0, 12), SegmentKind.DESIGNABLE_CDS, "cds"),
                Segment(Interval(12, 18), SegmentKind.WHITELISTED_REPEAT, "LTR"),
            ),
        )
        assert c.exempt == (Interval(12, 18),)


class TestTranslationUnit:
    def test_table_id_must_be_positive(self) -> None:
        with pytest.raises(ValueError, match="positive NCBI table id"):
            TranslationUnit(0, (Interval(0, 3),), "M")

    def test_defaults_to_requiring_an_initiator(self) -> None:
        """A complete ORF must start with a valid initiator; cassette fragments
        opt out explicitly."""
        tu = TranslationUnit(11, (Interval(0, 3),), "M")
        assert tu.starts_at_initiator is True


def test_reverse_complement_round_trips() -> None:
    seq = "ATGCGGTTACA"
    assert reverse_complement(reverse_complement(seq)) == seq
    assert reverse_complement("ATGC") == "GCAT"
