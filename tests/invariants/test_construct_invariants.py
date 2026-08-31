"""Invariants that only exist because BT5 evaluates the ASSEMBLED construct."""

from __future__ import annotations

import pytest
from bt5.core.result import VerificationError
from bt5.core.types import (
    Construct,
    Interval,
    Segment,
    SegmentKind,
    Topology,
    TranslationUnit,
    reverse_complement,
)
from bt5.verify import find_motifs, gc_fraction, verify_construct

CDS = "ATGAAACCCTAA"
BACKBONE = "GGATCCAAGCTTGTCGACCTGCAG"


def build(
    seq: str = CDS + BACKBONE,
    backbone_kind: SegmentKind = SegmentKind.BACKBONE,
) -> Construct:
    return Construct(
        sequence=seq,
        topology=Topology.CIRCULAR,
        segments=(
            Segment(Interval(0, 12), SegmentKind.DESIGNABLE_CDS, "insert"),
            Segment(Interval(12, len(seq)), backbone_kind, "vector"),
        ),
        translation_units=(
            TranslationUnit(11, tuple(Interval(i, i + 3) for i in range(0, 12, 3)), "MKP", True),
        ),
    )


def test_baseline_verifies() -> None:
    c = build()
    verify_construct(c, protein="MKP", table_id=11, original_backbone=c)


def test_i9_catches_a_single_base_backbone_edit() -> None:
    """The worst possible bug in a vector-context tool: silently editing the
    user's own vector. I9 makes it an exception rather than a shipped plasmid."""
    original = build()
    seq = original.sequence
    tampered = build(seq[:20] + ("A" if seq[20] != "A" else "C") + seq[21:])
    with pytest.raises(VerificationError) as exc:
        verify_construct(tampered, protein="MKP", table_id=11, original_backbone=original)
    assert exc.value.invariant == "I9"


def test_i6_finds_a_motif_spanning_the_origin() -> None:
    """The construct ends ...CTGCAG and begins ATGAAA, so the circular junction
    reads CTGCAGATGAAA. CAGATG exists ONLY across that junction."""
    c = build()
    assert "CAGATG" not in c.sequence
    assert find_motifs(c, ["CAGATG"]), "origin-spanning motif must be found"


def test_i6_finds_a_motif_only_present_on_the_reverse_strand() -> None:
    """A motif absent from the forward strand but present as its reverse
    complement must still be found.

    Every classic six-cutter is palindromic, so proving this needs a
    non-palindromic probe: GGATCCA sits at the CDS/backbone junction, and its
    reverse complement TGGATCC does NOT occur in the forward sequence.
    """
    c = build()
    forward = "GGATCCA"
    probe = reverse_complement(forward)  # TGGATCC
    assert probe != forward, "probe must be non-palindromic for this test to mean anything"
    assert forward in c.sequence
    assert probe not in c.sequence, "probe must be absent from the forward strand"
    assert find_motifs(c, [probe]), "a motif present only as a revcomp must be found"


def test_i6_finds_a_motif_spanning_the_cds_backbone_junction() -> None:
    """CDS ends ...CCCTAA, backbone begins GGATCC -> junction reads CCCTAAGGATCC."""
    c = build()
    assert find_motifs(c, ["TAAGGA"]), "junction-spanning motif must be found"


def test_whitelisted_repeats_are_exempt_from_scanning() -> None:
    """LTRs and ITRs irreducibly violate the repeat rules. The answer is a strain
    and temperature protocol, not a redesign, so they are immutable AND unscanned."""
    plain = build()
    exempted = build(backbone_kind=SegmentKind.WHITELISTED_REPEAT)
    assert find_motifs(plain, ["GGATCC"]), "site in plain backbone must be flagged"
    assert not find_motifs(exempted, ["GGATCC"]), "site in a whitelisted LTR must be ignored"


def test_interior_stop_is_refused() -> None:
    c = Construct(
        sequence="ATGTAACCCTAA",
        topology=Topology.LINEAR,
        segments=(Segment(Interval(0, 12), SegmentKind.DESIGNABLE_CDS, "cds"),),
        translation_units=(
            TranslationUnit(11, tuple(Interval(i, i + 3) for i in range(0, 12, 3)), "M*P", True),
        ),
    )
    with pytest.raises(VerificationError) as exc:
        verify_construct(c, protein="M*P", table_id=11)
    assert exc.value.invariant == "I5"


class TestI7OrderedSpanScope:
    """I7 holds the ORDERED spans to the band, not the plasmid, and not windows.

    These replace `test_gc_band_window_wraps_the_origin`, which asserted that a
    20 bp window at 100% GC straddling the origin must make the oracle refuse a
    construct whose overall GC was exactly 0.50. That test fused two claims, and
    they part company.

    The first -- a single extreme window is a breach -- is FALSE as measured.
    Sixteen of sixteen ladder probes carrying one extreme 50 bp window (4-96%
    local GC) were accepted by both IDT and Twist, and a probe with a 100 bp
    window at 23% was accepted while one at 21% GLOBAL was refused
    (docs/design/vendor-gc-calibration.md). The oracle was refusing DNA both
    vendors manufacture. `test_an_extreme_window_inside_an_in_band_span_passes`
    now asserts the measured truth in that test's place, so the change is
    recorded as an assertion rather than an absence.

    The second -- a circular scan must handle the origin -- is still true, and
    now has a HARDER test: a designable span can itself cross the origin, so
    `test_an_ordered_span_that_wraps_the_origin_is_measured` builds one whose
    head and tail differ in GC. A wrap-blind slice reads 0.50 and passes; the
    real ordered DNA reads 0.75 and must be refused.
    """

    @staticmethod
    def _circular(seq: str, segments: tuple[Segment, ...]) -> Construct:
        return Construct(sequence=seq, topology=Topology.CIRCULAR, segments=segments)

    def test_a_gc_extreme_insert_is_refused_though_the_plasmid_is_in_band(self) -> None:
        """The motivating case. A near-neutral backbone dilutes a 90% insert to
        0.47 across the construct, so the old whole-construct check passed a
        fragment both vendors deny."""
        insert = "GGC" * 10  # 30 bp at 100% -> the ordered span
        backbone = "AT" * 35  # 70 bp at 0% -> nobody synthesises this
        seq = insert + backbone
        c = self._circular(
            seq,
            (
                Segment(Interval(0, 30), SegmentKind.DESIGNABLE_CDS, "insert"),
                Segment(Interval(30, 100), SegmentKind.BACKBONE, "vector"),
            ),
        )
        assert 0.28 <= gc_fraction(seq) <= 0.77, "the PLASMID is in band; the tube is not"
        with pytest.raises(VerificationError) as exc:
            verify_construct(c, protein="", table_id=11, gc_bounds=(0.28, 0.77))
        assert exc.value.invariant == "I7"
        assert "1.000" in str(exc.value)

    def test_an_extreme_window_inside_an_in_band_span_passes(self) -> None:
        """A window is never a breach. This is the behaviour change, stated as an
        assertion so a future reversion fails loudly."""
        seq = "GC" * 25 + "AT" * 25  # 100 bp, global 0.50, two 50 bp windows at 1.0 and 0.0
        assert gc_fraction(seq[:50]) == 1.0
        assert gc_fraction(seq[50:]) == 0.0
        c = self._circular(seq, (Segment(Interval(0, 100), SegmentKind.DESIGNABLE_CDS, "cds"),))
        verify_construct(c, protein="", table_id=11, gc_bounds=(0.40, 0.60))

    def test_an_ordered_span_that_wraps_the_origin_is_measured(self) -> None:
        """The direct descendant of the retired test. Fails under whole-construct
        scope AND under wrap-blind slicing."""
        seq = "GGGGGGGGGG" + ("AT" * 15 + "GC" * 15) + "ACACACACAC"  # 80 bp
        c = self._circular(
            seq,
            (
                Segment(Interval(70, 90), SegmentKind.DESIGNABLE_CDS, "insert"),  # wraps
                Segment(Interval(10, 70), SegmentKind.BACKBONE, "vector"),
            ),
        )
        assert 0.40 <= gc_fraction(seq) <= 0.60, "the plasmid is in band"
        assert gc_fraction(seq[70:90]) == 0.50, "a wrap-blind slice reads in band and passes"
        with pytest.raises(VerificationError) as exc:
            verify_construct(c, protein="", table_id=11, gc_bounds=(0.40, 0.60))
        assert exc.value.invariant == "I7"
        assert "0.750" in str(exc.value), "the WRAPPED bases are 75% GC"

    def test_every_designable_span_is_held_on_its_own(self) -> None:
        """Two tubes are two synthesis reactions. Averaging into band is not a
        defence for either of them."""
        seq = "GGGAAAAAAA" * 3 + "ATGC" * 10 + "GGGGGGGAAA" * 3  # 0.30 / 0.50 / 0.70
        c = self._circular(
            seq,
            (
                Segment(Interval(0, 30), SegmentKind.DESIGNABLE_CDS, "insert one"),
                Segment(Interval(30, 70), SegmentKind.BACKBONE, "vector"),
                Segment(Interval(70, 100), SegmentKind.DESIGNABLE_CDS, "insert two"),
            ),
        )
        assert gc_fraction(seq) == 0.50, "dead centre of the band, and both tubes are outside it"
        with pytest.raises(VerificationError) as exc:
            verify_construct(c, protein="", table_id=11, gc_bounds=(0.40, 0.60))
        assert exc.value.invariant == "I7"
        assert "0.300" in str(exc.value), "the first span in sorted order"

    def test_a_backbone_the_user_owns_is_not_held_to_the_band(self) -> None:
        """The deliberate relaxation, asserted rather than emergent. I9 already
        proves BT5 never touched these bases, and no vendor is asked to
        synthesise them, so refusing over them would refuse the user's plasmid
        for a property they cannot change and nobody measured."""
        seq = "ATGCATGCAT" * 3 + "GC" * 20 + "ATGCATGCAT" * 3  # spans 0.50, backbone 1.00
        c = self._circular(
            seq,
            (
                Segment(Interval(0, 30), SegmentKind.DESIGNABLE_CDS, "insert one"),
                Segment(Interval(30, 70), SegmentKind.BACKBONE, "vector"),
                Segment(Interval(70, 100), SegmentKind.DESIGNABLE_CDS, "insert two"),
            ),
        )
        assert gc_fraction(seq) > 0.60, "the PLASMID is out of band and that is not I7's business"
        verify_construct(c, protein="", table_id=11, gc_bounds=(0.40, 0.60))

    def test_a_band_on_a_construct_with_nothing_to_order_is_refused(self) -> None:
        """Not a silent skip. A CDS mislabelled BACKBONE would switch I7 off with
        no sound at all -- the same armed-and-vacuous shape this PR removes."""
        seq = "ATGC" * 25
        c = self._circular(seq, (Segment(Interval(0, 100), SegmentKind.BACKBONE, "vector"),))
        with pytest.raises(VerificationError) as exc:
            verify_construct(c, protein="", table_id=11, gc_bounds=(0.40, 0.60))
        assert exc.value.invariant == "I7"
        assert "no designable segment" in str(exc.value)

    def test_no_band_requested_means_no_check(self) -> None:
        seq = "ATGC" * 25
        c = self._circular(seq, (Segment(Interval(0, 100), SegmentKind.BACKBONE, "vector"),))
        verify_construct(c, protein="", table_id=11)


class TestTheCallersDeclarationIsChecked:
    """`protein` and `table_id` were required and read by nothing."""

    def test_a_declared_protein_with_no_translation_unit_is_refused(self) -> None:
        """Previously passed vacuously: I3/I4/I5 loop over `translation_units or
        ()`, so a construct with none was verified against nothing."""
        c = Construct(
            sequence=CDS + BACKBONE,
            topology=Topology.CIRCULAR,
            segments=(Segment(Interval(0, len(CDS)), SegmentKind.DESIGNABLE_CDS, "cds"),),
        )
        with pytest.raises(VerificationError) as exc:
            verify_construct(c, protein="MKP", table_id=11)
        assert exc.value.invariant == "I3"
        assert "no translation unit" in str(exc.value)

    def test_a_declared_table_no_unit_uses_is_refused(self) -> None:
        """The silent one. Asking for table 12 against a unit declaring 11 was
        verified under 11, while the caller believed CTG=Ser had been checked."""
        c = build()
        with pytest.raises(VerificationError) as exc:
            verify_construct(c, protein="MKP", table_id=12)
        assert exc.value.invariant == "I3"
        assert "table 12" in str(exc.value)
