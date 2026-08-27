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
from bt5.verify import find_motifs, verify_construct

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


def test_gc_band_window_wraps_the_origin() -> None:
    """A GC-extreme region straddling the origin must be caught on a circular
    construct. Half the run sits at each end of the linear string."""
    seq = "GGGGGGGGGG" + "ATATATATATATATATATAT" + "CCCCCCCCCC"
    c = Construct(
        sequence=seq,
        topology=Topology.CIRCULAR,
        segments=(Segment(Interval(0, len(seq)), SegmentKind.DESIGNABLE_CDS, "cds"),),
    )
    with pytest.raises(VerificationError) as exc:
        verify_construct(c, protein="", table_id=11, gc_bounds=(0.40, 0.60), gc_window=20)
    assert exc.value.invariant == "I7"
