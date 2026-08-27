"""Assembly: the CDS/backbone boundary, made structural and then proved.

These are the lane's gate-G2 tests. A wrong coordinate model invalidates
everything built on top of it, so the checks here are deliberately the
adversarial ones: a motif planted so it exists ONLY across the CDS/backbone
junction, another planted so it exists ONLY across the origin, and a backbone
base edited on purpose to confirm I9 actually fires.
"""

from __future__ import annotations

import io
from dataclasses import replace

import pytest
from bt5.core.result import VerificationError
from bt5.core.types import (
    Construct,
    Feature,
    Interval,
    SegmentKind,
    Topology,
    reverse_complement,
)
from bt5.vector import (
    VectorBackbone,
    VectorError,
    assemble,
    construct_to_record,
    read_genbank,
    write_genbank,
)
from bt5.verify import verify_construct
from conftest import make_cds

LONGER = 150  # codons; the fixture's own CDS is 140
SHORTER = 120


def occurs(sequence: str, motif: str) -> bool:
    """Present on either strand of a LINEAR reading of `sequence`."""
    return motif in sequence or reverse_complement(motif) in sequence


def build(
    backbone: VectorBackbone, n_codons: int = LONGER, **kw: object
) -> tuple[object, str, str]:
    cds, protein = make_cds(n_codons)
    assembly = assemble(backbone, cds, protein=protein, table_id=1, **kw)  # type: ignore[arg-type]
    return assembly, cds, protein


class TestBackboneIsUntouched:
    def test_flanking_bases_are_byte_identical(self, backbone: VectorBackbone) -> None:
        """Checked directly, not only through I9, so a bad splice cannot hide."""
        assembly, cds, _ = build(backbone)
        site = assembly.site  # type: ignore[attr-defined]
        construct = assembly.construct  # type: ignore[attr-defined]
        assert construct.sequence[: site.interval.start] == backbone.sequence[: site.interval.start]
        assert (
            construct.sequence[site.interval.start :][len(cds) :]
            == (backbone.sequence[site.interval.end :])
        )

    def test_verify_construct_passes_with_i9_armed(self, backbone: VectorBackbone) -> None:
        assembly, _, protein = build(backbone)
        verify_construct(
            assembly.construct,  # type: ignore[attr-defined]
            protein=protein,
            table_id=1,
            original_backbone=assembly.reference,  # type: ignore[attr-defined]
        )

    def test_i9_catches_a_deliberate_backbone_edit(self, backbone: VectorBackbone) -> None:
        """The worst bug a vector-aware tool can have, made into an exception."""
        assembly, _, protein = build(backbone)
        construct: Construct = assembly.construct  # type: ignore[attr-defined]
        target = 1900  # inside a backbone segment, well clear of the CDS
        edited = "T" if construct.sequence[target] != "T" else "G"
        tampered = replace(
            construct,
            sequence=construct.sequence[:target] + edited + construct.sequence[target + 1 :],
        )
        with pytest.raises(VerificationError) as exc:
            verify_construct(
                tampered,
                protein=protein,
                table_id=1,
                original_backbone=assembly.reference,  # type: ignore[attr-defined]
            )
        assert exc.value.invariant == "I9"
        assert "1900" in str(exc.value)

    def test_the_reference_carries_the_input_backbone(self, backbone: VectorBackbone) -> None:
        assembly, _, _ = build(backbone)
        reference: Construct = assembly.reference  # type: ignore[attr-defined]
        site = assembly.site  # type: ignore[attr-defined]
        assert reference.sequence[: site.interval.start] == backbone.sequence[: site.interval.start]


class TestMotifsAcrossBoundaries:
    def test_a_motif_only_at_the_cds_backbone_junction_is_caught(
        self, backbone: VectorBackbone
    ) -> None:
        assembly, cds, protein = build(backbone)
        construct: Construct = assembly.construct  # type: ignore[attr-defined]
        end = assembly.cds_interval.end  # type: ignore[attr-defined]
        motif = construct.sequence[end - 6 : end + 6]

        assert not occurs(cds, motif), "the motif must not exist inside the CDS alone"
        assert not occurs(backbone.sequence, motif), "nor inside the backbone alone"

        with pytest.raises(VerificationError) as exc:
            verify_construct(construct, protein=protein, table_id=1, forbidden=[motif])
        assert exc.value.invariant == "I6"

    def test_a_motif_only_across_the_origin_is_caught(self, backbone: VectorBackbone) -> None:
        assembly, _, protein = build(backbone)
        construct: Construct = assembly.construct  # type: ignore[attr-defined]
        motif = construct.sequence[-6:] + construct.sequence[:6]

        assert not occurs(construct.sequence, motif), "the motif exists only across the origin"

        with pytest.raises(VerificationError) as exc:
            verify_construct(construct, protein=protein, table_id=1, forbidden=[motif])
        assert exc.value.invariant == "I6"

    def test_a_linear_construct_has_no_origin_junction(self, backbone: VectorBackbone) -> None:
        """The same motif must NOT be reported once the topology says linear."""
        assembly, _, protein = build(backbone)
        construct: Construct = assembly.construct  # type: ignore[attr-defined]
        motif = construct.sequence[-6:] + construct.sequence[:6]
        assert not occurs(construct.sequence, motif)
        linear = replace(construct, topology=Topology.LINEAR)
        verify_construct(linear, protein=protein, table_id=1, forbidden=[motif])


class TestSegments:
    def test_segments_tile_the_construct_exactly_once(self, backbone: VectorBackbone) -> None:
        assembly, _, _ = build(backbone)
        construct: Construct = assembly.construct  # type: ignore[attr-defined]
        cursor = 0
        for segment in construct.segments:
            assert segment.interval.start == cursor, "segments must be contiguous and ordered"
            cursor = segment.interval.end
        assert cursor == construct.length

    def test_the_cds_is_the_only_editable_region(self, backbone: VectorBackbone) -> None:
        assembly, cds, _ = build(backbone)
        construct: Construct = assembly.construct  # type: ignore[attr-defined]
        assert construct.editable == (assembly.cds_interval,)  # type: ignore[attr-defined]
        assert construct.slice(construct.editable[0]) == cds

    def test_ltrs_become_whitelisted_repeats(self, backbone: VectorBackbone) -> None:
        assembly, _, _ = build(backbone)
        construct: Construct = assembly.construct  # type: ignore[attr-defined]
        repeats = [s for s in construct.segments if s.kind is SegmentKind.WHITELISTED_REPEAT]
        assert len(repeats) == 2
        assert {s.label for s in repeats} == {"5' LTR", "3' LTR"}
        assert construct.slice(repeats[0].interval) == construct.slice(repeats[1].interval), (
            "the fixture's LTRs are a genuine 250 bp perfect direct repeat"
        )

    def test_the_annotated_intron_is_exempt(self, backbone: VectorBackbone) -> None:
        assembly, _, _ = build(backbone)
        construct: Construct = assembly.construct  # type: ignore[attr-defined]
        introns = [s for s in construct.segments if s.kind is SegmentKind.ANNOTATED_INTRON]
        assert [s.interval for s in introns] == [Interval(720, 800)]
        assert Interval(720, 800) in construct.exempt

    def test_an_intron_inside_the_cds_is_refused(self, backbone: VectorBackbone) -> None:
        """There are no introns in a BT5 CDS; silently exempting one would hide it."""
        bad = replace(
            backbone,
            features=tuple(
                replace(f, interval=Interval(900, 1000)) if f.kind == "intron" else f
                for f in backbone.features
            ),
        )
        cds, protein = make_cds(LONGER)
        with pytest.raises(VectorError, match="overlaps the designable CDS"):
            assemble(bad, cds, protein=protein, table_id=1)


class TestFeatureRemapping:
    def test_upstream_features_do_not_move(self, backbone: VectorBackbone) -> None:
        assembly, _, _ = build(backbone)
        construct: Construct = assembly.construct  # type: ignore[attr-defined]
        promoter = next(f for f in construct.features if f.kind == "promoter")
        assert promoter.interval == Interval(500, 720, 1)

    def test_downstream_features_shift_by_delta(self, backbone: VectorBackbone) -> None:
        assembly, _, _ = build(backbone)
        delta = assembly.remapper.delta  # type: ignore[attr-defined]
        construct: Construct = assembly.construct  # type: ignore[attr-defined]
        wpre = next(f for f in construct.features if f.qualifiers.get("label") == ("WPRE",))
        assert wpre.interval == Interval(1280 + delta, 1520 + delta, 1)

    def test_the_origin_spanning_feature_still_wraps(self, backbone: VectorBackbone) -> None:
        assembly, _, _ = build(backbone)
        construct: Construct = assembly.construct  # type: ignore[attr-defined]
        ori = next(f for f in construct.features if f.kind == "rep_origin")
        assert ori.interval.end > construct.length
        assert ori.interval.length == 280, "the origin must not change size"
        assert construct.slice(ori.interval) == (
            construct.sequence[ori.interval.start :] + construct.sequence[:80]
        )

    def test_the_source_feature_spans_the_new_length(self, backbone: VectorBackbone) -> None:
        assembly, _, _ = build(backbone)
        construct: Construct = assembly.construct  # type: ignore[attr-defined]
        source = next(f for f in construct.features if f.kind == "source")
        assert source.interval == Interval(0, construct.length)

    def test_the_new_cds_feature_records_the_table(self, backbone: VectorBackbone) -> None:
        assembly, _, _ = build(backbone)
        construct: Construct = assembly.construct  # type: ignore[attr-defined]
        cds_features = [f for f in construct.features if f.kind == "CDS"]
        assert len(cds_features) == 1, "the replaced CDS feature must not survive"
        assert cds_features[0].qualifiers["transl_table"] == ("1",)
        assert cds_features[0].interval == assembly.cds_interval  # type: ignore[attr-defined]

    def test_a_feature_overlapping_the_insert_is_dropped_and_reported(
        self, backbone: VectorBackbone
    ) -> None:
        overlapping = replace(backbone.features[3], interval=Interval(1200, 1400), uid="overlap")
        bb = replace(backbone, features=(*backbone.features, overlapping))
        assembly, _, _ = build(bb)
        construct: Construct = assembly.construct  # type: ignore[attr-defined]
        assert not any(f.uid == "overlap" for f in construct.features)
        assert any(  # type: ignore[attr-defined]
            n.kind == "change" and "overlaps the replaced CDS" in n.summary for n in assembly.notes
        )


class TestLengthChanges:
    @pytest.mark.parametrize("n_codons", [SHORTER, 140, LONGER])
    def test_any_length_assembles_and_verifies(
        self, backbone: VectorBackbone, n_codons: int
    ) -> None:
        assembly, cds, protein = build(backbone, n_codons)
        construct: Construct = assembly.construct  # type: ignore[attr-defined]
        assert construct.length == backbone.length + len(cds) - 420
        verify_construct(
            construct,
            protein=protein,
            table_id=1,
            original_backbone=assembly.reference,  # type: ignore[attr-defined]
        )


class TestOriginStraddlingSite:
    def test_the_vector_is_rotated_so_the_cds_is_contiguous(self, backbone: VectorBackbone) -> None:
        rotated = backbone.rotated(1000)
        assert rotated.find_insertion_site().interval.end > rotated.length
        assembly, cds, protein = build(rotated)
        construct: Construct = assembly.construct  # type: ignore[attr-defined]
        assert assembly.rotation == 1860  # type: ignore[attr-defined]
        assert assembly.cds_interval.end <= construct.length  # type: ignore[attr-defined]
        assert construct.slice(assembly.cds_interval) == cds  # type: ignore[attr-defined]
        verify_construct(
            construct,
            protein=protein,
            table_id=1,
            original_backbone=assembly.reference,  # type: ignore[attr-defined]
        )

    def test_the_rotation_is_reported(self, backbone: VectorBackbone) -> None:
        assembly, _, _ = build(backbone.rotated(1000))
        assert any("rotated" in n.summary for n in assembly.notes)  # type: ignore[attr-defined]


class TestValidation:
    def test_a_table_mismatch_is_refused(self, backbone: VectorBackbone) -> None:
        """Table 12 reassigns CTG to Ser; letting one side win is a silent bug."""
        cds, protein = make_cds(LONGER)
        with pytest.raises(VectorError, match="transl_table=1 .* asked for table 12"):
            assemble(backbone, cds, protein=protein, table_id=12)

    def test_a_frame_shifted_cds_is_refused(self, backbone: VectorBackbone) -> None:
        cds, protein = make_cds(LONGER)
        with pytest.raises(VectorError, match="not a multiple of 3"):
            assemble(backbone, cds[:-1], protein=protein, table_id=1)

    def test_an_ambiguous_cds_is_refused(self, backbone: VectorBackbone) -> None:
        cds, protein = make_cds(LONGER)
        with pytest.raises(VectorError, match="non-ACGT"):
            assemble(backbone, "N" + cds[1:], protein=protein, table_id=1)


class TestExport:
    def test_the_assembled_construct_round_trips_as_genbank(self, backbone: VectorBackbone) -> None:
        assembly, _, _ = build(backbone)
        construct: Construct = assembly.construct  # type: ignore[attr-defined]
        text = write_genbank(
            construct_to_record(construct, name="designed", compound_parts=assembly.compound_parts)  # type: ignore[attr-defined]
        )
        again = read_genbank(io.StringIO(text))
        assert again.sequence == construct.sequence
        assert again.topology == construct.topology
        assert {(f.kind, f.interval) for f in again.features} == {
            (f.kind, f.interval) for f in construct.features
        }


class TestReverseOrientedCassette:
    """A reverse-oriented cassette is the case a strand-blind assembler gets
    exactly backwards, and it is the common one in lentiviral vectors.

    Inserting the coding sequence verbatim and walking codons in ascending
    coordinates yields the reverse-complement protein. That still passes a naive
    per-codon check, so the real assertion here is that the independent oracle
    accepts the round trip.
    """

    def backbone(self) -> VectorBackbone:
        return VectorBackbone(
            sequence="ACGT" * 150,
            topology=Topology.CIRCULAR,
            features=(
                Feature(
                    interval=Interval(200, 290, -1),
                    kind="CDS",
                    qualifiers={"label": ("reverse transgene",), "transl_table": ("1",)},
                    uid="f0",
                ),
            ),
            name="reverse",
        )

    def test_the_plus_strand_carries_the_reverse_complement(self) -> None:
        cds, protein = make_cds(30)
        assembly = assemble(self.backbone(), cds, protein=protein, table_id=1)
        assert assembly.construct.sequence[200 : 200 + len(cds)] == reverse_complement(cds)

    def test_the_editable_region_slices_back_to_the_coding_sequence(self) -> None:
        cds, protein = make_cds(30)
        assembly = assemble(self.backbone(), cds, protein=protein, table_id=1)
        construct = assembly.construct
        assert construct.editable[0].strand == -1
        assert construct.slice(construct.editable[0]) == cds

    def test_the_codon_map_runs_high_to_low(self) -> None:
        cds, protein = make_cds(30)
        assembly = assemble(self.backbone(), cds, protein=protein, table_id=1)
        codons = assembly.construct.translation_units[0].codon_map
        assert codons[0] == Interval(287, 290, -1), "the first codon is at the HIGH end"
        assert all(c.strand == -1 for c in codons)
        assert [c.start for c in codons] == sorted((c.start for c in codons), reverse=True)

    def test_the_oracle_accepts_the_round_trip(self) -> None:
        cds, protein = make_cds(30)
        assembly = assemble(self.backbone(), cds, protein=protein, table_id=1)
        verify_construct(
            assembly.construct,
            protein=protein,
            table_id=1,
            original_backbone=assembly.reference,
        )

    def test_the_exported_cds_feature_is_on_the_complement(self) -> None:
        cds, protein = make_cds(30)
        assembly = assemble(self.backbone(), cds, protein=protein, table_id=1)
        text = write_genbank(construct_to_record(assembly.construct))
        assert "complement(201..290)" in text
