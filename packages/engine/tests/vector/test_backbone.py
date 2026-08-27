"""Insertion-site and untranslated-region detection.

The 5'UTR feeds the highest-weight objective in BT5, and the window that carries
it spans the UTR/CDS junction. So the two things that matter here are that the
detector is strand-aware, and that it says "absent" out loud instead of quietly
handing back something that is not a 5'UTR.
"""

from __future__ import annotations

from dataclasses import replace

import pytest
from bt5.core.types import Feature, Interval, Topology
from bt5.vector import VectorBackbone, VectorError, insertion_site_from_interval
from bt5.vector.backbone import rotate_interval


def feature(kind: str, start: int, end: int, strand: int = 1, label: str = "") -> Feature:
    return Feature(
        interval=Interval(start, end, strand),  # type: ignore[arg-type]
        kind=kind,
        qualifiers={"label": (label or kind,)},
        uid=f"{kind}{start}",
    )


def without(backbone: VectorBackbone, *kinds: str) -> VectorBackbone:
    drop = {k.lower() for k in kinds}
    return replace(
        backbone, features=tuple(f for f in backbone.features if f.kind.lower() not in drop)
    )


class TestInsertionSite:
    def test_finds_the_single_cds(self, backbone: VectorBackbone) -> None:
        site = backbone.find_insertion_site()
        assert site.interval == Interval(860, 1280, 1)
        assert site.label == "transgene"
        assert site.source == "annotated_cds"

    def test_reads_the_declared_translation_table(self, backbone: VectorBackbone) -> None:
        """The table comes from the file, and is never defaulted here."""
        assert backbone.find_insertion_site().detected_table_id == 1

    def test_an_unlabelled_table_stays_none(self, backbone: VectorBackbone) -> None:
        stripped = tuple(
            replace(f, qualifiers={k: v for k, v in f.qualifiers.items() if k != "transl_table"})
            if f.kind == "CDS"
            else f
            for f in backbone.features
        )
        site = replace(backbone, features=stripped).find_insertion_site()
        assert site.detected_table_id is None, "silence in the file must not become a guess"

    def test_no_cds_is_an_error_with_a_way_out(self, backbone: VectorBackbone) -> None:
        with pytest.raises(VectorError, match="no CDS feature"):
            without(backbone, "CDS").find_insertion_site()

    def test_two_cds_features_refuse_to_pick_one(self, backbone: VectorBackbone) -> None:
        """Choosing wrong here yields a plausible plasmid expressing the wrong gene."""
        second = feature("CDS", 100, 200, label="other")
        two = replace(backbone, features=(*backbone.features, second))
        with pytest.raises(VectorError, match="2 CDS features"):
            two.find_insertion_site()

    def test_a_label_disambiguates(self, backbone: VectorBackbone) -> None:
        second = feature("CDS", 100, 200, label="other")
        two = replace(backbone, features=(*backbone.features, second))
        assert two.find_insertion_site(label="transgene").interval == Interval(860, 1280, 1)


class TestFivePrimeUtr:
    def test_prefers_the_annotated_feature(self, backbone: VectorBackbone) -> None:
        utr = backbone.utr_context(backbone.find_insertion_site())
        assert utr.five_prime == Interval(800, 860, 1)
        assert utr.five_prime_source == "annotated_feature"

    def test_falls_back_to_the_promoter(self, backbone: VectorBackbone) -> None:
        stripped = without(backbone, "5'UTR")
        utr = stripped.utr_context(stripped.find_insertion_site())
        assert utr.five_prime == Interval(720, 860, 1)
        assert utr.five_prime_source == "derived_from_promoter"
        assumed = [n for n in utr.notes if n.kind == "assumption"]
        assert assumed, "inferring a transcription start must be stated, not silent"
        assert all(n.bears_on == "protein expression" for n in assumed)

    def test_a_derived_utr_containing_an_intron_says_so(self, backbone: VectorBackbone) -> None:
        """The mature 5'UTR after splicing is shorter than the span used."""
        stripped = without(backbone, "5'UTR")
        utr = stripped.utr_context(stripped.find_insertion_site())
        assert any("intron" in n.summary and "shorter" in n.summary for n in utr.notes)

    def test_absent_utr_degrades_instead_of_guessing(self, backbone: VectorBackbone) -> None:
        """G6: with no UTR the objective is unavailable, not silently CDS-only."""
        stripped = without(backbone, "5'UTR", "promoter")
        utr = stripped.utr_context(stripped.find_insertion_site())
        assert utr.five_prime is None
        assert not utr.has_five_prime
        assert utr.five_prime_source == "absent"
        unavailable = [n for n in utr.notes if n.kind == "unavailable"]
        assert any("5' folding objective" in n.summary for n in unavailable)
        assert all(n.action for n in unavailable if "5' folding" in n.summary), (
            "telling the user an objective is unavailable is only useful with a way out"
        )

    def test_an_implausibly_long_derived_utr_is_rejected(self, backbone: VectorBackbone) -> None:
        stripped = without(backbone, "5'UTR")
        utr = stripped.utr_context(stripped.find_insertion_site(), max_derived_utr=10)
        assert utr.five_prime is None
        assert any("max_derived_utr=10" in n.summary for n in utr.notes)


class TestReverseStrandCassette:
    """For a reverse-oriented cassette the 5' side is at HIGHER coordinates.

    A strand-blind detector picks the decoy below -- which is downstream of the
    CDS -- and calls it the 5'UTR, silently, on exactly the lentiviral layouts
    this tool exists to serve.
    """

    def backbone(self) -> VectorBackbone:
        return VectorBackbone(
            sequence="ACGT" * 150,
            topology=Topology.CIRCULAR,
            features=(
                feature("5'UTR", 140, 200, -1, "decoy downstream UTR"),
                feature("CDS", 200, 500, -1, "reverse transgene"),
                feature("5'UTR", 500, 560, -1, "real 5' UTR"),
                feature("promoter", 560, 600, -1, "reverse promoter"),
            ),
            name="reverse",
        )

    def test_upstream_is_at_higher_coordinates(self) -> None:
        bb = self.backbone()
        site = bb.find_insertion_site()
        assert site.strand == -1
        utr = bb.utr_context(site)
        assert utr.five_prime == Interval(500, 560, -1), "picked the decoy at lower coordinates"

    def test_a_forward_feature_is_not_used_for_a_reverse_cassette(self) -> None:
        bb = self.backbone()
        forward_only = replace(
            bb, features=(bb.features[1], feature("5'UTR", 500, 560, 1, "forward utr"))
        )
        utr = forward_only.utr_context(forward_only.find_insertion_site())
        assert utr.five_prime is None


class TestCircularAdjacency:
    def test_a_utr_abutting_across_the_origin_is_found(self, backbone: VectorBackbone) -> None:
        """Rotating so the CDS straddles the origin must not lose the UTR."""
        rotated = backbone.rotated(1000)
        site = rotated.find_insertion_site()
        assert site.interval.end > rotated.length, "the site now wraps"
        utr = rotated.utr_context(site)
        assert utr.five_prime_source == "annotated_feature"
        assert utr.five_prime == Interval(1800, 1860, 1)


class TestRotation:
    def test_rotation_preserves_every_span(self, backbone: VectorBackbone) -> None:
        rotated = backbone.rotated(1000)
        assert rotated.length == backbone.length
        for before, after in zip(backbone.features, rotated.features, strict=True):
            assert after.interval.length == before.interval.length
            assert rotated.slice(after.interval) == backbone.slice(before.interval)

    def test_rotation_is_recorded_not_silent(self, backbone: VectorBackbone) -> None:
        notes = backbone.rotated(1000).notes
        assert any(n.kind == "change" and "rotated" in n.summary for n in notes)

    def test_rotation_moves_existing_note_intervals(self) -> None:
        """A located note must not end up pointing at unrelated sequence."""
        from bt5.vector.notes import DesignNote

        bb = VectorBackbone(
            sequence="ACGT" * 150,
            topology=Topology.CIRCULAR,
            features=(feature("CDS", 200, 500),),
            notes=(DesignNote(kind="liability", summary="x", interval=Interval(10, 40)),),
        )
        moved = next(n for n in bb.rotated(100).notes if n.summary == "x")
        assert moved.interval is not None
        assert bb.rotated(100).slice(moved.interval) == bb.slice(Interval(10, 40))

    def test_rotating_by_zero_is_the_same_object(self, backbone: VectorBackbone) -> None:
        assert backbone.rotated(0) is backbone

    def test_a_linear_vector_cannot_be_rotated(self, backbone: VectorBackbone) -> None:
        linear = replace(backbone, topology=Topology.LINEAR)
        with pytest.raises(VectorError, match="only a circular vector"):
            linear.rotated(10)

    def test_rotate_interval_wraps_rather_than_going_negative(self) -> None:
        assert rotate_interval(Interval(10, 40), by=20, length=100) == Interval(90, 120)


class TestExplicitSite:
    def test_marking_by_hand_carries_no_detected_table(self) -> None:
        site = insertion_site_from_interval(Interval(10, 40), label="manual")
        assert site.source == "explicit"
        assert site.detected_table_id is None
