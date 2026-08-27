"""Insertion-site and untranslated-region detection.

The 5'UTR feeds the highest-weight objective in BT5, and the window that carries
it spans the UTR/CDS junction. So the two things that matter here are that the
detector is strand-aware, and that it says "absent" out loud instead of quietly
handing back something that is not a 5'UTR.
"""

from __future__ import annotations

from dataclasses import replace

import pytest
from bt5.core.types import Feature, Interval, Topology, reverse_complement
from bt5.vector import VectorBackbone, VectorError, insertion_site_from_interval
from bt5.vector.backbone import rotate_interval


def feature(
    kind: str,
    start: int,
    end: int,
    strand: int = 1,
    label: str = "",
    **extra: tuple[str, ...],
) -> Feature:
    return Feature(
        interval=Interval(start, end, strand),  # type: ignore[arg-type]
        kind=kind,
        qualifiers={"label": (label or kind,), **extra},
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


class TestRibosomeBindingSite:
    """Bacteria do not scan, so the promoter is the wrong 5' anchor for them.

    Cap-dependent scanning makes the transcription start the thing that defines a
    eukaryotic leader. A bacterial ribosome is recruited by the Shine-Dalgarno,
    and pET vectors annotate it. A detector that knows only about promoters
    reports "the transcription start is assumed" on a vector whose initiation
    element is annotated 30 nt from the start codon -- true, and beside the point.
    """

    #: T7 gene 10 leader, as pET28a carries it. The SD is the trailing AAGGAGA.
    SD = "TTTGTTTAACTTTAAGAAGGAGA"
    PROMOTER = "TAATACGACTCACTATAGG"

    def backbone(
        self,
        *,
        spacer: int = 6,
        with_promoter: bool = True,
        strand: int = 1,
        rbs_kind: str = "RBS",
        rbs_class: str | None = None,
    ) -> tuple[VectorBackbone, dict[str, Interval]]:
        """A minimal bacterial cassette, with every offset derived not counted."""
        lead = "AC" * 30
        gap = "CTAG" * 20
        cds = "ATG" + "GCT" * 40 + "TAA"
        parts = [lead, self.PROMOTER, gap, self.SD, "T" * spacer, cds]
        forward = "".join(parts) + "GC" * 40

        starts: list[int] = []
        at = 0
        for part in parts:
            starts.append(at)
            at += len(part)
        promoter_iv = Interval(starts[1], starts[1] + len(self.PROMOTER), strand)
        rbs_iv = Interval(starts[3], starts[3] + len(self.SD), strand)
        cds_iv = Interval(starts[5], starts[5] + len(cds), strand)

        if strand == -1:
            # Rebuild on the other strand: every span mirrors, so a detector that
            # is right by accident on the forward layout is caught here.
            n = len(forward)
            sequence = reverse_complement(forward)
            flip = lambda iv: Interval(n - iv.end, n - iv.start, -1)  # noqa: E731
            promoter_iv, rbs_iv, cds_iv = flip(promoter_iv), flip(rbs_iv), flip(cds_iv)
        else:
            sequence = forward

        extra = {"regulatory_class": (rbs_class,)} if rbs_class else {}
        features = [
            feature(rbs_kind, rbs_iv.start, rbs_iv.end, strand, "T7 g10 RBS", **extra),
            feature("CDS", cds_iv.start, cds_iv.end, strand, "transgene"),
        ]
        if with_promoter:
            features.insert(
                0,
                feature("promoter", promoter_iv.start, promoter_iv.end, strand, "T7 promoter"),
            )
        bb = VectorBackbone(
            sequence=sequence,
            topology=Topology.CIRCULAR,
            features=tuple(features),
            name="synthetic pET",
        )
        return bb, {"promoter": promoter_iv, "rbs": rbs_iv, "cds": cds_iv}

    def test_an_annotated_rbs_is_reported_with_its_spacing(self) -> None:
        bb, spans = self.backbone(spacer=6)
        utr = bb.utr_context(bb.find_insertion_site())
        assert utr.ribosome_binding_site == spans["rbs"]
        assert utr.has_ribosome_binding_site
        assert utr.rbs_spacing == 6
        assert bb.slice(spans["rbs"]).endswith("AAGGAGA")

    @pytest.mark.parametrize(
        ("kind", "regulatory_class"),
        [
            # SnapGene, and every Addgene deposit that came through it.
            ("RBS", None),
            # INSDC, which deprecated the bare key in favour of a class.
            ("regulatory", "ribosome_binding_site"),
            ("regulatory", "shine_dalgarno_sequence"),
            # Sequence Ontology type name, as a GFF3-derived import carries it.
            ("ribosome_binding_site", None),
        ],
    )
    def test_every_spelling_of_an_rbs_is_equivalent(
        self, kind: str, regulatory_class: str | None
    ) -> None:
        """Which spelling a file uses is a fact about its exporter, not its biology."""
        snapgene, _ = self.backbone(rbs_kind="RBS")
        other, _ = self.backbone(rbs_kind=kind, rbs_class=regulatory_class)
        expected = snapgene.utr_context(snapgene.find_insertion_site())
        got = other.utr_context(other.find_insertion_site())
        assert got.ribosome_binding_site == expected.ribosome_binding_site
        assert got.rbs_spacing == expected.rbs_spacing

    def test_an_untyped_regulatory_feature_is_not_an_rbs(self) -> None:
        """`regulatory` alone covers terminators and operators too."""
        bb, _ = self.backbone(rbs_kind="regulatory", rbs_class=None)
        assert bb.utr_context(bb.find_insertion_site()).ribosome_binding_site is None

    def test_the_rbs_anchors_a_leader_when_no_promoter_is_annotated(self) -> None:
        """Without this the whole objective is dropped over a missing promoter."""
        bb, spans = self.backbone(with_promoter=False, spacer=6)
        utr = bb.utr_context(bb.find_insertion_site())
        assert utr.five_prime_source == "derived_from_rbs"
        assert utr.five_prime is not None
        assert bb.slice(utr.five_prime).startswith(self.SD), (
            "an anchored leader that starts AFTER the SD drops the 7 bases that "
            "recruit the ribosome -- the whole reason for anchoring on it"
        )
        assert utr.five_prime.length == len(self.SD) + 6

    def test_an_anchored_leader_says_it_stops_at_the_sd(self) -> None:
        bb, _ = self.backbone(with_promoter=False)
        notes = bb.utr_context(bb.find_insertion_site()).notes
        assert any(
            n.kind == "assumption"
            and "anchored on the annotated ribosome binding site" in n.summary
            for n in notes
        )

    def test_no_rbs_and_no_promoter_still_degrades(self) -> None:
        """The fallback must not paper over a vector that really has nothing."""
        bb, _ = self.backbone(with_promoter=False)
        bare = without(bb, "RBS")
        utr = bare.utr_context(bare.find_insertion_site())
        assert utr.five_prime is None
        assert utr.five_prime_source == "absent"
        assert utr.ribosome_binding_site is None

    def test_a_promoter_leader_wins_but_names_the_rbs_inside_it(self) -> None:
        """The longer real leader is better; the flat 'assumed' wording was not."""
        bb, spans = self.backbone(with_promoter=True)
        utr = bb.utr_context(bb.find_insertion_site())
        assert utr.five_prime_source == "derived_from_promoter"
        assert utr.five_prime is not None
        assert utr.five_prime.start == spans["promoter"].end
        assert utr.ribosome_binding_site == spans["rbs"]
        assert any("ribosome binding site lies inside this span" in n.summary for n in utr.notes), (
            "reporting only 'the transcription start is assumed' buries the known element"
        )

    def test_an_rbs_out_of_reach_of_the_start_codon_is_a_liability(self) -> None:
        """pET28a annotates a CDS 66 nt past its SD: that CDS is not the real ORF."""
        bb, _ = self.backbone(spacer=80)
        utr = bb.utr_context(bb.find_insertion_site())
        assert utr.rbs_spacing == 80
        liabilities = [n for n in utr.notes if n.kind == "liability"]
        assert any("80 nt" in n.summary and "initiates" in n.summary for n in liabilities)
        assert all(n.action for n in liabilities)

    def test_an_unreachable_rbs_does_not_corroborate_the_start_codon(self) -> None:
        """It still falls inside the derived span, which is not the same thing."""
        bb, _ = self.backbone(spacer=80)
        utr = bb.utr_context(bb.find_insertion_site())
        assert not any("lies inside this span" in n.summary for n in utr.notes)

    def test_an_rbs_outside_the_leader_is_a_liability(self) -> None:
        """Two 5' annotations describing different transcripts is worth saying."""
        bb, spans = self.backbone(with_promoter=False)
        annotated = replace(
            bb,
            features=(
                *bb.features,
                feature("5'UTR", spans["rbs"].end, spans["cds"].start, 1, "short leader"),
            ),
        )
        utr = annotated.utr_context(annotated.find_insertion_site())
        assert utr.five_prime_source == "annotated_feature"
        assert any(n.kind == "liability" and "outside the 5'UTR" in n.summary for n in utr.notes)

    def test_a_reverse_cassette_finds_its_rbs_at_higher_coordinates(self) -> None:
        """Mirror the layout: a strand-blind search finds nothing, or the wrong thing."""
        bb, spans = self.backbone(strand=-1, spacer=6)
        site = bb.find_insertion_site()
        assert site.strand == -1
        utr = bb.utr_context(site)
        assert utr.ribosome_binding_site == spans["rbs"]
        assert utr.rbs_spacing == 6
        assert utr.ribosome_binding_site is not None
        assert utr.ribosome_binding_site.start >= site.interval.end
        assert bb.slice(utr.ribosome_binding_site).endswith("AAGGAGA")

    def test_a_forward_rbs_is_not_used_for_a_reverse_cassette(self) -> None:
        bb, spans = self.backbone(strand=-1)
        forward_rbs = replace(
            bb,
            features=tuple(
                feature("RBS", f.interval.start, f.interval.end, 1) if f.kind == "RBS" else f
                for f in bb.features
            ),
        )
        utr = forward_rbs.utr_context(forward_rbs.find_insertion_site())
        assert utr.ribosome_binding_site is None

    def test_an_rbs_across_the_origin_is_found(self) -> None:
        """Circular adjacency: the SD may sit the other side of position 0."""
        bb, spans = self.backbone(with_promoter=False, spacer=6)
        # Rotate so the origin lands between the SD and the start codon.
        rotated = bb.rotated(spans["cds"].start - 2)
        site = rotated.find_insertion_site()
        utr = rotated.utr_context(site)
        assert utr.rbs_spacing == 6, "linear arithmetic loses an SD across the origin"
        assert utr.five_prime_source == "derived_from_rbs"
        assert utr.five_prime is not None
        assert utr.five_prime.end > rotated.length, "the anchored leader wraps"
        assert rotated.slice(utr.five_prime).startswith(self.SD)


class TestDownstreamRegulatoryClass:
    """`regulatory` is a catch-all, so matching it on the key alone reads
    whichever of terminator, attenuator or operator happens to be nearest as the
    polyA signal -- and then hands back a 3'UTR that is not one."""

    def backbone(self, downstream_class: str | None) -> VectorBackbone:
        extra = {"regulatory_class": (downstream_class,)} if downstream_class else {}
        return VectorBackbone(
            sequence="ACGT" * 200,
            topology=Topology.CIRCULAR,
            features=(
                feature("CDS", 100, 400, 1, "transgene"),
                feature("regulatory", 500, 540, 1, "downstream element", **extra),
            ),
            name="regulatory",
        )

    def test_a_polya_class_regulatory_feature_bounds_the_three_prime_utr(self) -> None:
        bb = self.backbone("polyA_signal_sequence")
        utr = bb.utr_context(bb.find_insertion_site())
        assert utr.three_prime == Interval(400, 500, 1)

    def test_a_terminator_is_not_read_as_a_polya_signal(self) -> None:
        bb = self.backbone("terminator")
        assert bb.utr_context(bb.find_insertion_site()).three_prime is None

    def test_an_untyped_regulatory_feature_is_not_read_as_a_polya_signal(self) -> None:
        bb = self.backbone(None)
        assert bb.utr_context(bb.find_insertion_site()).three_prime is None
