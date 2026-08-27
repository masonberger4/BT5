"""Finding the insertion site on a vector somebody else annotated.

The rule these replace -- "take the one CDS feature" -- fails on every real
transfer vector tried, because the transgene is routinely not annotated as a CDS
while the selection marker always is. So the cases that matter are the ones where
a naive detector picks the wrong thing.
"""

from __future__ import annotations

from dataclasses import replace

import pytest
from bt5.core.types import Feature, Interval, Topology, reverse_complement
from bt5.vector import VectorBackbone, cloning_sites, find_orfs, suggest_insertion_sites
from bt5.vector.candidates import CONFIDENT_SCORE, score_candidate
from bt5.vector.markers import is_marker, is_recombination_site
from conftest import make_cds

FILLER = "TTATTATTAT"  # stop-rich, so it contributes no spurious ORFs


def feature(kind: str, start: int, end: int, strand: int = 1, label: str = "", **q: str) -> Feature:
    quals = {"label": (label or kind,)}
    quals.update({k: (v,) for k, v in q.items()})
    return Feature(Interval(start, end, strand), kind, quals, f"{kind}{start}")  # type: ignore[arg-type]


def backbone_with(cds: str, at: int, *, total: int = 3000, strand: int = 1, **kw: object):  # type: ignore[no-untyped-def]
    body = (FILLER * (total // 10 + 1))[:total]
    placed = cds if strand == 1 else reverse_complement(cds)
    seq = body[:at] + placed + body[at + len(cds) :]
    return VectorBackbone(sequence=seq[:total], topology=Topology.CIRCULAR, name="synthetic", **kw)  # type: ignore[arg-type]


class TestFindOrfs:
    def test_finds_a_forward_orf(self) -> None:
        cds, _ = make_cds(120)
        bb = backbone_with(cds, 300)
        assert Interval(300, 300 + len(cds), 1) in find_orfs(bb, table_id=1)

    def test_finds_a_reverse_orf_at_the_right_coordinates(self) -> None:
        cds, _ = make_cds(120)
        bb = backbone_with(cds, 300, strand=-1)
        found = [o for o in find_orfs(bb, table_id=1) if o.strand == -1]
        assert Interval(300, 300 + len(cds), -1) in found

    def test_a_reverse_orf_slices_back_to_the_coding_sequence(self) -> None:
        cds, _ = make_cds(120)
        bb = backbone_with(cds, 300, strand=-1)
        orf = next(o for o in find_orfs(bb, table_id=1) if o.strand == -1 and o.start == 300)
        assert bb.slice(orf) == cds

    def test_keeps_only_the_longest_orf_per_stop(self) -> None:
        """A downstream in-frame ATG describes the same gene."""
        cds, _ = make_cds(120)
        bb = backbone_with(cds, 300)
        ends = [o.end for o in find_orfs(bb, table_id=1) if o.strand == 1]
        assert len(ends) == len(set(ends))

    def test_respects_the_minimum_length(self) -> None:
        cds, _ = make_cds(40)
        bb = backbone_with(cds, 300)
        assert not [o for o in find_orfs(bb, table_id=1, min_aa=100) if o.start == 300]
        assert [o for o in find_orfs(bb, table_id=1, min_aa=20) if o.start == 300]

    def test_finds_an_orf_spanning_the_origin(self) -> None:
        """A circular vector has no privileged position; the transgene may straddle 0."""
        cds, _ = make_cds(120)
        bb = backbone_with(cds, 300)
        rotated = bb.rotated(360)  # drags the ORF across the origin
        site = next(o for o in find_orfs(rotated, table_id=1) if o.length == len(cds))
        assert site.end > rotated.length, "the ORF must be reported as wrapping"
        assert rotated.slice(site) == cds


class TestScoring:
    def orf_backbone(self, **kw: object):  # type: ignore[no-untyped-def]
        cds, _ = make_cds(200)
        return backbone_with(cds, 500, **kw), Interval(500, 500 + len(cds), 1)

    def test_a_selection_marker_is_penalised_below_confidence(self) -> None:
        """The failure that matters: an empty backbone's longest ORF is its marker."""
        bb, orf = self.orf_backbone()
        bb = replace(bb, features=(feature("CDS", 520, 900, label="PuroR"),))
        score, reasons = score_candidate(bb, orf, table_id=1)
        assert score < CONFIDENT_SCORE
        assert any("selection marker" in r for r in reasons)

    def test_a_non_marker_cds_inside_is_a_boost(self) -> None:
        bb, orf = self.orf_backbone()
        bb = replace(bb, features=(feature("CDS", 520, 560, label="Myc"),))
        score, reasons = score_candidate(bb, orf, table_id=1)
        assert score > 0
        assert any("Myc" in r for r in reasons)

    def test_a_kozak_at_the_start_codon_is_a_boost(self) -> None:
        bb, orf = self.orf_backbone()
        plain, _ = score_candidate(bb, orf, table_id=1)
        marked = replace(bb, features=(feature("regulatory", 494, 503, label="Kozak sequence"),))
        boosted, reasons = score_candidate(marked, orf, table_id=1)
        assert boosted > plain
        assert any("Kozak" in r for r in reasons)

    def test_a_signal_peptide_at_the_start_codon_is_a_boost(self) -> None:
        bb, orf = self.orf_backbone()
        plain, _ = score_candidate(bb, orf, table_id=1)
        marked = replace(bb, features=(feature("sig_peptide", 500, 560, label="Ig-kappa leader"),))
        boosted, _ = score_candidate(marked, orf, table_id=1)
        assert boosted > plain

    def test_a_marker_promoter_does_not_count_as_this_transcripts_promoter(self) -> None:
        bb, orf = self.orf_backbone()
        amp = replace(bb, features=(feature("promoter", 200, 480, label="AmpR promoter"),))
        real = replace(bb, features=(feature("promoter", 200, 480, label="CMV promoter"),))
        assert score_candidate(amp, orf, table_id=1)[0] < score_candidate(real, orf, table_id=1)[0]

    def test_every_candidate_explains_itself(self) -> None:
        bb, orf = self.orf_backbone()
        _, reasons = score_candidate(bb, orf, table_id=1)
        assert reasons, "a score with no reasons is not reviewable"


class TestCloningSites:
    def backbone(self) -> VectorBackbone:
        return VectorBackbone(
            sequence=(FILLER * 300)[:3000],
            topology=Topology.CIRCULAR,
            features=(
                feature("promoter", 100, 119, label="T7 promoter"),
                feature("promoter", 300, 520, label="CMV promoter"),
                feature("CDS", 900, 1400, label="PuroR"),
            ),
            name="empty",
        )

    def test_a_phage_promoter_does_not_open_a_cloning_site(self) -> None:
        """T7 and T3 are ~19 bp and drive sequencing, not expression."""
        sites = cloning_sites(self.backbone())
        assert all("T7" not in s.reasons[0] for s in sites)

    def test_an_expression_promoter_does(self) -> None:
        sites = cloning_sites(self.backbone())
        assert sites
        assert sites[0].interval == Interval(520, 900, 1)
        assert "CMV promoter" in sites[0].reasons[0]

    def test_a_marker_promoter_does_not(self) -> None:
        bb = self.backbone()
        bb = replace(bb, features=(feature("promoter", 300, 520, label="AmpR promoter"),))
        assert cloning_sites(bb) == ()


class TestSuggest:
    def test_a_transgene_outranks_the_selection_marker(self) -> None:
        cds, _ = make_cds(300)
        bb = backbone_with(cds, 700)
        bb = replace(
            bb,
            features=(
                feature("promoter", 400, 690, label="CMV promoter"),
                feature("polyA_signal", 1700, 1800, label="SV40 poly(A)"),
                feature("CDS", 2200, 2900, label="AmpR"),
            ),
        )
        best = suggest_insertion_sites(bb, table_id=1)[0]
        assert best.interval == Interval(700, 700 + len(cds), 1)

    def test_position_alone_is_not_enough_to_be_confident(self) -> None:
        """Promoter, polyA and length are a guess; the app should still ask."""
        cds, _ = make_cds(300)
        bb = replace(
            backbone_with(cds, 700),
            features=(
                feature("promoter", 400, 690, label="CMV promoter"),
                feature("polyA_signal", 1700, 1800, label="SV40 poly(A)"),
            ),
        )
        assert not suggest_insertion_sites(bb, table_id=1)[0].confident

    def test_annotation_that_pins_the_start_codon_makes_it_confident(self) -> None:
        cds, _ = make_cds(300)
        bb = replace(
            backbone_with(cds, 700),
            features=(
                feature("promoter", 400, 690, label="CMV promoter"),
                feature("regulatory", 694, 703, label="Kozak sequence"),
                feature("sig_peptide", 700, 760, label="Ig-kappa leader"),
                feature("polyA_signal", 1700, 1800, label="SV40 poly(A)"),
            ),
        )
        best = suggest_insertion_sites(bb, table_id=1)[0]
        assert best.interval == Interval(700, 700 + len(cds), 1)
        assert best.confident

    def test_when_nothing_is_confident_the_cloning_site_is_offered(self) -> None:
        """An empty backbone has no transgene; saying so is the honest answer."""
        cds, _ = make_cds(200)
        bb = backbone_with(cds, 700)
        bb = replace(
            bb,
            features=(
                feature("promoter", 400, 590, label="CMV promoter"),
                feature("CDS", 700, 1300, label="PuroR"),
            ),
        )
        suggestions = suggest_insertion_sites(bb, table_id=1)
        assert not any(c.confident for c in suggestions)
        assert any(c.kind == "cloning_site" for c in suggestions)


class TestMarkerVocabulary:
    @pytest.mark.parametrize(
        "text", ["AmpR", "PuroR", "confers resistance to puromycin", "beta-lactamase"]
    )
    def test_markers_are_recognised(self, text: str) -> None:
        assert is_marker(text)

    @pytest.mark.parametrize("text", ["EGFP", "Ig-kappa leader", "WPRE", "Myc"])
    def test_ordinary_features_are_not(self, text: str) -> None:
        assert not is_marker(text)

    @pytest.mark.parametrize("text", ["attB1", "loxP", "FRT site", "Gateway attP"])
    def test_recombination_sites_are_recognised(self, text: str) -> None:
        assert is_recombination_site(text)
