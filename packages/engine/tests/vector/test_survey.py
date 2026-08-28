"""Reading a plasmid BT5 was not asked to redesign.

The case this exists for is a construct with nothing to back-translate -- an
sgRNA under a U6 promoter, an empty backbone -- where the useful output is
still "here is what will make this hard to propagate". So the tests that matter
are that it produces a construct with NO designable region, that it says so
rather than inventing an insertion site, and that its findings agree with the
design path rather than drifting from it.
"""

from __future__ import annotations

import numpy as np
import pytest
from bt5.core.types import (
    Feature,
    Interval,
    SegmentKind,
    Topology,
    reverse_complement,
)
from bt5.vector import VectorBackbone
from bt5.vector.findings import repeat_liabilities
from bt5.vector.survey import construct_from_backbone, survey


def dna(n: int, *, seed: int = 5) -> str:
    rng = np.random.default_rng(seed)
    return "".join("ACGT"[i] for i in rng.integers(0, 4, n))


def orf(n_codons: int, *, seed: int = 13) -> str:
    """A real ATG..stop reading frame with no internal stop codon."""
    from Bio.Data import CodonTable

    table = CodonTable.unambiguous_dna_by_id[1]
    sense = sorted(table.forward_table)
    rng = np.random.default_rng(seed)
    picks = rng.integers(0, len(sense), n_codons - 2)
    return "ATG" + "".join(sense[i] for i in picks) + "TAA"


def feature(kind: str, start: int, end: int, strand: int = 1, label: str = "") -> Feature:
    return Feature(
        interval=Interval(start, end, strand),  # type: ignore[arg-type]
        kind=kind,
        qualifiers={"label": (label or kind,)},
        uid=f"{kind}{start}",
    )


def vector(
    sequence: str, *features: Feature, circular: bool = True, name: str = "plasmid"
) -> VectorBackbone:
    return VectorBackbone(
        sequence=sequence,
        topology=Topology.CIRCULAR if circular else Topology.LINEAR,
        features=features,
        name=name,
    )


class TestNothingIsDesignable:
    def test_the_construct_has_no_editable_region(self) -> None:
        """Not a degenerate case to work around: every base came from the user,
        so every finding against it is unfixable by codon choice."""
        c = construct_from_backbone(vector(dna(1200)))
        assert c.editable == ()
        assert not c.is_editable(Interval(10, 40))

    def test_no_translation_unit_is_invented(self) -> None:
        """A translation unit carries a genetic code, and guessing one is the
        failure mode this project refuses outright."""
        c = construct_from_backbone(vector(dna(1200), feature("CDS", 100, 400)))
        assert c.translation_units == ()

    def test_the_sequence_and_features_are_carried_through_unchanged(self) -> None:
        seq = dna(1200)
        bb = vector(seq, feature("promoter", 10, 60), feature("rep_origin", 700, 900))
        c = construct_from_backbone(bb)
        assert c.sequence == seq
        assert len(c.features) == 2

    def test_segments_tile_the_construct_exactly_once(self) -> None:
        bb = vector(dna(1200), feature("LTR", 100, 300, 1, "5' LTR"))
        c = construct_from_backbone(bb)
        covered = sorted((s.interval.start, s.interval.end) for s in c.segments)
        assert covered[0][0] == 0
        assert covered[-1][1] == c.length
        for (_, end), (start, _) in zip(covered, covered[1:], strict=False):
            assert end == start, "segments must abut without gaps or overlaps"

    def test_topology_is_carried_not_assumed(self) -> None:
        assert construct_from_backbone(vector(dna(600), circular=False)).topology is Topology.LINEAR


class TestExemptRegions:
    def test_an_ltr_becomes_a_whitelisted_repeat(self) -> None:
        bb = vector(dna(1200), feature("LTR", 100, 300, 1, "5' LTR"))
        kinds = {s.kind for s in construct_from_backbone(bb).segments}
        assert SegmentKind.WHITELISTED_REPEAT in kinds

    def test_an_itr_annotated_as_misc_feature_is_still_exempt(self) -> None:
        """AAV ITRs are routinely deposited as misc_feature; the label is the
        only signal, and it is the one the design path uses too."""
        bb = vector(dna(1200), feature("misc_feature", 100, 250, 1, "AAV2 ITR"))
        exempt = [s for s in construct_from_backbone(bb).segments if s.exempt_from_scanning]
        assert exempt
        assert exempt[0].label == "AAV2 ITR"

    def test_a_palindrome_inside_an_exempt_region_is_not_a_finding(self) -> None:
        """An AAV ITR IS a palindrome. Reporting it every run trains the reader
        to skip the section."""
        arm = dna(40, seed=17)
        itr = arm + reverse_complement(arm)
        seq = dna(300) + itr + dna(500, seed=19)
        bare = survey(vector(seq))
        marked = survey(
            vector(seq, feature("misc_feature", 290, 300 + len(itr) + 10, 1, "AAV2 ITR"))
        )
        assert any("inverted repeat" in n.summary for n in bare.notes)
        assert not any("inverted repeat" in n.summary for n in marked.notes)

    def test_the_report_says_when_nothing_was_exempted(self) -> None:
        """Which repeats are 'accepted' depends entirely on what the depositor
        annotated -- one AAV deposit in hand annotates its ITRs and another does
        not -- so a reader needs to know before judging a palindrome finding."""
        assert "no region is annotated as an ITR" in survey(vector(dna(900))).to_comment()

    def test_the_report_names_what_was_exempted(self) -> None:
        bb = vector(dna(1200), feature("LTR", 100, 300, 1, "5' LTR"))
        assert "5' LTR" in survey(bb).to_comment()


class TestFindingsAgreeWithTheDesignPath:
    def test_the_same_repeat_is_reported_either_way(self) -> None:
        """A survey that reimplemented the checks would drift from the design
        path exactly where it is hardest to notice."""
        unit = dna(40, seed=23)
        seq = dna(300) + unit + dna(60, seed=29) + unit + dna(300, seed=31)
        bb = vector(seq)
        direct = repeat_liabilities(construct_from_backbone(bb))
        assert direct
        assert [n.summary for n in survey(bb).notes] == [n.summary for n in direct]

    def test_a_hairpin_is_reported_with_its_own_action(self) -> None:
        """The only reason a user cares which kind of repeat they have is that
        the action differs; re-running the optimizer will not remove a stem."""
        arm = dna(40, seed=37)
        seq = dna(300) + arm + dna(10, seed=41) + reverse_complement(arm) + dna(300, seed=43)
        stems = [n for n in survey(vector(seq)).notes if "inverted repeat" in n.summary]
        assert stems
        assert all("sbcC" in n.action for n in stems)
        assert all("will not remove it" in n.action for n in stems)

    def test_a_clean_plasmid_reports_nothing(self) -> None:
        s = survey(vector(dna(3000, seed=47)))
        assert s.clean
        assert "nothing further to report" in s.to_comment()


class TestNoCodingSequence:
    def test_a_non_coding_plasmid_says_so_instead_of_guessing(self) -> None:
        """The pU6/sgRNA case: the only ORF is the selection marker, and
        offering to redesign it would be worse than offering nothing."""
        marker = orf(200, seed=53)
        seq = dna(400) + marker + dna(500, seed=59)
        scaffold = 400 + len(marker) + 100
        bb = vector(
            seq,
            feature("CDS", 400, 400 + len(marker), 1, "AmpR"),
            feature("misc_RNA", scaffold, scaffold + 76, 1, "gRNA scaffold"),
        )
        s = survey(bb)
        assert s.candidates, "the marker ORF is found -- and then rejected"
        assert not s.has_confident_site
        assert "no confident coding sequence" in s.to_comment()

    def test_a_plasmid_with_no_reading_frame_at_all_says_that_instead(self) -> None:
        """A different sentence, because it is a different situation: nothing
        was found, rather than something found and judged not worth redesigning."""
        s = survey(vector("AT" * 400))
        assert s.candidates == ()
        assert "no open reading frame found" in s.to_comment()

    def test_an_empty_backbone_still_gets_a_report(self) -> None:
        s = survey(vector(dna(2000, seed=61)))
        assert not s.has_confident_site
        assert "no sequence was designed" in s.to_comment()

    def test_the_survey_never_claims_to_have_generated_anything(self) -> None:
        comment = survey(vector(dna(2000, seed=67))).to_comment()
        generated = comment.split("== GENERATED BY BT5 ==")[1].split("== NOTED BY BT5 ==")[0]
        assert "nothing" in generated

    def test_no_note_predicts_an_outcome(self) -> None:
        banned = ("predict", "titer", "yield", "fold-improvement", "expression level")
        unit = dna(40, seed=71)
        seq = dna(300) + unit + dna(60, seed=73) + unit + dna(300, seed=79)
        for note in survey(vector(seq)).notes:
            assert not any(w in note.summary.lower() for w in banned)


class TestOriginSpanning:
    def test_a_finding_across_the_origin_is_reported_in_one_frame(self) -> None:
        arm = dna(40, seed=83)
        hairpin = arm + dna(10, seed=89) + reverse_complement(arm)
        cut = 45
        seq = hairpin[cut:] + dna(600, seed=97) + hairpin[:cut]
        stems = [n for n in survey(vector(seq)).notes if "inverted repeat" in n.summary]
        assert stems, "a stem-loop across the origin is still a stem-loop"

    def test_an_exempt_feature_wrapping_the_origin_is_still_exempt(self) -> None:
        seq = dna(1200)
        bb = vector(seq, feature("LTR", 1100, 1300, 1, "5' LTR"))
        exempt = [s for s in construct_from_backbone(bb).segments if s.exempt_from_scanning]
        assert exempt, "a wrapping feature must not lose its exemption"
        assert sum(s.interval.length for s in exempt) == 200


@pytest.mark.parametrize("length", [200, 1000, 5000])
def test_any_plasmid_surveys_without_raising(length: int) -> None:
    """Adversarial sizes: the survey path must degrade, never crash."""
    s = survey(vector(dna(length, seed=101)))
    assert s.construct.length == length
    assert s.to_comment()
