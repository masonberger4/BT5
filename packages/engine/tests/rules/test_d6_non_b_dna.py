"""D6: G-quadruplexes and telomere repeats, and the host-dependent enforcement.

The interesting property of this rule is that the evidence points in opposite
directions by host and both directions are well measured, so the tests are
mostly about the rule saying so rather than picking a side.
"""

from __future__ import annotations

import pytest
from bt5.core.context import Modality
from bt5.core.registry import discover, get
from bt5.core.spec import Enforcement, Evidence
from bt5.core.types import Construct, Interval, Segment, SegmentKind, Topology
from bt5.rules.catalog.d6_non_b_dna import (
    G4HUNTER_SEVERE,
    NonBDna,
    g4hunter_score,
    peak_g4hunter,
)
from conftest import construct, context, slot

discover()

#: A canonical four-tract G4: four G3 runs with short loops.
G4 = "GGGTTAGGGTTAGGGTTAGGG"
#: Its reverse complement -- a G4 on the OTHER strand, invisible to a forward scan.
C4 = "CCCTAACCCTAACCCTAACCC"
PAD = "ACTACTACTACTACTACTACT"


def hits(rule: NonBDna, c, ctx=None):
    return rule.evaluate(c, ctx or context(), None).breaches


class TestG4Hunter:
    def test_a_g_run_scores_positive_and_a_c_run_negative(self) -> None:
        """The sign says which strand, which is why callers take the absolute
        value rather than clamping at zero."""
        assert g4hunter_score("GGGG") == 4.0
        assert g4hunter_score("CCCC") == -4.0

    def test_a_run_scores_by_its_own_length_capped_at_four(self) -> None:
        assert g4hunter_score("G") == 1.0
        assert g4hunter_score("GG") == 2.0
        assert g4hunter_score("GGGGGGG") == 4.0, "capped at 4"

    def test_mixed_sequence_averages(self) -> None:
        assert g4hunter_score("GATC") == pytest.approx(0.0)

    def test_empty_is_zero_rather_than_an_error(self) -> None:
        assert g4hunter_score("") == 0.0

    def test_peak_takes_the_strongest_window(self) -> None:
        seq = "ACTACTACTACTACTACTACT" + "GGGG" * 8
        assert peak_g4hunter(seq) > peak_g4hunter("ACTACTACT" * 10)

    def test_a_sequence_shorter_than_the_window_scores_zero(self) -> None:
        """Not a partial window: a score over 8 nt is not comparable to one over
        25, and returning it anyway would put both on the same axis."""
        assert peak_g4hunter("GGGG") == 0.0


class TestDetection:
    def test_finds_a_canonical_g4(self) -> None:
        found = hits(NonBDna(), construct("ATG" + PAD + G4 + PAD + "TAA"))
        assert any(b.detail["motif"] == "G4" for b in found)

    def test_finds_a_g4_on_the_reverse_strand(self) -> None:
        """A C-rich forward strand IS a G4 -- on the strand the scan would miss."""
        found = hits(NonBDna(), construct("ATG" + PAD + C4 + PAD + "TAA"))
        g4s = [b for b in found if b.detail["motif"] == "G4"]
        assert g4s
        assert any(b.interval.strand == -1 for b in g4s)

    def test_the_reverse_hit_is_reported_at_its_forward_position(self) -> None:
        c = construct("ATG" + PAD + C4 + PAD + "TAA")
        breach = next(b for b in hits(NonBDna(), c) if b.detail["motif"] == "G4")
        assert c.sequence[breach.interval.start : breach.interval.end] == C4

    def test_ordinary_sequence_is_clean(self) -> None:
        assert not hits(NonBDna(), construct("ATG" + PAD * 3 + "TAA"))

    def test_finds_telomere_repeats(self) -> None:
        found = hits(NonBDna(), construct("ATG" + PAD + "TTAGGG" * 3 + PAD + "TAA"))
        assert any(b.detail["class"] == "telomere" for b in found)

    def test_a_single_telomere_unit_is_not_a_repeat(self) -> None:
        found = hits(NonBDna(), construct("ATG" + PAD + "TTAGGG" + PAD + "TAA"))
        assert not any(b.detail["class"] == "telomere" for b in found)

    def test_telomeres_can_be_switched_off(self) -> None:
        c = construct("ATG" + PAD + "TTTAGGG" * 3 + PAD + "TAA")
        assert any(b.detail["class"] == "telomere" for b in hits(NonBDna(telomeres=True), c))
        assert not any(b.detail["class"] == "telomere" for b in hits(NonBDna(telomeres=False), c))

    def test_a_g4_spanning_the_origin_is_found(self) -> None:
        """The plasmid does not know where the file's first character is."""
        seq = G4[10:] + PAD + G4[:10]
        c = Construct(
            seq,
            Topology.CIRCULAR,
            (Segment(Interval(0, len(seq)), SegmentKind.DESIGNABLE_CDS, "cds"),),
        )
        found = [b for b in hits(NonBDna(), c) if b.detail["motif"] == "G4"]
        assert found, "a G4 assembled across position 0 is still a G4"

    def test_the_same_sequence_linear_is_clean(self) -> None:
        seq = G4[10:] + PAD + G4[:10]
        c = Construct(
            seq,
            Topology.LINEAR,
            (Segment(Interval(0, len(seq)), SegmentKind.DESIGNABLE_CDS, "cds"),),
        )
        assert not [b for b in hits(NonBDna(), c) if b.detail["motif"] == "G4"]


class TestSeverity:
    def test_a_strong_g4_outranks_a_weak_match(self) -> None:
        strong = hits(NonBDna(), construct("ATG" + PAD + "GGGG" * 8 + PAD + "TAA"))
        assert any(b.magnitude >= 2.0 for b in strong), "G4Hunter severe"

    def test_a_regex_match_scoring_low_is_reported_not_dropped(self) -> None:
        """The regex and G4Hunter disagree often and neither is a census, so a
        low-scoring match is reported at low magnitude rather than discarded."""
        c = construct("ATG" + PAD + G4 + PAD + "TAA")
        found = [b for b in hits(NonBDna(), c) if b.detail["motif"] == "G4"]
        assert found
        assert all(b.magnitude > 0 for b in found)

    def test_the_g4hunter_score_travels_with_the_breach(self) -> None:
        c = construct("ATG" + PAD + "GGGG" * 8 + PAD + "TAA")
        breach = next(b for b in hits(NonBDna(), c) if b.detail["motif"] == "G4")
        assert float(breach.detail["g4hunter"]) >= G4HUNTER_SEVERE

    def test_a_backbone_g4_is_reported_but_unfixable(self) -> None:
        c = construct("ATGCTGTAA", PAD + G4 + PAD)
        found = [b for b in hits(NonBDna(), c) if b.detail["motif"] == "G4"]
        assert found
        assert not any(b.fixable_by_codon_choice for b in found)


class TestEnforcement:
    def test_hard_in_e_coli_where_g4s_are_measured_to_fold(self) -> None:
        assert NonBDna().enforcement_for(slot(modality=Modality.BACTERIAL_EXPRESSION)).is_hard

    def test_soft_in_mammalian_contexts_where_they_are_measured_not_to(self) -> None:
        """rG4s are globally unfolded in mammalian cells: median folding score
        0.06. The sequence is there and the structure mostly is not."""
        rule = NonBDna()
        for modality in (Modality.LENTIVIRAL, Modality.AAV, Modality.PLASMID_TRANSIENT):
            assert rule.enforcement_for(slot(modality=modality)) is Enforcement.SOFT

    def test_the_citations_carry_opposite_signs(self) -> None:
        """One URL would make the badge dishonest on precisely the rule where
        the disagreement is the point."""
        signs = {c.sign for c in NonBDna.citations}
        assert "supports" in signs
        assert "refutes" in signs

    def test_it_is_evidence_backed_not_contested(self) -> None:
        """Both directions are well measured; they simply apply to different
        hosts. That is a job for enforcement_for, not for a weaker badge."""
        assert NonBDna.evidence is Evidence.EVIDENCE_BACKED

    def test_the_default_weight_is_modest(self) -> None:
        """Every base spent removing a motif that does not fold is a base not
        spent on repeats, which are the measured top predictor of failure."""
        assert 0.0 < NonBDna.default_weight < 0.5
        assert NonBDna.weight_provenance.strip()

    def test_it_is_not_a_lattice_rule(self) -> None:
        """A G4 is variable-length with 1-7 nt loops and unbounded G-runs, so it
        is not a finite motif set the automaton can make unreachable."""
        assert NonBDna().lattice_terms(None) is None

    def test_it_does_not_apply_to_ivt_mrna(self) -> None:
        assert not NonBDna().gate(slot(modality=Modality.IVT_MRNA))


def test_a_zero_flag_threshold_is_refused() -> None:
    with pytest.raises(ValueError, match="must be positive"):
        NonBDna(g4hunter_flag=0.0)


def test_it_is_registered_under_its_brief_row() -> None:
    assert get("d6_non_b_dna") is NonBDna
    assert NonBDna.brief_ref == "2.D6"
