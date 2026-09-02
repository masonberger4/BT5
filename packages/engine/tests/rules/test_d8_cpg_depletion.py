"""D8: three CpG metrics that must not become one slider."""

from __future__ import annotations

import pytest
from bt5.core.context import HostId, Modality
from bt5.core.registry import discover, get
from bt5.core.spec import Enforcement, Evidence
from bt5.core.types import Construct, Interval, Segment, SegmentKind, Topology
from bt5.rules.catalog.d8_cpg_depletion import (
    MAG_TLR9_GENERAL,
    MAG_TLR9_SPECIFIC,
    TLR9_CONSENSUS_ONLY,
    ZAP_MIN_CPG,
    CpGDepletion,
    obs_over_exp,
)
from conftest import context, slot

discover()


def whole(seq: str, topology: Topology = Topology.CIRCULAR) -> Construct:
    return Construct(
        sequence=seq,
        topology=topology,
        segments=(Segment(Interval(0, len(seq)), SegmentKind.DESIGNABLE_CDS, "cds"),),
    )


def run(rule: CpGDepletion, c: Construct, host: HostId = HostId.HUMAN) -> list:
    ev = rule.evaluate(c, context(slot(host=host)), None)  # type: ignore[arg-type]
    return list(ev.breaches)


def metrics(breaches: list) -> set[str]:
    return {str(b.detail["metric"]) for b in breaches}


def hexamers(breaches: list) -> list:
    """TLR9 hexamer findings only -- not the total-count reading, which shares the
    metric label but is a measurement rather than a finding."""
    return [
        b
        for b in breaches
        if b.detail.get("metric") == "tlr9" and b.detail.get("reading") != "total_cpg"
    ]


class TestObsOverExp:
    """brief.md:131: obs/exp = (N_CpG x L) / (N_C x N_G)."""

    def test_alternating_cg(self) -> None:
        # CGCGCG: 3 C, 3 G, 3 CG, L=6 -> (3*6)/(3*3) = 2.0
        assert obs_over_exp("CGCGCG") == pytest.approx(2.0)

    def test_no_c_or_no_g_is_zero_not_a_division_error(self) -> None:
        """Undefined, and a window with no C cannot be an island under any
        definition -- so zero rather than an exception or a NaN."""
        assert obs_over_exp("AAAA") == 0.0
        assert obs_over_exp("CCCC") == 0.0
        assert obs_over_exp("GGGG") == 0.0

    def test_separated_c_and_g_score_below_one(self) -> None:
        assert obs_over_exp("CATCATGATGAT") < 1.0


class TestSeparateMetrics:
    """brief.md:128: three separate, separately-toggleable metrics -- do not collapse
    into one slider. Each toggle has to be independently effective or the instruction
    has not been honoured."""

    DENSE = "CGATCGTTACGCGATCGCTACG" * 20

    def test_each_metric_can_be_switched_off_alone(self) -> None:
        dense = whole(self.DENSE)
        assert "zap" not in metrics(run(CpGDepletion(zap=False), dense))
        assert "methylation" not in metrics(run(CpGDepletion(methylation=False), dense))
        assert "tlr9" not in metrics(run(CpGDepletion(tlr9=False), dense))

    def test_all_three_off_reports_nothing(self) -> None:
        c = whole(self.DENSE)
        assert run(CpGDepletion(tlr9=False, zap=False, methylation=False), c) == []

    def test_a_cpg_free_construct_is_clean(self) -> None:
        c = whole("ATTAATTAATTAATTAATTA" * 20)
        assert run(CpGDepletion(), c) == []

    def test_findings_carry_which_metric_produced_them(self) -> None:
        """The report cannot present three mechanisms separately if the breaches do
        not say which one each came from."""
        breaches = run(CpGDepletion(), whole(self.DENSE))
        assert breaches
        assert all("metric" in b.detail for b in breaches)


class TestZap:
    """brief.md:130: >=14 CpGs at mean inter-CpG spacing <=14 nt, worst window only."""

    def test_the_worst_window_is_reported_once_not_every_window(self) -> None:
        """ "Report the worst window, never the global count." A 400 nt CpG-dense
        construct contains hundreds of qualifying 200 nt windows."""
        breaches = [
            b for b in run(CpGDepletion(), whole("CGATCG" * 70)) if b.detail["metric"] == "zap"
        ]
        assert len(breaches) == 1

    def test_it_reports_spacing_and_not_only_a_count(self) -> None:
        breaches = [
            b for b in run(CpGDepletion(), whole("CGATCG" * 70)) if b.detail["metric"] == "zap"
        ]
        assert "mean_spacing" in breaches[0].detail
        assert float(breaches[0].detail["mean_spacing"]) <= 14.0

    def test_widely_spaced_cpgs_are_not_restricted(self) -> None:
        """brief.md:130: spacing >=32 nt is explicitly NOT restricted. One CpG every
        40 nt clears the count over a 200 nt window only if spacing is ignored."""
        c = whole(("CG" + "ATATATATATATATATATAT" * 2) * 30)
        assert "zap" not in metrics(run(CpGDepletion(), c))

    def test_too_few_cpgs_does_not_fire_however_tight(self) -> None:
        spaced = "CGCGCGCGCGCG" + "AT" * 200  # 6 CpGs, well under the threshold
        assert ZAP_MIN_CPG > 6
        assert "zap" not in metrics(run(CpGDepletion(), whole(spaced)))


class TestMethylationIslands:
    def test_takai_jones_is_stricter_than_gardiner_garden(self) -> None:
        """A 300 bp island at ~55% GC satisfies Gardiner-Garden's 200 bp length but
        not Takai-Jones' 500 bp."""
        c = whole("GCGCGATCGCGATCGCGATCGCGATCGCGA" * 10 + "ATATATATAT" * 30)
        assert "methylation" in metrics(run(CpGDepletion(island_criterion="gardiner_garden"), c))
        assert "methylation" not in metrics(run(CpGDepletion(island_criterion="takai_jones"), c))

    def test_an_unknown_criterion_is_refused_rather_than_defaulted(self) -> None:
        with pytest.raises(ValueError, match="gardiner_garden or takai_jones"):
            CpGDepletion(island_criterion="illumina")

    def test_a_low_gc_region_is_not_an_island_however_many_cpgs(self) -> None:
        assert "methylation" not in metrics(run(CpGDepletion(), whole("CGATATATAT" * 40)))


class TestTlr9:
    """brief.md:129: RRCGYY, with GTCGTT/TTCGTT (human) and GACGTT (mouse)."""

    def test_a_named_human_hexamer_escalates_in_a_human_host(self) -> None:
        c = whole("AAAA" + "GTCGTT" + "AAAA" * 10)
        breaches = hexamers(run(CpGDepletion(), c, HostId.HUMAN))
        assert breaches
        assert breaches[0].magnitude == MAG_TLR9_SPECIFIC
        assert breaches[0].detail["species"] == "human"

    def test_the_same_hexamer_does_not_escalate_in_a_mouse_host(self) -> None:
        """A human-attributed hexamer in a mouse cassette is a weaker finding, and
        brief.md:129 attributes each hexamer to one species."""
        c = whole("AAAA" + "GTCGTT" + "AAAA" * 10)
        breaches = hexamers(run(CpGDepletion(), c, HostId.MOUSE))
        assert breaches
        assert breaches[0].magnitude == MAG_TLR9_GENERAL

    def test_an_unmapped_host_gets_no_species_escalation_rather_than_a_guess(self) -> None:
        """brief.md:129 names human and mouse. CHO is neither, and the rule does not
        invent an attribution for it."""
        c = whole("AAAA" + "GTCGTT" + "AAAA" * 10)
        breaches = hexamers(run(CpGDepletion(), c, HostId.CHO))
        assert breaches
        assert breaches[0].magnitude == MAG_TLR9_GENERAL

    def test_a_general_rrcgyy_hexamer_is_reported_unattributed(self) -> None:
        c = whole("AAAA" + "AACGCC" + "AAAA" * 10)
        breaches = hexamers(run(CpGDepletion(), c, HostId.HUMAN))
        assert breaches
        assert breaches[0].detail["species"] == "unattributed"

    def test_a_non_rrcgyy_cpg_is_not_a_hexamer_finding(self) -> None:
        c = whole("AAAA" + "CCCGAA" + "AAAA" * 10)
        assert hexamers(run(CpGDepletion(), c, HostId.HUMAN)) == []

    def test_the_named_human_hexamers_are_not_rrcgyy(self) -> None:
        """The finding this rule would silently get wrong. brief.md:129 lists RRCGYY
        and then "specifically" GTCGTT/TTCGTT/GACGTT -- but position 2 of GTCGTT and
        position 1 of TTCGTT are T, and R is A or G. Reading the named hexamers as
        examples OF the consensus, and scanning only the consensus, produces a rule
        that never fires on either human motif."""
        assert not TLR9_CONSENSUS_ONLY.fullmatch("GTCGTT")
        assert not TLR9_CONSENSUS_ONLY.fullmatch("TTCGTT")
        assert TLR9_CONSENSUS_ONLY.fullmatch("GACGTT"), "the mouse motif does match"

    def test_the_total_cpg_count_is_reported_as_a_measurement(self) -> None:
        """brief.md:129 makes metric (a) "total CpG count + stimulatory hexamers", so
        the count is half the metric. Magnitude 0 and unfixable: a measurement, not a
        finding, and never a solver target."""
        c = whole("AAAA" + "CCCGAA" + "AAAA" * 10)
        totals = [b for b in run(CpGDepletion(), c) if b.detail.get("reading") == "total_cpg"]
        assert len(totals) == 1
        assert totals[0].magnitude == 0.0
        assert not totals[0].fixable_by_codon_choice
        assert totals[0].detail["cpg_total"] == 1.0

    def test_no_count_is_reported_when_there_are_no_cpgs(self) -> None:
        c = whole("ATTAATTAATTAATTAATTA" * 20)
        assert [b for b in run(CpGDepletion(), c) if b.detail.get("reading") == "total_cpg"] == []


class TestCircular:
    def test_a_cpg_across_the_origin_is_counted(self) -> None:
        """The construct ends C and begins G, so the only CpG in it exists solely
        because the molecule is circular. A rule evaluating a bare string misses it
        entirely -- which is why rules take a Construct (CLAUDE.md 3.3)."""
        seq = "G" + "ATATATATAT" * 4 + "C"
        assert "CG" not in seq, "the dinucleotide exists only across the origin"

        circular = _totals(run(CpGDepletion(), whole(seq, Topology.CIRCULAR)))
        linear = _totals(run(CpGDepletion(), whole(seq, Topology.LINEAR)))

        assert len(circular) == 1
        assert circular[0].detail["cpg_total"] == 1.0
        assert linear == [], "the same sequence, linear, has no CpG at all"


def _totals(breaches: list) -> list:
    return [b for b in breaches if b.detail.get("reading") == "total_cpg"]


class TestSpecShape:
    def test_it_is_soft_and_never_refuses(self) -> None:
        """brief.md:128's header carries no H/S marker, unlike D3 and D4."""
        assert CpGDepletion.enforcement is Enforcement.SOFT
        assert CpGDepletion().enforcement_for(slot()) is Enforcement.SOFT
        assert CpGDepletion.default_weight > 0.0

    def test_the_evidence_badge_is_contested(self) -> None:
        """The primary source for the ZAP arm reports that inhibition did NOT
        correlate with CpG number. An EVIDENCE_BACKED badge here would mark a
        contested quantity settled."""
        assert CpGDepletion.evidence is Evidence.CONTESTED

    def test_it_carries_a_citation_that_refutes(self) -> None:
        """Citation.sign exists for rules resting on sources with opposite signs; ZAP
        CpG is one of the three the type's own docstring names."""
        assert any(c.sign == "refutes" for c in CpGDepletion.citations)
        assert any(c.sign == "supports" for c in CpGDepletion.citations)

    def test_it_declares_the_gc_conflict(self) -> None:
        """brief.md:133: full depletion forces AGA/AGG for Arg and can drop GC below
        vendor floors."""
        assert "e2_gc_band" in CpGDepletion.conflicts_with

    def test_it_is_not_a_lattice_rule(self) -> None:
        """Forbidding CG outright would make most of the Arg/Ala/Pro/Ser/Thr codon
        space unreachable to express a soft preference."""
        assert CpGDepletion().lattice_terms(context()) is None  # type: ignore[arg-type]

    def test_a_soft_rule_explains_its_weight(self) -> None:
        assert CpGDepletion.weight_provenance.strip()

    def test_it_is_registered_under_its_brief_row(self) -> None:
        assert get("d8_cpg_depletion").brief_ref == "2.D8"

    @pytest.mark.parametrize(
        "modality", [Modality.LENTIVIRAL, Modality.AAV, Modality.PLASMID_TRANSIENT]
    )
    def test_it_applies_across_modalities(self, modality: Modality) -> None:
        assert CpGDepletion().gate(slot(modality=modality))
