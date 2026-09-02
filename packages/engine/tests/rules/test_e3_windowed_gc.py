"""E3: windowed GC, and the published thresholds 18 measured probes refuted.

The tests that matter most here are the negative ones. A rule that reported windowed GC
would look correct while quietly refusing constructs both vendors manufacture; what makes
this rule right is what it declines to enforce.
"""

from __future__ import annotations

import pytest
from bt5.core.context import Modality
from bt5.core.registry import discover, get
from bt5.core.spec import Direction, Enforcement, Evidence
from bt5.core.types import Construct, Interval, Segment, SegmentKind, Topology
from bt5.rules.catalog.e3_windowed_gc import (
    GENSCRIPT_HI,
    GENSCRIPT_LO,
    IDT_FLOOR,
    IDT_WINDOW,
    MAG_FLOOR,
    TWIST_TRIGGER_HI,
    TWIST_TRIGGER_LO,
    TWIST_WINDOW,
    WindowedGC,
)
from conftest import context, slot

discover()

#: 50% GC, comfortably inside every band in the rule.
CLEAN = "ACGT" * 150


def whole(seq: str, topology: Topology = Topology.CIRCULAR) -> Construct:
    return Construct(
        sequence=seq,
        topology=topology,
        segments=(Segment(Interval(0, len(seq)), SegmentKind.DESIGNABLE_CDS, "cds"),),
    )


def ev(rule: WindowedGC, seq: str, topology: Topology = Topology.CIRCULAR):
    return rule.evaluate(whole(seq, topology), context(), None)  # type: ignore[arg-type]


def tiers(result) -> list[str]:
    return [str(b.detail["tier"]) for b in result.breaches]


class TestNothingWindowedIsEverHard:
    """The load-bearing property. docs/design/vendor-gc-calibration.md measured 18
    probes at two vendors and found no windowed GC threshold that refused anything --
    including IDT's own stated floor, which its LOC probes tripped while staying green.
    """

    def test_the_rule_is_soft(self) -> None:
        assert WindowedGC.enforcement is Enforcement.SOFT
        assert WindowedGC().enforcement_for(slot()) is Enforcement.SOFT

    def test_a_soft_rule_carries_a_weight_and_explains_it(self) -> None:
        assert WindowedGC.default_weight > 0.0
        assert WindowedGC.weight_provenance.strip()

    @pytest.mark.parametrize("block", ["AT" * 25, "GC" * 25])
    def test_an_extreme_fifty_bp_window_is_only_ever_a_note(self, block: str) -> None:
        """The calibration's LOC ladder: one 50 bp window from 4% to 96% GC on a 50%
        background, 16 accepting verdicts out of 16. brief.md:140 would hard-fail every
        one of them, so the 50 bp reading must never rise above the note tier."""
        result = ev(WindowedGC(), CLEAN + block + CLEAN)
        fifty = [b for b in result.breaches if b.detail["window"] == float(TWIST_WINDOW)]
        assert fifty, "the trigger is still reported -- silence would lose the provenance"
        assert all(b.detail["tier"] == "note" for b in fifty)
        assert all(b.magnitude < MAG_FLOOR for b in fifty)

    def test_a_gc_rich_construct_produces_no_floor_finding(self) -> None:
        """The asymmetry the calibration measured: GLB_gc80 produced no window finding
        at all despite being 80% GC throughout. There is no windowed ceiling."""
        assert "floor" not in tiers(ev(WindowedGC(), "GC" * 300))


class TestIdtFloor:
    """The one windowed rule a vendor was observed to score against at all."""

    def test_a_low_gc_region_is_reported_at_the_floor_tier(self) -> None:
        result = ev(WindowedGC(), "AT" * 100 + CLEAN)
        assert "floor" in tiers(result)
        assert not result.passes, "the rule's own verdict, not a refusal"

    def test_the_finding_names_the_measured_geometry(self) -> None:
        breach = next(
            b for b in ev(WindowedGC(), "AT" * 100 + CLEAN).breaches if b.detail["tier"] == "floor"
        )
        assert breach.detail["window"] == float(IDT_WINDOW)
        assert f"{IDT_FLOOR:.0%}" in breach.message

    def test_the_floor_is_a_gradient_not_a_flag(self) -> None:
        """A partial improvement has to be visible in the score."""
        mild = ev(WindowedGC(), "ATATATATAC" * 10 + CLEAN)
        severe = ev(WindowedGC(), "AT" * 100 + CLEAN)
        assert severe.raw_score >= mild.raw_score

    def test_the_floor_is_configurable_and_validated(self) -> None:
        assert "floor" not in tiers(ev(WindowedGC(floor=0.0), "AT" * 100 + CLEAN))
        with pytest.raises(ValueError, match="GC fraction"):
            WindowedGC(floor=1.5)

    def test_a_balanced_construct_has_no_findings_at_all(self) -> None:
        result = ev(WindowedGC(), CLEAN)
        assert result.breaches == ()
        assert result.passes


class TestReportedNotEnforced:
    def test_the_twist_trigger_is_reported_as_a_note(self) -> None:
        result = ev(WindowedGC(), CLEAN + "GC" * 30 + CLEAN)
        notes = [b for b in result.breaches if b.detail["tier"] == "note"]
        assert notes
        assert all(b.magnitude < MAG_FLOOR for b in notes)

    def test_the_note_says_it_is_not_enforced(self) -> None:
        result = ev(WindowedGC(), CLEAN + "GC" * 30 + CLEAN)
        note = next(b for b in result.breaches if b.detail["tier"] == "note")
        assert "NOT ENFORCED" in note.message
        assert "16 verdicts of 16" in note.message

    def test_the_fifty_bp_report_can_be_switched_off(self) -> None:
        result = ev(WindowedGC(report_fifty_bp=False), CLEAN + "GC" * 30 + CLEAN)
        assert all(b.detail["window"] != float(TWIST_WINDOW) for b in result.breaches)

    def test_the_genscript_band_is_reported_at_its_own_tier(self) -> None:
        result = ev(WindowedGC(), "GC" * 100 + CLEAN)
        assert "band" in tiers(result)

    def test_the_published_trigger_bounds_are_carried_verbatim(self) -> None:
        """brief.md:140. The numbers are reported even though they do not gate --
        removing them would lose the provenance that makes the refusal auditable."""
        assert (TWIST_TRIGGER_LO, TWIST_TRIGGER_HI) == (0.10, 0.90)
        assert (GENSCRIPT_LO, GENSCRIPT_HI) == (0.25, 0.65)


class TestRegions:
    def test_one_contiguous_low_stretch_is_one_floor_breach(self) -> None:
        result = ev(WindowedGC(), "AT" * 100 + CLEAN)
        assert len([b for b in result.breaches if b.detail["tier"] == "floor"]) == 1

    def test_two_separated_stretches_are_two(self) -> None:
        result = ev(WindowedGC(), "AT" * 100 + CLEAN + "AT" * 100 + CLEAN)
        assert len([b for b in result.breaches if b.detail["tier"] == "floor"]) == 2

    def test_a_region_crossing_the_origin_is_one_breach(self) -> None:
        """merge_regions is shared with f5_at_window precisely so this case has one
        implementation: a stretch spanning the origin arrives as two runs."""
        half = 100
        seq = "AT" * half + CLEAN + "AT" * half
        result = ev(WindowedGC(), seq, Topology.CIRCULAR)
        assert len([b for b in result.breaches if b.detail["tier"] == "floor"]) == 1


class TestSpecShape:
    def test_the_evidence_badge_is_contested_despite_the_row_grade(self) -> None:
        """brief.md:140 grades the row A, but that grade attaches to published numbers
        the calibration contradicts. EVIDENCE_BACKED would let the refuted half of the
        row wear the badge."""
        assert WindowedGC.evidence is Evidence.CONTESTED

    def test_it_carries_the_citation_that_refutes_its_own_row(self) -> None:
        refuting = [c for c in WindowedGC.citations if c.sign == "refutes"]
        assert len(refuting) == 1
        assert "16 verdicts of 16" in refuting[0].label

    def test_it_is_a_band_rule_with_a_non_inverted_band(self) -> None:
        assert WindowedGC.direction is Direction.BAND
        assert WindowedGC.band == (GENSCRIPT_LO, GENSCRIPT_HI)
        assert WindowedGC.band[0] < WindowedGC.band[1]

    def test_it_steers_less_than_the_rule_that_owns_global_gc(self) -> None:
        """Global GC is what gates and e2_gc_band owns it. Three full-strength GC
        steering terms would be one preference counted three times."""
        from bt5.rules.catalog.e2_gc_band import GCBand

        assert WindowedGC.steering_weight < GCBand.steering_weight

    def test_it_declares_the_overlapping_gc_rules(self) -> None:
        assert set(WindowedGC.conflicts_with) == {"e2_gc_band", "f5_at_window"}

    def test_it_is_not_a_lattice_rule(self) -> None:
        assert WindowedGC().lattice_terms(context()) is None  # type: ignore[arg-type]

    def test_it_is_registered_under_its_brief_row(self) -> None:
        assert get("e3_windowed_gc").brief_ref == "2.E3"

    @pytest.mark.parametrize(
        "modality", [Modality.LENTIVIRAL, Modality.PLASMID_TRANSIENT, Modality.IVT_MRNA]
    )
    def test_it_applies_everywhere(self, modality: Modality) -> None:
        assert WindowedGC().gate(slot(modality=modality))
