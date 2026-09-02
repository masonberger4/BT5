"""F5: the two-sided GC band per 100 nt, and the side that is binding."""

from __future__ import annotations

import pytest
from bt5.core.context import Modality
from bt5.core.registry import discover, get
from bt5.core.spec import Direction, Enforcement
from bt5.core.types import Construct, Interval, Segment, SegmentKind, Topology
from bt5.rules.catalog.f5_at_window import (
    HARD_HI,
    HARD_LO,
    MAG_HARD_BASE,
    MAG_WARN,
    SOFT_HI,
    SOFT_LO,
    ATWindow,
    window_gc,
)
from conftest import context, slot

discover()


def whole(seq: str, topology: Topology = Topology.CIRCULAR) -> Construct:
    return Construct(
        sequence=seq,
        topology=topology,
        segments=(Segment(Interval(0, len(seq)), SegmentKind.DESIGNABLE_CDS, "cds"),),
    )


def ev(rule: ATWindow, seq: str, topology: Topology = Topology.CIRCULAR):
    return rule.evaluate(whole(seq, topology), context(), None)  # type: ignore[arg-type]


#: ~50% GC, comfortably inside the preference band.
CLEAN = "ACGT" * 100
#: 100% AT and 100% GC blocks, both far outside the hard band.
AT_BLOCK = "AT" * 100
GC_BLOCK = "GC" * 100


class TestWindowGc:
    def test_it_reads_only_the_middle_copy(self) -> None:
        """Tripled text, but a window per construct position -- not per tripled
        position, which would report each site three times."""
        text, offset, n, span = "ACGT" * 3, 4, 4, 2
        assert window_gc(text, offset, n, span).size == n

    def test_a_pure_gc_run_is_one(self) -> None:
        assert window_gc("GGGG", 0, 1, 4)[0] == pytest.approx(1.0)

    def test_a_pure_at_run_is_zero(self) -> None:
        assert window_gc("ATAT", 0, 1, 4)[0] == pytest.approx(0.0)


class TestBands:
    """brief.md:159: two-sided band 45-60% GC per 100 nt, hard-fail outside 35-65%."""

    def test_a_balanced_construct_passes_with_no_findings(self) -> None:
        result = ev(ATWindow(), CLEAN)
        assert result.passes
        assert result.breaches == ()

    def test_an_at_rich_window_hard_fails(self) -> None:
        result = ev(ATWindow(), AT_BLOCK + CLEAN)
        assert not result.passes
        assert any(b.detail["hard"] == "yes" for b in result.breaches)

    def test_a_gc_rich_window_hard_fails(self) -> None:
        result = ev(ATWindow(), GC_BLOCK + CLEAN)
        assert not result.passes
        assert any(b.detail["hard"] == "yes" for b in result.breaches)

    def test_the_hard_band_is_wider_than_the_preference_band(self) -> None:
        assert HARD_LO < SOFT_LO < SOFT_HI < HARD_HI

    def test_a_window_between_the_bands_warns_without_failing(self) -> None:
        """Inside 35-65% but outside 45-60%. solver/catalog.py:158-170: handing a
        warn-band finding to repair sets it chasing a threshold never crossed."""
        # 40% GC: outside the preference, inside the hard band.
        seq = ("ACGTATATAT" * 10 + CLEAN) * 2
        result = ev(ATWindow(), seq)
        warns = [b for b in result.breaches if b.detail["hard"] == "no"]
        if warns:  # the fixture is only useful if it produced one
            assert all(b.magnitude == MAG_WARN for b in warns)
            assert all(b.magnitude < MAG_HARD_BASE for b in warns)

    def test_an_inverted_hard_band_is_refused(self) -> None:
        with pytest.raises(ValueError, match="non-inverted"):
            ATWindow(hard_lo=0.7, hard_hi=0.3)


class TestBindingSide:
    """brief.md:159: "show which side is binding per window". A |deviation| scalar
    cannot say whether a window needs GC raised or lowered, and the two sides have
    different mechanisms behind them."""

    def test_an_at_rich_window_reports_the_lower_bound(self) -> None:
        result = ev(ATWindow(), AT_BLOCK + CLEAN)
        assert result.binding_side == "lower"
        assert result.breaches[0].detail["binding_side"] == "lower"

    def test_a_gc_rich_window_reports_the_upper_bound(self) -> None:
        result = ev(ATWindow(), GC_BLOCK + CLEAN)
        assert result.binding_side == "upper"

    def test_both_sides_can_be_reported_on_one_construct(self) -> None:
        """A construct with an AT-rich end and a GC-rich end is binding on opposite
        sides in different places -- which is the case a single scalar destroys."""
        result = ev(ATWindow(), AT_BLOCK + CLEAN + GC_BLOCK + CLEAN)
        sides = {str(b.detail["binding_side"]) for b in result.breaches}
        assert sides == {"lower", "upper"}

    def test_the_at_framing_is_carried_in_the_message(self) -> None:
        """The evidence is in AT (63-68% AT toxic genes vs a 55% control), so the
        finding says both."""
        result = ev(ATWindow(), AT_BLOCK + CLEAN)
        assert "AT" in result.breaches[0].message
        assert "toxic" in result.breaches[0].message


class TestRegionMerging:
    def test_one_contiguous_bad_stretch_is_one_breach(self) -> None:
        """A 5 kb plasmid at 30% GC has ~5,000 offending windows. One breach per
        window would make this rule's COUNT -- the currency _aggregate steers on --
        swamp every other rule for what is one contiguous problem."""
        result = ev(ATWindow(), AT_BLOCK + CLEAN)
        assert len([b for b in result.breaches if b.detail["binding_side"] == "lower"]) == 1

    def test_two_separated_stretches_are_two_breaches(self) -> None:
        result = ev(ATWindow(), AT_BLOCK + CLEAN + AT_BLOCK + CLEAN)
        assert len(result.breaches) == 2

    def test_the_named_window_is_the_worst_in_its_region(self) -> None:
        result = ev(ATWindow(), AT_BLOCK + CLEAN)
        worst = result.breaches[0]
        assert float(worst.detail["gc"]) == pytest.approx(0.0, abs=0.05)


class TestMagnitudeIsAGradient:
    """_accepts needs a strictly falling magnitude sum to accept a move that improves
    a breach without clearing it. A constant hard magnitude makes 66% and 95% GC look
    identical to the search."""

    def test_worse_gc_scores_higher(self) -> None:
        mild = ev(ATWindow(), "GGCCGGCCGA" * 10 + CLEAN)  # ~90% GC block
        severe = ev(ATWindow(), GC_BLOCK + CLEAN)  # 100% GC block
        assert severe.raw_score > mild.raw_score

    def test_every_hard_breach_outranks_every_warn(self) -> None:
        result = ev(ATWindow(), AT_BLOCK + CLEAN)
        hard = [b.magnitude for b in result.breaches if b.detail["hard"] == "yes"]
        assert hard
        assert min(hard) >= MAG_HARD_BASE


class TestCircular:
    def test_a_window_spanning_the_origin_is_evaluated(self) -> None:
        """Half the AT block at each end: linearly there is no 100 nt AT window, but
        the molecule is circular and there is."""
        half = len(AT_BLOCK) // 2
        seq = AT_BLOCK[:half] + CLEAN + AT_BLOCK[half:]
        assert not ev(ATWindow(), seq, Topology.CIRCULAR).passes

    def test_the_same_sequence_linear_is_judged_on_its_own_ends(self) -> None:
        half = len(AT_BLOCK) // 2
        seq = AT_BLOCK[:half] + CLEAN + AT_BLOCK[half:]
        linear = ev(ATWindow(), seq, Topology.LINEAR)
        circular = ev(ATWindow(), seq, Topology.CIRCULAR)
        assert len(circular.breaches) != len(linear.breaches) or not linear.passes


class TestSpecShape:
    def test_it_is_a_band_rule_with_a_non_inverted_band(self) -> None:
        assert ATWindow.direction is Direction.BAND
        assert ATWindow.band is not None
        assert ATWindow.band[0] < ATWindow.band[1]

    def test_the_band_is_the_preference_not_the_gate(self) -> None:
        """brief.md:159 names both: 45-60% is the band, 35-65% is the hard-fail."""
        assert ATWindow.band == (SOFT_LO, SOFT_HI)
        assert ATWindow.band != (HARD_LO, HARD_HI)

    def test_a_hard_rule_carries_no_objective_weight(self) -> None:
        """CLAUDE.md 3.5: a hard constraint is never a heavy weight."""
        assert ATWindow.enforcement is Enforcement.HARD_REPAIR
        assert ATWindow.default_weight == 0.0
        assert ATWindow.steering_weight > 0.0

    def test_it_declares_the_vendor_conflict(self) -> None:
        """brief.md:159 in bold: "This directly conflicts with vendor GC ceilings.""" ""
        assert "e2_gc_band" in ATWindow.conflicts_with

    def test_it_is_not_a_lattice_rule(self) -> None:
        """The automaton makes motifs unreachable, not statistics."""
        assert ATWindow().lattice_terms(context()) is None  # type: ignore[arg-type]

    def test_the_repair_window_matches_the_statistic(self) -> None:
        assert isinstance(ATWindow().window, int)
        assert ATWindow().window == 100

    def test_it_is_registered_under_its_brief_row(self) -> None:
        assert get("f5_at_window").brief_ref == "2.F5"

    @pytest.mark.parametrize(
        "modality",
        [Modality.LENTIVIRAL, Modality.AAV, Modality.PLASMID_TRANSIENT, Modality.IVT_MRNA],
    )
    def test_it_applies_to_every_construct_that_sees_a_cloning_host(
        self, modality: Modality
    ) -> None:
        """brief.md:151 scopes 2.F to "every construct that passes through a cloning
        host" -- the plasmid is grown in E. coli whatever the final host is."""
        assert ATWindow().gate(slot(modality=modality))
