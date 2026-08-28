"""E4: GC VARIATION, which is a different question from E2's GC band.

If E4 were E2 with different constants the panel would report every GC problem
twice, so most of what matters here is the boundary. The load-bearing test is
`test_every_window_in_band_and_still_a_finding`: a fragment E2 passes cleanly
and E4 flags. If that ever stops holding, this rule has stopped earning its
place in the catalog.
"""

from __future__ import annotations

import numpy as np
import pytest
from bt5.core.registry import discover, get
from bt5.core.services import Services
from bt5.core.spec import Direction, Enforcement, Evidence
from bt5.core.types import Construct
from bt5.rules.catalog.e2_gc_band import GCBand
from bt5.rules.catalog.e4_gc_extent import (
    HARD_RATIO,
    WARN_RATIO,
    GCExtent,
    chance_dgc,
    dgc,
    dispersion_ratio,
    extremes,
    gc_windows,
)
from conftest import construct, context

ADAPTER_ON = "twist_gene_fragment_adapter_on"


def block(n: int, gc_fraction: float, seed: int = 3) -> str:
    """`n` bases at approximately the requested GC, shuffled deterministically."""
    rng = np.random.default_rng(seed)
    gc = rng.integers(0, 2, n)
    at = rng.integers(0, 2, n)
    return "".join(
        ("GC"[g] if rng.random() < gc_fraction else "AT"[a]) for g, a in zip(gc, at, strict=True)
    )


@pytest.fixture
def svc() -> Services:
    from bt5.vector.kmers import ConstructKmerIndex

    return Services(
        fold=None,
        kmer=ConstructKmerIndex,
        tables=None,  # type: ignore[arg-type]
        rng=np.random.default_rng(42),
    )


def run(cds: str, svc: Services, **kw: object):
    c: Construct = construct(cds, block(300, 0.5, 11))
    return GCExtent(**kw).evaluate(c, context(), svc)  # type: ignore[arg-type]


class TestTheBoundaryWithE2:
    """The reason this rule exists at all."""

    def test_every_window_in_band_and_still_a_finding(self, svc: Services) -> None:
        """A gradient that never leaves E2's band but swings across the molecule.

        E2 asks whether any window is out of band and this one is not, so E2 is
        silent. The two ends still anneal at different temperatures, which is
        what breaks assembly, and only a variation measure sees it.
        """
        cds = block(600, 0.35, 5) + block(600, 0.65, 6)
        c = construct(cds, block(300, 0.5, 11))
        ctx = context()

        e2 = GCBand(gc_min=0.10, gc_max=0.90).evaluate(c, ctx, svc)
        e4 = GCExtent().evaluate(c, ctx, svc)

        assert e2.passes, "E2 must be clean, or this proves nothing about the boundary"
        assert not e2.breaches
        assert e4.breaches, "E4 must see the excursion E2 structurally cannot"

    def test_uniform_random_sequence_is_clean(self, svc: Services) -> None:
        """The test that killed the first design.

        Random 50% GC DNA has a 50 bp extent of ~36 points at 1200 bp purely by
        chance, so the brief's "target <= 25" fired on every sequence. Against
        the binomial floor the same sequence reads ~1.0 and is silent.
        """
        cds = block(1200, 0.5, 7)
        c = construct(cds, block(300, 0.5, 11))
        ev = GCExtent().evaluate(c, context(), svc)
        assert ev.raw_score == pytest.approx(1.0, abs=0.25)
        assert not ev.breaches


class TestTheScalarIsDGC:
    def test_a_two_block_fragment_scores_above_a_uniform_one(self, svc: Services) -> None:
        uniform = run(block(1200, 0.5, 7), svc).raw_score
        split = run(block(600, 0.25, 5) + block(600, 0.75, 6), svc).raw_score
        assert uniform == pytest.approx(1.0, abs=0.25), "chance reads as chance"
        assert split > uniform * 2

    def test_the_floor_tracks_composition_not_just_window(self) -> None:
        """A 70% GC fragment has a genuinely lower binomial floor than a 50% one,
        which is why the scalar is a ratio: comparing raw SDs across synonymous
        variants would rank composition rather than dispersion."""
        assert chance_dgc("GC" * 250 + "AT" * 250) == pytest.approx(5.0, abs=0.01)
        assert chance_dgc("G" * 700 + "AT" * 150) < 4.7

    def test_a_constant_sequence_is_below_chance_not_above(self) -> None:
        """Zero dispersion is the best possible, and the ratio must say so."""
        assert dispersion_ratio("GC" * 400) == pytest.approx(0.0, abs=1e-9)

    def test_dgc_is_zero_for_a_constant_sequence(self) -> None:
        assert dgc("GC" * 400) == pytest.approx(0.0, abs=1e-9)

    def test_dgc_needs_more_than_one_window(self) -> None:
        """Under one window there is no dispersion to measure, and 0.0 is the
        honest answer rather than a number derived from a single point."""
        assert dgc("ACGT" * 10, window=100) == 0.0

    def test_the_declared_unit_is_the_scalar_actually_returned(self, svc: Services) -> None:
        spec = get("e4_gc_extent")
        assert "relative to chance" in spec.unit
        assert spec.direction is Direction.LOWER_IS_BETTER


class TestTheBreachIsExtent:
    def test_it_reports_both_ends_and_anchors_on_the_outlier(self, svc: Services) -> None:
        """Most of the molecule sits near 50%; one short stretch is GC-rich."""
        cds = block(500, 0.5, 5) + block(120, 0.95, 6) + block(500, 0.5, 8)
        breach = run(cds, svc).breaches[0]
        assert breach.detail["outlier_side"] == "high"
        assert breach.detail["max_gc_pct"] > breach.detail["median_gc_pct"]
        assert "median" in breach.message

    def test_and_anchors_low_when_the_low_end_is_the_outlier(self, svc: Services) -> None:
        cds = block(500, 0.5, 5) + block(120, 0.05, 6) + block(500, 0.5, 8)
        assert run(cds, svc).breaches[0].detail["outlier_side"] == "low"

    def test_magnitude_scales_with_dispersion_and_saturates(self, svc: Services) -> None:
        wide = run(block(600, 0.05, 5) + block(600, 0.95, 6), svc)
        assert wide.breaches[0].detail["dispersion_ratio"] >= HARD_RATIO
        assert wide.breaches[0].magnitude == 1.0

        modest = run(block(600, 0.44, 5) + block(600, 0.56, 6), svc)
        assert modest.breaches, "above the reporting threshold"
        assert 0.0 < modest.breaches[0].magnitude < 1.0

    def test_the_breach_reports_chance_alongside_the_measurement(self, svc: Services) -> None:
        """A dispersion number is meaningless without the floor it beat."""
        detail = run(block(600, 0.05, 5) + block(600, 0.95, 6), svc).breaches[0].detail
        assert detail["chance_dgc_pct"] > 0
        assert detail["dgc_pct"] > detail["chance_dgc_pct"]

    def test_findings_are_fixable_by_codon_choice(self, svc: Services) -> None:
        """Unlike e9, this one the solver CAN act on -- which is the whole reason
        it carries a steering weight."""
        for b in run(block(600, 0.05, 5) + block(600, 0.95, 6), svc).breaches:
            assert b.fixable_by_codon_choice is True


class TestWindowing:
    def test_windows_step_by_one(self) -> None:
        """A stride can step over the extreme window and under-report the range."""
        assert len(gc_windows("ACGT" * 25, 50)) == 100 - 50 + 1

    def test_a_short_sequence_yields_one_whole_sequence_value(self) -> None:
        assert gc_windows("GGCC", 50) == [100.0]
        assert gc_windows("", 50) == []

    def test_extremes_returns_both_positions(self) -> None:
        lo_i, lo, hi_i, hi = extremes([50.0, 10.0, 90.0, 50.0])
        assert (lo_i, lo, hi_i, hi) == (1, 10.0, 2, 90.0)

    def test_windows_under_20_bp_are_refused(self) -> None:
        """Below ~20 bp the number measures codon composition, not a gradient."""
        with pytest.raises(ValueError, match="composition gradient"):
            GCExtent(extent_window=12)

    def test_a_threshold_at_or_below_chance_is_refused(self) -> None:
        """1.0 IS chance, so it would report a perfectly uniform sequence."""
        with pytest.raises(ValueError, match="at or below chance"):
            GCExtent(warn_ratio=1.0)
        assert GCExtent(warn_ratio=WARN_RATIO).warn_ratio == WARN_RATIO


class TestAdapters:
    def test_an_adapter_can_never_be_the_reported_outlier(self, svc: Services) -> None:
        """Documents why the adapter-fallback branch in `_extent` is defensive.

        For an adapter window to set the extent, every insert window must sit
        closer to the median than the adapter does -- and a fragment that uniform
        does not reach the dispersion gate. A 200 bp insert at 97% GC with both
        adapters reads BELOW chance, so nothing is reported at all.
        """
        ev = run(block(200, 0.97, 5), svc, vendor=ADAPTER_ON, extent_window=20)
        assert ev.raw_score < 1.0
        assert not ev.breaches

    def test_every_reported_interval_is_a_real_construct_coordinate(self, svc: Services) -> None:
        cds = block(600, 0.05, 5) + block(600, 0.95, 6)
        for b in run(cds, svc, vendor=ADAPTER_ON).breaches:
            assert b.interval.start >= 0
            assert b.interval.end > b.interval.start

    def test_the_default_order_has_no_adapters_to_skew_gc(self, svc: Services) -> None:
        plain = run(block(900, 0.97, 5), svc).raw_score
        with_adapters = run(block(900, 0.97, 5), svc, vendor=ADAPTER_ON).raw_score
        assert with_adapters > plain, "44 bp of ~55% GC widens the spread"


class TestTheContract:
    def test_it_is_soft_because_no_vendor_publishes_this_threshold(self) -> None:
        discover()
        spec = get("e4_gc_extent")
        assert spec.enforcement is Enforcement.SOFT
        assert spec.default_weight > 0.0
        assert spec.weight_provenance.strip()

    def test_the_evidence_is_contested_and_both_signs_are_carried(self) -> None:
        """Two synthesis models rank repeats-vs-GC in opposite orders. A single
        badge or a single citation would hide exactly the disagreement that
        makes this weight a judgement call."""
        spec = get("e4_gc_extent")
        assert spec.evidence is Evidence.CONTESTED
        signs = {c.sign for c in spec.citations}
        assert "supports" in signs
        assert "refutes" in signs

    def test_it_declares_the_conflict_with_e2_and_f5(self) -> None:
        spec = get("e4_gc_extent")
        assert "e2_gc_band" in spec.conflicts_with
        assert "f5_at_window" in spec.conflicts_with
