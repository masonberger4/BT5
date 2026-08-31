"""E2: per-fragment GC gates, not a window — pinned to the vendor calibration.

The load-bearing test is `test_matches_the_measured_vendor_line`: under the
gBlocks default E2 must breach the two sequences both vendors refused and pass
30–75%, and it must SPLIT on 80% — refused by IDT, made by Twist — because the
ceiling is the selected vendor's (#43 V3b). When that holds, E2 is the
manufacturability gate the 18-probe ladder measured; if it ever stops, E2 has
drifted back toward the folklore band. See docs/design/vendor-gc-calibration.md.
"""

from __future__ import annotations

import numpy as np
import pytest
from bt5.core.services import Services
from bt5.core.spec import Direction, Enforcement
from bt5.core.types import Construct
from bt5.rules.catalog.e2_gc_band import GCBand
from bt5.rules.vendors import VendorSelection
from conftest import construct, context


@pytest.fixture
def svc() -> Services:
    from bt5.vector.kmers import ConstructKmerIndex

    return Services(
        fold=None,
        kmer=ConstructKmerIndex,
        tables=None,  # type: ignore[arg-type]
        rng=np.random.default_rng(42),
    )


def block(n: int, gc: float, seed: int) -> str:
    """Exactly round(gc*n) strong bases, shuffled -- global GC hits the target to
    the base, matching the vendor panel's exact-count construction. E2 reads only
    global GC, so homopolymers from the shuffle are irrelevant here."""
    rng = np.random.default_rng(seed)
    n_gc = round(gc * n)
    bases = ["GC"[rng.integers(0, 2)] for _ in range(n_gc)]
    bases += ["AT"[rng.integers(0, 2)] for _ in range(n - n_gc)]
    rng.shuffle(bases)
    return "".join(bases)


def hits(seq: str, svc: Services, **kw: object) -> bool:
    c: Construct = construct(seq)
    return bool(GCBand(**kw).evaluate(c, context(), svc).breaches)  # type: ignore[arg-type]


class TestGlobalGate:
    def test_it_is_global_not_windowed(self, svc: Services) -> None:
        """One extreme 50 bp window on a mid-GC background must NOT breach.

        This is the case the old band got wrong: a 5% window in an otherwise 50%
        sequence was refused, though both vendors make it. Global GC ~48%.
        """
        seq = block(225, 0.5, 1) + block(50, 0.05, 2) + block(225, 0.5, 3)
        assert not hits(seq, svc)

    def test_low_global_gc_breaches(self, svc: Services) -> None:
        assert hits(block(500, 0.20, 4), svc)
        assert hits(block(500, 0.25, 5), svc)

    def test_the_floor_sits_between_25_and_30(self, svc: Services) -> None:
        """Both vendors refuse 25% and accept 30%; the 0.28 floor lands between."""
        assert hits(block(500, 0.25, 6), svc)
        assert not hits(block(500, 0.30, 7), svc)

    def test_the_ceiling_is_the_selected_vendors(self, svc: Services) -> None:
        """Twist ships 80% as Standard; IDT denies above 77%. The ceiling is the
        SELECTED vendor's, so 80% breaches under the gBlocks default and passes
        for a Twist order. 90% is over both."""
        twist = VendorSelection.of("twist_gene_fragment")
        assert not hits(block(500, 0.80, 8), svc, vendors=twist)
        assert hits(block(500, 0.80, 8), svc), "the gBlocks default caps at 77%"
        assert hits(block(500, 0.90, 9), svc)
        assert hits(block(500, 0.90, 9), svc, vendors=twist)


class TestBindingSide:
    def test_it_reports_which_side_is_binding(self, svc: Services) -> None:
        low = GCBand().evaluate(construct(block(500, 0.20, 10)), context(), svc)
        assert low.binding_side == "lower"
        high = GCBand().evaluate(construct(block(500, 0.92, 11)), context(), svc)
        assert high.binding_side == "upper"

    def test_a_breach_is_fixable_when_there_is_a_cds(self, svc: Services) -> None:
        ev = GCBand().evaluate(construct(block(500, 0.20, 12)), context(), svc)
        assert ev.breaches[0].fixable_by_codon_choice is True


class TestWindowsAreDisplayOnly:
    def test_windows_are_populated_but_never_a_breach(self, svc: Services) -> None:
        """The GC landscape is still reported; a window just cannot fail the rule."""
        seq = block(225, 0.5, 1) + block(50, 0.02, 2) + block(225, 0.5, 3)
        ev = GCBand().evaluate(construct(seq), context(), svc)
        assert ev.windows, "the landscape is still there for the report"
        assert not ev.breaches, "but the 2% window is not a breach"


class TestFragmentScope:
    """The scope fix (#43 V3b): E2 measures the ordered fragment, not the whole
    construct. Before, a near-neutral backbone diluted a GC-extreme insert to
    inside the band, so the one HARD_REPAIR GC rule never fired on a vector design
    though both vendors deny the tube. No older test built an insert-in-backbone
    construct and asked E2 about it, which is why the inertness went unseen.
    """

    #: Near-50% GC, so it drags any insert's whole-construct GC toward the middle.
    @staticmethod
    def _backbone() -> str:
        return block(2000, 0.5075, 99)

    def test_a_gc_rich_insert_breaches_though_the_whole_construct_is_in_band(
        self, svc: Services
    ) -> None:
        insert = block(900, 0.90, 1)  # IDT denies a 90% fragment
        c = construct(insert, self._backbone())
        from bt5.verify import gc_fraction

        # The number the old construct-scope rule read -- comfortably in band.
        assert 0.60 < gc_fraction(c.sequence) < 0.66
        ev = GCBand().evaluate(c, context(), svc)
        assert ev.breaches, "the fragment is 90% GC; the gBlocks default caps at 77%"
        assert ev.binding_side == "upper"
        b = ev.breaches[0]
        assert b.detail["gc"] == pytest.approx(0.90, abs=0.005)
        assert b.detail["fragment_bp"] == 900.0
        assert b.detail["vendor"] == "idt_gblocks"

    def test_a_gc_poor_insert_breaches_on_the_floor(self, svc: Services) -> None:
        insert = block(900, 0.20, 3)  # too LOW at fragment scope
        c = construct(insert, self._backbone())
        ev = GCBand().evaluate(c, context(), svc)
        assert ev.binding_side == "lower"
        b = ev.breaches[0]
        # A real interval WINDOW_MINUS_1 can localize -- not Interval(0, n) over
        # the whole construct, which made codon_span return every codon.
        assert (b.interval.start, b.interval.end) == (0, 900)
        assert b.fixable_by_codon_choice is True

    def test_a_backbone_gc_excursion_is_never_e2s_to_fix(self, svc: Services) -> None:
        """An immutable GC-rich backbone with an in-band insert must pass. The
        vendor never synthesizes the backbone, and the old construct-scope rule
        set fixable True whenever a CDS existed, sending repair after a target the
        codons provably could not reach."""
        c = construct(block(600, 0.50, 8), block(2000, 0.85, 7))
        ev = GCBand().evaluate(c, context(), svc)
        assert ev.passes, "the insert is in band; the backbone is not ours to fix"


class TestTheBandIsPerSelection:
    """The gate is the intersection of the selected vendors' bands, not the
    ClassVar. Override and `none` are the two ways it is not vendor-derived."""

    def test_a_multi_vendor_selection_intersects_to_the_strictest(self, svc: Services) -> None:
        # 79% GC: inside Twist's 0.80, over IDT's 0.77. A design for BOTH must
        # satisfy both, so the intersected ceiling 0.77 refuses it, and the finding
        # names IDT as the vendor that binds.
        both = VendorSelection.of("twist_gene_fragment", "idt_gblocks")
        c = construct(block(500, 0.79, 21))
        ev = GCBand(vendors=both).evaluate(c, context(), svc)
        assert ev.breaches
        assert ev.breaches[0].detail["vendor"] == "idt_gblocks"
        assert ev.breaches[0].detail["band_hi"] == 0.77

    def test_none_falls_back_to_the_loosest_envelope(self, svc: Services) -> None:
        """`none` has no vendor band; E2 refuses only what no configuration makes,
        the (0.28, 0.80) ClassVar. 80% passes there; 90% does not."""
        none = VendorSelection.of("none")
        lo = GCBand(vendors=none).evaluate(construct(block(500, 0.80, 22)), context(), svc)
        assert not lo.breaches, "80% is inside the (0.28, 0.80) envelope"
        hi = GCBand(vendors=none).evaluate(construct(block(500, 0.90, 23)), context(), svc)
        assert hi.breaches
        assert hi.breaches[0].detail["limit_source"] == "envelope"
        assert hi.breaches[0].detail["vendor"] == ""

    def test_an_explicit_bound_overrides_and_is_attributed_to_nobody(self, svc: Services) -> None:
        c = construct(block(500, 0.75, 24))
        ev = GCBand(gc_min=0.10, gc_max=0.60).evaluate(c, context(), svc)
        assert ev.breaches
        b = ev.breaches[0]
        assert b.detail["limit_source"] == "override"
        assert b.detail["vendor"] == ""
        assert "the band you set" in b.message


class TestContract:
    def test_still_hard_repair_and_band(self) -> None:
        from bt5.core.registry import discover, get

        discover()
        spec = get("e2_gc_band")
        assert spec.enforcement is Enforcement.HARD_REPAIR
        assert spec.direction is Direction.BAND
        assert spec.band == (0.28, 0.80)
        assert spec.default_weight == 0.0

    def test_matches_the_measured_vendor_line(self, svc: Services) -> None:
        """The pin. Uniform GC from 20% to 75%: under the gBlocks default E2
        breaches 20 and 25 (both vendors refused) and passes 30-75. 80% is the
        vendor split the per-vendor band exists to capture -- above IDT's 77%
        ceiling, below Twist's 80% -- so it breaches by default and passes Twist.
        """
        for pct in (20, 25):
            assert hits(block(500, pct / 100, 100 + pct), svc), f"{pct}% must breach"
        for pct in (30, 40, 50, 60, 65, 70, 75):
            assert not hits(block(500, pct / 100, 100 + pct), svc), f"{pct}% must pass"
        twist = VendorSelection.of("twist_gene_fragment")
        assert hits(block(500, 0.80, 180), svc), "80% is over IDT's 77% ceiling"
        assert not hits(block(500, 0.80, 180), svc, vendors=twist), "80% is Standard at Twist"
