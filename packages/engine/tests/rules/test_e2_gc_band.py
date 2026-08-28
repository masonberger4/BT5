"""E2: GLOBAL GC gates, not a window — pinned to the vendor calibration.

The load-bearing test is `test_matches_the_measured_vendor_line`: E2 must breach
exactly the two sequences both vendors refused and pass the sixteen they made,
including eight carrying a single 4–96% GC window. When that holds, E2 is the
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

    def test_the_ceiling_admits_twists_80_percent(self, svc: Services) -> None:
        """Twist ships 80% as Standard, so the universal hard band must not refuse
        it; IDT's stricter 77% ceiling is a per-vendor concern, not this default."""
        assert not hits(block(500, 0.80, 8), svc)
        assert hits(block(500, 0.90, 9), svc)


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
        """The pin. Uniform GC from 20% to 80%: E2 breaches exactly 20 and 25,
        the two both vendors refused, and passes 30-80, which both made.
        """
        for pct in (20, 25):
            assert hits(block(500, pct / 100, 100 + pct), svc), f"{pct}% must breach"
        for pct in (30, 40, 50, 60, 65, 70, 75, 80):
            assert not hits(block(500, pct / 100, 100 + pct), svc), f"{pct}% must pass"
