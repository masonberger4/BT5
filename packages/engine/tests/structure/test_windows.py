"""Where to fold. No ViennaRNA in the room, on purpose.

Folding the right sequence at the wrong offset returns a number that looks
completely reasonable -- correct units, plausible magnitude, no exception. So the
arithmetic that decides the offset is tested on its own, where a wrong answer is
a wrong interval rather than a slightly-off energy nobody can eyeball.
"""

from __future__ import annotations

import pytest
from bt5.core.types import Interval
from bt5.structure import (
    KUDLA_DOWNSTREAM,
    KUDLA_UPSTREAM,
    five_prime_window,
    sliding_windows,
    windows_touching,
)

WIDTH = KUDLA_UPSTREAM + KUDLA_DOWNSTREAM  # 41 nt


class TestFivePrimeWindow:
    """Kudla's -4..+37 window spans the UTR/CDS junction, so it cannot be
    computed from the CDS alone. That is the whole reason it takes construct
    geometry rather than a coding sequence."""

    def test_it_straddles_the_start_codon(self) -> None:
        w = five_prime_window(Interval(500, 1400), length=3000, circular=True)
        assert w == Interval(500 - KUDLA_UPSTREAM, 500 + KUDLA_DOWNSTREAM, 1)
        assert w.start < 500 < w.end
        assert w.length == WIDTH

    def test_a_reverse_cassette_takes_its_leader_from_higher_coordinates(self) -> None:
        """Getting this backwards folds the 3' end and calls it the 5' term."""
        w = five_prime_window(Interval(500, 1400, -1), length=3000, circular=True)
        assert w is not None
        assert w.strand == -1
        assert w == Interval(1400 - KUDLA_DOWNSTREAM, 1400 + KUDLA_UPSTREAM, -1)
        assert w.start < 1400 < w.end, "the window straddles the CDS's 5' edge, which is `end`"

    def test_the_two_strands_are_mirror_images(self) -> None:
        fwd = five_prime_window(Interval(500, 1400), length=3000, circular=True)
        rev = five_prime_window(Interval(500, 1400, -1), length=3000, circular=True)
        assert fwd is not None
        assert rev is not None
        assert fwd.length == rev.length
        assert fwd.start != rev.start, "a strand-blind implementation returns the same span"

    def test_a_window_across_the_origin_wraps(self) -> None:
        w = five_prime_window(Interval(2, 900), length=3000, circular=True)
        assert w is not None
        assert w.start == 3000 - 2, "wraps rather than going negative"
        assert w.end > 3000
        assert w.length == WIDTH

    def test_a_linear_construct_with_no_room_upstream_returns_none(self) -> None:
        """None, not a truncated window. A window clamped to what fits is a
        DIFFERENT window, and comparing it to a full one compares two things."""
        assert five_prime_window(Interval(2, 900), length=3000, circular=False) is None

    def test_a_linear_construct_with_no_room_downstream_returns_none(self) -> None:
        assert five_prime_window(Interval(2980, 2999), length=3000, circular=False) is None

    def test_a_reverse_cassette_at_the_far_end_of_a_linear_construct_returns_none(self) -> None:
        """The only way to overrun the high end, and it takes a reverse cassette.

        On the forward strand a CDS at least `downstream` long guarantees the
        window fits, so the upper-bound check is unreachable there. Reversed, the
        window extends `upstream` bases PAST the CDS's end, which on a linear
        construct can run off the edge.
        """
        cds = Interval(2900, 2999, -1)
        assert cds.length >= KUDLA_DOWNSTREAM, "otherwise the length guard catches it first"
        assert five_prime_window(cds, length=3000, circular=False) is None
        wrapped = five_prime_window(cds, length=3000, circular=True)
        assert wrapped is not None, "the same layout is fine on a circular molecule"
        assert wrapped.end > 3000

    def test_a_cds_shorter_than_the_window_returns_none(self) -> None:
        assert five_prime_window(Interval(500, 520), length=3000, circular=True) is None

    def test_the_window_width_is_a_parameter_and_matters(self) -> None:
        """The Kudla window carries only 4 leader bases, so on real deposits the
        in-context and CDS-alone energies can differ by as little as 0.10
        kcal/mol at that width and by ~9 at a 53-base leader. How much leader the
        window includes decides whether the leader matters at all."""
        narrow = five_prime_window(Interval(500, 1400), length=3000, circular=True)
        wide = five_prime_window(Interval(500, 1400), length=3000, circular=True, upstream=53)
        assert narrow is not None
        assert wide is not None
        assert wide.length > narrow.length
        assert wide.start < narrow.start

    def test_a_window_that_does_not_reach_the_cds_is_refused(self) -> None:
        with pytest.raises(ValueError, match="must extend into the CDS"):
            five_prime_window(Interval(500, 1400), length=3000, circular=True, downstream=0)


class TestSlidingWindows:
    def test_a_linear_sweep_emits_no_partial_final_window(self) -> None:
        """A short last window is not comparable to the others."""
        plan = sliding_windows(1000, size=100, step=10)
        assert all(w.length == 100 for w in plan)
        assert max(w.end for w in plan) <= 1000

    def test_a_circular_sweep_continues_past_the_origin(self) -> None:
        plan = sliding_windows(1000, size=100, step=10, circular=True)
        assert any(w.end > 1000 for w in plan), "structure across position 0 must be covered"
        assert all(w.length == 100 for w in plan)

    def test_every_base_of_a_circular_construct_is_covered(self) -> None:
        plan = sliding_windows(300, size=100, step=10, circular=True)
        covered = set()
        for w in plan:
            covered.update(p % 300 for p in range(w.start, w.end))
        assert covered == set(range(300))

    def test_a_construct_shorter_than_one_window_yields_nothing_linear(self) -> None:
        assert sliding_windows(50, size=100, step=10) == ()

    @pytest.mark.parametrize(("size", "step"), [(0, 10), (100, 0), (-1, 10)])
    def test_degenerate_parameters_are_refused(self, size: int, step: int) -> None:
        with pytest.raises(ValueError, match="must both be positive"):
            sliding_windows(1000, size=size, step=step)


class TestIncrementalInvalidation:
    """A single codon change touches a handful of windows. Recomputing all of
    them is correct and about a hundred times slower, which is the difference
    between a sweep that fits in a search loop and one that does not."""

    def test_a_change_invalidates_only_nearby_windows(self) -> None:
        plan = sliding_windows(1000, size=100, step=10)
        hit = windows_touching(Interval(500, 503), plan, length=1000, circular=False)
        assert hit
        assert len(hit) < len(plan) // 4, "a 3 nt change must not invalidate the whole sweep"
        for i in hit:
            assert plan[i].start < 503
            assert plan[i].end > 500

    def test_untouched_windows_really_do_not_overlap_the_change(self) -> None:
        plan = sliding_windows(1000, size=100, step=10)
        hit = set(windows_touching(Interval(500, 503), plan, length=1000, circular=False))
        for i, w in enumerate(plan):
            if i not in hit:
                assert w.end <= 500 or w.start >= 503

    def test_a_change_at_the_origin_invalidates_wrapped_windows(self) -> None:
        plan = sliding_windows(1000, size=100, step=10, circular=True)
        hit = windows_touching(Interval(0, 3), plan, length=1000, circular=True)
        assert any(plan[i].end > 1000 for i in hit), (
            "a change at position 0 sits inside the windows that wrap the origin"
        )

    def test_nothing_is_invalidated_by_a_change_outside_every_window(self) -> None:
        plan = sliding_windows(200, size=100, step=10)
        assert windows_touching(Interval(150, 151), plan, length=1000, circular=False) != ()
        far = sliding_windows(1000, size=10, step=100)
        assert windows_touching(Interval(55, 56), far, length=1000, circular=False) == ()
