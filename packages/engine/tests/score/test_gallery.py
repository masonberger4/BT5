"""Gate G4: does the gallery actually contain different designs?

Its failure invalidates a PRODUCT decision rather than a technical one -- if the
sweep cannot produce distinct candidates, the gallery is not a feature and the UI
built on it is a lie. The plan names it as one of the two gates that should
scare you most, and unlike G3 it needs no folding engine to answer.
"""

from __future__ import annotations

from collections.abc import Mapping

import pytest
from bt5.score import (
    G4_MIN_PAIRWISE_DISTANCE,
    Gallery,
    SweepPoint,
    build_gallery,
    codon_distance,
    greedy_max_min,
    pairwise_minimum,
    simplex_weights,
    sweep,
)

#: Three designs of the same 6-codon protein, maximally spread, plus a
#: near-duplicate of the first. The near-duplicate is the trap: a selector that
#: does not look at sequence space will happily include it.
#: Twelve codons, not six: at six, a single differing codon is 1/6 = 0.167,
#: which is ABOVE G4's 0.15 floor -- so a "near-duplicate" built that way is not
#: actually near enough to fail the gate, and the failing test would not fail.
FAR_A = "CTT" * 12
FAR_B = "CTG" * 12
FAR_C = "CTC" * 12
NEAR_A = "CTT" * 11 + "CTG"  # one codon from FAR_A: 1/12 = 0.083


class TestSimplexWeights:
    def test_every_weight_vector_sums_to_one(self) -> None:
        for w in simplex_weights(["a", "b", "c"], 4):
            assert sum(w.values()) == pytest.approx(1.0)

    def test_the_corners_are_included(self) -> None:
        """The single-objective extremes are the designs a user most needs to
        see, so a sweep that misses them misses the ends of the trade-off."""
        got = simplex_weights(["a", "b"], 4)
        assert {"a": 1.0, "b": 0.0} in got
        assert {"a": 0.0, "b": 1.0} in got

    def test_density_grows_with_steps(self) -> None:
        assert len(simplex_weights(["a", "b", "c"], 8)) > len(simplex_weights(["a", "b", "c"], 4))

    def test_degenerate_input_is_refused(self) -> None:
        with pytest.raises(ValueError, match="no objectives"):
            simplex_weights([], 4)
        with pytest.raises(ValueError, match="steps must be positive"):
            simplex_weights(["a"], 0)


class TestSweep:
    def test_an_unsolvable_weight_vector_is_dropped_not_substituted(self) -> None:
        """An infeasible corner of the simplex is a real answer about the
        constraint set, not a gap to paper over."""

        def solve(w: Mapping[str, float]) -> str | None:
            return None if w["a"] == 1.0 else FAR_A

        got = sweep(["a", "b"], solve, steps=4)
        assert len(got) == 4, "one of the five lattice points was infeasible"
        assert all(p.cds == FAR_A for p in got)


class TestGreedyMaxMin:
    def points(self, *seqs: str) -> list[SweepPoint]:
        return [SweepPoint(weights={"a": float(i)}, cds=s) for i, s in enumerate(seqs)]

    def test_it_seeds_with_the_farthest_pair(self) -> None:
        """Greedy max-min is only as good as where it starts; seeding from the
        extremes is what keeps both ends of the trade-off in a small gallery."""
        picked = greedy_max_min(self.points(FAR_A, NEAR_A, FAR_B), 2)
        chosen = {p.cds for p in picked}
        assert chosen == {FAR_A, FAR_B}, "it must not pick the near-duplicate first"

    def test_it_prefers_a_distant_design_over_a_near_duplicate(self) -> None:
        picked = greedy_max_min(self.points(FAR_A, NEAR_A, FAR_B, FAR_C), 3)
        assert NEAR_A not in {p.cds for p in picked}

    def test_asking_for_more_than_exist_returns_everything(self) -> None:
        pts = self.points(FAR_A, FAR_B)
        assert len(greedy_max_min(pts, 5)) == 2

    def test_a_gallery_of_zero_is_refused(self) -> None:
        with pytest.raises(ValueError, match="not a gallery"):
            greedy_max_min(self.points(FAR_A), 0)


class TestBuildGallery:
    def solver(self, designs: list[str]):
        """Return a different design per weight vector, cycling through `designs`."""
        seen: list[int] = []

        def solve(w: Mapping[str, float]) -> str | None:
            seen.append(1)
            return designs[len(seen) % len(designs)]

        return solve

    def test_duplicates_are_collapsed_before_selection(self) -> None:
        """Many weight vectors land on the same design; a greedy pick over
        duplicates spends slots proving they are identical."""
        g = build_gallery(["a", "b"], lambda w: FAR_A, steps=8, k=3)
        assert g.swept == 9
        assert g.distinct == 1
        assert len(g.picks) == 1

    def test_a_diverse_sweep_meets_g4(self) -> None:
        g = build_gallery(["a", "b"], self.solver([FAR_A, FAR_B, FAR_C]), steps=8, k=3)
        assert g.distinct == 3
        assert g.min_pairwise_distance == 1.0
        assert g.meets_g4

    def test_a_clustered_sweep_fails_g4_loudly(self) -> None:
        """The gate must be able to FAIL. A gallery of near-duplicates has to
        report as one, or G4 cannot invalidate the product decision it exists
        to test."""
        near_b = "CTT" * 11 + "CTC"
        g = build_gallery(["a", "b"], self.solver([FAR_A, NEAR_A, near_b]), steps=8, k=3)
        assert g.distinct == 3
        assert g.min_pairwise_distance < G4_MIN_PAIRWISE_DISTANCE
        assert not g.meets_g4

    def test_gallery_size_is_bounded(self) -> None:
        with pytest.raises(ValueError, match="gallery size must be"):
            build_gallery(["a"], lambda w: FAR_A, k=2)
        with pytest.raises(ValueError, match="gallery size must be"):
            build_gallery(["a"], lambda w: FAR_A, k=9)

    def test_the_reported_minimum_covers_the_picks_not_the_whole_sweep(self) -> None:
        """It is the selected panel a user sees, so it is the selected panel
        that has to be diverse."""
        g = build_gallery(["a", "b"], self.solver([FAR_A, NEAR_A, FAR_B, FAR_C]), steps=8, k=3)
        assert g.min_pairwise_distance == pytest.approx(pairwise_minimum([p.cds for p in g.picks]))


class TestDasAndDennis:
    """Why the gallery is not selected by evenly spaced weights.

    Das & Dennis (1997) showed an even spread of weights gives an UNEVEN spread
    of solutions, because the weight tracks the local slope of the Pareto front
    rather than position along it. Here is that failure made concrete: a solver
    whose output responds non-linearly to the weight -- which is the normal case,
    not a contrived one -- clusters under even weighting and spreads under
    max-min selection over the same sweep.
    """

    def steep_solver(self, w: Mapping[str, float]) -> str:
        # Almost every weight vector lands on FAR_A; only the extreme corner
        # moves. This is a front whose slope is highly uneven.
        return FAR_B if w["a"] > 0.875 else (NEAR_A if w["a"] > 0.75 else FAR_A)

    def test_even_weights_give_a_clustered_gallery(self) -> None:
        points = sweep(["a", "b"], self.steep_solver, steps=8)
        every_other = [p.cds for p in points[::3]][:3]
        assert pairwise_minimum(every_other) < G4_MIN_PAIRWISE_DISTANCE, (
            "evenly spaced weights land on the same design over and over"
        )

    def test_max_min_selection_over_the_same_sweep_spreads(self) -> None:
        g = build_gallery(["a", "b"], self.steep_solver, steps=8, k=3)
        assert g.distinct == 3
        assert g.min_pairwise_distance > 0
        far = greedy_max_min([SweepPoint({"a": 0.0}, s) for s in (FAR_A, NEAR_A, FAR_B)], 2)
        assert {p.cds for p in far} == {FAR_A, FAR_B}
        assert codon_distance(FAR_A, FAR_B) == 1.0


def test_a_gallery_dataclass_reports_g4_honestly() -> None:
    assert Gallery((), 0, 0, 1.0).meets_g4 is False, "an empty gallery is not a passing gallery"
