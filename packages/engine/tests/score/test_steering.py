"""The sweep's axes, and the property that makes them worth sweeping.

A gallery is only a gallery if the axes actually move the answer. These tests
pin the direction of each term and, most importantly, that the GC lean is the
axis that pulls two ways -- which is what gets the panel past gate G4's 15%
codon distance where a CAI-versus-repeats mixture alone does not.
"""

from __future__ import annotations

import pytest
from bt5.score.steering import (
    REPEAT_STEERING_PENALTY,
    SWEEP_AXES,
    blended_scorer,
    gc_fraction,
)

#: A toy relative-adaptiveness table: two leucine codons, one preferred.
USAGE = {"CTG": 1.0, "TTA": 0.1, "GGC": 1.0, "GGA": 0.2}

#: A prefix that CTG -- and only CTG -- extends into a repeated 9-mer. Appending
#: CTG makes the trailing "CCCCCCCTG" a second copy of the prefix's own first
#: nine bases; appending any other codon here does not repeat anything.
REPEATING_PREFIX = "CCCCCC" + "CTG" + "CCCCCC"

#: The same length with no repeated 9-mer in reach, for the negative.
UNIQUE_PREFIX = "ACGACGACG" + "TAGCTA"


class TestGcFraction:
    @pytest.mark.parametrize(
        ("codon", "expected"),
        [("GCC", 1.0), ("ATA", 0.0), ("GAT", 1 / 3), ("gcc", 1.0), ("", 0.0)],
    )
    def test_it_counts_g_and_c(self, codon: str, expected: float) -> None:
        assert gc_fraction(codon) == pytest.approx(expected)


class TestBlendedScorer:
    def test_an_absent_axis_costs_nothing(self) -> None:
        """A caller sweeping two axes pays for two. `simplex_weights` puts a 0.0
        on every axis it is not currently pushing, and a scorer that charged for
        those would make every lattice point score the same constant."""
        score = blended_scorer({}, usage=USAGE)
        assert score(0, "CTG", "") == 0.0

    def test_codon_adaptation_prefers_the_preferred_codon(self) -> None:
        """A COST, so preferring means scoring LOWER."""
        score = blended_scorer({"codon_adaptation": 1.0}, usage=USAGE)
        assert score(0, "CTG", "") < score(0, "TTA", "")

    def test_repeat_avoidance_only_charges_a_codon_that_extends_a_repeat(self) -> None:
        score = blended_scorer({"repeat_avoidance": 1.0}, usage=USAGE)
        assert score(0, "CTG", REPEATING_PREFIX) == pytest.approx(REPEAT_STEERING_PENALTY)
        assert score(0, "GGC", REPEATING_PREFIX) == 0.0
        assert score(0, "TTA", UNIQUE_PREFIX) == 0.0

    def test_the_two_gc_leans_pull_opposite_ways(self) -> None:
        """This is the axis that makes the panel diverse. Two sequences can differ
        at 30% of their bases as pure wobble; leaning GC one way and the other
        moves the CODON choice, which is the space G4 measures."""
        toward_at = blended_scorer({"gc_lean_at": 1.0}, usage=USAGE)
        toward_gc = blended_scorer({"gc_lean_gc": 1.0}, usage=USAGE)
        assert toward_at(0, "TTA", "") < toward_at(0, "CTG", "")
        assert toward_gc(0, "CTG", "") < toward_gc(0, "TTA", "")

    def test_the_steering_penalty_is_not_the_enforcement_penalty(self) -> None:
        """`repeat_breaking_scorer` uses 100.0 so a repeat can never be preferred.
        At 100 inside a weighted blend the repeat term swamps every other term at
        any non-zero weight, most of the simplex collapses onto one design, and
        the sweep stops sweeping. 4.0 still outranks the whole [0, 1] range of a
        codon-preference difference, and E5/F1 are HARD_REPAIR regardless -- the
        guarantee against repeats was never this scorer's job."""
        assert REPEAT_STEERING_PENALTY == 4.0
        score = blended_scorer({"repeat_avoidance": 0.5, "codon_adaptation": 0.5}, usage=USAGE)
        # CTG is the PREFERRED codon (w 1.0) and repeats here; TTA is the rare one
        # (w 0.1) and does not. The rare-but-unique codon must still win.
        repeating_preferred = score(0, "CTG", REPEATING_PREFIX)
        unique_unpreferred = score(0, "TTA", UNIQUE_PREFIX)
        assert unique_unpreferred < repeating_preferred

    def test_every_declared_axis_is_read(self) -> None:
        """A name in `SWEEP_AXES` that the scorer ignores is a wasted dimension:
        `simplex_weights` would spend lattice points on it and every one of them
        would solve to the same design."""
        for axis in SWEEP_AXES:
            scorer = blended_scorer({axis: 1.0}, usage=USAGE)
            costs = {scorer(0, codon, REPEATING_PREFIX) for codon in ("CTG", "TTA", "GGC", "GGA")}
            assert len(costs) > 1, f"axis {axis} does not distinguish any two codons"
