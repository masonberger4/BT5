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
    lean_targets,
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
        moves the CODON choice, which is the space G4 measures.

        Stated as "which codon each lean PICKS", not as a ranking of two chosen
        codons. The leans aim at targets inside the band rather than at 0% and
        100% GC, so the comparison that matters is between the two argmins --
        the earlier version compared two hand-picked codons and only held while
        the term was monotonic.
        """
        band = (0.30, 0.70)  # targets land at 0.44 and 0.56
        toward_at = blended_scorer({"gc_lean_at": 1.0}, usage=USAGE, gc_bounds=band)
        toward_gc = blended_scorer({"gc_lean_gc": 1.0}, usage=USAGE, gc_bounds=band)
        # Must span 1/3 AND 2/3 GC. A codon's GC is quantised to {0, 1/3, 2/3, 1},
        # and targets inside a sane band both sit between 1/3 and 2/3 -- so those
        # two are the values the leans actually choose between. A set offering
        # only {0, 2/3, 1} sends both leans to 2/3 and proves nothing, which is
        # how the first version of this test failed.
        codons = ("TTA", "GAT", "CTG", "GGC")  # 0, 1/3, 2/3, 1
        pick_at = min(codons, key=lambda c: toward_at(0, c, ""))
        pick_gc = min(codons, key=lambda c: toward_gc(0, c, ""))
        assert gc_fraction(pick_at) < gc_fraction(pick_gc), (
            f"the leans chose {pick_at} ({gc_fraction(pick_at):.0%} GC) and "
            f"{pick_gc} ({gc_fraction(pick_gc):.0%} GC) -- they must separate"
        )

    def test_a_lean_aims_inside_the_band_never_at_its_edge(self) -> None:
        """The regression that hung CI for two hours.

        The leans were a flat per-codon preference for GC or for AT, so the DP --
        which minimises a sum -- drove every codon to the extreme it could reach.
        `oracle_bounds()` passes e2's wide band (0.28-0.77), so Tier A emitted a
        ~70% GC CDS while `f5_at_window` hard-fails any 100-nt window outside
        35-65% GC. Tier B was handed a sequence out of band along its whole
        length and could not converge: an unsteered solve took 0.6 s and either
        lean ran past 70 s.

        Both targets must sit strictly inside the band AND clear of its edges,
        which is what keeps them inside the tighter windowed rules that live
        there.
        """
        for band in ((0.28, 0.77), (0.30, 0.70), (0.40, 0.60), (0.0, 1.0)):
            low, high = lean_targets(band)
            assert band[0] < low < high < band[1], f"targets escaped {band}"
            margin = 0.3 * (band[1] - band[0])
            assert low >= band[0] + margin - 1e-9
            assert high <= band[1] - margin + 1e-9

    def test_the_shipped_band_aims_inside_f5s_window(self) -> None:
        """Concrete, on the band the design actually passes.

        `f5_at_window` hard-fails a 100-nt window outside 35-65% GC and prefers
        45-60%. The band `oracle_bounds()` hands the solver is e2's 0.28-0.77,
        which is far wider -- so a lean aimed at ITS edges lands outside f5. Both
        targets must land inside f5's preferred window.
        """
        low, high = lean_targets((0.28, 0.77))
        assert 0.45 <= low <= 0.60, f"low lean target {low:.3f} is outside f5's window"
        assert 0.45 <= high <= 0.60, f"high lean target {high:.3f} is outside f5's window"

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
