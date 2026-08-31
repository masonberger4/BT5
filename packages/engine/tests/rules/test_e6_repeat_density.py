"""E6: repetitiveness as a density, and why that is a different rule.

The pair rules (E5, F1, F2) answer "are these two spans the same?". This one
answers "how much of this fragment is not unique?", and the gap between those
questions is the case the only vendor-outcome-trained synthesis model was built
to catch: many short repeats, no long one. Several tests below exist to hold
that boundary, and several more to hold the claims the docstring makes about
length dependence -- because that is the claim justifying the absence of a hard
threshold, and an unchecked justification rots.

`conftest` is imported at module level; see the note in test_f1_direct_repeats.
"""

from __future__ import annotations

import numpy as np
import pytest
from bt5.core.context import Modality
from bt5.core.registry import discover, get
from bt5.core.services import Services
from bt5.core.spec import Direction, Enforcement
from bt5.rules.catalog.d4_internal_polya import InternalPolyA
from bt5.rules.catalog.e5_synthesis_repeats import SynthesisRepeats
from bt5.rules.catalog.e6_repeat_density import (
    KMER_BP,
    WINDOW_FLAG,
    RepeatDensity,
    repetitive_kmers_per_100bp,
)
from bt5.rules.catalog.f2_near_perfect_repeats import NearPerfectRepeats
from bt5.rules.vendors import VendorSelection
from bt5.vector.kmers import ConstructKmerIndex
from conftest import construct, context, slot

discover()

#: One codon per residue, the collapse max-CAI produces on a repetitive protein.
GGGGS_LINKER = "GGTGGTGGTGGTAGC" * 20


def dna(n: int, seed: int = 3) -> str:
    rng = np.random.default_rng(seed)
    return "".join("ACGT"[i] for i in rng.integers(0, 4, n))


@pytest.fixture
def svc() -> Services:
    return Services(
        fold=None,
        kmer=ConstructKmerIndex,
        tables=None,  # type: ignore[arg-type]
        rng=np.random.default_rng(1),
    )


def evaluate(cds: str, svc: Services, **kw: object):
    return RepeatDensity(**kw).evaluate(  # type: ignore[arg-type]
        construct(cds, dna(300, 11)), context(), svc
    )


def many_short_repeats(seed: int = 61) -> str:
    """Every 9-mer recurs; no 12-mer does.

    Ten distinct 9-mers laid down twice, with DIFFERENT 3 bp separators on the
    two passes. Any 12-mer necessarily crosses a separator, so nothing at 12 bp
    or longer is repeated -- which is precisely the fragment the pair rules call
    clean and the vendor calls repetitive.
    """
    rng = np.random.default_rng(seed)

    def draw(n: int) -> str:
        return "".join("ACGT"[i] for i in rng.integers(0, 4, n))

    units = [draw(9) for _ in range(10)]
    first = "".join(u + draw(3) for u in units)
    second = "".join(u + draw(3) for u in units)
    return dna(150, 63) + first + dna(60, 65) + second + dna(150, 67)


class TestTheDensityBoundary:
    """Many short repeats, no long one: the gap the pair rules leave."""

    def test_the_pair_rules_call_this_fragment_clean(self, svc: Services) -> None:
        c = construct(many_short_repeats(), dna(300, 11))
        assert not SynthesisRepeats().evaluate(c, context(), svc).breaches
        assert not NearPerfectRepeats().evaluate(c, context(), svc).breaches

    def test_and_this_rule_does_not(self, svc: Services) -> None:
        """If this ever passes, E6 has stopped earning its place in the catalog."""
        clean = repetitive_kmers_per_100bp(dna(len(many_short_repeats()), 71))
        assert repetitive_kmers_per_100bp(many_short_repeats()) > 5 * max(clean, 0.05)


class TestDetection:
    def test_a_one_codon_per_residue_linker_is_flagged(self, svc: Services) -> None:
        """`(GGGGS)n` back-translated with one codon per residue is a perfect
        nucleotide repeat -- the plan's named failure mode for repetitive
        proteins, and the reason max-CAI is not the objective."""
        result = evaluate(dna(300, 41) + GGGGS_LINKER + dna(300, 43), svc)
        assert result.breaches
        assert not result.passes
        assert result.breaches[0].detail["peak_repeated_fraction"] == 1.0

    def test_ordinary_coding_sequence_is_clean(self, svc: Services) -> None:
        result = evaluate(dna(900, 41), svc)
        assert result.passes
        assert not result.breaches
        assert result.windows, "a clean fragment still reports its windows"

    def test_a_duplicated_cassette_element_is_flagged(self, svc: Services) -> None:
        """Two copies of the same 2A peptide in one ORF is a perfect direct
        repeat; the literature requires <=85% nucleotide identity between them."""
        element = dna(66, 45)
        cds = dna(300, 41) + element + dna(200, 47) + element + dna(300, 43)
        assert evaluate(cds, svc).breaches

    def test_one_repetitive_region_is_one_finding(self, svc: Services) -> None:
        """Windows overlap by design, so a 300 bp repetitive stretch trips
        fifteen of them. Reporting each is the failure E1 had with homopolymers
        and F2 had with adjacent seeds."""
        assert len(evaluate(dna(300, 41) + GGGGS_LINKER + dna(300, 43), svc).breaches) == 1

    def test_the_reported_fraction_is_named_as_a_peak(self, svc: Services) -> None:
        """A merged region is not uniformly as bad as its worst window, and a
        message implying otherwise overstates every finding that merged."""
        breach = evaluate(dna(300, 41) + GGGGS_LINKER + dna(300, 43), svc).breaches[0]
        assert "at its worst" in breach.message

    def test_findings_land_on_the_repetitive_region(self, svc: Services) -> None:
        cds = dna(300, 41) + GGGGS_LINKER + dna(300, 43)
        breach = evaluate(cds, svc).breaches[0]
        assert breach.interval.start < 300 + len(GGGGS_LINKER)
        assert breach.interval.end > 300
        assert breach.fixable_by_codon_choice


class TestLengthDependence:
    """The claim that justifies having no hard threshold, held to account."""

    def test_identical_composition_scores_higher_when_longer(self) -> None:
        short = np.mean([repetitive_kmers_per_100bp(dna(900, 5000 + s)) for s in range(15)])
        long = np.mean([repetitive_kmers_per_100bp(dna(3000, 6000 + s)) for s in range(15)])
        assert long > 2 * short, (
            "the metric is expected to grow linearly in length on identical "
            "composition; if it does not, the docstring's argument against a "
            "fixed threshold no longer holds"
        )

    def test_it_tracks_the_closed_form(self) -> None:
        """Expected value is 100n / (2 * 4^k). Measuring it rather than
        asserting it is what makes the length argument a fact."""
        n = 3000
        predicted = 100 * n / (2 * 4**KMER_BP)
        measured = np.mean([repetitive_kmers_per_100bp(dna(n, 7000 + s)) for s in range(20)])
        assert measured == pytest.approx(predicted, rel=0.25)

    def test_the_rule_declares_no_band(self) -> None:
        assert RepeatDensity.direction is Direction.LOWER_IS_BETTER
        assert RepeatDensity.band is None


class TestScope:
    def test_a_repetitive_backbone_is_not_the_vendors_problem(self, svc: Services) -> None:
        """Nobody is synthesising the user's existing vector."""
        c = construct(dna(400, 41), GGGGS_LINKER + dna(200, 43))
        assert not RepeatDensity().evaluate(c, context(), svc).breaches

    def test_the_worst_fragment_sets_the_score(self, svc: Services) -> None:
        """An order succeeds only if every tube does, so the summary is a max
        across fragments and not a mean that a clean tube could dilute."""
        clean = evaluate(dna(900, 41), svc).raw_score
        mixed = evaluate(dna(300, 41) + GGGGS_LINKER + dna(300, 43), svc).raw_score
        assert mixed > clean


class TestContract:
    def test_it_is_soft_so_it_can_be_weighted(self) -> None:
        assert RepeatDensity.enforcement is Enforcement.SOFT
        assert RepeatDensity.default_weight > 0.0
        assert RepeatDensity.weight_provenance.strip()

    def test_it_outweighs_the_conventional_repeat_threshold(self) -> None:
        """f2's 90%-over-40bp is a convention; this is the highest-Gini feature
        of the only model trained on real vendor outcomes."""
        assert RepeatDensity.default_weight > NearPerfectRepeats.default_weight

    def test_but_not_the_rule_with_a_measured_effect_size(self) -> None:
        """Feature importance ranks predictors; it does not measure effect size,
        and it should not outrank a measured 8-9x functional titer loss."""
        assert RepeatDensity.default_weight < InternalPolyA.default_weight

    def test_it_does_not_stack_a_third_full_strength_repeat_steering_term(self) -> None:
        assert RepeatDensity.steering_weight < SynthesisRepeats.steering_weight

    def test_it_is_not_a_lattice_rule(self) -> None:
        assert RepeatDensity().lattice_terms(None) is None

    def test_it_applies_in_every_context(self) -> None:
        for modality in Modality:
            assert RepeatDensity().gate(slot(modality=modality))

    def test_absurd_parameters_are_refused(self) -> None:
        with pytest.raises(ValueError, match="recur by chance everywhere"):
            RepeatDensity(k=4)
        with pytest.raises(ValueError, match="at least two"):
            RepeatDensity(k=9, window=12)
        with pytest.raises(ValueError, match="fraction in"):
            RepeatDensity(window_flag=0.0)
        with pytest.raises(ValueError, match="unknown vendor"):
            RepeatDensity(vendors=VendorSelection.of("acme"))

    def test_the_flag_clears_ordinary_sequence_by_a_wide_margin(self) -> None:
        assert WINDOW_FLAG >= 0.2

    def test_it_is_registered_under_its_brief_row(self) -> None:
        assert get("e6_repeat_density") is RepeatDensity
        assert RepeatDensity.brief_ref == "2.E6"
