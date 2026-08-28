"""The null and the percentile -- the app's entire honesty claim.

BT5 never reports a predicted expression level; it reports where a design sits
against random synonymous variants of the same protein, in the same construct,
scored the same way. Every test here defends one of the three properties that
make that sentence true.
"""

from __future__ import annotations

import numpy as np
import pytest
from bt5.core.spec import Direction
from bt5.core.types import Construct, Interval, Segment, SegmentKind, Topology
from bt5.score import (
    NullDistribution,
    band_deviation,
    normalise,
    null_distribution,
    percentile_of,
    synonymous_variant,
)
from conftest import make_cds, translate


def construct_of(cds: str) -> Construct:
    return Construct(
        sequence=cds,
        topology=Topology.LINEAR,
        segments=(Segment(Interval(0, len(cds)), SegmentKind.DESIGNABLE_CDS),),
    )


def null_of(values: list[float]) -> NullDistribution:
    return NullDistribution(
        values=tuple(values), kind="uniform_synonymous", windowed_fold_only=True, seed=1
    )


class TestSynonymousVariant:
    def test_a_variant_encodes_the_same_protein(self, synonyms: dict) -> None:
        """If it did not, it would be a null for a different design."""
        cds = make_cds(120)
        rng = np.random.default_rng(7)
        for _ in range(20):
            assert translate(synonymous_variant(cds, synonyms, rng)) == translate(cds)

    def test_a_variant_keeps_its_length_and_frame(self, synonyms: dict) -> None:
        cds = make_cds(120)
        v = synonymous_variant(cds, synonyms, np.random.default_rng(7))
        assert len(v) == len(cds)

    def test_the_same_seed_gives_the_same_variant(self, synonyms: dict) -> None:
        """An unseeded null makes a percentile irreproducible, which makes the
        number it produced unfalsifiable."""
        cds = make_cds(60)
        a = synonymous_variant(cds, synonyms, np.random.default_rng(11))
        b = synonymous_variant(cds, synonyms, np.random.default_rng(11))
        assert a == b

    def test_different_seeds_give_different_variants(self, synonyms: dict) -> None:
        cds = make_cds(60)
        a = synonymous_variant(cds, synonyms, np.random.default_rng(11))
        b = synonymous_variant(cds, synonyms, np.random.default_rng(12))
        assert a != b

    def test_weights_actually_bias_the_sampling(self, synonyms: dict, usage: dict) -> None:
        """host_frequency must differ from uniform, or the null_kind field lies."""
        cds = make_cds(200)
        rng_u = np.random.default_rng(5)
        rng_w = np.random.default_rng(5)
        uniform = synonymous_variant(cds, synonyms, rng_u)
        weighted = synonymous_variant(cds, synonyms, rng_w, weights=usage)

        def fraction_ending_c(s: str) -> float:
            return sum(s[i + 2] == "C" for i in range(0, len(s), 3)) / (len(s) // 3)

        assert fraction_ending_c(weighted) > fraction_ending_c(uniform), (
            "a 10x weight on C-ending codons must show up in the sampling"
        )

    def test_a_codon_with_no_synonym_passes_through(self, synonyms: dict) -> None:
        """ATG and TGG have no alternatives; substituting one changes the protein."""
        cds = "ATGTGGATGTGGTAA"
        v = synonymous_variant(cds, synonyms, np.random.default_rng(3))
        assert v[:12] == cds[:12]


class TestNullDistribution:
    def scorer(self, c: Construct) -> float:
        return c.sequence.count("G") + c.sequence.count("C")

    def test_it_scores_the_assembled_construct_not_the_bare_cds(self, synonyms: dict) -> None:
        """The load-bearing property. A percentile against a distribution that
        never contained a backbone is measured against the wrong thing."""
        cds = make_cds(40)
        seen: list[str] = []

        def build(variant: str) -> Construct:
            seen.append(variant)
            return construct_of("GGGGCCCC" + variant + "AAAATTTT")

        null_distribution(
            cds,
            synonyms=synonyms,
            build=build,
            score=self.scorer,
            seed=1,
            n=5,
            kind="uniform_synonymous",
            windowed_fold_only=True,
        )
        assert len(seen) == 5
        assert all(translate(v) == translate(cds) for v in seen)

    def test_it_is_reproducible_from_its_seed(self, synonyms: dict) -> None:
        cds = make_cds(40)
        kw = {
            "synonyms": synonyms,
            "build": construct_of,
            "score": self.scorer,
            "n": 10,
            "kind": "uniform_synonymous",
            "windowed_fold_only": True,
        }
        assert (
            null_distribution(cds, seed=4, **kw).values
            == null_distribution(  # type: ignore[arg-type]
                cds,
                seed=4,
                **kw,  # type: ignore[arg-type]
            ).values
        )

    def test_a_null_too_small_to_support_a_percentile_is_refused(self, synonyms: dict) -> None:
        with pytest.raises(ValueError, match="cannot support a percentile"):
            null_distribution(
                make_cds(10),
                synonyms=synonyms,
                build=construct_of,
                score=self.scorer,
                seed=1,
                n=1,
                kind="uniform_synonymous",
                windowed_fold_only=True,
            )

    def test_host_frequency_without_weights_is_refused(self, synonyms: dict) -> None:
        """Declaring host_frequency and sampling uniformly would make null_kind
        a false label on every score derived from it."""
        with pytest.raises(ValueError, match="needs codon weights"):
            null_distribution(
                make_cds(10),
                synonyms=synonyms,
                build=construct_of,
                score=self.scorer,
                seed=1,
                n=5,
                kind="host_frequency",
                windowed_fold_only=True,
            )

    def test_mean_and_sd_are_the_sample_statistics(self) -> None:
        null = null_of([1.0, 2.0, 3.0, 4.0])
        assert null.n == 4
        assert null.mean == pytest.approx(2.5)
        assert null.sd == pytest.approx(np.std([1, 2, 3, 4], ddof=1))


class TestBandDeviation:
    def test_inside_the_band_is_zero(self) -> None:
        assert band_deviation(0.8, (0.70, 0.90)) == 0.0

    def test_below_and_above_measure_the_gap(self) -> None:
        assert band_deviation(0.60, (0.70, 0.90)) == pytest.approx(0.10)
        assert band_deviation(0.95, (0.70, 0.90)) == pytest.approx(0.05)

    def test_an_inverted_band_is_refused(self) -> None:
        with pytest.raises(ValueError, match="inverted"):
            band_deviation(0.8, (0.90, 0.70))


class TestPercentile:
    def test_higher_is_better(self) -> None:
        null = null_of([1.0, 2.0, 3.0, 4.0])
        assert percentile_of(5.0, null, Direction.HIGHER_IS_BETTER) == 1.0
        assert percentile_of(0.0, null, Direction.HIGHER_IS_BETTER) == 0.0
        assert percentile_of(2.5, null, Direction.HIGHER_IS_BETTER) == 0.5

    def test_lower_is_better_is_the_mirror(self) -> None:
        null = null_of([1.0, 2.0, 3.0, 4.0])
        assert percentile_of(0.0, null, Direction.LOWER_IS_BETTER) == 1.0
        assert percentile_of(5.0, null, Direction.LOWER_IS_BETTER) == 0.0

    def test_ties_count_half(self) -> None:
        """A variant scoring exactly what the design scores is neither beaten
        nor beating; awarding it wholly to one side turns a null containing the
        design itself into a 1.0 or a 0.0."""
        null = null_of([1.0, 2.0, 2.0, 3.0])
        assert percentile_of(2.0, null, Direction.HIGHER_IS_BETTER) == pytest.approx(0.5)

    def test_band_is_scored_on_distance_to_the_band(self) -> None:
        """The failure this exists to prevent: ranking CAI as higher-is-better
        drives it to 1.0, and one codon per amino acid is a perfect repeat."""
        null = null_of([0.55, 0.60, 0.98, 0.99])
        band = (0.70, 0.90)
        inside = percentile_of(0.80, null, Direction.BAND, band)
        maxed = percentile_of(1.00, null, Direction.BAND, band)
        assert inside == 1.0, "inside the band beats every out-of-band variant"
        assert maxed < inside, "maximising must NOT be rewarded"

    def test_a_band_objective_without_a_band_is_refused(self) -> None:
        with pytest.raises(ValueError, match="without its band"):
            percentile_of(0.8, null_of([1.0, 2.0]), Direction.BAND)


class TestNormalise:
    def test_it_carries_the_null_provenance_onto_the_score(self) -> None:
        null = null_of([1.0, 2.0, 3.0, 4.0])
        score = normalise(
            spec_id="b1_five_prime_structure",
            raw=5.0,
            unit="kcal/mol",
            direction=Direction.HIGHER_IS_BETTER,
            null=null,
        )
        assert score.percentile == 1.0
        assert score.null_n == 4
        assert score.null_mean == pytest.approx(2.5)
        assert score.null_kind == "uniform_synonymous"
        assert score.windowed_fold_only is True
        assert score.raw == 5.0
        assert score.unit == "kcal/mol"


class TestRandomnessIsSpentOnlyOnChoices:
    """An invariant codon consumes no draw, so the stream depends only on
    positions that have a choice. Without that, inserting a single Met into a
    protein reshuffles every codon after it and changes the entire null."""

    @pytest.mark.parametrize("weighted", [False, True])
    def test_an_all_invariant_cds_leaves_the_generator_untouched(
        self, synonyms: dict, usage: dict, weighted: bool
    ) -> None:
        """Both sampling paths, because only one of them can break.

        `rng.integers(0, 1)` does NOT advance the generator -- numpy
        short-circuits a single-value range -- so on the uniform path the guard
        is unobservable. `rng.choice(1, p=[1.0])` DOES advance it, so on the
        weighted path, which is what null_kind='host_frequency' uses by default,
        dropping the guard silently changes every variant in the null.
        """
        w = usage if weighted else None
        rng = np.random.default_rng(21)
        cds = "ATGTGGATGTGG"  # Met and Trp: no synonyms in any table
        assert synonymous_variant(cds, synonyms, rng, weights=w) == cds
        assert rng.bit_generator.state == np.random.default_rng(21).bit_generator.state, (
            "an invariant codon must not spend randomness"
        )

    @pytest.mark.parametrize("weighted", [False, True])
    def test_an_inserted_invariant_codon_does_not_shift_later_choices(
        self, synonyms: dict, usage: dict, weighted: bool
    ) -> None:
        w = usage if weighted else None
        without = synonymous_variant("CTTCTTCTT", synonyms, np.random.default_rng(33), weights=w)
        with_met = synonymous_variant(
            "ATGCTTCTTCTT", synonyms, np.random.default_rng(33), weights=w
        )
        assert with_met[3:] == without, "the Met must not consume a draw and shift the rest"

    def test_the_caller_decides_whether_the_stop_varies(self, synonyms: dict) -> None:
        """The map is where that decision belongs, so both behaviours are
        reachable and neither is baked into the sampler."""
        fixed = {**synonyms, "TAA": ["TAA"]}
        rng = np.random.default_rng(3)
        assert synonymous_variant("ATGTAA", fixed, rng).endswith("TAA")
