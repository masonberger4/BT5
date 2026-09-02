"""The percentile, and the four ways it could quietly become a lie.

Every test here is about a number that would still LOOK like a percentile if the
code were wrong: a NaN scored as 0.0, a null built on bare CDSs, a uniform null
reported as a host-frequency one, a weighted total that silently drops its
highest-weight term. None of those raises; all of them put a confident number on
a report about a quantity nobody measured.
"""

from __future__ import annotations

import math

import pytest
from bt5.codon.tables import FileTableProvider
from bt5.core.spec import Breach, Direction, Evaluation
from bt5.core.types import Interval
from bt5.design.ranking import (
    NULL_N_BY_COST,
    Nulls,
    null_size,
    score_candidate,
    synonymous_map,
    unavailability,
    weighted_total,
)
from bt5.score.null import DEFAULT_NULL_N, NullDistribution

CODE = FileTableProvider().genetic_code(1)


class _Objective:
    """A minimal stand-in for a scored Spec.

    Hand-rolled rather than a real catalog rule: these tests are about what
    `score_candidate` does with an evaluation, and pinning them to c1's or e4's
    current thresholds would make them fail when the rules lane retunes a number
    that has nothing to do with this file.
    """

    def __init__(
        self,
        spec_id: str,
        *,
        unit: str = "count",
        direction: Direction = Direction.LOWER_IS_BETTER,
        band: tuple[float, float] | None = None,
        weight: float = 1.0,
        cost_class: str = "cheap",
    ) -> None:
        self.id = spec_id
        self.unit = unit
        self.direction = direction
        self.band = band
        self.default_weight = weight
        self.cost_class = cost_class


def _evaluation(spec_id: str, raw: float, *, reason: str = "") -> Evaluation:
    breaches = (
        (
            Breach(
                spec_id=spec_id,
                interval=Interval(0, 3),
                magnitude=0.0,
                message=f"unavailable: {reason}",
                fixable_by_codon_choice=False,
                detail={"unavailable_reason": reason},
            ),
        )
        if reason
        else ()
    )
    return Evaluation(spec_id=spec_id, passes=True, raw_score=raw, breaches=breaches)


def _null(values: tuple[float, ...], *, kind: str = "host_frequency") -> NullDistribution:
    return NullDistribution(values=values, kind=kind, windowed_fold_only=True, seed=0)  # type: ignore[arg-type]


class TestSynonymousMap:
    def test_every_codon_maps_to_its_own_family(self) -> None:
        mapping = synonymous_map(CODE)
        assert "ATG" in mapping
        assert set(mapping["TTA"]) == set(CODE.synonymous_codons("L"))
        assert "TTA" in mapping["CTG"]

    def test_stop_codons_are_absent_so_the_terminator_never_varies(self) -> None:
        """A null variant that changes the stop codon is a variant of a different
        construct. `code.families()` already excludes stops -- tables 27 and 28
        make TGA both Trp and a stop -- and `synonymous_variant` passes a codon
        the map does not offer through untouched."""
        mapping = synonymous_map(CODE)
        for stop in CODE.stop_codons:
            assert stop not in mapping


class TestNullSize:
    def test_cost_class_drives_it(self) -> None:
        assert null_size(_Objective("a", cost_class="cheap")) == DEFAULT_NULL_N
        assert null_size(_Objective("b", cost_class="moderate")) == NULL_N_BY_COST["moderate"]
        assert null_size(_Objective("c", cost_class="expensive")) == 0

    def test_an_unknown_cost_class_gets_no_null_rather_than_the_shipped_one(self) -> None:
        """The failure of guessing the other way is a 10 s budget blown by a rule
        nobody sized, and it would surface as a timing gate rather than as the
        missing sizing it is."""
        assert null_size(_Objective("d", cost_class="astronomical")) == 0

    def test_an_explicit_override_wins(self) -> None:
        assert null_size(_Objective("a", cost_class="cheap"), {"a": 7}) == 7


class TestUnavailability:
    def test_a_finite_score_is_a_measurement(self) -> None:
        assert unavailability(_evaluation("a", 0.42)) is None

    def test_nan_carries_the_rule_s_own_sentence(self) -> None:
        reason = unavailability(_evaluation("c1_cai", float("nan"), reason="no reference set"))
        assert reason == "no reference set"

    def test_nan_without_a_stated_reason_still_reports_one(self) -> None:
        reason = unavailability(_evaluation("a", float("nan")))
        assert reason is not None
        assert "no measurable score" in reason


class TestScoreCandidate:
    """`score_candidate` is the last place a NaN can become a number."""

    def test_a_nan_raw_score_never_becomes_a_percentile(self) -> None:
        """The whole reason `unavailability` exists.

        Every comparison against NaN is False, so `percentile_of` would count
        `better = 0`, `ties = 0` and return a confident 0.0 -- "worse than every
        one of these variants" -- about a quantity the rule explicitly declined
        to measure. Nothing raises; the report just becomes false.
        """
        spec = _Objective("c1_cai", unit="CAI")
        nulls = Nulls(by_spec={"c1_cai": _null((0.1, 0.2, 0.3))}, unavailable={}, seed=0)
        card = score_candidate(
            {"c1_cai": _evaluation("c1_cai", float("nan"), reason="no reference set for hek293")},
            [spec],  # type: ignore[list-item]
            nulls=nulls,
        )
        assert card.available == ()
        assert len(card.unavailable) == 1
        assert card.unavailable[0].unavailable_reason == "no reference set for hek293"
        assert math.isnan(card.unavailable[0].percentile)

    def test_an_objective_with_no_null_is_named_not_dropped(self) -> None:
        spec = _Objective("e8_kmer_uniqueness")
        nulls = Nulls(by_spec={}, unavailable={"e8_kmer_uniqueness": "too expensive"}, seed=0)
        card = score_candidate(
            {"e8_kmer_uniqueness": _evaluation("e8_kmer_uniqueness", 3.0)},
            [spec],  # type: ignore[list-item]
            nulls=nulls,
        )
        assert {s.spec_id for s in card.unavailable} == {"e8_kmer_uniqueness"}
        assert card.unavailable[0].unavailable_reason == "too expensive"

    def test_a_measured_objective_gets_its_percentile_and_its_null(self) -> None:
        spec = _Objective("e4_gc_extent", unit="fraction")
        nulls = Nulls(by_spec={"e4_gc_extent": _null((1.0, 2.0, 3.0, 4.0))}, unavailable={}, seed=0)
        card = score_candidate(
            {"e4_gc_extent": _evaluation("e4_gc_extent", 1.5)},
            [spec],  # type: ignore[list-item]
            nulls=nulls,
        )
        (score,) = card.available
        # LOWER_IS_BETTER: 1.5 beats 2.0, 3.0 and 4.0 -- three of four.
        assert score.percentile == pytest.approx(0.75)
        assert score.null_n == 4
        assert score.null_kind == "host_frequency"

    def test_the_null_kind_travels_onto_the_score(self) -> None:
        """A uniform-synonymous null reported as a host-frequency one is a
        different claim about what the design was compared against."""
        spec = _Objective("e4_gc_extent")
        nulls = Nulls(
            by_spec={"e4_gc_extent": _null((1.0, 2.0), kind="uniform_synonymous")},
            unavailable={},
            seed=0,
        )
        card = score_candidate(
            {"e4_gc_extent": _evaluation("e4_gc_extent", 0.5)},
            [spec],  # type: ignore[list-item]
            nulls=nulls,
        )
        assert card.available[0].null_kind == "uniform_synonymous"

    def test_a_band_objective_is_scored_on_distance_to_the_band(self) -> None:
        """Ranking CAI as higher-is-better drives it to 1.0, and one codon per
        amino acid is the documented failure this project is organised around."""
        spec = _Objective("c1_cai", direction=Direction.BAND, band=(0.7, 0.9))
        # 0.80 is inside the band; 0.95 and 0.50 are outside it in both directions.
        nulls = Nulls(by_spec={"c1_cai": _null((0.95, 0.50, 0.99))}, unavailable={}, seed=0)
        card = score_candidate(
            {"c1_cai": _evaluation("c1_cai", 0.80)},
            [spec],  # type: ignore[list-item]
            nulls=nulls,
        )
        assert card.available[0].percentile == pytest.approx(1.0)

    def test_an_objective_the_rules_lane_dropped_is_still_on_the_card(self) -> None:
        """A rule `check_engine_calibration` removed is absent from the rule set
        entirely, and an absent objective looks exactly like one that was never
        configured."""
        nulls = Nulls(by_spec={}, unavailable={}, seed=0)
        card = score_candidate(
            {},
            [],
            nulls=nulls,
            extra_unavailable={"b1_five_prime": ("kcal/mol", "no folding engine")},
        )
        assert [s.spec_id for s in card.unavailable] == ["b1_five_prime"]
        assert card.unavailable[0].unit == "kcal/mol"


class TestWeightedTotal:
    def test_it_renormalises_over_what_was_actually_evaluated(self) -> None:
        """A run missing its highest-weight objective must not score lower than
        one that has it -- the two are not comparable, and pretending the missing
        term scored zero would make the gallery's order an artefact of what
        happened to be installed."""
        specs = [_Objective("a", weight=3.0), _Objective("b", weight=1.0)]
        nulls = Nulls(
            by_spec={"b": _null((0.0, 1.0))},
            unavailable={"a": "no folding engine"},
            seed=0,
        )
        card = score_candidate(
            {
                "a": _evaluation("a", float("nan"), reason="no folding engine"),
                "b": _evaluation("b", -1.0),
            },
            specs,  # type: ignore[arg-type]
            nulls=nulls,
        )
        # Only b was measured, and it beats both null values, so the total is b's
        # own percentile -- not b's percentile diluted by a zero for a.
        assert card.total == pytest.approx(1.0)

    def test_nothing_available_is_zero_rather_than_a_division_by_zero(self) -> None:
        assert weighted_total([], [_Objective("a")]) == 0.0  # type: ignore[list-item]
