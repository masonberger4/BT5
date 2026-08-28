"""The honest-degradation contract on a scorecard.

A scorecard missing its highest-weight objective looks exactly like a scorecard
where that objective was never configured, and the difference is whether the
ranking means anything. `unavailable_reason` is what makes those two states
distinguishable, so these tests are about the distinction rather than the field.
"""

from __future__ import annotations

import math

import pytest
from bt5.core.result import ObjectiveScore, ScoreCard


def scored(spec_id: str, percentile: float = 0.5) -> ObjectiveScore:
    return ObjectiveScore(
        spec_id=spec_id,
        raw=1.0,
        unit="kcal/mol",
        percentile=percentile,
        null_n=200,
        null_mean=0.0,
        null_sd=1.0,
    )


class TestUnavailableObjective:
    def test_an_evaluated_score_is_available(self) -> None:
        assert scored("b1").available
        assert scored("b1").unavailable_reason == ""

    def test_an_unavailable_score_states_why(self) -> None:
        s = ObjectiveScore.unavailable("b1", "kcal/mol", "no folding engine installed")
        assert not s.available
        assert "folding engine" in s.unavailable_reason

    def test_its_numbers_are_not_mistakable_for_measurements(self) -> None:
        """NaN rather than 0.0. A zero percentile is a real, terrible ranking;
        a zero that means 'not measured' is the exact confusion this prevents."""
        s = ObjectiveScore.unavailable("b1", "kcal/mol", "no folding engine")
        assert math.isnan(s.percentile)
        assert math.isnan(s.raw)
        assert s.null_n == 0

    def test_an_unavailable_objective_must_give_a_reason(self) -> None:
        with pytest.raises(ValueError, match="must say why"):
            ObjectiveScore.unavailable("b1", "kcal/mol", "")


class TestScoreCardPartitions:
    def test_it_separates_measured_from_unmeasured(self) -> None:
        card = ScoreCard(
            scores=(
                scored("d1"),
                ObjectiveScore.unavailable("b1", "kcal/mol", "no folding engine"),
                scored("e2"),
            )
        )
        assert [s.spec_id for s in card.available] == ["d1", "e2"]
        assert [s.spec_id for s in card.unavailable] == ["b1"]

    def test_a_fully_measured_card_reports_nothing_unavailable(self) -> None:
        card = ScoreCard(scores=(scored("d1"), scored("e2")))
        assert card.unavailable == ()
        assert len(card.available) == 2

    def test_the_weighted_sum_has_a_defined_input_set(self) -> None:
        """The point of the partition: `available` is what a sum may touch, and
        summing a NaN percentile would poison the total silently."""
        card = ScoreCard(
            scores=(scored("d1", 0.8), ObjectiveScore.unavailable("b1", "kcal/mol", "absent"))
        )
        total = sum(s.percentile for s in card.available)
        assert total == pytest.approx(0.8)
        assert not math.isnan(total)

    def test_dropping_an_objective_is_distinguishable_from_never_having_it(self) -> None:
        """The two cards below differ ONLY in whether the absence is stated.
        Before `unavailable_reason` they were the same object."""
        silent = ScoreCard(scores=(scored("d1"),))
        honest = ScoreCard(
            scores=(scored("d1"), ObjectiveScore.unavailable("b1", "kcal/mol", "no engine"))
        )
        assert len(silent.available) == len(honest.available)
        assert silent.unavailable == ()
        assert honest.unavailable != ()
