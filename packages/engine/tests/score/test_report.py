"""The QC report, tested mostly for what it refuses to say.

BT5's differentiator is trustworthiness, and the report is where that is either
kept or lost: it is the last artefact a human reads before spending money on
DNA. So the assertions here are largely negative -- no predicted number, no
silently dropped objective, no defaulted genetic code.
"""

from __future__ import annotations

import re

import pytest
from bt5.core.result import Candidate, Conflict, DesignResult, ObjectiveScore, ScoreCard
from bt5.core.types import Construct, Interval, Provenance, Segment, SegmentKind, Topology
from bt5.score.report import (
    QcReport,
    build_report,
    render,
)

CDS = "ATGAAACCCGGGTTTTAA"


def construct(seq: str = CDS + "GGGCCCAAATTT") -> Construct:
    return Construct(
        sequence=seq,
        topology=Topology.CIRCULAR,
        segments=(
            Segment(Interval(0, len(CDS)), SegmentKind.DESIGNABLE_CDS, "cds"),
            Segment(Interval(len(CDS), len(seq)), SegmentKind.BACKBONE, "vector"),
        ),
    )


def score(spec_id: str, percentile: float = 0.5) -> ObjectiveScore:
    return ObjectiveScore(
        spec_id=spec_id,
        raw=1.0,
        unit="au",
        percentile=percentile,
        null_n=200,
        null_mean=0.5,
        null_sd=0.1,
    )


def candidate(*scores: ObjectiveScore, label: str = "a") -> Candidate:
    return Candidate(
        label=label,
        construct=construct(),
        cds=CDS,
        scorecard=ScoreCard(scores=scores),
        design_hash="abc123def456",
    )


class TestBuildReport:
    def result(self, *, native: bool = False, degradations: tuple[str, ...] = ()) -> DesignResult:
        cand = candidate(score("c1_cai", 0.8))
        return DesignResult(
            candidates=(cand,),
            native_baseline=cand if native else None,
            provenance=Provenance(
                app_version="0.1.0",
                seed=7,
                engine_versions={"viennarna": "2.7.2"},
                degradations=degradations,
            ),
        )

    def test_carries_the_translation_table_through(self) -> None:
        result = self.result()
        report = build_report(result, result.candidates[0], translation_table_id=11)
        assert report.translation_table_id == 11

    def test_the_table_has_no_default(self) -> None:
        """Defaulting it to 1 would be the silent-wrong-protein bug the contract
        refuses everywhere else, one layer further out."""
        result = self.result()
        with pytest.raises(TypeError):
            build_report(result, result.candidates[0])  # type: ignore[call-arg]

    def test_degradations_travel_from_provenance(self) -> None:
        result = self.result(degradations=("no folding engine installed",))
        report = build_report(result, result.candidates[0], translation_table_id=1)
        assert "no folding engine installed" in report.degradations
        assert not report.is_complete

    def test_the_native_baseline_is_reported_as_present(self) -> None:
        result = self.result(native=True)
        report = build_report(result, result.candidates[0], translation_table_id=1)
        assert report.native_baseline_available

    def test_unavailable_objectives_are_separated_not_dropped(self) -> None:
        """A scorecard missing its highest-weight term looks exactly like one
        where that term was never configured."""
        cand = candidate(
            score("c1_cai", 0.8),
            ObjectiveScore.unavailable("b1_five_prime", "kcal/mol", "no folding engine"),
        )
        result = DesignResult(candidates=(cand,))
        report = build_report(result, cand, translation_table_id=1)
        assert [s.spec_id for s in report.scored] == ["c1_cai"]
        assert [s.spec_id for s in report.unavailable] == ["b1_five_prime"]
        assert not report.is_complete

    def test_a_run_with_everything_evaluated_is_complete(self) -> None:
        """Nothing unavailable and nothing degraded means complete.

        This needed a vendor named explicitly while the report computed a
        screening burden: BT5's own default had no error-free length on file, so
        every default run carried a degradation and no run was ever complete
        without overriding the vendor. The report no longer sizes colony
        picking, so completeness is once again about the objectives alone.
        """
        result = self.result()
        report = build_report(result, result.candidates[0], translation_table_id=1)
        assert report.is_complete


class TestRender:
    def report(self, **kw: object) -> QcReport:
        cand = candidate(
            score("c1_cai", 0.82),
            ObjectiveScore.unavailable("b1_five_prime", "kcal/mol", "no folding engine installed"),
        )
        result = DesignResult(
            candidates=(cand,),
            native_baseline=cand,
            conflicts=(
                Conflict(
                    interval=Interval(3, 12),
                    spec_ids=("b8_kozak", "d1_restriction_sites"),
                    kind="mutually_exclusive",
                    binding_spec_id="d1_restriction_sites",
                ),
            ),
            provenance=Provenance(
                app_version="0.1.0", seed=7, engine_versions={"viennarna": "2.7.2"}
            ),
        )
        return build_report(
            result,
            cand,
            translation_table_id=11,
            preset_id="ecoli_expression",
            advisories=("uAUG at 812 in your 5'UTR: no codon can move it",),
            strain_protocol=("Stbl3 at 30 C; covers the long repeats only",),
            **kw,  # type: ignore[arg-type]
        )

    def test_prints_the_genetic_code(self) -> None:
        assert "translation table 11" in render(self.report())

    def test_names_what_was_not_evaluated(self) -> None:
        text = render(self.report())
        assert "NOT evaluated" in text
        assert "no folding engine installed" in text

    def test_says_the_percentile_is_not_a_prediction(self) -> None:
        text = render(self.report())
        assert "not a prediction" in text
        assert "never a predicted expression level" in text

    def test_uses_no_banned_prediction_vocabulary(self) -> None:
        """The same grep CI runs over engine source, applied to what the report
        actually EMITS -- which is where a user would read it."""
        banned = re.compile(
            r"\b(predicted_expression|expression_score|titer_prediction|predicted_titer|"
            r"fold_improvement|expression_level|predicted_yield)\b"
        )
        text = render(self.report())
        offending = [line for line in text.splitlines() if banned.search(line)]
        # "never a predicted expression level" is the disclaimer itself.
        assert all("never a predicted" in line for line in offending), offending

    def test_states_the_native_baseline(self) -> None:
        assert "native sequence" in render(self.report())

    def test_shows_the_conflict_with_the_binding_rule(self) -> None:
        text = render(self.report())
        assert "b8_kozak" in text
        assert "binding: d1_restriction_sites" in text

    def test_shows_the_advisory_apart_from_the_conflicts(self) -> None:
        text = render(self.report())
        assert "Not fixable by codon choice" in text
        assert "uAUG at 812" in text

    def test_the_strain_protocol_keeps_its_caveat(self) -> None:
        assert "covers the long repeats only" in render(self.report())

    def test_the_design_hash_is_on_the_report(self) -> None:
        assert "abc123def456" in render(self.report())
