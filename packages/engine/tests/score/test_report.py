"""The QC report, tested mostly for what it refuses to say.

BT5's differentiator is trustworthiness, and the report is where that is either
kept or lost: it is the last artefact a human reads before spending money on
DNA. So the assertions here are largely negative -- no predicted number, no
silently dropped objective, no defaulted genetic code.
"""

from __future__ import annotations

import math
import re

import pytest
from bt5.core.result import Candidate, Conflict, DesignResult, ObjectiveScore, ScoreCard
from bt5.core.types import Construct, Interval, Provenance, Segment, SegmentKind, Topology
from bt5.rules.vendors import DEFAULT_VENDOR, PROFILES, orderable_keys
from bt5.score.report import (
    ERROR_FREE_BP,
    QcReport,
    build_report,
    render,
    screening_burden,
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


class TestScreeningBurden:
    def test_matches_the_published_formula(self) -> None:
        b = screening_burden(1500, vendor="idt_eblocks")
        assert b.error_free_bp == ERROR_FREE_BP["idt_eblocks"]
        assert b.p_perfect == pytest.approx(math.exp(-1500 / 5000))

    def test_colonies_reach_the_requested_confidence(self) -> None:
        """The number is only useful if picking that many actually gets you there."""
        for length in (300, 1500, 3000, 6000):
            b = screening_burden(length, vendor="idt_eblocks")
            reached = 1.0 - (1.0 - b.p_perfect) ** b.colonies_to_pick
            assert reached >= b.confidence
            one_fewer = 1.0 - (1.0 - b.p_perfect) ** (b.colonies_to_pick - 1)
            if b.colonies_to_pick > 1:
                assert one_fewer < b.confidence, "and it must not be more than needed"

    def test_a_longer_construct_costs_more_picking(self) -> None:
        short = screening_burden(1000, vendor="idt_eblocks")
        long = screening_burden(4000, vendor="idt_eblocks")
        assert long.colonies_to_pick > short.colonies_to_pick
        assert long.p_perfect < short.p_perfect

    def test_the_vendor_changes_the_answer(self) -> None:
        """Twist's error-free length is longer, so the same insert needs less
        screening. Reporting one vendor's number under another's name is the
        kind of quiet error a lab pays for in plates."""
        twist = screening_burden(6000, vendor="twist_gene_fragment")
        eblocks = screening_burden(6000, vendor="idt_eblocks")
        assert twist is not None
        assert eblocks is not None
        assert twist.colonies_to_pick < eblocks.colonies_to_pick

    def test_an_unknown_vendor_is_refused_rather_than_defaulted(self) -> None:
        with pytest.raises(ValueError, match="unknown vendor 'acme'"):
            screening_burden(1000, vendor="acme")

    def test_a_real_vendor_with_no_figure_on_file_returns_none_rather_than_raising(
        self,
    ) -> None:
        """Two different failures, and collapsing them is the bug this guards.

        "acme" is nobody's product and that is a caller error. `idt_gblocks` IS
        orderable -- it is BT5's own default -- and simply has no published
        error-free length in this repo. Raising for it would make a data gap
        look like a bad request; answering with eBlocks' 5000 would be worse.
        """
        assert "idt_gblocks" in PROFILES
        assert "idt_gblocks" not in ERROR_FREE_BP
        assert screening_burden(1000, vendor="idt_gblocks") is None

    def test_no_vendor_chosen_is_refused_rather_than_answered_with_a_shrug(self) -> None:
        """`none` is not a vendor whose fidelity we happen to lack.

        How many colonies to pick is a question only a vendor can answer, so
        returning None for `none` would read as "this product has no figure on
        file" when the truth is that nobody picked a product. Same refusal E1
        and E9 make, for the same reason.
        """
        with pytest.raises(ValueError, match="not orderable from anyone"):
            screening_burden(1000, vendor="none")

    def test_no_fourth_vendor_namespace(self) -> None:
        """Every key here must name a real orderable configuration.

        This table used to be keyed on "twist" and "idt_eblocks", overlapping
        the real registry by one key and by accident, so screening_burden raised
        for BT5's own DEFAULT_VENDOR and accepted "twist", which names nothing
        anyone can order. That is the third instance of the bug PR #53 fixed
        twice. See issue #54.
        """
        assert set(ERROR_FREE_BP) <= set(orderable_keys())

    def test_the_vendor_numbers_carry_a_verification_date(self) -> None:
        """Vendor figures drift -- Twist moved its homopolymer limit from 14 to
        30 bp between 2023 and 2026 -- so an undated one is unmaintainable."""
        b = screening_burden(1000, vendor="idt_eblocks")
        assert b is not None
        assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", b.last_verified)

    def test_rejects_nonsense_inputs(self) -> None:
        with pytest.raises(ValueError, match="length must be positive"):
            screening_burden(0)
        with pytest.raises(ValueError, match="confidence must be in"):
            screening_burden(1000, confidence=1.0)


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
        """Vendor named explicitly: completeness now includes the screening line,
        and BT5's default vendor has no error-free length on file."""
        result = self.result()
        report = build_report(
            result, result.candidates[0], translation_table_id=1, vendor="idt_eblocks"
        )
        assert report.is_complete

    def test_a_vendor_with_no_error_free_length_degrades_rather_than_going_quiet(
        self,
    ) -> None:
        """The default vendor is gBlocks and BT5 has no fidelity figure for it.

        The report must SAY that, not merely omit the line -- a report missing
        its screening burden looks exactly like one where nobody asked for it.
        So `burden` is None, a degradation names the vendor, and `is_complete`
        is False until the number exists. That is deliberate pressure: the
        alternative is answering with eBlocks' number under gBlocks' name.
        """
        result = self.result()
        report = build_report(result, result.candidates[0], translation_table_id=1)
        assert report.burden is None
        assert not report.is_complete
        assert any("no published error-free length" in d for d in report.degradations)
        assert any(DEFAULT_VENDOR in d for d in report.degradations)


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

    def test_shows_the_screening_burden_as_an_instruction(self) -> None:
        text = render(self.report())
        assert "pick" in text
        assert "colonies" in text

    def test_the_strain_protocol_keeps_its_caveat(self) -> None:
        assert "covers the long repeats only" in render(self.report())

    def test_the_design_hash_is_on_the_report(self) -> None:
        assert "abc123def456" in render(self.report())
