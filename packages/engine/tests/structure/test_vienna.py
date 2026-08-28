"""The folding engine, and the provenance that makes its numbers comparable.

NEVER a byte-exact assertion on a dG. Energy parameters are floating point and
version-dependent, and a snapshot of one would turn a ViennaRNA upgrade into a
diff nobody can read instead of the scientific change it is. Every energy here
is compared with a tolerance, or against another energy from the same run.

The availability tests run everywhere. The folding tests need the `fold` extra,
which CI does not currently install -- see the PR for that.
"""

from __future__ import annotations

import pytest
from bt5.core.services import FoldEnergy
from bt5.core.types import Interval, reverse_complement
from bt5.structure import (
    CALIBRATED_VERSION,
    ENGINE_NAME,
    PARAM_SET,
    ViennaFold,
    degradation_reason,
    installed_version,
    is_calibrated,
    load_fold_engine,
    slice_of,
    vienna_available,
)

SEQ = "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCAC"


class TestAvailabilityAndDegradation:
    """These run with or without ViennaRNA -- the honest-degradation path is the
    half that MUST work when the engine is missing, so it cannot be skipped
    alongside the engine it stands in for."""

    def test_availability_and_the_engine_factory_agree(self) -> None:
        assert (load_fold_engine() is not None) == vienna_available()

    def test_an_absent_engine_yields_none_not_a_stub(self) -> None:
        """A stub returning plausible kcal/mol would flow through the scorers,
        the null and the percentile unchallenged and emerge as a confident rank."""
        engine = load_fold_engine()
        if engine is None:
            assert degradation_reason() is not None
            assert "not installed" in degradation_reason()  # type: ignore[operator]

    def test_a_degradation_is_reported_whenever_the_result_is_not_calibrated(self) -> None:
        reason = degradation_reason()
        if vienna_available() and is_calibrated():
            assert reason is None
        else:
            assert reason, "an uncalibrated or absent engine must say so"

    def test_the_calibrated_version_is_the_one_pinned_in_pyproject(self) -> None:
        """CLAUDE.md section 6: a version bump is a scientific change. This is
        the constant a benchmark comparability guard reads."""
        import pathlib

        text = pathlib.Path("pyproject.toml").read_text()
        assert f'viennarna=={CALIBRATED_VERSION}"' in text

    def test_the_calibration_key_matches_the_spec_field_form(self) -> None:
        """Spec.engine_calibration documents 'viennarna:rna_turner2004'; a rule
        comparing its declared calibration against a result's key must match."""
        energy = FoldEnergy(-1.0, ENGINE_NAME, "2.7.2", PARAM_SET)
        assert energy.calibration_key == "viennarna:rna_turner2004"


class TestSliceOf:
    """Window extraction, tested without folding: a wrong slice is a wrong
    number that looks entirely reasonable."""

    def test_a_plain_window(self) -> None:
        assert slice_of(SEQ, Interval(3, 9)) == SEQ[3:9]

    def test_a_wrapping_window(self) -> None:
        w = slice_of(SEQ, Interval(len(SEQ) - 5, len(SEQ) + 4))
        assert w == SEQ[-5:] + SEQ[:4]
        assert len(w) == 9

    def test_a_reverse_strand_window_is_reverse_complemented(self) -> None:
        assert slice_of(SEQ, Interval(3, 9, -1)) == reverse_complement(SEQ[3:9])

    def test_a_wrapping_reverse_window(self) -> None:
        iv = Interval(len(SEQ) - 5, len(SEQ) + 4, -1)
        assert slice_of(SEQ, iv) == reverse_complement(SEQ[-5:] + SEQ[:4])


class TestFolding:
    @pytest.fixture(autouse=True)
    def _needs_vienna(self) -> None:
        pytest.importorskip("RNA", reason="the `fold` extra is not installed")

    def engine(self, **kw: float) -> ViennaFold:
        return ViennaFold(**kw)  # type: ignore[arg-type]

    def test_an_energy_carries_everything_needed_to_compare_it(self) -> None:
        e = self.engine().mfe_window(SEQ, Interval(0, 40))
        assert e.engine == ENGINE_NAME
        assert e.param_set == PARAM_SET
        assert e.engine_version == installed_version()
        assert e.temperature_c == 37.0
        assert e.dangles == 2
        assert e.dg_kcal_mol < 0

    def test_temperature_changes_the_answer(self) -> None:
        """Measured 23-33% swings. A dG without its temperature is not a number."""
        cold = self.engine(temperature_c=30.0).mfe_window(SEQ, Interval(0, 60)).dg_kcal_mol
        warm = self.engine(temperature_c=42.0).mfe_window(SEQ, Interval(0, 60)).dg_kcal_mol
        assert cold < warm, "colder folds more stably"
        assert abs(cold - warm) > 0.5

    def test_dangles_changes_the_answer_and_is_recorded(self) -> None:
        a = self.engine(dangles=0).mfe_window(SEQ, Interval(0, 60))
        b = self.engine(dangles=2).mfe_window(SEQ, Interval(0, 60))
        assert a.dangles == 0
        assert b.dangles == 2
        assert a.dg_kcal_mol != pytest.approx(b.dg_kcal_mol, abs=1e-6)

    def test_dna_and_rna_spelling_agree(self) -> None:
        """ViennaRNA reads T as U today. This test is what catches a future
        version that stops -- otherwise it would surface as a shifted baseline."""
        dna = self.engine().mfe_window(SEQ, Interval(0, 60)).dg_kcal_mol
        rna = self.engine().mfe_window(SEQ.replace("T", "U"), Interval(0, 60)).dg_kcal_mol
        assert dna == pytest.approx(rna, abs=1e-6)

    def test_a_window_folds_exactly_the_sequence_the_slice_names(self) -> None:
        iv = Interval(10, 55)
        assert self.engine().mfe_window(SEQ, iv).dg_kcal_mol == pytest.approx(
            self.engine().mfe(slice_of(SEQ, iv)).dg_kcal_mol, abs=1e-6
        )

    def test_a_wrapping_window_folds_across_the_origin(self) -> None:
        iv = Interval(len(SEQ) - 20, len(SEQ) + 20)
        got = self.engine().mfe_window(SEQ, iv).dg_kcal_mol
        expected = self.engine().mfe(SEQ[-20:] + SEQ[:20]).dg_kcal_mol
        assert got == pytest.approx(expected, abs=1e-6)

    def test_a_reverse_window_folds_the_reverse_complement(self) -> None:
        iv = Interval(10, 55, -1)
        got = self.engine().mfe_window(SEQ, iv).dg_kcal_mol
        expected = self.engine().mfe(reverse_complement(SEQ[10:55])).dg_kcal_mol
        assert got == pytest.approx(expected, abs=1e-6)

    def test_a_hairpin_folds_more_stably_than_random_sequence(self) -> None:
        """A sanity check with a known sign, not a magnitude assertion."""
        stem = "GGGCCCAAAGGGCCC"
        hairpin = stem + "TTTT" + reverse_complement(stem)
        loose = "ATATATATATATATATATATATATATATATATAT"
        assert self.engine().mfe(hairpin).dg_kcal_mol < self.engine().mfe(loose).dg_kcal_mol

    def test_non_acgu_input_is_refused(self) -> None:
        with pytest.raises(ValueError, match="non-ACGU"):
            self.engine().mfe("ACGTNNNN")

    def test_accessibility_is_a_probability(self) -> None:
        p = self.engine().accessibility(SEQ, Interval(0, 60), 5)
        assert p is not None
        assert 0.0 <= p <= 1.0

    def test_accessibility_is_none_when_the_window_is_shorter_than_the_stretch(self) -> None:
        """A question with no answer, rather than a zero."""
        assert self.engine().accessibility(SEQ, Interval(0, 4), 10) is None

    def test_a_longer_stretch_is_less_likely_to_be_unpaired(self) -> None:
        short = self.engine().accessibility(SEQ, Interval(0, 70), 1)
        long = self.engine().accessibility(SEQ, Interval(0, 70), 10)
        assert short is not None
        assert long is not None
        assert long <= short

    def test_a_zero_length_stretch_is_refused(self) -> None:
        with pytest.raises(ValueError, match="must be positive"):
            self.engine().accessibility(SEQ, Interval(0, 60), 0)


class TestCost:
    """The windowed fold is the null-model and interactive-loop primitive, so
    its cost is a property worth pinning. Whole-sequence MFE is O(n^3) and
    measured SLOWER than docs/PLAN.md states -- 0.77 s at 1 kb against the
    plan's ~0.24 s -- which is exactly why `mfe` is report-time only.

    Generous bounds: this catches an accidental whole-sequence fold in a loop,
    not milliseconds on a shared runner.
    """

    @pytest.fixture(autouse=True)
    def _needs_vienna(self) -> None:
        pytest.importorskip("RNA", reason="the `fold` extra is not installed")

    def test_a_kudla_window_is_fast_enough_for_a_two_hundred_variant_null(self) -> None:
        import time

        engine = ViennaFold()
        seq = SEQ * 3
        start = time.perf_counter()
        for _ in range(200):
            engine.mfe_window(seq, Interval(0, 41))
        elapsed = time.perf_counter() - start
        assert elapsed < 5.0, (
            f"200 windowed folds took {elapsed:.2f}s; gate G3 allows 2s for a "
            f"200-variant null, and a full sweep instead of a 5' window costs ~148s"
        )
