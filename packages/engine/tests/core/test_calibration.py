"""`engine_calibration` had to start meaning something.

It is a required ClassVar on every Spec, and until now nothing read it. Both
reference rules set it to None, so the first rule to carry a ViennaRNA-calibrated
-39 kcal/mol threshold would have merged green against whatever engine happened
to be installed -- which the research brief names as the single most likely
correctness bug in the folding feature. It is silent by construction: the
comparison succeeds, the rule fires or does not, and nothing looks wrong.
"""

from __future__ import annotations

import pytest
from bt5.core.registry import (
    CalibrationMismatchError,
    _validate,
    check_engine_calibration,
)
from bt5.core.spec import Direction, Enforcement, Evidence


class FakeEngine:
    name = "viennarna"
    param_set = "rna_turner2004"


class OtherEngine:
    name = "seqfold"
    param_set = "default"


def rule(rule_id: str, calibration: str | None):
    return type(
        f"Rule_{rule_id}",
        (),
        {
            "id": rule_id,
            "engine_calibration": calibration,
            "enforcement": Enforcement.SOFT,
            "evidence": Evidence.EVIDENCE_BACKED,
            "direction": Direction.LOWER_IS_BETTER,
        },
    )


class TestCalibrationGate:
    def test_an_uncalibrated_rule_runs_against_anything(self) -> None:
        """A motif or GC rule does not depend on an energy model at all."""
        specs = [rule("d1", None)]
        assert check_engine_calibration(specs, FakeEngine()) == tuple(specs)  # type: ignore[arg-type]

    def test_a_matching_rule_runs(self) -> None:
        specs = [rule("b1", "viennarna:rna_turner2004")]
        assert check_engine_calibration(specs, FakeEngine()) == tuple(specs)  # type: ignore[arg-type]

    def test_a_mismatched_rule_raises_rather_than_being_skipped(self) -> None:
        """Skipping would drop a constraint silently. This is a configuration
        error, and the run should stop rather than quietly enforce less."""
        specs = [rule("b1", "viennarna:rna_turner2004")]
        with pytest.raises(CalibrationMismatchError, match="calibrated"):
            check_engine_calibration(specs, OtherEngine())  # type: ignore[arg-type]

    def test_the_error_names_the_rule_and_both_calibrations(self) -> None:
        specs = [rule("b1", "viennarna:rna_turner2004")]
        with pytest.raises(CalibrationMismatchError) as exc:
            check_engine_calibration(specs, OtherEngine())  # type: ignore[arg-type]
        assert "b1" in str(exc.value)
        assert "seqfold:default" in str(exc.value)
        assert "viennarna:rna_turner2004" in str(exc.value)

    def test_with_no_engine_a_calibrated_rule_is_unrunnable_not_an_error(self) -> None:
        """Absence is a degradation to report, not a misconfiguration to crash
        on -- see ObjectiveScore.unavailable."""
        specs = [rule("b1", "viennarna:rna_turner2004"), rule("d1", None)]
        runnable = check_engine_calibration(specs, None)  # type: ignore[arg-type]
        assert [c.id for c in runnable] == ["d1"]  # type: ignore[attr-defined]

    def test_a_real_engine_satisfies_a_real_rule_calibration(self) -> None:
        """End to end against the shipped engine rather than a fake, so the two
        halves cannot drift into agreeing only with each other."""
        pytest.importorskip("RNA", reason="the `fold` extra is not installed")
        from bt5.structure import ViennaFold

        specs = [rule("b1", "viennarna:rna_turner2004")]
        assert check_engine_calibration(specs, ViennaFold()) == tuple(specs)  # type: ignore[arg-type]


class TestWellFormedAtImportTime:
    """A typo in a calibration string matches nothing, so the rule is never
    refused -- exactly the failure the gate exists to prevent. Caught where it is
    cheapest, at registration.

    `_validate` is called directly rather than through `register`, which would
    mutate the process-wide registry. An earlier version of this file cleared
    and restored it and broke four unrelated tests; the registry is shared
    state, and a test that reaches for it has to put it back exactly.
    """

    def probe(self, calibration: str | None):
        cls = rule("zz_probe", calibration)
        for attr, value in (
            ("citations", ("x",)),
            ("last_verified", "2026-08-28"),
            ("weight_provenance", "test"),
        ):
            setattr(cls, attr, value)
        return cls

    @pytest.mark.parametrize("bad", ["viennarna", "viennarna:", ":turner", "vienna rna:t", ""])
    def test_a_malformed_calibration_is_refused(self, bad: str) -> None:
        with pytest.raises(ValueError, match="engine:param_set"):
            _validate(self.probe(bad))  # type: ignore[arg-type]

    def test_a_well_formed_calibration_passes(self) -> None:
        _validate(self.probe("viennarna:rna_turner2004"))  # type: ignore[arg-type]

    def test_no_calibration_at_all_is_allowed(self) -> None:
        """Most rules need no engine, and requiring a value would force every
        motif rule to declare a folding calibration it does not use."""
        _validate(self.probe(None))  # type: ignore[arg-type]
