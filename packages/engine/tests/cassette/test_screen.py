"""The biosecurity screen, and the one property it exists to guarantee.

The load-bearing tests are `TestNeverClaimsClear`. Everything else checks
behaviour; those check the negative -- that no path which did not obtain an
explicit commec "Clear" ever reports `status="clear"`. That is the bug this
whole module exists to prevent: a screen that reads as clean when it never ran
passes every type check and every green test, and hands a user a false
assurance about a hazard.

commec is absent in CI (the bootstrap install is `dev,fold,export`, not
`screen`), so the degraded path -- `load_screen` returning a `NullScreen` --
is the one actually exercised here. The commec-backed paths are driven through
an injected fake `runner`, never a real commec process.
"""

from __future__ import annotations

import subprocess

import pytest
from bt5.cassette import screen as screen_mod
from bt5.cassette.screen import (
    BiosecurityBlockedError,
    CommecOutcome,
    CommecScreen,
    NullScreen,
    _status_for,
    commec_available,
    guard_emission,
    load_screen,
    screen_degradation_reason,
)
from bt5.core.context import BiosecurityVerdict
from hypothesis import given
from hypothesis import strategies as st

PROTEIN = "MKVLA"


def _runner_returning(outcome: str, version: str | None = "commec-db-2026.01"):
    def _run(protein: str, database_dir, timeout_s: float) -> CommecOutcome:
        return CommecOutcome(outcome=outcome, database_version=version)

    return _run


def _runner_raising(exc: Exception):
    def _run(protein: str, database_dir, timeout_s: float) -> CommecOutcome:
        raise exc

    return _run


def _commec(runner, *, timeout_s: float = 900.0) -> CommecScreen:
    return CommecScreen(database_dir="/nonexistent/db", timeout_s=timeout_s, runner=runner)


class TestStatusMapping:
    """commec's Clear / Warning / Flag mapped to BT5's status vocabulary."""

    @pytest.mark.parametrize(
        ("outcome", "status"),
        [
            ("Clear", "clear"),
            ("clear", "clear"),
            ("Warning", "flag"),
            ("warn", "flag"),
            ("Flag", "block"),
            ("FLAG", "block"),
        ],
    )
    def test_known_outcomes(self, outcome: str, status: str) -> None:
        assert _status_for(outcome) == status

    @pytest.mark.parametrize(
        "outcome", ["", "   ", "error", "skipped", "cleared_flag", "pass", "unknown"]
    )
    def test_unrecognised_outcome_is_not_run(self, outcome: str) -> None:
        assert _status_for(outcome) == "not_run"


class TestNullScreen:
    def test_reports_not_run(self) -> None:
        verdict = NullScreen().screen(PROTEIN)
        assert verdict.status == "not_run"
        assert verdict.database_version is None
        assert verdict.may_proceed is True  # not_run does not itself refuse

    def test_carries_its_reason(self) -> None:
        verdict = NullScreen(reason="because reasons").screen(PROTEIN)
        assert verdict.detail == "because reasons"

    def test_empty_protein_raises(self) -> None:
        with pytest.raises(ValueError, match="empty protein"):
            NullScreen().screen("")


class TestCommecScreen:
    def test_clear_records_database_version(self) -> None:
        verdict = _commec(_runner_returning("Clear")).screen(PROTEIN)
        assert verdict.status == "clear"
        assert verdict.database_version == "commec-db-2026.01"

    def test_warning_flags_but_proceeds(self) -> None:
        verdict = _commec(_runner_returning("Warning")).screen(PROTEIN)
        assert verdict.status == "flag"
        assert verdict.may_proceed is True
        assert verdict.database_version == "commec-db-2026.01"

    def test_flag_blocks(self) -> None:
        verdict = _commec(_runner_returning("Flag")).screen(PROTEIN)
        assert verdict.status == "block"
        assert verdict.may_proceed is False
        assert verdict.database_version == "commec-db-2026.01"

    def test_real_verdict_always_carries_database_version(self) -> None:
        for outcome in ("Clear", "Warning", "Flag"):
            verdict = _commec(_runner_returning(outcome, "v-42")).screen(PROTEIN)
            assert verdict.database_version == "v-42"

    def test_empty_protein_raises(self) -> None:
        with pytest.raises(ValueError, match="empty protein"):
            _commec(_runner_returning("Clear")).screen("")


class TestDegradationNeverLies:
    """Every way commec can fail becomes not_run -- never clear."""

    def test_timeout_is_not_run(self) -> None:
        runner = _runner_raising(subprocess.TimeoutExpired(cmd="commec", timeout=900.0))
        verdict = _commec(runner).screen(PROTEIN)
        assert verdict.status == "not_run"
        assert "unfinished" in verdict.detail

    def test_nonzero_exit_is_not_run(self) -> None:
        runner = _runner_raising(subprocess.CalledProcessError(1, "commec"))
        verdict = _commec(runner).screen(PROTEIN)
        assert verdict.status == "not_run"

    def test_missing_binary_is_not_run(self) -> None:
        runner = _runner_raising(screen_mod.ScreenUnavailableError("commec not on PATH"))
        verdict = _commec(runner).screen(PROTEIN)
        assert verdict.status == "not_run"
        assert "PATH" in verdict.detail

    def test_os_error_is_not_run(self) -> None:
        verdict = _commec(_runner_raising(OSError("boom"))).screen(PROTEIN)
        assert verdict.status == "not_run"

    def test_unparseable_result_is_not_run(self) -> None:
        verdict = _commec(_runner_raising(ValueError("no recommendation"))).screen(PROTEIN)
        assert verdict.status == "not_run"

    def test_unrecognised_outcome_is_not_run_not_clear(self) -> None:
        verdict = _commec(_runner_returning("Errored")).screen(PROTEIN)
        assert verdict.status == "not_run"
        assert "does not recognise" in verdict.detail


class TestNeverClaimsClear:
    """The one property this module exists to guarantee."""

    def test_absent_commec_load_screen_is_not_clear(self, monkeypatch) -> None:
        # The CI-exercised path: commec is genuinely absent in the bootstrap
        # install, so this exercises the real degradation without mocking.
        assert commec_available() is False
        verdict = load_screen().screen(PROTEIN)
        assert verdict.status != "clear"
        assert verdict.status == "not_run"

    def test_commec_present_but_no_database_is_not_clear(self, monkeypatch) -> None:
        monkeypatch.setattr(screen_mod, "commec_available", lambda: True)
        verdict = load_screen(database_dir=None).screen(PROTEIN)
        assert verdict.status != "clear"
        assert verdict.status == "not_run"

    @given(st.text())
    def test_only_explicit_clear_reads_clear(self, outcome: str) -> None:
        # No arbitrary outcome word may map to "clear" except commec's own
        # "Clear" (case- and whitespace-insensitive). This is the invariant the
        # default-to-not_run in `_status_for` enforces.
        if _status_for(outcome) == "clear":
            assert outcome.strip().lower() == "clear"

    @given(
        st.sampled_from(
            [
                subprocess.TimeoutExpired(cmd="commec", timeout=1.0),
                subprocess.CalledProcessError(1, "commec"),
                OSError("io"),
                ValueError("parse"),
                screen_mod.ScreenUnavailableError("gone"),
            ]
        )
    )
    def test_no_failure_mode_reads_clear(self, exc: Exception) -> None:
        verdict = _commec(_runner_raising(exc)).screen(PROTEIN)
        assert verdict.status != "clear"


class TestLoadScreen:
    def test_absent_commec_returns_null_screen(self) -> None:
        assert isinstance(load_screen(), NullScreen)

    def test_present_commec_with_database_returns_commec_screen(self, monkeypatch) -> None:
        monkeypatch.setattr(screen_mod, "commec_available", lambda: True)
        assert isinstance(load_screen(database_dir="/some/db"), CommecScreen)

    def test_present_commec_without_database_returns_null_screen(self, monkeypatch) -> None:
        monkeypatch.setattr(screen_mod, "commec_available", lambda: True)
        assert isinstance(load_screen(database_dir=None), NullScreen)


class TestDegradationReason:
    def test_absent_commec_reports_a_reason(self) -> None:
        reason = screen_degradation_reason()
        assert reason is not None
        assert "commec is not installed" in reason

    def test_present_commec_no_database_reports_a_reason(self, monkeypatch) -> None:
        monkeypatch.setattr(screen_mod, "commec_available", lambda: True)
        reason = screen_degradation_reason(database_dir=None)
        assert reason is not None
        assert "no reference database" in reason

    def test_present_commec_with_database_is_none(self, monkeypatch) -> None:
        monkeypatch.setattr(screen_mod, "commec_available", lambda: True)
        assert screen_degradation_reason(database_dir="/some/db") is None


class TestGuardEmission:
    def test_block_refuses(self) -> None:
        verdict = BiosecurityVerdict("block", "v1", "regulated pathogen match")
        with pytest.raises(BiosecurityBlockedError) as excinfo:
            guard_emission(verdict)
        assert excinfo.value.verdict is verdict

    @pytest.mark.parametrize("status", ["clear", "flag", "not_run"])
    def test_non_block_proceeds(self, status: str) -> None:
        guard_emission(BiosecurityVerdict(status, None, ""))  # does not raise

    def test_guard_never_downgrades_a_block(self) -> None:
        # There is deliberately no argument that turns a block into a pass.
        verdict = BiosecurityVerdict("block", None, "")
        with pytest.raises(BiosecurityBlockedError):
            guard_emission(verdict)
