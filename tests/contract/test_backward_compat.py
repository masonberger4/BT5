"""Every recorded fixture must still build against today's contract.

docs/PLAN.md's amendment protocol names this test by name: a MAJOR amendment
ships a deprecation shim and `test_backward_compat` re-parsing every recorded
fixture inside the amendment PR. That is what makes the two-window rule real --
the shim is not "written", it is "written and demonstrated on the values that
existed before it".

The manifest and these fixtures answer different questions and neither
subsumes the other. The manifest asks "did the declared shape change?" and can
be regenerated into agreement with anything. A fixture asks "does a value
recorded under the old contract still construct?", and regenerating it does not
make an old caller work again -- it only stops recording that it doesn't. That
is why the fixtures are checked for staleness too.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import fixtures  # noqa: E402

RECORDED = fixtures.load_all()
REGENERATE = "python tests/contract/regenerate.py"


def test_there_are_fixtures_at_all() -> None:
    """A silently empty fixture directory turns every test below into a no-op."""
    assert RECORDED, f"no recorded fixtures; run {REGENERATE}"
    assert len(RECORDED) >= 15, f"only {len(RECORDED)} fixtures; the contract has more shapes"


@pytest.mark.parametrize("name", sorted(RECORDED))
def test_a_recorded_value_still_constructs(name: str) -> None:
    """The whole point. Fails on a renamed field, a removed field, or a new
    required one -- the three changes that break every existing caller -- and
    passes on a new defaulted field, which breaks none."""
    record = RECORDED[name]
    try:
        revived = fixtures.decode(record["value"])
    except TypeError as exc:
        pytest.fail(
            f"fixture {name!r} ({record['recorded_as']}) no longer constructs: {exc}\n"
            f"This is a MAJOR contract change. It needs a deprecation shim that "
            f"accepts the recorded form, not a regenerated fixture -- regenerating "
            f"stops recording the breakage without fixing it."
        )
    assert type(revived).__module__ + "." + type(revived).__qualname__ == record["recorded_as"]


@pytest.mark.parametrize("name", sorted(RECORDED))
def test_a_recorded_value_round_trips_unchanged(name: str) -> None:
    """Construct, re-encode, compare. Catches a field that still ACCEPTS the
    recorded value but no longer stores it -- a __post_init__ that normalises,
    a default that silently wins over what was passed."""
    record = RECORDED[name]
    assert fixtures.encode(fixtures.decode(record["value"])) == record["value"]


def test_the_fixtures_are_not_stale() -> None:
    """Recorded fixtures must match what the recorder produces today.

    Without this, a fixture could be edited by hand into agreement with
    whatever broke it, which is the same failure mode CLAUDE.md bans for
    snapshots: `--snapshot-update` is not a fix.
    """
    live = {name: fixtures.encode(value) for name, value in fixtures.specimens().items()}
    recorded = {name: record["value"] for name, record in RECORDED.items()}
    assert live.keys() == recorded.keys(), (
        f"fixture set changed; run {REGENERATE}. "
        f"Missing: {sorted(live.keys() - recorded.keys())}; "
        f"orphaned: {sorted(recorded.keys() - live.keys())}"
    )
    differing = sorted(k for k in live if live[k] != recorded[k])
    assert not differing, (
        f"recorded fixtures differ from live values: {differing}. Run {REGENERATE}"
    )


def test_the_codec_carries_nan_rather_than_rounding_it_away() -> None:
    """ObjectiveScore.unavailable() fills raw and percentile with NaN
    deliberately: 0.0 would read as a real, terrible score. JSON has no NaN, so
    a codec that quietly dropped it would record the exact lie the contract
    added that field to prevent."""
    import math

    from bt5.core.result import ObjectiveScore

    revived = fixtures.decode(RECORDED["objective_unavailable"]["value"])
    assert isinstance(revived, ObjectiveScore)
    assert not revived.available
    assert math.isnan(revived.raw)
    assert math.isnan(revived.percentile)


def test_a_wrapping_interval_survives_the_round_trip() -> None:
    """end > length is BT5's one representation of an origin wrap. A codec that
    normalised it would erase the case the coordinate model exists for."""
    from bt5.core.types import Construct

    revived = fixtures.decode(RECORDED["construct_circular"]["value"])
    assert isinstance(revived, Construct)
    cds = revived.editable[0]
    assert cds.wraps(revived.length), "the recorded CDS spans the origin"
    assert revived.is_editable(cds)
