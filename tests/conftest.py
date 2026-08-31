from __future__ import annotations

import os
import sys
from pathlib import Path

from hypothesis import HealthCheck, settings

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "packages" / "engine" / "src"))
sys.path.insert(0, str(ROOT / "tests"))

settings.register_profile(
    "ci", max_examples=200, deadline=None, suppress_health_check=[HealthCheck.too_slow]
)
settings.register_profile("dev", max_examples=50, deadline=None)
settings.register_profile(
    "nightly", max_examples=2000, deadline=None, suppress_health_check=[HealthCheck.too_slow]
)
#: `dev` locally, `ci` in CI, `nightly` when asked for. Read from the
#: environment, because THE ENVIRONMENT VARIABLE ALONE DOES NOTHING: Hypothesis
#: has no HYPOTHESIS_PROFILE support of its own -- its only env-driven selection
#: is the built-in `ci` profile via `is_in_ci()`, and this file overwrites that
#: name and then unconditionally loaded `dev` on top. So `invariants`, which
#: docs/PLAN.md calls the single most valuable gate, ran at 50 examples instead
#: of 200 for its whole life, and `nightly` was unreachable code.
#:
#: `--hypothesis-profile` on the pytest command line still wins: the plugin
#: applies it after conftest import.
_PROFILE = os.environ.get("HYPOTHESIS_PROFILE", "dev")
if _PROFILE not in ("dev", "ci", "nightly"):
    raise RuntimeError(
        f"HYPOTHESIS_PROFILE={_PROFILE!r} is not a registered profile; "
        f"expected one of dev, ci, nightly"
    )
settings.load_profile(_PROFILE)
