from __future__ import annotations

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
settings.load_profile("dev")
