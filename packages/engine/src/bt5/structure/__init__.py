"""M6 structure lane: RNA secondary structure behind the `FoldEngine` protocol.

Two halves, deliberately separable. `windows` decides WHERE to fold and holds no
dependency on ViennaRNA, so the arithmetic that actually goes wrong -- origin
wraps, strand, whether a window fits -- is tested without an engine present.
`vienna` produces the numbers, and never produces one without the temperature,
dangles, engine version and parameter set that make it comparable to anything.
"""

from bt5.structure.vienna import (
    ACCESSIBILITY_MAX_BP_SPAN,
    ACCESSIBILITY_WINDOW,
    CALIBRATED_VERSION,
    ENGINE_NAME,
    PARAM_SET,
    FoldUnavailableError,
    ViennaFold,
    degradation_reason,
    installed_version,
    is_calibrated,
    load_fold_engine,
    slice_of,
    vienna_available,
)
from bt5.structure.windows import (
    KUDLA_DOWNSTREAM,
    KUDLA_UPSTREAM,
    SWEEP_SIZE,
    SWEEP_STEP,
    five_prime_window,
    sliding_windows,
    windows_touching,
)

__all__ = [
    "ACCESSIBILITY_MAX_BP_SPAN",
    "ACCESSIBILITY_WINDOW",
    "CALIBRATED_VERSION",
    "ENGINE_NAME",
    "KUDLA_DOWNSTREAM",
    "KUDLA_UPSTREAM",
    "PARAM_SET",
    "SWEEP_SIZE",
    "SWEEP_STEP",
    "FoldUnavailableError",
    "ViennaFold",
    "degradation_reason",
    "five_prime_window",
    "installed_version",
    "is_calibrated",
    "load_fold_engine",
    "slice_of",
    "sliding_windows",
    "vienna_available",
    "windows_touching",
]
