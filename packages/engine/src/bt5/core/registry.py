"""Rule autodiscovery.

There is deliberately NO committed catalog file and NO hand-maintained
`RULES = [...]` list. Both were tried in competing designs and both collide on
every pull request at the highest-volume lane -- entries sort by rule id, so hunks
interleave and every merge invalidates every other open rules PR.

Instead the registry walks the package tree and collects `@register`-decorated
classes. Adding a rule edits ZERO shared files, not even an `__init__.py`.
The UI catalogue is served from the live registry at request time.
"""

from __future__ import annotations

import importlib
import pkgutil
import re
from collections.abc import Iterator, Sequence

from bt5.core.services import FoldEngine
from bt5.core.spec import Enforcement, Evidence, Spec

_REGISTRY: dict[str, type[Spec]] = {}


class CalibrationMismatchError(RuntimeError):
    """A rule's thresholds were measured on an engine that is not the one running.

    Every structure threshold in BT5 is a kcal/mol number calibrated against one
    engine and one parameter set -- the Boel -39 dual gate, the cap-proximal
    -30/-50/-60 ladder. Applying a ViennaRNA-calibrated number to another
    engine's output is named in the research brief as the single most likely
    correctness bug in the folding feature, and it is silent: the comparison
    succeeds, the rule fires or does not, and nothing looks wrong.
    """


def register(cls: type[Spec]) -> type[Spec]:
    """Class decorator. Validates the contract at import time."""
    rule_id = getattr(cls, "id", None)
    if not rule_id:
        raise ValueError(f"{cls.__name__} must define a class-level `id`")
    if rule_id in _REGISTRY and _REGISTRY[rule_id] is not cls:
        raise ValueError(
            f"duplicate rule id {rule_id!r}: {_REGISTRY[rule_id].__name__} and {cls.__name__}"
        )
    _validate(cls)
    _REGISTRY[rule_id] = cls
    return cls


def _validate(cls: type[Spec]) -> None:
    """The contract every rule must satisfy. Mirrored by a CI test over all rules."""
    name = cls.__name__
    if not getattr(cls, "citations", ()):
        raise ValueError(f"{name}: at least one citation is required")
    if (
        getattr(cls, "enforcement", None) is Enforcement.SOFT
        and not getattr(cls, "weight_provenance", "").strip()
    ):
        raise ValueError(
            f"{name}: SOFT rules must carry a non-empty weight_provenance. "
            f"The default weight vector is what 90% of users actually get."
        )
    is_folklore = getattr(cls, "evidence", None) is Evidence.FOLKLORE
    if is_folklore and getattr(cls, "default_enabled", False):
        raise ValueError(f"{name}: FOLKLORE rules must ship default_enabled=False")
    if not getattr(cls, "last_verified", ""):
        raise ValueError(f"{name}: last_verified (ISO date) is required")
    calibration = getattr(cls, "engine_calibration", None)
    if calibration is not None and not _WELL_FORMED_CALIBRATION.match(calibration):
        raise ValueError(
            f"{name}: engine_calibration must read 'engine:param_set', e.g. "
            f"'viennarna:rna_turner2004', got {calibration!r}. It is compared "
            f"against FoldEnergy.calibration_key, so a typo silently matches "
            f"nothing and the rule would never be refused."
        )


#: `engine:param_set`, the form `FoldEnergy.calibration_key` produces.
_WELL_FORMED_CALIBRATION = re.compile(r"^[A-Za-z0-9_.+-]+:[A-Za-z0-9_.+-]+$")


def check_engine_calibration(
    specs: Sequence[type[Spec]], engine: FoldEngine | None
) -> tuple[type[Spec], ...]:
    """Refuse rules calibrated against a different engine; return the runnable ones.

    Called when the engine is chosen, which is the only moment both facts are
    known: a rule declares its calibration at import time, and the engine that
    will run is injected later. That is why this is not part of `_validate`.

    Three outcomes, deliberately different:
      - a rule declaring no calibration runs against anything, because it does
        not depend on an engine at all;
      - a rule whose calibration matches runs;
      - a rule whose calibration DIFFERS raises, rather than being skipped. A
        skipped rule is a missing constraint that nobody sees; the mismatch is a
        configuration error and the run should stop.

    With no engine at all, calibrated rules are returned as unrunnable rather
    than refused -- absence is a degradation to report (see
    `ObjectiveScore.unavailable`), not a misconfiguration to crash on.
    """
    active = f"{engine.name}:{engine.param_set}" if engine is not None else None
    runnable: list[type[Spec]] = []
    wrong: list[str] = []
    for cls in specs:
        declared = getattr(cls, "engine_calibration", None)
        if declared is None:
            runnable.append(cls)
        elif active is None:
            continue  # unrunnable, but not an error: no engine is available
        elif declared == active:
            runnable.append(cls)
        else:
            wrong.append(f"{getattr(cls, 'id', cls.__name__)} wants {declared}")
    if wrong:
        raise CalibrationMismatchError(
            f"active folding engine is {active!r}, but these rules were calibrated "
            f"against another: {', '.join(sorted(wrong))}. Their thresholds are in "
            f"kcal/mol measured on a different energy model, so the comparison "
            f"would succeed while meaning nothing."
        )
    return tuple(runnable)


def discover(package: str = "bt5.rules.catalog") -> None:
    """Import every module under `package` so decorators run."""
    try:
        mod = importlib.import_module(package)
    except ModuleNotFoundError:
        return
    for info in pkgutil.walk_packages(mod.__path__, prefix=f"{package}."):
        importlib.import_module(info.name)


def all_specs() -> tuple[type[Spec], ...]:
    return tuple(_REGISTRY[k] for k in sorted(_REGISTRY))


def get(rule_id: str) -> type[Spec]:
    return _REGISTRY[rule_id]


def iter_specs() -> Iterator[type[Spec]]:
    yield from all_specs()


def clear() -> None:
    """Test-only."""
    _REGISTRY.clear()
