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
from collections.abc import Iterator

from bt5.core.spec import Enforcement, Evidence, Spec

_REGISTRY: dict[str, type[Spec]] = {}


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
