"""A machine-readable record of what `bt5.core` promises, and how it may change.

docs/PLAN.md freezes `core/` once M1, M2 and M3 land, and gives the amendment
protocol: MINOR on a fast path, MAJOR requires an RFC plus a deprecation shim,
a two-window rule, and a backward-compatibility test over recorded fixtures.
This module is the machinery that makes "frozen" a check rather than an
intention.

Two ideas do all the work.

FIRST, the surface is EXTRACTED, not hand-written. The names come from the
module's own AST -- what the file actually defines at top level -- rather than
from `dir()`, which cannot tell `Strand = Literal[1, -1]` (a contract type) from
`import re` (not). A hand-maintained list of what core promises is a list that
goes stale exactly when it matters.

SECOND, MINOR and MAJOR are decided by asking who breaks. A field added WITH a
default breaks nobody: every existing constructor call still works. The same
field added WITHOUT one breaks every caller at once. That is the whole
classification, applied uniformly:

    MINOR   a new type, a new defaulted field, a new enum member, a field that
            gains a default
    MAJOR   a removed or renamed anything, a changed annotation, a field that
            LOSES its default, a changed signature, a new protocol method

The last two are not hypothetical. PR #22 added a required `duplex()` to the
FoldEngine protocol and PR #26 made `Breach.fixable_by_codon_choice` required;
both were right, and both are exactly the shape of change this gate exists to
make visible and deliberate rather than incidental.
"""

from __future__ import annotations

import ast
import dataclasses
import enum
import importlib
import inspect
from collections.abc import Mapping
from pathlib import Path
from typing import Any

#: Every module in the frozen contract. Listed explicitly: a new core module is
#: itself a contract change, and discovering them by globbing would let one
#: appear without anyone deciding it should.
CORE_MODULES: tuple[str, ...] = (
    "bt5.core.types",
    "bt5.core.context",
    "bt5.core.spec",
    "bt5.core.result",
    "bt5.core.services",
    "bt5.core.registry",
)

MANIFEST_PATH = Path(__file__).resolve().parent / "manifest.json"
RFC_DIR = Path(__file__).resolve().parents[2] / "docs" / "rfcs"


# --- extraction ----------------------------------------------------------


def _defined_names(module: Any) -> tuple[str, ...]:
    """Public top-level names the module's SOURCE defines, in source order.

    From the AST rather than `dir()`, because `dir()` cannot distinguish a type
    alias the contract promises from a module the file happens to import.
    """
    source = Path(inspect.getfile(module)).read_text(encoding="utf-8")
    out: list[str] = []
    for node in ast.parse(source).body:
        if isinstance(node, ast.ClassDef | ast.FunctionDef):
            out.append(node.name)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            out.append(node.target.id)
        elif isinstance(node, ast.Assign):
            out.extend(t.id for t in node.targets if isinstance(t, ast.Name))
    return tuple(n for n in out if not n.startswith("_"))


def _canonical(value: Any) -> str:
    """A repr that does not move between runs.

    Set iteration order is hash-dependent, so `repr(frozenset("ACGT"))` is not
    stable across processes and would make the manifest churn for no reason.
    """
    if isinstance(value, frozenset | set):
        inner = ", ".join(_canonical(v) for v in sorted(value, key=repr))
        return f"{type(value).__name__}({{{inner}}})"
    if isinstance(value, Mapping):
        items = ", ".join(
            f"{_canonical(k)}: {_canonical(v)}"
            for k, v in sorted(value.items(), key=lambda kv: repr(kv[0]))
        )
        return f"{{{items}}}"
    return repr(value)


def _signature(func: Any) -> str:
    try:
        return str(inspect.signature(func))
    except (TypeError, ValueError):  # pragma: no cover - builtins only
        return "(...)"


def _members(cls: type) -> dict[str, str]:
    """Public methods and properties declared on the class itself.

    Inherited members are excluded: they belong to the base's entry, and
    recording them twice would report one change as two.
    """
    out: dict[str, str] = {}
    for name, value in sorted(vars(cls).items()):
        if name.startswith("_"):
            continue
        if isinstance(value, property):
            fget = value.fget
            out[name] = f"property{_signature(fget)}" if fget else "property"
        elif isinstance(value, classmethod | staticmethod):
            kind = "classmethod" if isinstance(value, classmethod) else "staticmethod"
            out[name] = f"{kind}{_signature(value.__func__)}"
        elif callable(value):
            out[name] = _signature(value)
    return out


def _dataclass_entry(cls: type) -> dict[str, Any]:
    fields = []
    for f in dataclasses.fields(cls):
        has_default = f.default is not dataclasses.MISSING
        has_factory = f.default_factory is not dataclasses.MISSING  # type: ignore[misc]
        entry: dict[str, Any] = {
            "name": f.name,
            "annotation": str(f.type),
            "required": not (has_default or has_factory),
        }
        if has_default:
            entry["default"] = _canonical(f.default)
        elif has_factory:
            entry["default"] = "<factory>"
        fields.append(entry)
    params = cls.__dataclass_params__  # type: ignore[attr-defined]
    return {
        "kind": "dataclass",
        "frozen": bool(params.frozen),
        "slots": bool(getattr(cls, "__slots__", None) is not None),
        "fields": fields,
        "members": _members(cls),
    }


def _class_entry(cls: type) -> dict[str, Any]:
    if dataclasses.is_dataclass(cls):
        return _dataclass_entry(cls)
    if isinstance(cls, type) and issubclass(cls, enum.Enum):
        return {
            "kind": "enum",
            "base": cls.__mro__[1].__name__,
            "members": {m.name: _canonical(m.value) for m in cls},
            "methods": _members(cls),
        }
    if issubclass(cls, BaseException):
        return {
            "kind": "exception",
            "base": cls.__mro__[1].__name__,
            "members": _members(cls),
        }
    # A Protocol, or a plain class. ClassVar annotations are part of the
    # promise for Spec, whose whole surface is class-level declarations.
    return {
        "kind": "protocol" if getattr(cls, "_is_protocol", False) else "class",
        "annotations": {k: str(v) for k, v in sorted(getattr(cls, "__annotations__", {}).items())},
        "members": _members(cls),
    }


def extract() -> dict[str, Any]:
    """The current contract surface, as plain JSON-able data."""
    surface: dict[str, Any] = {}
    for module_name in CORE_MODULES:
        module = importlib.import_module(module_name)
        entries: dict[str, Any] = {}
        for name in _defined_names(module):
            value = getattr(module, name, None)
            if value is None and name not in vars(module):
                continue
            if inspect.isclass(value):
                entries[name] = _class_entry(value)
            elif inspect.isfunction(value):
                entries[name] = {"kind": "function", "signature": _signature(value)}
            else:
                entries[name] = {"kind": "constant", "value": _canonical(value)}
        surface[module_name] = dict(sorted(entries.items()))
    return surface


# --- classification ------------------------------------------------------


@dataclasses.dataclass(frozen=True, slots=True)
class Change:
    """One difference between two surfaces, and who it breaks."""

    severity: str  # "MINOR" | "MAJOR"
    path: str
    detail: str

    def __str__(self) -> str:
        return f"{self.severity}  {self.path}: {self.detail}"


def _field_index(entry: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    return {f["name"]: f for f in entry.get("fields", [])}


def _diff_dataclass(path: str, old: Mapping[str, Any], new: Mapping[str, Any]) -> list[Change]:
    out: list[Change] = []
    before, after = _field_index(old), _field_index(new)

    for name in sorted(set(before) - set(after)):
        out.append(Change("MAJOR", f"{path}.{name}", "field removed; every reader of it breaks"))

    for name in sorted(set(after) - set(before)):
        field = after[name]
        if field["required"]:
            out.append(
                Change(
                    "MAJOR",
                    f"{path}.{name}",
                    "field added with no default; every constructor call breaks at once",
                )
            )
        else:
            out.append(Change("MINOR", f"{path}.{name}", "field added with a default"))

    for name in sorted(set(before) & set(after)):
        was, now = before[name], after[name]
        if was["annotation"] != now["annotation"]:
            out.append(
                Change(
                    "MAJOR",
                    f"{path}.{name}",
                    f"annotation changed: {was['annotation']} -> {now['annotation']}",
                )
            )
        if was["required"] and not now["required"]:
            out.append(Change("MINOR", f"{path}.{name}", "field gained a default"))
        elif not was["required"] and now["required"]:
            out.append(
                Change(
                    "MAJOR",
                    f"{path}.{name}",
                    "field lost its default; callers that omitted it now fail",
                )
            )
        elif was.get("default") != now.get("default"):
            out.append(
                Change(
                    "MAJOR",
                    f"{path}.{name}",
                    f"default changed: {was.get('default')} -> {now.get('default')}. "
                    f"A silently different default is a behaviour change nobody reads.",
                )
            )

    for key in ("frozen", "slots"):
        if old.get(key) != new.get(key):
            out.append(
                Change("MAJOR", f"{path}", f"{key} changed: {old.get(key)} -> {new.get(key)}")
            )
    return out


def _diff_mapping(
    path: str, old: Mapping[str, str], new: Mapping[str, str], *, noun: str, added: str
) -> list[Change]:
    out: list[Change] = []
    for name in sorted(set(old) - set(new)):
        out.append(Change("MAJOR", f"{path}.{name}", f"{noun} removed"))
    for name in sorted(set(new) - set(old)):
        out.append(Change(added, f"{path}.{name}", f"{noun} added"))
    for name in sorted(set(old) & set(new)):
        if old[name] != new[name]:
            out.append(
                Change("MAJOR", f"{path}.{name}", f"{noun} changed: {old[name]} -> {new[name]}")
            )
    return out


def _diff_entry(path: str, old: Mapping[str, Any], new: Mapping[str, Any]) -> list[Change]:
    if old.get("kind") != new.get("kind"):
        return [Change("MAJOR", path, f"kind changed: {old.get('kind')} -> {new.get('kind')}")]

    kind = new["kind"]
    out: list[Change] = []
    if kind == "dataclass":
        out += _diff_dataclass(path, old, new)
        out += _diff_mapping(path, old["members"], new["members"], noun="member", added="MINOR")
    elif kind == "enum":
        # A new member is MINOR: nothing that existed stops working. An
        # exhaustive match downstream is the caller's problem to keep total,
        # and refusing new members would freeze every vocabulary in BT5.
        out += _diff_mapping(path, old["members"], new["members"], noun="member", added="MINOR")
        out += _diff_mapping(path, old["methods"], new["methods"], noun="method", added="MINOR")
    elif kind == "protocol":
        # A new protocol method is MAJOR and a new dataclass field is not,
        # because the roles are reversed: BT5 CONSTRUCTS dataclasses and
        # IMPLEMENTS protocols, so an added requirement lands on every
        # implementer -- which for FoldEngine means every lane's fake.
        out += _diff_mapping(
            path, old["annotations"], new["annotations"], noun="declaration", added="MAJOR"
        )
        out += _diff_mapping(path, old["members"], new["members"], noun="method", added="MAJOR")
    elif kind in ("class", "exception"):
        out += _diff_mapping(path, old["members"], new["members"], noun="member", added="MINOR")
        if old.get("base") != new.get("base"):
            out.append(
                Change("MAJOR", path, f"base changed: {old.get('base')} -> {new.get('base')}")
            )
    elif kind == "function":
        if old["signature"] != new["signature"]:
            out.append(
                Change(
                    "MAJOR", path, f"signature changed: {old['signature']} -> {new['signature']}"
                )
            )
    elif kind == "constant" and old["value"] != new["value"]:
        out.append(Change("MAJOR", path, f"value changed: {old['value']} -> {new['value']}"))
    return out


def diff(old: Mapping[str, Any], new: Mapping[str, Any]) -> tuple[Change, ...]:
    """Every difference between two surfaces, each classified MINOR or MAJOR."""
    out: list[Change] = []
    for module in sorted(set(old) | set(new)):
        before = old.get(module, {})
        after = new.get(module, {})
        if module not in new:
            out.append(Change("MAJOR", module, "module removed from the contract"))
            continue
        if module not in old:
            out.append(Change("MINOR", module, "module added to the contract"))
            continue
        for name in sorted(set(before) - set(after)):
            out.append(Change("MAJOR", f"{module}.{name}", "removed from the contract"))
        for name in sorted(set(after) - set(before)):
            out.append(Change("MINOR", f"{module}.{name}", "added to the contract"))
        for name in sorted(set(before) & set(after)):
            out += _diff_entry(f"{module}.{name}", before[name], after[name])
    return tuple(out)


def majors(changes: tuple[Change, ...]) -> tuple[Change, ...]:
    return tuple(c for c in changes if c.severity == "MAJOR")


# --- the amendment decision ----------------------------------------------


@dataclasses.dataclass(frozen=True, slots=True)
class Verdict:
    """Whether this branch's contract changes are properly amended."""

    changes: tuple[Change, ...]
    problems: tuple[str, ...]
    amendment: Mapping[str, Any] | None = None

    @property
    def ok(self) -> bool:
        return not self.problems

    @property
    def major(self) -> tuple[Change, ...]:
        return majors(self.changes)


def _rfc_exists(name: str) -> bool:
    return (RFC_DIR / Path(name).name).exists()


def review(
    baseline: Mapping[str, Any] | None,
    current: Mapping[str, Any],
    live: Mapping[str, Any],
    *,
    rfc_exists: Any = _rfc_exists,
) -> Verdict:
    """Decide whether `current` amends `baseline` legitimately.

    Pure, and takes `rfc_exists` as a predicate, so the decision is testable
    without a baseline branch or real RFC files on disk. A gate whose logic can
    only be exercised by pushing to CI is a gate nobody has tested against the
    case it exists to reject.
    """
    if baseline is None:
        return Verdict(changes=(), problems=())

    changes = diff(baseline["surface"], live)
    if not majors(changes):
        return Verdict(changes=changes, problems=())

    was = baseline.get("contract_version", 0)
    now = current.get("contract_version", 0)
    amendment = next((a for a in current.get("amendments", []) if a.get("version") == now), None)

    problems: list[str] = []
    if now <= was:
        problems.append(
            f"contract_version is still {now}; a MAJOR change must bump it (was {was} on the base)"
        )
    if amendment is None:
        problems.append(f"no entry in `amendments` for contract version {now}")
    else:
        if not amendment.get("summary"):
            problems.append(f"amendment {now} has no summary")
        if not amendment.get("rfc"):
            problems.append("the amendment names no RFC")
        elif not rfc_exists(amendment["rfc"]):
            problems.append(f"the amendment names an RFC that does not exist: {amendment['rfc']}")

    return Verdict(changes=changes, problems=tuple(problems), amendment=amendment)
