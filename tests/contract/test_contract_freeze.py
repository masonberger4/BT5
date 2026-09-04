"""`core/` is frozen: the live contract must match what was recorded.

This is the local half of the gate, and it answers one question -- did anyone
change `bt5.core` without recording it? A drifted manifest means the freeze is
decorative, because every later comparison is against a baseline that already
moved.

The MINOR-versus-MAJOR half runs in CI, against main's manifest
(.github/scripts/check-contract-freeze.sh). It cannot run here: once you
regenerate, the local manifest and the live code agree by construction, and the
severity of what you just did is only visible against the version you started
from.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import surface  # noqa: E402

REGENERATE = "python tests/contract/regenerate.py"


@pytest.fixture(scope="module")
def manifest() -> dict:
    assert surface.MANIFEST_PATH.exists(), f"no recorded contract; run {REGENERATE}"
    return json.loads(surface.MANIFEST_PATH.read_text(encoding="utf-8"))


def test_the_live_contract_matches_the_recorded_one(manifest: dict) -> None:
    changes = surface.diff(manifest["surface"], surface.extract())
    if not changes:
        return
    major = surface.majors(changes)
    listing = "\n".join(f"  {c}" for c in changes)
    pytest.fail(
        f"bt5.core has drifted from the recorded contract "
        f"({len(changes)} change(s), {len(major)} MAJOR):\n{listing}\n\n"
        f"Run `{REGENERATE}` to record it. If any change above is MAJOR you also "
        f"owe an amendment: bump contract_version, add an entry to `amendments` "
        f"naming an RFC under docs/rfcs/, and ship the deprecation shim."
    )


def test_the_manifest_covers_every_core_module(manifest: dict) -> None:
    """A new module inside core/ is a contract change like any other, and the
    surface list is explicit so one cannot appear without a decision."""
    core = Path(__file__).resolve().parents[2] / "packages/engine/src/bt5/core"
    on_disk = {f"bt5.core.{p.stem}" for p in core.glob("*.py") if not p.stem.startswith("_")}
    assert on_disk == set(surface.CORE_MODULES), (
        f"core modules on disk do not match the frozen list. "
        f"Unrecorded: {sorted(on_disk - set(surface.CORE_MODULES))}; "
        f"recorded but missing: {sorted(set(surface.CORE_MODULES) - on_disk)}"
    )
    assert set(manifest["surface"]) == set(surface.CORE_MODULES)


def test_every_amendment_names_an_rfc_that_exists(manifest: dict) -> None:
    """An amendment log entry pointing at a missing RFC is worse than no log:
    it reads as a documented decision and is not one."""
    for entry in manifest["amendments"]:
        for key in ("version", "rfc", "summary"):
            assert entry.get(key), f"amendment entry missing {key!r}: {entry}"
        rfc = surface.RFC_DIR / Path(entry["rfc"]).name
        assert rfc.exists(), f"amendment {entry['version']} names a missing RFC: {entry['rfc']}"


def test_the_contract_version_matches_the_amendment_log(manifest: dict) -> None:
    """Version 1 is the freeze itself; each MAJOR amendment adds exactly one."""
    assert manifest["contract_version"] == len(manifest["amendments"]) + 1, (
        "contract_version must be one more than the number of recorded "
        "amendments; a bump with no entry is an undocumented breaking change"
    )


class TestClassification:
    """The rules the gate applies, stated as tests so they are reviewable.

    Each case is a two-surface diff, so the classifier is exercised on the
    shapes it will actually meet rather than on the current contract, which by
    definition has no differences to classify.
    """

    def surface_with(self, entry: dict) -> dict:
        return {"bt5.core.types": {"Thing": entry}}

    def dataclass(self, *fields: dict) -> dict:
        return {
            "kind": "dataclass",
            "frozen": True,
            "slots": True,
            "fields": list(fields),
            "members": {},
        }

    def field(self, name: str, annotation: str = "int", **kw: object) -> dict:
        return {"name": name, "annotation": annotation, "required": True, **kw}

    def test_a_defaulted_field_is_minor(self) -> None:
        """It breaks nobody: every existing constructor call still works."""
        old = self.surface_with(self.dataclass(self.field("a")))
        new = self.surface_with(
            self.dataclass(self.field("a"), self.field("b", required=False, default="0"))
        )
        changes = surface.diff(old, new)
        assert [c.severity for c in changes] == ["MINOR"]

    def test_the_same_field_without_a_default_is_major(self) -> None:
        """Every constructor call breaks at once."""
        old = self.surface_with(self.dataclass(self.field("a")))
        new = self.surface_with(self.dataclass(self.field("a"), self.field("b")))
        assert surface.majors(surface.diff(old, new))

    def test_a_field_losing_its_default_is_major(self) -> None:
        """The real case: Breach.fixable_by_codon_choice in PR #26."""
        old = self.surface_with(self.dataclass(self.field("a", required=False, default="True")))
        new = self.surface_with(self.dataclass(self.field("a")))
        major = surface.majors(surface.diff(old, new))
        assert major
        assert "lost its default" in major[0].detail

    def test_a_field_gaining_a_default_is_minor(self) -> None:
        old = self.surface_with(self.dataclass(self.field("a")))
        new = self.surface_with(self.dataclass(self.field("a", required=False, default="0")))
        assert not surface.majors(surface.diff(old, new))

    def test_a_changed_annotation_is_major(self) -> None:
        old = self.surface_with(self.dataclass(self.field("a", "int")))
        new = self.surface_with(self.dataclass(self.field("a", "str")))
        assert surface.majors(surface.diff(old, new))

    def test_a_changed_default_is_major(self) -> None:
        """A silently different default is a behaviour change nobody reads."""
        old = self.surface_with(self.dataclass(self.field("a", required=False, default="1")))
        new = self.surface_with(self.dataclass(self.field("a", required=False, default="2")))
        assert surface.majors(surface.diff(old, new))

    def test_a_removed_field_is_major(self) -> None:
        old = self.surface_with(self.dataclass(self.field("a"), self.field("b")))
        new = self.surface_with(self.dataclass(self.field("a")))
        assert surface.majors(surface.diff(old, new))

    def test_a_new_protocol_method_is_major(self) -> None:
        """Reversed from a dataclass field, because the roles are reversed: BT5
        IMPLEMENTS protocols, so a new requirement lands on every implementer.
        The real case is FoldEngine.duplex in PR #22."""
        old = {"m": {"P": {"kind": "protocol", "annotations": {}, "members": {"a": "(self)"}}}}
        new = {
            "m": {
                "P": {
                    "kind": "protocol",
                    "annotations": {},
                    "members": {"a": "(self)", "b": "(self)"},
                }
            }
        }
        major = surface.majors(surface.diff(old, new))
        assert major
        assert "method added" in major[0].detail

    def test_a_new_enum_member_is_minor(self) -> None:
        """Nothing that existed stops working; refusing this would freeze every
        vocabulary in BT5."""
        old = {
            "m": {"E": {"kind": "enum", "base": "StrEnum", "members": {"A": "'a'"}, "methods": {}}}
        }
        new = {
            "m": {
                "E": {
                    "kind": "enum",
                    "base": "StrEnum",
                    "members": {"A": "'a'", "B": "'b'"},
                    "methods": {},
                }
            }
        }
        assert not surface.majors(surface.diff(old, new))

    def test_a_changed_enum_value_is_major(self) -> None:
        """The value is what gets serialised into a GenBank note and a report."""
        old = {
            "m": {"E": {"kind": "enum", "base": "StrEnum", "members": {"A": "'a'"}, "methods": {}}}
        }
        new = {
            "m": {"E": {"kind": "enum", "base": "StrEnum", "members": {"A": "'z'"}, "methods": {}}}
        }
        assert surface.majors(surface.diff(old, new))

    def test_a_changed_function_signature_is_major(self) -> None:
        old = {"m": {"f": {"kind": "function", "signature": "(a: int) -> int"}}}
        new = {"m": {"f": {"kind": "function", "signature": "(a: int, b: int) -> int"}}}
        assert surface.majors(surface.diff(old, new))

    def test_a_changed_constant_is_major(self) -> None:
        """A threshold that moves silently is the quietest breaking change there is."""
        old = {"m": {"K": {"kind": "constant", "value": "30"}}}
        new = {"m": {"K": {"kind": "constant", "value": "40"}}}
        assert surface.majors(surface.diff(old, new))

    def test_a_removed_name_is_major_and_a_new_one_is_minor(self) -> None:
        old = {"m": {"A": {"kind": "constant", "value": "1"}}}
        new = {"m": {"B": {"kind": "constant", "value": "1"}}}
        changes = surface.diff(old, new)
        assert {c.severity for c in changes} == {"MAJOR", "MINOR"}

    def test_an_identical_surface_produces_no_changes(self) -> None:
        live = surface.extract()
        assert surface.diff(live, live) == ()


class TestAmendmentReview:
    """The decision CI makes, exercised against the cases it exists to reject.

    A gate whose logic can only be run by pushing to CI is a gate nobody has
    tested. `surface.review` is pure and takes an `rfc_exists` predicate for
    exactly this reason.
    """

    def surfaces(self) -> tuple[dict, dict]:
        """A baseline and a live surface differing by one MAJOR change: a field
        that lost its default. The shape of the real PR #26 change."""
        base = {
            "m": {
                "T": {
                    "kind": "dataclass",
                    "frozen": True,
                    "slots": True,
                    "fields": [
                        {"name": "a", "annotation": "bool", "required": False, "default": "True"}
                    ],
                    "members": {},
                }
            }
        }
        live = {
            "m": {
                "T": {
                    "kind": "dataclass",
                    "frozen": True,
                    "slots": True,
                    "fields": [{"name": "a", "annotation": "bool", "required": True}],
                    "members": {},
                }
            }
        }
        return base, live

    def manifest(self, version: int, amendments: list[dict], live: dict) -> dict:
        return {"contract_version": version, "amendments": amendments, "surface": live}

    def test_no_baseline_is_the_freeze_itself(self) -> None:
        _, live = self.surfaces()
        assert surface.review(None, self.manifest(1, [], live), live).ok

    def test_a_minor_only_change_needs_no_amendment(self) -> None:
        base, _ = self.surfaces()
        live = json.loads(json.dumps(base))
        live["m"]["T"]["fields"].append(
            {"name": "b", "annotation": "int", "required": False, "default": "0"}
        )
        verdict = surface.review(
            {"surface": base, "contract_version": 1}, self.manifest(1, [], live), live
        )
        assert verdict.changes
        assert not verdict.major
        assert verdict.ok

    def test_a_major_change_with_no_amendment_is_rejected(self) -> None:
        base, live = self.surfaces()
        verdict = surface.review(
            {"surface": base, "contract_version": 1}, self.manifest(1, [], live), live
        )
        assert verdict.major
        assert not verdict.ok
        assert any("contract_version is still 1" in p for p in verdict.problems)

    def test_a_major_change_with_a_complete_amendment_passes(self) -> None:
        base, live = self.surfaces()
        amendment = {"version": 2, "rfc": "docs/rfcs/0001-example.md", "summary": "why"}
        verdict = surface.review(
            {"surface": base, "contract_version": 1},
            self.manifest(2, [amendment], live),
            live,
            rfc_exists=lambda _: True,
        )
        assert verdict.major
        assert verdict.ok, verdict.problems

    def test_an_amendment_naming_a_missing_rfc_is_rejected(self) -> None:
        """A log entry pointing at a missing RFC reads as a documented decision
        and is not one."""
        base, live = self.surfaces()
        amendment = {"version": 2, "rfc": "docs/rfcs/0001-nope.md", "summary": "why"}
        verdict = surface.review(
            {"surface": base, "contract_version": 1},
            self.manifest(2, [amendment], live),
            live,
            rfc_exists=lambda _: False,
        )
        assert not verdict.ok
        assert any("does not exist" in p for p in verdict.problems)

    def test_a_bumped_version_with_no_matching_entry_is_rejected(self) -> None:
        base, live = self.surfaces()
        verdict = surface.review(
            {"surface": base, "contract_version": 1},
            self.manifest(2, [], live),
            live,
            rfc_exists=lambda _: True,
        )
        assert not verdict.ok
        assert any("no entry in `amendments`" in p for p in verdict.problems)

    def test_an_amendment_without_a_summary_is_rejected(self) -> None:
        """The RFC is the argument; the summary is what a reader sees first."""
        base, live = self.surfaces()
        amendment = {"version": 2, "rfc": "docs/rfcs/0001-example.md", "summary": ""}
        verdict = surface.review(
            {"surface": base, "contract_version": 1},
            self.manifest(2, [amendment], live),
            live,
            rfc_exists=lambda _: True,
        )
        assert not verdict.ok
        assert any("no summary" in p for p in verdict.problems)

    def test_the_real_pr22_and_pr26_changes_are_caught(self) -> None:
        """Both landed before the freeze and both were right. Against a
        pre-#22 baseline they are exactly what this gate must flag: a new
        protocol method every implementer must supply, and a field that lost
        its default."""
        live = surface.extract()
        base = json.loads(json.dumps(live))
        for field in base["bt5.core.spec"]["Breach"]["fields"]:
            if field["name"] == "fixable_by_codon_choice":
                field["required"] = False
                field["default"] = "True"
        del base["bt5.core.services"]["FoldEngine"]["members"]["duplex"]

        major = surface.majors(surface.diff(base, live))
        details = {c.path: c.detail for c in major}
        assert "bt5.core.spec.Breach.fixable_by_codon_choice" in details
        assert "lost its default" in details["bt5.core.spec.Breach.fixable_by_codon_choice"]
        assert "bt5.core.services.FoldEngine.duplex" in details
