"""The contract every registered rule must satisfy.

A rule author gets this entire suite for free and cannot ship a rule missing an
evidence badge, a citation, or (for a scored rule) a written justification for its
default weight. The contract test catches MISSING provenance; it cannot catch
WRONG provenance, which is why `last_verified` plus a standing review rotation
exists alongside it.
"""

from __future__ import annotations

import datetime as dt

import pytest
from bt5.core.registry import all_specs, discover
from bt5.core.spec import Direction, Enforcement, Evidence

discover()
SPECS = all_specs()


def test_at_least_one_rule_is_registered() -> None:
    assert SPECS, "autodiscovery found no rules"


@pytest.mark.parametrize("spec", SPECS, ids=lambda s: s.id)
class TestRuleContract:
    def test_has_citations(self, spec: type) -> None:
        assert spec.citations, f"{spec.id}: at least one citation is required"
        for c in spec.citations:
            assert c.url.startswith("https://"), f"{spec.id}: citation {c.url} must be https"
            assert c.label.strip(), f"{spec.id}: citation needs a label"

    def test_soft_rules_justify_their_default_weight(self, spec: type) -> None:
        """The default weight vector is what 90% of users actually get."""
        if spec.enforcement is Enforcement.SOFT:
            assert spec.weight_provenance.strip(), (
                f"{spec.id}: a SOFT rule must explain, in prose, where its default weight came from"
            )

    def test_folklore_ships_disabled(self, spec: type) -> None:
        if spec.evidence is Evidence.FOLKLORE:
            assert not spec.default_enabled, f"{spec.id}: folklore rules must default off"

    def test_vendor_asserted_rules_are_not_disabled_by_default(self, spec: type) -> None:
        """Vendors enforce these at order time; shipping them off is backwards."""
        if spec.evidence is Evidence.VENDOR_ASSERTED:
            assert spec.default_enabled, f"{spec.id}: vendor-asserted rules should default on"

    def test_last_verified_is_a_real_date(self, spec: type) -> None:
        dt.date.fromisoformat(spec.last_verified)

    def test_band_rules_declare_a_band(self, spec: type) -> None:
        if spec.direction is Direction.BAND:
            assert spec.band is not None, f"{spec.id}: BAND direction requires a band"
            lo, hi = spec.band
            assert lo < hi, f"{spec.id}: band {spec.band} is inverted"

    def test_hard_rules_carry_no_objective_weight(self, spec: type) -> None:
        """A hard constraint is guaranteed by construction, or by repair plus the
        independent validator -- NEVER by a penalty weight.

        A hard rule may still carry a `steering_weight` to nudge the Tier-A DP
        (windowed GC does), but that term is not the objective function and the
        weighted sum never sees it.
        """
        if spec.enforcement.is_hard:
            assert spec.default_weight == 0.0, (
                f"{spec.id}: hard rules must not carry an OBJECTIVE weight; use "
                f"steering_weight if the DP needs nudging"
            )

    def test_only_lattice_rules_skip_steering(self, spec: type) -> None:
        """HARD_LATTICE is unreachable by construction, so steering is pointless."""
        if spec.enforcement is Enforcement.HARD_LATTICE:
            assert spec.steering_weight == 0.0, (
                f"{spec.id}: HARD_LATTICE is already unreachable; steering is a no-op"
            )

    def test_lattice_rules_declare_motifs(self, spec: type) -> None:
        if spec.enforcement is Enforcement.HARD_LATTICE:
            terms = spec().lattice_terms(None)
            assert terms is not None and terms.forbidden, (
                f"{spec.id}: HARD_LATTICE means the automaton makes it unreachable, "
                f"so it must declare forbidden motifs"
            )

    def test_declares_a_brief_reference(self, spec: type) -> None:
        assert spec.brief_ref, f"{spec.id}: must cite its section in docs/research/brief.md"

    def test_param_schema_is_a_json_schema_object(self, spec: type) -> None:
        assert spec.param_schema.get("type") == "object", (
            f"{spec.id}: param_schema drives the UI controls and must be a JSON Schema object"
        )

    def test_is_probeable_per_slot(self, spec: type) -> None:
        """Every spec must be no-arg constructible AND expose the per-slot probe.

        `score/presets.py::_unscored_enforcement` is the guard that stops a
        preset weighting a rule that is HARD in the modality the preset is
        pinned to. It probes by calling `spec()` and then `enforcement_for(slot)`,
        and it has two fallbacks that both resolve to "treat this spec as
        scored", i.e. ALLOW the weight: a `TypeError` from the constructor, and a
        missing `enforcement_for`. Either one silently drops the guard back to
        the class-level `enforcement` ClassVar.

        That ClassVar is only the FLOOR. Reading it alone is exactly the defect
        #72 fixed: `d4_internal_polya` declares SOFT and its `enforcement_for`
        returns HARD_REPAIR on every packaged modality, and two shipped presets
        came to weight it. A future rule that took a required constructor
        argument would sail straight back into that hole.

        CAUGHT HERE RATHER THAN THERE, deliberately (#82). At the preset the
        symptom is a weight that silently appears; at the rule it is a plain fact
        about the rule, with a message naming it. Same argument
        `core/registry.py::check_engine_calibration` makes for raising rather
        than skipping: a skipped rule is a missing constraint nobody sees.

        Unreachable today -- every catalog rule has all-default `__init__`
        parameters and implements the full Protocol -- but that is a coincidence
        of the current catalog, not a property of it, and nothing else gates it.
        """
        try:
            rule = spec()
        except TypeError as exc:  # pragma: no cover - the failure this exists for
            pytest.fail(
                f"{spec.id}: must be constructible with no arguments. "
                f"`_unscored_enforcement` probes rules by calling `spec()`, and on "
                f"TypeError it falls back to the `enforcement` ClassVar -- the floor, "
                f"not what the solver routes on -- which is how a HARD rule silently "
                f"regains an objective weight. Give every parameter a default. ({exc})"
            )
        for member in ("enforcement_for", "gate"):
            assert callable(getattr(rule, member, None)), (
                f"{spec.id}: must expose a callable `{member}`; it is a member of the "
                f"Spec Protocol, and `_unscored_enforcement` treats its absence as "
                f"permission to weight the rule"
            )


def test_rule_ids_are_unique_and_sorted_stably() -> None:
    ids = [s.id for s in SPECS]
    assert len(ids) == len(set(ids))
    assert ids == sorted(ids)


def test_no_rule_declares_a_conflict_with_itself() -> None:
    for s in SPECS:
        assert s.id not in s.conflicts_with, f"{s.id} conflicts with itself"
