"""The default weight vector is the product, so it is tested like one.

docs/PLAN.md: 90% of users never move a slider. Everything here is really one
question asked four ways -- can a number reach the objective function without an
argument attached to it?
"""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

import pytest
from bt5.core.context import ContextSlot, HostId, Modality, SlotRole
from bt5.core.registry import all_specs, discover
from bt5.core.spec import (
    Citation,
    Direction,
    Enforcement,
    Evidence,
    LocalizationPolicy,
    RepairPolicy,
)
from bt5.rules.catalog.d4_internal_polya import InternalPolyA
from bt5.score import PRESETS, Preset, PresetError, WeightEntry, preset_for, resolve
from bt5.score.presets import BACTERIAL, LENTIVIRAL, _slots_admitted_by, get


def lentiviral_slot(role: SlotRole = "producer") -> ContextSlot:
    """A real packaging slot. HEK293 is locked to NCBI table 1."""
    return ContextSlot(role=role, host=HostId.HEK293, modality=Modality.LENTIVIRAL, table_id=1)


def a_preset(modality: Modality, *entries: WeightEntry) -> Preset:
    return Preset(
        id="p",
        title="p",
        rationale="a rationale long enough to be an argument rather than a label",
        modality=modality,
        entries=entries,
    )


def fake_spec(
    spec_id: str,
    brief_ref: str,
    *,
    enforcement: Enforcement = Enforcement.SOFT,
    default_weight: float = 1.0,
    conflicts_with: tuple[str, ...] = (),
    gate_returns: bool = True,
    escalates_to: Enforcement | None = None,
) -> type:
    """A Spec-shaped class, built without touching the process-wide registry.

    Registering real classes here would leak into every other test in the
    session, which has bitten this suite before.

    `gate_returns` and `escalates_to` exist to separate two things that
    COINCIDE in every real catalog rule, and therefore cannot be told apart by
    a test built on one. `escalates_to` sets what `enforcement_for` returns,
    independently of the class-level `enforcement` floor; `gate_returns` sets
    whether the rule applies at all. Without them a double's `enforcement_for`
    just echoes its ClassVar, so a test cannot distinguish "the guard read the
    ClassVar" from "the guard asked per slot" -- which is the entire question
    this file exists to settle.
    """

    class _Fake:
        id: ClassVar[str] = spec_id
        version: ClassVar[str] = "1.0.0"
        title: ClassVar[str] = spec_id
        enforcement: ClassVar[Enforcement] = Enforcement.SOFT
        evidence: ClassVar[Evidence] = Evidence.EVIDENCE_BACKED
        direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
        unit: ClassVar[str] = "au"
        citations: ClassVar[tuple[Citation, ...]] = (Citation("x", "https://example.org"),)
        last_verified: ClassVar[str] = "2026-08-28"
        weight_provenance: ClassVar[str] = "test"
        default_enabled: ClassVar[bool] = True
        default_weight: ClassVar[float] = 1.0
        steering_weight: ClassVar[float] = 0.0
        band: ClassVar[tuple[float, float] | None] = None
        localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.WHOLE_SCOPE
        repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
        cost_class: ClassVar[str] = "cheap"
        conflicts_with: ClassVar[tuple[str, ...]] = ()
        param_schema: ClassVar[dict[str, object]] = {}
        # Assigned below: a class body cannot read an enclosing local whose name
        # it also binds.
        brief_ref: ClassVar[str] = ""
        engine_calibration: ClassVar[str | None] = None

        # `Spec` requires both, and `resolve()` now asks them per slot. A double
        # missing them would let the resolver quietly fall back to the ClassVar
        # -- the very read this test file exists to stop it making.
        def gate(self, slot: ContextSlot) -> bool:
            return gate_returns

        def enforcement_for(self, slot: ContextSlot) -> Enforcement:
            # Deliberately NOT `type(self).enforcement` unless asked: a double
            # that echoes its own ClassVar cannot prove the guard consulted
            # the slot rather than the class.
            return escalates_to if escalates_to is not None else type(self).enforcement

    _Fake.brief_ref = brief_ref
    _Fake.enforcement = enforcement
    _Fake.default_weight = default_weight
    _Fake.conflicts_with = conflicts_with
    _Fake.__name__ = spec_id
    return _Fake


class TestPresetContract:
    def test_every_shipped_preset_carries_a_rationale(self) -> None:
        for preset in PRESETS:
            assert preset.rationale.strip(), f"{preset.id} has no rationale"
            assert len(preset.rationale) > 80, (
                f"{preset.id}: a rationale short enough to be a label is not an "
                f"argument, and for most users this vector is the only claim "
                f"BT5 makes about what matters"
            )

    def test_a_preset_without_a_rationale_cannot_be_built(self) -> None:
        with pytest.raises(PresetError, match="rationale"):
            Preset(id="x", title="x", rationale="  ", modality=Modality.AAV)

    def test_shipped_preset_ids_are_unique(self) -> None:
        ids = [p.id for p in PRESETS]
        assert len(set(ids)) == len(ids)

    def test_a_ref_cannot_be_weighted_twice_in_one_preset(self) -> None:
        with pytest.raises(PresetError, match="more than once"):
            Preset(
                id="x",
                title="x",
                rationale="a rationale long enough to be an argument rather than a label",
                modality=Modality.AAV,
                entries=(WeightEntry("2.C1", 0.1, "a"), WeightEntry("2.C1", 0.2, "b")),
            )

    def test_a_negative_weight_is_refused(self) -> None:
        """Inverting an objective's direction is a different rule, not a weight."""
        with pytest.raises(PresetError, match="must be >= 0"):
            WeightEntry("2.C1", -1.0, "note")

    def test_lookup_by_id_and_by_modality(self) -> None:
        assert get("lentiviral_hek293") is LENTIVIRAL
        assert preset_for(Modality.BACTERIAL_EXPRESSION) is BACTERIAL
        with pytest.raises(PresetError, match="no preset"):
            get("nope")

    def test_an_uncurated_modality_returns_none_rather_than_the_nearest_thing(self) -> None:
        """Handing a lentiviral vector to an IVT mRNA design because it was on
        the shelf is worse than saying there is no curated default."""
        assert preset_for(Modality.IVT_MRNA) is None


class TestResolve:
    def test_binds_refs_to_the_spec_ids_present_in_this_build(self) -> None:
        preset = Preset(
            id="p",
            title="p",
            rationale="a rationale long enough to be an argument rather than a label",
            modality=Modality.AAV,
            entries=(WeightEntry("2.C1", 1.0),),
        )
        resolved = resolve(preset, [fake_spec("c1_cai", "2.C1")])
        assert resolved.weights == {"c1_cai": 1.0}
        assert resolved.unimplemented == ()
        assert resolved.degradations == ()

    def test_an_unimplemented_objective_is_reported_not_swallowed(self) -> None:
        """A user who believes an objective is being optimised, and whose rule
        for it does not exist, is reading a ranking that does not mean what they
        think it means."""
        preset = Preset(
            id="p",
            title="p",
            rationale="a rationale long enough to be an argument rather than a label",
            modality=Modality.AAV,
            entries=(WeightEntry("2.B1", 1.0),),
        )
        resolved = resolve(preset, [fake_spec("c1_cai", "2.C1")])
        assert resolved.weights == {}
        assert resolved.unimplemented == ("2.B1",)
        assert resolved.degradations, "an absent objective must reach Provenance.degradations"
        assert "2.B1" in resolved.degradations[0]

    @pytest.mark.parametrize(
        "enforcement",
        [Enforcement.HARD_LATTICE, Enforcement.HARD_REPAIR, Enforcement.HARD_CHECK],
    )
    def test_weighting_a_hard_rule_is_refused(self, enforcement: Enforcement) -> None:
        """The whole point of the Enforcement enum: a hard constraint is never
        enforced by a penalty weight, and giving one a weight both adds a term
        for something already guaranteed and implies it is a trade-off."""
        preset = Preset(
            id="p",
            title="p",
            rationale="a rationale long enough to be an argument rather than a label",
            modality=Modality.AAV,
            entries=(WeightEntry("2.D1", 0.5, "note"),),
        )
        spec = fake_spec("d1_sites", "2.D1", enforcement=enforcement, default_weight=0.0)
        with pytest.raises(PresetError, match="never by a penalty weight"):
            resolve(preset, [spec])

    def test_zero_weight_on_a_hard_rule_is_allowed(self) -> None:
        """Zero is not a penalty. A preset may name a hard rule to say out loud
        that it contributes nothing to the sum."""
        preset = Preset(
            id="p",
            title="p",
            rationale="a rationale long enough to be an argument rather than a label",
            modality=Modality.AAV,
            entries=(WeightEntry("2.D1", 0.0),),
        )
        spec = fake_spec(
            "d1_sites", "2.D1", enforcement=Enforcement.HARD_LATTICE, default_weight=0.0
        )
        assert resolve(preset, [spec]).weights == {"d1_sites": 0.0}

    def test_overriding_a_default_weight_without_saying_why_is_refused(self) -> None:
        preset = Preset(
            id="p",
            title="p",
            rationale="a rationale long enough to be an argument rather than a label",
            modality=Modality.AAV,
            entries=(WeightEntry("2.C1", 0.25),),  # no note
        )
        with pytest.raises(PresetError, match="no note saying"):
            resolve(preset, [fake_spec("c1_cai", "2.C1", default_weight=1.0)])

    def test_the_same_override_with_a_note_is_accepted(self) -> None:
        preset = Preset(
            id="p",
            title="p",
            rationale="a rationale long enough to be an argument rather than a label",
            modality=Modality.AAV,
            entries=(WeightEntry("2.C1", 0.25, "the evidence for CAI does not support more"),),
        )
        assert resolve(preset, [fake_spec("c1_cai", "2.C1", default_weight=1.0)]).weights == {
            "c1_cai": 0.25
        }

    def test_two_rules_claiming_one_brief_ref_is_an_error(self) -> None:
        """Presets key on brief_ref, so a duplicate makes every preset ambiguous."""
        preset = Preset(
            id="p",
            title="p",
            rationale="a rationale long enough to be an argument rather than a label",
            modality=Modality.AAV,
            entries=(WeightEntry("2.C1", 1.0),),
        )
        with pytest.raises(PresetError, match="two rules claim"):
            resolve(preset, [fake_spec("c1_a", "2.C1"), fake_spec("c1_b", "2.C1")])


class TestShippedWeights:
    def test_every_entry_that_departs_from_a_default_carries_a_note(self) -> None:
        """Enforced at resolve time against real rules; asserted here against the
        preset text itself, so the argument exists before the rule does."""
        for preset in PRESETS:
            for entry in preset.entries:
                assert entry.note.strip(), (
                    f"{preset.id}/{entry.brief_ref}: every shipped weight states its case"
                )

    def test_the_packaging_presets_weight_repeats_above_the_expression_proxies(self) -> None:
        """Q4a: recA- strains cover the LONG repeats only. The 15-100 bp repeats
        codon choice controls sit in the RecA-INDEPENDENT regime the strain does
        not suppress, so they are BT5's job and are weighted like it."""
        for preset in (get("lentiviral_hek293"), get("aav_hek293")):
            by_ref = preset.by_ref
            repeat = by_ref["2.F2"].weight
            for ref in ("2.C1", "2.C3"):
                assert repeat > by_ref[ref].weight, (
                    f"{preset.id}: repeats must outweigh {ref}; the repeat rules are "
                    f"mechanically real and the codon-composition proxies are the "
                    f"ones nine benchmarked optimizers could not beat native with"
                )

    def test_the_strain_protocol_says_what_it_does_not_cover(self) -> None:
        """The report must not tell a user that recA- covers the short repeats."""
        for preset in (get("lentiviral_hek293"), get("aav_hek293")):
            assert preset.strain_protocol, f"{preset.id} recommends a strain"
            text = " ".join(preset.strain_protocol).lower()
            assert "long repeats" in text or "does not cover" in text, (
                f"{preset.id}: a strain recommendation that does not state its "
                f"limit reads as blanket protection it does not provide"
            )

    def test_the_bacterial_preset_gives_the_five_prime_window_top_weight(self) -> None:
        """Kudla 2009: the -4..+37 window explains 44-59% of expression variance;
        CAI's r = 0.14 was not significant. Nothing else in BT5 has that support."""
        by_ref = BACTERIAL.by_ref
        top = max(by_ref.values(), key=lambda e: e.weight)
        assert top.brief_ref == "2.B1"
        assert by_ref["2.B1"].weight > by_ref["2.C1"].weight

    def test_every_weighted_ref_is_a_real_row_in_the_brief(self) -> None:
        """The guard for the bug this test was written after.

        The packaging presets originally weighted a `2.G1` that does not exist:
        section 2.G of the brief is prose, not a numbered table, so the ref was
        invented. It failed silently in the worst way -- `resolve` reported it
        as an unimplemented objective, which reads exactly like a rule nobody
        has written yet rather than a row nobody will ever write.

        Rows appear either as a table row `| D4 |` or as bolded prose `**D4 `,
        because the brief uses both forms.
        """
        brief = (Path(__file__).resolve().parents[4] / "docs" / "research" / "brief.md").read_text()
        missing: list[str] = []
        for preset in PRESETS:
            for entry in preset.entries:
                row = entry.brief_ref.split(".")[-1]
                if f"| {row} |" not in brief and f"**{row} " not in brief:
                    missing.append(f"{preset.id} -> {entry.brief_ref}")
        assert not missing, f"presets weight brief rows that do not exist: {missing}"

    def test_the_packaging_presets_do_not_weight_internal_polya(self) -> None:
        """This asserted the opposite until issue #72, and was wrong.

        The measured 8-9x functional titer loss is exactly why 2.D4 must NOT be
        weighted here. `enforcement_for` makes d4 HARD_REPAIR in both packaged
        modalities, so it is removed by Tier-B repair and proven by the
        independent validator, which refuses to emit. A weight on top of that
        adds a term to the sum for something already guaranteed and tells the
        user the guarantee is a slider.
        """
        for preset in (get("lentiviral_hek293"), get("aav_hek293")):
            assert "2.D4" not in preset.by_ref

    def test_no_shipped_preset_weights_a_hard_rule_in_this_build(self) -> None:
        """Runs resolve() against the live registry, which is what the pipeline
        will do. It is allowed to report unimplemented objectives; it is not
        allowed to put a hard rule in the weighted sum."""
        for preset in PRESETS:
            resolve(preset)  # raises if any bound rule is not SOFT


class TestAHardRuleNeverCarriesWeightInTheSlotItIsHardIn:
    """The invariant issue #72 was about, pinned on its own.

    `Enforcement` is per SLOT, not per class. `resolve()` used to read the
    class-level `enforcement` ClassVar with no slot in scope, so it asked
    whether a rule is hard EVERYWHERE instead of whether it is hard HERE --
    and passed a preset weighting a rule that is hard in the very modality the
    preset is pinned to.

    d4_internal_polya is the case, and it is not academic: both packaging
    presets shipped `WeightEntry("2.D4", 1.0)` past the old guard.
    """

    def test_d4_disagrees_with_its_own_classvar_in_a_lentiviral_slot(self) -> None:
        """The disagreement the resolver has to see. If this ever stops being
        true the rest of this class proves nothing, so assert it first."""
        rule = InternalPolyA()
        slot = lentiviral_slot()
        assert InternalPolyA.enforcement is Enforcement.SOFT, "the ClassVar is the FLOOR"
        assert rule.gate(slot), "d4 applies in a lentiviral slot"
        assert rule.enforcement_for(slot) is Enforcement.HARD_REPAIR

    def test_resolve_refuses_to_weight_d4_in_a_lentiviral_preset(self) -> None:
        """The regression. This preset resolved cleanly before #72."""
        preset = a_preset(Modality.LENTIVIRAL, WeightEntry("2.D4", 1.0, "a note"))
        with pytest.raises(PresetError, match="never by a penalty weight"):
            resolve(preset, [InternalPolyA])

    def test_the_refusal_names_the_per_slot_enforcement_not_the_classvar(self) -> None:
        """Reporting `soft` here would send a reader looking at the ClassVar,
        which is not the thing that made the weight illegal."""
        preset = a_preset(Modality.AAV, WeightEntry("2.D4", 1.0, "a note"))
        with pytest.raises(PresetError) as excinfo:
            resolve(preset, [InternalPolyA])
        assert "hard_repair" in str(excinfo.value)
        assert "aav" in str(excinfo.value)

    def test_the_same_rule_is_still_weightable_where_it_is_genuinely_soft(self) -> None:
        """Not a blanket ban on d4. In a plasmid modality the same hexamer costs
        a little expression and nothing else, `enforcement_for` returns SOFT,
        and refusing there would be the mirror image of the bug -- a guard
        that is also not asking about the slot in front of it.
        """
        preset = a_preset(Modality.PLASMID_TRANSIENT, WeightEntry("2.D4", 0.7))
        assert resolve(preset, [InternalPolyA]).weights == {"d4_internal_polya": 0.7}

    def test_the_guard_scopes_to_the_presets_own_modality(self) -> None:
        """d4 is HARD_REPAIR under lentiviral yet weightable in a BACTERIAL
        preset, because the guard asks about the modality in front of it.

        Named for what it actually proves. It was called
        `test_a_rule_gated_off_in_every_slot_is_not_treated_as_hard`, which it
        did NOT establish: d4 both gates off AND is soft under bacterial, so
        the assertion below holds whether or not the guard consults `gate` at
        all. The gate skip is pinned separately, on a double where the two
        come apart.
        """
        rule = InternalPolyA()
        bacterial = ContextSlot(
            role="propagation",
            host=HostId.E_COLI_K12,
            modality=Modality.BACTERIAL_EXPRESSION,
            table_id=11,
        )
        assert not rule.gate(bacterial)
        preset = a_preset(Modality.BACTERIAL_EXPRESSION, WeightEntry("2.D4", 0.7))
        assert resolve(preset, [InternalPolyA]).weights == {"d4_internal_polya": 0.7}

    def test_a_rule_that_gates_off_keeps_its_weight_even_where_it_would_be_hard(
        self,
    ) -> None:
        """Pins the gate skip itself -- `if not gate(slot): continue`.

        No real rule can pin it: every catalog rule that gates off in a slot is
        also SOFT there, so the branch is invisible to a test built on one.
        This double gates OFF while `enforcement_for` says HARD_REPAIR, so the
        two conditions come apart and deleting the skip changes the answer.

        The behaviour: a rule that never fires in this modality contributes
        nothing to its sum, so there is no penalty weight to refuse. Refusing
        anyway would block a legitimate weight.
        """
        spec = fake_spec(
            "gates_off",
            "2.Z1",
            enforcement=Enforcement.SOFT,
            default_weight=0.5,
            gate_returns=False,
            escalates_to=Enforcement.HARD_REPAIR,
        )
        preset = a_preset(Modality.LENTIVIRAL, WeightEntry("2.Z1", 0.5))
        assert resolve(preset, [spec]).weights == {"gates_off": 0.5}

    def test_a_soft_floor_that_escalates_in_an_active_slot_is_refused(self) -> None:
        """The same shape as d4, on a double: SOFT ClassVar, gate ON,
        HARD_REPAIR per slot. Independent of whether d4 keeps its current
        HARD_MODALITIES, so the invariant survives a change to that rule."""
        spec = fake_spec(
            "escalates",
            "2.Z2",
            enforcement=Enforcement.SOFT,
            default_weight=0.5,
            escalates_to=Enforcement.HARD_REPAIR,
        )
        preset = a_preset(Modality.LENTIVIRAL, WeightEntry("2.Z2", 0.5))
        with pytest.raises(PresetError, match="never by a penalty weight"):
            resolve(preset, [spec])

    def test_a_hard_classvar_is_refused_even_when_enforcement_for_says_soft(
        self,
    ) -> None:
        """Pins the ClassVar floor check, which nothing else did.

        The ClassVar is the FLOOR: `enforcement_for` escalates from it and must
        never de-escalate below it. A rule declaring HARD_LATTICE keeps its
        guarantee no matter what it answers per slot -- so dropping the floor
        check must not pass silently.
        """
        spec = fake_spec(
            "hard_floor",
            "2.Z3",
            enforcement=Enforcement.HARD_LATTICE,
            default_weight=0.5,
            escalates_to=Enforcement.SOFT,
        )
        preset = a_preset(Modality.LENTIVIRAL, WeightEntry("2.Z3", 0.5))
        with pytest.raises(PresetError, match="never by a penalty weight"):
            resolve(preset, [spec])

    def test_every_shipped_weight_is_scored_in_every_slot_its_modality_admits(
        self,
    ) -> None:
        """The invariant over the live registry, computed WITHOUT going through
        `resolve()` -- a guard cannot be its own witness. This is the assertion
        that would have caught #72 before it shipped.
        """
        discover()
        by_ref = {spec.brief_ref: spec for spec in all_specs()}
        for preset in PRESETS:
            for entry in preset.entries:
                spec = by_ref.get(entry.brief_ref)
                if spec is None or not entry.weight:
                    continue
                rule = spec()
                for slot in _slots_admitted_by(preset.modality):
                    if not rule.gate(slot):
                        continue
                    level = rule.enforcement_for(slot)
                    assert level.is_scored, (
                        f"{preset.id} weights {spec.id} ({entry.brief_ref}) at "
                        f"{entry.weight}, but it is {level.value} in a "
                        f"{slot.role}/{slot.host}/{slot.modality} slot. A hard "
                        f"constraint is enforced by repair plus the independent "
                        f"validator, never by a penalty weight."
                    )
