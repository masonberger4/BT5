"""The default weight vector is the product, so it is tested like one.

docs/PLAN.md: 90% of users never move a slider. Everything here is really one
question asked four ways -- can a number reach the objective function without an
argument attached to it?
"""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

import pytest
from bt5.core.context import Modality
from bt5.core.spec import (
    Citation,
    Direction,
    Enforcement,
    Evidence,
    LocalizationPolicy,
    RepairPolicy,
)
from bt5.score import PRESETS, Preset, PresetError, WeightEntry, preset_for, resolve
from bt5.score.presets import BACTERIAL, LENTIVIRAL, get


def fake_spec(
    spec_id: str,
    brief_ref: str,
    *,
    enforcement: Enforcement = Enforcement.SOFT,
    default_weight: float = 1.0,
    conflicts_with: tuple[str, ...] = (),
) -> type:
    """A Spec-shaped class, built without touching the process-wide registry.

    Registering real classes here would leak into every other test in the
    session, which has bitten this suite before.
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

    def test_the_packaging_presets_weight_internal_polya(self) -> None:
        """The measured 8-9x functional titer loss, which is the whole reason
        these two presets differ from the plasmid case."""
        for preset in (get("lentiviral_hek293"), get("aav_hek293")):
            assert "2.D4" in preset.by_ref

    def test_no_shipped_preset_weights_a_hard_rule_in_this_build(self) -> None:
        """Runs resolve() against the live registry, which is what the pipeline
        will do. It is allowed to report unimplemented objectives; it is not
        allowed to put a hard rule in the weighted sum."""
        for preset in PRESETS:
            resolve(preset)  # raises if any bound rule is not SOFT
