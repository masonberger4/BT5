"""The wiring that makes the rule catalog run.

Before this module's subject existed, `grep -rn "\\.evaluate(" src/` returned
nothing. Most of what is tested here is not "does the code work" but "does the
code route": every hazard is a way for the wiring to look like it works while
quietly enforcing less than it claims. The Tier-B mechanics themselves --
round-robin target selection, the fixable/advisory partition, per-rule repair
discipline -- live in `repair.py` and are tested in `test_repair_seam.py`; this
file tests what the catalog HANDS to that machinery.
"""

from __future__ import annotations

from typing import ClassVar

import pytest
from bt5.core.context import DesignContext, HostId, Modality
from bt5.core.services import Services
from bt5.core.spec import (
    Breach,
    Citation,
    Direction,
    Enforcement,
    Evaluation,
    Evidence,
    LatticeTerms,
    LocalizationPolicy,
    RepairPolicy,
    Spec,
)
from bt5.core.types import Construct, Interval
from bt5.rules.catalog.e1_homopolymers import Homopolymers
from bt5.rules.catalog.e2_gc_band import GCBand
from bt5.rules.vendors import VendorSelection
from bt5.solver.catalog import RuleSet, build_rule_set, default_services
from conftest import context, linear_cds, slot


class Fake:
    """A minimal Spec. Registering a real rule for a test would pollute the
    process-global registry that every other test discovers from."""

    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "fake"
    evidence: ClassVar[Evidence] = Evidence.EVIDENCE_BACKED
    direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
    unit: ClassVar[str] = "count"
    citations: ClassVar[tuple[Citation, ...]] = (Citation("t", "https://example.org"),)
    last_verified: ClassVar[str] = "2026-08-31"
    weight_provenance: ClassVar[str] = "test double"
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.0
    steering_weight: ClassVar[float] = 0.0
    band: ClassVar[tuple[float, float] | None] = None
    cost_class: ClassVar[str] = "cheap"
    conflicts_with: ClassVar[tuple[str, ...]] = ()
    param_schema: ClassVar[dict[str, object]] = {"type": "object"}
    brief_ref: ClassVar[str] = "test"
    engine_calibration: ClassVar[str | None] = None

    def __init__(
        self,
        spec_id: str,
        *,
        enforcement: Enforcement = Enforcement.HARD_REPAIR,
        localization: LocalizationPolicy = LocalizationPolicy.WINDOW_MINUS_1,
        repair: RepairPolicy = RepairPolicy.SINGLE_PASS,
        window: int | None = None,
        breaches: tuple[Breach, ...] = (),
        passes: bool = False,
        gates: bool = True,
        hard_in: Modality | None = None,
    ) -> None:
        self.id = spec_id
        self.enforcement = enforcement
        self.localization = localization
        self.repair = repair
        if window is not None:
            self.window = window
        self._breaches = breaches
        self._passes = passes
        self._gates = gates
        self._hard_in = hard_in
        self.calls = 0

    def gate(self, slot_: object) -> bool:
        return self._gates

    def enforcement_for(self, slot_: object) -> Enforcement:
        if self._hard_in is not None:
            modality = getattr(slot_, "modality", None)
            return Enforcement.HARD_REPAIR if modality is self._hard_in else Enforcement.SOFT
        return self.enforcement

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        self.calls += 1
        return Evaluation(self.id, self._passes, 0.0, self._breaches)

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        return None


def breach(spec_id: str, start: int, magnitude: float, *, fixable: bool = True) -> Breach:
    return Breach(
        spec_id=spec_id,
        interval=Interval(start, start + 6),
        magnitude=magnitude,
        message=f"{spec_id} at {start}",
        fixable_by_codon_choice=fixable,
    )


def ruleset(*specs: Spec, ctx: DesignContext | None = None, svc: Services | None = None) -> RuleSet:
    services = svc if svc is not None else default_services(seed=1, autoload_fold=False)
    return RuleSet(specs=tuple(specs), ctx=ctx or context(), svc=services)


CDS = "ATG" + "GCT" * 20 + "TAA"
PROTEIN = "M" + "A" * 20


class TestTheSlotLoopIsTheRulesJob:
    """Hazard 3. Rules iterate `ctx.active_slots` themselves, so a caller that
    also looped would report every finding once per slot."""

    def test_a_spec_is_evaluated_once_per_construct_not_once_per_slot(self) -> None:
        spy = Fake("spy")
        rules = ruleset(
            spy,
            ctx=context(
                slot(),
                slot("producer", HostId.HEK293, Modality.LENTIVIRAL, 1),
                slot("target", HostId.HUMAN, Modality.GENOME_INTEGRATED, 1),
            ),
        )
        rules.findings(linear_cds(CDS, PROTEIN))
        assert spy.calls == 1

    def test_a_spec_no_active_slot_gates_is_recorded_not_dropped_in_silence(self) -> None:
        rules = build_rule_set(
            context(slot("producer", HostId.HEK293, Modality.IVT_MRNA, 1)),
            default_services(seed=1, autoload_fold=False),
        )
        assert "e2_gc_band" not in {s.id for s in rules.specs}
        assert "e2_gc_band" in rules.gated_out

    def test_enforcement_is_resolved_per_slot_not_read_off_the_classvar(self) -> None:
        """D4 and D6 return HARD_REPAIR in some modalities and SOFT in others.
        Reading the ClassVar would demote a lentiviral polyA signal to a
        weighted preference."""
        d4 = Fake(
            "d4_like",
            enforcement=Enforcement.SOFT,
            hard_in=Modality.LENTIVIRAL,
            breaches=(breach("d4_like", 3, 1.0),),
        )
        c = linear_cds(CDS, PROTEIN)

        bacterial = ruleset(d4, ctx=context(slot()))
        assert bacterial.findings(c).repairable == ()

        lentiviral = ruleset(
            d4, ctx=context(slot("producer", HostId.HEK293, Modality.LENTIVIRAL, 1))
        )
        assert len(lentiviral.findings(c).repairable) == 1


class TestRoutingToTheSolver:
    """What the breach finder must and must not hand to `repair()`."""

    def test_a_passing_rule_contributes_nothing_to_the_solver(self) -> None:
        """`Evaluation.passes` is NOT `not breaches`: E5 passes on a warn-band
        finding, E7 under its hard tract, F1 under its hard length. Chasing those
        stagnates the search and reports a design the catalog accepts as
        infeasible."""
        warn = Fake("e5_like", passes=True, breaches=(breach("e5_like", 6, 0.5),))
        found = ruleset(warn).findings(linear_cds(CDS, PROTEIN))
        assert found.repairable == ()
        assert found.hard_check == ()

    def test_an_unfixable_hard_repair_breach_reaches_the_finder_for_partitioning(self) -> None:
        """It belongs in the finder output so `repair()._partition` can carry it
        on `RepairOutcome.advisory` -- a polyA hexamer in the user's own LTR is
        reported, not chased. Keeping it OUT of the finder is the old bug where an
        unfixable backbone breach aborted the whole pass."""
        rules = ruleset(Fake("d4_like", breaches=(breach("d4_like", 0, 2.0, fixable=False),)))
        found = rules.findings(linear_cds(CDS, PROTEIN))
        assert len(found.repairable) == 1
        assert found.repairable[0].fixable_by_codon_choice is False

    def test_hard_check_findings_never_reach_the_finder(self) -> None:
        """HARD_CHECK means real, reported, never chased -- an over-length
        fragment, an ITR palindrome. It routes to `hard_check`, not the solver."""
        e9 = Fake(
            "e9_like",
            enforcement=Enforcement.HARD_CHECK,
            breaches=(breach("e9_like", 0, 1.0, fixable=False),),
        )
        rules = ruleset(e9)
        c = linear_cds(CDS, PROTEIN)
        assert rules.breach_finder()(c) == ()
        assert len(rules.advise()(c)) == 1

    def test_soft_findings_are_neither_repaired_nor_advised(self) -> None:
        soft = Fake("soft", enforcement=Enforcement.SOFT, breaches=(breach("soft", 0, 9.0),))
        found = ruleset(soft).findings(linear_cds(CDS, PROTEIN))
        assert found.repairable == ()
        assert found.hard_check == ()
        assert len(found.evaluations) == 1


class TestThePolicies:
    """Hazard 4. `repair()` solves it per rule; the catalog supplies the map."""

    def test_each_hard_repair_rule_gets_its_own_localization(self) -> None:
        rules = ruleset(
            Fake("wide", localization=LocalizationPolicy.WHOLE_SCOPE),
            Fake("paired", localization=LocalizationPolicy.PAIRED_SEGMENTS),
        )
        pols = rules.policies()
        assert pols["wide"].localization is LocalizationPolicy.WHOLE_SCOPE
        assert pols["paired"].localization is LocalizationPolicy.PAIRED_SEGMENTS

    def test_a_rules_instance_window_travels_in_its_policy(self) -> None:
        rules = ruleset(Fake("windowed", window=80))
        assert rules.policies(default_window=50)["windowed"].window == 80

    def test_a_rule_without_a_window_falls_back_to_the_default(self) -> None:
        rules = ruleset(Fake("plain"))
        assert rules.policies(default_window=42)["plain"].window == 42

    def test_each_rules_repair_discipline_travels_per_rule(self) -> None:
        """No global escalation: a FIXED_POINT rule declares it on itself, which
        is where CLAUDE.md 3.6 puts the splice-removal requirement."""
        rules = ruleset(
            Fake("single", repair=RepairPolicy.SINGLE_PASS),
            Fake("fixed", repair=RepairPolicy.FIXED_POINT),
        )
        pols = rules.policies()
        assert pols["single"].repair is RepairPolicy.SINGLE_PASS
        assert pols["fixed"].repair is RepairPolicy.FIXED_POINT

    def test_only_hard_repair_rules_get_a_policy(self) -> None:
        rules = ruleset(
            Fake("hard", enforcement=Enforcement.HARD_REPAIR),
            Fake("soft", enforcement=Enforcement.SOFT),
            Fake("check", enforcement=Enforcement.HARD_CHECK),
        )
        assert set(rules.policies()) == {"hard"}

    def test_the_shipped_hard_repair_rules_declare_three_localizations(self) -> None:
        """The reason a per-rule map is needed rather than one global policy: one
        value gives at least two of these the wrong repair window."""
        rules = build_rule_set(context(slot()), default_services(seed=1, autoload_fold=False))
        declared = {p.localization for p in rules.policies().values()}
        assert len(declared) >= 3


class TestTheOracleGetsTheSameNumbers:
    """The #59 failure was the oracle and a rule enforcing different contracts."""

    def test_gc_bounds_come_from_the_e2_instance_not_its_classvar(self) -> None:
        rules = build_rule_set(context(slot()), default_services(seed=1, autoload_fold=False))
        e2 = next(s for s in rules.specs if s.id == "e2_gc_band")
        assert rules.oracle_bounds().gc_bounds == (e2.gc_min, e2.gc_max)  # type: ignore[attr-defined]
        assert rules.oracle_bounds().gc_bounds != GCBand.band

    def test_the_band_narrows_with_the_selection(self) -> None:
        svc = default_services(seed=1, autoload_fold=False)
        gblocks = build_rule_set(context(slot()), svc)
        twist = build_rule_set(
            context(slot()), svc, vendors=VendorSelection.of("twist_gene_fragment")
        )
        assert gblocks.oracle_bounds().gc_bounds != twist.oracle_bounds().gc_bounds

    def test_an_adapter_on_selection_disarms_i7_rather_than_half_arming_it(self) -> None:
        rules = build_rule_set(
            context(slot()),
            default_services(seed=1, autoload_fold=False),
            vendors=VendorSelection.of("twist_gene_fragment_adapter_on"),
        )
        assert rules.oracle_bounds().gc_bounds is None

    def test_no_e2_means_no_band(self) -> None:
        assert ruleset(Fake("x")).oracle_bounds().gc_bounds is None

    def test_forbidden_motifs_come_from_the_hard_lattice_rules_only(self) -> None:
        rules = build_rule_set(context(slot()), default_services(seed=1, autoload_fold=False))
        motifs = rules.forbidden()
        assert "GAATTC" in motifs  # D1, EcoRI
        assert "A" * (Homopolymers().max_at_run + 1) in motifs
        assert all(len(m) >= 6 for m in motifs)


class TestBuildingTheSet:
    def test_the_default_selection_reaches_every_rule_that_accepts_one(self) -> None:
        rules = build_rule_set(
            context(slot()),
            default_services(seed=1, autoload_fold=False),
            vendors=VendorSelection.of("twist_gene_fragment"),
        )
        chosen = {
            s.vendors.label  # type: ignore[attr-defined]
            for s in rules.specs
            if hasattr(s, "vendors")
        }
        assert chosen == {"twist_gene_fragment"}

    def test_an_override_naming_no_parameter_says_which_rule_and_what_exists(self) -> None:
        with pytest.raises(ValueError, match="e2_gc_band"):
            build_rule_set(
                context(slot()),
                default_services(seed=1, autoload_fold=False),
                overrides={"e2_gc_band": {"nonsense": 1}},
            )

    def test_an_override_is_applied(self) -> None:
        rules = build_rule_set(
            context(slot()),
            default_services(seed=1, autoload_fold=False),
            overrides={"e2_gc_band": {"gc_min": 0.35, "gc_max": 0.65}},
        )
        assert rules.oracle_bounds().gc_bounds == (0.35, 0.65)

    def test_a_calibrated_rule_with_no_engine_is_recorded_rather_than_vanishing(self) -> None:
        rules = build_rule_set(context(slot()), default_services(seed=1, autoload_fold=False))
        assert "b1_five_prime" in rules.unrunnable

    def test_services_never_fabricate_a_folding_engine(self) -> None:
        assert default_services(seed=1, autoload_fold=False).fold is None
