"""The wiring that makes the rule catalog run.

Before this module's subject existed, `grep -rn "\\.evaluate(" src/` returned
nothing. Most of what is tested here is therefore not "does the code work" but
"does the code route", because every one of these four hazards is a way for the
wiring to look like it works while quietly enforcing less than it claims.
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
from bt5.solver.repair import BreachCost
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
        breaches: tuple[Breach, ...] = (),
        passes: bool = False,
        gates: bool = True,
        hard_in: Modality | None = None,
    ) -> None:
        self.id = spec_id
        self.enforcement = enforcement
        self.localization = localization
        self.repair = repair
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
        """'Did not apply' and 'passed' are different claims, and only the
        report can tell them apart."""
        rules = build_rule_set(
            context(slot("producer", HostId.HEK293, Modality.IVT_MRNA, 1)),
            default_services(seed=1, autoload_fold=False),
        )
        assert "e2_gc_band" not in {s.id for s in rules.specs}
        assert "e2_gc_band" in rules.gated_out

    def test_enforcement_is_resolved_per_slot_not_read_off_the_classvar(self) -> None:
        """D4 and D6 return HARD_REPAIR in some modalities and SOFT in others.
        Reading the ClassVar would demote a lentiviral polyA signal -- which cut
        functional titer 8-9x -- to a weighted preference."""
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


class TestRoutingAwayFromTheSolver:
    """Hazard 7 and the `passes` trap. What the solver must NOT be handed."""

    def test_a_passing_rules_sub_threshold_breaches_are_reported_not_repaired(self) -> None:
        """`Evaluation.passes` is NOT `not breaches`: E5 passes on a warn-band
        finding, E7 under its hard tract, F1 under its hard length. Chasing those
        stagnates the search and reports a design the catalog accepts as
        infeasible."""
        warn = Fake("e5_like", passes=True, breaches=(breach("e5_like", 6, 0.5),))
        found = ruleset(warn).findings(linear_cds(CDS, PROTEIN))
        assert found.repairable == ()
        assert len(found.advisory) == 1

    def test_hard_check_findings_never_reach_the_finder(self) -> None:
        e9 = Fake(
            "e9_like",
            enforcement=Enforcement.HARD_CHECK,
            breaches=(breach("e9_like", 0, 1.0, fixable=False),),
        )
        rules = ruleset(e9)
        c = linear_cds(CDS, PROTEIN)
        assert rules.breach_finder()(c) == ()
        assert len(rules.advise()(c)) == 1

    def test_an_unfixable_breach_never_reaches_the_finder(self) -> None:
        """Routing one into the solver exhausts the mutation space chasing a fix
        that does not exist, and reports a fine design infeasible."""
        rules = ruleset(Fake("f", breaches=(breach("f", 0, 2.0, fixable=False),)))
        found = rules.findings(linear_cds(CDS, PROTEIN))
        assert found.repairable == ()
        assert len(found.advisory) == 1

    def test_soft_findings_are_neither_repaired_nor_advised(self) -> None:
        soft = Fake("soft", enforcement=Enforcement.SOFT, breaches=(breach("soft", 0, 9.0),))
        found = ruleset(soft).findings(linear_cds(CDS, PROTEIN))
        assert found.repairable == ()
        assert found.advisory == ()
        assert len(found.evaluations) == 1

    def test_the_finder_skips_rules_that_cannot_contribute(self) -> None:
        """The finder runs once per candidate, up to 256 per iteration. A SOFT
        rule cannot return a repairable breach, so evaluating it here is a
        quarter of a million wasted k-mer indexes."""
        soft = Fake("soft", enforcement=Enforcement.SOFT)
        hard = Fake("hard")
        rules = ruleset(soft, hard)
        rules.breach_finder()(linear_cds(CDS, PROTEIN))
        assert soft.calls == 0
        assert hard.calls == 1


class TestTheCostReduction:
    """Hazard 2. `sum(b.magnitude)` mixed nucleotides with GC fractions."""

    def test_clearing_a_breach_beats_shrinking_one_in_another_rules_units(self) -> None:
        cost = BreachCost()
        cost([breach("nt", 0, 4.0), breach("gc", 9, 0.05)])  # freeze the scale
        cleared = cost([breach("nt", 0, 4.0)])
        shrunk = cost([breach("nt", 0, 3.0), breach("gc", 9, 0.04)])
        assert cleared < shrunk

    def test_a_small_magnitude_in_its_own_units_is_not_starved_by_a_large_one(self) -> None:
        """Under `sum(magnitude)` a GC delta of 0.05 is invisible beside a
        repeat finding of 4.0, so the search never works on it."""
        cost = BreachCost()
        cost([breach("nt", 0, 4.0), breach("gc", 9, 0.05)])
        assert cost.normalised(breach("gc", 9, 0.05)) == pytest.approx(1.0)
        assert cost.normalised(breach("nt", 0, 4.0)) == pytest.approx(1.0)
        assert cost.normalised(breach("gc", 9, 0.025)) == pytest.approx(0.5)

    def test_the_scale_does_not_move_between_iterations(self) -> None:
        """A moving normaliser makes the objective non-stationary: the same two
        states compare one way early and the other way late, so `best_cost`
        stops being a bound the search can descend."""
        cost = BreachCost()
        start = [breach("a", 0, 10.0), breach("b", 9, 1.0)]
        cost(start)
        later = [breach("a", 0, 1.0)]
        cost(later)  # 'a' is now the only rule; a moving scale would rescale it
        assert cost(later) == cost(later)
        assert cost.normalised(breach("a", 0, 1.0)) == pytest.approx(0.1)

    def test_a_rule_first_seen_late_is_recorded_rather_than_rescaling_the_run(self) -> None:
        cost = BreachCost()
        cost([breach("a", 0, 2.0)])
        cost([breach("a", 0, 2.0), breach("late", 9, 1.0)])
        assert cost.late == ("late",)

    def test_an_empty_breach_set_is_the_minimum(self) -> None:
        cost = BreachCost()
        cost([breach("a", 0, 1.0)])
        assert cost(()) < cost([breach("a", 0, 0.001)])


class TestThePolicies:
    """Hazard 4. Two policies, joined in opposite directions."""

    def test_localization_is_looked_up_per_rule(self) -> None:
        rules = ruleset(
            Fake("wide", localization=LocalizationPolicy.WHOLE_SCOPE),
            Fake("paired", localization=LocalizationPolicy.PAIRED_SEGMENTS),
        )
        assert rules.localization_for("wide") is LocalizationPolicy.WHOLE_SCOPE
        assert rules.localization_for("paired") is LocalizationPolicy.PAIRED_SEGMENTS

    def test_an_unknown_spec_id_falls_back_to_the_generic_policy(self) -> None:
        assert ruleset().localization_for("nobody") is LocalizationPolicy.WINDOW_MINUS_1

    def test_the_shipped_hard_repair_rules_declare_three_different_policies(self) -> None:
        """The reason per-breach lookup is needed rather than tidy: one global
        value gives at least two of these the wrong repair window."""
        rules = build_rule_set(context(slot()), default_services(seed=1, autoload_fold=False))
        declared = {s.id: s.localization for s in rules.repair_specs()}
        assert len(set(declared.values())) >= 3

    def test_an_all_single_pass_catalog_does_not_downgrade_a_fixed_point_request(self) -> None:
        """Every rule shipped today declares SINGLE_PASS while `repair()`
        defaults to FIXED_POINT, so a plain min-join would silently weaken Tier B
        from iterate-to-convergence to stop-on-first-stall."""
        rules = ruleset(Fake("a", repair=RepairPolicy.SINGLE_PASS))
        assert rules.repair_policy(RepairPolicy.FIXED_POINT) is RepairPolicy.FIXED_POINT

    def test_one_rule_asking_for_fixed_point_escalates_a_single_pass_request(self) -> None:
        rules = ruleset(Fake("splice", repair=RepairPolicy.FIXED_POINT))
        assert rules.repair_policy(RepairPolicy.SINGLE_PASS) is RepairPolicy.FIXED_POINT

    def test_two_contributors_escalate_even_when_both_claim_single_pass(self) -> None:
        """SINGLE_PASS is a rule's claim about ITSELF. It says nothing about
        whether recoding a window for E2 creates a repeat E5 refuses."""
        rules = ruleset(
            Fake("a", repair=RepairPolicy.SINGLE_PASS),
            Fake("b", repair=RepairPolicy.SINGLE_PASS),
        )
        assert rules.repair_policy(RepairPolicy.SINGLE_PASS) is RepairPolicy.FIXED_POINT

    def test_a_lone_single_pass_rule_stays_single_pass_when_asked(self) -> None:
        rules = ruleset(Fake("a", repair=RepairPolicy.SINGLE_PASS))
        assert rules.repair_policy(RepairPolicy.SINGLE_PASS) is RepairPolicy.SINGLE_PASS


class TestTheOracleGetsTheSameNumbers:
    """The #59 failure was the oracle and a rule enforcing different contracts."""

    def test_gc_bounds_come_from_the_e2_instance_not_its_classvar(self) -> None:
        """E2 says in as many words that its ClassVar is the loosest demonstrated
        envelope and the gate is the selected vendors' intersection."""
        rules = build_rule_set(context(slot()), default_services(seed=1, autoload_fold=False))
        e2 = next(s for s in rules.specs if s.id == "e2_gc_band")
        assert rules.oracle_bounds().gc_bounds == (e2.gc_min, e2.gc_max)  # type: ignore[attr-defined]
        assert rules.oracle_bounds().gc_bounds != GCBand.band

    def test_the_band_narrows_with_the_selection(self) -> None:
        svc = default_services(seed=1, autoload_fold=False)
        gblocks = build_rule_set(context(slot()), svc)
        twist = build_rule_set(
            context(slot()),
            svc,
            vendors=VendorSelection.of("twist_gene_fragment"),
        )
        assert gblocks.oracle_bounds().gc_bounds != twist.oracle_bounds().gc_bounds

    def test_an_adapter_on_selection_disarms_i7_rather_than_half_arming_it(self) -> None:
        """E2 measures the fragment including adapters; I7 measures the
        designable span alone. With adapters those are different bases, so a
        shared number would be #59 from a third side."""
        rules = build_rule_set(
            context(slot()),
            default_services(seed=1, autoload_fold=False),
            vendors=VendorSelection.of("twist_gene_fragment_adapter_on"),
        )
        assert rules.oracle_bounds().gc_bounds is None

    def test_no_e2_means_no_band(self) -> None:
        assert ruleset(Fake("x")).oracle_bounds().gc_bounds is None

    def test_forbidden_motifs_come_from_the_hard_lattice_rules_only(self) -> None:
        """Letting a SOFT rule contribute here would turn a weighted preference
        into an absolute guarantee the user cannot trade away."""
        rules = build_rule_set(context(slot()), default_services(seed=1, autoload_fold=False))
        motifs = rules.forbidden()
        assert "GAATTC" in motifs  # D1, EcoRI
        e1 = Homopolymers()
        assert "A" * (e1.max_at_run + 1) in motifs
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
        """`check_engine_calibration` skips these silently; a constraint that
        disappears without comment is the same failure as a repair that never
        ran."""
        rules = build_rule_set(context(slot()), default_services(seed=1, autoload_fold=False))
        assert "b1_five_prime" in rules.unrunnable

    def test_services_never_fabricate_a_folding_engine(self) -> None:
        assert default_services(seed=1, autoload_fold=False).fold is None
