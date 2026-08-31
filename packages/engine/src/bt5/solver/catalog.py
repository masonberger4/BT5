"""The catalog, wired -- the one place a rule is actually run.

Until this module existed, `grep -rn "\\.evaluate(" packages/engine/src/` returned
nothing: fifteen rules and every vendor profile were exercised only from tests,
and the `BreachFinder` seam in `repair.py` had no constructor. A rule that no
production code path calls is a rule that cannot refuse anything.

FOUR THINGS MUST COME FROM ONE PLACE or they can disagree, and one of them
already did (#59): the Tier-A forbidden set, the Tier-B breach finder, the
repair policies, and the parameters handed to the independent validator. E2
gates the GC of an ordered fragment against the selected vendor's band while
invariant I7 gates the same span against `gc_bounds`; supply those from two
sources and the optimizer and its own oracle enforce different contracts, in
both directions. So this module is one object built from one `VendorSelection`,
not a bag of functions each resolving its own defaults.

WHAT THIS MODULE DOES NOT DO. It does not loop over context slots -- rules do
that themselves (`b1:219`, `d4:237`, `f3:174`), so a caller that also looped
would double-count every finding. It does not decide weights; the weighted sum
is M3's and only ever sees SOFT rules. And it does not send an unfixable
finding to the solver: `Breach.fixable_by_codon_choice` exists precisely so a
uAUG in the user's own 5'UTR or an over-length fragment is reported rather than
chased through a mutation space that cannot contain a fix.
"""

from __future__ import annotations

import inspect
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, cast

from bt5.core.registry import all_specs, check_engine_calibration, discover
from bt5.core.services import Services, TableProvider
from bt5.core.spec import (
    Breach,
    Enforcement,
    Evaluation,
    LocalizationPolicy,
    RepairPolicy,
    Spec,
)
from bt5.rules.vendors import DEFAULT_SELECTION, VendorSelection
from bt5.solver.pipeline import OptimizeResult, optimize
from bt5.solver.reference import CodonScorer
from bt5.solver.repair import (
    NO_RULES,
    Assembler,
    BreachCost,
    BreachFinder,
    no_rules,
)

if TYPE_CHECKING:
    from bt5.codon.tables import NcbiGeneticCode
    from bt5.core.context import ContextSlot, DesignContext
    from bt5.core.services import FoldEngine
    from bt5.core.types import Construct

#: Which enforcement wins when a breach names no slot. A finding that is hard in
#: ANY applicable context is hard for the design, because the three contexts are
#: never collapsed -- one sequence has to survive propagation, packaging and the
#: target cell, so taking the weakest reading would let a lentiviral polyA breach
#: be scored away by the E. coli slot's opinion of it.
_SEVERITY: Mapping[Enforcement, int] = {
    Enforcement.HARD_LATTICE: 4,
    Enforcement.HARD_REPAIR: 3,
    Enforcement.HARD_CHECK: 2,
    Enforcement.SOFT: 1,
    Enforcement.REPORT_ONLY: 0,
}


@dataclass(frozen=True, slots=True)
class OracleBounds:
    """What `verify_construct` must be told so it checks what the rules enforced.

    ONE FIELD, AND THE TWO ABSENCES ARE THE INTERESTING PART.

    `max_homopolymer` is not derived, because arming it would be strictly worse
    than what already happens. E1 is HARD_LATTICE and publishes
    `A*(max_at+1)` / `G*(max_gc+1)` through `lattice_terms` (`e1:181-190`); those
    reach the oracle as `forbidden`, and I6 proves them absent on both strands
    over the circular construct WHILE HONOURING `Construct.exempt`
    (`verify.py:148`). I8's homopolymer scan (`verify.py:304-310`) does not honour
    exemptions and takes a single ceiling for every base, where E1 bands A/T and
    G/C separately -- 9 and 5 under the gBlocks default. The looser number
    catches nothing I6 has not already caught; the tighter one refuses runs E1
    permits. There is no honest single value, and no gap for it to close.

    `max_repeat` is not derived either, and that one is measurable rather than
    arguable: I8's repeat scan (`verify.py:311-319`) walks every k-mer of the
    whole construct with no reference to `exempt` at all, so on any real
    lentiviral or AAV vector it fires on the identical LTRs or ITRs -- the
    `WHITELISTED_REPEAT` segments that exist precisely to be skipped. F1 gets
    this right by calling `both_arms_exempt` (`f1:206`); I8 has no equivalent.
    Teaching it one is an oracle change and carries its own label.
    """

    gc_bounds: tuple[float, float] | None = None


@dataclass(frozen=True, slots=True)
class Findings:
    """One evaluation of the whole catalog, already routed.

    `repairable` is what Tier B may chase. `advisory` is real and reported and
    NOT chased -- routing it into the solver is what makes a design with a
    perfectly good CDS report itself infeasible because the vector is 200 bp
    too long for the chosen vendor.
    """

    repairable: tuple[Breach, ...] = ()
    advisory: tuple[Breach, ...] = ()
    evaluations: tuple[Evaluation, ...] = ()


@dataclass(frozen=True)
class RuleSet:
    """The catalog, instantiated once and bound to one design context.

    Instances, not classes: a rule's thresholds come from its constructor (E2's
    band from the selection, E1's run limits, E5's length), so `type[Spec]` is
    not enough to evaluate anything, and building a fresh instance per repair
    iteration would recompute the vendor intersection a thousand times.
    """

    specs: tuple[Spec, ...]
    ctx: DesignContext
    svc: Services
    vendors: VendorSelection = DEFAULT_SELECTION
    #: spec ids no active slot gates. Recorded, not dropped in silence: an
    #: objective the user believes is on and that simply did not apply is a
    #: different thing from one that passed, and only the report can say which.
    gated_out: tuple[str, ...] = ()
    #: spec ids whose thresholds are calibrated against a folding engine that is
    #: not available. `check_engine_calibration` returns these rather than
    #: raising, because absence is a degradation to report.
    unrunnable: tuple[str, ...] = ()

    # -- Tier A -------------------------------------------------------------

    def forbidden(self) -> tuple[str, ...]:
        """Forward motifs for the automaton, from the HARD_LATTICE rules only.

        The solver closes this set under reverse complement, so a rule lists one
        strand and gets both. Restricted to HARD_LATTICE deliberately: that
        enforcement class MEANS "guaranteed by construction in the Tier-A DP",
        and letting a SOFT rule contribute here would turn a weighted preference
        into an absolute guarantee with no way for the user to trade it away.
        """
        out: list[str] = []
        for spec in self.specs:
            if spec.enforcement is not Enforcement.HARD_LATTICE:
                continue
            terms = spec.lattice_terms(self.ctx)
            if terms is not None:
                out.extend(terms.forbidden)
        return tuple(sorted(set(out)))

    # -- Tier B -------------------------------------------------------------

    def findings(self, c: Construct, specs: Sequence[Spec] | None = None) -> Findings:
        """Run rules against the assembled construct, once each, and route them.

        THE SOLVER CHASES ONLY BREACHES FROM A RULE THAT SAYS IT FAILED, and
        `Evaluation.passes` is emphatically not `not breaches`. Three of the four
        always-HARD_REPAIR rules deliberately emit sub-threshold findings while
        passing: E5 passes on `worst == "warn"` (`e5:259`), E7 on
        `worst <= hard_tract` (`e7:264`), F1 on `worst < hard_len` (`f1:251`).
        Handing those to `repair()` sets it chasing findings the rule itself
        calls acceptable -- and since no codon choice can clear a threshold that
        was never crossed, the search stagnates and raises
        `InfeasibleConstraints` on a design the catalog accepts.

        A passing rule's breaches are still REPORTED. "Not a failure" and "not
        worth mentioning" are different claims, and a warn-band repeat is exactly
        what a user wants to see before paying for a tube.
        """
        repairable: list[Breach] = []
        advisory: list[Breach] = []
        evaluations: list[Evaluation] = []
        for spec in self.specs if specs is None else specs:
            ev = spec.evaluate(c, self.ctx, self.svc)
            evaluations.append(ev)
            for breach in ev.breaches:
                enforcement = self._enforcement_of(spec, breach)
                if not enforcement.is_hard:
                    continue  # SOFT and REPORT_ONLY are M3's; see CLAUDE.md 3.5
                if (
                    enforcement is Enforcement.HARD_REPAIR
                    and breach.fixable_by_codon_choice
                    and not ev.passes
                ):
                    repairable.append(breach)
                else:
                    advisory.append(breach)
        return Findings(tuple(repairable), tuple(advisory), tuple(evaluations))

    def breach_finder(self) -> BreachFinder:
        """The HARD_REPAIR, codon-fixable breaches. This is what Tier B repairs.

        SCOPED TO `repair_specs()`, AND THAT IS A PERFORMANCE DECISION WITH A
        CORRECTNESS ARGUMENT UNDER IT. This callable runs once per candidate --
        up to 256 per iteration in the guided-random branch, for up to a thousand
        iterations -- so every rule it touches is paid for a quarter of a million
        times. E8 builds a k-mer index over the whole construct; F2 extends seeds
        through mismatches; B1 folds. None of them can ever return a repairable
        breach, because a SOFT rule's findings are the weighted sum's business
        and the weighted sum is not in this loop. Evaluating them here buys
        nothing and costs the interactive budget: with all thirteen rules a
        repetitive 500 aa design did not converge in two minutes; with the five
        that can contribute, it is seconds.

        `advise()` still runs the full set, once, on the construct that actually
        ships -- so nothing is dropped from the report, only from the inner loop.
        """
        specs = self.repair_specs()
        return lambda c: self.findings(c, specs).repairable

    def advise(self) -> Callable[[Construct], tuple[Breach, ...]]:
        """The hard findings no codon can fix. Reported, never chased.

        Every rule, not just the repairable ones, and evaluated on the final
        construct: an advisory describes what shipped.
        """
        return lambda c: self.findings(c).advisory

    def cost(self) -> BreachCost:
        """A FRESH cost function: the scale freezes per repair run, not per RuleSet."""
        return BreachCost()

    def localization_for(self, spec_id: str) -> LocalizationPolicy:
        """The rule's own declared policy, or the generic default.

        Per breach rather than per run because the four HARD_REPAIR rules
        currently declare three different policies -- E2 `WINDOW_MINUS_1`, E5 and
        F1 `PAIRED_SEGMENTS`, E7 `WHOLE_SCOPE`. One global value gives at least
        two of them the wrong repair window on every run: a repeat pair is only
        fixable by editing one copy, and widening it by `window - 1` instead
        searches a region that does not contain the fix.
        """
        for spec in self.specs:
            if spec.id == spec_id:
                return spec.localization
        return LocalizationPolicy.WINDOW_MINUS_1

    def repair_policy(self, requested: RepairPolicy = RepairPolicy.FIXED_POINT) -> RepairPolicy:
        """Escalate to FIXED_POINT if any kept rule needs it. NEVER de-escalate.

        Unlike localisation this cannot be per breach: it controls when the one
        shared search loop stops, so there is exactly one value and it has to be
        the safe one. CLAUDE.md section 3.6 fixes the direction -- a single pass
        over a set containing a FIXED_POINT rule ships a construct whose cryptic
        donors were removed INTO NEW DONORS, and the validator passes it because
        the specific 9-mer it was told about is gone. The reverse mistake, running
        a SINGLE_PASS rule to convergence, costs iterations and nothing else.

        This is also why the join only ever goes up. Every rule shipped today
        declares SINGLE_PASS while `repair()` defaults to FIXED_POINT, so a plain
        "join over the contributing rules" would compute SINGLE_PASS and quietly
        downgrade Tier B from iterate-to-convergence to stop-on-first-stall -- a
        regression introduced by wiring the catalog in, which is the opposite of
        the point.
        """
        if requested is RepairPolicy.FIXED_POINT:
            return RepairPolicy.FIXED_POINT
        repairing = self.repair_specs()
        if any(spec.repair is RepairPolicy.FIXED_POINT for spec in repairing):
            return RepairPolicy.FIXED_POINT
        # A rule declaring SINGLE_PASS is making a claim about ITSELF: "fixing one
        # of my breaches cannot create another of mine." It says nothing about
        # whether recoding a window to bring GC into E2's band creates a 20 bp
        # repeat E5 refuses. Two such rules in one loop invalidate each other's
        # claim, so more than one contributor promotes unconditionally.
        if len(repairing) > 1:
            return RepairPolicy.FIXED_POINT
        return requested

    def repair_specs(self) -> tuple[Spec, ...]:
        """The rules that can put work into Tier B, in any active slot."""
        return tuple(
            spec
            for spec in self.specs
            if any(
                spec.enforcement_for(slot) is Enforcement.HARD_REPAIR
                for slot in self.ctx.active_slots
            )
        )

    def demands_fixed_point(self) -> tuple[str, ...]:
        """Which rules force the escalation, for the report."""
        return tuple(
            sorted(spec.id for spec in self.specs if spec.repair is RepairPolicy.FIXED_POINT)
        )

    # -- The validator ------------------------------------------------------

    def oracle_bounds(self) -> OracleBounds:
        """Numbers for `verify_construct`, taken from the rules it backstops.

        Passed as parameters, never resolved by the oracle itself: `verify.py`
        may not import `bt5.rules` (an AST check enforces it), and that
        independence is what stops the validator from sharing a code path with
        the scorer that guided the search.
        """
        return OracleBounds(gc_bounds=self._gc_bounds())

    def _gc_bounds(self) -> tuple[float, float] | None:
        """E2's own resolved band, so E2 and I7 cannot enforce different contracts.

        Read off the INSTANCE, never `Spec.band`: E2 says in as many words that
        its ClassVar (0.28, 0.80) is the loosest demonstrated envelope and that
        the gate is the selected vendors' intersection, computed per instance.

        Three ways this returns None, each of them a refusal to half-arm the
        oracle rather than an oversight:

        - no E2 at all (disabled, excluded, or gated off in every active slot --
          it gates off for IVT mRNA), so there is no band to enforce;
        - E2 could not resolve numbers;
        - an ADAPTER-ON selection. E2 measures the fragment the vendor
          synthesises, adapters included (`fragments()` splices them on); I7
          measures the designable span alone, because the oracle has no vendor
          data and its documented extension point for adapters is a parameter
          nobody passes yet. With adapters those are different bases and so
          different numbers, which is #59's contradiction from a third side.
          `VendorSelection.of()` guarantees an adapter-on selection is
          single-key, so this is exactly one configuration today.
        """
        for spec in self.specs:
            if spec.id != "e2_gc_band":
                continue
            lo = getattr(spec, "gc_min", None)
            hi = getattr(spec, "gc_max", None)
            if not (isinstance(lo, float) and isinstance(hi, float)):
                return None
            # getattr rather than an import: the solver stays free of a
            # compile-time dependency on the vendor catalogue. Pinned by test, so
            # a rename in M4 fails loudly instead of silently reading 0.
            adapters = getattr(getattr(spec, "vendors", None), "adapters", None)
            if getattr(adapters, "total", 0):
                return None
            return (lo, hi)
        return None

    # -- internals ----------------------------------------------------------

    def _enforcement_of(self, spec: Spec, breach: Breach) -> Enforcement:
        """How hard is this finding, in this design's contexts?

        `enforcement_for` is not decoration: D4 (internal polyA) and D6 (non-B
        DNA) both return HARD_REPAIR in some modalities and SOFT in others, so
        reading the ClassVar would silently demote a lentiviral polyA signal --
        which raised expression 3-6.5x while cutting functional titer 8-9x -- to
        a weighted preference.
        """
        slots = self._gating_slots()
        if breach.slot_role is not None:
            for slot in slots:
                if slot.role == breach.slot_role:
                    return spec.enforcement_for(slot)
        if not slots:
            return spec.enforcement
        return max((spec.enforcement_for(s) for s in slots), key=lambda e: _SEVERITY[e])

    def _gating_slots(self) -> tuple[ContextSlot, ...]:
        return self.ctx.active_slots


def build_rule_set(
    ctx: DesignContext,
    svc: Services,
    *,
    vendors: VendorSelection = DEFAULT_SELECTION,
    overrides: Mapping[str, Mapping[str, object]] | None = None,
    include: Callable[[type[Spec]], bool] | None = None,
) -> RuleSet:
    """Discover, filter, calibrate, instantiate and gate the catalog.

    The order matters and each step earns its place:

    1. `discover()` walks `bt5.rules.catalog` so a new rule file needs no edit
       here. There is deliberately no committed catalog list.
    2. `default_enabled` drops FOLKLORE rules, which the contract forces off.
    3. `check_engine_calibration` RAISES when a rule's kcal/mol thresholds were
       measured on an engine other than the one running, and reports rather than
       runs them when no engine is available at all.
    4. Instantiation passes the one `VendorSelection` to every rule that accepts
       one -- detected from the signature, so adding an E-rule does not mean
       editing a list here.
    5. Gating is hoisted OUT of the evaluation loop: `gate()` reads only the
       slot, so it cannot change between repair iterations.
    """
    discover()
    catalog = [cls for cls in all_specs() if cls.default_enabled]
    if include is not None:
        catalog = [cls for cls in catalog if include(cls)]
    runnable = check_engine_calibration(catalog, svc.fold)
    unrunnable = tuple(sorted({cls.id for cls in catalog} - {cls.id for cls in runnable}))

    specs: list[Spec] = []
    gated_out: list[str] = []
    slots = ctx.active_slots
    for cls in runnable:
        spec = _instantiate(cls, vendors, (overrides or {}).get(cls.id))
        if any(spec.gate(slot) for slot in slots):
            specs.append(spec)
        else:
            gated_out.append(spec.id)

    return RuleSet(
        specs=tuple(specs),
        ctx=ctx,
        svc=svc,
        vendors=vendors,
        gated_out=tuple(sorted(gated_out)),
        unrunnable=unrunnable,
    )


def _instantiate(
    cls: type[Spec], vendors: VendorSelection, overrides: Mapping[str, object] | None
) -> Spec:
    """`cls(**params)`, with the shared vendor selection injected where accepted."""
    accepted = inspect.signature(cls.__init__).parameters
    kwargs: dict[str, object] = {}
    if "vendors" in accepted:
        kwargs["vendors"] = vendors
    for name, value in (overrides or {}).items():
        if name not in accepted:
            raise ValueError(
                f"{cls.id}: no parameter {name!r}; its param_schema advertises "
                f"{sorted(n for n in accepted if n != 'self')}"
            )
        # A selection arrives from JSON as a list of vendor keys, and
        # `require_selection` refuses a bare string on purpose.
        kwargs[name] = (
            VendorSelection.of(*(str(k) for k in value))
            if name == "vendors" and isinstance(value, list | tuple)
            else value
        )
    # `Spec` is a Protocol, so the class object is not statically instantiable
    # even though every registered rule is a concrete class.
    factory = cast("Callable[..., Spec]", cls)
    return factory(**kwargs)


def default_services(
    *,
    seed: int,
    fold: FoldEngine | None = None,
    autoload_fold: bool = True,
) -> Services:
    """Resolve the concrete providers a rule needs.

    THE IMPORTS ARE INSIDE THE FUNCTION, and that is not style. Importing this
    module must not drag the vector, codon and structure lanes into the solver's
    module graph: the whole point of `Services` is that a rule RECEIVES its
    providers rather than importing them, and a top-level import here would make
    the solver depend on M2, M5 and M6 at import time for a convenience it does
    not itself use. Every provider stays injectable.

    `fold` stays None when ViennaRNA is absent -- never a stub. Every threshold
    in BT5 is a kcal/mol number, so a stub returning plausible energies would
    flow through the scorers, the null and the percentile unchallenged and come
    out the far end as a confident rank.

    THE THREE CASTS ARE REAL DEFECTS, not type-checker appeasement, and they are
    narrow on purpose so that fixing one removes exactly one. Nothing in `src/`
    had ever constructed a `Services` before this module, so no lane's concrete
    provider had ever been checked against the protocol it claims to implement,
    and all three fail:

      - `ViennaFold.version` is a read-only property where `FoldEngine.version`
        is a settable `ClassVar[str]` (M6);
      - `NcbiGeneticCode.table_id` is a read-only property where
        `GeneticCode.table_id` is a plain attribute (M5);
      - `FileTableProvider.usage` returns a `CodonUsage`, not the
        `Mapping[str, float]` the protocol declares, and is `@cache`-wrapped on
        top (M5). `packages/engine/tests/rules/conftest.py` already papered over
        this one with a test double returning `dict[str, float]`, so the two
        implementations of `usage()` do not even agree with each other.

    All three are other lanes' to fix, and correcting the protocols instead
    would be a `core/` amendment. They are runtime-harmless for every rule
    shipped today -- nothing calls `svc.tables.usage()` -- but the third would
    break the first rule that does, which is why it is written down here rather
    than smoothed over.
    """
    import numpy as np

    from bt5.codon.tables import FileTableProvider
    from bt5.vector.kmers import ConstructKmerIndex

    engine = fold
    if engine is None and autoload_fold:
        from bt5.structure.vienna import load_fold_engine

        engine = cast("FoldEngine | None", load_fold_engine())

    return Services(
        fold=engine,
        kmer=ConstructKmerIndex,
        tables=cast("TableProvider", FileTableProvider()),
        rng=np.random.default_rng(seed),
    )


def optimize_with(
    rules: RuleSet,
    protein: str,
    code: NcbiGeneticCode,
    *,
    assemble: Assembler,
    score: CodonScorer | None = None,
    original_backbone: Construct | None = None,
    repair_policy: RepairPolicy = RepairPolicy.FIXED_POINT,
    gc_window: int = 50,
    left_flank: str = "",
    right_flank: str = "",
    seed: int = 0,
    table_id: int | None = None,
) -> OptimizeResult:
    """`optimize()` with every argument derived from ONE rule set.

    This is the call that makes the catalog real, and the reason it exists here
    rather than in `pipeline` is that all seven derived arguments have to come
    from the same instantiated rules. Supply the motifs from one place and the
    validator's GC band from another and you get #59 again: E2 gating a fragment
    against the selected vendor's band while I7 gates the same span against a
    different one, refusing constructs the rules pass and passing constructs the
    rules refuse.

    Pass `original_backbone` -- `vector.assemble()` returns it as
    `Assembly.reference` -- to arm I9. Without it the invariant that proves BT5
    did not touch a single backbone base does not run.
    """
    bounds = rules.oracle_bounds()
    return optimize(
        protein,
        code,
        assemble=assemble,
        find_breaches=rules.breach_finder(),
        forbidden=rules.forbidden(),
        score=score,
        gc_bounds=bounds.gc_bounds,
        gc_window=gc_window,
        localization=rules.localization_for,
        repair_policy=rules.repair_policy(repair_policy),
        cost=rules.cost(),
        advise=rules.advise(),
        original_backbone=original_backbone,
        left_flank=left_flank,
        right_flank=right_flank,
        seed=seed,
        table_id=table_id,
    )


__all__ = [
    "NO_RULES",
    "BreachCost",
    "Findings",
    "OracleBounds",
    "RuleSet",
    "build_rule_set",
    "default_services",
    "no_rules",
    "optimize_with",
]
