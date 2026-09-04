"""The catalog, wired -- the one place a rule is actually run.

Until this module existed, `grep -rn "\\.evaluate(" packages/engine/src/` returned
nothing: fifteen rules and every vendor profile were exercised only from tests,
and the `BreachFinder` seam in `repair.py` had no constructor. A rule that no
production code path calls is a rule that cannot refuse anything.

THIS IS THE CALLER `repair.py` WAS BUILT FOR. `RulePolicy`'s own docstring says
the per-rule policies must be *caller-supplied* -- because importing the registry
to read them would pull the entire rule catalog, Biopython and the vendor
registry into the solver, and because the values live on rule INSTANCES
(`e2.window` is `self.window`), not the classes. This module is that caller: the
one place #58 sanctions the solver lane reaching the rules, and it hands
`repair()` and `optimize()` a `forbidden` set, a `find_breaches`, a
`policies` map and a `gc_bounds` -- all derived from ONE `VendorSelection` so
they cannot disagree, which is how #59 happened when E2 and the oracle resolved
the same band from two places.

WHAT THIS MODULE DOES NOT DO. It does not loop over context slots -- rules do
that themselves (`b1:219`, `d4:237`, `f3:174`), so a caller that also looped
would double-count every finding. It does not decide weights; the weighted sum
is M3's and only ever sees SOFT rules. And it does not hand the solver a finding
no codon can fix: a HARD_CHECK (an over-length fragment, an ITR palindrome)
never enters the breach finder, and an unfixable HARD_REPAIR (a polyA hexamer in
the user's own LTR) is routed by `repair()`'s own `_partition` onto
`RepairOutcome.advisory`, not chased.
"""

from __future__ import annotations

import inspect
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, cast

from bt5.core.registry import all_specs, check_engine_calibration, discover
from bt5.core.services import Services, TableProvider
from bt5.core.spec import Breach, Enforcement, Evaluation, Spec
from bt5.rules.vendors import DEFAULT_SELECTION, VendorSelection
from bt5.solver.pipeline import OptimizeResult, optimize
from bt5.solver.reference import CodonScorer
from bt5.solver.repair import Assembler, BreachFinder, RulePolicy

if TYPE_CHECKING:
    from bt5.codon.tables import NcbiGeneticCode
    from bt5.core.context import DesignContext
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

#: The hexamer default for MOTIF_LEN_MINUS_1 localisation. Only rules that
#: localise by motif length read it (D4's polyA is a hexamer), and `repair()`'s
#: junction guard scans by the longest forbidden *pattern* regardless, so this is
#: a localisation width, never a correctness bound.
_DEFAULT_MOTIF_LEN = 6


@dataclass(frozen=True, slots=True)
class OracleBounds:
    """What `verify_construct` must be told so it checks what the rules enforced.

    ONE FIELD, AND THE TWO ABSENCES ARE THE INTERESTING PART.

    `max_homopolymer` is not derived, because arming it would be strictly worse
    than what already happens. E1 is HARD_LATTICE and publishes
    `A*(max_at+1)` / `G*(max_gc+1)` through `lattice_terms` (`e1:181-190`); those
    reach the oracle as `forbidden`, and I6 proves them absent on both strands
    over the circular construct WHILE HONOURING `Construct.exempt`
    (`verify.py:148`). I8's homopolymer scan does not honour exemptions and takes
    a single ceiling for every base, where E1 bands A/T and G/C separately. The
    looser number catches nothing I6 has not; the tighter one refuses runs E1
    permits.

    `max_repeat` is not derived either, and that one is measurable rather than
    arguable: I8's repeat scan walks every k-mer of the whole construct with no
    reference to `exempt`, so on any real lentiviral or AAV vector it fires on the
    identical LTRs or ITRs -- the `WHITELISTED_REPEAT` segments that exist
    precisely to be skipped. Filed as an oracle follow-up.
    """

    gc_bounds: tuple[float, float] | None = None


@dataclass(frozen=True, slots=True)
class Findings:
    """One evaluation of the whole catalog, already routed.

    `repairable` is what the breach finder returns -- every HARD_REPAIR breach
    from a rule that FAILED, fixable or not; `repair()` partitions those, chasing
    the fixable ones and carrying the rest on `RepairOutcome.advisory`.
    `hard_check` is the HARD_CHECK family (over-length fragment, ITR palindrome):
    real, reported, and NEVER handed to the solver, which would exhaust the
    mutation space chasing a fix no codon can make.
    """

    repairable: tuple[Breach, ...] = ()
    hard_check: tuple[Breach, ...] = ()
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

    def findings(self, c: Construct) -> Findings:
        """Run every kept rule against the assembled construct, once each, routed.

        THE SOLVER CHASES ONLY BREACHES FROM A RULE THAT SAYS IT FAILED, and
        `Evaluation.passes` is emphatically not `not breaches`. Three of the four
        always-HARD_REPAIR rules deliberately emit sub-threshold findings while
        passing: E5 passes on `worst == "warn"` (`e5:259`), E7 on
        `worst <= hard_tract` (`e7:264`), F1 on `worst < hard_len` (`f1:251`).
        Handing those to `repair()` sets it chasing findings the rule itself
        calls acceptable -- and since no codon choice can clear a threshold that
        was never crossed, the search stagnates and raises `InfeasibleConstraints`
        on a design the catalog accepts.

        HARD_REPAIR breaches -- fixable or not -- go to `repairable`; `repair()`'s
        own `_partition` splits them. HARD_CHECK breaches go to `hard_check` and
        never near the solver. SOFT and REPORT_ONLY are M3's business.
        """
        repairable: list[Breach] = []
        hard_check: list[Breach] = []
        evaluations: list[Evaluation] = []
        for spec in self.specs:
            ev = spec.evaluate(c, self.ctx, self.svc)
            evaluations.append(ev)
            if ev.passes:
                continue  # the rule's own verdict; a warn-band finding is not a failure
            for breach in ev.breaches:
                enforcement = self._enforcement_of(spec, breach)
                if enforcement is Enforcement.HARD_REPAIR:
                    repairable.append(breach)
                elif enforcement is Enforcement.HARD_CHECK:
                    hard_check.append(breach)
                # SOFT / REPORT_ONLY: not the solver's business (CLAUDE.md 3.5)
        return Findings(tuple(repairable), tuple(hard_check), tuple(evaluations))

    def breach_finder(self) -> BreachFinder:
        """The HARD_REPAIR breaches Tier B repairs.

        Both fixable and unfixable -- `repair()._partition` routes the fixable
        ones into the search and carries the rest on `RepairOutcome.advisory`,
        which is where a polyA hexamer in the user's own LTR belongs: reported,
        never chased.

        SCOPED TO THE HARD_REPAIR RULES, which is a performance decision with a
        correctness argument under it. This callable runs once per candidate --
        up to `max_candidates` per iteration -- so every rule it touches is paid
        for hundreds of times per iteration. A SOFT rule can never return a
        repairable breach (the weighted sum is not in this loop), and E8's k-mer
        index or B1's fold evaluated here would be pure waste.

        That paragraph described the intent and not the code until now: this
        returned `self.findings(c).repairable`, and `findings()` walks every
        spec. Measured on a 500 aa protein assembled into the 3.1 kb synthetic
        lentiviral backbone against the post-#92/#93 catalog: 65 ms per repair
        candidate, against 7 ms for the eight HARD_REPAIR specs actually
        consumed -- a ~9x waste. Two rules were most of it,
        `f3_inverted_repeats` at ~29 ms and `e8_kmer_uniqueness` (the rule this
        docstring already named) at ~20 ms, and every result was discarded.

        The narrowing is behaviour-preserving, and the argument is short enough
        to check: `repair_specs()` selects the specs where
        `enforcement_for(slot)` is HARD_REPAIR for some active slot, and
        `_enforcement_of` returns only values drawn from that same set -- so a
        spec outside `repair_specs()` cannot produce a breach that reaches
        `repairable`.

        With ONE exception, handled explicitly below rather than inherited:
        `_enforcement_of` falls back to the `enforcement` ClassVar when there
        are no active slots, and `repair_specs()` is empty there because its
        `any()` is over nothing. `DesignContext` requires a slot but does not
        require one to be ENABLED, so a fully disabled context would silently
        change behaviour. It keeps the full walk.

        `findings()` is deliberately left alone -- `advise()` reads its
        `hard_check` and the design lane reads its `evaluations` for the
        scorecard, so both need every spec.
        """
        if not self.ctx.active_slots:
            # No enabled slot, so `repair_specs()` is empty while
            # `_enforcement_of` still reads each rule's ClassVar. Narrowing here
            # would drop breaches the unscoped path returns.
            return lambda c: self.findings(c).repairable

        # Hoisted out of the per-candidate loop: `gate()` and `enforcement_for`
        # read only the slot, so neither can change between repair iterations.
        specs = self.repair_specs()

        def find(c: Construct) -> tuple[Breach, ...]:
            out: list[Breach] = []
            for spec in specs:
                ev = spec.evaluate(c, self.ctx, self.svc)
                if ev.passes:
                    continue  # the rule's own verdict; a warn-band finding is not a failure
                out.extend(
                    b
                    for b in ev.breaches
                    if self._enforcement_of(spec, b) is Enforcement.HARD_REPAIR
                )
            return tuple(out)

        return find

    def advise(self) -> Callable[[Construct], tuple[Breach, ...]]:
        """The HARD_CHECK findings no codon can fix. Reported, never chased.

        Distinct from `RepairOutcome.advisory` (unfixable HARD_REPAIR, carried by
        the search): these never enter `repair()` at all. Surfaced for the report
        layer; there is no consumer in the solver.
        """
        return lambda c: self.findings(c).hard_check

    def policies(self, default_window: int = 50) -> dict[str, RulePolicy]:
        """A `RulePolicy` per HARD_REPAIR rule, for `repair()`.

        This is hazard 4 of #58, and `repair()`'s per-rule machinery is where it
        is solved: the four HARD_REPAIR rules declare three different
        localisations (E2 `WINDOW_MINUS_1`, E5 and F1 `PAIRED_SEGMENTS`, E7
        `WHOLE_SCOPE`), so one global `policy=` gives at least two of them the
        wrong repair window. Each rule's own `localization`, `repair` discipline
        and instance `window` travel in its policy.

        `repair` is read per rule and never globally escalated: `repair()` applies
        each rule's discipline to its own breaches (SINGLE_PASS retires a breach
        after one attempt, FIXED_POINT re-targets until the rule stops producing
        it), so a splice-removal rule that needs FIXED_POINT declares it on itself
        -- the escalation is per rule, where CLAUDE.md 3.6 puts it.
        """
        out: dict[str, RulePolicy] = {}
        for spec in self.repair_specs():
            window = getattr(spec, "window", default_window)
            out[spec.id] = RulePolicy(
                localization=spec.localization,
                repair=spec.repair,
                window=int(window) if isinstance(window, int) else default_window,
                motif_len=_DEFAULT_MOTIF_LEN,
                priority=0,  # pure round-robin; repair() already prevents starvation
            )
        return out

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

        Three ways this returns None, each a refusal to half-arm the oracle:

        - no E2 at all (disabled, excluded, or gated off in every active slot --
          it gates off for IVT mRNA), so there is no band to enforce;
        - E2 could not resolve numbers;
        - an ADAPTER-ON selection. E2 measures the fragment the vendor
          synthesises, adapters included; I7 measures the designable span alone,
          because the oracle has no vendor data. With adapters those are different
          bases and different numbers, which is #59's contradiction from a third
          side. `VendorSelection.of()` guarantees an adapter-on selection is
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
            # compile-time dependency on the vendor catalogue. Pinned by test.
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
        slots = self.ctx.active_slots
        if breach.slot_role is not None:
            for slot in slots:
                if slot.role == breach.slot_role:
                    return spec.enforcement_for(slot)
        if not slots:
            return spec.enforcement
        return max((spec.enforcement_for(s) for s in slots), key=lambda e: _SEVERITY[e])


def build_rule_set(
    ctx: DesignContext,
    svc: Services,
    *,
    vendors: VendorSelection = DEFAULT_SELECTION,
    overrides: Mapping[str, Mapping[str, object]] | None = None,
    include: Callable[[type[Spec]], bool] | None = None,
) -> RuleSet:
    """Discover, filter, calibrate, instantiate and gate the catalog.

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
    providers rather than importing them. Every provider stays injectable.

    `fold` stays None when ViennaRNA is absent -- never a stub. Every threshold
    in BT5 is a kcal/mol number, so a stub returning plausible energies would
    flow through the scorers, the null and the percentile unchallenged and come
    out the far end as a confident rank.

    THE THREE CASTS ARE REAL DEFECTS, not type-checker appeasement, and they are
    narrow on purpose so that fixing one removes exactly one. Nothing in `src/`
    had ever constructed a `Services` before this module, so no lane's concrete
    provider had ever been checked against the protocol it claims to implement,
    and all three fail: `ViennaFold.version` is a read-only property where the
    protocol wants a settable ClassVar (M6); `NcbiGeneticCode.table_id` is a
    read-only property where the protocol wants a plain attribute (M5); and
    `FileTableProvider.usage` returns a `CodonUsage`, not the `Mapping[str,float]`
    the protocol declares, `@cache`-wrapped on top (M5). All three are other
    lanes' to fix, and correcting the protocols would be a `core/` amendment.
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
    gc_window: int = 50,
    max_candidates: int = 256,
    left_flank: str = "",
    right_flank: str = "",
    seed: int = 0,
    table_id: int | None = None,
) -> OptimizeResult:
    """`optimize()` with every argument derived from ONE rule set.

    This is the call that makes the catalog real, and the reason it exists here
    rather than in `pipeline` is that all four derived arguments have to come from
    the same instantiated rules. Supply the motifs from one place and the
    validator's GC band from another and you get #59 again: E2 gating a fragment
    against the selected vendor's band while I7 gates the same span against a
    different one, refusing constructs the rules pass and passing constructs the
    rules refuse.

    Pass `original_backbone` -- `vector.assemble()` returns it as
    `Assembly.reference` -- to arm I9. Without it the invariant that proves BT5
    did not touch a single backbone base does not run.

    HARD_REPAIR advisories (a polyA in the user's LTR) come back on
    `result.repair_outcome.advisory`; HARD_CHECK findings (an over-length
    fragment) are read separately with `rules.advise()(result.construct)`, since
    `OptimizeResult` carries only the search's own outcome.
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
        policies=rules.policies(gc_window),
        max_candidates=max_candidates,
        original_backbone=original_backbone,
        left_flank=left_flank,
        right_flank=right_flank,
        seed=seed,
        table_id=table_id,
    )


__all__ = [
    "Findings",
    "OracleBounds",
    "RuleSet",
    "build_rule_set",
    "default_services",
    "optimize_with",
]
