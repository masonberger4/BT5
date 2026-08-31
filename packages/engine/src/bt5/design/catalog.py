"""Run the whole rule catalog against one design's context, and hand the solver
exactly the rules it should chase -- no more.

This is where the enforcement enum stops being metadata and starts steering the
search. Three facts drive everything here:

- **Enforcement is per slot, not per class.** `d4_internal_polya` is SOFT by its
  ClassVar and HARD_REPAIR in a lentiviral slot. So selection reads
  `enforcement_for(slot)`, never `.enforcement` -- the opposite of what the
  preset resolver does, and getting it wrong puts a titer-killing polyA in the
  weighted sum instead of the repair set.
- **The slot loop is the rule's job.** A rule's `evaluate` already iterates
  `ctx.active_slots` internally, so the finder calls each rule ONCE. Looping over
  slots here would double-count every breach.
- **The backbone already carries forbidden motifs.** A vector's own 5'UTR can
  hold an XbaI site; `LatticeTerms` would forbid it, and then Tier A refuses to
  place a flank that contains it, or the validator refuses the whole construct.
  So the lattice forbidden set is partitioned: motifs the backbone carries become
  advisories and are excluded from both the automaton and the validator, and only
  the motifs a codon choice could actually create are enforced.
"""

from __future__ import annotations

import inspect
from collections.abc import Sequence
from dataclasses import dataclass

from bt5.core.context import DesignContext
from bt5.core.services import Services
from bt5.core.spec import Breach, Enforcement, Spec
from bt5.core.types import DNA_ALPHABET, Construct, reverse_complement
from bt5.design.errors import DesignError
from bt5.rules.vendors import VendorSelection
from bt5.solver.repair import BreachFinder, RulePolicy
from bt5.vector.backbone import InsertionSite, VectorBackbone

#: motif_len used when localizing a motif rule's breach. The polyA hexamers d4
#: works are 6 nt; the other HARD_REPAIR rule (e2) localizes by window, not motif,
#: so this value is only read for MOTIF_LEN_MINUS_1 rules and 6 is correct there.
_MOTIF_LEN = 6


@dataclass(frozen=True)
class Catalog:
    """The rules of one design, sorted into what each stage of the pipeline sees."""

    rules: tuple[Spec, ...]  # every rule gated on at least one active slot
    hard_repair: tuple[Spec, ...]  # the solver's BreachFinder chases these
    hard_check: tuple[Spec, ...]  # reported and blocking, never chased (advisories)
    scored: tuple[Spec, ...]  # SOFT and not hard for ANY active slot (the 3.5 guard)
    usable_forbidden: tuple[str, ...]  # lattice motifs a codon choice could create
    carried_forbidden: tuple[str, ...]  # motifs the backbone already carries
    policies: dict[str, RulePolicy]  # per-rule dispatch for the HARD_REPAIR set


def _instantiate(cls: type[Spec], vendors: VendorSelection) -> Spec:
    """Build a rule, threading the design's vendor selection into the rules that
    take one so E2's band and the solver's gc_bounds come from the SAME
    selection. Every catalog rule has an all-defaulted constructor."""
    if "vendors" in inspect.signature(cls).parameters:
        return cls(vendors=vendors)  # type: ignore[call-arg]
    return cls()


def _gated_slots(rule: Spec, ctx: DesignContext) -> list[object]:
    return [slot for slot in ctx.active_slots if rule.gate(slot)]


def _immutable_backbone(backbone: VectorBackbone, site: InsertionSite) -> str:
    """The backbone with the insertion span masked out.

    A forbidden motif inside the span that gets REPLACED is the solver's to
    remove; only a motif outside it is one the backbone truly carries. Masking
    with a non-ACGT character means those replaced bases can never match a motif,
    so the scan sees exactly the immutable sequence.
    """
    seq = list(backbone.sequence)
    n = len(seq)
    for p in range(site.interval.start, site.interval.end):
        seq[p % n] = "N"
    return "".join(seq)


def _carries(motif: str, immutable: str, *, circular: bool) -> bool:
    """Does the immutable backbone contain this motif or its reverse complement?

    Both strands, because the solver closes the forbidden set under reverse
    complement -- a motif the backbone carries on either strand is one no codon
    choice can remove. Circular-aware: a motif can straddle the origin.
    """
    hay = immutable
    if circular and len(motif) > 1:
        hay = immutable + immutable[: len(motif) - 1]
    return motif in hay or reverse_complement(motif) in hay


def _lattice_forbidden(rules: Sequence[Spec], ctx: DesignContext) -> tuple[str, ...]:
    """The union of every active rule's forbidden motifs, deduplicated and sorted.

    Every motif is asserted ACGT-only: `LatticeTerms.forbidden` is documented
    IUPAC, but the automaton `KeyError`s on a degenerate base, so a rule that ever
    emits one must be caught here rather than deep in Tier A. (Filed as an M1 fix;
    today only d1 and e1 populate the field and both are pure ACGT.)
    """
    motifs: set[str] = set()
    for rule in rules:
        terms = rule.lattice_terms(ctx)
        if terms is None:
            continue
        for motif in terms.forbidden:
            bad = sorted(set(motif.upper()) - DNA_ALPHABET)
            if bad:
                raise DesignError(
                    f"rule {rule.id} forbids {motif!r}, which contains non-ACGT "
                    f"characters {bad}; the automaton cannot consume IUPAC codes"
                )
            motifs.add(motif.upper())
    return tuple(sorted(motifs))


def build_catalog(
    ctx: DesignContext,
    *,
    backbone: VectorBackbone,
    site: InsertionSite,
    vendors: VendorSelection,
    default_window: int = 50,
) -> Catalog:
    """Instantiate the catalog, gate it to `ctx`, and sort it for the pipeline."""
    from bt5.core.registry import all_specs, discover

    discover()
    instantiated = [_instantiate(cls, vendors) for cls in all_specs()]
    rules = tuple(r for r in instantiated if _gated_slots(r, ctx))

    hard_repair: list[Spec] = []
    hard_check: list[Spec] = []
    scored: list[Spec] = []
    for rule in rules:
        slots = _gated_slots(rule, ctx)
        resolved = [rule.enforcement_for(slot) for slot in slots]  # type: ignore[arg-type]
        if any(e is Enforcement.HARD_REPAIR for e in resolved):
            hard_repair.append(rule)
        if any(e is Enforcement.HARD_CHECK for e in resolved):
            hard_check.append(rule)
        # The 3.5 guard the preset resolver misses: a rule that is SOFT by its
        # ClassVar but hard in ANY active slot must NOT enter the weighted sum.
        if rule.enforcement.is_scored and not any(e.is_hard for e in resolved):
            scored.append(rule)

    all_forbidden = _lattice_forbidden(rules, ctx)
    immutable = _immutable_backbone(backbone, site)
    usable: list[str] = []
    carried: list[str] = []
    for motif in all_forbidden:
        if _carries(motif, immutable, circular=backbone.is_circular):
            carried.append(motif)
        else:
            usable.append(motif)

    policies = {
        rule.id: RulePolicy(
            localization=rule.localization,
            repair=rule.repair,
            window=int(getattr(rule, "window", default_window)),
            motif_len=_MOTIF_LEN,
            priority=0,
        )
        for rule in hard_repair
    }

    return Catalog(
        rules=rules,
        hard_repair=tuple(hard_repair),
        hard_check=tuple(hard_check),
        scored=tuple(scored),
        usable_forbidden=tuple(usable),
        carried_forbidden=tuple(carried),
        policies=policies,
    )


def breach_finder(
    hard_repair: Sequence[Spec], ctx: DesignContext, services: Services
) -> BreachFinder:
    """A BreachFinder over exactly the HARD_REPAIR rules.

    Each rule is evaluated ONCE -- its own `evaluate` loops the gated slots -- so
    a breach is attributed to its slot without being counted per slot here.
    HARD_CHECK rules are deliberately absent: no codon choice moves their
    findings, so they belong in the advisories, not the search.
    """

    def find(construct: Construct) -> tuple[Breach, ...]:
        out: list[Breach] = []
        for rule in hard_repair:
            out.extend(rule.evaluate(construct, ctx, services).breaches)
        return tuple(out)

    return find
