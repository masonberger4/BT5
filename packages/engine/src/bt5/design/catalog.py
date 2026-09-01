"""The design lane's slice of catalog work -- what `bt5.solver.catalog` cannot know.

`bt5.solver.catalog.build_rule_set` runs the whole catalog against a context and
hands the pipeline a forbidden set, a breach finder, per-rule policies and the
oracle's GC band, all from one `VendorSelection` so they cannot disagree. This
module does NOT re-implement any of that. Two things the solver cannot do,
because it has no backbone, live here:

- The forbidden set it derives can include a motif the USER'S OWN BACKBONE
  already carries -- an XbaI site in a 5'UTR. Handed to Tier A that makes it
  refuse a flank it cannot recode; handed to the validator it makes I6 refuse the
  whole construct. Only the design lane, which holds the backbone and the
  insertion site, can partition those out.
- The scored objectives -- SOFT, and not hard in ANY active slot (the §3.5
  guard) -- are the report's to name `unavailable`. The solver leaves weights to
  M3 and never assembles that set.
"""

from __future__ import annotations

from collections.abc import Sequence

from bt5.core.spec import Spec
from bt5.core.types import DNA_ALPHABET, reverse_complement
from bt5.design.errors import DesignError
from bt5.solver.catalog import RuleSet
from bt5.vector.backbone import InsertionSite, VectorBackbone


def _immutable_backbone(backbone: VectorBackbone, site: InsertionSite) -> str:
    """The backbone with the insertion span masked out.

    A forbidden motif inside the span that gets REPLACED is the solver's to
    remove; only a motif outside it is one the backbone truly carries. Masking
    with a non-ACGT character means those replaced bases can never match a motif.
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


def partition_forbidden(
    forbidden: Sequence[str], backbone: VectorBackbone, site: InsertionSite
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Split the solver's forbidden set into usable vs. backbone-carried.

    `usable` motifs a codon choice could actually create -- fed to the automaton
    and the validator. `carried` motifs the backbone already holds outside the
    insert -- excluded from both, and surfaced as advisories, because no codon can
    remove them.

    Every motif is asserted ACGT-only: `LatticeTerms.forbidden` is documented
    IUPAC but the automaton `KeyError`s on a degenerate base (filed as an M1 fix,
    #73); today only d1 and e1 populate it and both are pure ACGT.
    """
    immutable = _immutable_backbone(backbone, site)
    usable: list[str] = []
    carried: list[str] = []
    for motif in forbidden:
        bad = sorted(set(motif.upper()) - DNA_ALPHABET)
        if bad:
            raise DesignError(
                f"forbidden motif {motif!r} contains non-ACGT characters {bad}; the "
                f"automaton cannot consume IUPAC codes"
            )
        if _carries(motif, immutable, circular=backbone.is_circular):
            carried.append(motif)
        else:
            usable.append(motif)
    return tuple(usable), tuple(carried)


def scored_objectives(rule_set: RuleSet) -> tuple[Spec, ...]:
    """The SOFT rules that stay soft in every active slot -- the §3.5 guard.

    A rule SOFT by its ClassVar but hard in any active slot (d4 under lentiviral)
    must NOT enter the weighted sum; the preset resolver misses this because it
    reads the ClassVar (filed as #72). These become `unavailable` objectives in
    the scorecard, named rather than dropped.
    """
    slots = rule_set.ctx.active_slots
    return tuple(
        spec
        for spec in rule_set.specs
        if spec.enforcement.is_scored
        and not any(spec.enforcement_for(slot).is_hard for slot in slots if spec.gate(slot))
    )
