"""E9 -- can this fragment be ordered at all.

Every other manufacturability rule asks whether the vendor will *succeed* at
building the sequence. This one asks whether they will accept the order, and it
is the only E-rule whose answer no codon can change.

**The minimum is the bound that surprises people.** Everyone remembers that a
gene fragment has a length ceiling. Almost nobody remembers there is a floor:
Twist Gene Fragments and IDT eBlocks both start at 300 bp, which is 100 codons,
so a protein shorter than about 100 residues cannot be ordered as a fragment
from either. Only gBlocks reach lower, to 125 bp, and below that no gene
fragment product exists at all. A short peptide, a tag-only construct, a small
domain -- exactly the things someone back-translates without a second thought --
land under the floor, and without this rule BT5 would design one, rank it
against a null, hash it, and write it into the order file. The user would find
out at checkout.

**Why HARD_CHECK rather than HARD_REPAIR.** Length is fixed by the protein and
the genetic code: every synonymous codon is three bases, so the mutation space
the solver searches does not contain a single sequence of a different length.
Routing this to the solver would send it hunting through a space where no
solution exists and report infeasibility on a design that is fine -- precisely
the failure `fixable_by_codon_choice` exists to prevent. The fix is a different
product, a different vendor, or splitting the order, and all three are decisions
for a person.

That also makes the ADVISORY content the rule's real output. A finding that says
only "too short" is useless when nothing in the sequence can respond to it, so
every breach names the configurations that WOULD accept this length. "125-3000
bp as a gBlock, nowhere else" is the sentence the user can act on.

**One unresolved question, reported rather than guessed.** These ranges are
published against "gene fragments" without saying whether the bound applies to
the ordered insert or to the insert plus its adapters. It only matters within 44
bp of a boundary, and only for an adapter-on order -- but inside that margin the
two readings disagree about whether the order is accepted. Rather than pick one,
the rule computes both and reports a distinct, non-blocking finding when they
disagree, so a design near a boundary is flagged as unresolved instead of being
silently passed or silently failed.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import ClassVar

from bt5.core.context import ContextSlot, DesignContext
from bt5.core.registry import register
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
)
from bt5.core.types import Construct
from bt5.rules.fragment import VENDOR_ADAPTERS, VENDOR_LENGTHS, fragments

#: A fragment that cannot be ordered under either reading of the bound.
UNORDERABLE = 1.0
#: A fragment whose acceptance depends on the insert-vs-total question above.
#: Below 1.0 deliberately: it is a real finding and not a proven rejection, so it
#: is reported without failing the rule.
AMBIGUOUS = 0.5

DEFAULT_VENDOR = "twist_gene_fragment"


def alternatives(length: int, exclude: str) -> tuple[str, ...]:
    """Every other configuration whose range contains `length`."""
    return tuple(
        sorted(
            name
            for name, (lo, hi) in VENDOR_LENGTHS.items()
            if name != exclude and lo <= length <= hi
        )
    )


@register
class LengthTiers:
    id: ClassVar[str] = "e9_length_tiers"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "Orderable length range for the chosen vendor configuration"
    enforcement: ClassVar[Enforcement] = Enforcement.HARD_CHECK
    evidence: ClassVar[Evidence] = Evidence.VENDOR_ASSERTED
    direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
    unit: ClassVar[str] = "fragments outside the orderable length range"
    band: ClassVar[tuple[float, float] | None] = None
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "Twist gene fragments run 300 bp to 5 kb; 300 bp is a floor, not a "
            "recommendation, and clonal genes rather than fragments are what reach 7 kb",
            "https://www.twistbioscience.com/faq/gene-synthesis/are-there-any-sequence-limitationsdesign-guidelines-genes-which-i-should-follow",
            2026,
            sign="supports",
        ),
        Citation(
            "IDT eBlocks Gene Fragments are 300-1500 bp",
            "https://www.idtdna.com/pages/products/genes-and-gene-fragments/double-stranded-dna-fragments/eblocks-gene-fragments",
            2026,
            sign="supports",
        ),
        Citation(
            "IDT synthesises gBlocks Gene Fragments between 125 bp and 3 kb -- the only "
            "configuration here that accepts an insert under 300 bp",
            "https://www.idtdna.com/pages/support/faqs/what-is-the-length-of-gblocks-gene-fragments-that-idt-can-synthesize",
            2026,
            sign="supports",
        ),
    )
    last_verified: ClassVar[str] = "2026-08-28"
    weight_provenance: ClassVar[str] = ""  # hard rule; never weighted
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.0
    #: Zero, and not because nobody got round to tuning it. Steering biases the
    #: DP toward some codons over others, and every codon is three bases: there
    #: is no synonymous choice that moves a length. A non-zero value here would
    #: be a term that cannot affect its own objective.
    steering_weight: ClassVar[float] = 0.0
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.WHOLE_SCOPE
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "cheap"
    conflicts_with: ClassVar[tuple[str, ...]] = ()
    brief_ref: ClassVar[str] = "2.E9"
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "vendor": {
                "type": "string",
                "default": DEFAULT_VENDOR,
                "enum": sorted(VENDOR_LENGTHS),
                "description": (
                    "Which vendor configuration the fragment is ordered as. Each "
                    "has its own orderable length range, and the minimum differs "
                    "between them by more than the maximum does."
                ),
            },
        },
    }

    def __init__(self, vendor: str = DEFAULT_VENDOR) -> None:
        if vendor not in VENDOR_LENGTHS:
            raise ValueError(
                f"unknown vendor {vendor!r}; have {sorted(VENDOR_LENGTHS)}. "
                f"'none' is not orderable from anyone and has no length range."
            )
        self.vendor = vendor

    def gate(self, slot: ContextSlot) -> bool:
        # Every construct BT5 designs is ordered as DNA, in every context.
        return True

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """None. The automaton decides from a bounded suffix of bases; this rule
        is a property of the whole molecule's length and of nothing else."""
        return None

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        adapters = VENDOR_ADAPTERS[self.vendor]
        lo, hi = VENDOR_LENGTHS[self.vendor]
        breaches: list[Breach] = []
        frags = fragments(c, adapters)

        for frag in frags:
            length = frag.ordered.length
            total = length + adapters.total
            ordered_ok = lo <= length <= hi
            total_ok = lo <= total <= hi
            if ordered_ok and total_ok:
                continue

            alts = alternatives(length, self.vendor)
            where = (
                "orderable instead as " + ", ".join(alts)
                if alts
                else "no configured vendor accepts this length"
            )
            if ordered_ok != total_ok:
                message = (
                    f"{length} bp of ordered DNA plus {adapters.total} bp of "
                    f"{self.vendor} adapters is {total} bp, and the published "
                    f"{lo}-{hi} bp range does not say whether it applies to the "
                    f"insert or to the total. One reading accepts this order and "
                    f"the other rejects it; confirm with the vendor before building"
                )
                magnitude = AMBIGUOUS
            else:
                side = "below the" if length < lo else "above the"
                message = (
                    f"{length} bp is {side} {lo}-{hi} bp {self.vendor} range, so "
                    f"this fragment cannot be ordered as designed. No codon choice "
                    f"changes a length -- {where}"
                )
                magnitude = UNORDERABLE

            breaches.append(
                Breach(
                    spec_id=self.id,
                    interval=frag.origin,
                    magnitude=magnitude,
                    message=message,
                    fixable_by_codon_choice=False,
                    detail={
                        "vendor": self.vendor,
                        "ordered_bp": float(length),
                        "with_adapters_bp": float(total),
                        "min_bp": float(lo),
                        "max_bp": float(hi),
                        "alternatives": ", ".join(alts) if alts else "none",
                    },
                )
            )

        return Evaluation(
            spec_id=self.id,
            passes=not any(b.magnitude >= UNORDERABLE for b in breaches),
            raw_score=float(sum(b.magnitude for b in breaches)),
            breaches=tuple(breaches),
            n_evaluated=len(frags),
        )
