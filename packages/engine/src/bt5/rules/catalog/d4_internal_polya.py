"""D4 -- internal polyadenylation signals, on the strand that actually gets packaged.

This is the rule where getting the strand wrong is worst, and where it is
easiest. An internal polyA signal in a lentiviral transfer vector raised
expression 3-6.5x and cut FUNCTIONAL TITER 8-9x with CMV or EF1a promoters: the
genome is truncated before packaging, so the thing you measure goes up while the
thing you need goes down. A scan that reads the wrong strand returns clean on
exactly that construct.

So this rule is NOT a lattice rule, and that is deliberate. `LatticeTerms.forbidden`
is closed under reverse complement by the solver, which is right for a
restriction site and wrong here: AATAAA matters on the packaged strand and its
reverse complement TTTATT does not. Forbidding both would refuse designs for a
signal that cannot fire. So it reads `strand_for(ctx, slot)` per slot and
repairs, rather than being guaranteed by construction.

Enforcement is per-modality, which is what `enforcement_for` exists for. The
class attribute is the FLOOR: soft everywhere, escalated to HARD_REPAIR on the
lentiviral and genome-integrated modalities where the titer cost is measured.
A weighted sum is per-slot, so in a slot where this rule is hard it must be
excluded from the sum rather than weighted -- the same requirement `pipeline`
asserts for any hard rule.

The downstream element is what separates a hexamer that fires from one that does
not. A canonical hexamer alone is common; a hexamer WITH a GU-rich downstream
element 10-40 nt 3' of it is a functioning cleavage site, so that pairing
escalates rather than adding a second independent finding.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import ClassVar

from bt5.core.context import ContextSlot, DesignContext, Modality
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
    strand_for,
)
from bt5.core.types import Construct, Interval, Strand, reverse_complement

#: The two dominant hexamers: ~61.6% and ~15% of mammalian poly(A) sites.
CANONICAL: tuple[str, ...] = ("AATAAA", "ATTAAA")

#: Documented variants, individually much weaker. Reported, never hard.
VARIANTS: tuple[str, ...] = (
    "AGTAAA",
    "TATAAA",
    "CATAAA",
    "GATAAA",
    "AATATA",
    "AATACA",
    "AATAGA",
    "ACTAAA",
    "AAGAAA",
    "AATGAA",
)

#: Where the downstream element sits relative to the END of the hexamer.
DSE_FROM = 10
DSE_TO = 40

#: Modalities where an internal polyA is a measured titer loss rather than a
#: preference: the genome is truncated before packaging, or before integration.
HARD_MODALITIES: frozenset[Modality] = frozenset(
    {Modality.LENTIVIRAL, Modality.AAV, Modality.GENOME_INTEGRATED}
)


def _has_downstream_element(window: str) -> bool:
    """A GU-rich element: GTGT/TGTG, or a U-run of 4+ T within any 6 nt."""
    if "GTGT" in window or "TGTG" in window:
        return True
    return any(window[i : i + 6].count("T") >= 4 for i in range(max(0, len(window) - 5)))


@register
class InternalPolyA:
    id: ClassVar[str] = "d4_internal_polya"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "Internal polyadenylation signals on the packaged strand"
    #: The FLOOR. `enforcement_for` escalates on the packaged modalities.
    enforcement: ClassVar[Enforcement] = Enforcement.SOFT
    evidence: ClassVar[Evidence] = Evidence.EVIDENCE_BACKED
    direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
    unit: ClassVar[str] = "weighted signals"
    band: ClassVar[tuple[float, float] | None] = None
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "Internal polyA in a lentiviral transfer vector raised expression 3-6.5x but "
            "cut FUNCTIONAL TITER 8-9x with CMV or EF1a (and not at all with beta-actin "
            "or PSA/Pb) -- the measured quantity moves opposite to the needed one",
            "https://pubmed.ncbi.nlm.nih.gov/18627247/",
            2008,
            sign="supports",
        ),
        Citation(
            "AATAAA and ATTAAA account for ~61.6% and ~15% of dominant mammalian poly(A) "
            "hexamers; the remaining variants are individually far weaker",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC7100133/",
            2020,
            sign="supports",
        ),
    )
    last_verified: ClassVar[str] = "2026-08-28"
    weight_provenance: ClassVar[str] = (
        "High for a soft rule, because the effect size is measured rather than "
        "inferred and it is directional in the expensive way: 8-9x functional titer "
        "loss, on a construct whose expression assay looks BETTER. It is not weighted "
        "higher still because the promoter dependence is real -- the same study saw no "
        "effect with beta-actin or PSA/Pb -- so a hexamer is a liability, not a defect. "
        "In the packaged modalities enforcement_for escalates it out of the weighted "
        "sum entirely, which is where the 8-9x actually bites."
    )
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.7
    steering_weight: ClassVar[float] = 0.3
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.MOTIF_LEN_MINUS_1
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "cheap"
    conflicts_with: ClassVar[tuple[str, ...]] = ()
    brief_ref: ClassVar[str] = "2.D4"
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "include_variants": {
                "type": "boolean",
                "default": True,
                "description": "Also report the ten weak variant hexamers.",
            }
        },
    }

    def __init__(self, include_variants: bool = True) -> None:
        self.include_variants = include_variants

    def gate(self, slot: ContextSlot) -> bool:
        # Bacterial hosts do not 3'-process this way, and an IVT mRNA is
        # transcribed from a linear template with a defined end.
        return slot.modality not in (
            Modality.BACTERIAL_EXPRESSION,
            Modality.IVT_MRNA,
        )

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        """Hard where a truncated genome costs titer, soft where it costs a little
        expression. One frozen level would be wrong in every job."""
        if slot.modality in HARD_MODALITIES:
            return Enforcement.HARD_REPAIR
        return Enforcement.SOFT

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """Deliberately None. `forbidden` is closed under reverse complement by
        the solver, which would also forbid TTTATT -- a sequence that cannot
        polyadenylate anything -- and refuse designs for a signal on the strand
        nobody packages."""
        return None

    def _scan(self, c: Construct, strand: Strand, role: str) -> list[Breach]:
        """Every hexamer on ONE strand, in construct coordinates."""
        n = c.length
        forward = c.sequence
        # Read the requested strand, but report positions in construct
        # coordinates: a breach the user cannot locate on their own map is not
        # actionable, whichever strand found it.
        seq = forward if strand == 1 else reverse_complement(forward)
        scan = seq + seq[:DSE_TO] if c.is_circular else seq

        hexamers = (*CANONICAL, *VARIANTS) if self.include_variants else CANONICAL
        breaches: list[Breach] = []
        for hexamer in hexamers:
            canonical = hexamer in CANONICAL
            start = scan.find(hexamer)
            while start != -1:
                if start < n:
                    tail = scan[start + 6 + DSE_FROM : start + 6 + DSE_TO]
                    escalated = canonical and _has_downstream_element(tail)
                    lo = start if strand == 1 else n - start - 6
                    iv = Interval(max(0, lo), max(0, lo) + 6, strand)
                    breaches.append(
                        Breach(
                            spec_id=self.id,
                            interval=iv,
                            magnitude=3.0 if escalated else (1.0 if canonical else 0.2),
                            message=(
                                f"polyA hexamer {hexamer!r} on the packaged strand "
                                f"({'+' if strand == 1 else '-'}) at {iv.start}"
                                + (
                                    "; a GU-rich downstream element sits 10-40 nt 3' of it, "
                                    "so this is a functioning cleavage site rather than a "
                                    "hexamer that happens to be present"
                                    if escalated
                                    else ""
                                )
                            ),
                            slot_role=role,
                            fixable_by_codon_choice=c.overlaps_editable(iv),
                            detail={
                                "hexamer": hexamer,
                                "canonical": "yes" if canonical else "no",
                                "downstream_element": "yes" if escalated else "no",
                                "strand": float(strand),
                            },
                        )
                    )
                start = scan.find(hexamer, start + 1)
        return breaches

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        """One pass per gated slot, on that slot's own packaged strand.

        Two slots can disagree about which strand matters -- a reverse-oriented
        cassette packaged forward is the case -- so the breaches carry
        `slot_role` and the conflict detector sees them as what they are.
        """
        breaches: list[Breach] = []
        seen: set[tuple[str, int, int]] = set()
        for slot in ctx.active_slots:
            if not self.gate(slot):
                continue
            strand = strand_for(ctx, slot)
            for breach in self._scan(c, strand, slot.role):
                key = (breach.spec_id, breach.interval.start, int(strand))
                if key in seen:
                    continue
                seen.add(key)
                breaches.append(breach)

        return Evaluation(
            spec_id=self.id,
            passes=not any(b.magnitude >= 1.0 for b in breaches),
            raw_score=float(sum(b.magnitude for b in breaches)),
            breaches=tuple(breaches),
            n_evaluated=c.length,
        )
