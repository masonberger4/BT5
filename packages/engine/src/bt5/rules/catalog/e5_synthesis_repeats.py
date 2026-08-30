"""E5 -- repeats in the fragment the vendor actually builds.

Not F1 with different numbers. F1 asks whether the assembled PLASMID will delete
a chunk in E. coli; this asks whether the vendor's assembly PCR of the ORDERED
FRAGMENT will mis-prime. Different molecule, different mechanism, different
scope -- see `bt5.rules.fragment` for why the scope is the load-bearing part.

Two consequences a reader should expect and not mistake for a bug:

  * A repeat spanning the origin is an F1 finding and NOT an E5 finding. The
    ordered fragment is linear; its two ends are never in the same tube.
  * A repeat between the insert and the backbone is an F1 finding and NOT an E5
    finding. The vendor never receives the backbone.

**Why Tm and not just length.** The published vendor thresholds are lengths --
Twist and GenScript both name 20 bp -- but the physics is duplex stability at
the annealing step, and length is a poor proxy for it. Measured with the
parameters pinned below:

    GCGCGCGCGCGCGC          14 bp   65.2 C   mis-primes; under every length rule
    ACGTACGTACGTACGTACGT    20 bp   53.1 C   flagged by length; stable at neither

So a pure length rule flags the harmless one and misses the dangerous one. E5
carries both criteria and treats either as sufficient.

**Why the floor is 12 bp.** Not a convention. A pure-GC 12-mer melts at 58.7 C
and a pure-GC 10-mer at 49.7 C, so 12 bp is the shortest repeat that can reach a
60 C anneal AT ALL, whatever its composition. Scanning below it cannot find a
mis-priming repeat, and would only add chance hits: an 8-mer has roughly 17
chance occurrences per 1.5 kb. The brief's "warn >12-16 bp" and the physics
agree, which is the only reason the number is worth pinning.

The brief's "flag repeat clusters at 8-9 bp" clause is deliberately NOT here.
That is a density statement, not a pair statement, and E6 measures it as one --
implementing it here as well would report the same repetitiveness twice.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import ClassVar

from Bio.SeqUtils import MeltingTemp

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
from bt5.core.types import Construct, Interval
from bt5.rules.fragment import Fragment, fragments
from bt5.rules.vendors import (
    DEFAULT_SELECTION,
    DEFAULT_VENDOR,
    VendorSelection,
    all_keys,
    require_selection,
)

#: The shortest repeat that can reach a 60 C duplex at any composition
#: (pure-GC 12-mer = 58.7 C, pure-GC 10-mer = 49.7 C). Also the brief's warn floor.
MIN_LENGTH_BP = 12
#: Twist and GenScript both publish 20 bp as the point of rejection.
HARD_LENGTH_BP = 20
#: Above this a repeat is not a complexity surcharge, it is an unbuildable order.
SEVERE_LENGTH_BP = 200
#: Assembly PCR anneal. A repeat whose duplex is stable here can mis-prime.
ANNEAL_C = 60.0
#: A tandem array must not flood the panel.
MAX_FINDINGS = 200

#: The ruler. A melting temperature is meaningless without the conditions it was
#: computed under, exactly as a dG is meaningless without its energy parameters,
#: so the conditions are pinned here, reported in `engine_calibration`, and
#: travel with every finding rather than being inherited from a library default
#: that can move underneath us.
TM_NN_TABLE = MeltingTemp.DNA_NN3  # Allawi & SantaLucia 1997
TM_NA_MM = 50.0
TM_OLIGO_NM = 25.0
TM_CONDITIONS = "Tm_NN / DNA_NN3 (Allawi & SantaLucia 1997), Na+ 50 mM, 25 nM each strand"


def duplex_tm(seq: str) -> float:
    """Nearest-neighbour duplex Tm under the pinned conditions above."""
    return float(
        # Biopython ships no stubs for MeltingTemp; the float() is the real
        # boundary, and it is what keeps the untyped return out of the rule.
        MeltingTemp.Tm_NN(  # type: ignore[no-untyped-call]
            seq,
            nn_table=TM_NN_TABLE,
            Na=TM_NA_MM,
            dnac1=TM_OLIGO_NM,
            dnac2=TM_OLIGO_NM,
        )
    )


def severity(length: int, tm: float, *, hard_len: int, anneal_c: float) -> str:
    """Three bands. Either criterion alone is sufficient for `hard`.

    Length because that is what the vendors publish and enforce at order time;
    Tm because that is the actual mechanism and it catches the short GC-rich
    repeats the published length rules miss.
    """
    if length >= SEVERE_LENGTH_BP:
        return "severe"
    if length >= hard_len or tm >= anneal_c:
        return "hard"
    return "warn"


@register
class SynthesisRepeats:
    id: ClassVar[str] = "e5_synthesis_repeats"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "Direct repeats in the synthesized fragment, by length and duplex Tm"
    enforcement: ClassVar[Enforcement] = Enforcement.HARD_REPAIR
    evidence: ClassVar[Evidence] = Evidence.VENDOR_ASSERTED
    direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
    unit: ClassVar[str] = "weighted repeat findings in the ordered fragment"
    band: ClassVar[tuple[float, float] | None] = None
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "Twist gene synthesis FAQ: repeated sequences are a named complexity "
            "trigger, with 20 bp the published limit. Adapters are an OPTION and "
            "not the default -- Twist states that adapter sequences are not added "
            "by default to Gene Fragments, and adapter-on is chosen at checkout",
            "https://www.twistbioscience.com/faq/gene-synthesis/are-adapter-sequences-appended-ends-my-sequences",
            2026,
            sign="supports",
        ),
        Citation(
            "Repetitive 9-mers per 100 bp and longest repetitive sequence are the two "
            "highest-importance features of the only synthesis-success model trained "
            "on real vendor outcomes (random forest, 1,076 orders, F1 0.928)",
            "https://pubs.acs.org/doi/10.1021/acssynbio.9b00460",
            2020,
            sign="supports",
        ),
        Citation(
            "Repeats and high GC are the documented cause of Gibson misassembly, and "
            "no two junctions may share a homology sequence -- the same duplex "
            "stability argument as assembly PCR mis-priming",
            "https://blog.addgene.org/plasmids-101-gibson-assembly",
            2023,
            sign="supports",
        ),
        Citation(
            "Nearest-neighbour thermodynamics used for the duplex Tm criterion; the "
            "parameter set is pinned so a Tm never travels without its conditions",
            "https://pubmed.ncbi.nlm.nih.gov/9236978/",
            1997,
            sign="supports",
        ),
    )
    last_verified: ClassVar[str] = "2026-08-28"
    weight_provenance: ClassVar[str] = ""  # hard rule; never weighted
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.0
    #: Below F1's 1.0. F1 steers on the whole plasmid and therefore already
    #: suppresses most of what E5 would; the increment here is the short GC-rich
    #: repeat and the adapter collision, which F1 cannot see at all.
    steering_weight: ClassVar[float] = 0.7
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.PAIRED_SEGMENTS
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "moderate"
    conflicts_with: ClassVar[tuple[str, ...]] = ()
    brief_ref: ClassVar[str] = "2.E5"
    #: None, and NOT the Tm parameter set, however tempting. This field is
    #: compared against the FOLD engine's calibration key by
    #: `check_engine_calibration`, so declaring 'biopython_tm_nn:...' here would
    #: raise CalibrationMismatchError against ViennaRNA on every run -- it names
    #: the engine a rule's THRESHOLDS are measured on, and this rule's are in bp
    #: and degrees C, not kcal/mol. The Tm conditions travel with each finding
    #: in `detail` instead, which is where the ruler actually needs to be.
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "min_len": {"type": "integer", "default": MIN_LENGTH_BP, "minimum": MIN_LENGTH_BP},
            "hard_len": {"type": "integer", "default": HARD_LENGTH_BP},
            "anneal_c": {"type": "number", "default": ANNEAL_C},
            "vendors": {
                "type": "array",
                "minItems": 1,
                "uniqueItems": True,
                "items": {"type": "string", "enum": list(all_keys())},
                "default": [DEFAULT_VENDOR],
                "description": (
                    "Which vendor product the fragment is ordered as. Only the "
                    "adapter-on options carry adapters; a plain Gene Fragment "
                    "order is the ordered DNA and nothing else."
                ),
            },
        },
    }

    def __init__(
        self,
        min_len: int = MIN_LENGTH_BP,
        hard_len: int = HARD_LENGTH_BP,
        anneal_c: float = ANNEAL_C,
        vendors: VendorSelection = DEFAULT_SELECTION,
    ) -> None:
        if min_len < MIN_LENGTH_BP:
            raise ValueError(
                f"min_len {min_len} is below the shortest repeat that can reach a "
                f"60 C duplex at any composition (a pure-GC 12-mer is 58.7 C), so it "
                f"can only add chance hits"
            )
        if hard_len < min_len:
            raise ValueError(f"hard_len {hard_len} must not be below min_len {min_len}")
        if not 30.0 <= anneal_c <= 80.0:
            raise ValueError(f"anneal_c {anneal_c} is outside any real PCR protocol")
        self.min_len = min_len
        self.hard_len = hard_len
        self.anneal_c = anneal_c
        self.vendors = require_selection(vendors)

    def gate(self, slot: ContextSlot) -> bool:
        # Every construct BT5 designs is ordered as DNA, including an IVT mRNA
        # template, so this applies in every context.
        return True

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """None. Whether a codon completes a repeat depends on the whole
        fragment, not on a bounded suffix."""
        return None

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        adapters = self.vendors.adapters
        breaches: list[Breach] = []
        worst = "warn"
        scanned = 0

        for frag in fragments(c, adapters):
            scanned += frag.ordered.length
            for breach, band in self._scan(frag, svc):
                if len(breaches) >= MAX_FINDINGS:
                    break
                breaches.append(breach)
                if band == "severe" or (band == "hard" and worst != "severe"):
                    worst = band

        return Evaluation(
            spec_id=self.id,
            passes=worst == "warn",
            raw_score=float(sum(b.magnitude for b in breaches)),
            breaches=tuple(breaches),
            n_evaluated=scanned,
        )

    def _scan(self, frag: Fragment, svc: Services) -> list[tuple[Breach, str]]:
        index = svc.kmer.of(frag.construct, self.min_len)
        out: list[tuple[Breach, str]] = []

        for first, second in index.duplicates(self.min_len):
            span = Interval(first.start, max(first.end, second.end))
            mapped = frag.to_construct(span)
            if mapped is None:
                continue  # the vendor's adapters colliding with themselves
            unit = frag.sequence[first.start : first.end]
            tm = duplex_tm(unit)
            length = first.length
            band = severity(length, tm, hard_len=self.hard_len, anneal_c=self.anneal_c)
            adapter = frag.touches_adapter(span)

            out.append(
                (
                    Breach(
                        spec_id=self.id,
                        interval=mapped,
                        magnitude={"warn": 0.5, "hard": 2.0, "severe": 4.0}[band],
                        message=(
                            f"{length} bp repeat in the ordered fragment at "
                            f"{first.start} and {second.start} (fragment coordinates), "
                            f"duplex Tm {tm:.1f} C. "
                            + (
                                f"Stable at the {self.anneal_c:.0f} C anneal, so it can "
                                f"mis-prime during assembly PCR"
                                if tm >= self.anneal_c
                                else f"Below the {self.anneal_c:.0f} C anneal"
                            )
                            + (
                                f"; at or over the {self.hard_len} bp vendor limit"
                                if length >= self.hard_len
                                else ""
                            )
                            + (
                                ". One copy is vendor adapter sequence, which BT5 cannot "
                                "change -- recode the insert side"
                                if adapter
                                else ""
                            )
                        ),
                        # The fragment IS the designable region, so anything that
                        # maps back is recodable. `overlaps_editable` on the
                        # parent is still the honest test: it is what the rest of
                        # the catalog uses and it keeps the answer consistent if
                        # the parent's segmentation ever disagrees.
                        fixable_by_codon_choice=frag.construct.overlaps_editable(span),
                        detail={
                            "length": float(length),
                            "duplex_tm_c": round(tm, 1),
                            "severity": band,
                            "tm_conditions": TM_CONDITIONS,
                            "vendor": frag.vendor,
                            "involves_adapter": "yes" if adapter else "no",
                        },
                    ),
                    band,
                )
            )
        return out
