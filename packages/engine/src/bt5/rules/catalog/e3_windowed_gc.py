"""E3 -- windowed GC at 50 bp and 100 bp, with enforcement set by measurement.

`brief.md:140`, grade **A**: *"Windowed GC 50 bp: hard-fail any window <10% or >90%
(Twist High-Complexity trigger); warn <25% or >75%. Windowed GC 100 bp: warn outside
25-65% (GenScript GenTitan -- the only vendor publishing a windowed rule)."*

**The 50 bp hard-fail is not encoded as a gate, and that is the whole point of this
rule.** `docs/design/vendor-gc-calibration.md` submitted 18 probes to two vendors on
2026-08-28. Eight of them carried a single 50 bp window swept from **4% to 96% GC** on a
byte-identical 50% background. All eight were **green at IDT and Standard at Twist** --
sixteen of sixteen verdicts, spanning a range far wider than the 10%/90% trigger the
brief quotes. Its conclusion is flat: *"A single extreme 50 bp window is irrelevant to
both vendors."*

So the published trigger is a threshold no vendor was observed to apply. Encoding it at
`HARD_REPAIR` would make BT5 refuse to emit constructs both vendors manufacture without
comment -- the same failure `docs/design/repeats.md` §2 records for Twist's published
repeat trigger, and the same shape as the E4 row `brief.md:141` struck through on
2026-08-28 for being below the chance floor. This rule reports it and does not gate on it.

**What this rule DOES carry is the one windowed rule a vendor was observed to apply.**
IDT's own remediation text, quoted in the calibration: *"a window of **100 bases**
starting at base 117 with a GC content of 15%. Redesign this region to have a GC content
**greater than 30%**."* That is a 100 bp **floor** with **no ceiling** -- and the
asymmetry is confirmed from the other side, since `GLB_gc80` produced no window finding at
all despite being 80% GC throughout. The calibration's own recommendation is *"Keep a
windowed floor, at IDT's geometry."* This rule is that, at the weight the next paragraph
explains.

**Why a two-sided windowed band is the wrong shape**, in the calibration's words: *"IDT
enforces a windowed floor and no windowed ceiling at all."* It also independently matches
SCP4ssd's finding that *"local fragments with low GC content might have a more important
impact than fragments with high GC content"* -- an ML model trained on synthesis outcomes
and a production order checker agreeing from opposite directions.

**And even IDT's floor is SOFT, which is why this whole rule is.** The floor is a
scoring contributor, not a gate. The calibration says so while explaining the LOC
results: *"a 50 bp window at 4% GC, diluted by 50 bp of 50% background, is ~27% across
100 bp -- barely under the 30% floor, **scoring a few points and staying green**."* Those
probes tripped IDT's windowed rule and were accepted anyway, because what denies an order
is the TOTAL complexity score reaching 24, and for GC that total is driven by GLOBAL GC
(`score = 1.40 x GC% - 83.8`, denial at 77%). So no windowed GC threshold was measured to
refuse anything at either vendor. A `HARD_REPAIR` windowed rule -- including one built on
IDT's own stated number -- would refuse constructs IDT scores green, which is the same
mistake as the 50 bp trigger one step further in. Global GC is what gates, and that is
`e2_gc_band`'s job, not this rule's.

**Overlap with `e2_gc_band`.** E2 carries `brief_ref = "2.E2/2.E3"` and already implements
a 50 bp windowed band. This rule does not touch it: E2's rework is X7 in the calibration
(*"E2 as shipped ... is wrong in every dimension at once"*), an owner-level change to what
the app refuses to build. See `docs/decisions/2026-09-02-e3-windowed-gc-measured.md`.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import ClassVar

import numpy as np

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

# `window_gc` and `merge_regions` live in f5_at_window. Both are S4-lane rules doing
# the same windowed-GC geometry, and the wrap case in `merge_regions` is subtle enough
# that two copies would drift apart -- which is exactly how a circular construct ends up
# reporting one contiguous region as two.
from bt5.rules.catalog.f5_at_window import merge_regions, window_gc

#: IDT's geometry, from its own remediation text in
#: `docs/design/vendor-gc-calibration.md`: "a window of 100 bases ... Redesign this
#: region to have a GC content greater than 30%". A FLOOR, with no ceiling.
IDT_WINDOW = 100
IDT_FLOOR = 0.30

#: `brief.md:140`, the GenScript GenTitan 100 bp warn band.
GENSCRIPT_WINDOW = 100
GENSCRIPT_LO, GENSCRIPT_HI = 0.25, 0.65

#: `brief.md:140`, Twist's published 50 bp trigger. REPORTED, NEVER GATED -- see the
#: module docstring: 8 probes from 4% to 96% GC in a 50 bp window were accepted by both
#: vendors, 16 verdicts out of 16.
TWIST_WINDOW = 50
TWIST_TRIGGER_LO, TWIST_TRIGGER_HI = 0.10, 0.90
TWIST_WARN_LO, TWIST_WARN_HI = 0.25, 0.75

#: Three tiers, all soft. The floor is the strongest reading because it is the only
#: windowed threshold a vendor was observed to score against at all.
MAG_FLOOR = 1.0
MAG_BAND = 0.3
MAG_NOTE = 0.1


@register
class WindowedGC:
    id: ClassVar[str] = "e3_windowed_gc"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "Windowed GC: IDT's 100 bp floor, with the 50 bp trigger reported"
    #: SOFT, not HARD_REPAIR. No windowed GC threshold was measured to refuse anything
    #: at either vendor -- IDT's own floor is a scored contributor and its LOC probes
    #: tripped it while staying green. See the module docstring.
    enforcement: ClassVar[Enforcement] = Enforcement.SOFT
    #: CONTESTED, not EVIDENCE_BACKED, even though brief.md:140 grades the row A. The
    #: grade attaches to the published numbers, and 18 measured probes contradict the
    #: 50 bp half of them. A badge saying "evidence backed" on a rule whose own
    #: calibration says the headline threshold does not gate would be the dishonest
    #: half of the row winning.
    evidence: ClassVar[Evidence] = Evidence.CONTESTED
    direction: ClassVar[Direction] = Direction.BAND
    unit: ClassVar[str] = "GC fraction per window"
    #: The GenScript 100 bp band from brief.md:140. Nothing here gates; this is the
    #: band the UI renders, and BAND rules must declare one.
    band: ClassVar[tuple[float, float] | None] = (GENSCRIPT_LO, GENSCRIPT_HI)
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "IDT's gBlocks complexity checker states its own windowed rule in its "
            "remediation text: a 100-base window at 15% GC, 'Redesign this region to "
            "have a GC content greater than 30%' -- a floor, with no windowed ceiling",
            "https://www.idtdna.com/pages/products/genes-and-gene-fragments/double-stranded-dna-fragments/gblocks-gene-fragments",
            2026,
            sign="supports",
        ),
        Citation(
            "GenScript GenTitan is the only vendor publishing a windowed rule, at "
            "100 bp outside 25-65% -- the source of this rule's reported band",
            "https://www.genscript.com/gentitan-gene-fragments.html",
            2026,
            sign="supports",
        ),
        Citation(
            "Twist publishes a 50 bp trigger at <10% or >90% GC, but 8 probes carrying "
            "one 50 bp window swept 4-96% GC on an identical background were accepted "
            "by BOTH vendors -- 16 verdicts of 16. The published trigger does not gate, "
            "so this rule reports it and refuses to enforce it",
            "https://www.twistbioscience.com/faq/gene-synthesis",
            2026,
            sign="refutes",
        ),
    )
    last_verified: ClassVar[str] = "2026-09-02"
    weight_provenance: ClassVar[str] = (
        "Mid. The quantity is real and measured -- IDT states a 100 bp / 30% floor and "
        "was observed to score against it -- but every windowed GC signal at both "
        "vendors was a scoring CONTRIBUTOR rather than a gate: the 18-probe calibration "
        "found 16 of 16 accepting verdicts on 50 bp windows from 4% to 96% GC, and the "
        "LOC probes tripped IDT's 100 bp floor and stayed green. So the weight carries "
        "the whole of this rule's influence, and it is set below e2_gc_band's steering "
        "because global GC is what actually gates and e2 owns it. Raising this weight "
        "would recreate, through the objective, the over-enforcement the calibration "
        "was run to remove."
    )
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.4
    #: Deliberately small. e2_gc_band already steers GC at 0.5 and f5_at_window at 0.25;
    #: a third full-strength GC steering term is one preference counted three times.
    steering_weight: ClassVar[float] = 0.1
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.WINDOW_MINUS_1
    #: A windowed statistic is monotone in the quantity it measures: raising GC in one
    #: window cannot manufacture a new below-floor window elsewhere, so a pass cannot
    #: recreate what it removed and FIXED_POINT would describe a loop that never runs
    #: (CLAUDE.md 3.6 asks for this justification when SINGLE_PASS is chosen).
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "cheap"
    #: Both are windowed GC bands on overlapping geometry; e2 additionally hard-enforces
    #: a 50 bp band this rule declines to gate on. f5_at_window pulls the same quantity
    #: the other way, from propagation toxicity rather than synthesis.
    conflicts_with: ClassVar[tuple[str, ...]] = ("e2_gc_band", "f5_at_window")
    brief_ref: ClassVar[str] = "2.E3"
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "floor": {
                "type": "number",
                "minimum": 0.0,
                "maximum": 1.0,
                "default": IDT_FLOOR,
                "description": (
                    "GC fraction below which a 100 bp window hard-fails. IDT's stated "
                    "remediation threshold, and the only windowed rule any vendor was "
                    "measured to apply."
                ),
            },
            "report_fifty_bp": {
                "type": "boolean",
                "default": True,
                "description": (
                    "Report Twist's published 50 bp trigger. Reported only -- 16 of 16 "
                    "vendor verdicts on 50 bp windows from 4% to 96% GC were accepting."
                ),
            },
        },
    }

    def __init__(self, floor: float = IDT_FLOOR, report_fifty_bp: bool = True) -> None:
        if not 0.0 <= floor <= 1.0:
            raise ValueError(f"floor must be a GC fraction in [0, 1], got {floor}")
        self.floor = floor
        self.report_fifty_bp = report_fifty_bp
        #: Read by `solver/catalog.py:236` into `RulePolicy.window`: a 100 nt statistic
        #: can only be changed by bases within 99 of it.
        self.window: int = IDT_WINDOW

    def gate(self, slot: ContextSlot) -> bool:
        """Every context: this is about the fragment a vendor synthesises, which every
        construct needs regardless of host or modality."""
        return True

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """None. The automaton makes motifs unreachable, not windowed statistics."""
        return None

    def _emit(
        self,
        c: Construct,
        fracs: np.ndarray,
        mask: np.ndarray,
        span: int,
        magnitude: float,
        describe: object,
    ) -> list[Breach]:
        """One breach per contiguous offending region, at its most extreme window."""
        out: list[Breach] = []
        for region in merge_regions(mask, circular=c.is_circular):
            block = fracs[region]
            # Furthest from the middle of the reported band is the window to name.
            local = int(np.argmax(np.abs(block - (GENSCRIPT_LO + GENSCRIPT_HI) / 2)))
            start, gc = int(region[local]), float(block[local])
            iv = Interval(start, start + span)
            tier = (
                "floor" if magnitude >= MAG_FLOOR else ("band" if magnitude >= MAG_BAND else "note")
            )
            out.append(
                Breach(
                    spec_id=self.id,
                    interval=iv,
                    # A gradient below the floor, so a partial improvement is visible
                    # in the score rather than flat until it clears.
                    magnitude=(
                        magnitude + max(0.0, self.floor - gc) * 10.0
                        if tier == "floor"
                        else magnitude
                    ),
                    message=describe(start, gc, span),  # type: ignore[operator]
                    fixable_by_codon_choice=c.overlaps_editable(iv),
                    detail={
                        "gc": gc,
                        "window": float(span),
                        "tier": tier,
                    },
                )
            )
        return out

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        text, offset = c.tripled()
        n = c.length
        breaches: list[Breach] = []

        span100 = min(IDT_WINDOW, n)
        gc100 = window_gc(text, offset, n, span100)
        if gc100.size:
            breaches += self._emit(
                c,
                gc100,
                gc100 < self.floor,
                span100,
                MAG_FLOOR,
                lambda s, gc, w: (
                    f"{w} nt window at {s}: GC {gc:.0%}, below IDT's stated {self.floor:.0%} "
                    "floor. This is the only windowed GC rule any vendor was measured to "
                    "apply, and it is a floor with no ceiling"
                ),
            )
            # The GenScript band, reported. Only the part not already hard-failed.
            warn100 = ((gc100 < GENSCRIPT_LO) | (gc100 > GENSCRIPT_HI)) & (gc100 >= self.floor)
            breaches += self._emit(
                c,
                gc100,
                warn100,
                span100,
                MAG_BAND,
                lambda s, gc, w: (
                    f"{w} nt window at {s}: GC {gc:.0%}, outside GenScript GenTitan's "
                    f"{GENSCRIPT_LO:.0%}-{GENSCRIPT_HI:.0%} published band"
                ),
            )

        if self.report_fifty_bp and n >= 1:
            span50 = min(TWIST_WINDOW, n)
            gc50 = window_gc(text, offset, n, span50)
            if gc50.size:
                trigger = (gc50 < TWIST_TRIGGER_LO) | (gc50 > TWIST_TRIGGER_HI)
                breaches += self._emit(
                    c,
                    gc50,
                    trigger,
                    span50,
                    MAG_NOTE,
                    lambda s, gc, w: (
                        f"{w} nt window at {s}: GC {gc:.0%}, past Twist's published "
                        f"{TWIST_TRIGGER_LO:.0%}/{TWIST_TRIGGER_HI:.0%} trigger. REPORTED, "
                        "NOT ENFORCED: 8 probes carrying one 50 bp window from 4% to 96% GC "
                        "were accepted by both vendors, 16 verdicts of 16"
                    ),
                )
                warn50 = ((gc50 < TWIST_WARN_LO) | (gc50 > TWIST_WARN_HI)) & ~trigger
                breaches += self._emit(
                    c,
                    gc50,
                    warn50,
                    span50,
                    MAG_NOTE,
                    lambda s, gc, w: (
                        f"{w} nt window at {s}: GC {gc:.0%}, outside the "
                        f"{TWIST_WARN_LO:.0%}-{TWIST_WARN_HI:.0%} warn band. Reported only"
                    ),
                )

        return Evaluation(
            spec_id=self.id,
            # ONLY the measured IDT floor refuses. Everything else is a reading.
            # The rule's own verdict, not a refusal: enforcement is SOFT, so
            # catalog.py never routes these breaches to the solver either way.
            passes=not any(b.detail["tier"] == "floor" for b in breaches),
            raw_score=max((b.magnitude for b in breaches), default=0.0),
            breaches=tuple(breaches),
            n_evaluated=n,
            binding_side=("lower" if any(b.detail["tier"] == "floor" for b in breaches) else None),
        )
