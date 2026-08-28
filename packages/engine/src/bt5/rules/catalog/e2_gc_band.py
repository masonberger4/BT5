"""E2 -- GC content within a manufacturability band, measured against both vendors.

REFERENCE RULE. Copy this file's shape when adding a scored rule.

**This rule gates on GLOBAL GC, and that is a measured decision.** An 18-sequence
ladder run through the IDT and Twist order-entry complexity checkers
(docs/design/vendor-gc-calibration.md) settled two things the earlier 50 bp
windowed band got wrong:

1. **No vendor gates on a single window.** Sixteen of sixteen probes carrying one
   extreme 50 bp window -- 4% to 96% GC on an otherwise neutral background --
   were Standard at Twist and green at IDT. A windowed band is the wrong
   instrument: it refuses sequences both vendors manufacture. LOC_gc05, accepted
   by both, has a 100 bp window at 23% -- LOWER than GLB_gc25 which was refused at
   21%. Only the OVERALL GC separates accepted from refused.

2. **The gate is global, two-sided, and asymmetric between vendors.** Both refuse
   <=25% and accept >=30%, so the floor is 0.28. At the top Twist ships 80% as
   Standard while IDT denies above ~77%, so the ceiling is vendor-specific: 0.80
   is the loosest demonstrated acceptance, and IDT's stricter ceiling belongs in
   the per-vendor profile (issue #43), not in this default.

So the band is (0.28, 0.80) on OVERALL GC. The old 40-60% is not a
manufacturability bound at all -- it is a design preference, and it survives only
as the steering nudge (`steering_weight`), never as a gate. `windows` is still
populated for the report's GC landscape, but a window is never a breach.

This also corrects an earlier proposal of mine: a HARD windowed FLOOR (IDT's
stated 100 bp / 30%) is wrong for the same reason -- it would refuse LOC_gc05
(23% window) while a lower-global sequence squeaks past. IDT's 30% figure is a
soft SCORED penalty, not a refuse-to-emit gate. Windowed GC belongs in a soft
rule; E4 (dGC) already carries part of it.

Two things kept deliberately:

- `Direction.BAND`, not LOWER_IS_BETTER. GC has an optimum in the middle;
  collapsing to |deviation| throws away which side is binding, which the conflict
  panel needs because manufacturability pulls GC down while suppressing cryptic
  E. coli promoters (F5) pulls it up.
- `Enforcement.HARD_REPAIR`, not HARD_LATTICE. Global GC is not decidable from a
  bounded codon suffix, so the Tier-A DP is STEERED by a Lagrangian term, Tier-B
  REPAIRS what remains, and the independent validator REFUSES TO EMIT if overall
  GC is still out of band.
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
)
from bt5.core.types import Construct, Interval
from bt5.verify import gc_fraction


@register
class GCBand:
    id: ClassVar[str] = "e2_gc_band"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "Global GC within a two-sided manufacturability band"
    enforcement: ClassVar[Enforcement] = Enforcement.HARD_REPAIR
    evidence: ClassVar[Evidence] = Evidence.VENDOR_ASSERTED
    direction: ClassVar[Direction] = Direction.BAND
    unit: ClassVar[str] = "fraction GC"
    band: ClassVar[tuple[float, float] | None] = (0.28, 0.80)
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "Measured 2026-08-28 against the IDT and Twist order-entry complexity "
            "checkers, 18-sequence ladder: both refuse global GC <=25% and accept "
            ">=30%; Twist ships 80% as Standard while IDT denies above ~77%. The "
            "floor is 0.28 and the ceiling is the loosest vendor's demonstrated "
            "acceptance. See docs/design/vendor-gc-calibration.md",
            "https://www.twistbioscience.com/faq/gene-synthesis/what-do-scoring-results-my-gene-mean",
            2026,
            sign="supports",
        ),
        Citation(
            "Twist's PUBLISHED trigger is local GC <10% or >90% over a 50 bp window, "
            "but the calibration above shows no single window drives acceptance, so "
            "the 50 bp windowed band is NOT encoded as a gate here",
            "https://www.twistbioscience.com/faq/gene-synthesis",
            2026,
            sign="qualifies",
        ),
        Citation(
            "Toxic horizontally-acquired E. coli genes run 63-68% AT vs a 55% non-toxic "
            "control -- suppressing cryptic AT-rich promoters pushes GC UP, opposing the "
            "vendor floor, which is why this is a two-sided band and F5 conflicts",
            "https://www.nature.com/articles/nmicrobiol2016249",
            2016,
            sign="qualifies",
        ),
    )
    last_verified: ClassVar[str] = "2026-08-28"
    weight_provenance: ClassVar[str] = (
        "Hard rule, so default_weight is 0 and this explains the BAND, not a weight. "
        "(0.28, 0.80) is the measured manufacturability envelope: both vendors refuse "
        "global GC <=25% and accept >=30% (floor 0.28), and the loosest vendor, Twist, "
        "ships 80% as Standard (ceiling 0.80). The 40-60% target is a design "
        "preference carried by steering_weight, not a gate; IDT's stricter ~77% "
        "ceiling is vendor-specific and deferred to the per-vendor profile (#43)."
    )
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.0  # hard rule: never in the objective
    steering_weight: ClassVar[float] = 0.5  # nudges the Tier-A DP toward the 40-60% preference
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.WINDOW_MINUS_1
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "cheap"
    conflicts_with: ClassVar[tuple[str, ...]] = ("f5_at_window", "d8_cpg_depletion")
    brief_ref: ClassVar[str] = "2.E2/2.E3"
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "gc_min": {"type": "number", "default": 0.28, "minimum": 0.0, "maximum": 1.0},
            "gc_max": {"type": "number", "default": 0.80, "minimum": 0.0, "maximum": 1.0},
            "window": {
                "type": "integer",
                "default": 50,
                "minimum": 10,
                "description": "Resolution of the reported GC landscape only; not a gate.",
            },
        },
    }

    def __init__(self, gc_min: float = 0.28, gc_max: float = 0.80, window: int = 50) -> None:
        if gc_min >= gc_max:
            raise ValueError(f"gc_min {gc_min} must be below gc_max {gc_max}")
        self.gc_min, self.gc_max, self.window = gc_min, gc_max, window

    def gate(self, slot: ContextSlot) -> bool:
        # IVT mRNA is not synthesised as a plasmid insert in the same way; the
        # manufacturability argument does not transfer unchanged.
        return slot.modality is not Modality.IVT_MRNA

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        # Deliberately None: windowed GC is not expressible in the automaton state.
        # Tier A steers via a Lagrangian term supplied by the solver, not here.
        return None

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        seq = c.sequence
        n = len(seq)
        step = max(1, self.window // 5)
        scan = seq + seq[: self.window - 1] if c.is_circular else seq

        # The windows are for the report's GC landscape ONLY. Measurement (see the
        # module docstring) shows no vendor refuses on a single window, so a window
        # is never a breach -- the gate below is on OVERALL GC.
        windows: list[tuple[Interval, float]] = []
        for start in range(0, n if c.is_circular else max(1, n - self.window + 1), step):
            chunk = scan[start : start + self.window]
            if len(chunk) < self.window:
                continue
            windows.append((Interval(start, start + self.window), gc_fraction(chunk)))

        overall = gc_fraction(seq)
        breaches: list[Breach] = []
        binding: str | None = None
        if overall < self.gc_min:
            dev, binding = self.gc_min - overall, "lower"
        elif overall > self.gc_max:
            dev, binding = overall - self.gc_max, "upper"
        if binding is not None:
            breaches.append(
                Breach(
                    spec_id=self.id,
                    interval=Interval(0, n),
                    magnitude=dev,
                    message=(
                        f"overall GC {overall:.1%} is outside the manufacturability "
                        f"band [{self.gc_min:.0%}, {self.gc_max:.0%}] ({binding} bound "
                        f"binding). A single window is not a vendor gate; global GC is"
                    ),
                    # Global GC is movable wherever the construct has a CDS to
                    # recode; a construct that is all backbone is the user's vector.
                    fixable_by_codon_choice=bool(c.editable),
                    detail={"gc": overall, "binding_side": binding},
                )
            )

        return Evaluation(
            spec_id=self.id,
            passes=not breaches,
            raw_score=overall,
            breaches=tuple(breaches),
            windows=tuple(windows),
            n_evaluated=len(windows),
            binding_side=binding,
        )
