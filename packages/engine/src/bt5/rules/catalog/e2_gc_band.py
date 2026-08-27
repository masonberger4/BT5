"""E2 -- windowed GC content, as a two-sided band.

REFERENCE RULE. Copy this file's shape when adding a scored rule.

Two things here are easy to get wrong and both are deliberate:

1. `Direction.BAND`, not LOWER_IS_BETTER. GC has an optimum in the middle, and
   collapsing a band to |deviation| throws away WHICH SIDE is binding -- which is
   exactly what the conflict panel needs, because manufacturability pushes GC down
   while suppressing cryptic E. coli promoters pushes it up. The report must be
   able to say which constraint is binding in which window.

2. `Enforcement.HARD_REPAIR`, not HARD_LATTICE. Deciding whether a codon pushes a
   50 bp window past its GC bound requires the G+C count over the previous ~17
   codons, and enumerating that history in the DP state is combinatorially
   intractable with no literature precedent. So: the Tier-A DP is STEERED by a
   Lagrangian term, Tier-B REPAIRS what remains, and the independent validator
   REFUSES TO EMIT if a window is still out of band. A violating sequence never
   reaches the user -- worst case the app declines and reports the conflict.

Numbers: Twist's PUBLISHED 50 bp window bound is 10-90%. The widely repeated
"35-65%" has no vendor source behind it and is not encoded here.
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
    title: ClassVar[str] = "Windowed GC content within a two-sided band"
    enforcement: ClassVar[Enforcement] = Enforcement.HARD_REPAIR
    evidence: ClassVar[Evidence] = Evidence.VENDOR_ASSERTED
    direction: ClassVar[Direction] = Direction.BAND
    unit: ClassVar[str] = "fraction GC"
    band: ClassVar[tuple[float, float] | None] = (0.40, 0.60)
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "Twist gene synthesis FAQ: High-Complexity trigger is local GC below 10% "
            "or above 90% in a 50 bp window (NOT the widely repeated 35-65%)",
            "https://www.twistbioscience.com/faq/gene-synthesis",
            2026,
            sign="supports",
        ),
        Citation(
            "Toxic horizontally-acquired E. coli genes run 63-68% AT vs a 55% non-toxic "
            "control -- suppressing cryptic AT-rich promoters pushes GC UP, opposing the "
            "vendor ceiling",
            "https://www.nature.com/articles/nmicrobiol2016249",
            2016,
            sign="qualifies",
        ),
    )
    last_verified: ClassVar[str] = "2026-08-27"
    weight_provenance: ClassVar[str] = (
        "Moderate default weight. The vendor acceptance envelope is genuinely wide "
        "(Twist rejects only <10% or >90% over 50 bp), so a narrow target is a "
        "preference rather than a hard requirement; but GC extremes correlate with "
        "synthesis surcharges and assembly failure, and repeats -- not GC -- are the "
        "top predictor of synthesis failure, so GC should not dominate the objective."
    )
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.0  # hard rule: never in the objective
    steering_weight: ClassVar[float] = 0.5  # nudges the Tier-A DP toward the band
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.WINDOW_MINUS_1
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "cheap"
    conflicts_with: ClassVar[tuple[str, ...]] = ("f5_at_window", "d8_cpg_depletion")
    brief_ref: ClassVar[str] = "2.E2/2.E3"
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "gc_min": {"type": "number", "default": 0.40, "minimum": 0.0, "maximum": 1.0},
            "gc_max": {"type": "number", "default": 0.60, "minimum": 0.0, "maximum": 1.0},
            "window": {"type": "integer", "default": 50, "minimum": 10},
        },
    }

    def __init__(self, gc_min: float = 0.40, gc_max: float = 0.60, window: int = 50) -> None:
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

        windows: list[tuple[Interval, float]] = []
        breaches: list[Breach] = []
        worst_dev, binding = 0.0, None

        for start in range(0, n if c.is_circular else max(1, n - self.window + 1), step):
            chunk = scan[start : start + self.window]
            if len(chunk) < self.window:
                continue
            frac = gc_fraction(chunk)
            iv = Interval(start, start + self.window)
            windows.append((iv, frac))

            if frac < self.gc_min:
                dev, side = self.gc_min - frac, "lower"
            elif frac > self.gc_max:
                dev, side = frac - self.gc_max, "upper"
            else:
                continue

            breaches.append(
                Breach(
                    spec_id=self.id,
                    interval=iv,
                    magnitude=dev,
                    message=(
                        f"GC {frac:.1%} in the {self.window}bp window at {start} is outside "
                        f"[{self.gc_min:.0%}, {self.gc_max:.0%}] ({side} bound binding)"
                    ),
                    detail={"gc": frac, "binding_side": side},
                )
            )
            if dev > worst_dev:
                worst_dev, binding = dev, side

        return Evaluation(
            spec_id=self.id,
            passes=not breaches,
            raw_score=gc_fraction(seq),
            breaches=tuple(breaches),
            windows=tuple(windows),
            n_evaluated=len(windows),
            binding_side=binding,
        )
