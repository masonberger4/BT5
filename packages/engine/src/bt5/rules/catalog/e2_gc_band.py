"""E2 -- GC content within a manufacturability band, measured on the ordered fragment.

REFERENCE RULE. Copy this file's shape when adding a scored rule.

**This rule gates on the GC of the SYNTHESISED FRAGMENT, and that is a measured
decision.** An 18-sequence ladder run through the IDT and Twist order-entry
complexity checkers (docs/design/vendor-gc-calibration.md) settled two things the
earlier 50 bp windowed band got wrong:

1. **No vendor gates on a single window.** Sixteen of sixteen probes carrying one
   extreme 50 bp window -- 4% to 96% GC on an otherwise neutral background --
   were Standard at Twist and green at IDT. A windowed band is the wrong
   instrument: it refuses sequences both vendors manufacture. LOC_gc05, accepted
   by both, has a 100 bp window at 23% -- LOWER than GLB_gc25 which was refused at
   21%. Only the OVERALL GC OF THE FRAGMENT separates accepted from refused.

2. **The gate is per-fragment, two-sided, and asymmetric between vendors.** Both
   refuse <=25% and accept >=30%, so the floor is 0.28. At the top Twist ships 80%
   as Standard while IDT denies above ~77%, so the ceiling is vendor-specific, and
   the effective ceiling is the SELECTED vendor's: 0.77 under the gBlocks default,
   0.80 for a Twist order, and the intersection across a multi-vendor selection.

**Scope is the fragment, not the whole construct.** The vendor's complexity
checker measures the tube they synthesise -- the ordered insert plus any adapters
-- not the assembled plasmid. E2 therefore evaluates `fragments(c, adapters)` one
tube at a time, exactly as E4-E7 do, and an order succeeds only if every tube
does. Measuring `c.sequence` instead (as this rule did before #43 V3b) dragged a
GC-extreme insert toward 0.5 against a near-neutral backbone, so the one
HARD_REPAIR rule whose job is refusing unbuildable DNA never fired on a vector
design: a 900 bp / 90% insert reads 0.629 across a 2 kb / 50% backbone and passed,
though both vendors deny the fragment. The band is per-vendor; the scope is the
fragment; the two only mean anything together.

The `band` ClassVar is the LOOSEST demonstrated envelope (0.28, 0.80) -- the outer
bound of what any single configuration accepts, kept non-None and `lo < hi` for
the contract. It is not the gate. The gate is `self.gc_min`/`self.gc_max`, the
intersection of the SELECTED vendors' bands (or an explicit override), computed
per instance. A single advertised number cannot be the gate once the ceiling is
vendor-specific.

The old 40-60% is not a manufacturability bound at all -- it is a design
preference, and it survives only as the steering nudge (`steering_weight`), never
as a gate. `windows` is still populated for the report's GC landscape, now per
fragment and mapped back to parent coordinates, but a window is never a breach.

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
- `Enforcement.HARD_REPAIR`, not HARD_LATTICE. A fragment's GC is not decidable
  from a bounded codon suffix, so the Tier-A DP is STEERED by a Lagrangian term,
  Tier-B REPAIRS what remains, and the independent validator REFUSES TO EMIT if a
  fragment is still out of band. Infeasibility from the protein side cannot
  happen: no amino acid forces a per-codon GC above 2/3 (poly-Ala/Gly/Pro reach
  0.656-0.989), so an intersected band with a 0.28 floor and a >=0.77 ceiling is
  always reachable by recoding.
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
from bt5.rules.fragment import fragments
from bt5.rules.vendors import (
    DEFAULT_SELECTION,
    DEFAULT_VENDOR,
    VendorSelection,
    all_keys,
    require_selection,
)
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
        "Hard rule, so default_weight is 0 and this explains the BAND ClassVar, not a "
        "weight. (0.28, 0.80) is the LOOSEST demonstrated envelope, not the gate: both "
        "vendors refuse GC <=25% and accept >=30% (floor 0.28), and the loosest vendor, "
        "Twist, ships 80% as Standard (ceiling 0.80). The actual gate is the SELECTED "
        "vendors' intersected band (IDT's 0.77 under the gBlocks default), computed per "
        "instance. The 40-60% target is a design preference carried by steering_weight."
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
            # No numeric default: the gate is the selected vendors' intersected
            # band, so a single advertised number would lie once the ceiling is
            # vendor-specific. A value here overrides that side of the band.
            "gc_min": {"type": "number", "minimum": 0.0, "maximum": 1.0},
            "gc_max": {"type": "number", "minimum": 0.0, "maximum": 1.0},
            "vendors": {
                "type": "array",
                "minItems": 1,
                "uniqueItems": True,
                # all_keys: like E4-E7, E2 answers the scope question even for
                # `none` (the fragment is linear, the backbone is not synthesised),
                # falling back to the loosest envelope when no vendor band exists.
                "items": {"type": "string", "enum": list(all_keys())},
                "default": [DEFAULT_VENDOR],
                "description": (
                    "Which vendor configurations the fragment is ordered as. The "
                    "band is the intersection of their published GC bands; explicit "
                    "gc_min/gc_max override it."
                ),
            },
            "window": {
                "type": "integer",
                "default": 50,
                "minimum": 10,
                "description": "Resolution of the reported GC landscape only; not a gate.",
            },
        },
    }

    def __init__(
        self,
        vendors: VendorSelection = DEFAULT_SELECTION,
        gc_min: float | None = None,
        gc_max: float | None = None,
        window: int = 50,
    ) -> None:
        # `require_selection`, not `orderable_only`: E2 has an answer even with no
        # vendor chosen (the scope is the fragment either way), so it accepts
        # `none` like E4-E7 and falls back to the loosest envelope for its band.
        self.vendors = require_selection(vendors)
        band = self.vendors.gc_band()
        lo_keys: tuple[str, ...] = ()
        hi_keys: tuple[str, ...] = ()
        if band is None:
            # `none`: no vendor publishes a band, so refuse only what no
            # configuration would make -- the ClassVar's loosest envelope, which
            # is a concrete tuple for this rule.
            assert self.band is not None
            lo, hi = self.band
        else:
            (lo, lo_keys), (hi, hi_keys) = band
        # An explicit bound is the user's own number: enforced, attributed to
        # nobody. Tracked per side so a gc_min override does not silence gc_max's
        # vendor, the same discipline E1 applies per axis.
        self._lo_override = gc_min is not None
        self._hi_override = gc_max is not None
        self._floor_binders = () if self._lo_override else lo_keys
        self._ceil_binders = () if self._hi_override else hi_keys
        self._floor_from_envelope = band is None and not self._lo_override
        self._ceil_from_envelope = band is None and not self._hi_override
        self.gc_min = lo if gc_min is None else gc_min
        self.gc_max = hi if gc_max is None else gc_max
        self.window = window
        if self.gc_min >= self.gc_max:
            raise ValueError(f"gc_min {self.gc_min} must be below gc_max {self.gc_max}")

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
        adapters = self.vendors.adapters
        step = max(1, self.window // 5)
        mid = (self.gc_min + self.gc_max) / 2

        windows: list[tuple[Interval, float]] = []
        breaches: list[Breach] = []
        # The representative GC for the report and the steering term: the fragment
        # furthest outside the band, or -- if all are in band -- the one furthest
        # from its middle. A breach always outranks an in-band fragment.
        rep_key = float("-inf")
        rep_gc: float | None = None
        rep_side: str | None = None

        for frag in fragments(c, adapters):
            seq = frag.sequence
            n = len(seq)
            gc = gc_fraction(seq)

            # Windows are the report's GC landscape ONLY, mapped back to parent
            # coordinates; a window is never a breach. A window wholly inside a
            # vendor adapter has no parent coordinate and is dropped.
            for start in range(0, max(1, n - self.window + 1), step):
                chunk = seq[start : start + self.window]
                if len(chunk) < self.window:
                    continue
                parent = frag.to_construct(Interval(start, start + self.window))
                if parent is None:
                    continue
                windows.append((parent, gc_fraction(chunk)))

            side: str | None = None
            dev = 0.0
            binders: tuple[str, ...] = ()
            from_envelope = overridden = False
            edge = ""
            if gc < self.gc_min:
                dev, side = self.gc_min - gc, "lower"
                binders, from_envelope = self._floor_binders, self._floor_from_envelope
                overridden, edge = self._lo_override, f"under the {self.gc_min:.0%} floor"
            elif gc > self.gc_max:
                dev, side = gc - self.gc_max, "upper"
                binders, from_envelope = self._ceil_binders, self._ceil_from_envelope
                overridden, edge = self._hi_override, f"over the {self.gc_max:.0%} ceiling"

            key = dev if side is not None else -abs(gc - mid) - 1.0
            if key > rep_key:
                rep_key, rep_gc, rep_side = key, gc, side

            if side is None:
                continue

            if overridden:
                source = "the band you set"
            elif from_envelope:
                source = "the loosest envelope, no vendor chosen"
            else:
                source = f"the {', '.join(binders)} band"
            breaches.append(
                Breach(
                    spec_id=self.id,
                    # The parent coordinates of the ordered DNA -- a real interval
                    # WINDOW_MINUS_1 can localize and codon_span can walk, unlike
                    # the whole-construct Interval(0, n) this rule emitted before.
                    interval=frag.origin,
                    magnitude=dev,
                    message=(
                        f"ordered fragment GC {gc:.1%} is outside {source} "
                        f"[{self.gc_min:.0%}, {self.gc_max:.0%}] ({side} bound binding, "
                        f"{edge}). A single window is not a vendor gate; the whole "
                        f"synthesized fragment is"
                    ),
                    # The fragment is recodeable CDS and no amino acid forces GC
                    # outside the band, so this is always honestly fixable -- the
                    # old construct-scope True-whenever-a-CDS-exists could chase a
                    # GC-rich immutable backbone it never had the freedom to reach.
                    fixable_by_codon_choice=frag.construct.overlaps_editable(frag.ordered),
                    detail={
                        "gc": gc,
                        "binding_side": side,
                        "vendor": ", ".join(binders),
                        "band_lo": self.gc_min,
                        "band_hi": self.gc_max,
                        "limit_source": "override"
                        if overridden
                        else ("envelope" if from_envelope else "vendor"),
                        "fragment_bp": float(frag.ordered.length),
                    },
                )
            )

        return Evaluation(
            spec_id=self.id,
            passes=not breaches,
            raw_score=rep_gc if rep_gc is not None else mid,
            breaches=tuple(breaches),
            windows=tuple(windows),
            n_evaluated=len(windows),
            binding_side=rep_side,
        )
