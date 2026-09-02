"""F5 -- the AT-window rule, as a two-sided band the brief itself asks for.

`brief.md:159`: *"every 100-nt window <=55% AT (>=45% GC); none above 60% AT. Toxic
horizontally-acquired E. coli genes are 63-68% AT vs a non-toxic 55% control ...
**This directly conflicts with vendor GC ceilings** -- resolve as a two-sided band
45-60% GC per 100 nt, hard-fail outside 35-65%, and show which side is binding per
window."* Evidence grade **A**.

Three things that follow from that sentence and are easy to get wrong:

**It is a GC band, stated in AT.** The row is written in AT because the toxicity
evidence is (63-68% AT genes vs a 55% control), but the brief converts it itself:
<=55% AT is >=45% GC. Everything here is stored as GC so it can be compared with
`e2_gc_band` without a mental flip on every line, and the AT figure is echoed in the
message because that is the form the evidence is in.

**The conflict with the vendor ceiling is real, declared, and not resolved here.**
Suppressing cryptic AT-rich promoters pushes GC *up*; the vendor floor and IDT's
denial above ~77% global GC push it *down*. `e2_gc_band` already names this rule in its
`conflicts_with` and carries the same Nature citation with `sign="qualifies"`; this rule
names it back. Neither silently wins -- the conflict panel is what surfaces it.

**"Show which side is binding per window" is a requirement, not a nicety.** A
`|deviation|` scalar cannot say whether a window needs GC raised or lowered, and on a
two-sided band with two different mechanisms behind the two sides, that is the only
actionable part of the finding. Each breach carries `detail["binding_side"]`, and
`Evaluation.binding_side` reports the worst window's side.

**This is a propagation rule, not a synthesis rule.** Section `brief.md:151` scopes 2.F
to *"every construct that passes through a cloning host"*, so it gates everywhere: the
plasmid is grown in E. coli whatever the final host is. That also means it is NOT
governed by `docs/design/vendor-gc-calibration.md`, which measured what vendors refuse to
SYNTHESISE. The two are different questions about the same number.
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

#: `brief.md:159`: "every 100-nt window".
WINDOW = 100

#: `brief.md:159`: "resolve as a two-sided band 45-60% GC per 100 nt". The preference.
SOFT_LO, SOFT_HI = 0.45, 0.60

#: `brief.md:159`: "hard-fail outside 35-65%". The gate.
HARD_LO, HARD_HI = 0.35, 0.65

#: The same band expressed as AT, which is the form the toxicity evidence is in.
#: 45% GC == 55% AT; 40% GC == 60% AT.
AT_SOFT_HI = 1.0 - SOFT_LO
AT_HARD_HI = 1.0 - HARD_LO

MAG_WARN = 0.3
MAG_HARD_BASE = 1.0
#: A gradient, not a flag. `_accepts` (`solver/repair.py:270-312`) needs a strictly
#: falling magnitude sum to accept a move that does not clear the breach outright, so a
#: constant hard magnitude would make "62% GC" and "95% GC" look identical to the search.
MAG_PER_GC_POINT = 10.0


def window_gc(text: str, offset: int, n: int, span: int) -> np.ndarray:
    """GC fraction of every `span`-nt window starting in the middle copy.

    Vectorised, and that is not premature. HARD_REPAIR means `breach_finder` calls this
    once per candidate, up to 256 per repair iteration: a Python loop over 5,000 window
    starts costs ~4.7 ms here, which is 1.2 s per iteration before any other rule runs.
    That is the same shape of cost that turned a 9-second design test into one that had
    not finished in 500 seconds earlier in this branch. A cumulative sum over the
    tripled text makes every window O(1) and the whole scan one pass.
    """
    arr = np.frombuffer(text.encode("ascii"), dtype=np.uint8)
    is_gc = (arr == ord("G")) | (arr == ord("C"))
    prefix = np.concatenate(([0], np.cumsum(is_gc, dtype=np.int32)))
    starts = np.arange(offset, min(offset + n, len(text) - span + 1), dtype=np.int64)
    if starts.size == 0:
        return np.empty(0, dtype=np.float64)
    return (prefix[starts + span] - prefix[starts]) / span


@register
class ATWindow:
    id: ClassVar[str] = "f5_at_window"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "AT-window rule: two-sided GC band per 100 nt"
    enforcement: ClassVar[Enforcement] = Enforcement.HARD_REPAIR
    evidence: ClassVar[Evidence] = Evidence.EVIDENCE_BACKED  # brief.md:159 grades it A
    direction: ClassVar[Direction] = Direction.BAND
    unit: ClassVar[str] = "GC fraction per 100 nt window"
    #: The PREFERENCE band brief.md:159 names, not the gate. The gate is
    #: HARD_LO/HARD_HI; `band` is what the UI renders and what BAND rules declare.
    band: ClassVar[tuple[float, float] | None] = (SOFT_LO, SOFT_HI)
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "Toxic horizontally-acquired E. coli genes run 63-68% AT against a 55% AT "
            "non-toxic control -- the measurement the 55% AT (45% GC) ceiling comes from",
            "https://www.nature.com/articles/nmicrobiol2016249",
            2016,
            sign="supports",
        ),
        Citation(
            "The same band read from the other side: suppressing cryptic AT-rich "
            "promoters pushes GC UP, against the vendor ceiling, which is why this is a "
            "two-sided band and why it conflicts with e2_gc_band rather than refining it",
            "https://www.nature.com/articles/nmicrobiol2016249",
            2016,
            sign="qualifies",
        ),
    )
    last_verified: ClassVar[str] = "2026-09-02"
    weight_provenance: ClassVar[str] = ""  # hard rule; the band is the gate, not a weight
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.0  # hard rule: never in the objective
    #: Nudges the Tier-A DP toward the middle of the preference band. Lower than
    #: e2_gc_band's 0.5 deliberately: two GC steering terms pulling at full strength on
    #: overlapping windows is one preference counted twice.
    steering_weight: ClassVar[float] = 0.25
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.WINDOW_MINUS_1
    #: A windowed statistic cannot recreate itself the way a splice donor can: raising
    #: GC in one window does not manufacture a new out-of-band window elsewhere, because
    #: the change is monotone in the quantity being measured. SINGLE_PASS is sufficient
    #: (CLAUDE.md 3.6 asks for the justification when it is chosen).
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "cheap"
    #: brief.md:159 says so in bold: "This directly conflicts with vendor GC ceilings."
    conflicts_with: ClassVar[tuple[str, ...]] = ("e2_gc_band",)
    brief_ref: ClassVar[str] = "2.F5"
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "window": {
                "type": "integer",
                "minimum": 20,
                "default": WINDOW,
                "description": "Window size in nt. brief.md:159 states 100.",
            },
            "hard_lo": {
                "type": "number",
                "minimum": 0.0,
                "maximum": 1.0,
                "default": HARD_LO,
                "description": "GC fraction below which a window hard-fails.",
            },
            "hard_hi": {
                "type": "number",
                "minimum": 0.0,
                "maximum": 1.0,
                "default": HARD_HI,
                "description": "GC fraction above which a window hard-fails.",
            },
        },
    }

    def __init__(
        self,
        window: int = WINDOW,
        hard_lo: float = HARD_LO,
        hard_hi: float = HARD_HI,
    ) -> None:
        if hard_lo >= hard_hi:
            raise ValueError(f"hard band must be non-inverted, got ({hard_lo}, {hard_hi})")
        if window < 1:
            raise ValueError(f"window must be positive, got {window}")
        #: Read by `solver/catalog.py:236` into `RulePolicy.window`, which is why
        #: localization is WINDOW_MINUS_1: a 100 nt statistic can only be changed by
        #: bases within 99 of it, and catalog.py hard-codes motif_len to 6 for everyone.
        self.window = window
        self.hard_lo = hard_lo
        self.hard_hi = hard_hi

    def gate(self, slot: ContextSlot) -> bool:
        """Everywhere. brief.md:151 scopes 2.F to "every construct that passes through
        a cloning host", which is every construct BT5 builds -- the plasmid is grown in
        E. coli whatever the final host is."""
        return True

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """None. A windowed compositional band is not a finite motif set; the automaton
        makes motifs unreachable, not statistics. Steering plus Tier-B repair is the
        mechanism, which is what `steering_weight` above is for."""
        return None

    def _windows(self, c: Construct) -> np.ndarray:
        """GC fraction per window, indexed by start in construct coordinates."""
        text, offset = c.tripled()
        return window_gc(text, offset, c.length, min(self.window, c.length))

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        n = c.length
        span = min(self.window, n)
        fracs = self._windows(c)
        if fracs.size == 0:
            return Evaluation(spec_id=self.id, passes=True, raw_score=0.0, n_evaluated=0)

        # Merge consecutive offending windows into regions and report one breach each.
        # A 5 kb plasmid at 30% GC has ~5,000 offending windows; emitting one breach
        # per window would make this rule's breach COUNT -- the currency `_aggregate`
        # steers on (`solver/repair.py:253-267`) -- swamp every other rule in the
        # catalog for what is one contiguous problem.
        breaches: list[Breach] = []
        worst_side: str | None = None
        worst_dev = 0.0

        offending = (fracs < SOFT_LO) | (fracs > SOFT_HI)
        # Run boundaries from the mask itself rather than a Python sweep: `np.diff` on
        # the padded mask marks every rising and falling edge at once.
        edges = np.diff(np.concatenate(([0], offending.view(np.int8), [0])))
        run_starts = np.flatnonzero(edges == 1)
        run_ends = np.flatnonzero(edges == -1)
        regions = [np.arange(a, b) for a, b in zip(run_starts, run_ends, strict=True)]

        # An offending stretch that crosses the origin arrives as two runs -- one
        # ending at the last window start, one beginning at the first -- because the
        # mask is a line and the molecule is not. Reporting it as two regions would
        # double this rule's breach count on exactly the construct where the finding
        # is one contiguous problem, and would name two "worst" windows for it.
        if c.is_circular and len(regions) > 1 and offending[0] and offending[-1]:
            regions[0] = np.concatenate((regions[-1], regions[0]))
            regions.pop()

        def flush(idx: np.ndarray) -> None:
            nonlocal worst_side, worst_dev
            block = fracs[idx]
            # The window that sits furthest outside the hard band is the one to name.
            headroom = np.minimum(block - self.hard_lo, self.hard_hi - block)
            local = int(np.argmin(headroom))
            start, gc = int(idx[local]), float(block[local])
            below = gc < self.hard_lo
            hard = below or gc > self.hard_hi
            side = "lower" if gc < (self.hard_lo + self.hard_hi) / 2 else "upper"
            excess = (self.hard_lo - gc) if below else max(0.0, gc - self.hard_hi)
            if not hard:
                excess = 0.0
            deviation = abs(gc - (SOFT_LO if side == "lower" else SOFT_HI))
            if deviation > worst_dev:
                worst_dev, worst_side = deviation, side
            iv = Interval(start, start + span)
            breaches.append(
                Breach(
                    spec_id=self.id,
                    interval=iv,
                    magnitude=(
                        MAG_HARD_BASE + excess * 100.0 / MAG_PER_GC_POINT if hard else MAG_WARN
                    ),
                    message=(
                        f"{span} nt window at {start}: GC {gc:.0%} ({1 - gc:.0%} AT), "
                        f"{'outside' if hard else 'inside'} the hard band "
                        f"{self.hard_lo:.0%}-{self.hard_hi:.0%} and outside the "
                        f"{SOFT_LO:.0%}-{SOFT_HI:.0%} preference. The {side} bound is "
                        "binding"
                        + (
                            f"; AT above {AT_HARD_HI:.0%} is the range toxic "
                            "horizontally-acquired E. coli genes occupy"
                            if side == "lower" and hard
                            else ""
                        )
                    ),
                    fixable_by_codon_choice=c.overlaps_editable(iv),
                    detail={
                        "gc": gc,
                        "at": 1.0 - gc,
                        "binding_side": side,
                        "hard": "yes" if hard else "no",
                    },
                )
            )

        for region in regions:
            flush(region)

        return Evaluation(
            spec_id=self.id,
            # Only the HARD band refuses. A window inside 35-65% but outside the 45-60%
            # preference is a warn-band finding, and `solver/catalog.py:158-170` is
            # explicit that handing those to repair sets it chasing a threshold that was
            # never crossed until the search stagnates on a design the catalog accepts.
            passes=not any(b.magnitude >= MAG_HARD_BASE for b in breaches),
            raw_score=max((abs(b.magnitude) for b in breaches), default=0.0),
            breaches=tuple(breaches),
            n_evaluated=n,
            binding_side=worst_side,
        )
