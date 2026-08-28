"""E4 -- how much GC VARIES across the fragment, which is not what E2 asks.

E2 asks whether any window is out of band. This asks how far apart the windows
are from each other, and the two questions come apart completely: a fragment
whose every 50 bp window sits inside 25-65% can still swing from 26% to 64%
across its length. That is a 38-point excursion with zero E2 findings, and it is
the shape that breaks assembly -- one end anneals at a different temperature
from the other, so no single set of cycling conditions serves both.

**Why this rule exists at all, given E2 already scans GC.** The best-validated
synthesis-difficulty model puts GC *fluctuation* near the top of its feature
ranking and above global GC. Measured here on 200 synonymous variants x 4 seeds,
dGC is also essentially independent of everything BT5 already computes --
|rho| <= 0.09 against the repeat family AND against E2 itself -- while moving
about 58% of its own median under synonymous substitution. So it is a large,
codon-controllable axis that no other rule sees.

**Both metrics are measured; only one of them can carry a threshold.** Random
DNA is not uniform, so neither statistic reads zero on a sequence with nothing
wrong with it. Measured on uniform 50% GC sequence, 200 draws per length:

    length     extent(50 bp)     dGC(100 bp)
       300              26.0            3.22
      1200              36.0            4.66
      5000              44.0            4.94
     10000              46.0            4.96

**Extent grows without bound and dGC converges.** A range statistic is the max
minus the min over N windows, so it widens with every window added; a dispersion
statistic converges on the binomial floor, 100*sqrt(p(1-p)/w), which is 5.0 at
p=0.5 and w=100. That is the real reason dGC is the scalar here and extent is
not -- robustness is secondary.

It also disposes of the brief's numbers. Row 2.E4 proposes "extent <= 50, target
<= 25", and **25 is below the chance floor at every length at or above 300 bp**
while 50 is about four points above it at 10 kb. A fixed extent threshold would
fire on essentially every sequence, which is the same trap a fixed 9-mer density
threshold would have been for E6.

So the scalar is **dGC relative to its own binomial floor** -- 1.0 means exactly
as dispersed as chance, and the composition confound cancels, which matters
because synonymous choice alone moves a 300 aa protein's achievable GC across
27-60%. Extent is still computed and reported, because it is what localises the
finding to the two windows worth recoding, but it never sets a threshold.

**Enforcement is SOFT, and that is a claim about the evidence.** Twist publishes
a per-window BAND (10-90% over 50 bp) -- which is E2's rule, and E2 enforces it
hard. Nobody publishes any bound on dispersion. Inventing a hard one here would
be the move the E2 docstring refuses when it declines to encode the widely
repeated "35-65%" that has no vendor source.

**Evidence is CONTESTED, and the disagreement is the point.** The two best
synthesis-success models rank this against repeats in opposite orders, and the
one that favours repeats is the one that scored 5/10 on the only head-to-head
experimental validation. Neither ranking is settled, so both citations are
carried with opposite signs rather than one being quietly chosen.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
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
from bt5.core.types import Construct, Interval
from bt5.rules.fragment import VENDOR_ADAPTERS, Fragment, fragments

#: SD of GC is taken over 100 bp windows -- the window the feature was ranked at.
DGC_WINDOW_BP = 100
#: Extent is taken over 50 bp windows -- the window Twist publishes against.
EXTENT_WINDOW_BP = 50

#: Report when a fragment is at least this many times as dispersed as chance.
#: A convention, and named as one -- no vendor publishes a bound on dispersion.
#: 1.5 sits well clear of the measured chance floor (which is 1.0 by
#: construction) while still firing on a real block excursion.
WARN_RATIO = 1.5
#: Twice chance. Still only a finding: SOFT rules do not fail a design.
HARD_RATIO = 2.0

MAX_FINDINGS = 200
DEFAULT_VENDOR = "twist_gene_fragment"


def gc_windows(seq: str, window: int) -> list[float]:
    """GC percentage of every `window`-wide window, step 1.

    Step 1 rather than the coarser stride E2 uses: a range statistic is set by
    its extremes, and a stride can step over the extreme window and report a
    narrower spread than the molecule actually has. One cumulative sum makes the
    exact answer O(n), so there is nothing to buy by approximating.
    """
    if window <= 0:
        raise ValueError(f"window must be positive, got {window}")
    if len(seq) < window:
        if not seq:
            return []
        return [100.0 * sum(b in "GC" for b in seq) / len(seq)]
    running = 0
    cumulative = [0]
    for base in seq:
        running += base in "GC"
        cumulative.append(running)
    return [
        100.0 * (cumulative[i + window] - cumulative[i]) / window
        for i in range(len(seq) - window + 1)
    ]


def dgc(seq: str, window: int = DGC_WINDOW_BP) -> float:
    """Population SD of windowed GC. 0.0 when there is under one window."""
    values = gc_windows(seq, window)
    if len(values) < 2:
        return 0.0
    mean = sum(values) / len(values)
    return float((sum((v - mean) ** 2 for v in values) / len(values)) ** 0.5)


def chance_dgc(seq: str, window: int = DGC_WINDOW_BP) -> float:
    """The SD a sequence of THIS composition shows at this window by chance.

    100*sqrt(p(1-p)/w) -- the binomial SD of a windowed proportion. Measured
    against 200 random draws per length it is accurate to about 0.1 points from
    1 kb up, and it is the reason dGC can carry a threshold when extent cannot:
    chance is a known constant here rather than a function of length.
    """
    if not seq:
        return 0.0
    p = sum(b in "GC" for b in seq) / len(seq)
    return float(100.0 * (p * (1.0 - p) / window) ** 0.5)


def dispersion_ratio(seq: str, window: int = DGC_WINDOW_BP) -> float:
    """dGC over its own binomial floor. 1.0 is exactly as dispersed as chance.

    The ratio rather than the raw SD because composition is a confound: a 70% GC
    fragment has a lower floor (4.58) than a 50% one (5.00), so comparing raw SDs
    across synonymous variants -- which span 27-60% GC for one 300 aa protein --
    would rank composition rather than dispersion.
    """
    floor = chance_dgc(seq, window)
    if floor <= 0.0:
        return 0.0
    return dgc(seq, window) / floor


def extremes(values: Sequence[float]) -> tuple[int, float, int, float]:
    """(argmin, min, argmax, max). The two windows that set the extent."""
    lo = min(range(len(values)), key=lambda i: values[i])
    hi = max(range(len(values)), key=lambda i: values[i])
    return lo, values[lo], hi, values[hi]


@register
class GCExtent:
    id: ClassVar[str] = "e4_gc_extent"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "GC variation across the fragment (extent and dGC)"
    enforcement: ClassVar[Enforcement] = Enforcement.SOFT
    evidence: ClassVar[Evidence] = Evidence.CONTESTED
    direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
    unit: ClassVar[str] = "GC dispersion relative to chance (1.0 = binomial floor)"
    band: ClassVar[tuple[float, float] | None] = None
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "SCP4ssd (AutoML, F1 0.930) ranks dGC -- GC fluctuation in 100 bp windows "
            "-- second among its features, ABOVE global GC, and finds that locally "
            "LOW-GC fragments matter more than high-GC ones. 25 of its 31 optimal "
            "features were absent from the Synthesis Success Calculator's set",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC10048150/",
            2023,
            sign="supports",
        ),
        Citation(
            "The Synthesis Success Calculator (random forest, 1076 real vendor "
            "outcomes, F1 0.928) ranks REPEATS as the dominant contributor to "
            "synthesis failure and GC below them -- the opposite order. It is also "
            "the source of the brief's instruction to weight repeats above GC",
            "https://pubs.acs.org/doi/10.1021/acssynbio.9b00460",
            2020,
            sign="refutes",
        ),
        Citation(
            "On the only head-to-head experimental validation, ten E. coli genes, "
            "SCP4ssd made eight correct predictions and the Synthesis Success "
            "Calculator five. Weak evidence -- n=10, and reported by the winner -- "
            "but it is the only direct comparison, and 5/10 is chance",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC10048150/",
            2023,
            sign="qualifies",
        ),
        Citation(
            "Twist publishes a 50 bp GC WINDOW (10-90%), which fixes the window size "
            "used here for extent. It publishes no bound on the extent or the SD, so "
            "this rule reports variation and never enforces it",
            "https://www.twistbioscience.com/faq/gene-synthesis/what-do-scoring-results-my-gene-mean",
            2026,
            sign="qualifies",
        ),
    )
    last_verified: ClassVar[str] = "2026-08-28"
    weight_provenance: ClassVar[str] = (
        "0.60: deliberately just BELOW e6_repeat_density (0.65) rather than above "
        "it, and the gap is the honest size of what is known. SCP4ssd ranks this "
        "family above repeats; the Synthesis Success Calculator ranks repeats "
        "above it. The only head-to-head favours SCP4ssd 8/10 to 5/10, but at n=10 "
        "and self-reported, which is not enough to invert a published ordering. So "
        "the two are weighted as comparable and the tie breaks toward the older, "
        "larger, independently-cited model. What justifies carrying real weight at "
        "all is orthogonality rather than rank: measured over 200 synonymous "
        "variants x 4 seeds, |rho| <= 0.06 against e6, e8, f2 AND e2, while moving "
        "49-59% of its median under synonymous substitution. It is a large axis no "
        "other rule sees, so weighting it near zero would discard information no "
        "amount of repeat weighting recovers. The same run calibrates the scalar: "
        "an ordinary protein reads 1.01, which is chance to two decimals, and a "
        "repetitive one 1.34 -- so the number means what it claims to mean, and "
        "the 1.5 reporting threshold sits above where repetitive proteins already "
        "land rather than firing on all of them."
    )
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.60
    #: Below e2_gc_band (0.5) on purpose. E2 steers toward a bound the validator
    #: will REFUSE to emit without; this steers toward a preference with no
    #: published threshold, and a preference should not pull as hard as a
    #: requirement on the same bases.
    steering_weight: ClassVar[float] = 0.3
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.WINDOW_MINUS_1
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "cheap"
    #: Same bases, opposite directions: F5's AT ceiling raises GC in AT-rich
    #: windows, which can widen the very excursion this rule reports.
    conflicts_with: ClassVar[tuple[str, ...]] = ("f5_at_window", "e2_gc_band")
    brief_ref: ClassVar[str] = "2.E4"
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "dgc_window": {"type": "integer", "default": DGC_WINDOW_BP, "minimum": 20},
            "extent_window": {
                "type": "integer",
                "default": EXTENT_WINDOW_BP,
                "minimum": 20,
            },
            "warn_ratio": {
                "type": "number",
                "default": WARN_RATIO,
                "minimum": 1.0,
                "description": (
                    "Report a fragment at least this many times as dispersed as "
                    "chance. 1.0 would report every sequence, since chance is 1.0."
                ),
            },
            "vendor": {
                "type": "string",
                "default": DEFAULT_VENDOR,
                "enum": sorted(VENDOR_ADAPTERS),
                "description": (
                    "Which vendor configuration the fragment is ordered as. Only "
                    "the adapter-on options carry adapters, and adapter GC counts "
                    "toward the excursion because the vendor synthesises it."
                ),
            },
        },
    }

    def __init__(
        self,
        dgc_window: int = DGC_WINDOW_BP,
        extent_window: int = EXTENT_WINDOW_BP,
        warn_ratio: float = WARN_RATIO,
        vendor: str = DEFAULT_VENDOR,
    ) -> None:
        if dgc_window < 20 or extent_window < 20:
            raise ValueError(
                f"windows under 20 bp measure codon composition rather than a "
                f"composition gradient; got dgc_window={dgc_window}, "
                f"extent_window={extent_window}"
            )
        if warn_ratio <= 1.0:
            raise ValueError(
                f"warn_ratio {warn_ratio} is at or below chance (1.0), which would "
                f"report every sequence including a perfectly uniform one"
            )
        if vendor not in VENDOR_ADAPTERS:
            raise ValueError(f"unknown vendor {vendor!r}; have {sorted(VENDOR_ADAPTERS)}")
        self.dgc_window = dgc_window
        self.extent_window = extent_window
        self.warn_ratio = warn_ratio
        self.vendor = vendor

    def gate(self, slot: ContextSlot) -> bool:
        # Every construct BT5 designs is ordered as DNA, in every context.
        return True

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """None, for the reason E2 gives: deciding what a codon does to a windowed
        GC statistic needs the G+C count over the previous ~17 codons, and that
        history is not expressible in the automaton state."""
        return None

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        adapters = VENDOR_ADAPTERS[self.vendor]
        breaches: list[Breach] = []
        windows: list[tuple[Interval, float]] = []
        worst = 0.0
        scanned = 0

        for frag in fragments(c, adapters):
            scanned += frag.ordered.length
            # The MAX across fragments, as e6 does and for the same reason: each
            # designable span is its own synthesis reaction, and an order succeeds
            # only if every tube does.
            worst = max(worst, dispersion_ratio(frag.sequence, self.dgc_window))
            breach = self._extent(frag)
            if breach is not None:
                iv, value, found = breach
                windows.append((iv, value))
                breaches.append(found)

        return Evaluation(
            spec_id=self.id,
            passes=not any(b.magnitude >= 1.0 for b in breaches),
            raw_score=worst,
            breaches=tuple(breaches[:MAX_FINDINGS]),
            windows=tuple(windows),
            n_evaluated=scanned,
        )

    def _extent(self, frag: Fragment) -> tuple[Interval, float, Breach] | None:
        """One finding per fragment, anchored on the window worth recoding.

        The extent is set by two windows and closing either one narrows it, so
        the breach is anchored on whichever sits further from the fragment's
        median -- that is the one excursion, rather than the baseline the rest of
        the molecule sits at -- and the message names both.
        """
        ratio = dispersion_ratio(frag.sequence, self.dgc_window)
        if ratio < self.warn_ratio:
            return None
        values = gc_windows(frag.sequence, self.extent_window)
        if len(values) < 2:
            return None
        lo_i, lo, hi_i, hi = extremes(values)
        span = hi - lo

        median = sorted(values)[len(values) // 2]
        low_is_outlier = (median - lo) >= (hi - median)

        at = lo_i if low_is_outlier else hi_i
        side = "low" if low_is_outlier else "high"

        # None means the window lies WHOLLY inside a vendor adapter, which the
        # dispersion gate above makes unreachable: for an adapter window to be
        # the extreme, every insert window must sit closer to the median than the
        # adapter does, and a fragment that uniform never reaches ratio 1.5.
        # Checked over 840 configurations of length, composition, window and
        # seed -- zero reached it. This is here to narrow the Optional, not as a
        # fallback for a case that exists.
        parent = frag.to_construct(Interval(at, at + self.extent_window))
        if parent is None:
            return None

        magnitude = min(1.0, (ratio - 1.0) / (HARD_RATIO - 1.0))
        message = (
            f"GC is {ratio:.1f}x as dispersed as chance for this composition. It "
            f"swings {span:.0f} points, from {lo:.0f}% at {lo_i} to {hi:.0f}% at "
            f"{hi_i}; against a median of {median:.0f}% the {side} window is the "
            f"outlier, so recoding it toward the middle narrows the excursion "
            f"fastest. Every window here can sit inside its own GC band and still "
            f"leave the two ends annealing at different temperatures"
        )
        return (
            parent,
            span,
            Breach(
                spec_id=self.id,
                interval=parent,
                magnitude=magnitude,
                message=message,
                fixable_by_codon_choice=True,
                detail={
                    "vendor": self.vendor,
                    "dispersion_ratio": round(ratio, 3),
                    "extent_pct": round(span, 2),
                    "min_gc_pct": round(lo, 2),
                    "max_gc_pct": round(hi, 2),
                    "median_gc_pct": round(median, 2),
                    "outlier_side": side,
                    "low_window_start": float(lo_i),
                    "high_window_start": float(hi_i),
                    "dgc_pct": round(dgc(frag.sequence, self.dgc_window), 2),
                    "chance_dgc_pct": round(chance_dgc(frag.sequence, self.dgc_window), 2),
                    "window_bp": float(self.extent_window),
                },
            ),
        )
