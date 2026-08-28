"""E6 -- how repetitive the fragment is overall, as one number.

Every other repeat rule in this catalog reports PAIRS: these two spans are
identical, here is the geometry. This one reports a DENSITY, and the difference
is not cosmetic. A fragment carrying twenty 9 bp repeats and no 20 bp one passes
E5, F1 and F2 and still fails at the vendor -- which is exactly the case the
only published synthesis-success model was trained to catch. Repetitive 9-mers
per 100 bp is its highest-Gini feature, over a random forest on 1,076 real
vendor outcomes at F1 0.928 (https://pubs.acs.org/doi/10.1021/acssynbio.9b00460).

**Why there is no hard threshold on this number, and why that is not timidity.**
The metric is length-dependent, and steeply. Chance collisions grow as n^2 while
the per-100-bp denominator grows as n, so the expected value rises linearly in
length at 100n / (2 * 4^9). Measured on uniform random sequence, which tracks
that closed form:

    300 bp  0.075    900 bp  0.203    1500 bp  0.307    3000 bp  0.611
                                                        5000 bp  0.984

Identical composition, an order of magnitude of spread. A fixed cutoff would
therefore fail every long fragment and pass every short one regardless of design
quality. The cited model uses the feature alongside length, so the model absorbs
this; a threshold cannot. BT5 already owns the right correction -- every
objective is normalised to a percentile against 200-500 random synonymous
variants OF THIS PROTEIN, which are the same length by construction -- so this
rule is SOFT, reports its native units, and leaves the comparison to the null.
A length correction invented here would be a second, worse null competing with
the real one.

**The window flag is measured, not chosen.** Fraction of a 100 bp window's
9-mers that occur again anywhere in the fragment:

    uniform random        p50 0.000, p99 0.043, max 0.120
    codon-biased CDS      p50 0.076, p90 0.208
    (GGGGS)x20 linker     1.000 everywhere

0.25 clears the p99 of both baselines several times over and sits far under the
linker, which is the case this rule exists to catch.

**A known limitation, in the feature itself.** "Repetitive 9-mers per 100 bp"
counts DISTINCT 9-mers that recur, which is the cited feature and is kept
faithful to it. That count rises with the PERIOD of a tandem array: 6 copies of
CAT contain 3 recurring 9-mers while the same 18 bp as 3 copies of CACCAT
contain 4, so the metric ranks the less repetitive encoding marginally worse.
It mis-orders only short-period tandem arrays, which are E1's (homopolymers) and
E7's (short tandem repeats) to band, and both encodings are far under this
rule's resolution anyway -- 18 bp of repetition is not a density finding. For
dispersed repeats, the case this rule is for, distinct-count and position-count
agree. The per-window fraction reported alongside is position-based and does not
have the artefact.

**What this rule deliberately does not do.** The cited model's other top feature,
longest repetitive sequence, is E5's length band, reported there with the pair
geometry and the duplex Tm. Computing it here as well would put the same base
pairs in the panel twice under two names. Locally, per window, the longest run
of repeated 9-mers IS reported here, because "this window is 60% repeated" and
"these two spans are identical" are different statements about different objects.
"""

from __future__ import annotations

from collections import Counter
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
from bt5.core.types import Construct, Interval
from bt5.rules.fragment import VENDOR_ADAPTERS, Fragment, fragments

#: The model's k. Not tunable downward without leaving the cited feature behind.
KMER_BP = 9
#: The model's denominator, and the window this rule localises on.
WINDOW_BP = 100
#: Window step. Fine enough to place a repetitive region, coarse enough that a
#: 1.5 kb fragment costs 70 windows rather than 1,400.
STEP_BP = 20
#: Above the 99th percentile of ordinary coding sequence (0.076) with margin.
WINDOW_FLAG = 0.25
MAX_FINDINGS = 200


def repetitive_kmers_per_100bp(seq: str, k: int = KMER_BP) -> float:
    """The cited feature: distinct k-mers occurring more than once, per 100 bp."""
    if len(seq) < k:
        return 0.0
    counts = Counter(seq[i : i + k] for i in range(len(seq) - k + 1))
    recurring = sum(1 for v in counts.values() if v >= 2)
    return recurring / (len(seq) / 100.0)


def _duplicate_flags(seq: str, k: int) -> list[bool]:
    """Per k-mer start position: does this k-mer occur again in the fragment?

    Anywhere in the fragment, not just in the window. A 9-mer repeated 500 bp
    away is the liability; one that happens to sit twice inside one window is
    the same liability seen closer up.
    """
    if len(seq) < k:
        return []
    counts = Counter(seq[i : i + k] for i in range(len(seq) - k + 1))
    return [counts[seq[i : i + k]] >= 2 for i in range(len(seq) - k + 1)]


def _longest_run(flags: list[bool]) -> int:
    """Longest consecutive stretch of repeated k-mer positions."""
    best = run = 0
    for flag in flags:
        run = run + 1 if flag else 0
        best = max(best, run)
    return best


@register
class RepeatDensity:
    id: ClassVar[str] = "e6_repeat_density"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "Repetitive 9-mer density in the synthesized fragment"
    enforcement: ClassVar[Enforcement] = Enforcement.SOFT
    evidence: ClassVar[Evidence] = Evidence.EVIDENCE_BACKED
    direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
    unit: ClassVar[str] = "repetitive 9-mers per 100 bp"
    band: ClassVar[tuple[float, float] | None] = None
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "Repetitive 9-mers per 100 bp and longest repetitive sequence are the two "
            "highest-importance features of the Synthesis Success Calculator -- a "
            "random forest over 1,076 real vendor outcomes, F1 0.928 -- ranking above "
            "every GC feature",
            "https://pubs.acs.org/doi/10.1021/acssynbio.9b00460",
            2020,
            sign="supports",
        ),
        Citation(
            "The AAV liability is repetitiveness, not GC: a designed stuffer at GC "
            "43.5-44.8% cut yield up to 68% while a LOWER-GC natural stuffer of "
            "identical length cost neither yield nor bioactivity",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC12207685/",
            2025,
            sign="supports",
        ),
        Citation(
            "Twist names repeated sequence as a complexity trigger in its own right, "
            "separately from length and GC -- the vendor-side confirmation that "
            "repetitiveness is priced and rejected on independently of the pair rules",
            "https://www.twistbioscience.com/faq/gene-synthesis",
            2026,
            sign="supports",
        ),
    )
    last_verified: ClassVar[str] = "2026-08-28"
    weight_provenance: ClassVar[str] = (
        "0.65: above f2_near_perfect_repeats (0.60), whose 90%/40bp threshold is a "
        "convention, because this is the single highest-Gini feature of the only "
        "synthesis model trained on real vendor outcomes. Below d4_internal_polya "
        "(0.70) on purpose -- feature importance RANKS predictors, it does not "
        "measure effect size, and it should not out-weight a rule backed by a "
        "measured 8-9x functional titer loss. The brief's instruction to weight "
        "repeats above GC is satisfied structurally rather than numerically: BT5's "
        "GC rules are hard constraints and carry no objective weight at all."
    )
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.65
    #: Deliberately modest. This is a smooth scalar and the DP could chase it
    #: hard, but f1 already steers repeats at 1.0 and e5 at 0.7 against largely
    #: the same bases; stacking a third full-strength repeat term would spend the
    #: sequence's freedom on one family. The increment here is the density
    #: signal specifically -- the many-short-repeats case the pair rules miss.
    steering_weight: ClassVar[float] = 0.4
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.WINDOW_MINUS_1
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "cheap"
    conflicts_with: ClassVar[tuple[str, ...]] = ()
    brief_ref: ClassVar[str] = "2.E6"
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "k": {"type": "integer", "default": KMER_BP, "minimum": 6},
            "window": {"type": "integer", "default": WINDOW_BP, "minimum": 20},
            "window_flag": {
                "type": "number",
                "default": WINDOW_FLAG,
                "minimum": 0.0,
                "maximum": 1.0,
            },
            "vendor": {
                "type": "string",
                "default": "twist_gene_fragment",
                "enum": sorted(VENDOR_ADAPTERS),
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
        k: int = KMER_BP,
        window: int = WINDOW_BP,
        window_flag: float = WINDOW_FLAG,
        vendor: str = "twist_gene_fragment",
    ) -> None:
        if k < 6:
            raise ValueError(
                f"k {k} is short enough to recur by chance everywhere: a 6-mer has "
                f"~0.4 expected occurrences per 1.5 kb even in random sequence"
            )
        if window < 2 * k:
            raise ValueError(f"window {window} must hold at least two {k}-mers")
        if not 0.0 < window_flag <= 1.0:
            raise ValueError(f"window_flag must be a fraction in (0, 1], got {window_flag}")
        if vendor not in VENDOR_ADAPTERS:
            raise ValueError(f"unknown vendor {vendor!r}; have {sorted(VENDOR_ADAPTERS)}")
        self.k = k
        self.window = window
        self.window_flag = window_flag
        self.vendor = vendor

    def gate(self, slot: ContextSlot) -> bool:
        return True

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """None. A k-mer's multiplicity depends on the whole fragment."""
        return None

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        adapters = VENDOR_ADAPTERS[self.vendor]
        breaches: list[Breach] = []
        windows: list[tuple[Interval, float]] = []
        worst_density = 0.0
        scanned = 0

        for frag in fragments(c, adapters):
            scanned += frag.ordered.length
            # The MAX across fragments, not the mean. Each designable span is a
            # separate synthesis reaction, and an order succeeds only if every
            # tube does, so the worst fragment is the decision-relevant number.
            worst_density = max(worst_density, repetitive_kmers_per_100bp(frag.sequence, self.k))
            local, flagged = self._windows(frag)
            windows.extend(local)
            breaches.extend(flagged)

        return Evaluation(
            spec_id=self.id,
            passes=not breaches,
            raw_score=worst_density,
            breaches=tuple(breaches[:MAX_FINDINGS]),
            windows=tuple(windows),
            n_evaluated=scanned,
        )

    def _windows(self, frag: Fragment) -> tuple[list[tuple[Interval, float]], list[Breach]]:
        flags = _duplicate_flags(frag.sequence, self.k)
        if not flags:
            return [], []
        span = self.window - self.k + 1
        starts = range(0, max(1, len(flags) - span + 1), STEP_BP)

        measured: list[tuple[Interval, float, int]] = []
        for start in starts:
            seg = flags[start : start + span]
            if not seg:
                continue
            iv = Interval(start, min(start + self.window, len(frag.sequence)))
            measured.append((iv, sum(seg) / len(seg), _longest_run(seg)))

        windows = [
            (mapped, frac)
            for iv, frac, _ in measured
            if (mapped := frag.to_construct(iv)) is not None
        ]
        return windows, self._merge(frag, measured)

    def _merge(self, frag: Fragment, measured: list[tuple[Interval, float, int]]) -> list[Breach]:
        """One breach per contiguous repetitive REGION, not per window.

        The windows overlap by design, so a 300 bp repetitive stretch trips
        fifteen of them. Reporting each is the same failure E1 had with
        homopolymers and F2 had with adjacent seeds: one physical problem, many
        findings, and a conflict panel that reads as fifteen defects.
        """
        out: list[Breach] = []
        current: tuple[int, int, float, int] | None = None

        for iv, frac, run in measured:
            if frac < self.window_flag:
                continue
            if current is not None and iv.start <= current[1]:
                current = (
                    current[0],
                    max(current[1], iv.end),
                    max(current[2], frac),
                    max(current[3], run),
                )
                continue
            if current is not None:
                out.append(self._breach(frag, *current))
            current = (iv.start, iv.end, frac, run)

        if current is not None:
            out.append(self._breach(frag, *current))
        return out

    def _breach(self, frag: Fragment, start: int, end: int, frac: float, run: int) -> Breach:
        iv = Interval(start, end)
        mapped = frag.to_construct(iv)
        target = mapped if mapped is not None else iv
        return Breach(
            spec_id=self.id,
            interval=target,
            magnitude=frac,
            message=(
                # The fraction is the PEAK over the windows merged into this
                # region, not the region's own average. Attributing the peak to
                # the whole span would overstate every finding that merged.
                f"repetitive region of {end - start} bp at {target.start}: at its "
                f"worst, {frac:.0%} of the {self.k}-mers in a {self.window} bp window "
                f"occur again in the ordered fragment, with a repeated stretch of "
                f"{run + self.k - 1} bp. Repetitive {self.k}-mers per 100 bp is the "
                f"highest-importance feature of the only synthesis-success model "
                f"trained on real vendor outcomes"
            ),
            fixable_by_codon_choice=frag.construct.overlaps_editable(iv),
            detail={
                "peak_repeated_fraction": round(frac, 3),
                "longest_repeated_bp": float(run + self.k - 1),
                "k": float(self.k),
            },
        )
