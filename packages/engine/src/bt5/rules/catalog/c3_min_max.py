"""C3 -- %MinMax, the windowed rare/common codon profile, as a soft band with a
ceiling at the metric's own neutral point.

%MinMax (Clarke & Clark 2008) slides a window along the CDS and asks, per window,
how far the codons actually chosen sit between the host's AVERAGE synonymous usage
and its most extreme. The paper's own scale is the whole argument for this rule's
band: +100 is "a sequence window encoded using only the most common codons", 0 is
"codon usage equal to the mean of all possible codon choices", and -100 is "a
sequence window encoded using only the most rare codons".

Two decisions carry this file.

**The ceiling is 0, and it is the operative half.** brief.md:73 heads section 2.C
"all S, soft bands, never maximized", and 0 is not a threshold anybody picked -- it
is the point Clarke & Clark define as average host usage. A ceiling there says
exactly one thing: do not push the CDS above what the host already does. That is
the failure mode with the evidence behind it. Ranaghan 2021 benchmarked nine
commercial and academic optimizers and found "a roughly equivalent chance that an
algorithm-optimized CDS will increase or diminish recombinant yields"; Welch 2009's
deliberately high-usage control expressed at ~15% of their best variant. Neither
result says a below-average CDS is better -- they say chasing the top codon is not
better, which is a ceiling, not a target.

**The floor is -100 and is non-binding by construction.** It is the metric's
definitional minimum, so it can never be breached, and that is deliberate rather
than an oversight: the evidence affirmatively says the low side should NOT be
penalised. Clarke & Clark's own finding is that rare codons CLUSTER non-randomly --
their title -- and brief.md:84 (C8) asks for >=80% of native rare-codon clusters to
be RETAINED. A floor that punished %Min excursions would fight both. Declaring the
band one-sided and saying so beats inventing a floor to look symmetric.

`Direction.BAND` rather than `LOWER_IS_BETTER` for the same reason C1 is a band: a
monotone "minimise %Max" objective has its optimum at an all-rare-codon sequence,
which is a catastrophic design and which nothing in the evidence base recommends.
A band cannot run away.

**Scope: this rule computes nothing in this build, and says so per host.** %MinMax
is defined on raw codon usage FREQUENCIES -- X_ij, and the per-family X_avg,i,
X_max,i, X_min,i. `data/codon_usage/` ships exactly one table, Sharp & Li's
relative adaptiveness w-index, and w is not a frequency: w_ij = (count_ij + 0.5) /
max synonymous (count + 0.5), so each family is renormalised to its own peak and
the peak is then discarded. Per family the differences survive the rescaling --
w_ij - w_avg,i = (X_ij - X_avg,i) / K_i -- but %MinMax sums those differences
ACROSS families, and recovering the per-family K_i that would make the sum
comparable is exactly what normalising to 1.0 threw away. So %MinMax computed on w
is a different statistic wearing this one's name and citation. Substituting it
because it is the table on disk is the same move as scoring a mammalian CDS against
the E. coli w-index, and `docs/decisions/2026-09-01-c1-cai-soft-band.md` already
refused that one: "the fallback would always succeed and always be wrong -- a
plausible-looking number measuring nothing".

So `MINMAX_REFERENCE_SET` is empty and every host reports the objective unavailable
with the reason. The preset machinery gets what it was actually waiting for: 2.C3
stops reading as an objective NOBODY IMPLEMENTED and starts reading as one whose
reference data this build does not ship, which is a materially different statement
to a user weighing a ranking. The arithmetic below is real, tested and complete; it
needs one entry in the map on the day a frequency table lands.

The gate is `slot.role != "propagation"`, and for C1's reason: a lentiviral job
propagates in E. coli and expresses in HEK293, so a rule keyed on "is any host
E. coli" would compute a confident number for the one host with a table and report
it as the objective for a protein made somewhere else.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import ClassVar

from bt5.core.context import ContextSlot, DesignContext, HostId
from bt5.core.registry import register
from bt5.core.services import GeneticCode, Services
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
from bt5.core.types import Construct, Interval, Strand

#: Clarke & Clark 2008: "The resulting values are typically averaged over an
#: 18-codon window", and "All results shown in Figures 2-4 used a window size of
#: 18." brief.md:79 gives the range the literature spans ("z = 17-18 codons
#: (CHARMING default 10)"); docs/PLAN.md:661 resolves it by decision -- "Pin 18
#: (CodonTransformer's value; CHARMING defaults to 10, other sources say 17).
#: Expose in `param_schema`, print in the report." 18 is therefore both the pinned
#: value and the one its primary source measured on.
WINDOW_CODONS = 18

#: The band. Both edges come from the metric's definition rather than from a
#: threshold anyone chose -- see the module docstring.
BAND_LO = -100.0
BAND_HI = 0.0

#: HostId -> the reference-set key a per-codon usage FREQUENCY table would be
#: looked up under, mirroring `c1_cai.CAI_REFERENCE_SET`.
#:
#: Deliberately empty. No shipped reference set carries frequencies: the one file
#: in `data/codon_usage/` is a relative-adaptiveness w-index, which %MinMax cannot
#: be computed from (module docstring). This is the whole of the wiring a frequency
#: table needs -- add the host and its reference-set stem here and the rule is live
#: for that host, with no other change to this file.
MINMAX_REFERENCE_SET: Mapping[HostId, str] = {}

#: What `TableProvider` implementations raise when a host's table is absent.
#: `FileTableProvider` raises `FileNotFoundError` (an `OSError`); `NotImplementedError`
#: is what it raises for the kinds that are a later lane.
_MISSING_TABLE = (OSError, LookupError, NotImplementedError, ValueError)


def family_statistics(
    code: GeneticCode, freq: Mapping[str, float]
) -> Mapping[str, tuple[float, float, float, float]]:
    """codon -> (X_ij, X_avg,i, X_max,i, X_min,i), the four terms %MinMax needs.

    The last three are properties of the family the codon belongs to, so every
    codon of one amino acid shares them; the first is the codon's own usage
    frequency. Positions carrying no codon-choice information are simply absent
    from the result and drop out of every window.

    Excluded, and read from the INJECTED table rather than hard-coded: stops, and
    any amino acid with a single non-stop codon under this table. NCBI table 4
    makes TGA a second Trp codon, so Trp does carry information there, and tables
    27/28 make TGA both Trp and a stop (CLAUDE.md 3.1, 3.2). A family the reference
    set does not fully cover is dropped whole rather than averaged over the codons
    it happens to list -- a partial family's mean is not the family's mean, and the
    difference would ride silently into every window overlapping it.
    """
    stats: dict[str, tuple[float, float, float, float]] = {}
    seen: set[str] = set()
    for codon in freq:
        if code.is_stop(codon):
            continue
        aa = code.translate(codon)
        if not aa or aa == "*" or aa in seen:
            continue
        seen.add(aa)
        try:
            family = code.synonymous_codons(aa)
        except (ValueError, KeyError):
            continue
        if len(family) < 2:
            continue
        values = [freq.get(member) for member in family]
        if any(v is None for v in values):
            continue
        present = [float(v) for v in values if v is not None]
        mean, most, least = sum(present) / len(present), max(present), min(present)
        for member in family:
            stats[member] = (float(freq[member]), mean, most, least)
    return stats


def min_max_profile(
    codons: Sequence[str],
    stats: Mapping[str, tuple[float, float, float, float]],
    window: int = WINDOW_CODONS,
) -> list[float]:
    """The signed %MinMax value centred on every codon position.

    Clarke & Clark 2008, verbatim: "%Max = 100 * sum(X_ij - X_avg,i) /
    sum(X_max,i - X_avg,i)" over the window's codons when the window's usage is
    above average, and the analogous %Min below it. Both of the paper's equations
    return a positive number by definition; the profile is signed here (+%Max,
    -%Min) because that is the form the paper plots and reads -- "clusters of
    predominantly rare codons appear as negative (%Min) peaks" -- and an unsigned
    profile could not distinguish the two failure directions this rule's band
    treats so differently.

    Windows SHRINK at the termini rather than dropping positions (brief.md:79), so
    the profile is exactly as long as `codons` and the first and last codons are
    scored rather than silently omitted. A window whose informative positions are
    all single-codon families scores 0.0 -- correct rather than undefined: nothing
    in it could have been chosen differently, so it is neither common nor rare.
    """
    if window < 1:
        raise ValueError(f"window must be at least 1 codon, got {window}")
    n = len(codons)
    behind = (window - 1) // 2
    ahead = window // 2
    profile: list[float] = []
    for centre in range(n):
        lo = max(0, centre - behind)
        hi = min(n, centre + ahead + 1)
        actual = average = peak = trough = 0.0
        for codon in codons[lo:hi]:
            row = stats.get(codon)
            if row is None:
                continue
            used, mean, most, least = row
            actual += used
            average += mean
            peak += most
            trough += least
        if actual >= average:
            span = peak - average
            profile.append(100.0 * (actual - average) / span if span > 0.0 else 0.0)
        else:
            span = average - trough
            profile.append(-100.0 * (average - actual) / span if span > 0.0 else 0.0)
    return profile


@register
class MinMax:
    id: ClassVar[str] = "c3_min_max"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "%MinMax windowed codon usage profile, capped at host average"
    enforcement: ClassVar[Enforcement] = Enforcement.SOFT
    #: brief.md:79 grades C3 "B" -- the brief's own legend (brief.md:46) reads
    #: "B = replicated but contested or single-lab". CONTESTED is the badge that
    #: says exactly that; EVIDENCE_BACKED would upgrade the brief's own grading.
    evidence: ClassVar[Evidence] = Evidence.CONTESTED
    direction: ClassVar[Direction] = Direction.BAND
    unit: ClassVar[str] = "%MinMax (+100 all-common, 0 host-average, -100 all-rare)"
    band: ClassVar[tuple[float, float] | None] = (BAND_LO, BAND_HI)
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "Clarke & Clark 2008, the metric itself and the source of every number "
            "in this file: %Max = 100*sum(X_ij - X_avg,i)/sum(X_max,i - X_avg,i) "
            "with %Min analogous below average; -100/0/+100 are 'only the most "
            "rare codons' / 'codon usage equal to the mean of all possible codon "
            "choices' / 'only the most common codons'; and 'The resulting values "
            "are typically averaged over an 18-codon window', with 'All results "
            "shown in Figures 2-4 used a window size of 18'. Defined on raw usage "
            "FREQUENCIES, which is why this rule refuses to compute it from a "
            "relative-adaptiveness table",
            "https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0003412",
            2008,
            sign="supports",
        ),
        Citation(
            "Clarke & Clark's actual finding -- rare codons cluster non-randomly, "
            "enriched at the 5' and 3' termini of E. coli genes. The direct reason "
            "the band's floor is non-binding: %Min excursions are the structure the "
            "paper reports, not a defect to optimise away",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC2565806/",
            2008,
            sign="qualifies",
        ),
        Citation(
            "Ranaghan 2021 benchmarked nine commercial and academic optimizers and "
            "found 'a roughly equivalent chance that an algorithm-optimized CDS will "
            "increase or diminish recombinant yields', with three tools "
            "non-deterministic. The evidence for a CEILING specifically: pushing "
            "codon usage up is a coin flip, so it is not a target",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC7893858/",
            2021,
            sign="supports",
        ),
        Citation(
            "Welch 2009: their deliberately high-usage control expressed at ~15% of "
            "the best variant, and the paper states flatly that 'CAI has no value in "
            "predicting gene expression'. Max-usage is a MEASURED failure mode, "
            "which is what makes 0 a ceiling rather than a midpoint to aim at",
            "https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0007002",
            2009,
            sign="supports",
        ),
        Citation(
            "Cambray 2018: all computable design features together explain 5-31% "
            "(mean ~14%) of protein-level variance. The ceiling every "
            "codon-composition objective sits under, and why C3 reports a band "
            "position and never a predicted expression level",
            "https://www.nature.com/articles/nbt.4238",
            2018,
            sign="qualifies",
        ),
    )
    last_verified: ClassVar[str] = "2026-09-02"
    weight_provenance: ClassVar[str] = (
        "0.3, matching what all three shipped presets already assign 2.C3, and set "
        "against C1's 0.2 rather than derived independently: %MinMax and CAI measure "
        "the same underlying quantity -- how far synonymous choice has been pushed "
        "toward the host's common codons -- and %MinMax is the better-behaved of the "
        "two, because it is windowed and so localises WHERE the pushing happened "
        "instead of collapsing the ORF to one geometric mean. Half a notch above C1 "
        "reproduces that ordering without claiming the evidence is stronger; it is "
        "the same evidence, graded B by the brief, and it sits far below B1's 1.0 "
        "(r = 0.66 over 154 measured variants) because nothing here has an effect "
        "size of that kind. steering_weight is 0.0 and must stay there for C1's "
        "reason with the sign reversed: a Lagrangian term on codon usage steers the "
        "Tier-A DP monotonically, and toward this band's interior means toward RARE "
        "codons, which is a worse design than the one it was correcting."
    )
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.3
    #: Zero, and structurally so -- see weight_provenance and lattice_terms.
    steering_weight: ClassVar[float] = 0.0
    #: WHOLE_SCOPE: the breached quantity is the mean of the profile over the whole
    #: coding scope, and no single window is "the" reason a mean left the band. The
    #: per-window profile still travels in `Evaluation.windows`, which is where the
    #: report reads it from (docs/PLAN.md:661 -- "print in the report").
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.WHOLE_SCOPE
    #: SINGLE_PASS, and safe here in the way it is not for splice donors
    #: (CLAUDE.md 3.6). Repair lowers %Max by swapping common codons toward their
    #: family average; each swap lowers that codon's contribution to every window
    #: containing it and raises none of them, so a pass cannot manufacture a new
    #: above-ceiling window the way removing one donor creates another.
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "cheap"
    #: None declared. The obvious candidates are not conflicts: C1's band and this
    #: ceiling both push AWAY from max-usage, and the repeat rules (E5/E6/F1) are
    #: what max-usage manufactures, so C3 relieves them rather than fighting them.
    #: C8 (rare-cluster preservation) would conflict with a floor -- which is one
    #: more reason this band does not have an operative one.
    conflicts_with: ClassVar[tuple[str, ...]] = ()
    brief_ref: ClassVar[str] = "2.C3"
    #: Arithmetic over a usage table, not a folding energy: no engine's parameter
    #: set is involved and none needs to be matched.
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "window": {
                "type": "integer",
                "default": WINDOW_CODONS,
                "minimum": 1,
                "description": (
                    "Sliding window in CODONS. 18 is Clarke & Clark's own value and "
                    "the value pinned by decision (PLAN.md:661); the literature "
                    "spans 10 (CHARMING) to 18. The paper reports window sizes from "
                    "5 to 30 giving similar distributions with more noise at the "
                    "small end, so this is a resolution knob, not a threshold."
                ),
            },
            "min_max_max": {
                "type": "number",
                "default": BAND_HI,
                "minimum": -100.0,
                "maximum": 100.0,
                "description": (
                    "Ceiling of the target band, and the operative half. 0 is the "
                    "metric's own neutral point -- codon usage equal to the mean of "
                    "all possible synonymous choices -- so the default says 'do not "
                    "push above what the host already does' rather than encoding a "
                    "threshold anybody picked."
                ),
            },
            "min_max_min": {
                "type": "number",
                "default": BAND_LO,
                "minimum": -100.0,
                "maximum": 100.0,
                "description": (
                    "Floor of the target band. -100 is the metric's definitional "
                    "minimum, so the default floor can never be breached: rare-codon "
                    "clusters are structure this rule declines to penalise. Raise it "
                    "only with evidence that a specific %Min depth is harmful."
                ),
            },
        },
    }

    def __init__(
        self,
        window: int = WINDOW_CODONS,
        min_max_min: float = BAND_LO,
        min_max_max: float = BAND_HI,
    ) -> None:
        if window < 1:
            raise ValueError(f"window must be at least 1 codon, got {window}")
        if not -100.0 <= min_max_min < min_max_max <= 100.0:
            raise ValueError(
                f"%MinMax band must satisfy -100 <= min_max_min {min_max_min} < "
                f"min_max_max {min_max_max} <= 100"
            )
        self.window = window
        self.min_max_min = min_max_min
        self.min_max_max = min_max_max

    def gate(self, slot: ContextSlot) -> bool:
        """Every slot whose protein is actually translated.

        Not `host`, for C1's reason: the E. coli slot of a lentiviral job is a
        PROPAGATION slot, so a host-keyed rule would report a confident number for
        the one host likely to have a table and attribute it to a protein made in
        HEK293.
        """
        return slot.role != "propagation"

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """Deliberately None, for C1's argument in both directions.

        A `codon_weights` term would steer the Tier-A DP monotonically, and this
        rule's band has an operative ceiling, so the steer would run toward RARE
        codons -- trading the failure mode the brief warns about for a worse one it
        does not. The band is not expressible in the automaton regardless: the
        automaton decides from a bounded codon suffix, and while one %MinMax window
        is bounded, the scored quantity is the MEAN over every window in the ORF, so
        no bounded state can know whether the next codon takes it above the ceiling.
        `forbidden` is empty because no individual codon is illegal.
        """
        return None

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        slots = [s for s in ctx.active_slots if self.gate(s)]
        if not slots:
            return self._unavailable(
                c,
                "no translating slot in this context (propagation-only), so there is "
                "no host whose codon usage the CDS should be measured against",
            )

        editable = sorted(c.editable)
        if not editable:
            return self._unavailable(c, "no designable CDS to compute a %MinMax profile over")

        per_slot: list[tuple[ContextSlot, str, float, list[float], int]] = []
        for slot in slots:
            reference_set = MINMAX_REFERENCE_SET.get(slot.host)
            if reference_set is None:
                return self._unavailable(
                    c,
                    f"no codon usage FREQUENCY reference set for host {slot.host}. "
                    f"%MinMax is defined on raw usage frequencies, and this build "
                    f"ships one codon usage table -- Sharp & Li's E. coli relative "
                    f"adaptiveness w-index -- which normalises every synonymous "
                    f"family to its own peak and discards that peak, so the "
                    f"per-family scale %MinMax sums across families is not "
                    f"recoverable from it. Computing %MinMax from w would produce a "
                    f"different statistic under this one's name and citation",
                )
            try:
                freq = svc.tables.usage(reference_set)
                code = svc.tables.genetic_code(slot.table_id)
            except _MISSING_TABLE as exc:
                return self._unavailable(
                    c,
                    f"codon usage data for host {slot.host} could not be loaded "
                    f"(reference set {reference_set!r}, NCBI table {slot.table_id}): {exc}",
                )
            numeric = self._as_frequencies(freq)
            if numeric is None:
                return self._unavailable(
                    c,
                    f"the {reference_set!r} reference set did not supply a mapping of "
                    f"codon to usage frequency, so there are no X_ij to take %MinMax "
                    f"over",
                )

            codons = self._codons(c, editable, strand_for(ctx, slot))
            if codons is None:
                return self._unavailable(
                    c,
                    "the designable CDS is not a whole number of codons, so it cannot "
                    "be read in frame; an out-of-frame %MinMax profile describes a "
                    "protein that is not the one being made",
                )
            stats = family_statistics(code, numeric)
            informative = sum(1 for codon in codons if codon in stats)
            if not informative:
                return self._unavailable(
                    c,
                    "the CDS has no informative codons -- every position is a stop or "
                    "a single-codon family under this table, which carry no "
                    "codon-choice information and are excluded from %MinMax by "
                    "definition",
                )
            profile = min_max_profile(codons, stats, self.window)
            mean = sum(profile) / len(profile)
            per_slot.append((slot, reference_set, mean, profile, informative))

        # The slot furthest OUTSIDE the band binds, and among in-band slots the one
        # furthest from the band's middle. Averaging across slots would let a
        # comfortable slot hide one at +80, and +80 is the finding that matters.
        # The -1.0 offset keeps every in-band slot strictly below every breaching
        # one, so a breach always outranks an in-band slot however central it is.
        mid = (self.min_max_min + self.min_max_max) / 2

        def rank(row: tuple[ContextSlot, str, float, list[float], int]) -> float:
            dev = self._deviation(row[2])
            return dev if dev > 0.0 else -abs(row[2] - mid) - 1.0

        bound, reference_set, mean, profile, informative = max(per_slot, key=rank)

        side = self._side(mean)
        breaches: list[Breach] = []
        if side is not None:
            edge = self.min_max_min if side == "lower" else self.min_max_max
            why = (
                "a CDS pushed below the host's average synonymous usage across the "
                "ORF, which no evidence in the catalog recommends"
                if side == "lower"
                else "a CDS pushed toward the host's most common codons, which "
                "Ranaghan 2021 measured as a coin flip on yield and Welch 2009 "
                "measured expressing at ~15% of their best variant"
            )
            breaches.append(
                Breach(
                    spec_id=self.id,
                    # WHOLE_SCOPE: no sub-interval is "the" reason a mean over every
                    # window in the ORF left the band, so the finding is the CDS.
                    interval=Interval(editable[0].start, editable[-1].end),
                    magnitude=self._deviation(mean),
                    message=(
                        f"mean %MinMax {mean:+.1f} for the {bound.role} slot "
                        f"({bound.host}, {reference_set}) is outside the target band "
                        f"[{self.min_max_min:+.0f}, {self.min_max_max:+.0f}] "
                        f"({side} bound binding, "
                        f"{'under' if side == 'lower' else 'over'} {edge:+.0f}), "
                        f"indicating {why}. %MinMax is a soft band and a descriptive "
                        f"profile, never a target to maximize"
                    ),
                    # Codon choice is the only thing that moves %MinMax, and there is
                    # designable CDS here by construction (`editable` is non-empty).
                    fixable_by_codon_choice=True,
                    slot_role=bound.role,
                    detail={
                        "min_max_mean": mean,
                        "min_max_peak": max(profile),
                        "min_max_trough": min(profile),
                        "binding_side": side,
                        "band_lo": self.min_max_min,
                        "band_hi": self.min_max_max,
                        "host": str(bound.host),
                        # Part of the metric's definition: the same CDS scored
                        # against another reference set is a different number, so a
                        # %MinMax that travels without it is not reproducible.
                        "reference_set": reference_set,
                        "table_id": float(bound.table_id),
                        "window_codons": float(self.window),
                        "informative_codons": float(informative),
                    },
                )
            )

        return Evaluation(
            spec_id=self.id,
            # SOFT, so this verdict is a report signal only -- the solver never
            # routes a SOFT breach to repair (solver/catalog.py:189).
            passes=not breaches,
            raw_score=mean,
            breaches=tuple(breaches),
            windows=self._windows(c, editable, profile, strand_for(ctx, bound)),
            # The informative codons, not the profile's length: the profile has an
            # entry per codon including the ones that carry no choice, and the
            # honest denominator for "how much evidence is behind this number" is
            # the count that could have been chosen differently.
            n_evaluated=informative,
            binding_side=side,
        )

    # -- internals ----------------------------------------------------------

    def _side(self, mean: float) -> str | None:
        if mean < self.min_max_min:
            return "lower"
        if mean > self.min_max_max:
            return "upper"
        return None

    def _deviation(self, mean: float) -> float:
        """Distance outside the band; 0.0 when in band. Rule-native magnitude."""
        if mean < self.min_max_min:
            return self.min_max_min - mean
        if mean > self.min_max_max:
            return mean - self.min_max_max
        return 0.0

    def _codons(
        self, c: Construct, editable: Sequence[Interval], strand: Strand
    ) -> tuple[str, ...] | None:
        """The designable CDS as codons, in reading order for `strand`.

        `Construct.slice` is wrap- and strand-aware, so an insert cloned across the
        origin is read as one span rather than two truncated ones. On the reverse
        strand the segments are visited in reverse genomic order, because the last
        segment on the plus strand is the FIRST one a ribosome reading the minus
        strand meets. Returns None when the scope is not a whole number of codons.
        """
        spans = list(editable) if strand == 1 else list(reversed(editable))
        seq = "".join(c.slice(Interval(iv.start, iv.end, strand)) for iv in spans)
        if len(seq) % 3:
            return None
        return tuple(seq[i : i + 3] for i in range(0, len(seq), 3))

    def _windows(
        self,
        c: Construct,
        editable: Sequence[Interval],
        profile: Sequence[float],
        strand: Strand,
    ) -> tuple[tuple[Interval, float], ...]:
        """The profile mapped back to parent coordinates, one entry per codon.

        This is the report's %MinMax landscape and nothing else -- a window is
        never a breach (PLAN.md:661 asks for the profile to be printed). Codons are
        laid back down over the editable spans in the same reading order `_codons`
        read them, so a reverse-strand CDS gets its profile on the right bases.

        A codon whose three bases are not contiguous in construct coordinates -- it
        wraps the origin, or bridges two editable spans -- is DROPPED rather than
        given a plausible-looking interval, the same discipline `fragments` applies
        to a window with no parent coordinate. The profile itself keeps every
        position; only its projection onto the report loses the seam.
        """
        spans = list(editable) if strand == 1 else list(reversed(editable))
        offsets: list[int] = []
        for span in spans:
            positions = [p % c.length for p in range(span.start, span.end)]
            offsets.extend(positions if strand == 1 else list(reversed(positions)))
        out: list[tuple[Interval, float]] = []
        for index, value in enumerate(profile):
            base = index * 3
            if base + 2 >= len(offsets):
                break
            trio = sorted(offsets[base : base + 3])
            if trio[2] - trio[0] != 2:
                continue
            out.append((Interval(trio[0], trio[0] + 3), value))
        return tuple(out)

    def _as_frequencies(self, table: object) -> Mapping[str, float] | None:
        """`table` as a codon -> frequency mapping, or None if it is not one.

        `TableProvider.usage` is declared to return `Mapping[str, float]`
        (core/services.py:137), which is the frequency channel this rule needs.
        `FileTableProvider.usage` currently returns a `CodonUsage` dataclass
        instead -- a divergence from the protocol it implements, and one that only
        bites a caller who takes the protocol at its word, as this rule does. It is
        checked rather than assumed so the mismatch surfaces as a stated
        unavailability with a reason instead of an AttributeError from inside a
        rule, and so the day a frequency table lands the failure is loud.
        """
        if not isinstance(table, Mapping) or not table:
            return None
        out: dict[str, float] = {}
        for codon, value in table.items():
            # `bool` is an `int` subclass, so it would pass the numeric check and
            # silently become 0.0/1.0 usage.
            if not isinstance(codon, str) or isinstance(value, bool):
                return None
            if not isinstance(value, int | float):
                return None
            out[codon] = float(value)
        return out

    def _unavailable(self, c: Construct, reason: str) -> Evaluation:
        """NaN plus a breach carrying the reason -- B1's pattern, C1's argument.

        `Evaluation` has no "could not be computed" field (`ObjectiveScore.
        unavailable` is M3's type and is built from this downstream), so the reason
        travels in the one string channel a rule has. NaN rather than a number
        because every value in [-100, +100] is a real, meaningful %MinMax: 0.0 would
        read as a CDS sitting exactly at host-average usage, which is this band's
        ideal, and reporting the ideal for a quantity nobody measured is an
        affirmative false claim.

        `passes=True`: nothing about this construct failed. The objective was not
        computed, which is a different statement and must not read as a breach of a
        band that was never evaluated.
        """
        where = sorted(c.editable)[0] if c.editable else Interval(0, min(1, c.length) or 1)
        return Evaluation(
            spec_id=self.id,
            passes=True,
            raw_score=float("nan"),
            breaches=(
                Breach(
                    spec_id=self.id,
                    interval=where,
                    magnitude=0.0,
                    message=f"%MinMax objective unavailable: {reason}",
                    # Nothing about the coding sequence causes this and no codon
                    # choice fixes it -- the missing input is a reference set or a
                    # translating slot, not the DNA.
                    fixable_by_codon_choice=False,
                    detail={"unavailable_reason": reason},
                ),
            ),
            n_evaluated=0,
        )
