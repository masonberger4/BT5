"""C1 -- CAI as a soft band around the host's own usage, never a maximization target.

The brief's row is `2.C1`, under a section header that is itself the specification:
`### 2.C Codon composition (all S, soft bands, never maximized)` (brief.md:73). The
row asks for Sharp & Li's index -- `w_i = count_i / max_synonymous_count` over a
HIGHLY-EXPRESSED reference set, pseudocount 0.5 before the family max, and
`CAI = exp(mean ln w)` excluding stops and single-codon families -- scored against
a **band**, with the instruction "Never 1.0" (brief.md:77).

**The band is per host, and (0.70-0.90) is E. coli's.** brief.md:77 offers those
numbers as an example -- "e.g. 0.70-0.90, or +/-0.1 of host median" -- and they were
calibrated on a strong-bias organism. Measured against a composition-neutral random
synonymous encoding, E. coli's chance CAI is 0.238, so the 0.70 floor sits 0.46 above
chance; on the three mammalian w-tables chance is 0.656/0.633/0.660 and the same floor
sits 0.04-0.07 above chance, where it no longer discriminates. `CAI_BAND` therefore
carries one band per host: E. coli's pair unchanged, and for the weak-bias hosts an
inert 0.0 floor plus a ceiling scaled to that host's own chance-to-1.0 headroom. The
full argument -- including the rescaled floor that was rejected because it would flag
NATIVE human sequence as a finding -- is on `CAI_BAND` below, and the `band` ClassVar
is only the loosest envelope over it.

**Why a band, and why the direction is not HIGHER_IS_BETTER.** `Direction`'s own
docstring names this rule as the reason BAND exists: "a monotone weighted sum over a
collapsed |deviation| drives CAI toward 1.0, which is the precise failure the evidence
refutes (max-CAI collapses to one codon per amino acid and produces perfect nucleotide
repeats)". That is not a hypothetical. In Kudla 2009's 154-variant GFP series CAI
gave r = 0.14, not significant, against the 5' folding window's r = 0.66; Welch 2009
states outright that CAI "has no value in predicting gene expression" and their
deliberately high-CAI control expressed at ~15% of the best variant. The ceiling --
0.90 on E. coli, and its own on every other host -- is the operative half of this
rule, and the half that transfers: it is what stops the optimizer buying an
unmeasurable codon gain with a measurable repeat problem in E5/E6/F1.

**The evidence badge is CONTESTED, and the brief agrees in an unusual way.** The row's
own grade is `A (that it's weak)` -- grade-A evidence FOR THE CLAIM THAT CAI IS A WEAK
PREDICTOR, not grade-A evidence that CAI predicts anything. Meanwhile Boel 2016 measured
codon content as 3-5x MORE influential than structure on a soluble-protein readout under
T7, on 6,348 genes. Those two results are not reconciled and BT5 cannot adjudicate them,
which is what `Evidence.CONTESTED` is for. `EVIDENCE_BACKED` would badge a disputed
objective as settled; `FOLKLORE` would ship it disabled and silently drop an objective
all three presets weight.

**Codon data arrives through `Services`, never by importing M5.** `svc.tables.weights()`
and `svc.tables.genetic_code()` are the whole interface. `exp(mean ln w)` is recomputed
here rather than delegated to `CodonUsage.cai` precisely so the rules lane keeps no
import edge into the codon lane -- `Services` is what decouples M4 from M5, and a
single convenience import would erase that.

**Single-codon families are found from the TABLE, not from a hard-coded "M and W".**
Under NCBI table 4 TGA is Trp, so Trp has two codons and DOES carry information; under
tables 27 and 28 TGA is both Trp and a stop. `len(code.synonymous_codons(aa)) > 1` is
the only spelling that stays correct across tables, and `synonymous_codons` already
excludes stops from every family (CLAUDE.md 3.1, 3.2).

**Honest unavailability, because the data really is missing.** `data/codon_usage/`
now ships four reference sets -- Sharp & Li's E. coli w-index plus highly-expressed
human, mouse and CHO sets (S6, #90) -- against nine `HostId` values. So C1 computes a
number for E. coli K-12/BL21, HUMAN, HEK293, MOUSE and CHO, and reports the objective
UNAVAILABLE for S_CEREVISIAE, P_PASTORIS and SF9, which S6 deferred deliberately. All
three shipped presets now yield a CAI, `lentiviral_hek293` and `aav_hek293` included.
The three that remain unavailable are still the correct outcome and not a gap to paper
over: inventing a w-table for them -- or, worse, scoring an insect or yeast CDS against
whichever table happened to be on disk -- would produce a number that looks entirely
reasonable and measures nothing. A human highly-expressed
reference set is an evidence-bearing decision with its own provenance burden, and
`data/codon_usage/**` is a protected path; it belongs in its own labelled change.

**Why the gate is `role`, not `modality` or `host`.** CAI is a statement about
translation. The E. coli slot of a lentiviral job is a PROPAGATION slot -- the plasmid
is maintained there and the transgene is never expressed -- so a rule keyed on "is any
host E. coli" would find the one host BT5 has a table for, compute a confident CAI, and
report it as an objective for a design whose protein is made in HEK293. Excluding
propagation slots is what makes the mammalian presets report unavailable instead of
reporting the wrong number.
"""

from __future__ import annotations

import math
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

#: brief.md:77 -- "Target a **band** (e.g. 0.70-0.90, or +/-0.1 of host median).
#: Never 1.0." The ceiling is the operative half; see the module docstring.
#:
#: These are **E. coli's** bounds, not universal ones, and they are not the
#: constructor's defaults either -- `cai_min`/`cai_max` default to None so that a
#: caller passing 0.90 explicitly is distinguishable from one who passed nothing
#: (see `__init__`). They are the base `CAI_BAND` is derived from and the pair
#: E. coli keeps exactly.
BAND_LO = 0.70
BAND_HI = 0.90

#: HostId -> the band C1 actually scores against, because (0.70, 0.90) does not
#: mean the same thing on every host.
#:
#: brief.md:77 offers its numbers as an EXAMPLE -- "Target a band (e.g. 0.70-0.90,
#: or +/-0.1 of host median)" -- and they were calibrated on E. coli, a strong-bias
#: organism. Measured on the shipped tables, with a composition-neutral random
#: synonymous encoding as the chance baseline:
#:
#:     E. coli  chance CAI 0.238  -> the 0.70 floor sits 0.46 ABOVE chance
#:     human    chance CAI 0.656  -> 0.70 sits 0.04 above chance
#:     mouse    chance CAI 0.633  -> 0.07 above
#:     CHO      chance CAI 0.660  -> 0.04 above
#:
#: Mammalian codon bias is weak (brief.md:206: "isochore GC, not selection"), so
#: the w-tables are nearly flat and the same floor loses almost all of its
#: discriminating power. Two halves, and they do not transfer the same way:
#:
#: **The floor does not transfer at all, and is inert for the weak-bias hosts.**
#: Rescaling it to keep E. coli's headroom would put human's floor at 0.864 --
#: ABOVE where a native human CDS sits -- so C1 would flag native sequence as
#: "rare codons across the ORF" and hand the optimizer pressure to raise its CAI.
#: That is precisely what the evidence forbids: brief.md:206 marks the CAI weight
#: "very low" for CHO/HEK and the default mode "Native or harmonize", and
#: brief.md:13's 18-glycoprotein / 90-screen Expi293F benchmark concluded "codon
#: optimization to make human proteins, in a human cell line, did not generate
#: increased yields", with native constructs most consistent. There is no evidence
#: that a low-CAI mammalian CDS is worse, so nothing here claims one. The floor is
#: 0.0 and can never bind -- declared, like c3_min_max's, rather than left to look
#: symmetric.
#:
#: **The ceiling does transfer, and is the operative half.** Max-CAI collapse is a
#: MECHANICAL failure -- it drives each amino acid onto one codon and manufactures
#: perfect direct repeats -- and that is true of any organism. It is scaled to each
#: host's own chance-to-1.0 headroom so it means the same thing: E. coli's 0.90 is
#: 0.8687 of that headroom, and every ceiling below is that same fraction of its
#: own host's. E. coli's pair is therefore unchanged at exactly (0.70, 0.90).
#:
#: The ceilings are transcribed rather than computed at import: `chance_cai` needs a
#: w-table and a genetic code, and those arrive through `Services`, never by importing
#: M5 at module scope. `TestBandCalibration` re-derives each of them -- and the
#: rejected rescaled floor -- from the shipped tables through that same helper, so a
#: transcription that drifts from the data fails a test rather than shipping.
#:
#: Keyed on `HostId` alone, and the constants assume each host's default genetic code
#: (11 for E. coli, 1 for the mammals): the chance baseline is a function of the
#: (w-table, code) pair, so a slot pairing one of these hosts with a non-standard
#: `table_id` gets a ceiling calibrated for the standard one. The keys must stay in
#: step with `CAI_REFERENCE_SET` -- a host with a table but no band has no scoreable
#: band, and `_band_for` reports it unavailable rather than falling back to E. coli's.
CAI_BAND: Mapping[HostId, tuple[float, float]] = {
    HostId.E_COLI_K12: (BAND_LO, BAND_HI),
    HostId.E_COLI_BL21: (BAND_LO, BAND_HI),
    HostId.HUMAN: (0.0, 0.9548),
    HostId.HEK293: (0.0, 0.9548),
    HostId.MOUSE: (0.0, 0.9519),
    HostId.CHO: (0.0, 0.9553),
}

#: The fraction of a host's chance-to-1.0 headroom that its ceiling sits at,
#: taken from E. coli's published 0.90. The one number the table above is
#: derived from; the test re-derives it too.
CEILING_FRACTION_OF_HEADROOM = 0.8687

#: HostId -> the reference-set key `TableProvider.weights` is looked up under.
#:
#: The provider keys on the REFERENCE SET, not on the host, so this map is what
#: turns a slot's host into a lookup key; each value is the stem of a file in
#: `data/codon_usage/`. A host absent from this map has no CAI reference set in
#: this build and C1 reports its objective unavailable rather than borrowing one.
#:
#: THREE entries are approximations, and all three are stated rather than
#: hidden: the reference set used travels in every breach's
#: `detail["reference_set"]`, so a report can always name what it scored against.
#:
#: - BL21 shares K-12's entry. Sharp & Li's w-index was computed from
#:   highly-expressed E. coli K-12 genes; applying it to the B-strain is
#:   universal practice for the standard T7 expression host.
#: - CHO maps to a whole-organism *Cricetulus griseus* table, and this is the
#:   WEAKEST of the three: CHO is an aneuploid, heavily rearranged line, and
#:   78% of that table's contributing transcripts (71 of 91) are RefSeq
#:   PREDICTED `XM_` models rather than curated `NM_` ones, because C. griseus
#:   RefSeq is largely gene-prediction based. A predicted CDS can carry
#:   model-derived frame or boundary error, and codon counts are exactly what
#:   that corrupts. Disclosed by S6 in its own `_provenance`; repeated here
#:   because `TableProvider.weights()` returns only the bare `w` map, so a
#:   reader of a CHO breach never sees the provenance otherwise.
#: - HEK293 shares HUMAN's entry. HEK293 is a Homo sapiens cell line, so this is
#:   the same move one taxon down, and it is the mapping the data lane shipped
#:   the human set FOR: `docs/decisions/2026-09-02-s6-host-data-and-real-backbone.md`
#:   records "HEK293 -> human is the load-bearing one: both shipped mammalian
#:   presets (lentiviral_hek293, aav_hek293) key on HEK293". Codon usage is a
#:   property of the organism's translational machinery, not of the cell line.
#:
#: S_CEREVISIAE, P_PASTORIS and SF9 remain absent and therefore still report
#: unavailable -- S6 deferred them deliberately (no shipped preset consumes them,
#: and their RefSeq coverage is materially messier). That is the correct state,
#: not an oversight: a host absent from this map has no reference set in this
#: build, and C1 says so rather than borrowing one.
CAI_REFERENCE_SET: Mapping[HostId, str] = {
    HostId.E_COLI_K12: "sharp_li_1987_ecoli_w",
    HostId.E_COLI_BL21: "sharp_li_1987_ecoli_w",
    HostId.HUMAN: "human_highly_expressed_refseq_w",
    HostId.HEK293: "human_highly_expressed_refseq_w",
    HostId.MOUSE: "mouse_highly_expressed_refseq_w",
    HostId.CHO: "cho_highly_expressed_refseq_w",
}


def chance_cai(w: Mapping[str, float], code: GeneticCode) -> float:
    """Expected CAI of a random synonymous encoding of any protein. Exact, no RNG.

    The baseline the per-host ceilings are scaled against: what CAI a sequence scores
    for no reason other than the table's shape. Per informative family the expectation
    of ln w is the family mean of ln w, and CAI is the geometric mean over families,
    so the expected CAI is the geometric mean of the per-family geometric means
    outright -- no sampling.

    **For a protein of uniform amino-acid composition.** The mean over families here
    is unweighted, so this is the expectation for a protein carrying equal counts of
    all 18 informative amino acids, not for any particular one. A real proteome-like
    composition moves it by roughly 0.05 (human ~0.61 rather than ~0.656, E. coli
    ~0.20 rather than ~0.238). The ceilings are near-insensitive to that, because
    both terms of `chance + f * (1 - chance)` move together and `f` is itself derived
    from a chance value: recomputing the whole chain from composition-weighted
    baselines moves human's ceiling from 0.9548 to ~0.9513. The floor-above-chance
    figures quoted around this file are the uniform ones and would widen by ~0.05
    under a weighted composition, which strengthens rather than weakens the argument
    the mammalian floor is inert.

    Shares `_informative`'s family predicate deliberately: family size from the TABLE
    first, positive weight second. A baseline that admitted families the metric drops
    (or dropped families the metric admits) would calibrate a ceiling for a slightly
    different index than the one being scored, and nothing would say so.
    """
    logs: list[float] = []
    seen: set[str] = set()
    for codon in sorted(w):
        # The guard covers the whole table lookup, not only the family call: a
        # w-table carrying a non-canonical key (lowercase, `NNN`, a degenerate
        # codon) makes `is_stop`/`translate` raise KeyError, and a baseline that
        # died on one odd key would take the ceilings with it.
        try:
            if code.is_stop(codon):
                continue
            aa = code.translate(codon)
            if not aa or aa == "*" or aa in seen:
                continue
            seen.add(aa)
            family = code.synonymous_codons(aa)
        except (ValueError, KeyError):
            continue
        if len(family) < 2:
            continue
        weights = [w[c] for c in family if w.get(c, 0.0) > 0.0]
        if not weights:
            continue
        logs.append(sum(math.log(x) for x in weights) / len(weights))
    if not logs:
        return float("nan")
    return math.exp(sum(logs) / len(logs))


#: What `TableProvider` implementations raise when a host's table is absent.
#: `FileTableProvider` raises `FileNotFoundError` (an `OSError`); `NotImplementedError`
#: is what it raises for the tAI/stAI/CSC kinds that are a later lane.
_MISSING_TABLE = (OSError, LookupError, NotImplementedError, ValueError)


@register
class CodonAdaptationIndex:
    id: ClassVar[str] = "c1_cai"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "Codon adaptation index within a soft band, never maximized"
    enforcement: ClassVar[Enforcement] = Enforcement.SOFT
    #: See the module docstring: the brief grades C1 "A (that it's weak)" while
    #: Boel 2016 measured the opposite ordering on a different readout. A disputed
    #: objective badged EVIDENCE_BACKED would be asserting the dispute settled.
    evidence: ClassVar[Evidence] = Evidence.CONTESTED
    direction: ClassVar[Direction] = Direction.BAND
    unit: ClassVar[str] = "CAI (geometric mean of relative adaptiveness)"
    #: The LOOSEST envelope over `CAI_BAND`, not the gate -- e2_gc_band's
    #: convention, and `solver/catalog.py:272` states it: "Read off the INSTANCE,
    #: never `Spec.band`." Computed rather than transcribed so it cannot drift when
    #: a host is added. The band actually scored against is `_band_for(slot.host)`.
    band: ClassVar[tuple[float, float] | None] = (
        min(lo for lo, _ in CAI_BAND.values()),
        max(hi for _, hi in CAI_BAND.values()),
    )
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "Sharp & Li 1987, the index itself: w_i = count_i / max synonymous count "
            "over a HIGHLY-EXPRESSED reference set, and CAI as the geometric mean of w "
            "over sense codons. The reference set is part of the definition -- swapping "
            "it moves every CAI value, which a benchmark reads as an improvement",
            "https://pubmed.ncbi.nlm.nih.gov/3547335/",
            1987,
            sign="supports",
        ),
        Citation(
            "Kudla 2009, 154 synonymous GFP variants over a 250-fold expression range: "
            "CAI gave r = 0.14, NOT significant, while the -4..+37 folding window "
            "explained 44% of the variance (r = 0.66). This is why C1 carries 0.2 "
            "against B1's 1.0",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC3902468/",
            2009,
            sign="refutes",
        ),
        Citation(
            "Welch 2009: 'CAI has no value in predicting gene expression'; their "
            "deliberately high-CAI control expressed at ~15% of the best variant. The "
            "direct evidence for the 0.90 CEILING -- max-CAI is a measured failure mode, "
            "not a theoretical one",
            "https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0007002",
            2009,
            sign="refutes",
        ),
        Citation(
            "Boel 2016, 6,348 genes: codon content 3-5x MORE influential than structure "
            "on a soluble-protein readout under T7. The unresolved counter-result that "
            "keeps C1 in the objective at all rather than at weight zero, and the reason "
            "this rule is badged CONTESTED rather than EVIDENCE_BACKED",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC5054687/",
            2016,
            sign="supports",
        ),
        Citation(
            "The evidence that scopes this rule BY HOST, and the reason the "
            "mammalian floor is inert. brief.md:206 grades the CAI weight 'very "
            "low (isochore GC, not selection)' for CHO/HEK with default mode "
            "'Native or harmonize', and brief.md:215 marks the per-host evidence "
            "strength 'low for human, mouse, CHO, Sf9, Tni' against high for "
            "E. coli. An 18-glycoprotein / 90-screen Expi293F benchmark concluded "
            "'codon optimization to make human proteins, in a human cell line, did "
            "not generate increased yields', with native and harmonized constructs "
            "most consistent -- so nothing here claims a low-CAI mammalian CDS is "
            "worse than a native one",
            "https://proteininnovation.org/2026/03/codon-optimization-native-codon-mammalian-protein-expression/",
            2026,
            sign="refutes",
        ),
        Citation(
            "A 2026 Pichia study found CAI NEGATIVELY correlated with titer "
            "(coefficient -0.81 for trastuzumab). Cited as the sharpest statement "
            "of the direction this rule must not assume: higher CAI is not "
            "reliably better, which is why the band has a ceiling and why the "
            "weak-bias hosts have no operative floor",
            "https://europepmc.org/article/MED/41701818",
            2026,
            sign="refutes",
        ),
        Citation(
            "Ranaghan 2021 benchmarked nine commercial and academic optimizers and "
            "found 'a roughly equivalent chance that an algorithm-optimized CDS will "
            "increase or diminish recombinant yields', with three tools "
            "non-deterministic. Half of the ceiling C1 sits under -- and cited here "
            "rather than under Cambray below, which is 2018 and could not have "
            "reported a 2021 benchmark",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC7893858/",
            2021,
            sign="qualifies",
        ),
        Citation(
            "Cambray 2018: all computable design features together explain 5-31% "
            "(mean ~14%) of protein-level variance. The other half of that ceiling, "
            "and the reason C1's output is a band position rather than a predicted "
            "expression level",
            "https://www.nature.com/articles/nbt.4238",
            2018,
            sign="qualifies",
        ),
    )
    last_verified: ClassVar[str] = "2026-09-02"
    weight_provenance: ClassVar[str] = (
        "0.2, matching what all three shipped presets already assign 2.C1, and low on "
        "purpose. The brief grades C1's evidence 'A (that it's weak)' -- grade A for the "
        "claim that CAI is a WEAK predictor, not for any predictive power. The number "
        "behind that grade is Kudla 2009's r = 0.14 (not significant) for CAI against "
        "B1's r = 0.66 on the same 154 variants, so C1 is weighted at a fifth of B1 to "
        "reproduce that ordering rather than a guess. It is not zero because Boel 2016 "
        "measured codon content as 3-5x more influential than structure on a different "
        "readout, and BT5 cannot adjudicate that disagreement; a weight of zero would "
        "settle it by omission. steering_weight is 0.0 and must stay there: a Lagrangian "
        "term on codon weights would steer the Tier-A DP toward MAXIMUM CAI, which is "
        "the one thing this rule's own band exists to forbid. The `band` ClassVar is "
        "not a weight and not the gate: it is the LOOSEST envelope over `CAI_BAND`, "
        "whose bounds are per host because (0.70, 0.90) was calibrated on E. coli and "
        "sits 0.46 above that organism's chance CAI while the same floor sits 0.04-0.07 "
        "above chance on the mammalian tables. The gate is the slot host's own pair."
    )
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.2
    #: Zero, and structurally so -- not merely "not needed yet". See lattice_terms.
    steering_weight: ClassVar[float] = 0.0
    #: CAI is one number over the whole coding scope; there is no sub-interval to
    #: widen, because no single codon is the reason the aggregate left the band.
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.WHOLE_SCOPE
    #: SINGLE_PASS, and the downgrade is safe here in a way it is not for splice
    #: donors (CLAUDE.md 3.6): CAI is a monotone aggregate of independent codon
    #: choices, so moving codons toward the band cannot manufacture a new
    #: out-of-band region the way removing one donor creates another.
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "cheap"
    #: Raising CAI collapses each family toward its single top codon, which is how
    #: a max-CAI sequence manufactures perfect direct repeats and low-complexity
    #: tracts. `Direction`'s own docstring names this consequence.
    conflicts_with: ClassVar[tuple[str, ...]] = (
        "e5_synthesis_repeats",
        "e6_repeat_density",
        "f1_direct_repeats",
    )
    brief_ref: ClassVar[str] = "2.C1"
    #: CAI is arithmetic over a w-table, not a folding energy; there is no engine
    #: whose parameter set the band was calibrated against.
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        # No numeric defaults, e2_gc_band's reasoning applied per host instead of
        # per vendor: the gate is the slot host's own band, so a single advertised
        # number would lie -- a form showing 0.70/0.90 for a HEK293 job would be
        # describing E. coli's band. A value here overrides that side of it.
        "properties": {
            "cai_min": {
                "type": "number",
                "minimum": 0.0,
                "maximum": 1.0,
                "description": (
                    "Overrides the floor of the host's band (CAI_BAND); omitted, the "
                    "host's own floor applies -- 0.70 for E. coli, inert 0.0 for the "
                    "weak-bias mammalian hosts. Below the floor the CDS uses rare codons."
                ),
            },
            "cai_max": {
                "type": "number",
                "minimum": 0.0,
                "maximum": 1.0,
                "description": (
                    "Overrides the ceiling of the host's band (CAI_BAND); omitted, the "
                    "host's own ceiling applies. The operative half: it is what stops "
                    "the optimizer chasing CAI toward 1.0, which collapses each amino "
                    "acid to one codon and manufactures perfect repeats."
                ),
            },
        },
    }

    def __init__(self, cai_min: float | None = None, cai_max: float | None = None) -> None:
        """`None` means "unset", and it has to be a real sentinel.

        Not `cai_min=BAND_LO` compared by value: BAND_LO/BAND_HI are E. coli's
        bounds, so `cai_max=0.90` on a HEK293 job is a caller asking for the
        TIGHTER published ceiling, and reading that as "unset" would hand them the
        host's looser 0.9548 instead -- silently permitting the higher CAI this
        rule exists to refuse, through the rule's own parameter. `e2_gc_band.py`
        uses the same sentinel per vendor for the same reason.
        """
        for name, value in (("cai_min", cai_min), ("cai_max", cai_max)):
            if value is not None and not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must lie in [0, 1], got {value}")
        if cai_max is not None and cai_max >= 1.0:
            # "Never 1.0" (brief.md:77). A ceiling AT 1.0 is not a band -- it
            # permits the single-codon-per-amino-acid collapse this rule exists
            # to refuse, while still looking like a two-sided constraint.
            raise ValueError(
                "cai_max must be below 1.0: a ceiling at 1.0 permits the max-CAI "
                "collapse to one codon per amino acid, which is the failure mode "
                "the band exists to forbid"
            )
        if cai_min is not None and cai_max is not None and cai_min >= cai_max:
            raise ValueError(f"CAI band must satisfy cai_min {cai_min} < cai_max {cai_max}")
        # An explicit bound is the user's own number and wins over the per-host
        # band; tracked per side so a cai_min override does not silence the
        # host's ceiling, the discipline e2_gc_band applies per vendor. The
        # attributes below are only meaningful on the side whose flag is set --
        # `_band_for` composes the effective pair; None here means "this side is
        # the host's", and an instance built with no arguments must not report
        # `cai_min == 0.70` for a HEK293 job -- that is the exact confusion the
        # sentinel was introduced to remove.
        self._lo_override = cai_min is not None
        self._hi_override = cai_max is not None
        self.cai_min: float | None = cai_min
        self.cai_max: float | None = cai_max

    def gate(self, slot: ContextSlot) -> bool:
        """Every slot whose protein is actually translated.

        Not `modality`, and emphatically not `host`: see the module docstring.
        A propagation slot maintains the plasmid and never makes the protein, so
        scoring its codon adaptation would score a translation event that does
        not happen -- and, since E. coli is the one host with a w-table in this
        build, it is exactly the slot that would light up wrongly.
        """
        return slot.role != "propagation"

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """Deliberately None, and this is the load-bearing decision in the file.

        The tempting term is `LatticeTerms.codon_weights = {codon: log w}`. It
        would work, and it would be wrong: the Tier-A DP maximizes what it is
        given, so a codon-weight term steers monotonically toward MAXIMUM CAI --
        the exact behaviour "never maximized" (brief.md:73) forbids and that
        Welch 2009 measured expressing at ~15% of the best variant. A band is not
        expressible as a lattice term either: the automaton decides from a bounded
        codon suffix, while CAI is a geometric mean over the whole ORF, so no
        bounded state can know whether the next codon takes the aggregate above
        0.90. `forbidden` is likewise empty -- CAI is not a motif constraint and
        no individual codon is illegal.

        Tier B is the right and only place: the aggregate is measurable there,
        and being out of band is a report finding a human weighs, not a
        constraint the search chases (SOFT breaches never reach `repair()`;
        solver/catalog.py:189).
        """
        return None

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        slots = [s for s in ctx.active_slots if self.gate(s)]
        if not slots:
            return self._unavailable(
                c,
                "no translating slot in this context (propagation-only), so there is "
                "no host whose codon usage the CDS should be adapted to",
            )

        editable = sorted(c.editable)
        if not editable:
            return self._unavailable(c, "no designable CDS to compute codon adaptation over")

        # Every band composed BEFORE any per-slot work. `_band_for` raises when an
        # explicit bound inverts against a host's band, and resolving lazily inside
        # the loop made which failure you got depend on slot order: the same rule and
        # construct reported an honest `unavailable` for a deferred host when that
        # slot came first, and raised when it came second. A contradictory parameter
        # is a contradictory parameter whichever slot is examined first.
        bands = [self._band_for(s.host) for s in slots]

        # Each slot carries its OWN band: the hosts do not share one.
        per_slot: list[tuple[ContextSlot, str, float, int, tuple[float, float]]] = []
        for slot, band in zip(slots, bands, strict=True):
            reference_set = CAI_REFERENCE_SET.get(slot.host)
            if reference_set is None:
                return self._unavailable(
                    c,
                    f"no CAI reference set for host {slot.host}. This build ships "
                    f"reference sets for E. coli, human, mouse and CHO; "
                    f"{slot.host} is not among them. Scoring a {slot.host} CDS "
                    f"against another organism's weights would produce a "
                    f"plausible-looking number measuring nothing",
                )
            if band is None:
                return self._unavailable(
                    c,
                    f"no calibrated CAI band for host {slot.host}. A band is a "
                    f"statement about that host's own codon bias -- (0.70, 0.90) is "
                    f"E. coli's and sits 0.46 above its chance CAI, where the same "
                    f"floor sits 0.04-0.07 above chance on a weak-bias host -- so scoring "
                    f"{slot.host} against another organism's band would report a "
                    f"threshold nobody calibrated for it",
                )
            try:
                w = svc.tables.weights(reference_set, "cai")
                code = svc.tables.genetic_code(slot.table_id)
            except _MISSING_TABLE as exc:
                return self._unavailable(
                    c,
                    f"codon data for host {slot.host} could not be loaded "
                    f"(reference set {reference_set!r}, NCBI table {slot.table_id}): {exc}",
                )
            if not w:
                return self._unavailable(
                    c,
                    f"the {reference_set!r} reference set supplied no relative "
                    f"adaptiveness weights, so CAI has no w to take a geometric mean of",
                )

            codons = self._codons(c, editable, strand_for(ctx, slot))
            if codons is None:
                return self._unavailable(
                    c,
                    "the designable CDS is not a whole number of codons, so it cannot "
                    "be read in frame; an out-of-frame CAI is a number about a protein "
                    "that is not the one being made",
                )
            logs = [math.log(weight) for weight in self._informative(codons, code, w)]
            if not logs:
                return self._unavailable(
                    c,
                    "the CDS has no informative codons -- every position is a stop or a "
                    "single-codon family (ATG/TGG under most tables), which carry no "
                    "codon-choice information and are excluded from CAI by definition",
                )
            per_slot.append(
                (
                    slot,
                    reference_set,
                    math.exp(sum(logs) / len(logs)),
                    len(logs),
                    band,
                )
            )

        # The slot furthest OUTSIDE ITS OWN band binds, and among in-band slots the
        # one CLOSEST to a binding edge. Averaging CAI across slots would let a
        # comfortable slot hide one at 0.97, and 0.97 is the finding that matters.
        # Comparing raw CAI across slots would be worse still now that the bands
        # differ per host: 0.85 is comfortably in band for human and a breach for
        # E. coli. The -1.0 offset keeps every in-band slot strictly below every
        # breaching one (deviations are >= 0 and the normalized distance below is in
        # [0, 1]), so a breach always outranks an in-band slot whichever band it
        # breached.
        #
        # Distance to the nearest BINDING EDGE, as a fraction of the band's width,
        # and that normalization is load-bearing rather than cosmetic. Two earlier
        # spellings both fail:
        #
        # - Distance to the band's MIDDLE. With an inert 0.0 floor the middle is
        #   0.477, below where any real mammalian CDS ever sits, so every mammalian
        #   slot ranks as maximally uninteresting and the report always goes to the
        #   E. coli slot -- discarding the CAI the job is actually about.
        # - RAW distance to an edge. E. coli's band is 0.20 wide, so an in-band
        #   E. coli slot is always within 0.10 of an edge, while a mammalian band is
        #   ~0.955 wide. A HEK293 CDS 0.028 below its max-CAI ceiling still loses to
        #   an E. coli slot merely near its floor -- the same pathology, narrower.
        #
        # An INERT edge is excluded outright: `lo <= 0.0` can never bind (CAI is a
        # geometric mean of positive weights, so it is strictly positive), and
        # "0.2 above a floor that cannot fire" is not a distance to anything.
        #
        # Deviations ARE compared across hosts here, and they are in raw CAI units,
        # unnormalized. Human's chance-to-1.0 headroom is 0.344 against E. coli's
        # 0.762, so a full max-CAI collapse yields |1 - hi| = 0.045 on HEK293 against
        # 0.100 on E. coli -- a factor of 2.2, which is exactly that headroom ratio.
        # A headroom-normalized magnitude would be more comparable within C1, but
        # `Breach.magnitude` is unnormalized repo-wide and one rule normalizing
        # unilaterally would misrank against the others in `score/conflicts.py`.
        # Recorded so the ordering is not read as calibrated.

        def rank(row: tuple[ContextSlot, str, float, int, tuple[float, float]]) -> float:
            lo, hi = row[4]
            dev = self._deviation(row[2], row[4])
            if dev > 0.0:
                return dev
            edges = [hi - row[2]] if lo <= 0.0 else [row[2] - lo, hi - row[2]]
            return -min(edges) / (hi - lo) - 1.0

        bound, reference_set, cai, n, band = max(per_slot, key=rank)
        cai_min, cai_max = band

        side = self._side(cai, band)
        breaches: list[Breach] = []
        if side is not None:
            edge = cai_min if side == "lower" else cai_max
            why = (
                "rare codons across the ORF"
                if side == "lower"
                else "the max-CAI collapse toward one codon per amino acid, which "
                "manufactures perfect direct repeats"
            )
            breaches.append(
                Breach(
                    spec_id=self.id,
                    # WHOLE_SCOPE: no sub-interval is "the" reason a geometric mean
                    # over the whole ORF left the band, so the finding is the CDS.
                    interval=Interval(editable[0].start, editable[-1].end),
                    magnitude=self._deviation(cai, band),
                    message=(
                        f"CAI {cai:.3f} for the {bound.role} slot ({bound.host}, "
                        f"{reference_set}) is outside the target band "
                        f"[{cai_min:.3f}, {cai_max:.3f}] ({side} bound "
                        f"binding, {'under' if side == 'lower' else 'over'} "
                        f"{edge:.3f}), indicating {why}. CAI is a soft band and a "
                        f"descriptive statistic: r = 0.14 (not significant) against "
                        f"expression in Kudla 2009, so this is a finding to weigh, "
                        f"never a target to maximize"
                    ),
                    # Codon choice is the ONLY thing that moves CAI, and there is
                    # designable CDS here by construction (`editable` is non-empty).
                    fixable_by_codon_choice=True,
                    slot_role=bound.role,
                    detail={
                        "cai": cai,
                        "binding_side": side,
                        "band_lo": cai_min,
                        "band_hi": cai_max,
                        "host": str(bound.host),
                        # The reference set is part of CAI's definition: the same
                        # CDS scored against another set is a different number, so
                        # a CAI that travels without it is not reproducible.
                        "reference_set": reference_set,
                        "table_id": float(bound.table_id),
                        "informative_codons": float(n),
                    },
                )
            )

        return Evaluation(
            spec_id=self.id,
            # SOFT, so this verdict is a report signal only -- the solver never
            # routes a SOFT breach to repair (solver/catalog.py:189). Saying
            # `passes=False` for an out-of-band CAI is what puts it in front of a
            # human, which for an objective this weakly evidenced is the point.
            passes=not breaches,
            raw_score=cai,
            breaches=tuple(breaches),
            n_evaluated=n,
            binding_side=side,
        )

    # -- internals ----------------------------------------------------------

    def _band_for(self, host: HostId) -> tuple[float, float] | None:
        """The band this host is actually scored against, or None if it has none.

        `CAI_BAND` supplies each side, and an explicit constructor bound overrides
        it -- per side, so setting one does not discard the other.

        **No fallback to (BAND_LO, BAND_HI) for an unmapped host.** That would
        reinstate the exact bug this map exists to fix, silently and one host at a
        time: the next host to gain a reference set would be scored against
        E. coli's band with no error and no test failure. A host with a w-table but
        no calibrated band has no band, and the caller reports it unavailable.
        Both bounds supplied explicitly is the one exception -- then the band is
        the caller's own and needs no host calibration.
        """
        base = CAI_BAND.get(host)
        if base is None:
            if self.cai_min is None or self.cai_max is None:
                return None
            base = (BAND_LO, BAND_HI)  # both sides are about to be overridden
        lo = base[0] if self.cai_min is None else self.cai_min
        hi = base[1] if self.cai_max is None else self.cai_max
        if lo >= hi:
            # Reachable only by composition: one explicit bound against a host
            # bound, e.g. cai_min=0.93 on E. coli's 0.90 ceiling. `__init__` cannot
            # catch it because the host is not known until evaluation, and a band
            # that is silently inverted reports a max-CAI sequence as "rare codons".
            raise ValueError(
                f"cai_min {lo} is not below cai_max {hi} for host {host}: an "
                f"explicit bound composed against that host's band inverts it"
            )
        return (lo, hi)

    def _side(self, cai: float, band: tuple[float, float]) -> str | None:
        lo, hi = band
        if cai < lo:
            return "lower"
        if cai > hi:
            return "upper"
        return None

    def _deviation(self, cai: float, band: tuple[float, float]) -> float:
        """Distance outside the band; 0.0 when in band. Rule-native magnitude."""
        lo, hi = band
        if cai < lo:
            return lo - cai
        if cai > hi:
            return cai - hi
        return 0.0

    def _codons(
        self, c: Construct, editable: Sequence[Interval], strand: Strand
    ) -> tuple[str, ...] | None:
        """The designable CDS as codons, in reading order for `strand`.

        `Construct.slice` is wrap- and strand-aware, so an insert cloned across
        the origin is read as one span rather than as two truncated ones. On the
        reverse strand the segments are visited in reverse genomic order, because
        the last segment on the plus strand is the FIRST one a ribosome reading
        the minus strand meets -- concatenating them in plus-strand order would
        put the codons in the wrong sequence and give a CAI for a protein nobody
        makes. Returns None when the scope is not a whole number of codons.
        """
        spans = list(editable) if strand == 1 else list(reversed(editable))
        seq = "".join(c.slice(Interval(iv.start, iv.end, strand)) for iv in spans)
        if len(seq) % 3:
            return None
        return tuple(seq[i : i + 3] for i in range(0, len(seq), 3))

    def _informative(
        self, codons: Sequence[str], code: GeneticCode, w: Mapping[str, float]
    ) -> list[float]:
        """The w of every codon that carries codon-choice information.

        Excluded, per the brief and Sharp & Li's definition: stops, and any codon
        whose amino acid has exactly one non-stop codon under THIS table. The
        family size is read from the table rather than hard-coded to ATG/TGG
        because it is table-dependent -- NCBI table 4 makes TGA a second Trp
        codon, so Trp does carry information there, and tables 27/28 make TGA
        both Trp and a stop (CLAUDE.md 3.1, 3.2). A codon the reference set has
        no positive weight for is skipped rather than treated as zero, which
        would send ln w to -infinity and collapse the whole geometric mean.
        """
        out: list[float] = []
        for codon in codons:
            if code.is_stop(codon):
                continue
            aa = code.translate(codon)
            if not aa or aa == "*":
                continue
            try:
                if len(code.synonymous_codons(aa)) < 2:
                    continue
            except (ValueError, KeyError):
                continue
            weight = w.get(codon)
            if weight is None or weight <= 0.0:
                continue
            out.append(weight)
        return out

    def _unavailable(self, c: Construct, reason: str) -> Evaluation:
        """NaN plus a breach carrying the reason -- B1's pattern, same argument.

        `Evaluation` has no "could not be computed" field (`ObjectiveScore.
        unavailable` is M3's type and is built from this downstream), so the
        reason travels in the one string channel a rule has. NaN rather than a
        number because every value in [0, 1] is a real, meaningful CAI: reporting
        0.0 would read as a catastrophically rare-codon sequence, and reporting
        the band's midpoint would read as a design that is exactly on target.
        Both are affirmative false claims about a quantity nobody measured.

        `passes=True`: nothing about this construct failed. The objective was not
        computed, which is a different statement and must not read as a breach of
        a band that was never evaluated.
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
                    message=f"CAI objective unavailable: {reason}",
                    # Nothing about the coding sequence causes this and no codon
                    # choice fixes it -- the missing input is a reference set or
                    # a translating slot, not the DNA.
                    fixable_by_codon_choice=False,
                    detail={"unavailable_reason": reason},
                ),
            ),
            n_evaluated=0,
        )
