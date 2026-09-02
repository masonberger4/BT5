"""D8 -- CpG, as three separate metrics that must not become one slider.

`brief.md:128` is unusually prescriptive about the shape: *"three separate,
separately-toggleable metrics (do not collapse into one slider)"*. They are separate
because they disagree. TLR9 sensing counts CpGs; ZAP-mediated decay explicitly does
NOT -- *"the magnitude of ZAP-mediated inhibition was not correlated with the number of
CpGs introduced"* (`brief.md:130`) -- it is driven by SPACING; and methylation silencing
is about a compositional island, not a count at all. A single "CpG depletion" slider
would average three mechanisms that move in different directions and report a number
that means nothing for any of them.

So each metric is its own toggle, each contributes its own breaches with its own
`detail["metric"]`, and the rule reports the WORST ZAP window rather than a global CpG
count, exactly as `brief.md:130` instructs.

**Enforcement is SOFT.** `brief.md:128`'s header carries no H/S marker, unlike D3's
`(H/S, ...)` or D4's `(H for lentiviral sense strand, S elsewhere)`. Absent a hard grade
this stays a weighted objective and is never allowed to refuse a construct -- and CpG
content is one of the places over-enforcement has a measurable cost, since
`brief.md:133` warns that full depletion forces AGA/AGG for Arg and can drop GC below
vendor floors. That is a real conflict with `e2_gc_band`, which already names this rule
in its own `conflicts_with`.

**No strand handling, and that is not an oversight.** `CG` is its own reverse
complement, so a CpG dinucleotide is the same object on both strands. This rule has no
`strand_of_interest` to read and no reverse scan to do -- unlike D3 and D4, where
getting the strand wrong is the whole hazard.

**CpG arises across codon boundaries** (`brief.md:132`): within codons via CGN Arg,
GCG, CCG, TCG and ACG, *and* wherever one codon ends C and the next begins G. Evaluating
the dinucleotide on the assembled construct rather than per codon is what catches the
second kind, which is also why this is a `Construct` rule and not a string rule.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from typing import ClassVar

from bt5.core.context import ContextSlot, DesignContext, HostId
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

# -- (a) TLR9 -----------------------------------------------------------------

#: `brief.md:129`, with the species each is attributed to, verbatim.
TLR9_SPECIFIC: Mapping[str, str] = {
    "GTCGTT": "human",
    "TTCGTT": "human",
    "GACGTT": "mouse",
}

#: `brief.md:129`: stimulatory hexamers `RRCGYY`, *"specifically"* the three above.
#:
#: TWO OF THE THREE ARE NOT `RRCGYY`, and reading "specifically" as "these are examples
#: of it" produces a rule that silently never fires on either human hexamer:
#:
#:   GTCGTT -> G T C G T T; position 2 is T, and R is A or G.  NOT a match.
#:   TTCGTT -> T T C G T T; position 1 is T.                   NOT a match.
#:   GACGTT -> G A C G T T.                                    matches.
#:
#: That is not a defect in the brief -- the classic human TLR9 motifs really are not
#: the RRCGYY consensus -- so the named hexamers are scanned as their own alternatives
#: UNIONED with the consensus, never as a filter applied to it. Pinned by
#: `TestTlr9::test_the_named_human_hexamers_are_not_rrcgyy`.
TLR9_CONSENSUS = "[AG][AG]CG[CT][CT]"
TLR9_MOTIF = re.compile("|".join([*TLR9_SPECIFIC, TLR9_CONSENSUS]))
TLR9_CONSENSUS_ONLY = re.compile(TLR9_CONSENSUS)

#: Only the two species `brief.md:129` names are mapped. A host it does not name gets
#: the general RRCGYY finding and no species escalation, rather than a guess.
HOST_SPECIES: Mapping[HostId, str] = {
    HostId.HUMAN: "human",
    HostId.HEK293: "human",
    HostId.MOUSE: "mouse",
}

# -- (b) ZAP / KHNYN ----------------------------------------------------------

#: `brief.md:130`: "sliding 200-nt window; flag any window with >=14 CpGs at mean
#: inter-CpG spacing <=14 nt (peak sensitivity 6-14 nt). Spacing >=32 nt is not
#: restricted."
ZAP_WINDOW = 200
ZAP_MIN_CPG = 14
ZAP_MAX_MEAN_SPACING = 14.0
ZAP_UNRESTRICTED_SPACING = 32.0

# -- (c) Methylation islands --------------------------------------------------

#: `brief.md:131`: Gardiner-Garden (>=200 bp, GC >=50%, obs/exp >=0.6) or Takai-Jones
#: strict (>=500 bp, GC >=55%, >=0.65).
GARDINER_GARDEN = (200, 0.50, 0.60)
TAKAI_JONES = (500, 0.55, 0.65)

MAG_ISLAND = 1.0
MAG_ZAP = 1.0
MAG_TLR9_SPECIFIC = 0.6
MAG_TLR9_GENERAL = 0.2


def _cpg_positions(text: str) -> list[int]:
    """Every CpG start in `text`. `CG` is its own reverse complement, so one pass
    over the forward strand is the whole answer."""
    out: list[int] = []
    i = text.find("CG")
    while i != -1:
        out.append(i)
        i = text.find("CG", i + 1)
    return out


def obs_over_exp(chunk: str) -> float:
    """`brief.md:131`: obs/exp = (N_CpG x L) / (N_C x N_G).

    Zero when the window has no C or no G: the ratio is undefined there, and a
    window with no C cannot be a CpG island under any definition.
    """
    n_c, n_g = chunk.count("C"), chunk.count("G")
    if n_c == 0 or n_g == 0:
        return 0.0
    return (chunk.count("CG") * len(chunk)) / (n_c * n_g)


@register
class CpGDepletion:
    id: ClassVar[str] = "d8_cpg_depletion"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "CpG: TLR9 sensing, ZAP/KHNYN decay, and methylation islands"
    enforcement: ClassVar[Enforcement] = Enforcement.SOFT
    #: CONTESTED, and the ZAP citation is why: the same paper that establishes the
    #: mechanism reports that inhibition did NOT correlate with the number of CpGs
    #: introduced. Grading this EVIDENCE_BACKED would badge a contested quantity as
    #: settled in exactly the place the evidence badge is supposed to be load-bearing.
    evidence: ClassVar[Evidence] = Evidence.CONTESTED
    direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
    unit: ClassVar[str] = "weighted CpG findings"
    band: ClassVar[tuple[float, float] | None] = None
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "CpG-depleted AAVrh32.33 evaded immune detection -- the TLR9 arm",
            "https://www.jci.org/articles/view/68205",
            2013,
            sign="supports",
        ),
        Citation(
            "ZAP/KHNYN restriction peaks at 6-14 nt inter-CpG spacing and is not "
            "restricted at spacing >=32 nt, which is why this rule reports the worst "
            "200-nt window rather than a global count",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC9519448/",
            2022,
            sign="supports",
        ),
        Citation(
            "The same work reports that the magnitude of ZAP-mediated inhibition was NOT "
            "correlated with the number of CpGs introduced -- so a CpG COUNT is the wrong "
            "quantity for this arm, and depleting CpGs is not monotonically protective",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC9519448/",
            2022,
            sign="refutes",
        ),
    )
    last_verified: ClassVar[str] = "2026-09-02"
    weight_provenance: ClassVar[str] = (
        "Low-to-mid, and lower than d4_internal_polya's 0.7 on purpose. Two of the three "
        "metrics rest on mechanisms whose magnitude is contested in their own primary "
        "source -- ZAP inhibition did not correlate with CpG number -- and the third, "
        "methylation, is a compositional heuristic rather than a measured effect on any "
        "construct BT5 builds. The weight also has to price a KNOWN cost on the other "
        "side: brief.md:133 records that full depletion forces AGA/AGG for Arg and can "
        "drop GC below vendor floors, so a high weight here buys a contested benefit with "
        "a certain manufacturability loss. e2_gc_band already declares that conflict."
    )
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.35
    steering_weight: ClassVar[float] = 0.15
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.WINDOW_MINUS_1
    #: SOFT rules are not repaired, so the policy is inert. Declared rather than
    #: defaulted: removing a CpG cannot create a CpG elsewhere -- a synonymous codon
    #: swap either leaves the dinucleotide or does not -- so even under repair this
    #: would not need FIXED_POINT, unlike d3_splicing.
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "cheap"
    #: brief.md:133: full depletion forces AGA/AGG for Arg and can drop GC below vendor
    #: floors. e2_gc_band already names this rule from its side.
    conflicts_with: ClassVar[tuple[str, ...]] = ("e2_gc_band",)
    brief_ref: ClassVar[str] = "2.D8"
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "tlr9": {
                "type": "boolean",
                "default": True,
                "description": "Metric (a): CpG count and stimulatory RRCGYY hexamers.",
            },
            "zap": {
                "type": "boolean",
                "default": True,
                "description": (
                    "Metric (b): worst 200-nt window by CpG count and mean inter-CpG "
                    "spacing. Reports one window, never a global count."
                ),
            },
            "methylation": {
                "type": "boolean",
                "default": True,
                "description": (
                    "Metric (c): CpG islands by Gardiner-Garden or Takai-Jones strict."
                ),
            },
            "island_criterion": {
                "type": "string",
                "enum": ["gardiner_garden", "takai_jones"],
                "default": "gardiner_garden",
                "description": (
                    "Gardiner-Garden (>=200 bp, GC >=50%, obs/exp >=0.6) is the looser "
                    "and the default; Takai-Jones strict is >=500 bp, GC >=55%, >=0.65."
                ),
            },
        },
    }

    def __init__(
        self,
        tlr9: bool = True,
        zap: bool = True,
        methylation: bool = True,
        island_criterion: str = "gardiner_garden",
    ) -> None:
        if island_criterion not in ("gardiner_garden", "takai_jones"):
            raise ValueError(
                f"island_criterion must be gardiner_garden or takai_jones, got {island_criterion!r}"
            )
        self.tlr9 = tlr9
        self.zap = zap
        self.methylation = methylation
        self.island_criterion = island_criterion
        #: Read by solver/catalog.py:236 if this rule is ever escalated to repair.
        self.window: int = ZAP_WINDOW

    def gate(self, slot: ContextSlot) -> bool:
        """Every context. TLR9 and ZAP are host-immune mechanisms and methylation is a
        mammalian one, but a plasmid propagated in E. coli still carries the CpGs it
        will have in the target cell, so the finding is worth reporting wherever the
        construct is being built."""
        return True

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """None. Forbidding `CG` outright is not a constraint, it is a different
        organism: CpG appears in CGN Arg, GCG, CCG, TCG and ACG, so the automaton would
        make most of the Arg/Ala/Pro/Ser/Thr codon space unreachable to express a soft
        preference. Steering is the mechanism for that (CLAUDE.md 3.5)."""
        return None

    # -- metrics --------------------------------------------------------------

    def _tlr9(
        self, c: Construct, text: str, offset: int, n: int, species: str | None
    ) -> list[Breach]:
        breaches: list[Breach] = []
        # `brief.md:129` makes metric (a) "total CpG count + stimulatory hexamers", so
        # the count is half the metric and is reported even when no hexamer is present.
        # Magnitude 0.0 and not fixable: it is a measurement, not a finding, so it
        # stays out of `passes`, out of the weighted sum's meaning, and off the
        # solver's target list -- the same channel d3_splicing uses to say a scan did
        # not run.
        total = sum(1 for p in _cpg_positions(text) if offset <= p < offset + n)
        if total:
            where = sorted(c.editable)[0] if c.editable else Interval(0, min(1, n) or 1)
            breaches.append(
                Breach(
                    spec_id=self.id,
                    interval=where,
                    magnitude=0.0,
                    message=(
                        f"total CpG count {total} over {n} nt "
                        f"({total / n:.3f} per nt); reported because TLR9 sensing scales "
                        "with total CpG, unlike the ZAP arm below"
                    ),
                    fixable_by_codon_choice=False,
                    detail={
                        "metric": "tlr9",
                        "reading": "total_cpg",
                        "cpg_total": float(total),
                        "cpg_per_nt": total / n,
                    },
                )
            )
        for match in TLR9_MOTIF.finditer(text):
            if not offset <= match.start() < offset + n:
                continue
            hexamer = match.group()
            named = TLR9_SPECIFIC.get(hexamer)
            lo = (match.start() - offset) % n
            iv = Interval(lo, lo + len(hexamer))
            breaches.append(
                Breach(
                    spec_id=self.id,
                    interval=iv,
                    # Escalated only when the hexamer is one brief.md:129 names AND the
                    # slot's host is the species it names it for. A mouse-specific
                    # hexamer in a human cassette is a weaker finding, and the brief
                    # maps no other host, so no other host is guessed at.
                    magnitude=(
                        MAG_TLR9_SPECIFIC
                        if named is not None and named == species
                        else MAG_TLR9_GENERAL
                    ),
                    message=(
                        f"TLR9 stimulatory hexamer {hexamer!r} at {lo}"
                        + (f", documented for {named}" if named else ", matching RRCGYY")
                    ),
                    fixable_by_codon_choice=c.overlaps_editable(iv),
                    detail={
                        "metric": "tlr9",
                        "hexamer": hexamer,
                        "species": named or "unattributed",
                    },
                )
            )
        return breaches

    def _zap(self, c: Construct, text: str, offset: int, n: int) -> list[Breach]:
        """The single worst 200-nt window, never a global count (`brief.md:130`).

        Two pointers over the CpG positions rather than a slice per window. This rule
        is SOFT, but `RuleSet.findings` evaluates EVERY spec and filters afterwards --
        its "scoped to the HARD_REPAIR rules" docstring describes `breach_finder`'s
        intent, not what the loop does -- so a soft rule is still paid for once per
        candidate, up to 256 per repair iteration. Re-slicing 200 characters 5,000
        times per evaluation is minutes of wall clock at that multiplier.
        """
        worst: tuple[float, int, int, float] | None = None
        span = min(ZAP_WINDOW, n)
        cpg = _cpg_positions(text)
        lo_i = hi_i = 0
        for start in range(offset, offset + n):
            if start + span > len(text):
                break
            while lo_i < len(cpg) and cpg[lo_i] < start:
                lo_i += 1
            while hi_i < len(cpg) and cpg[hi_i] < start + span - 1:
                hi_i += 1
            count = hi_i - lo_i
            if count < ZAP_MIN_CPG:
                continue
            spacing = (cpg[hi_i - 1] - cpg[lo_i]) / (count - 1)
            if spacing > ZAP_MAX_MEAN_SPACING:
                continue
            score = count / max(spacing, 1e-9)
            if worst is None or score > worst[0]:
                worst = (score, start, count, spacing)
        if worst is None:
            return []
        _, start, count, spacing = worst
        lo = (start - offset) % n
        iv = Interval(lo, lo + span)
        return [
            Breach(
                spec_id=self.id,
                interval=iv,
                magnitude=MAG_ZAP,
                message=(
                    f"ZAP/KHNYN: worst {span}-nt window at {lo} carries {count} CpGs at "
                    f"{spacing:.1f} nt mean spacing (restriction peaks at 6-14 nt; "
                    f">={ZAP_UNRESTRICTED_SPACING:.0f} nt is not restricted). Note the "
                    "source reports inhibition did NOT scale with CpG number, so this "
                    "is a spacing finding and not a count"
                ),
                fixable_by_codon_choice=c.overlaps_editable(iv),
                detail={
                    "metric": "zap",
                    "cpg_count": float(count),
                    "mean_spacing": spacing,
                },
            )
        ]

    def _islands(self, c: Construct, text: str, offset: int, n: int) -> list[Breach]:
        min_len, min_gc, min_ratio = (
            GARDINER_GARDEN if self.island_criterion == "gardiner_garden" else TAKAI_JONES
        )
        if n < min_len:
            return []
        breaches: list[Breach] = []
        # Step by a fifth of the window, the same stride e2_gc_band uses: a full
        # per-base sweep of a 500 bp criterion over a plasmid buys resolution the
        # threshold cannot use.
        step = max(1, min_len // 5)
        best: tuple[float, int] | None = None
        for start in range(offset, offset + n, step):
            chunk = text[start : start + min_len]
            if len(chunk) < min_len:
                break
            gc = (chunk.count("G") + chunk.count("C")) / len(chunk)
            if gc < min_gc:
                continue
            ratio = obs_over_exp(chunk)
            if ratio < min_ratio:
                continue
            if best is None or ratio > best[0]:
                best = (ratio, start)
        if best is not None:
            ratio, start = best
            lo = (start - offset) % n
            iv = Interval(lo, lo + min_len)
            chunk = text[start : start + min_len]
            gc = (chunk.count("G") + chunk.count("C")) / len(chunk)
            breaches.append(
                Breach(
                    spec_id=self.id,
                    interval=iv,
                    magnitude=MAG_ISLAND,
                    message=(
                        f"CpG island ({self.island_criterion}) at {lo}: {min_len} bp, "
                        f"GC {gc:.0%}, obs/exp {ratio:.2f} -- a methylation-silencing "
                        "substrate in a mammalian host"
                    ),
                    fixable_by_codon_choice=c.overlaps_editable(iv),
                    detail={
                        "metric": "methylation",
                        "criterion": self.island_criterion,
                        "obs_exp": ratio,
                        "gc": gc,
                    },
                )
            )
        return breaches

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        text, offset = c.tripled()
        n = c.length
        species = next(
            (HOST_SPECIES[s.host] for s in ctx.active_slots if s.host in HOST_SPECIES), None
        )
        breaches: list[Breach] = []
        if self.tlr9:
            breaches += self._tlr9(c, text, offset, n, species)
        if self.zap:
            breaches += self._zap(c, text, offset, n)
        if self.methylation:
            breaches += self._islands(c, text, offset, n)

        return Evaluation(
            spec_id=self.id,
            # A ZAP window and a CpG island each cross a threshold brief.md STATES;
            # a TLR9 hexamer note does not, so it is a warn-band finding that must not
            # read as a failure (solver/catalog.py:158-170). Enforcement is SOFT either
            # way, so `passes` never refuses a construct here -- it only decides what
            # counts as this rule's own verdict.
            passes=not any(b.magnitude >= MAG_ISLAND for b in breaches),
            raw_score=float(sum(b.magnitude for b in breaches)),
            breaches=tuple(breaches),
            n_evaluated=n,
        )
