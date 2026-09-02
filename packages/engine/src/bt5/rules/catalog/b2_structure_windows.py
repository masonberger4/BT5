"""B2 -- 5' secondary structure in Cambray's two windows, scored with its sign.

brief.md:62, graded **A**, bacteria: "Penalize structure in **STR(-30:+30)**;
neutral-to-mildly-reward moderate structure in **STR(+31:+90)**." The windows are
Cambray 2018's, from 244,000 designed sequences in full factorial, where phenotype
was "dominated by secondary structures and their interactions" and that dominance
localised to exactly these two spans.

**Only the proximal window is scored, and that is the brief's instruction, not a
shortcut.** brief.md:333 is explicit: structure's sign flips between contexts and
these are "three different objectives on three different molecules -- model them
as separate, sign-switched, window-specific terms and never average them into one
'structure' slider." One `Evaluation` carries one `raw_score`, so averaging the
two windows is precisely the thing available and precisely the thing forbidden.
STR(-30:+30) has a sign the evidence supports -- less structure is better --
and is what `raw_score` reports. STR(+31:+90) is measured, travels in `windows`
and in `detail["distal_dg"]`, and is NOT folded into the number.

**The distal window is reported rather than rewarded** because "neutral-to-mildly-
reward moderate structure" has no published magnitude to encode. brief.md:331 says
the same thing from the other side: "we cannot principledly set the relative weight
of 5'-structure vs codon terms past codon 17. Expose both as weighted hypotheses."
A reward with an invented coefficient would look like evidence and be a guess.

**Overlap with B1 is real and is priced in.** B1 measures Kudla's -4..+37 window,
which sits inside STR(-30:+30), and both gate to bacterial expression, so unlike
B8 and B1 these two DO co-apply. B2 is weighted at half of B1's 1.0 to avoid
counting one phenomenon twice; see `weight_provenance`.

**Never put a dG in a byte-exact snapshot** (CLAUDE.md 6), and
`engine_calibration` names the engine every kcal/mol here was measured on.
`registry.check_engine_calibration` refuses the run outright if another engine is
active, because comparing a threshold measured on one energy model against another
model's output succeeds while meaning nothing.

That declaration is INHERITED rather than earned, and the distinction is worth
stating: this rule codes no kcal/mol threshold of its own -- it reports a raw MFE,
`band` is None and `passes` is always True -- so it has no number that could be
mis-calibrated. What it does have is a `raw_score` only comparable to B1's if both
came from one engine, and a weight argued from B1's r = 0.66. Declaring the
calibration therefore constrains more than this rule strictly needs. Left in
deliberately: the failure it prevents is silent, and a rule reporting energies
that a preset weights alongside B1's is exactly where a mixed-engine run would go
unnoticed.
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
    strand_for,
)
from bt5.core.types import Construct, Interval
from bt5.rules.catalog.b1_five_prime import (
    ENGINE_CALIBRATION,
    FIVE_PRIME_UTR_KINDS,
    leader_of,
)

#: STR(-30:+30): 30 bases of leader and 30 of CDS, the A of ATG being +1.
PROXIMAL_UPSTREAM = 30
PROXIMAL_DOWNSTREAM = 30

#: STR(+31:+90): the next 60 bases of CDS, starting 30 in.
DISTAL_START = 30
DISTAL_LENGTH = 60


@register
class StructureWindows:
    id: ClassVar[str] = "b2_structure_windows"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "5' secondary structure in Cambray's proximal window"
    enforcement: ClassVar[Enforcement] = Enforcement.SOFT
    #: brief.md:62 grades B2 "A" -- 244,000 sequences in full factorial is the
    #: largest controlled dataset the catalog cites.
    evidence: ClassVar[Evidence] = Evidence.EVIDENCE_BACKED
    #: Less negative is less structure is better, matching B1's framing of the
    #: same quantity. The DISTAL window's sign is the opposite and is why that
    #: window is not in this number at all.
    direction: ClassVar[Direction] = Direction.HIGHER_IS_BETTER
    unit: ClassVar[str] = "kcal/mol (MFE of STR(-30:+30))"
    band: ClassVar[tuple[float, float] | None] = None
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "Cambray 2018, 244,000 designed sequences in full factorial: phenotype "
            "'dominated by secondary structures and their interactions', localized "
            "to STR(-30:+30) and STR(+31:+90). The source of both windows, and the "
            "reason they are two terms rather than one",
            "https://www.nature.com/articles/nbt.4238",
            2018,
            sign="supports",
        ),
        Citation(
            "Kudla 2009, 154 synonymous GFP variants over a 250-fold range: the "
            "-4..+37 window MFE explained 44% of variance (r = 0.66), 59% in a "
            "second promoter system, while whole-mRNA MFE gave r = 0.16 (n.s.). "
            "That window sits INSIDE STR(-30:+30), which is why B2 is weighted "
            "below B1 rather than beside it",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC3902468/",
            2009,
            sign="supports",
        ),
        Citation(
            "The sign of a structure term is context-dependent and unreconciled: "
            "LinearDesign maximizes global structure for IVT mRNA half-life, and "
            "this benchmark of five schemes over 18 glycoproteins and 90 "
            "small-scale screens found it 'the worst performer ... which produced "
            "the lowest yields' for a DNA transgene in a human cell line. The "
            "direct evidence for refusing to average the two windows into one "
            "structure slider. Sourced here and NOT to Cambray above, which is "
            "2018 and could not have discussed LinearDesign or a 2026 benchmark. "
            "SCOPE: brief.md:333 adds two clauses this page does not carry -- "
            "'2-fold lower normalized MFE' and 'no correlation between MFE and "
            "yield' -- which belong to the preprint it links (bioRxiv "
            "10.64898/2026.03.18.712111) and are deliberately not asserted here, "
            "because that preprint has not been read against this claim",
            "https://proteininnovation.org/2026/03/codon-optimization-native-codon-mammalian-protein-expression/",
            2026,
            sign="qualifies",
        ),
    )
    last_verified: ClassVar[str] = "2026-09-02"
    weight_provenance: ClassVar[str] = (
        "0.5, deliberately HALF of b1_five_prime's 1.0 rather than equal to it. B1 "
        "and B2 both gate to bacterial expression, so unlike B8 they co-apply, and "
        "Kudla's -4..+37 window sits inside Cambray's STR(-30:+30): the two measure "
        "overlapping spans of one phenomenon. Weighting them equally would count 5' "
        "structure twice in the same sum and let it dominate every other objective "
        "by arithmetic rather than by evidence. B1 keeps the larger share because "
        "its window is the one with the published per-variant effect size (r = 0.66 "
        "over 154 measured variants); B2 earns its half from the larger and better "
        "controlled dataset behind the window's LOCATION (244,000 sequences, full "
        "factorial). NOTE no shipped preset weights 2.B2, so this default reaches a "
        "user only outside a preset; adding it belongs to the score lane. "
        "steering_weight is 0.0 because a folding free energy is not decidable from "
        "a bounded codon suffix, which is b1's reason and unchanged here."
    )
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.5
    steering_weight: ClassVar[float] = 0.0
    #: The finding is a whole 60 nt window, not a base within it.
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.WINDOW_MINUS_1
    #: SINGLE_PASS: unlike a splice donor, flattening structure in one window
    #: cannot manufacture structure in it again -- MFE is a property of the whole
    #: window, and a repair that raises it has raised it (CLAUDE.md 3.6).
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    #: Folding is the expensive one in this catalog; it drives null sampling.
    cost_class: ClassVar[str] = "expensive"
    conflicts_with: ClassVar[tuple[str, ...]] = ()
    brief_ref: ClassVar[str] = "2.B2"
    #: Every kcal/mol here is measured on this engine and parameter set. Applying
    #: it to another engine's output fails silently, so the registry refuses.
    engine_calibration: ClassVar[str | None] = ENGINE_CALIBRATION
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "report_distal": {
                "type": "boolean",
                "default": True,
                "description": (
                    "Measure STR(+31:+90) and report it alongside the score. It is "
                    "never added to the score: brief.md:333 forbids averaging the "
                    "two windows, and the distal window's 'mildly reward' has no "
                    "published magnitude to encode. Turn off to skip one fold."
                ),
            },
        },
    }

    def __init__(self, report_distal: bool = True) -> None:
        self.report_distal = report_distal

    def gate(self, slot: ContextSlot) -> bool:
        """Bacteria, per brief.md:62's Context column -- B1's gate exactly.

        The eukaryotic cap-proximal story is B11's row with its own thresholds and
        a weaker evidence grade; it is not this rule with a different host.
        """
        return slot.modality is Modality.BACTERIAL_EXPRESSION

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """None. A folding free energy is not decidable from a bounded suffix."""
        return None

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        slots = [s for s in ctx.active_slots if self.gate(s)]
        if not slots:
            return self._unavailable(
                c,
                "no bacterial expression slot in this context. Cambray's windows "
                "were measured in E. coli; the eukaryotic cap-proximal rule is "
                "B11's row, with its own thresholds",
            )
        if svc.fold is None:
            return self._unavailable(
                c,
                "no folding engine installed, so neither structure window can be "
                "measured; install ViennaRNA to enable it",
            )

        editable = sorted(c.editable)
        if not editable:
            return self._unavailable(c, "no designable CDS to take 5' windows from")

        measured: list[tuple[Interval, float]] = []
        proximal: list[float] = []
        distal_dg: float | None = None
        for slot in slots:
            strand = strand_for(ctx, slot)
            cds = editable[0] if strand == 1 else editable[-1]
            oriented = Interval(cds.start, cds.end, strand)
            window = self._proximal(oriented, c)
            if window is None:
                return self._unavailable(
                    c,
                    f"the CDS starts too close to the end of this linear construct "
                    f"for STR(-{PROXIMAL_UPSTREAM}:+{PROXIMAL_DOWNSTREAM}) to fit",
                )
            leader = leader_of(window, PROXIMAL_UPSTREAM)
            if not self._annotated_leader(c, leader):
                return self._unavailable(
                    c,
                    f"no annotated 5'UTR covering the {PROXIMAL_UPSTREAM} bases "
                    f"upstream of the start codon. STR(-{PROXIMAL_UPSTREAM}:"
                    f"+{PROXIMAL_DOWNSTREAM}) spans the UTR/CDS junction, and "
                    f"unannotated sequence there may be promoter rather than "
                    f"leader -- folding it would report a molecule that is never "
                    f"transcribed. Annotate the 5'UTR in your map and re-run",
                    interval=leader,
                )
            dg = svc.fold.mfe_window(c.sequence, window).dg_kcal_mol
            proximal.append(dg)
            measured.append((window, dg))

            if self.report_distal:
                far = self._distal(oriented, c)
                if far is not None:
                    far_dg = svc.fold.mfe_window(c.sequence, far).dg_kcal_mol
                    measured.append((far, far_dg))
                    distal_dg = far_dg if distal_dg is None else min(distal_dg, far_dg)

        # The MOST NEGATIVE proximal window across gated slots binds: the most
        # structured 5' end is the finding, and averaging would let a permissive
        # slot hide it. This is B1's aggregation and the same argument.
        worst = min(proximal)
        detail: dict[str, float | str] = {
            "proximal_dg": worst,
            "proximal_window": f"STR(-{PROXIMAL_UPSTREAM}:+{PROXIMAL_DOWNSTREAM})",
            "calibration": ENGINE_CALIBRATION,
        }
        if distal_dg is not None:
            # Reported, never scored. brief.md:333: never average the two windows.
            detail["distal_dg"] = distal_dg
            detail["distal_window"] = f"STR(+{DISTAL_START + 1}:+{DISTAL_START + DISTAL_LENGTH})"

        return Evaluation(
            spec_id=self.id,
            # An objective, not a constraint: there is no published cutoff for
            # these windows, only a monotone relationship over Cambray's factorial.
            passes=True,
            raw_score=worst,
            breaches=(),
            windows=tuple(measured),
            n_evaluated=PROXIMAL_UPSTREAM + PROXIMAL_DOWNSTREAM,
        )

    # -- internals ----------------------------------------------------------

    def _proximal(self, cds: Interval, c: Construct) -> Interval | None:
        """STR(-30:+30), or None when it does not fit a linear construct.

        None rather than a clamped window, for B1's reason: a window silently
        shortened to what fits is a DIFFERENT window, and its MFE is not
        comparable to one that had its whole leader.
        """
        return self._window(cds, c, PROXIMAL_UPSTREAM, PROXIMAL_DOWNSTREAM)

    def _distal(self, cds: Interval, c: Construct) -> Interval | None:
        """STR(+31:+90). Entirely inside the CDS, so no leader is needed.

        None when the CDS is shorter than 90 nt -- a short ORF simply has no
        distal window, which is not a failure of the rule and must not make the
        proximal measurement unavailable.
        """
        if cds.length < DISTAL_START + DISTAL_LENGTH:
            return None
        start = (
            cds.start + DISTAL_START if cds.strand == 1 else cds.end - DISTAL_START - DISTAL_LENGTH
        )
        # Same origin-spanning hazard as `_window`; same reason it is invisible.
        if start >= c.length or start < 0:
            if not c.is_circular:
                return None
            start %= c.length
        return Interval(start, start + DISTAL_LENGTH, cds.strand)

    def _window(
        self, cds: Interval, c: Construct, upstream: int, downstream: int
    ) -> Interval | None:
        """`bt5.structure.windows.five_prime_window`'s arithmetic, as b1 has it.

        Rules receive Services and never import another lane, and window geometry
        is not on the Services protocol.
        """
        if cds.length < downstream:
            return None
        span = upstream + downstream
        start = cds.start - upstream if cds.strand == 1 else cds.end - downstream
        if start < 0:
            if not c.is_circular:
                return None
            start += c.length
        elif start >= c.length:
            # A CDS spanning the origin has `cds.end > c.length`, so on the minus
            # strand the window can START past the end. `Construct.slice` handles
            # an end beyond the sequence but not a start beyond it -- `seq[start:]`
            # is silently "" and the window returns a DIFFERENT LENGTH, folded on
            # different bases, with a plausible dG and passes=True. Normalising is
            # what keeps that from being a wrong answer nobody can see.
            if not c.is_circular:
                return None
            start %= c.length
        if not c.is_circular and start + span > c.length:
            return None
        return Interval(start, start + span, cds.strand)

    def _annotated_leader(self, c: Construct, leader: Interval) -> bool:
        return any(
            f.interval.contains(leader, c.length, c.is_circular)
            for f in c.features
            if f.kind in FIVE_PRIME_UTR_KINDS
        )

    def _unavailable(
        self, c: Construct, reason: str, *, interval: Interval | None = None
    ) -> Evaluation:
        """NaN plus a breach carrying the reason -- B1's pattern and its argument.

        NaN rather than 0.0 because 0 kcal/mol is a real and maximally GOOD value
        for this quantity: a completely unstructured 5' end. Reporting it for "we
        could not measure" would put the objective at the top of the ranking
        exactly when it is unknown.
        """
        where = interval or (sorted(c.editable)[0] if c.editable else Interval(0, 1))
        return Evaluation(
            spec_id=self.id,
            passes=True,
            raw_score=float("nan"),
            breaches=(
                Breach(
                    spec_id=self.id,
                    interval=where,
                    magnitude=0.0,
                    message=f"5' structure window objective unavailable: {reason}",
                    # The missing input is an annotation, a host or the engine
                    # install -- never a codon the solver could have chosen.
                    fixable_by_codon_choice=False,
                    detail={
                        "unavailable_reason": reason,
                        "calibration": ENGINE_CALIBRATION,
                    },
                ),
            ),
            n_evaluated=0,
        )
