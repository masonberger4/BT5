"""D3 -- cryptic splice donors and acceptors, on the strand that is transcribed.

This is the rule `CLAUDE.md` §3.6 was written for, and the first in the catalog to
declare `RepairPolicy.FIXED_POINT`. The reason is in `RepairPolicy`'s own docstring
(`core/spec.py:106-117`) and in `brief.md:104`: point-mutating one cryptic donor
activates cryptic donors nearby (the A2UCOE case), so a single-pass "remove motifs"
step ships a construct whose donors were removed INTO new donors, and the validator
passes it because the specific 9-mer it was told to avoid is gone.

What this rule does NOT do, and why:

  - **It does not iterate.** `evaluate()` is a pure, idempotent detector; the solver
    owns the loop. `RuleSet.breach_finder` wraps the catalog into one `BreachFinder`
    and `repair()` re-calls it on every candidate, up to 256 per iteration. The fixed
    point is over the set of breaches this rule emits, and it is reached when the rule
    stops emitting them -- which is what declaring FIXED_POINT buys, and it is the
    whole of what it buys. A `while` loop in here would run inside every one of those
    256 evaluations and would still not be the fixed point.

  - **It is not a lattice rule, and returns no `forbidden`.** `LatticeTerms.forbidden`
    is closed under reverse complement by the solver, which is right for a restriction
    site and wrong here for the same reason it is wrong in `d4_internal_polya`:
    `GGTAAG` is a donor on the transcribed strand and its complement `CTTACC` is not a
    donor on anything anyone transcribes. Forbidding both refuses designs over a site
    that cannot fire. So this rule reads `strand_for(ctx, slot)` per slot and repairs,
    rather than being guaranteed by construction.

  - **It does not score with MaxEntScan unless it is given a model.** `brief.md:288`
    records MaxEntScan redistribution as ambiguous and recommends training a
    max-entropy model on GENCODE instead; no such model ships here. So the 9-mer/23-mer
    scored path reports its objective UNAVAILABLE rather than returning a plausible
    number, and the motif scan -- which needs no licensed data -- is what actually runs.
    A stub returning a believable bit score is the one thing the honesty apparatus
    exists to prevent.

What the evidence does and does not support (`brief.md:321`): no published dataset
cleanly quantifies the titer cost of a cryptic splice donor inside a therapeutic ORF --
every quantified case is confounded. The V5 sub-case is the exception and is graded
evidence A: the standard V5 nucleotide encoding contains `G|GTAAG` and spliced in 17/17
genes tested, with 13/17 randomly chosen genes showing aberrant splicing from
vector/tag context. That is why V5 escalates and a generic cryptic donor does not.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import ClassVar

from bt5.core.context import ContextSlot, DesignContext, HostId, Modality
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
from bt5.core.types import Construct, Interval, Strand, reverse_complement

#: The MaxEntScan 5' donor window: 3 exonic + 6 intronic (`brief.md:102`). This is
#: also the reported breach span, which is load-bearing rather than cosmetic --
#: `_breach_key` is (spec_id, start, end), so anchoring the span on the invariant GT
#: is what makes "the repair weakened this donor but did not remove it" the SAME
#: breach across iterations instead of a new one.
DONOR_EXONIC = 3
DONOR_INTRONIC = 6
DONOR_WINDOW = DONOR_EXONIC + DONOR_INTRONIC

#: `brief.md:102` literal blacklist. Both are 6-mers, and `brief.md:90`'s
#: expected-occurrence budget makes a 6-mer hard "only if strongly evidenced". The V5
#: result (17/17, evidence A, `brief.md:105`) is that evidence.
LITERAL_DONORS: tuple[str, ...] = ("GGTAAG", "GGTGAG")

#: `brief.md:102` `AN|GT(A/G)AG` -- the exon|intron consensus, written here without the
#: boundary bar.
#:
#: SEVEN CHARACTERS, BUT NOT A 7-MER, and `brief.md:90`'s budget is about expected
#: OCCURRENCES rather than about length. Two positions are degenerate -- N is 4-fold and
#: R is 2-fold -- so P = (1/4)^5 x 1 x (1/2) = 1/2048, against 1/16384 for a true 7-mer.
#: Effective k = log4(2048) = 5.5, and measured on random 50% GC DNA it fires
#: **2.85 times per 1.5 kb**, which is brief.md:90's own figure for a 5-mer (2.9), not
#: its 0.18 for a 7-mer.
#:
#: brief.md:90 puts a <=5-mer at SOFT ONLY, so a consensus-only match is reported at the
#: FLAG tier. Only the literal 6-mers (0.70/1.5 kb measured, against the brief's 0.73)
#: and V5-overlapping hits reach the hard tier. Pinned by
#: `TestDonorTiers::test_the_degenerate_consensus_is_not_a_seven_mer`.
DONOR_CONSENSUS = "ANGTRAG"

#: `brief.md:102` coarse screen `GTNNG`. 5 nt, and `brief.md:90` is explicit that a
#: <=5-mer is SOFT ONLY -- ~2.9 expected occurrences per 1.5 kb. Reported, never hard,
#: and never placed in `forbidden`.
DONOR_COARSE = "GTNNG"

#: `brief.md:103` acceptor context. A >=10-nt window >=80% pyrimidine within 5-40 nt
#: upstream of the AG, plus a branch-point-like `YTNAY` 18-40 nt upstream.
PPT_MIN_LEN = 10
PPT_MIN_PYRIMIDINE = 0.80
PPT_FROM, PPT_TO = 5, 40
BRANCH_POINT = "YTNAY"
BRANCH_FROM, BRANCH_TO = 18, 40

#: The V5 epitope tag, as PROTEIN (`brief.md:190`). Matched on the translation, never
#: on a nucleotide string: the liability is the *standard encoding*, so a nucleotide
#: match would both miss every other encoding and match nothing once repair has
#: recoded it -- the rule would report success exactly when it stopped working.
V5_PEPTIDE = "GKPIPNPLLGLDST"

#: Magnitudes. Graded on purpose: `_accepts` (`solver/repair.py:270-312`) uses breach
#: COUNT as the cross-rule currency and magnitude only to break ties within one rule,
#: so a repair that weakens a donor without removing it is only ever accepted -- and
#: only ever visible as progress -- if the tiers differ.
MAG_STRONG = 3.0  # literal blacklist, the 7-mer consensus, or inside a V5 tag
MAG_FLAGGED = 1.0  # MaxEntScan > 3 bits, when a model is available
MAG_COARSE = 0.2  # GTNNG only

#: `brief.md:102`: "flag >3 bits, hard-constrain >6-8". The hard band is stated as a
#: range; the rule takes its conservative end, so a site is only hard-constrained on a
#: score the brief unambiguously calls hard.
MAXENT_FLAG_BITS = 3.0
MAXENT_HARD_BITS = 8.0

#: `brief.md:223`, the cryptic-splice-donor row: warn | warn | HARD (titer + safety) |
#: warn | n/a | HARD (fusion transcripts).
HARD_MODALITIES: frozenset[Modality] = frozenset({Modality.LENTIVIRAL, Modality.GENOME_INTEGRATED})

#: `brief.md:101`: eukaryotic Pol II contexts only -- irrelevant for E. coli, yeast
#: heterologous CDS, and IVT mRNA.
OFF_MODALITIES: frozenset[Modality] = frozenset({Modality.BACTERIAL_EXPRESSION, Modality.IVT_MRNA})
OFF_HOSTS: frozenset[HostId] = frozenset(
    {HostId.E_COLI_K12, HostId.E_COLI_BL21, HostId.S_CEREVISIAE, HostId.P_PASTORIS}
)

_PYRIMIDINES = frozenset("CT")
_IUPAC: Mapping[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "N": frozenset("ACGT"),
}


def _matches(text: str, pattern: str) -> bool:
    """IUPAC match. Local to this rule on purpose -- `verify.py`'s expander is the
    ORACLE's, and `tests/data_integrity/test_oracle_independence.py` exists to keep it
    off every lane's code path (see `docs/decisions/2026-09-01-expand-forbidden-iupac.md`,
    which rejected sharing it for the solver for the same reason)."""
    if len(text) != len(pattern):
        return False
    return all(base in _IUPAC[code] for base, code in zip(text, pattern, strict=True))


class _Donor:
    """One GT dinucleotide plus the verdict on its context."""

    __slots__ = ("magnitude", "position", "reason")

    def __init__(self, position: int, magnitude: float, reason: str) -> None:
        self.position = position
        self.magnitude = magnitude
        self.reason = reason


@register
class Splicing:
    id: ClassVar[str] = "d3_splicing"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "Cryptic splice donors and acceptors on the transcribed strand"
    #: The FLOOR. `enforcement_for` escalates on the two modalities `brief.md:223`
    #: grades hard; everywhere else the row says "warn", which is SOFT.
    enforcement: ClassVar[Enforcement] = Enforcement.SOFT
    evidence: ClassVar[Evidence] = Evidence.EVIDENCE_BACKED
    direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
    unit: ClassVar[str] = "weighted splice sites"
    band: ClassVar[tuple[float, float] | None] = None
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "The standard V5 nucleotide encoding carries a G|GTAAG donor and spliced in "
            "17/17 genes tested; 13/17 randomly chosen genes showed aberrant splicing "
            "from vector/tag context -- the evidence-A case for recoding tags",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC9379414/",
            2022,
            sign="supports",
        ),
        Citation(
            "A2UCOE: point-mutating one cryptic donor activates cryptic donors nearby, "
            "which is why removal must iterate to a fixed point rather than pass once",
            "https://jvi.asm.org/content/86/9088",
            2012,
            sign="supports",
        ),
        Citation(
            "MaxEntScan (Yeo & Burge 2004) is the 9-mer/23-mer maximum-entropy model the "
            "bit thresholds are stated against; redistribution is ambiguous, so this rule "
            "scores only when a model is supplied and reports unavailable otherwise",
            "https://www.gencodegenes.org/",
            2004,
            sign="qualifies",
        ),
        Citation(
            "No published dataset cleanly quantifies the titer cost of a cryptic splice "
            "donor inside a therapeutic ORF -- every quantified case is confounded, so "
            "the generic donor tiers are a liability ranking and not a predicted cost",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC9379414/",
            2022,
            sign="qualifies",
        ),
    )
    last_verified: ClassVar[str] = "2026-09-02"
    weight_provenance: ClassVar[str] = (
        "Mid-high for a soft rule, and deliberately below d4_internal_polya's 0.7. The "
        "mechanism is not in doubt and the V5 sub-case is evidence A (17/17), but "
        "brief.md:321 is explicit that no clean dataset quantifies what a cryptic donor "
        "inside a therapeutic ORF actually costs -- every quantified case is confounded. "
        "A weight tuned as if the cost were measured would be a predicted quantity in "
        "disguise. Where the cost IS attested -- lentiviral titer and safety, and fusion "
        "transcripts from an integrated cassette -- enforcement_for takes the rule out of "
        "the weighted sum entirely rather than raising this number, which is the "
        "mechanism CLAUDE.md 3.5 requires: a hard constraint is never a heavy weight."
    )
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.5
    steering_weight: ClassVar[float] = 0.3
    #: WINDOW_MINUS_1, not MOTIF_LEN_MINUS_1, and the instance attribute `window`
    #: below is why. `solver/catalog.py:241` hard-codes `motif_len` to 6 for every
    #: rule, so MOTIF_LEN_MINUS_1 would widen the repair window by 5 when the donor
    #: geometry needs 8. `catalog.py:236` reads `window` off the spec, so this is the
    #: only per-rule width knob a rule author can reach without editing the solver.
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.WINDOW_MINUS_1
    #: CLAUDE.md 3.6. See the module docstring for what the fixed point is over.
    repair: ClassVar[RepairPolicy] = RepairPolicy.FIXED_POINT
    cost_class: ClassVar[str] = "cheap"
    conflicts_with: ClassVar[tuple[str, ...]] = ()
    brief_ref: ClassVar[str] = "2.D3"
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "scan_acceptors": {
                "type": "boolean",
                "default": False,
                "description": (
                    "Also report 3' acceptors carrying a pyrimidine tract and a branch "
                    "point. Off by default: brief.md:103 only FLAGS acceptors and states "
                    "no hard cutoff, but this rule's enforcement class is set by its "
                    "donor half, so an acceptor breach emitted here inherits HARD_REPAIR "
                    "and the solver starts chasing sites the brief never authorised "
                    "constraining. A REPORT_ONLY acceptor rule is the right home."
                ),
            },
            "report_coarse": {
                "type": "boolean",
                "default": False,
                "description": (
                    "Report the coarse GTNNG screen. Off by default: at 5 nt it is below "
                    "the chance floor -- P(GTNNG) = 1/64 per position is ~78 hits on a "
                    "5 kb plasmid of random sequence -- so it reports noise rather than "
                    "findings. brief.md:90 already grades a <=5-mer soft-only."
                ),
            },
        },
    }

    def __init__(self, scan_acceptors: bool = False, report_coarse: bool = False) -> None:
        self.scan_acceptors = scan_acceptors
        self.report_coarse = report_coarse
        #: Read by `solver/catalog.py:236` into `RulePolicy.window`. The donor window
        #: is the geometry that matters: a 9-mer can only be created by bases within 8
        #: of it. Must stay an `int` -- catalog.py falls back to 50 for anything else.
        self.window: int = DONOR_WINDOW
        #: No MaxEntScan model ships (see the module docstring). Injectable so the
        #: scored path is real code rather than a comment, and so the honest
        #: "unavailable" is what a caller without a model actually gets.
        self.maxent: object | None = None

    def gate(self, slot: ContextSlot) -> bool:
        """Eukaryotic Pol II contexts only (`brief.md:101`).

        Off for bacteria and IVT mRNA by modality, and off for the yeast hosts, whose
        heterologous CDS is the case `brief.md:208` marks "Off (heterologous CDS)".
        """
        return slot.modality not in OFF_MODALITIES and slot.host not in OFF_HOSTS

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        """`brief.md:223`: HARD for lentiviral (titer + safety) and genome-integrated
        (fusion transcripts); warn -- i.e. SOFT -- for the plasmid modalities and AAV."""
        if slot.modality in HARD_MODALITIES:
            return Enforcement.HARD_REPAIR
        return Enforcement.SOFT

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """Deliberately None; see the module docstring. The solver closes `forbidden`
        under reverse complement, which would forbid CTTACC because GGTAAG is a donor
        and refuse designs over a site nothing transcribes."""
        return None

    # -- geometry -------------------------------------------------------------

    def _linearized(self, c: Construct, strand: Strand) -> tuple[str, int]:
        """The construct as one string on `strand`, with the offset of its middle copy.

        A circular construct is tripled, exactly as `Construct.tripled()` does it, so
        an origin-spanning donor and its 40 nt of upstream acceptor context are
        ordinary substring lookups rather than a special case each site has to
        remember. `reverse_complement(S * 3) == reverse_complement(S) * 3`, so the
        middle copy stays the middle copy on either strand.
        """
        text, offset = c.tripled()
        if strand == -1:
            text = reverse_complement(text)
        return text, offset

    def _to_construct_coords(
        self, index: int, length: int, n: int, strand: Strand, text_len: int, offset: int
    ) -> int:
        """Map a hit in the linearized text back to a forward construct coordinate.

        The modulo is not defensive. Mapping a minus-strand hit back goes NEGATIVE
        when the site crosses the origin, and clamping it to zero would move the
        breach to position 0 -- an interval whose bases are not the site it names, on
        exactly the reverse-oriented cassette this rule exists for.
        """
        start = index - offset if strand == 1 else (text_len - index - length) - offset
        return start % n

    # -- donors ---------------------------------------------------------------

    @staticmethod
    def _at(text: str, lo: int, hi: int) -> str:
        """`text[lo:hi]`, but empty rather than wrong when it runs off either end.

        Plain slicing is the trap here: a linear construct puts the first GT at index
        0, and `text[p - 2 : p + 5]` then reads from the END of the string and matches
        a consensus assembled out of two ends that are not adjacent to each other.
        """
        if lo < 0 or hi > len(text):
            return ""
        return text[lo:hi]

    def _classify_donor(self, text: str, p: int) -> _Donor | None:
        """Verdict on the GT at `text[p:p+2]`, from context alone.

        `p` is the invariant GT; every pattern is expressed as an offset from it so
        the three shapes in `brief.md:102` cannot drift out of register with each
        other. Returns the STRONGEST verdict, so a site matching both the literal
        blacklist and the consensus is one breach, not two.
        """
        if self._at(text, p - 1, p + 5) in LITERAL_DONORS:
            literal = text[p - 1 : p + 5]
            return _Donor(p, MAG_STRONG, f"literal donor {literal!r}")

        consensus = self._at(text, p - 2, p + 5)
        if consensus and _matches(consensus, DONOR_CONSENSUS):
            # MAG_FLAGGED, not MAG_STRONG: the consensus is 7 characters but only a
            # 5.5-mer by occurrence, and brief.md:90 puts a <=5-mer at soft only.
            return _Donor(p, MAG_FLAGGED, f"donor consensus AN|GT(A/G)AG as {consensus!r}")

        if self.report_coarse:
            coarse = self._at(text, p, p + 5)
            if coarse and _matches(coarse, DONOR_COARSE):
                return _Donor(p, MAG_COARSE, f"coarse GTNNG screen as {coarse!r}")
        return None

    def _score_donor(self, text: str, p: int) -> float | None:
        """MaxEntScan bits for the 9-mer at this GT, or None when no model is loaded."""
        if self.maxent is None:
            return None
        nine = self._at(text, p - DONOR_EXONIC, p + DONOR_INTRONIC)
        if not nine:
            return None
        return float(self.maxent.score5(nine))  # type: ignore[attr-defined]

    # -- acceptors ------------------------------------------------------------

    def _acceptor_context(self, text: str, p: int) -> str | None:
        """`brief.md:103`: an AG at `p` with BOTH a pyrimidine tract and a branch point.

        Both, not either. A lone AG dinucleotide occurs every ~16 bases and reporting
        one would bury every real finding; the brief pairs the two context features
        for exactly that reason.
        """
        tract = text[max(0, p - PPT_TO) : max(0, p - PPT_FROM)]
        has_tract = any(
            sum(base in _PYRIMIDINES for base in tract[i : i + PPT_MIN_LEN])
            >= PPT_MIN_PYRIMIDINE * PPT_MIN_LEN
            for i in range(max(0, len(tract) - PPT_MIN_LEN + 1))
        )
        if not has_tract:
            return None

        upstream = text[max(0, p - BRANCH_TO) : max(0, p - BRANCH_FROM)]
        width = len(BRANCH_POINT)
        branch = next(
            (
                upstream[i : i + width]
                for i in range(max(0, len(upstream) - width + 1))
                if _matches(upstream[i : i + width], BRANCH_POINT)
            ),
            None,
        )
        if branch is None:
            return None
        return f"pyrimidine tract >=80% within {PPT_FROM}-{PPT_TO} nt and branch point {branch!r}"

    # -- V5 -------------------------------------------------------------------

    def _v5_spans(self, c: Construct) -> list[Interval]:
        """Construct intervals covering any V5 tag, found on the TRANSLATION.

        `brief.md:105` is about the standard V5 *encoding*, so matching a nucleotide
        string would miss every other encoding and -- worse -- would stop matching the
        moment repair recoded it, making the rule report success exactly when it
        stopped working. The peptide is the stable identity.
        """
        spans: list[Interval] = []
        for unit in c.translation_units:
            start = unit.protein.find(V5_PEPTIDE)
            while start != -1:
                codons = unit.codon_map[start : start + len(V5_PEPTIDE)]
                if codons:
                    spans.append(Interval(codons[0].start, codons[-1].end))
                start = unit.protein.find(V5_PEPTIDE, start + 1)
        return spans

    # -- evaluation -----------------------------------------------------------

    def _scan(self, c: Construct, strand: Strand, role: str, v5: list[Interval]) -> list[Breach]:
        n = c.length
        text, offset = self._linearized(c, strand)
        breaches: list[Breach] = []
        # `str.find`, not `range(len)`. This runs once per CANDIDATE inside the repair
        # search -- up to 256 per iteration (solver/catalog.py:200-206) -- so a
        # per-position Python loop over a tripled 5 kb plasmid is 15k iterations paid
        # tens of thousands of times, which is minutes of wall clock, not milliseconds.
        # Only sites whose anchor lies in the middle copy are real; the flanking copies
        # exist so their context is available, not to be reported twice.
        limit = offset + n
        p = text.find("GT", offset)
        while p != -1 and p < limit:
            # A donor needs its whole 9-mer to be real sequence. On a circular
            # construct the tripling guarantees that; on a linear one, a GT within 3
            # of the start or 6 of the end genuinely has no context to score, and
            # inventing it by clamping would report a site assembled from an edge.
            if p - DONOR_EXONIC < 0 or p + DONOR_INTRONIC > len(text):
                p = text.find("GT", p + 1)
                continue
            hit = self._classify_donor(text, p)
            bits = self._score_donor(text, p)
            if bits is not None and bits > MAXENT_FLAG_BITS:
                magnitude = MAG_STRONG if bits > MAXENT_HARD_BITS else MAG_FLAGGED
                if hit is None or magnitude > hit.magnitude:
                    hit = _Donor(p, magnitude, f"MaxEntScan 5' donor {bits:.2f} bits")
            if hit is None:
                p = text.find("GT", p + 1)
                continue

            lo = self._to_construct_coords(
                p - DONOR_EXONIC, DONOR_WINDOW, n, strand, len(text), offset
            )
            iv = Interval(lo, lo + DONOR_WINDOW, strand)
            in_v5 = any(self._overlaps(iv, span, n) for span in v5)
            magnitude = MAG_STRONG if in_v5 else hit.magnitude
            breaches.append(
                Breach(
                    spec_id=self.id,
                    interval=iv,
                    magnitude=magnitude,
                    message=(
                        f"5' splice donor: {hit.reason} on the transcribed strand "
                        f"({'+' if strand == 1 else '-'}) at {iv.start}"
                        + (
                            "; inside a V5 tag, whose standard encoding spliced in 17/17 "
                            "genes tested -- recode the tag rather than accepting it"
                            if in_v5
                            else ""
                        )
                    ),
                    slot_role=role,
                    fixable_by_codon_choice=c.overlaps_editable(iv),
                    detail={
                        "kind": "donor",
                        "reason": hit.reason,
                        "strand": float(strand),
                        "v5_tag": "yes" if in_v5 else "no",
                        **({"maxent_bits": bits} if bits is not None else {}),
                    },
                )
            )
            p = text.find("GT", p + 1)

        if self.scan_acceptors:
            q = text.find("AG", offset)
            while q != -1 and q < limit:
                p, q = q, text.find("AG", q + 1)
                if p + 2 > len(text):
                    continue
                reason = self._acceptor_context(text, p)
                if reason is None:
                    continue
                lo = self._to_construct_coords(p, 2, n, strand, len(text), offset)
                iv = Interval(lo, lo + 2, strand)
                breaches.append(
                    Breach(
                        spec_id=self.id,
                        interval=iv,
                        # Never MAG_STRONG: brief.md:103 gives a flag threshold for
                        # acceptors and no hard one, and this rule does not invent the
                        # cutoff the brief declined to state.
                        magnitude=MAG_FLAGGED,
                        message=(
                            f"3' splice acceptor: AG at {iv.start} with {reason}, on the "
                            f"transcribed strand ({'+' if strand == 1 else '-'})"
                        ),
                        slot_role=role,
                        fixable_by_codon_choice=c.overlaps_editable(iv),
                        detail={"kind": "acceptor", "reason": reason, "strand": float(strand)},
                    )
                )
        return breaches

    @staticmethod
    def _overlaps(a: Interval, b: Interval, n: int) -> bool:
        """Overlap of two possibly-wrapping intervals, by expanding to residues.

        `Interval.overlaps` requires the construct length and circularity and is the
        right tool, but both of these are construct-coordinate spans of bounded width,
        so the explicit residue sets are cheaper to read and cannot get the wrap wrong.
        """
        ra = {i % n for i in range(a.start, a.end)}
        rb = {i % n for i in range(b.start, b.end)}
        return not ra.isdisjoint(rb)

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        """One pass per gated slot, on that slot's own transcribed strand.

        Two slots can disagree about which strand is transcribed -- a reverse-oriented
        cassette is the case `brief.md:244` names -- so breaches carry `slot_role` and
        are deduplicated on (span, strand) rather than collapsed.
        """
        v5 = self._v5_spans(c)
        breaches: list[Breach] = []
        seen: set[tuple[int, int, int, str]] = set()
        gated = False
        for slot in ctx.active_slots:
            if not self.gate(slot):
                continue
            gated = True
            strand = strand_for(ctx, slot)
            for breach in self._scan(c, strand, slot.role, v5):
                key = (
                    breach.interval.start,
                    breach.interval.end,
                    int(strand),
                    str(breach.detail.get("kind", "")),
                )
                if key in seen:
                    continue
                seen.add(key)
                breaches.append(breach)

        if gated and self.maxent is None:
            breaches.append(self._scored_path_unavailable(c))

        return Evaluation(
            spec_id=self.id,
            passes=not any(b.magnitude >= MAG_STRONG for b in breaches),
            raw_score=float(sum(b.magnitude for b in breaches)),
            breaches=tuple(breaches),
            n_evaluated=c.length if gated else 0,
        )

    def _scored_path_unavailable(self, c: Construct) -> Breach:
        """The MaxEntScan half did not run, said out loud.

        c1_cai's pattern (`c1_cai.py:490-525`) and the same argument: a rule that
        silently reports only what it could compute reads as a clean scan of
        everything it claims to cover. Magnitude 0.0 and not fixable, so it lands on
        `RepairOutcome.advisory` and never becomes a target the solver chases.
        """
        where = sorted(c.editable)[0] if c.editable else Interval(0, min(1, c.length) or 1)
        reason = (
            "no MaxEntScan model is loaded, so the 9-mer donor and 23-mer acceptor bit "
            "scores were not computed; the motif and context scans below did run"
        )
        return Breach(
            spec_id=self.id,
            interval=where,
            magnitude=0.0,
            message=f"splice scoring unavailable: {reason}",
            fixable_by_codon_choice=False,
            detail={"unavailable_reason": reason},
        )
