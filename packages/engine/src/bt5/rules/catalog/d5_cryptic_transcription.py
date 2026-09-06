"""D5 -- bacterial cryptic transcription, and the antisense hexamer nobody else forbids.

`brief.md:110-112`. A codon optimizer that is blind to promoter architecture can
build a working *E. coli* promoter inside a CDS by synonymous substitution alone.

The dengue-2 case (PMC3069047) is why the hazard is taken seriously: a cryptic
promoter in a flavivirus cDNA drove ~10^6.6 mRNA copies/ug uninduced. It is NOT
a test vector for this rule, and the docstring says so where the rule cannot see
it (`_scored_path_unavailable` and the `AT_TRACT_GAP` note below). The brief's
summary of it also does not survive the primary source -- see the citation label.

Five parts, and they do not all have the same shape:

  (a) hexamer within 1 mismatch of `TTGACA`, a 15-19 bp spacer, then a hexamer
      within 1 mismatch of `TATAAT`
  (b) the extended -10, `TGnTATAAT`, which needs no -35 at all
  (c) the AT-tract rule: `TATTTAT` or `AATTT` sitting at promoter positions
      -23..-17, i.e. with EXACTLY FOUR bases between the tract and a -10-like
      hexamer -- 103/103 randomized appropriately positioned AT-tracts were
      active. See `AT_TRACT_GAP`: the brief's "5 bp upstream" is a position
      offset, not a gap width, and reading it as a gap is off by one.
  (d) the sigma-38 variant -10, `TATACT`
  (e) the antisense promoter, novel and implemented by no competitor: forbid
      `ATTATA` (the reverse complement of `TATAAT`) on the sense strand. Only
      4.76% of 484,741 natural *E. coli* CDSs contain it -- 2-3x
      under-represented -- yet 77.28% of clean CDSs can silently acquire it by
      synonymous substitution. Nature avoids it; naive optimization walks in.

WHY THIS RULE CARRIES TWO ENFORCEMENT MECHANISMS
------------------------------------------------
(b), (d) and (e) are finite string sets, so the Tier-A automaton can make
them unreachable by construction. (a) and (c) are GEOMETRY -- two elements at a
declared distance -- and no finite forbidden-string set expresses "15 to 19 bases
apart". Enumerating it would be 19 x 5 x 19 x 4^spacer patterns.

So the ClassVar `enforcement` is HARD_LATTICE and `enforcement_for` returns
HARD_REPAIR, and those are two independent switches rather than a contradiction:
`solver/catalog.py:149` reads the CLASSVAR when it collects forbidden patterns
for the automaton (and those same patterns reach the independent validator),
while `solver/catalog.py:302` and `_enforcement_of` read `enforcement_for` when
they route breaches to Tier-B. `d4_internal_polya` already ships the same split
in the other direction -- a SOFT ClassVar with a HARD_REPAIR `enforcement_for`.

The cost of that choice, stated: `test_only_lattice_rules_skip_steering` forces
`steering_weight = 0.0` on any HARD_LATTICE rule, so the DP gets NO nudge away
from the (a)/(c) geometry and repair does all of that work. Paying it buys the
automaton guarantee plus the oracle's own proof for the three motif parts, which
is the half with the disaster anchor and the novel finding.

WHAT LISTING `ATTATA` ALSO FORBIDS -- deliberate, not an accident
-----------------------------------------------------------------
`ATTATA` IS the reverse complement of `TATAAT`, and both the solver and
`verify.find_motifs` close the pattern set under reverse complement. So listing
`ATTATA` necessarily forbids a bare sense-strand `TATAAT` too. That is a real
constraint beyond the literal words of (e), and it is kept rather than worked
around: a sense-strand `TATAAT` is the -10 box that (a), (b) and (c) are each
built around, and the expected-occurrence budget is small -- a 6-mer appears
~0.73x per 1.5 kb, so the pair costs about one site per 1.5 kb of design space.
The 1-mismatch variants that (a) and (c) turn on are NOT forbidden outright, so
the geometry halves of this rule still do real work.

THE QUANTITATIVE PATH IS NOT BUILT, AND SAYS SO
------------------------------------------------
`brief.md:111` names the Salis Promoter Calculator -- 346 parameters, validated
on 17,396 promoters -- as the quantitative scorer. It is a MODEL, not a
threshold, and reimplementing it here would be a second, unvalidated model
wearing a validated one's citation. `_scored_path_unavailable` reports its
absence with magnitude 0.0, exactly as `d3_splicing` does for MaxEntScan and for
the same reason c1_cai returns NaN rather than 0.0: a rule that silently reports
only what it could compute reads as a clean scan of everything it claims to
cover.

REVERSE STRAND (CLAUDE.md 3.4)
-------------------------------
The motif parts list FORWARD patterns only and let the solver close the set. The
geometry scan is a directional model -- a promoter reads one way -- so it reads
`strand_for(ctx, slot)` per slot rather than scanning both strands itself. The
antisense hazard is not handled by scanning the other strand; it is handled by
(e), which is a forward motif on the sense strand. That is the whole point of the
brief giving a specific hexamer rather than asking for a reverse-strand promoter
architecture scan.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from itertools import combinations, product
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
    strand_for,
)
from bt5.core.types import Construct, Interval, Strand, reverse_complement
from bt5.verify import find_motifs

#: `brief.md:111`. The two classical sigma-70 elements, each matched within one
#: mismatch. Held as consensus strings and expanded once, because 19 C-level
#: `str.find` scans beat a Python loop comparing six characters per position --
#: this rule is evaluated once per repair candidate.
MINUS_35 = "TTGACA"
MINUS_10 = "TATAAT"

#: `brief.md:111`, inclusive on both ends.
SPACER_MIN = 15
SPACER_MAX = 19

#: `brief.md:111` (b) and (d). The brief writes the extended -10 as `TGnTATAAT`,
#: and that IUPAC form is kept here because it is what the row says.
EXTENDED_MINUS_10 = "TGNTATAAT"
SIGMA38_MINUS_10 = "TATACT"

#: The extended -10 PRE-EXPANDED to its four ACGT forms, which is what
#: `lattice_terms` actually declares.
#:
#: The solver itself handles IUPAC: PR #77 closed issue #73 by expanding
#: degenerate bases before the reverse-complement closure
#: (docs/decisions/2026-09-01-expand-forbidden-iupac.md), and one `N` is four
#: patterns, far under `MAX_PATTERN_EXPANSION = 1024`. But the DESIGN lane still
#: carries the interim guard #73 itself calls "a guard, not a fix":
#: `design/catalog.py:77` raises `DesignError` on any non-ACGT character in the
#: forbidden set, and it was never retired when #77 landed. It runs before the
#: solver sees the pattern, so a bare `TGNTATAAT` fails every end-to-end design.
#:
#: Expanding here rather than removing the guard there is a lane call, not a
#: preference: `design/` is M11 and this rule owns only `rules/catalog/`. The
#: expansion is exactly what the solver would compute, so the constraint is
#: unchanged -- each of the four is still closed under reverse complement by the
#: solver. When the stale guard goes, this collapses back to the one IUPAC string.
EXTENDED_MINUS_10_PATTERNS: tuple[str, ...] = tuple(f"TG{base}TATAAT" for base in "ACGT")

#: `brief.md:112` (e). The novel part.
ANTISENSE_MINUS_10 = "ATTATA"

#: `brief.md:111` (c). Bases BETWEEN the AT-tract and the -10-like hexamer.
#: "Appropriately positioned" is the entire finding, so the offset is exact and
#: not a window.
#:
#: FOUR, not five, and the difference is a real defect this rule shipped with
#: until the citation was audited against its own source. `brief.md:111` says the
#: tract's "3' end sits exactly 5 bp upstream of a -10-like hexamer", which reads
#: naturally as five intervening bases. It is a POSITION offset: Warman & Grainger
#: place the tract at promoter positions -23..-17 and the -10 at -12..-7, and -17
#: to -12 is five positions but FOUR bases in between. Their own Table 1 primers
#: settle it -- `tatttat` + `tgac` + `tataat`, a 4 nt gap.
#:
#: Implemented as five, the rule scored ZERO breaches on both constructs the
#: 103/103 result was measured on, and fired instead on a one-base-shifted class
#: nobody has tested. Pinned by `test_the_papers_own_constructs_are_caught`.
AT_TRACTS: tuple[str, ...] = ("TATTTAT", "AATTT")
AT_TRACT_GAP = 4

#: What each declared pattern is called in a finding, keyed by the pattern as
#: `lattice_terms` declares it -- so the four expansions of the extended -10 all
#: carry the same name.
MOTIF_LABELS: dict[str, str] = {
    ANTISENSE_MINUS_10: "antisense -10 (revcomp of TATAAT) on the sense strand",
    SIGMA38_MINUS_10: "sigma-38 -10 variant",
    **dict.fromkeys(EXTENDED_MINUS_10_PATTERNS, "extended -10 (needs no -35)"),
}

#: The widest promoter architecture this rule recognises: -35 (6) + longest
#: spacer (19) + -10 (6). Read by `solver/catalog.py:286` into
#: `RulePolicy.window`, which is what gives Tier-B a repair window able to reach
#: either element. Must stay an `int`; catalog.py falls back to 50 otherwise.
ARCHITECTURE_SPAN = len(MINUS_35) + SPACER_MAX + len(MINUS_10)

#: Any breach at or above this magnitude means the rule failed. The
#: "scored path unavailable" note is 0.0 and must not trip it.
MAG_PROMOTER = 1.0


def _within_one_mismatch(consensus: str, budget: int) -> tuple[str, ...]:
    """Every ACGT string within `budget` substitutions of `consensus`.

    `TTGACA` at budget 1 is 19 strings (the exact match plus 6 positions x 3
    substitutions). Enumerated once at construction so the hot path is
    `str.find`, never a per-position character comparison.
    """
    out: set[str] = {consensus}
    for k in range(1, budget + 1):
        for where in combinations(range(len(consensus)), k):
            alternatives = [[b for b in "ACGT" if b != consensus[i]] for i in where]
            for choice in product(*alternatives):
                variant = list(consensus)
                for i, base in zip(where, choice, strict=True):
                    variant[i] = base
                out.add("".join(variant))
    return tuple(sorted(out))


def _substring(c: Construct, start: int, length: int) -> str:
    """The bases at [start, start+length), wrapping the origin when it must.

    `find_motifs` reports the pattern it was GIVEN, which for an IUPAC pattern or
    a reverse-complement hit is not the sequence actually present. `Breach.message`
    is contracted to name the exact offending substring, so read it back.
    """
    seq = c.sequence
    n = len(seq)
    if n == 0:
        return ""
    if start + length <= n:
        return seq[start : start + length]
    return seq[start:] + seq[: (start + length) % n]


def _starts(haystack: str, needles: Sequence[str]) -> set[int]:
    """Every start index in `haystack` at which any needle occurs."""
    found: set[int] = set()
    for needle in needles:
        at = haystack.find(needle)
        while at != -1:
            found.add(at)
            at = haystack.find(needle, at + 1)
    return found


@register
class CrypticTranscription:
    id: ClassVar[str] = "d5_cryptic_transcription"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "Bacterial cryptic transcription"
    enforcement: ClassVar[Enforcement] = Enforcement.HARD_LATTICE
    #: CONTESTED, not EVIDENCE_BACKED. `brief.md:110` marks D5 "H/S" and
    #: attaches no A/B/C evidence letter (the legend is at `brief.md:46`), so the
    #: badge is this rule's own call. The parts are not equally supported: the
    #: AT-tract is single-lab experimental, the dengue anchor is one construct
    #: whose architecture was partly disproved by its own authors, and `ATTATA`
    #: -- the part enforced most aggressively, unreachable by construction --
    #: rests on a computational study with no transcription assay. Badging the
    #: whole rule to its strongest part would be the overclaim.
    evidence: ClassVar[Evidence] = Evidence.CONTESTED
    direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
    unit: ClassVar[str] = "promoters"
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "AT-tract rule: 103/103 randomized AT-tracts at promoter positions "
            "-23..-17 were transcriptionally active (44 with, 59 without a partial "
            "-35). Also the source for the sigma-38 -10 consensus TATACT, and for "
            "the TGn motif -- which this paper places upstream of TATACT, not of "
            "TATAAT",
            "https://academic.oup.com/nar/article/48/9/4891/5820884",
        ),
        Citation(
            "Dengue-2 cryptic promoter: a BPROM-PREDICTED -35 TCAACG at nt 53 and "
            "-10 TTTTTAAT at nt 72 (a 13 bp spacer, not the 17 brief.md:111 states) "
            "gave ~10^6.6 mRNA copies/ug uninduced. Deletion showed the -10 "
            "essential and the -35 NOT, so this QUALIFIES the -35/-10 architecture "
            "rather than supporting it; the paper reports cloning and propagation "
            "difficulty, and does not use the word uncloneable",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC3069047/",
            sign="qualifies",
        ),
        Citation(
            "Antisense hexamer ATTATA is 2-3x under-represented (4.76% of 484,741 "
            "natural E. coli CDSs) yet 77.28% of clean CDSs can acquire it by "
            "synonymous substitution. COMPUTATIONAL only -- the study demonstrates "
            "the asymmetry, not that the hexamer transcribes",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC13029128/",
        ),
    )
    last_verified: ClassVar[str] = "2026-09-06"
    weight_provenance: ClassVar[str] = ""  # hard rule; never enters the weighted sum
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.0
    #: Forced to 0.0 by the HARD_LATTICE contract assertion. The consequence for
    #: the geometry half is argued in the module docstring.
    steering_weight: ClassVar[float] = 0.0
    band: ClassVar[tuple[float, float] | None] = None
    #: The architecture spans up to 31 bp, so a motif-width window cannot reach
    #: the other element of the pair.
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.WINDOW_MINUS_1
    #: MANDATORY here, not inherited (CLAUDE.md 3.6). Mutating one element of a
    #: promoter creates new -10-like hexamers nearby, and a new hexamer 5 bp
    #: downstream of an AT-tract that was previously harmless is a NEW promoter
    #: the single pass would never re-scan for. This rule can create instances of
    #: what it removes, which is exactly the condition FIXED_POINT exists for.
    repair: ClassVar[RepairPolicy] = RepairPolicy.FIXED_POINT
    cost_class: ClassVar[str] = "cheap"
    #: Both oppose this rule by pushing sequence the way it must not go. D5
    #: removes AT-rich hexamers, which raises local GC against f5's 60% ceiling;
    #: d8 depletes CpG, whose synonymous fixes lower GC and walk back toward the
    #: AT-rich hexamers D5 forbids. `e2_gc_band` is NOT declared: its band is
    #: (0.28, 0.80) over the whole construct, far too wide for a 6-mer
    #: substitution to bind against.
    conflicts_with: ClassVar[tuple[str, ...]] = ("d8_cpg_depletion", "f5_at_window")
    brief_ref: ClassVar[str] = "2.D5"
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "mismatch_budget": {
                "type": "integer",
                "minimum": 0,
                "maximum": 2,
                "default": 1,
                "description": (
                    "Substitutions tolerated when matching the -35 and -10 consensus "
                    "hexamers. brief.md:111 specifies 1; 2 is offered for a paranoid "
                    "scan and costs 154 variants per consensus instead of 19."
                ),
            },
            "spacer_bp": {
                "type": "array",
                "items": {"type": "integer"},
                "default": [SPACER_MIN, SPACER_MAX],
                "description": ("Inclusive -35/-10 spacer bounds. brief.md:111 gives 15-19."),
            },
        },
    }

    def __init__(
        self,
        mismatch_budget: int = 1,
        spacer_bp: tuple[int, int] | None = None,
    ) -> None:
        self.mismatch_budget = mismatch_budget
        self.spacer_min, self.spacer_max = spacer_bp or (SPACER_MIN, SPACER_MAX)
        self.minus_35 = _within_one_mismatch(MINUS_35, mismatch_budget)
        self.minus_10 = _within_one_mismatch(MINUS_10, mismatch_budget)
        #: Read by `solver/catalog.py:286` into `RulePolicy.window`.
        self.window: int = len(MINUS_35) + self.spacer_max + len(MINUS_10)
        #: No promoter-strength model ships. Injectable so the scored path is
        #: real code rather than a comment, and so the honest "unavailable" is
        #: what a caller without a model actually gets.
        self.promoter_calculator: object | None = None

    def gate(self, slot: ContextSlot) -> bool:
        """Every slot, including IVT mRNA -- and that exclusion is the trap.

        `d3`, `d4` and `d6` all gate OFF `IVT_MRNA`, and copying them here would
        be wrong. Those three read the TRANSCRIPT: an in-vitro transcript is not
        replicated and has no Pol II, so their hazards cannot fire. D5 reads the
        TEMPLATE PLASMID during *E. coli* propagation, which an IVT mRNA
        construct still goes through. `brief.md:110` says so without
        qualification -- "applies to ALL constructs regardless of final host".

        That one line carries it alone. `brief.md:209` was also cited here and
        should not have been: it sits in section 3.1, "by expression host", whose
        columns are E. coli / yeast / CHO-HEK / insect / plant and which has no
        IVT mRNA column at all. The table that does (3.2, `brief.md:219-231`)
        carries no cryptic-promoter row -- which is the gap that makes d3, d4 and
        d6 gate IVT off, and is not evidence about this rule either way.
        """
        return True

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        """HARD_REPAIR, which is what routes the geometry to Tier-B.

        Deliberately NOT the HARD_LATTICE ClassVar. Returning that would leave
        (a) and (c) matching neither arm of `findings()`'s router
        (`solver/catalog.py:184-189`), so every geometry breach would be
        discovered, reported into the scorecard, and then silently dropped by the
        solver -- the rule's own findings ignored by the thing that could fix
        them. The motif half does not need routing: it is already unreachable.
        """
        return Enforcement.HARD_REPAIR

    def lattice_terms(self, ctx: DesignContext | None) -> LatticeTerms:
        """The three parts a finite pattern set can express. Forward only.

        The solver closes these under reverse complement at automaton
        construction, so `ATTATA` also removes sense-strand `TATAAT` -- argued in
        the module docstring, and kept.
        """
        return LatticeTerms(
            forbidden=(ANTISENSE_MINUS_10, *EXTENDED_MINUS_10_PATTERNS, SIGMA38_MINUS_10)
        )

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        """Motif parts on the closed set; geometry once per slot, on its strand."""
        breaches: list[Breach] = list(self._motif_breaches(c))

        seen: set[tuple[int, int, int, str]] = set()
        for slot in ctx.active_slots:
            if not self.gate(slot):
                continue
            strand = strand_for(ctx, slot)
            for breach in self._geometry_breaches(c, strand, slot.role):
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

        if self.promoter_calculator is None:
            breaches.append(self._scored_path_unavailable(c))

        return Evaluation(
            spec_id=self.id,
            passes=not any(b.magnitude >= MAG_PROMOTER for b in breaches),
            raw_score=float(sum(b.magnitude for b in breaches)),
            breaches=tuple(breaches),
            n_evaluated=c.length,
        )

    def _motif_breaches(self, c: Construct) -> list[Breach]:
        """(b), (d) and (e): the parts the automaton guarantees.

        Still evaluated, for two reasons. A backbone-carried instance is real,
        unfixable and belongs in the report; and the scorecard reads
        `Evaluation`, so a rule that reported nothing here would look like a
        clean scan of the parts it never checked.
        """
        labels = MOTIF_LABELS
        out: list[Breach] = []
        for motif, pos in find_motifs(c, list(labels)):
            iv = Interval(pos, pos + len(motif))
            found = _substring(c, pos, len(motif))
            out.append(
                Breach(
                    spec_id=self.id,
                    interval=iv,
                    magnitude=MAG_PROMOTER,
                    message=f"{labels[motif]}: {found!r} at {pos} (pattern {motif})",
                    fixable_by_codon_choice=c.overlaps_editable(iv),
                    detail={"kind": "motif", "motif": motif, "found": found},
                )
            )
        return out

    def _geometry_breaches(self, c: Construct, strand: Strand, role: str) -> list[Breach]:
        """(a) and (c): two elements at a declared distance, on ONE strand.

        Positions are reported in construct coordinates whichever strand found
        them -- a breach the user cannot find on their own map is not actionable.
        The `% n` mapping is the one from `d4_internal_polya:197`: a minus-strand
        hit crossing the origin goes negative, and clamping it to zero would move
        the breach onto bases that are not the motif it names.
        """
        n = c.length
        if n == 0:
            return []
        forward = c.sequence
        seq = forward if strand == 1 else reverse_complement(forward)
        span = len(MINUS_35) + self.spacer_max + len(MINUS_10)
        scan = seq + seq[: span - 1] if c.is_circular else seq

        m35 = _starts(scan, self.minus_35)
        m10 = _starts(scan, self.minus_10)
        out: list[Breach] = []

        def fwd(start: int, length: int) -> int:
            """A scan-strand start, in CONSTRUCT coordinates.

            Every position that reaches a user goes through here, the interval
            and the message alike. They used to disagree: `emit` remapped the
            interval while the messages were formatted from raw scan indices, so
            on a reverse-oriented cassette the -35 and -10 were reported at each
            other's coordinates -- a breach whose prose contradicted its own
            interval.
            """
            return (start if strand == 1 else n - start - length) % n

        def emit(start: int, length: int, kind: str, message: str) -> None:
            if start >= n:
                return  # the wrapped copy of a hit already reported at its true start
            lo = fwd(start, length)
            iv = Interval(lo, lo + length, strand)
            out.append(
                Breach(
                    spec_id=self.id,
                    interval=iv,
                    magnitude=MAG_PROMOTER,
                    message=message,
                    slot_role=role,
                    fixable_by_codon_choice=c.overlaps_editable(iv),
                    detail={"kind": kind, "strand": float(strand)},
                )
            )

        # (a) -35 + 15-19 bp spacer + -10.
        for i in sorted(m35):
            for spacer in range(self.spacer_min, self.spacer_max + 1):
                j = i + len(MINUS_35) + spacer
                if j in m10:
                    emit(
                        i,
                        len(MINUS_35) + spacer + len(MINUS_10),
                        "sigma70_architecture",
                        f"sigma-70 promoter architecture on the {'+' if strand == 1 else '-'} "
                        f"strand: a -35-like hexamer at {fwd(i, len(MINUS_35))}, a {spacer} bp "
                        f"spacer, then a -10-like hexamer at {fwd(j, len(MINUS_10))}",
                    )

        # (c) AT-tract whose 3' end sits exactly 5 bp upstream of a -10-like hexamer.
        for j in sorted(m10):
            for tract in AT_TRACTS:
                end = j - AT_TRACT_GAP
                start = end - len(tract)
                if start < 0 or scan[start:end] != tract:
                    continue
                emit(
                    start,
                    len(tract) + AT_TRACT_GAP + len(MINUS_10),
                    "at_tract",
                    f"AT-tract {tract!r} at {fwd(start, len(tract))} with exactly "
                    f"{AT_TRACT_GAP} bases between it and a -10-like hexamer at "
                    f"{fwd(j, len(MINUS_10))}; 103/103 randomized tracts in this "
                    f"position were active",
                )

        return out

    def _scored_path_unavailable(self, c: Construct) -> Breach:
        """The Salis Promoter Calculator did not run, said out loud.

        Magnitude 0.0 and not fixable, so it lands on `RepairOutcome.advisory`
        and never becomes a target the solver chases -- `d3_splicing`'s pattern
        for its absent MaxEntScan model, and the same argument.
        """
        where = sorted(c.editable)[0] if c.editable else Interval(0, min(1, c.length) or 1)
        reason = (
            "the Salis Promoter Calculator (346 parameters, validated on 17,396 "
            "promoters) is not implemented, so no promoter STRENGTH was computed; "
            "the motif and spacer-geometry scans above did run, and they report "
            "presence rather than strength"
        )
        return Breach(
            spec_id=self.id,
            interval=where,
            magnitude=0.0,
            message=f"cryptic promoter strength unavailable: {reason}",
            fixable_by_codon_choice=False,
            detail={"kind": "unavailable", "unavailable_reason": reason},
        )
