"""B9 -- additional and out-of-frame ATGs near the start, and across the junction.

brief.md:69, graded **A**, type **H**: "No additional ATG in the first 50 nt of
CDS; penalize any out-of-frame ATG anywhere that has BOTH -3 purine and +4 G.
Ensure optimization does not create an upstream out-of-frame AUG at the UTR/CDS
junction." Three clauses, and this file implements all three against one scan.

The effect size is on B10, B9's non-editable twin (brief.md:70): uORFs occur in
~half of human transcripts and typically cut protein 30-80%. B10 covers the
uAUGs already in the user's 5'UTR, which no codon choice can reach; B9 covers the
ones the optimizer would otherwise CREATE, which is the half BT5 is responsible
for.

**Why HARD_REPAIR and not HARD_LATTICE.** `forbidden` motifs are the right answer
for a motif rule (CLAUDE.md 3.4), and they are the wrong answer here: `ATG` as a
forbidden string would ban the start codon and every in-frame Met. What makes an
ATG a defect is its FRAME and its POSITION, and the Aho-Corasick automaton knows
neither -- it decides from a bounded codon suffix with no notion of offset from
the CDS start. So the constraint is enforced by repair plus the independent
validator, which refuses to emit.

**`RepairPolicy.FIXED_POINT`, and it is mandatory.** This is the splice-donor case
CLAUDE.md 3.6 names. Recoding one codon to destroy an ATG shifts the bases around
it and can manufacture a new ATG at a neighbouring offset -- `...ATGA...` repaired
at one codon becomes `...GATG...` just as easily. A single pass would ship a
construct whose ATGs were removed INTO new ATGs and the validator would pass it,
because the specific offset it was told about is clean.

**Directional, and deliberately scanned by hand.** CLAUDE.md 3.4 bans hand-rolling
a reverse-strand scan for MOTIF rules, because the solver closes `forbidden` under
reverse complement. B9 is the exception the same rule names: an out-of-frame ATG
matters only on the SENSE strand of the transcript being made. The reverse
complement of the transgene is a different molecule that no ribosome scans here,
and closing this set under reverse complement would flag every CAT triplet in the
construct. So the strand comes from `strand_for(ctx, slot)` and the scan runs in
reading order, never from a hard-coded strand 1.

**The first-50-nt clause is literal, including in-frame Met.** brief.md:69 says "No
ADDITIONAL ATG", not "no out-of-frame ATG", for that window: an in-frame ATG at
codon 4 is a perfectly good alternative start that produces an N-terminally
truncated protein. Where the residue is genuinely Met the codon is forced, so the
breach ships `fixable_by_codon_choice=False` and goes to the advisor rather than
driving the search into an infeasibility it cannot resolve (docs/PLAN.md:372).
Whether Met is forced is read from the INJECTED table, never hard-coded: NCBI
tables differ on which codons a residue has (CLAUDE.md 3.1).

**Scope boundary: an upstream IN-FRAME AUG is not reported.** brief.md:69's third
clause asks only about an upstream OUT-OF-FRAME AUG. An in-frame one makes an
N-terminally EXTENDED product rather than a truncated or frameshifted one -- a
real design concern, and a different one the brief does not put on this row. It
is skipped deliberately rather than folded in under a message that would describe
it wrongly.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import ClassVar

from bt5.core.context import ContextSlot, DesignContext
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
from bt5.rules.catalog.b1_five_prime import FIVE_PRIME_UTR_KINDS
from bt5.rules.catalog.b8_kozak import EUKARYOTIC_HOSTS, PURINES

#: brief.md:69. Offsets 0..49 of the CDS; an ATG counts as inside when it STARTS
#: inside, so one at offset 48 is caught rather than half-ignored.
SCAN_NT = 50

#: How much annotated 5'UTR to pull in so the junction clause can run. Six bases
#: is the smallest window that lets an ATG starting two bases before the junction
#: still have its own -3 evaluated, and it matches B8's context width.
JUNCTION_UPSTREAM = 6

#: The magnitude floor for the hard clause, kept above every context-only finding
#: so a first-50 hit always outranks a distant strong-context one in a sort.
HARD_MAGNITUDE = 10.0


@register
class OutOfFrameAtg:
    id: ClassVar[str] = "b9_out_of_frame_atg"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "Additional and strong-context out-of-frame ATGs"
    #: Repair plus the independent validator. See the module docstring for why
    #: HARD_LATTICE cannot express a frame- and position-dependent constraint.
    enforcement: ClassVar[Enforcement] = Enforcement.HARD_REPAIR
    #: brief.md:69 grades B9 "A", and B10 carries the measured number.
    evidence: ClassVar[Evidence] = Evidence.EVIDENCE_BACKED
    direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
    unit: ClassVar[str] = "count of additional / strong-context out-of-frame ATGs"
    band: ClassVar[tuple[float, float] | None] = None
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "uORFs occur in ~half of human transcripts and typically cut protein "
            "30-80%. The effect size behind B9, carried on its non-editable twin "
            "B10 (brief.md:70): B10 reports the uAUGs already in the user's 5'UTR, "
            "which no codon choice reaches, and B9 stops the optimizer creating "
            "new ones",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC7100133/",
            2020,
            sign="supports",
        ),
        Citation(
            "Noderer 2014's -3 purine / +4 G conjunction is the same initiation "
            "context B8 scores, applied here to the WRONG start: an out-of-frame "
            "ATG in strong context is the one a scanning ribosome is most likely "
            "to initiate at, which is why the brief penalises those specifically "
            "rather than every out-of-frame ATG in the sequence",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC4299517/",
            2014,
            sign="supports",
        ),
    )
    last_verified: ClassVar[str] = "2026-09-02"
    #: Hard rule: default_weight is 0.0 and this explains the enforcement class,
    #: not a weight. The weighted sum only ever sees SOFT (CLAUDE.md 3.5).
    weight_provenance: ClassVar[str] = (
        "Hard, so default_weight is 0.0 and no preset may weight 2.B9. "
        "steering_weight is 0.0 as well, and that is a limitation rather than a "
        "judgement: the Tier-A DP would benefit from a nudge away from codons that "
        "complete an ATG across a codon boundary, and there is no channel for it. "
        "solver/catalog.py:151-153 reads LatticeTerms.forbidden and nothing else, "
        "so codon_weights, codon_pair_weights and positional are all inert; a "
        "non-zero steering_weight would claim a nudge the engine does not perform. "
        "Enforcement is by repair plus the independent validator instead."
    )
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.0
    steering_weight: ClassVar[float] = 0.0
    #: Each finding is one ATG triplet, and the message names its offset and frame.
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.MOTIF_LEN_MINUS_1
    #: MANDATORY, not a preference: removing one ATG can create another at a
    #: neighbouring offset. CLAUDE.md 3.6, and the module docstring.
    repair: ClassVar[RepairPolicy] = RepairPolicy.FIXED_POINT
    cost_class: ClassVar[str] = "cheap"
    #: Raising +4 to G for B8 strengthens the context of anything downstream that
    #: happens to be an out-of-frame ATG, and removing an ATG can spoil a Kozak.
    conflicts_with: ClassVar[tuple[str, ...]] = ("b8_kozak",)
    brief_ref: ClassVar[str] = "2.B9"
    #: A triplet scan, not a folding energy.
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "scan_nt": {
                "type": "integer",
                "default": SCAN_NT,
                "minimum": 3,
                "description": (
                    "Length of the CDS 5' window in which ANY additional ATG is a "
                    "breach, in frame or not. brief.md:69's 50 nt."
                ),
            },
            "strong_context_only": {
                "type": "boolean",
                "default": True,
                "description": (
                    "Beyond the 5' window, report only out-of-frame ATGs carrying "
                    "BOTH -3 purine and +4 G, which is what brief.md:69 asks for. "
                    "Turn off to report every out-of-frame ATG in the transcript; "
                    "expect many findings on any real CDS."
                ),
            },
        },
    }

    def __init__(self, scan_nt: int = SCAN_NT, strong_context_only: bool = True) -> None:
        if scan_nt < 3:
            raise ValueError(f"scan_nt must be at least one codon, got {scan_nt}")
        self.scan_nt = scan_nt
        self.strong_context_only = strong_context_only

    def gate(self, slot: ContextSlot) -> bool:
        """Eukaryotic translating slots, for B8's reason.

        Scanning initiation is the mechanism the rule is about; a bacterial
        ribosome finds its start by Shine-Dalgarno pairing, and internal
        initiation there is B7's row with its own TIR model (brief.md:67).
        """
        return slot.host in EUKARYOTIC_HOSTS and slot.role != "propagation"

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """None, and `forbidden` would be actively wrong here.

        Listing `ATG` would ban the start codon and every in-frame Met, and the
        solver closes `forbidden` under reverse complement (CLAUDE.md 3.4), so it
        would additionally ban every `CAT`. What makes an ATG a defect is frame
        and offset, which the automaton cannot see.
        """
        return None

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        slots = [s for s in ctx.active_slots if self.gate(s)]
        if not slots:
            return self._unavailable(
                c,
                "no eukaryotic translating slot in this context. Scanning initiation "
                "is the mechanism this rule is about; bacterial internal initiation "
                "is B7's row and uses a TIR model, not an ATG scan",
            )

        editable = sorted(c.editable)
        if not editable:
            return self._unavailable(c, "no designable CDS to scan for additional ATGs")

        breaches: list[Breach] = []
        seen: set[tuple[int, int, int]] = set()
        scanned = 0
        junction_checked = False  # rebound per slot inside the loop
        for slot in slots:
            strand = strand_for(ctx, slot)
            cds = editable[0] if strand == 1 else editable[-1]
            try:
                code = svc.tables.genetic_code(slot.table_id)
            except (OSError, LookupError, NotImplementedError, ValueError) as exc:
                return self._unavailable(
                    c, f"NCBI table {slot.table_id} could not be loaded: {exc}"
                )
            lead = self._leader(c, cds, strand)
            # Per slot, NOT a running OR across slots: a plus-strand slot with an
            # annotated leader must not make a minus-strand slot (whose leader is
            # the unannotated trailer) report that its junction was scanned.
            junction_checked = lead > 0
            span = self._span(c, cds, strand, lead)
            if span is None:
                return self._unavailable(
                    c, "the CDS and its junction context could not be read from the construct"
                )
            transcript = c.slice(span)
            scanned = max(scanned, len(transcript))
            for breach in self._scan(c, transcript, lead, span, slot, code, junction_checked):
                # One ATG is one finding, however many slots read it. Two
                # eukaryotic translating slots (producer HEK293 + target CHO) are
                # a routine lentiviral job and resolve to the same strand and the
                # same sequence, so unioning would put the same interval in the
                # conflict panel twice and double a raw_score whose declared unit
                # is a COUNT. B8 takes `min` over slots and C3 picks one binding
                # slot; this is the same choice, made by interval.
                key = (breach.interval.start, breach.interval.end, breach.interval.strand)
                if key not in seen:
                    seen.add(key)
                    breaches.append(breach)

        return Evaluation(
            spec_id=self.id,
            passes=not breaches,
            raw_score=float(len(breaches)),
            breaches=tuple(breaches),
            n_evaluated=scanned,
        )

    # -- internals ----------------------------------------------------------

    def _scan(
        self,
        c: Construct,
        transcript: str,
        lead: int,
        span: Interval,
        slot: ContextSlot,
        code: GeneticCode,
        junction_checked: bool,
    ) -> list[Breach]:
        """Every ATG in reading order, classified by frame and offset.

        `lead` is how many bases of annotated 5'UTR precede the CDS in
        `transcript`, so a CDS-relative offset is `index - lead` and is negative
        for an ATG starting upstream of the CDS. Nothing can literally straddle
        the junction while the CDS begins with ATG -- a triplet ending at the
        junction would need the CDS's own A as its G -- so a negative offset means
        "upstream", which is the uAUG the third clause of brief.md:69 is about.
        """
        out: list[Breach] = []
        for index in range(len(transcript) - 2):
            if transcript[index : index + 3] != "ATG":
                continue
            offset = index - lead
            if offset == 0:
                continue  # the start codon itself
            # `offset % 3` is already the right frame class for negatives in
            # Python (-1 % 3 == 2), so the old `offset >= 0` guard did not dodge
            # a hazard -- it mislabelled every upstream in-frame AUG as
            # out-of-frame while `detail["frame"]` said 0.
            in_frame = offset % 3 == 0
            near_start = 0 <= offset < self.scan_nt
            strong = self._strong_context(transcript, index)
            if near_start:
                reason = f"additional ATG in the first {self.scan_nt} nt of the CDS"
            elif not in_frame and strong and self.strong_context_only:
                reason = "out-of-frame ATG in strong initiation context"
            elif in_frame and offset < 0:
                # An upstream IN-frame AUG makes an N-terminally EXTENDED product,
                # not a truncated or frameshifted one, so brief.md:69's third
                # clause (out-of-frame uAUGs) does not cover it. Nothing else in
                # the catalog does either, and dropping it silently was the wrong
                # half of the choice between a wrong message and no message.
                reason = "upstream in-frame ATG, an N-terminal extension"
            elif not in_frame and not self.strong_context_only:
                reason = "out-of-frame ATG"
            else:
                continue
            out.append(
                self._breach(
                    c,
                    index,
                    offset,
                    span,
                    slot,
                    code,
                    in_frame,
                    near_start,
                    strong,
                    reason,
                    junction_checked,
                    transcript,
                )
            )
        return out

    def _strong_context(self, transcript: str, index: int) -> bool:
        """B8's conjunction applied to this ATG: -3 purine AND +4 G.

        An ATG without three readable bases in front of it cannot be shown to be
        in strong context, so it is not claimed to be -- the first-50 clause is
        what catches those, and it does not depend on context.
        """
        if index < 3 or index + 3 >= len(transcript):
            return False
        return transcript[index - 3] in PURINES and transcript[index + 3] == "G"

    def _breach(
        self,
        c: Construct,
        index: int,
        offset: int,
        span: Interval,
        slot: ContextSlot,
        code: GeneticCode,
        in_frame: bool,
        near_start: bool,
        strong: bool,
        reason: str,
        junction_checked: bool,
        transcript: str,
    ) -> Breach:
        where = self._locate(span, index, c.length)
        forced = in_frame and self._met_is_forced(code)
        upstream_of_cds = offset < 0
        # Built outside the f-string: a multi-line expression inside one is
        # PEP 701, i.e. 3.12+, and this project is 3.11 only.
        note = (
            ""
            if junction_checked
            else " (no annotated 5'UTR, so the UTR/CDS junction was not scanned)"
        )
        if forced:
            note = (
                "; the residue is Met, so this ATG is forced and must be moved or "
                "accepted rather than recoded"
            ) + note
        outcome = (
            "an N-terminally extended product"
            if in_frame and upstream_of_cds
            else "a truncated or frameshifted product"
        )
        frame = "in frame" if in_frame else f"frame +{offset % 3}"
        if upstream_of_cds:
            frame += ", starts upstream of the CDS"
        if strong:
            frame += ", -3 purine and +4 G"
        return Breach(
            spec_id=self.id,
            interval=where,
            magnitude=HARD_MAGNITUDE if near_start else 1.0,
            message=(
                f"{reason} at CDS offset {offset:+d} ({frame}) for the "
                f"{slot.role} slot ({slot.host}). A scanning ribosome can initiate "
                f"here and make {outcome}{note}"
            ),
            # Forced Met cannot be recoded, and an ATG outside the designable
            # region is not the solver's to touch.
            fixable_by_codon_choice=not forced and c.overlaps_editable(where),
            slot_role=slot.role,
            detail={
                "cds_offset": float(offset),
                "frame": float(offset % 3),
                "in_frame": str(in_frame),
                "near_start": str(near_start),
                "strong_context": str(strong),
                "upstream_of_cds": str(upstream_of_cds),
                "forced_met": str(forced),
                "junction_checked": str(junction_checked),
                "context": transcript[max(0, index - 6) : index + 5],
                "host": str(slot.host),
            },
        )

    @staticmethod
    def _met_is_forced(code: GeneticCode) -> bool:
        """Does the injected table give Met any codon other than ATG?

        Read from the table rather than assumed, per CLAUDE.md 3.1: which codons a
        residue has is table-dependent, and a rule that hard-codes "Met is ATG
        only" is wrong on the tables where it is not.
        """
        try:
            return len(code.synonymous_codons("M")) < 2
        except (ValueError, KeyError):
            return True

    def _leader(self, c: Construct, cds: Interval, strand: Strand) -> int:
        """Bases of ANNOTATED 5'UTR to pull in ahead of the CDS, 0 if none.

        Unannotated upstream sequence is not assumed to be leader: it may be
        promoter or backbone, and an ATG there is not one a ribosome scanning this
        transcript would ever meet. B1 refuses outright for the same reason; B9
        degrades instead, because its first-50 clause does not need the leader and
        silently skipping the whole rule would be worse than skipping one clause.
        """
        probe = self._span(c, cds, strand, JUNCTION_UPSTREAM)
        if probe is None:
            return 0
        # On the minus strand the transcript's 5' end is at HIGHER coordinates, so
        # the leader is the TAIL of the span -- b1's `leader_of`, same argument.
        if strand == 1:
            leader = Interval(probe.start, probe.start + JUNCTION_UPSTREAM, strand)
        else:
            leader = Interval(probe.end - JUNCTION_UPSTREAM, probe.end, strand)
        annotated = any(
            f.interval.contains(leader, c.length, c.is_circular)
            for f in c.features
            if f.kind in FIVE_PRIME_UTR_KINDS
        )
        return JUNCTION_UPSTREAM if annotated else 0

    def _span(self, c: Construct, cds: Interval, strand: Strand, lead: int) -> Interval | None:
        """`lead` bases of leader plus the whole CDS, in reading order.

        The leader is prepended in CONSTRUCT coordinates on the plus strand and
        appended on the minus strand, because that is where the transcript's 5'
        end is in each case. Returns None rather than a clamped span when the
        leader would run off a linear end.
        """
        length = cds.length + lead
        start = cds.start - lead if strand == 1 else cds.start
        if start < 0:
            if not c.is_circular:
                return None
            start += c.length
        if not c.is_circular and start + length > c.length:
            return None
        return Interval(start, start + length, strand)

    def _locate(self, span: Interval, index: int, construct_length: int) -> Interval:
        """The ATG's own three bases, in construct coordinates.

        On the minus strand reading order runs from high coordinates to low, so
        the triplet's construct start is `span.end - index - 3`.
        """
        if span.strand == 1:
            start = (span.start + index) % construct_length
        else:
            start = (span.end - index - 3) % construct_length
        # Carrying `span.strand`: a minus-strand hit reported on a plus-strand
        # interval points at bases that read CAT, and `Interval.contains` tells
        # callers to compare `.strand` themselves.
        return Interval(start, start + 3, span.strand)

    def _unavailable(self, c: Construct, reason: str) -> Evaluation:
        """NaN plus a breach carrying the reason -- B1's pattern, B1's argument.

        NaN rather than 0.0 because 0.0 is a real and good score here: it is "no
        additional ATGs found". Reporting it for a scan that never ran would tell
        a user the one thing this rule exists to deny them.
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
                    message=f"out-of-frame ATG scan unavailable: {reason}",
                    fixable_by_codon_choice=False,
                    detail={"unavailable_reason": reason},
                ),
            ),
            n_evaluated=0,
        )
