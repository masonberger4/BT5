"""B8 -- Kozak initiation context, scored on the two positions the brief defines.

brief.md:68 states the target as `GCCRCCATGG` (equivalently `gccRccATGG`), asks
that +4 = G and +5 = C be set by codon-2 choice where the residue permits, and
grades the row **A** -- the brief's own legend (brief.md:46) reads "A = large
controlled dataset or physical necessity". The dataset is Noderer 2014, which
measured all 65,536 variants of the -6..+5 context and found a 12-fold range,
-3 purine worth +58% over -3U, and +4G/+5C 24.8% better than +4G/+5A with
cooperativity between the two.

**What this rule scores, and what it does not.** The brief defines exactly one
tier: "Score strong (-3 purine AND +4 G) / adequate / weak". It names the
conjunction for *strong* and leaves the other two words undefined. Rather than
invent cutoffs, the ordinal here is the minimal completion of a three-level scale
whose top is a two-term conjunction -- strong = both determinants, adequate =
exactly one, weak = neither. That is a reading of the brief, not a measurement,
and it is flagged as such in `weight_provenance` so a later reviewer can replace
it with Noderer's own efficiency table rather than discovering a number nobody
sourced. +5 is reported but NOT scored: its 24.8% is measured relative to +4G, so
it is a refinement of a context that already has +4 G, and folding it into the
same ordinal would double-count the +4 term.

**The -3 base is usually not yours to change.** +4 and +5 are codon 2 and are
designable; -3 is in the user's 5'UTR. So a weak context caused only by -3 sets
`fixable_by_codon_choice=False` and goes to the advisor, while one caused by +4
is a real instruction to the solver. Reporting both as fixable would send the
search after a base it cannot edit (docs/PLAN.md:372).

**Gating on host is correct here, and it is the one place it is.**
`docs/decisions/2026-09-01-c1-cai-soft-band.md` warns that host-keyed rules are
dangerous because a lentiviral job propagates in E. coli and expresses in HEK293.
That warning is about inferring the TRANSLATING host from the presence of an
E. coli slot. Kozak is different in kind: brief.md:205 keys the rule family on
host explicitly -- "Kozak (B8/B9) | No (SD instead) | Yes | Yes, +4/+5 designable"
-- because a bacterial ribosome does not read a Kozak context at all; it reads a
Shine-Dalgarno, which is B6/B7's row. Scoring Kozak on an E. coli slot would not
be a number measured on the wrong host, it would be a number measuring nothing.
The propagation slot is excluded on top, for C1's reason.

**The NcoI conflict is surfaced, not silently broken.** brief.md:96 asks for
exactly that: "Surface the NcoI CCATGG subset-of Kozak GCCACCATGG conflict rather
than silently breaking one." A strong Kozak ending ...CCATGG contains NcoI, so
raising +4 to G is what manufactures the site. D1 is HARD_LATTICE, so when NcoI is
in the selected enzyme set the automaton simply will not emit that codon and this
rule's preference loses -- which is the correct precedence, a hard constraint over
a soft objective. `conflicts_with` names it so the conflict panel can show it.

**No lattice term, and no steering weight.** The tempting term is
`LatticeTerms.positional`, which is shaped exactly for "prefer a G-initial codon at
index 1". It is declared on the dataclass and read by nothing: `solver/catalog.py`
consumes `terms.forbidden` and only that, so `positional`, `codon_weights` and
`codon_pair_weights` are inert. A non-zero `steering_weight` here would claim a
nudge the engine does not perform.
"""

from __future__ import annotations

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
    strand_for,
)
from bt5.core.types import Construct, Interval
from bt5.rules.catalog.b1_five_prime import FIVE_PRIME_UTR_KINDS, leader_of

#: Noderer 2014 measured -6..+5. brief.md:68's target `GCCRCCATGG` spans -6..+4
#: and the +5 = C instruction extends it by one, so -6..+5 is both the measured
#: window and the reported context.
UPSTREAM = 6
DOWNSTREAM = 5

#: Purines. The -3 determinant is "purine", not "A": brief.md:68 writes the target
#: with `R` at -3 precisely because A and G both qualify, and Noderer's +58% is
#: measured for purine against U.
PURINES = frozenset("AG")

#: Hosts whose ribosomes scan for a Kozak context. brief.md:205 keys the B8/B9 row
#: on host and reads "No (SD instead)" for E. coli -- the bacterial path is B6/B7,
#: not a weaker version of this one. Listed positively rather than as "not E. coli"
#: so a host added later has to be classified deliberately.
EUKARYOTIC_HOSTS = frozenset(
    {
        HostId.HUMAN,
        HostId.HEK293,
        HostId.CHO,
        HostId.S_CEREVISIAE,
        HostId.P_PASTORIS,
        HostId.SF9,
        HostId.MOUSE,
    }
)

#: What `TableProvider` implementations raise when a table is absent.
_MISSING_TABLE = (OSError, LookupError, NotImplementedError, ValueError)

#: The three tiers of brief.md:68, as the raw score. HIGHER_IS_BETTER.
WEAK, ADEQUATE, STRONG = 0.0, 1.0, 2.0


@register
class KozakContext:
    id: ClassVar[str] = "b8_kozak"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "Kozak initiation context strength at -3 and +4"
    enforcement: ClassVar[Enforcement] = Enforcement.SOFT
    #: brief.md:68 grades B8 "A" = large controlled dataset. Noderer 2014 is a
    #: complete 65,536-variant enumeration of the context, which is as close to
    #: exhaustive as this catalog's evidence gets.
    evidence: ClassVar[Evidence] = Evidence.EVIDENCE_BACKED
    direction: ClassVar[Direction] = Direction.HIGHER_IS_BETTER
    unit: ClassVar[str] = "Kozak tier (0 weak, 1 adequate, 2 strong)"
    band: ClassVar[tuple[float, float] | None] = None
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "Noderer 2014 measured all 65,536 variants of the -6..+5 initiation "
            "context and found a 12-fold range: -3 purine is worth +58% over -3U, "
            "and +4G/+5C is 24.8% better than +4G/+5A with cooperativity between "
            "+4 and +5. The source of both determinants this rule scores, and the "
            "reason +5 is reported but not scored separately -- its effect is "
            "measured relative to a context that already has +4 G",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC4299517/",
            2014,
            sign="supports",
        ),
        Citation(
            "Cambray 2018: all computable design features together explain 5-31% "
            "(mean ~14%) of protein-level variance. A strong Kozak is a real and "
            "well-measured lever and still sits under this ceiling, which is why "
            "this rule reports a tier and never a predicted expression level",
            "https://www.nature.com/articles/nbt.4238",
            2018,
            sign="qualifies",
        ),
    )
    last_verified: ClassVar[str] = "2026-09-02"
    weight_provenance: ClassVar[str] = (
        "0.8 -- the highest in the catalog AFTER b1_five_prime, and strictly below "
        "it. An earlier draft set 1.0, reasoning that B1 gates to the "
        "BACTERIAL_EXPRESSION modality and B8 to eukaryotic hosts, so on any sane "
        "context they do not both fire and the number is not a ranking of one "
        "against the other. `test_b1_five_prime.py::TestContract::"
        "test_it_is_soft_and_carries_the_highest_weight_in_the_catalog` says "
        "otherwise, and it is right: it asserts B1 outranks EVERY other spec "
        "catalog-wide, because B1 is 'the only objective in BT5 justified by a "
        "measured effect size (r = 0.66, 44% of variance) rather than a feature "
        "ranking'. Noderer's numbers are within-context percentages on individual "
        "determinants, not variance explained on the objective, so they do not "
        "meet that bar however large the dataset. 0.8 keeps B8 second only to B1, "
        "which the evidence does support. Note also that the disjointness argument "
        "was overstated regardless: `ContextSlot.__post_init__` validates host "
        "against table_id only, so a (HEK293, BACTERIAL_EXPRESSION) slot is "
        "constructible and would fire both rules. "
        "The evidence is otherwise strong: the brief grades both rows A, "
        "B1's is Kudla's r = 0.66 over 154 variants and B8's is Noderer's complete "
        "65,536-variant enumeration with a 12-fold range. NOTE two things a "
        "reviewer should not have to discover: no shipped preset weights 2.B8 at "
        "all, so this default only reaches a user who enables the rule outside a "
        "preset (adding it belongs to the score lane, not this one); and the "
        "weak/adequate/strong ordinal is a READING of brief.md:68, which defines "
        "only 'strong (-3 purine AND +4 G)' and leaves the lower two words "
        "undefined. Replacing the ordinal with Noderer's measured efficiency "
        "table would be strictly better and needs no change to the band or gate."
    )
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.8
    #: Zero because the engine honours no steering here, not because none is
    #: wanted: `solver/catalog.py:151-153` reads `terms.forbidden` and nothing
    #: else, so the `positional` term this rule would use is inert.
    steering_weight: ClassVar[float] = 0.0
    #: The finding is two specific bases, and both are named in the message with
    #: their own coordinates, so the interval is the context window itself.
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.WINDOW_MINUS_1
    #: SINGLE_PASS is safe: the only editable determinant is +4, one base of one
    #: codon, and setting it cannot create a second weak Kozak elsewhere. This is
    #: not the splice-donor case where a repair manufactures new instances of what
    #: it removed (CLAUDE.md 3.6).
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "cheap"
    #: brief.md:96 asks for this to be surfaced rather than silently broken: a
    #: strong Kozak ending ...CCATGG contains the NcoI site CCATGG, so raising +4
    #: to G is what creates it. D1 is HARD_LATTICE and wins when NcoI is selected.
    #: b9 declares the reverse edge already; the interaction is bidirectional --
    #: raising +4 to G strengthens the context of any downstream out-of-frame ATG,
    #: and removing one can spoil a Kozak -- so both sides say so, as d1 does.
    conflicts_with: ClassVar[tuple[str, ...]] = ("d1_restriction_sites", "b9_out_of_frame_atg")
    brief_ref: ClassVar[str] = "2.B8"
    #: A base-identity check, not a folding energy: no engine calibration applies.
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "require_strong": {
                "type": "boolean",
                "default": False,
                "description": (
                    "Report anything below 'strong' as a finding rather than only "
                    "'weak'. Off by default: an adequate context already carries "
                    "one of the two measured determinants, and flagging every "
                    "adequate start would put a finding on most real vectors."
                ),
            },
        },
    }

    def __init__(self, require_strong: bool = False) -> None:
        self.require_strong = require_strong

    def gate(self, slot: ContextSlot) -> bool:
        """Eukaryotic translating slots only -- see the module docstring.

        Both halves matter. A bacterial ribosome reads a Shine-Dalgarno, not a
        Kozak (brief.md:205), and a propagation slot never translates the
        transgene at all.
        """
        return slot.host in EUKARYOTIC_HOSTS and slot.role != "propagation"

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """None. `positional` would fit and the solver does not read it.

        `LatticeTerms.positional` is `h(codon_index, codon)` and is exactly the
        shape of "prefer a G-initial codon at index 1". `solver/catalog.py`
        consumes `terms.forbidden` alone, so supplying one would look like
        steering and do nothing. `forbidden` is not the answer either: no codon is
        illegal here, and a weak Kozak is a finding, never an unreachable state.
        """
        return None

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        slots = [s for s in ctx.active_slots if self.gate(s)]
        if not slots:
            return self._unavailable(
                c,
                "no eukaryotic translating slot in this context. A bacterial "
                "ribosome reads a Shine-Dalgarno rather than a Kozak context "
                "(brief.md:205), so there is nothing here a Kozak score describes",
            )

        editable = sorted(c.editable)
        if not editable:
            return self._unavailable(c, "no designable CDS whose start codon to read")

        scored: list[tuple[ContextSlot, Interval, str, float]] = []
        no_initiator: str | None = None
        no_initiator_at: Interval | None = None
        for slot in slots:
            strand = strand_for(ctx, slot)
            cds = editable[0] if strand == 1 else editable[-1]
            window = self._window(Interval(cds.start, cds.end, strand), c)
            if window is None:
                return self._unavailable(
                    c,
                    f"the CDS starts too close to the end of this linear construct "
                    f"for the -{UPSTREAM}..+{DOWNSTREAM} Kozak context to fit",
                )
            leader = leader_of(window, UPSTREAM)
            if not self._annotated_leader(c, leader):
                return self._unavailable(
                    c,
                    f"no annotated 5'UTR covering the {UPSTREAM} bases upstream of "
                    f"the start codon. The -3 base is only a Kozak determinant if it "
                    f"is really in the transcript leader; unannotated, it may be "
                    f"promoter or vector backbone, and reporting it as context would "
                    f"describe a transcript that is never made. Annotate the 5'UTR "
                    f"in your map and re-run",
                    interval=leader,
                )
            context = c.slice(window)
            if len(context) != UPSTREAM + DOWNSTREAM:
                return self._unavailable(
                    c, "the Kozak context window could not be read from the construct"
                )
            # A Kozak tier is a statement about an initiation event, so there has
            # to be an initiator. Without this the rule happily scores the context
            # around any triplet -- its own reverse-strand fixture lands on TTA --
            # and reports a confident "adequate" for a start that is not there.
            # `is_start` from the INJECTED table, never a hard-coded "ATG": which
            # triplets initiate is table-dependent (CLAUDE.md 3.1).
            try:
                code = svc.tables.genetic_code(slot.table_id)
            except _MISSING_TABLE as exc:
                return self._unavailable(
                    c, f"NCBI table {slot.table_id} could not be loaded: {exc}"
                )
            initiator = context[UPSTREAM : UPSTREAM + 3]
            if not code.is_start(initiator):
                # SKIP this slot; do not abandon the rule. An antisense slot reads
                # a window that begins mid-codon on its own strand, and returning
                # unavailable for the whole Evaluation would throw away a
                # perfectly readable tier for the slot that DOES initiate.
                # Unavailable only if no slot has an initiator -- see below.
                no_initiator = (
                    f"the designable CDS begins {initiator!r} for the {slot.role} "
                    f"slot ({slot.host}), which is not a start codon under NCBI "
                    f"table {slot.table_id}. A Kozak tier describes how well a "
                    f"scanning ribosome initiates AT a start codon, so scoring one "
                    f"here would describe an initiation event that cannot happen"
                )
                no_initiator_at = window
                continue
            scored.append((slot, window, context, self._tier(context)))

        if not scored:
            # `slots` was non-empty and every other per-slot failure returns, so
            # the only way to arrive here is the initiator skip above.
            assert no_initiator is not None
            return self._unavailable(c, no_initiator, interval=no_initiator_at)

        # The WEAKEST slot binds. Averaging tiers across slots would let a strong
        # context in one hide a weak one in another, and the weak one is the finding.
        bound, window, context, tier = min(scored, key=lambda row: row[3])

        breaches: list[Breach] = []
        if tier == WEAK or (self.require_strong and tier < STRONG):
            breaches.append(self._breach(bound, window, context, tier, c))

        return Evaluation(
            spec_id=self.id,
            passes=not breaches,
            raw_score=tier,
            breaches=tuple(breaches),
            windows=tuple((row[1], row[3]) for row in scored),
            # Two bases carry the score, whatever the context's length.
            n_evaluated=2,
        )

    # -- internals ----------------------------------------------------------

    def _tier(self, context: str) -> float:
        """strong = both determinants, adequate = one, weak = neither.

        brief.md:68 defines only the top tier, "(-3 purine AND +4 G)". This is the
        minimal completion of a three-level scale whose top is a conjunction of
        two terms; it is a reading, not a measurement, and `weight_provenance`
        says so.
        """
        return float(self.purine_at_minus_3(context)) + float(self.g_at_plus_4(context))

    @staticmethod
    def purine_at_minus_3(context: str) -> bool:
        """-3 is the 4th of the six leader bases in a -6..+5 context string."""
        return context[UPSTREAM - 3] in PURINES

    @staticmethod
    def g_at_plus_4(context: str) -> bool:
        """+4 is the base after ATG, i.e. the first base of codon 2."""
        return context[UPSTREAM + 3] == "G"

    @staticmethod
    def c_at_plus_5(context: str) -> bool:
        """Reported, never scored: its effect is measured relative to +4 G."""
        return context[UPSTREAM + 4] == "C"

    def _breach(
        self, slot: ContextSlot, window: Interval, context: str, tier: float, c: Construct
    ) -> Breach:
        purine, g4 = self.purine_at_minus_3(context), self.g_at_plus_4(context)
        missing = []
        if not purine:
            missing.append(f"-3 is {context[UPSTREAM - 3]}, not a purine (A or G)")
        if not g4:
            missing.append(f"+4 is {context[UPSTREAM + 3]}, not G")
        name = {STRONG: "strong", ADEQUATE: "adequate"}.get(tier, "weak")
        return Breach(
            spec_id=self.id,
            interval=window,
            magnitude=STRONG - tier,
            message=(
                f"Kozak context {context[:UPSTREAM].lower()}{context[UPSTREAM:]} for "
                f"the {slot.role} slot ({slot.host}) is {name}: {'; '.join(missing)}. "
                f"The target is gccRccATGG with +5 = C (currently "
                f"{context[UPSTREAM + 4]}). "
                f"Noderer 2014 measured a 12-fold range across this context, -3 "
                f"purine at +58% over -3U; this is a finding to weigh, "
                f"never a target to maximize"
            ),
            # +4 is codon 2 and is designable; -3 is in the user's 5'UTR and no
            # codon choice reaches it. A context weak ONLY at -3 is advice, not a
            # target for the search (PLAN.md:372).
            fixable_by_codon_choice=not g4 and c.overlaps_editable(window),
            slot_role=slot.role,
            detail={
                "tier": tier,
                "tier_name": name,
                "context": context,
                "minus_3": context[UPSTREAM - 3],
                "plus_4": context[UPSTREAM + 3],
                "plus_5": context[UPSTREAM + 4],
                "c_at_plus_5": str(self.c_at_plus_5(context)),
                "purine_at_minus_3": str(purine),
                "g_at_plus_4": str(g4),
                "host": str(slot.host),
            },
        )

    def _window(self, cds: Interval, c: Construct) -> Interval | None:
        """The -6..+5 context, or None when it does not fit.

        B1's `_window` with this rule's offsets, and deliberately not imported
        from it: sharing the method would couple two rules' geometry so that
        changing B1's Kudla window silently moved B8's Kozak context. Returns None
        rather than a clamped window, for B1's reason -- a context shortened to
        what fits is a DIFFERENT context, and -3 would no longer be at -3.
        """
        if cds.length < DOWNSTREAM:
            return None
        span = UPSTREAM + DOWNSTREAM
        start = cds.start - UPSTREAM if cds.strand == 1 else cds.end - DOWNSTREAM
        if start < 0:
            if not c.is_circular:
                return None
            start += c.length
        elif start >= c.length:
            # An origin-spanning CDS puts a minus-strand window start past the
            # end, where `Construct.slice` yields a short string. The length
            # check in `evaluate` catches it, but degrading to "unavailable" for
            # a window that exists is still the wrong answer.
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
        """NaN plus a breach carrying the reason -- B1's pattern, B1's argument.

        NaN rather than 0.0 because 0.0 is a real tier here: it is `WEAK`, the
        worst context this rule can report. Returning it for "we could not read
        the context" would put a finding in front of a user about a start codon
        nobody measured.
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
                    message=f"Kozak objective unavailable: {reason}",
                    # The missing input is an annotation, a host or a construct
                    # geometry -- never a codon the solver could have chosen.
                    fixable_by_codon_choice=False,
                    detail={"unavailable_reason": reason},
                ),
            ),
            n_evaluated=0,
        )
