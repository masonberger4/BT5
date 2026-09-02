"""B8: Kozak context, and the two determinants brief.md:68 actually defines.

The tier is checked against the brief's own conjunction -- "strong (-3 purine AND
+4 G)" -- and against Noderer 2014's finding that -3 is a PURINE class rather than
an A, because "-3 must be A" is the plausible misreading that would pass every
mechanical check and silently downgrade half the strong contexts in the wild.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
from bt5.codon.tables import FileTableProvider
from bt5.core.context import (
    BiosecurityVerdict,
    ContextSlot,
    DesignContext,
    HostId,
    Modality,
)
from bt5.core.registry import discover, get
from bt5.core.services import Services
from bt5.core.spec import Direction, Enforcement, Evaluation, Evidence
from bt5.core.types import (
    Construct,
    Feature,
    Interval,
    Segment,
    SegmentKind,
    Strand,
    Topology,
)
from bt5.rules.catalog.b8_kozak import (
    ADEQUATE,
    DOWNSTREAM,
    EUKARYOTIC_HOSTS,
    STRONG,
    UPSTREAM,
    WEAK,
    KozakContext,
)
from bt5.vector.kmers import ConstructKmerIndex

discover()

#: HEK293's locked NCBI table. A mismatched (host, table) pair raises in
#: ContextSlot.__post_init__, so these travel together everywhere below.
HEK_TABLE = 1

#: A strong context: -3 = A (purine), +4 = G, +5 = C. This is brief.md:68's
#: GCCRCCATGG target with the +5 = C instruction applied.
STRONG_LEADER = "GCCACC"
STRONG_CDS = "ATGGCTGAAGGTATCAAATAA"

#: Ends in CAT, so the MINUS-strand transcript begins ATG. 18 nt, placed at
#: [6, 24) by `kozak_construct` with the default 6 nt leader.
REVERSE_START_CDS = "ATGGCTGAAGGTATCCAT"


# -- local helpers ------------------------------------------------------------
# Defined here, not in conftest: that file is shared with the other rules session
# and is read-only for both of us (docs/buildout/README.md).


def kozak_construct(
    leader: str = STRONG_LEADER,
    cds: str = STRONG_CDS,
    *,
    circular: bool = True,
    annotate: bool = True,
    trailer: str = "",
    extra_utr: Interval | None = None,
) -> Construct:
    """`leader` + `cds` + `trailer`, with the leader annotated as 5'UTR.

    conftest's `construct()` builds no features, and B8 reports unavailable
    without an annotated leader, so every test here needs this instead.
    """
    seq = leader + cds + trailer
    segments = [
        Segment(Interval(len(leader), len(leader) + len(cds)), SegmentKind.DESIGNABLE_CDS, "cds")
    ]
    if leader:
        segments.insert(0, Segment(Interval(0, len(leader)), SegmentKind.BACKBONE, "utr"))
    if trailer:
        segments.append(
            Segment(Interval(len(leader) + len(cds), len(seq)), SegmentKind.BACKBONE, "vector")
        )
    features: list[Feature] = []
    if annotate and leader:
        features.append(Feature(Interval(0, len(leader)), "5'UTR"))
    if extra_utr is not None:
        features.append(Feature(extra_utr, "five_prime_UTR"))
    return Construct(
        sequence=seq,
        topology=Topology.CIRCULAR if circular else Topology.LINEAR,
        segments=tuple(segments),
        features=tuple(features),
    )


def services() -> Services:
    # A real provider: B8 asks the INJECTED table whether +1..+3 initiates, so a
    # None here would only ever exercise the error path.
    return Services(
        fold=None,
        kmer=ConstructKmerIndex,
        tables=FileTableProvider(),  # type: ignore[arg-type]
        rng=np.random.default_rng(42),
    )


def slot(
    role: str = "producer",
    host: HostId = HostId.HEK293,
    modality: Modality = Modality.LENTIVIRAL,
    table: int = HEK_TABLE,
    strand_of_interest: Strand = 1,
) -> ContextSlot:
    return ContextSlot(role, host, modality, table, strand_of_interest)  # type: ignore[arg-type]


def context(*slots: ContextSlot, orientation: Strand = 1) -> DesignContext:
    return DesignContext(
        slots=slots or (slot(),),
        cassette_orientation=orientation,
        seed=42,
        screen=BiosecurityVerdict("not_run"),
    )


def evaluate(
    c: Construct | None = None,
    ctx: DesignContext | None = None,
    rule: KozakContext | None = None,
) -> Evaluation:
    return (rule or KozakContext()).evaluate(
        c if c is not None else kozak_construct(), ctx or context(), services()
    )


def is_unavailable(ev: Evaluation) -> bool:
    return (
        math.isnan(ev.raw_score)
        and len(ev.breaches) == 1
        and ev.breaches[0].message.startswith("Kozak objective unavailable:")
        and not ev.breaches[0].fixable_by_codon_choice
        and ev.n_evaluated == 0
        and ev.passes
    )


# -- the tiers ----------------------------------------------------------------


class TestTier:
    @pytest.mark.parametrize("base", ["A", "G"])
    def test_minus_3_is_a_purine_class_not_the_letter_a(self, base: str) -> None:
        """brief.md:68 writes `R` at -3, and Noderer's +58% is purine vs U.

        "-3 must be A" is the misreading this test exists to stop: it would
        downgrade every G-at-minus-3 context, which is half of all strong ones.
        """
        ev = evaluate(kozak_construct(leader=f"GCC{base}CC"))
        assert ev.raw_score == STRONG

    @pytest.mark.parametrize("base", ["C", "T"])
    def test_a_pyrimidine_at_minus_3_loses_that_determinant(self, base: str) -> None:
        ev = evaluate(kozak_construct(leader=f"GCC{base}CC"))
        assert ev.raw_score == ADEQUATE

    def test_strong_needs_both_determinants(self) -> None:
        assert evaluate(kozak_construct()).raw_score == STRONG

    def test_losing_plus_4_drops_a_tier(self) -> None:
        """+4 is the first base of codon 2, so ATG then C rather than G."""
        ev = evaluate(kozak_construct(cds="ATGCCTGAAGGTATCAAATAA"))
        assert ev.raw_score == ADEQUATE

    def test_losing_both_is_weak(self) -> None:
        ev = evaluate(kozak_construct(leader="GCCTCC", cds="ATGCCTGAAGGTATCAAATAA"))
        assert ev.raw_score == WEAK
        assert not ev.passes

    def test_plus_5_is_reported_but_never_scored(self) -> None:
        """Noderer's 24.8% is measured relative to +4G, so scoring it double-counts.

        Two constructs differing only at +5 must land on the same tier.
        """
        with_c = evaluate(kozak_construct(cds="ATGGCTGAAGGTATCAAATAA"))
        with_a = evaluate(kozak_construct(cds="ATGGATGAAGGTATCAAATAA"))
        assert with_c.raw_score == with_a.raw_score == STRONG

    def test_an_adequate_context_is_not_a_finding_by_default(self) -> None:
        """Flagging every adequate start would put a finding on most real vectors."""
        ev = evaluate(kozak_construct(cds="ATGCCTGAAGGTATCAAATAA"))
        assert ev.passes
        assert ev.breaches == ()

    def test_require_strong_makes_adequate_a_finding(self) -> None:
        ev = evaluate(
            kozak_construct(cds="ATGCCTGAAGGTATCAAATAA"), rule=KozakContext(require_strong=True)
        )
        assert not ev.passes
        assert ev.breaches[0].magnitude == STRONG - ADEQUATE


class TestFixability:
    def test_a_missing_plus_4_is_the_solvers_to_fix(self) -> None:
        """+4 is codon 2 and designable."""
        ev = evaluate(kozak_construct(leader="GCCTCC", cds="ATGCCTGAAGGTATCAAATAA"))
        assert ev.breaches[0].fixable_by_codon_choice

    def test_a_context_weak_only_at_minus_3_is_advice_not_a_target(self) -> None:
        """-3 lives in the user's 5'UTR; no codon choice reaches it.

        Sending the search after it is how a solver chases spurious infeasibility
        (docs/PLAN.md:372).
        """
        ev = evaluate(kozak_construct(leader="GCCTCC"), rule=KozakContext(require_strong=True))
        assert ev.raw_score == ADEQUATE
        assert not ev.breaches[0].fixable_by_codon_choice

    def test_the_breach_names_what_is_wrong_and_where(self) -> None:
        ev = evaluate(kozak_construct(leader="GCCTCC", cds="ATGCCTGAAGGTATCAAATAA"))
        message = ev.breaches[0].message
        assert "-3 is T, not a purine" in message
        assert "+4 is C, not G" in message
        assert ev.breaches[0].detail["minus_3"] == "T"
        assert ev.breaches[0].detail["plus_4"] == "C"

    def test_the_breach_never_predicts_expression(self) -> None:
        """CLAUDE.md: BT5 never reports a predicted expression level."""
        message = (
            evaluate(kozak_construct(leader="GCCTCC", cds="ATGCCTGAAGGTATCAAATAA"))
            .breaches[0]
            .message.lower()
        )
        for banned in ("will increase", "predicted expression", "fold-improvement"):
            assert banned not in message


# -- gating -------------------------------------------------------------------


class TestGating:
    @pytest.mark.parametrize("host", sorted(EUKARYOTIC_HOSTS))
    def test_every_eukaryotic_host_is_scored(self, host: HostId) -> None:
        from bt5.core.context import LOCKED_TRANSLATION_TABLE

        assert KozakContext().gate(slot(host=host, table=LOCKED_TRANSLATION_TABLE[host]))

    @pytest.mark.parametrize("host", [HostId.E_COLI_K12, HostId.E_COLI_BL21])
    def test_bacteria_read_a_shine_dalgarno_not_a_kozak(self, host: HostId) -> None:
        """brief.md:205: "Kozak (B8/B9) | No (SD instead)". B6/B7 own that path.

        Scoring Kozak here would not be a number measured on the wrong host -- it
        would be a number measuring nothing.
        """
        assert not KozakContext().gate(
            slot(host=host, modality=Modality.BACTERIAL_EXPRESSION, table=11)
        )

    def test_a_propagation_slot_never_translates_the_transgene(self) -> None:
        assert not KozakContext().gate(slot(role="propagation"))

    def test_a_bacterial_only_context_is_unavailable_not_weak(self) -> None:
        """Reporting WEAK (0.0) would claim a measurement nobody made."""
        ctx = context(
            slot(host=HostId.E_COLI_K12, modality=Modality.BACTERIAL_EXPRESSION, table=11)
        )
        ev = evaluate(ctx=ctx)
        assert is_unavailable(ev)
        assert "Shine-Dalgarno" in ev.breaches[0].message

    def test_a_second_eukaryotic_slot_does_not_mask_the_finding(self) -> None:
        """`min` across slots, not a mean: a strong slot must not average a weak one away.

        Both slots read the same context here, so this pins the aggregation shape
        rather than the tie-break; the strand case in TestStrand is where the two
        actually diverge.
        """
        ev = evaluate(
            kozak_construct(leader="GCCTCC", cds="ATGCCTGAAGGTATCAAATAA"),
            # Distinct roles: DesignContext.__post_init__ rejects duplicates.
            context(
                slot(role="producer", host=HostId.HEK293),
                slot(role="target", host=HostId.CHO, table=HEK_TABLE),
            ),
        )
        assert ev.raw_score == WEAK
        assert len(ev.windows) == 2


# -- unavailability -----------------------------------------------------------


class TestUnavailable:
    def test_an_unannotated_leader_is_unavailable(self) -> None:
        """-3 is only a determinant if it is really in the transcript leader.

        Unannotated it may be promoter or backbone, and scoring it would describe
        a transcript that is never made -- B1's argument for the same refusal.
        """
        ev = evaluate(kozak_construct(annotate=False))
        assert is_unavailable(ev)
        assert "no annotated 5'UTR" in ev.breaches[0].message

    def test_a_cds_too_close_to_a_linear_end_is_unavailable(self) -> None:
        """A context clamped to what fits is a DIFFERENT context: -3 moves."""
        ev = evaluate(kozak_construct(leader="GC", circular=False))
        assert is_unavailable(ev)

    def test_the_same_cds_wraps_fine_when_circular(self) -> None:
        """A plasmid has no end, so the leader is upstream across the origin."""
        c = kozak_construct(leader="GC", circular=True)
        # The two-base leader cannot cover -6..-1, so this is still unavailable --
        # but for the ANNOTATION reason, not the geometry one.
        ev = evaluate(c)
        assert is_unavailable(ev)
        assert "too close to the end" not in ev.breaches[0].message

    def test_no_designable_cds_is_unavailable(self) -> None:
        c = Construct(sequence="ACGT" * 10, topology=Topology.CIRCULAR, segments=())
        ev = evaluate(c)
        assert is_unavailable(ev)
        assert "no designable CDS" in ev.breaches[0].message

    def test_unavailable_is_not_a_weak_context(self) -> None:
        """NaN, not 0.0. WEAK is a real tier -- the worst this rule can report --
        so returning it for "we could not read the context" would put a finding in
        front of a user about a start codon nobody measured."""
        ev = evaluate(kozak_construct(annotate=False))
        assert math.isnan(ev.raw_score)
        assert ev.raw_score != WEAK
        assert ev.passes is True


class TestStrand:
    """CLAUDE.md 3.4: a directional model is not revcomp-symmetric.

    On the minus strand the transcript's 5' end is at HIGHER construct
    coordinates, so the leader is the TAIL of the context window. Getting that
    backwards checks the wrong six bases for UTR annotation and then scores a
    -3 that is really the far end of the CDS -- b1's `leader_of` exists for
    exactly this and B8 reuses it rather than re-deriving it.
    """

    #: Linear, so the reverse-strand leader lands in the trailer instead of
    #: wrapping round to the annotated 5' leader and hiding the bug.
    def _flanked(self, **kw: object) -> Construct:
        return kozak_construct(trailer="TTTTTTTTTT", circular=False, **kw)  # type: ignore[arg-type]

    def test_the_forward_strand_reads_the_annotated_leader(self) -> None:
        ev = evaluate(self._flanked(), context(slot()))
        assert ev.raw_score == STRONG

    def test_the_reverse_strand_looks_at_the_other_end(self) -> None:
        """The 5'UTR annotation at the low end is the wrong end for a minus-strand
        transcript, so the context is unavailable rather than silently scored."""
        ev = evaluate(self._flanked(), context(slot(), orientation=-1))
        assert is_unavailable(ev)
        assert "no annotated 5'UTR" in ev.breaches[0].message

    def test_the_reverse_strand_reads_a_real_start_codon_on_its_own_strand(self) -> None:
        """The proof the rule reads the TAIL, with a falsifiable value.

        This used to assert only `raw_score in (WEAK, ADEQUATE, STRONG)`, which
        `_tier` cannot violate -- it proved the call did not fail, not that the
        right bases were read. Worse, its construct put the minus-strand "start"
        on TTA, so it was asserting a Kozak tier for a start codon that is not
        one; the rule now refuses that outright.

        Here the CDS ENDS in CAT, so the minus-strand transcript begins ATG. The
        CDS is [6, 24), the window is [19, 30) and the minus-strand leader is
        [24, 30) -- inside the trailer, which is where a minus-strand
        transcript's 5' end sits. The context reads AAAAAAATGGA: -3 is A (purine)
        and +4 is G, so STRONG, and nothing but reading the correct tail in the
        correct direction produces it.
        """
        c = kozak_construct(
            cds=REVERSE_START_CDS,
            trailer="T" * 10,
            circular=False,
            extra_utr=Interval(24, 30),
        )
        ev = evaluate(c, context(slot(), orientation=-1))
        assert not is_unavailable(ev)
        assert ev.raw_score == STRONG
        assert ev.breaches == ()

    def test_a_window_whose_start_is_not_a_start_codon_is_unavailable(self) -> None:
        """A Kozak tier is a claim about initiation, so there must be an initiator.

        Without this the rule scores the context around any triplet -- the old
        reverse-strand fixture landed on TTA -- and reports a confident tier for
        a start that is not there.
        """
        # Leader annotated at the HIGH end, so the annotation path is satisfied
        # and only the initiator is wrong: the default CDS ends TAA, whose
        # minus-strand read is TTA. This is the exact fixture the old
        # tautological test scored as "adequate".
        c = kozak_construct(trailer="T" * 10, circular=False, extra_utr=Interval(27, 33))
        ev = evaluate(c, context(slot(), orientation=-1))
        assert is_unavailable(ev)
        assert "not a start codon" in ev.breaches[0].message

    def test_the_cassette_orientation_composes_with_the_slot(self) -> None:
        """`strand_for` multiplies the two; two negatives are a forward read."""
        both = context(slot(strand_of_interest=-1), orientation=-1)
        assert evaluate(self._flanked(), both).raw_score == STRONG


# -- metadata -----------------------------------------------------------------


class TestSpecMetadata:
    def test_it_is_a_soft_objective_that_never_steers(self) -> None:
        """solver/catalog.py reads terms.forbidden and nothing else.

        `LatticeTerms.positional` is exactly the right shape for "prefer a
        G-initial codon at index 1" and is consumed by nothing, so a non-zero
        steering_weight would claim a nudge the engine does not perform.
        """
        assert KozakContext.enforcement is Enforcement.SOFT
        assert KozakContext.direction is Direction.HIGHER_IS_BETTER
        assert KozakContext.steering_weight == 0.0
        assert KozakContext().lattice_terms(context()) is None

    def test_it_declares_no_band(self) -> None:
        """Only Direction.BAND requires one, and this is an ordinal."""
        assert KozakContext.band is None

    def test_the_evidence_badge_matches_the_briefs_grade(self) -> None:
        """brief.md:68 grades B8 'A'; Noderer enumerated all 65,536 variants."""
        assert KozakContext.evidence is Evidence.EVIDENCE_BACKED

    def test_it_surfaces_the_ncoi_conflict_rather_than_breaking_it(self) -> None:
        """brief.md:96 asks for exactly this.

        A strong Kozak ending ...CCATGG contains NcoI, so raising +4 to G creates
        the site. D1 is HARD_LATTICE and wins; the conflict panel needs to know.
        """
        assert "d1_restriction_sites" in KozakContext.conflicts_with

    def test_the_b9_conflict_is_declared_from_both_sides(self) -> None:
        """b9 already declared the reverse edge; d1/b8 is symmetric too.

        Raising +4 to G strengthens the context of any downstream out-of-frame
        ATG, and removing one can spoil a Kozak -- the interaction runs both ways
        and the conflict panel reads both tuples.
        """
        from bt5.rules.catalog.b9_out_of_frame_atg import OutOfFrameAtg

        assert "b9_out_of_frame_atg" in KozakContext.conflicts_with
        assert "b8_kozak" in OutOfFrameAtg.conflicts_with

    def test_it_is_second_in_the_catalog_and_never_ties_b1(self) -> None:
        """The ordering the weight was lowered from 1.0 to establish.

        `test_b1_five_prime.py` asserts B1 outranks everything; nothing asserted
        where B8 then sits, so the "second only to B1" claim in
        `weight_provenance` could have drifted silently. Both halves are pinned
        here: strictly below B1, and strictly above every other spec.
        """
        from bt5.core.registry import all_specs
        from bt5.rules.catalog.b1_five_prime import FivePrimeFolding

        assert KozakContext.default_weight == 0.8
        assert KozakContext.default_weight < FivePrimeFolding.default_weight
        rest = [
            s.default_weight
            for s in all_specs()
            if s.id not in (KozakContext.id, FivePrimeFolding.id)
        ]
        assert KozakContext.default_weight > max(rest)

    def test_the_window_offsets_match_noderers_measured_context(self) -> None:
        assert (UPSTREAM, DOWNSTREAM) == (6, 5)

    def test_the_context_window_is_reported_for_every_scored_slot(self) -> None:
        ev = evaluate()
        assert len(ev.windows) == 1
        assert ev.windows[0][0].length == UPSTREAM + DOWNSTREAM
        assert ev.n_evaluated == 2


def test_it_is_registered_under_its_brief_row() -> None:
    assert get("b8_kozak") is KozakContext
    assert KozakContext.brief_ref == "2.B8"


def test_no_preset_weights_it_yet() -> None:
    """Recorded, not asserted as desirable: adding it belongs to the score lane.

    If this starts failing, a preset has taken 2.B8 on and this rule's
    default_weight is suddenly load-bearing -- which is the moment to re-read
    weight_provenance rather than the moment to delete the test.
    """
    from bt5.score.presets import PRESETS

    assert not [
        entry for preset in PRESETS for entry in preset.entries if entry.brief_ref == "2.B8"
    ]


def test_the_source_matches_its_filename() -> None:
    assert KozakContext.id == Path(__file__).stem.removeprefix("test_")
