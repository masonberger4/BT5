"""B9: the three clauses of brief.md:69, and the ATGs that are NOT findings.

Every sequence here is built from a GCT filler that cannot produce an ATG at any
offset, so each construct contains exactly the ATGs it is named for. A test that
accidentally carried a second ATG would pass while measuring the wrong one.
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
from bt5.core.spec import Direction, Enforcement, Evaluation, Evidence, RepairPolicy
from bt5.core.types import (
    Construct,
    Feature,
    Interval,
    Segment,
    SegmentKind,
    Strand,
    Topology,
)
from bt5.rules.catalog.b9_out_of_frame_atg import (
    HARD_MAGNITUDE,
    JUNCTION_UPSTREAM,
    SCAN_NT,
    OutOfFrameAtg,
)
from bt5.vector.kmers import ConstructKmerIndex

discover()

HEK_TABLE = 1

#: Ala. GCT repeated contains no ATG at any offset, and neither does its junction
#: with ATG at the start or TAA at the end -- so every construct below holds
#: exactly the ATGs it is named for.
FILLER = "GCT"

#: Start codon only.
CLEAN = "ATG" + FILLER * 20 + "TAA"

#: An out-of-frame ATG at CDS offset 4, inside the 50 nt window.
#: A(0)T(1)G(2) A(3)A(4)T(5) G(6)... -> s[4:7] == "ATG".
NEAR_OUT_OF_FRAME = "ATG" + "AAT" + FILLER * 19 + "TAA"

#: An IN-frame ATG at offset 3 -- a real Met at residue 2, and still an
#: alternative start that makes an N-terminally truncated protein.
NEAR_IN_FRAME_MET = "ATG" + "ATG" + FILLER * 19 + "TAA"

#: An out-of-frame ATG at offset 55, past the window, in STRONG context:
#: s[52] == "A" (-3 purine) and s[58] == "G" (+4 G).
FAR_STRONG = "ATG" + FILLER * 16 + "GAT" + "AAT" + "GGT" + "TAA"

#: The same offset 55, but s[52] == "C" and s[58] == "C": weak context.
FAR_WEAK = "ATG" + FILLER * 17 + "AAT" + FILLER + "TAA"


# -- local helpers ------------------------------------------------------------
# Defined here, not in conftest: that file is shared with the other rules session
# and is read-only for both of us (docs/buildout/README.md).


def atg_construct(
    cds: str = CLEAN,
    leader: str = "",
    *,
    circular: bool = True,
    annotate: bool = True,
) -> Construct:
    """`leader` + `cds`, with the leader annotated as 5'UTR unless told not to."""
    seq = leader + cds
    segments = [Segment(Interval(len(leader), len(seq)), SegmentKind.DESIGNABLE_CDS, "cds")]
    features: list[Feature] = []
    if leader:
        segments.insert(0, Segment(Interval(0, len(leader)), SegmentKind.BACKBONE, "utr"))
        if annotate:
            features.append(Feature(Interval(0, len(leader)), "5'UTR"))
    return Construct(
        sequence=seq,
        topology=Topology.CIRCULAR if circular else Topology.LINEAR,
        segments=tuple(segments),
        features=tuple(features),
    )


def services() -> Services:
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
    rule: OutOfFrameAtg | None = None,
) -> Evaluation:
    return (rule or OutOfFrameAtg()).evaluate(
        c if c is not None else atg_construct(), ctx or context(), services()
    )


def is_unavailable(ev: Evaluation) -> bool:
    return (
        math.isnan(ev.raw_score)
        and len(ev.breaches) == 1
        and ev.breaches[0].message.startswith("out-of-frame ATG scan unavailable:")
        and ev.n_evaluated == 0
        and ev.passes
    )


# -- the fixtures are what they claim -----------------------------------------


class TestTheFixturesThemselves:
    """If a construct carries an ATG nobody meant, every test below lies."""

    @pytest.mark.parametrize(
        ("name", "cds", "expected"),
        [
            ("CLEAN", CLEAN, [0]),
            ("NEAR_OUT_OF_FRAME", NEAR_OUT_OF_FRAME, [0, 4]),
            ("NEAR_IN_FRAME_MET", NEAR_IN_FRAME_MET, [0, 3]),
            ("FAR_STRONG", FAR_STRONG, [0, 55]),
            ("FAR_WEAK", FAR_WEAK, [0, 55]),
        ],
    )
    def test_each_construct_holds_exactly_the_atgs_it_is_named_for(
        self, name: str, cds: str, expected: list[int]
    ) -> None:
        found = [i for i in range(len(cds) - 2) if cds[i : i + 3] == "ATG"]
        assert found == expected, name

    def test_the_far_contexts_differ_only_in_strength(self) -> None:
        """Both hold an out-of-frame ATG at 55; only FAR_STRONG has -3R and +4G."""
        assert (FAR_STRONG[52], FAR_STRONG[58]) == ("A", "G")
        assert (FAR_WEAK[52], FAR_WEAK[58]) == ("C", "C")
        assert 55 % 3 != 0


# -- the clauses --------------------------------------------------------------


class TestFirstFiftyNucleotides:
    def test_a_clean_cds_has_no_findings(self) -> None:
        ev = evaluate(atg_construct(CLEAN))
        assert ev.passes
        assert ev.raw_score == 0.0

    def test_the_start_codon_is_not_an_additional_atg(self) -> None:
        """Offset 0 is the thing the rule protects, not a breach."""
        assert evaluate(atg_construct(CLEAN)).breaches == ()

    def test_an_out_of_frame_atg_in_the_window_is_a_breach(self) -> None:
        ev = evaluate(atg_construct(NEAR_OUT_OF_FRAME))
        assert not ev.passes
        assert ev.raw_score == 1.0
        assert ev.breaches[0].detail["cds_offset"] == 4.0
        assert ev.breaches[0].magnitude == HARD_MAGNITUDE

    def test_an_in_frame_atg_in_the_window_is_also_a_breach(self) -> None:
        """brief.md:69 says "No ADDITIONAL ATG", not "no out-of-frame ATG".

        An in-frame ATG at codon 2 is a perfectly good alternative start that
        makes an N-terminally truncated protein.
        """
        ev = evaluate(atg_construct(NEAR_IN_FRAME_MET))
        assert not ev.passes
        assert ev.breaches[0].detail["in_frame"] == "True"

    def test_a_forced_met_is_advice_not_a_target(self) -> None:
        """Met has one codon under table 1, so no recoding removes this ATG.

        Sending the solver after it drives a search into an infeasibility it
        cannot resolve (docs/PLAN.md:372).
        """
        ev = evaluate(atg_construct(NEAR_IN_FRAME_MET))
        assert ev.breaches[0].detail["forced_met"] == "True"
        assert not ev.breaches[0].fixable_by_codon_choice

    def test_an_out_of_frame_atg_is_the_solvers_to_fix(self) -> None:
        ev = evaluate(atg_construct(NEAR_OUT_OF_FRAME))
        assert ev.breaches[0].detail["forced_met"] == "False"
        assert ev.breaches[0].fixable_by_codon_choice

    def test_the_window_is_configurable_and_defaults_to_the_briefs_50(self) -> None:
        assert SCAN_NT == 50
        assert OutOfFrameAtg().scan_nt == 50
        # Shrink it below the offending offset and the near-start clause stops
        # firing; the ATG at 55 is then weak-context and not reported either.
        assert evaluate(atg_construct(FAR_WEAK), rule=OutOfFrameAtg(scan_nt=6)).passes

    def test_scan_nt_must_be_at_least_one_codon(self) -> None:
        with pytest.raises(ValueError, match="at least one codon"):
            OutOfFrameAtg(scan_nt=2)


class TestStrongContextBeyondTheWindow:
    def test_a_strong_context_out_of_frame_atg_is_reported(self) -> None:
        """brief.md:69: "penalize any out-of-frame ATG anywhere that has BOTH
        -3 purine and +4 G"."""
        ev = evaluate(atg_construct(FAR_STRONG))
        assert not ev.passes
        assert ev.breaches[0].detail["strong_context"] == "True"
        assert ev.breaches[0].detail["cds_offset"] == 55.0

    def test_a_weak_context_out_of_frame_atg_is_not(self) -> None:
        """The brief penalises strong-context ones specifically.

        Reporting every out-of-frame ATG would put dozens of findings on any real
        CDS and bury the ones a ribosome would actually use.
        """
        ev = evaluate(atg_construct(FAR_WEAK))
        assert ev.passes
        assert ev.breaches == ()

    def test_it_carries_less_weight_than_a_first_window_hit(self) -> None:
        near = evaluate(atg_construct(NEAR_OUT_OF_FRAME)).breaches[0].magnitude
        far = evaluate(atg_construct(FAR_STRONG)).breaches[0].magnitude
        assert far < near

    def test_turning_off_the_context_filter_reports_every_one(self) -> None:
        rule = OutOfFrameAtg(strong_context_only=False)
        assert not evaluate(atg_construct(FAR_WEAK), rule=rule).passes

    def test_an_in_frame_atg_beyond_the_window_is_never_a_finding(self) -> None:
        """It is an ordinary internal Met, not an alternative reading frame."""
        cds = "ATG" + FILLER * 17 + "ATG" + FILLER + "TAA"
        assert cds[54:57] == "ATG"
        assert evaluate(atg_construct(cds)).passes


class TestTheJunctionClause:
    #: A leader whose last three bases are ATG, ending exactly at the junction.
    LEADER = "GCCATG"

    def test_an_upstream_atg_is_seen_when_the_leader_is_annotated(self) -> None:
        """brief.md:69's third clause. Rules take a Construct precisely so a
        junction-region finding cannot be missed (CLAUDE.md 3.3)."""
        c = atg_construct(CLEAN, leader=self.LEADER)
        ev = evaluate(c, rule=OutOfFrameAtg(strong_context_only=False))
        upstream = [b for b in ev.breaches if b.detail["upstream_of_cds"] == "True"]
        assert len(upstream) == 1
        assert upstream[0].detail["cds_offset"] == -3.0

    def test_the_cds_start_codon_is_still_not_a_finding(self) -> None:
        c = atg_construct(CLEAN, leader=self.LEADER)
        ev = evaluate(c, rule=OutOfFrameAtg(strong_context_only=False))
        assert all(b.detail["cds_offset"] != 0.0 for b in ev.breaches)

    def test_an_unannotated_leader_is_not_scanned_and_says_so(self) -> None:
        """Unannotated upstream sequence may be promoter or backbone.

        An ATG there is not one a ribosome scanning THIS transcript would meet,
        so it is not reported -- but the gap is stated rather than silent.
        """
        c = atg_construct(NEAR_OUT_OF_FRAME, leader=self.LEADER, annotate=False)
        ev = evaluate(c, rule=OutOfFrameAtg(strong_context_only=False))
        assert all(b.detail["upstream_of_cds"] == "False" for b in ev.breaches)
        assert ev.breaches[0].detail["junction_checked"] == "False"
        assert "junction was not scanned" in ev.breaches[0].message

    def test_the_cds_clauses_still_run_without_a_leader(self) -> None:
        """Degrading to a partial scan beats skipping the rule outright."""
        c = atg_construct(NEAR_OUT_OF_FRAME, leader=self.LEADER, annotate=False)
        ev = evaluate(c)
        assert not ev.passes
        assert not is_unavailable(ev)

    def test_the_leader_width_matches_b8s_context(self) -> None:
        assert JUNCTION_UPSTREAM == 6


# -- gating and strand --------------------------------------------------------


class TestGating:
    def test_bacteria_use_a_tir_model_not_an_atg_scan(self) -> None:
        """brief.md:67 (B7) owns internal initiation in bacteria."""
        assert not OutOfFrameAtg().gate(
            slot(host=HostId.E_COLI_K12, modality=Modality.BACTERIAL_EXPRESSION, table=11)
        )

    def test_a_propagation_slot_never_translates_the_transgene(self) -> None:
        assert not OutOfFrameAtg().gate(slot(role="propagation"))

    def test_a_bacterial_only_context_is_unavailable_not_clean(self) -> None:
        """raw_score 0.0 means "scanned, found nothing" -- the one claim a
        skipped scan must not make."""
        ctx = context(
            slot(host=HostId.E_COLI_K12, modality=Modality.BACTERIAL_EXPRESSION, table=11)
        )
        ev = evaluate(atg_construct(NEAR_OUT_OF_FRAME), ctx)
        assert is_unavailable(ev)
        assert ev.raw_score != 0.0


class TestStrand:
    """CLAUDE.md 3.4: directional, so the scan follows strand_for, not strand 1."""

    def test_the_reverse_strand_reads_the_other_sequence(self) -> None:
        """The reverse complement of a clean CDS is a different molecule.

        GCT filler reverse-complements to AGC, and the TAA stop to TTA, so the
        minus-strand reading has its own ATG population -- which is exactly why
        this rule must not be closed under reverse complement.
        """
        forward = evaluate(atg_construct(NEAR_OUT_OF_FRAME), context(slot()))
        reverse = evaluate(atg_construct(NEAR_OUT_OF_FRAME), context(slot(), orientation=-1))
        # The minus strand of this CDS reads TTAAGCAGCAGC... and holds no ATG at
        # all, so a rule closed under reverse complement would report a finding
        # here that no ribosome in this construct could ever act on.
        assert forward.raw_score == 1.0
        assert not forward.passes
        assert reverse.raw_score == 0.0
        assert reverse.passes

    def test_two_negatives_compose_to_a_forward_read(self) -> None:
        both = context(slot(strand_of_interest=-1), orientation=-1)
        forward = evaluate(atg_construct(NEAR_OUT_OF_FRAME), context(slot()))
        assert evaluate(atg_construct(NEAR_OUT_OF_FRAME), both).raw_score == forward.raw_score


# -- metadata -----------------------------------------------------------------


class TestSpecMetadata:
    def test_it_is_hard_and_carries_no_objective_weight(self) -> None:
        """CLAUDE.md 3.5: hard constraints are never enforced by a penalty weight."""
        assert OutOfFrameAtg.enforcement is Enforcement.HARD_REPAIR
        assert OutOfFrameAtg.enforcement.is_hard
        assert OutOfFrameAtg.default_weight == 0.0

    def test_repair_is_fixed_point_and_that_is_mandatory(self) -> None:
        """CLAUDE.md 3.6. Removing one ATG can create another at a neighbouring
        offset, and a single pass would ship a construct whose ATGs were removed
        INTO new ATGs with the validator passing it."""
        assert OutOfFrameAtg.repair is RepairPolicy.FIXED_POINT

    def test_it_declares_no_forbidden_motifs(self) -> None:
        """`ATG` as a forbidden string would ban the start codon and every
        in-frame Met, and the solver closes `forbidden` under reverse complement,
        so it would ban every CAT as well."""
        assert OutOfFrameAtg().lattice_terms(context()) is None

    def test_it_is_lower_is_better_with_no_band(self) -> None:
        assert OutOfFrameAtg.direction is Direction.LOWER_IS_BETTER
        assert OutOfFrameAtg.band is None

    def test_the_evidence_badge_matches_the_briefs_grade(self) -> None:
        assert OutOfFrameAtg.evidence is Evidence.EVIDENCE_BACKED

    def test_the_scan_length_is_reported_honestly(self) -> None:
        ev = evaluate(atg_construct(CLEAN))
        assert ev.n_evaluated == len(CLEAN)


def test_it_is_registered_under_its_brief_row() -> None:
    assert get("b9_out_of_frame_atg") is OutOfFrameAtg
    assert OutOfFrameAtg.brief_ref == "2.B9"
    assert OutOfFrameAtg.id == Path(__file__).stem.removeprefix("test_")


def test_no_preset_may_weight_a_hard_rule() -> None:
    """presets.resolve raises PresetError if one does; this catches it earlier."""
    from bt5.score.presets import PRESETS

    assert not [
        entry for preset in PRESETS for entry in preset.entries if entry.brief_ref == "2.B9"
    ]
