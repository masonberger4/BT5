"""B2: Cambray's two windows, and the refusal to average them into one number.

The fold engine here is a deterministic stub, not ViennaRNA. That is deliberate
twice over: these tests are about WINDOW GEOMETRY and aggregation, which are the
parts a real engine would obscure rather than check, and CLAUDE.md 6 forbids
putting a real dG in a byte-exact expectation -- a ViennaRNA bump is a scientific
change, and a test asserting its output would turn that into a red build in the
rules lane instead of the calibration one.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
from bt5.core.context import (
    BiosecurityVerdict,
    ContextSlot,
    DesignContext,
    HostId,
    Modality,
)
from bt5.core.registry import discover, get
from bt5.core.services import FoldEnergy, Services
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
from bt5.rules.catalog.b1_five_prime import ENGINE_CALIBRATION
from bt5.rules.catalog.b2_structure_windows import (
    DISTAL_LENGTH,
    DISTAL_START,
    PROXIMAL_DOWNSTREAM,
    PROXIMAL_UPSTREAM,
    StructureWindows,
)
from bt5.vector.kmers import ConstructKmerIndex

discover()

ECOLI_TABLE = 11

#: Long enough that STR(+31:+90) fits: 30 nt of leader plus 120 nt of CDS.
LEADER = "GC" * 15
CDS = "ATG" + "GCT" * 38 + "TAA"


class StubFold:
    """A fold engine whose dG is a function of the window, not of physics.

    Returns `-1.0 * iv.start - iv.length / 100`, so every window has a distinct,
    predictable energy and a test can prove WHICH interval was folded. Its
    identity matches ENGINE_CALIBRATION so `check_engine_calibration` is happy.
    """

    name = "viennarna"
    version = "2.7.2"
    param_set = "rna_turner2004"

    def __init__(self) -> None:
        self.folded: list[Interval] = []

    def _energy(self, iv: Interval) -> FoldEnergy:
        self.folded.append(iv)
        return FoldEnergy(
            dg_kcal_mol=-1.0 * iv.start - iv.length / 100,
            engine=self.name,
            engine_version=self.version,
            param_set=self.param_set,
        )

    def mfe(self, seq: str) -> FoldEnergy:
        return self._energy(Interval(0, len(seq)))

    def mfe_window(self, seq: str, iv: Interval) -> FoldEnergy:
        return self._energy(iv)

    def accessibility(self, seq: str, iv: Interval, u: int) -> float | None:
        return None

    def duplex(self, a: str, b: str) -> FoldEnergy:
        return self._energy(Interval(0, max(1, len(a))))


# -- local helpers ------------------------------------------------------------
# Defined here, not in conftest: that file is shared with the other rules session
# and is read-only for both of us (docs/buildout/README.md).


def folded_construct(
    cds: str = CDS,
    leader: str = LEADER,
    *,
    circular: bool = True,
    annotate: bool = True,
) -> Construct:
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


def services(fold: StubFold | None = None) -> Services:
    return Services(
        fold=fold if fold is not None else StubFold(),  # type: ignore[arg-type]
        kmer=ConstructKmerIndex,
        tables=None,  # type: ignore[arg-type]
        rng=np.random.default_rng(42),
    )


def slot(
    role: str = "producer",
    host: HostId = HostId.E_COLI_K12,
    modality: Modality = Modality.BACTERIAL_EXPRESSION,
    table: int = ECOLI_TABLE,
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
    svc: Services | None = None,
    rule: StructureWindows | None = None,
) -> Evaluation:
    return (rule or StructureWindows()).evaluate(
        c if c is not None else folded_construct(), ctx or context(), svc or services()
    )


def is_unavailable(ev: Evaluation) -> bool:
    return (
        math.isnan(ev.raw_score)
        and len(ev.breaches) == 1
        and ev.breaches[0].message.startswith("5' structure window objective unavailable:")
        and not ev.breaches[0].fixable_by_codon_choice
        and ev.n_evaluated == 0
        and ev.passes
    )


# -- window geometry ----------------------------------------------------------


class TestWindows:
    def test_the_proximal_window_spans_the_junction(self) -> None:
        """STR(-30:+30): 30 nt of leader and 30 of CDS, so it crosses the junction.

        This is why the rule needs a Construct and an annotated 5'UTR rather than
        a CDS string (CLAUDE.md 3.3).
        """
        fold = StubFold()
        evaluate(svc=services(fold))
        proximal = fold.folded[0]
        assert proximal.start == len(LEADER) - PROXIMAL_UPSTREAM
        assert proximal.length == PROXIMAL_UPSTREAM + PROXIMAL_DOWNSTREAM == 60

    def test_the_distal_window_is_entirely_inside_the_cds(self) -> None:
        """STR(+31:+90) needs no leader: it starts 30 nt into the coding sequence."""
        fold = StubFold()
        evaluate(svc=services(fold))
        distal = fold.folded[1]
        assert distal.start == len(LEADER) + DISTAL_START
        assert distal.length == DISTAL_LENGTH == 60

    def test_the_two_windows_do_not_overlap(self) -> None:
        """Cambray localised the effect to two spans, not one 150 nt smear."""
        fold = StubFold()
        evaluate(svc=services(fold))
        proximal, distal = fold.folded[0], fold.folded[1]
        assert proximal.end == distal.start

    def test_a_short_cds_has_no_distal_window_and_that_is_not_a_failure(self) -> None:
        """An ORF under 90 nt simply has nowhere for STR(+31:+90) to sit.

        The proximal measurement must still be reported -- making the whole
        objective unavailable because the far window did not fit would lose a
        number that was measured perfectly well.
        """
        fold = StubFold()
        short = folded_construct(cds="ATG" + "GCT" * 15 + "TAA")
        ev = evaluate(short, svc=services(fold))
        assert not is_unavailable(ev)
        assert not math.isnan(ev.raw_score)
        assert len(fold.folded) == 1
        assert len(ev.windows) == 1

    def test_report_distal_off_skips_the_second_fold(self) -> None:
        fold = StubFold()
        evaluate(svc=services(fold), rule=StructureWindows(report_distal=False))
        assert len(fold.folded) == 1


# -- the refusal to average ---------------------------------------------------


class TestTheTwoWindowsAreNotAveraged:
    """brief.md:333: "never average them into one 'structure' slider"."""

    def test_the_score_is_the_proximal_window_alone(self) -> None:
        fold = StubFold()
        ev = evaluate(svc=services(fold))
        proximal_dg = -1.0 * fold.folded[0].start - fold.folded[0].length / 100
        assert ev.raw_score == pytest.approx(proximal_dg)

    def test_the_distal_energy_is_reported_but_not_in_the_score(self) -> None:
        fold = StubFold()
        ev = evaluate(svc=services(fold))
        distal_dg = -1.0 * fold.folded[1].start - fold.folded[1].length / 100
        assert ev.raw_score != pytest.approx(distal_dg)
        assert ev.raw_score != pytest.approx((ev.raw_score + distal_dg) / 2)
        assert distal_dg in [dg for _, dg in ev.windows]

    def test_both_windows_travel_for_the_report(self) -> None:
        ev = evaluate()
        assert len(ev.windows) == 2

    def test_measuring_the_distal_window_at_all_does_not_move_the_score(self) -> None:
        """The decisive one: if the distal energy contributed, folding it or not
        would change raw_score. The stub gives the two windows very different
        energies (-0.6 vs -60.6), so any averaging would be unmissable."""
        with_distal = evaluate(rule=StructureWindows(report_distal=True))
        without = evaluate(rule=StructureWindows(report_distal=False))
        assert with_distal.raw_score == pytest.approx(without.raw_score)
        assert len(with_distal.windows) == 2
        assert len(without.windows) == 1


# -- aggregation and availability ---------------------------------------------


class TestAggregation:
    def test_a_second_slot_does_not_dilute_the_score(self) -> None:
        """`min` over gated slots, not a mean: the most structured 5' end is the
        finding. Both slots read the same window here, so this pins the shape --
        a mean would give the same answer only because the inputs are equal, and
        that is what the two-window test above is for."""
        ev = evaluate(ctx=context(slot(role="producer"), slot(role="target")))
        single = evaluate(ctx=context(slot(role="producer")))
        assert ev.raw_score == pytest.approx(single.raw_score)

    def test_it_is_an_objective_not_a_constraint(self) -> None:
        """No published cutoff for these windows, only a monotone relationship."""
        ev = evaluate()
        assert ev.passes
        assert ev.breaches == ()


class TestUnavailable:
    def test_no_folding_engine_is_unavailable_not_zero(self) -> None:
        """0 kcal/mol is a real and maximally GOOD value -- an unstructured 5' end.

        Reporting it for "we could not measure" puts the objective at the top of
        the ranking exactly when it is unknown.
        """
        engineless = Services(
            fold=None,
            kmer=ConstructKmerIndex,
            tables=None,  # type: ignore[arg-type]
            rng=np.random.default_rng(42),
        )
        ev = evaluate(svc=engineless)
        assert is_unavailable(ev)
        assert "no folding engine" in ev.breaches[0].message

    def test_an_unannotated_leader_is_unavailable(self) -> None:
        """The window spans the junction, and unannotated upstream sequence may be
        promoter -- folding it reports a molecule that is never transcribed."""
        ev = evaluate(folded_construct(annotate=False))
        assert is_unavailable(ev)
        assert "no annotated 5'UTR" in ev.breaches[0].message

    def test_a_eukaryotic_context_is_unavailable(self) -> None:
        """brief.md:62's Context column is Bacteria; B11 owns the cap-proximal rule."""
        ev = evaluate(ctx=context(slot(host=HostId.HEK293, modality=Modality.LENTIVIRAL, table=1)))
        assert is_unavailable(ev)
        assert "bacterial expression slot" in ev.breaches[0].message

    def test_a_cds_too_close_to_a_linear_end_is_unavailable(self) -> None:
        """A window clamped to what fits is a DIFFERENT window, not comparable."""
        ev = evaluate(folded_construct(leader="GCGC", circular=False))
        assert is_unavailable(ev)

    def test_the_reason_carries_the_calibration(self) -> None:
        ev = evaluate(folded_construct(annotate=False))
        assert ev.breaches[0].detail["calibration"] == ENGINE_CALIBRATION


# -- metadata -----------------------------------------------------------------


class TestSpecMetadata:
    def test_it_names_the_engine_its_numbers_were_measured_on(self) -> None:
        """registry.check_engine_calibration refuses the run under another engine,
        because the comparison would succeed while meaning nothing (CLAUDE.md 6)."""
        assert StructureWindows.engine_calibration == ENGINE_CALIBRATION
        assert StructureWindows.engine_calibration == "viennarna:rna_turner2004"

    def test_it_is_soft_higher_is_better_and_never_steers(self) -> None:
        assert StructureWindows.enforcement is Enforcement.SOFT
        assert StructureWindows.direction is Direction.HIGHER_IS_BETTER
        assert StructureWindows.steering_weight == 0.0
        assert StructureWindows().lattice_terms(context()) is None

    def test_it_is_weighted_below_b1_because_the_windows_overlap(self) -> None:
        """Kudla's -4..+37 sits inside STR(-30:+30) and both gate to bacteria, so
        equal weights would count one phenomenon twice."""
        from bt5.rules.catalog.b1_five_prime import FivePrimeFolding

        assert StructureWindows.default_weight == 0.5
        assert StructureWindows.default_weight < FivePrimeFolding.default_weight

    def test_folding_is_declared_expensive(self) -> None:
        """cost_class drives null sampling, and this is the slow one."""
        assert StructureWindows.cost_class == "expensive"

    def test_the_evidence_badge_matches_the_briefs_grade(self) -> None:
        assert StructureWindows.evidence is Evidence.EVIDENCE_BACKED

    def test_the_window_constants_are_the_briefs(self) -> None:
        assert (PROXIMAL_UPSTREAM, PROXIMAL_DOWNSTREAM) == (30, 30)
        assert (DISTAL_START + 1, DISTAL_START + DISTAL_LENGTH) == (31, 90)


def test_it_is_registered_under_its_brief_row() -> None:
    assert get("b2_structure_windows") is StructureWindows
    assert StructureWindows.brief_ref == "2.B2"
    assert StructureWindows.id == Path(__file__).stem.removeprefix("test_")


def test_no_preset_weights_it_yet() -> None:
    """Adding 2.B2 to a preset belongs to the score lane, not this one."""
    from bt5.score.presets import PRESETS

    assert not [
        entry for preset in PRESETS for entry in preset.entries if entry.brief_ref == "2.B2"
    ]
