"""B1: the highest-weight objective, and the three ways it declines to answer.

Most of this file is about gate G6. The window Kudla measured spans the UTR/CDS
junction, so B1 is the one rule in the catalog that can be handed everything it
needs except one annotation and still have to say "I cannot compute this". The
tempting failure -- fold the backbone that happens to precede the CDS and call it
the leader -- produces a number that looks entirely reasonable and describes a
molecule that is never transcribed.

Measured, on 300 random 4 nt leaders spliced onto random 37 nt windows: the
leader changes the number 53% of the time and by >= 1 kcal/mol a third of the
time, up to 7 kcal/mol, and never in the positive direction. So you cannot tell
from the number whether the leader was included -- which is exactly why the rule
has to KNOW, and why an unannotated 5'UTR is a refusal rather than a guess.

No dG is asserted byte-exactly anywhere here; the assertions are relations.
Energy parameters determine every dG, so a pinned literal would turn a ViennaRNA
upgrade into a mystery failure rather than a reviewable one.

`conftest` is imported at module level; see the note in test_f1_direct_repeats.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from bt5.core.context import (
    BiosecurityVerdict,
    ContextSlot,
    DesignContext,
    HostId,
    Modality,
)
from bt5.core.registry import all_specs, discover, get
from bt5.core.services import Services
from bt5.core.spec import Direction, Enforcement
from bt5.core.types import (
    Construct,
    Feature,
    Interval,
    Segment,
    SegmentKind,
    Strand,
    Topology,
)
from bt5.rules.catalog.b1_five_prime import (
    ENGINE_CALIBRATION,
    KUDLA_DOWNSTREAM,
    KUDLA_UPSTREAM,
    FivePrimeFolding,
    leader_of,
)
from bt5.structure.vienna import load_fold_engine
from bt5.vector.kmers import ConstructKmerIndex

discover()

LEAD = 200


def dna(n: int, seed: int = 3) -> str:
    rng = np.random.default_rng(seed)
    return "".join("ACGT"[i] for i in rng.integers(0, 4, n))


@pytest.fixture
def svc() -> Services:
    """A REAL folding engine. ViennaRNA is a declared, pinned dependency, so a
    skip here would hide a broken install rather than tolerate one."""
    return Services(
        fold=load_fold_engine(),
        kmer=ConstructKmerIndex,
        tables=None,  # type: ignore[arg-type]
        rng=np.random.default_rng(1),
    )


@pytest.fixture
def no_fold() -> Services:
    return Services(
        fold=None,
        kmer=ConstructKmerIndex,
        tables=None,  # type: ignore[arg-type]
        rng=np.random.default_rng(1),
    )


def build(
    utr: str,
    cds: str,
    *,
    annotate: str | None = "upstream",
    circular: bool = True,
    lead: int = LEAD,
) -> Construct:
    """`lead` backbone, then `utr`, then `cds`, then trailing backbone.

    `annotate` places the 5'UTR feature: "upstream" of the CDS (the forward
    cassette's leader), "downstream" of it (a reverse cassette's leader, which
    sits at HIGHER construct coordinates), or None for no annotation at all.
    """
    head = dna(lead, 41) if lead else ""
    tail = dna(200, 43)
    seq = head + utr + cds + tail
    utr_start, cds_start = len(head), len(head) + len(utr)
    cds_end = cds_start + len(cds)

    features: tuple[Feature, ...] = ()
    if annotate == "upstream":
        features = (Feature(Interval(utr_start, cds_start), "5'UTR"),)
    elif annotate == "downstream":
        features = (Feature(Interval(cds_end, cds_end + len(utr)), "5'UTR"),)

    segments = [Segment(Interval(cds_start, cds_end), SegmentKind.DESIGNABLE_CDS, "cds")]
    if cds_start:
        segments.append(Segment(Interval(0, cds_start), SegmentKind.BACKBONE, "v"))
    if cds_end < len(seq):
        segments.append(Segment(Interval(cds_end, len(seq)), SegmentKind.BACKBONE, "v2"))
    return Construct(
        sequence=seq,
        topology=Topology.CIRCULAR if circular else Topology.LINEAR,
        segments=tuple(segments),
        features=features,
    )


def context(
    modality: Modality = Modality.BACTERIAL_EXPRESSION,
    *,
    orientation: Strand = 1,
    host: HostId = HostId.E_COLI_K12,
    table: int = 11,
) -> DesignContext:
    return DesignContext(
        slots=(ContextSlot("producer", host, modality, table),),
        cassette_orientation=orientation,
        seed=42,
        screen=BiosecurityVerdict("not_run"),
    )


#: A leader/CDS pair measured to shift the window by 4.0 kcal/mol, so the
#: difference is far outside any plausible parameter jitter.
STRONG_LEADER = dna(4, 903)
STRONG_CDS = "ATG" + dna(34, 0)


class TestG6TheLeaderIsPartOfTheMeasurement:
    """The window spans the UTR/CDS junction. That is the whole rule."""

    def test_the_leader_changes_the_number(self, svc: Services) -> None:
        """Folded with the real leader spliced on, against the same 37 coding
        bases alone. If these ever agree, B1 has stopped measuring Kudla's
        window and the 5'UTR requirement below is protecting nothing."""
        with_leader = (
            FivePrimeFolding().evaluate(build(STRONG_LEADER, STRONG_CDS), context(), svc).raw_score
        )
        assert svc.fold is not None
        cds_alone = svc.fold.mfe_window(STRONG_CDS, Interval(0, KUDLA_DOWNSTREAM)).dg_kcal_mol
        assert with_leader < cds_alone - 2.0, (
            f"leader-inclusive {with_leader:.2f} vs CDS-alone {cds_alone:.2f}: the "
            f"window is not picking up the leader"
        )

    def test_adding_a_leader_never_makes_the_window_less_structured(self, svc: Services) -> None:
        """Four more bases can only add pairing options, never remove them, so
        the leader-inclusive dG is bounded above by the CDS-alone one. A
        violation means the window is landing somewhere other than intended."""
        assert svc.fold is not None
        for seed in range(6):
            cds = "ATG" + dna(34, seed)
            got = FivePrimeFolding().evaluate(build(dna(4, 900 + seed), cds), context(), svc)
            alone = svc.fold.mfe_window(cds, Interval(0, KUDLA_DOWNSTREAM)).dg_kcal_mol
            assert got.raw_score <= alone + 1e-6

    def test_the_window_is_the_published_one(self, svc: Services) -> None:
        got = FivePrimeFolding().evaluate(build(dna(12, 5), STRONG_CDS), context(), svc)
        window, _ = got.windows[0]
        assert window.length == KUDLA_UPSTREAM + KUDLA_DOWNSTREAM == 41
        assert got.n_evaluated == 41


class TestG6DegradesInsteadOfGuessing:
    """Three refusals, each a number BT5 declines to invent."""

    def unavailable(self, result) -> str:
        assert math.isnan(result.raw_score), (
            "0.0 kcal/mol is a real value for this quantity -- an unstructured 5' "
            "end, the best possible score -- so it must never stand in for 'unknown'"
        )
        assert result.breaches
        return str(result.breaches[0].detail["unavailable_reason"])

    def test_no_annotated_utr(self, svc: Services) -> None:
        """The sequence upstream of a CDS is often promoter, which is not
        transcribed. Folding it would report a molecule that never exists."""
        got = FivePrimeFolding().evaluate(
            build(STRONG_LEADER, STRONG_CDS, annotate=None), context(), svc
        )
        reason = self.unavailable(got)
        assert "no annotated 5'UTR" in reason
        assert "Annotate the 5'UTR" in reason, "a refusal is only useful with a way out"

    def test_no_folding_engine(self, no_fold: Services) -> None:
        assert "no folding engine" in self.unavailable(
            FivePrimeFolding().evaluate(build(STRONG_LEADER, STRONG_CDS), context(), no_fold)
        )

    def test_the_window_does_not_fit(self, svc: Services) -> None:
        """A linear construct whose CDS starts at position 0 has no leader, and a
        window clamped to what fits is a different window."""
        got = FivePrimeFolding().evaluate(
            build("", STRONG_CDS, annotate=None, circular=False, lead=0), context(), svc
        )
        assert "too close to the end" in self.unavailable(got)

    def test_a_non_bacterial_context(self, svc: Services) -> None:
        """Eukaryotic scanning initiation makes cap-proximal structure the
        analogous term (B11's), with its own position-dependent ladder."""
        got = FivePrimeFolding().evaluate(
            build(STRONG_LEADER, STRONG_CDS), context(Modality.LENTIVIRAL), svc
        )
        assert "bacterial expression" in self.unavailable(got)

    def test_a_refusal_is_not_the_codon_choices_fault(self, svc: Services) -> None:
        got = FivePrimeFolding().evaluate(
            build(STRONG_LEADER, STRONG_CDS, annotate=None), context(), svc
        )
        assert not got.breaches[0].fixable_by_codon_choice

    def test_the_e_coli_slot_of_a_viral_job_is_not_an_expression_event(self, svc: Services) -> None:
        """The gate reads modality, not host. A lentiviral plasmid is propagated
        in E. coli and never expressed there, so scoring its 5' structure would
        be scoring a translation event that does not happen."""
        assert not FivePrimeFolding().gate(
            ContextSlot("propagation", HostId.E_COLI_K12, Modality.LENTIVIRAL, 11)
        )
        assert FivePrimeFolding().gate(
            ContextSlot("producer", HostId.E_COLI_K12, Modality.BACTERIAL_EXPRESSION, 11)
        )


class TestStrand:
    """For a reverse-oriented cassette the 5' end is at HIGHER coordinates."""

    def test_the_leader_of_a_forward_window_is_its_head(self) -> None:
        assert leader_of(Interval(100, 141, 1), 4) == Interval(100, 104, 1)

    def test_the_leader_of_a_reverse_window_is_its_tail(self) -> None:
        assert leader_of(Interval(100, 141, -1), 4) == Interval(137, 141, -1)

    def test_a_reverse_cassette_reads_the_utr_at_the_other_end(self, svc: Services) -> None:
        got = FivePrimeFolding().evaluate(
            build(STRONG_LEADER, STRONG_CDS, annotate="downstream"),
            context(orientation=-1),
            svc,
        )
        assert not math.isnan(got.raw_score)
        window, _ = got.windows[0]
        assert window.strand == -1

    def test_a_reverse_cassette_is_not_satisfied_by_the_forward_utr(self, svc: Services) -> None:
        """The test that catches a hard-coded strand 1: annotate the leader where
        a FORWARD cassette would have it, then ask about a reverse one. The
        annotation is real and in the wrong place, and the rule must decline."""
        got = FivePrimeFolding().evaluate(
            build(STRONG_LEADER, STRONG_CDS, annotate="upstream"),
            context(orientation=-1),
            svc,
        )
        assert math.isnan(got.raw_score)
        assert "no annotated 5'UTR" in str(got.breaches[0].detail["unavailable_reason"])

    def test_the_two_orientations_fold_different_sequence(self, svc: Services) -> None:
        forward = FivePrimeFolding().evaluate(
            build(STRONG_LEADER, STRONG_CDS, annotate="upstream"), context(), svc
        )
        reverse = FivePrimeFolding().evaluate(
            build(STRONG_LEADER, STRONG_CDS, annotate="downstream"),
            context(orientation=-1),
            svc,
        )
        assert forward.windows[0][0] != reverse.windows[0][0]


class TestContract:
    def test_it_is_soft_and_carries_the_highest_weight_in_the_catalog(self) -> None:
        assert FivePrimeFolding.enforcement is Enforcement.SOFT
        assert FivePrimeFolding.default_weight == 1.0
        others = [s.default_weight for s in all_specs() if s.id != FivePrimeFolding.id]
        assert FivePrimeFolding.default_weight > max(others), (
            "B1 is the only objective in BT5 justified by a measured effect size "
            "(r = 0.66, 44% of variance) rather than a feature ranking"
        )

    def test_higher_is_better_and_not_a_band(self) -> None:
        """Less negative is less structure is better initiation, monotone across
        the whole 250-fold range Kudla observed -- unlike CAI or GC, which have
        an optimum in the middle and would be misrepresented by a monotone term."""
        assert FivePrimeFolding.direction is Direction.HIGHER_IS_BETTER
        assert FivePrimeFolding.band is None

    def test_it_declares_the_engine_its_thresholds_are_measured_on(self) -> None:
        """Unlike e5, whose Tm is not a folding energy, B1's quantity really is a
        ViennaRNA kcal/mol -- so a mismatched engine must refuse the run rather
        than silently compare numbers taken with two different rulers."""
        assert FivePrimeFolding.engine_calibration == ENGINE_CALIBRATION
        assert ENGINE_CALIBRATION == "viennarna:rna_turner2004"

    def test_it_does_not_steer_the_dp(self) -> None:
        """Tier A decides from a bounded suffix; a fold is not that. A cheap
        proxy would steer on a quantity this rule does not score."""
        assert FivePrimeFolding.steering_weight == 0.0

    def test_it_is_not_a_lattice_rule(self) -> None:
        assert FivePrimeFolding().lattice_terms(None) is None

    def test_a_cds_only_window_is_refused(self) -> None:
        with pytest.raises(ValueError, match="fold the CDS alone"):
            FivePrimeFolding(upstream=0)
        with pytest.raises(ValueError, match="reach into the CDS"):
            FivePrimeFolding(downstream=0)

    def test_the_shipped_bacterial_preset_now_binds_it(self) -> None:
        from bt5.score.presets import BACTERIAL, resolve

        resolved = resolve(BACTERIAL)
        assert "b1_five_prime" in resolved.weights
        assert "2.B1" not in resolved.unimplemented

    def test_it_is_registered_under_its_brief_row(self) -> None:
        assert get("b1_five_prime") is FivePrimeFolding
        assert FivePrimeFolding.brief_ref == "2.B1"
