"""F2: the repeats F1 cannot see, and the boundary between the two rules.

The tests that matter most are the ones with a known ground truth: a repeat
planted at a chosen length and a chosen number of mismatches, where the rule
must report exactly that and not an artefact of its own stopping rule.
"""

from __future__ import annotations

import numpy as np
import pytest
from bt5.core.context import Modality
from bt5.core.registry import discover, get
from bt5.core.services import Services
from bt5.core.spec import Enforcement
from bt5.core.types import Construct, Interval, Segment, SegmentKind, Topology
from bt5.rules.catalog.f1_direct_repeats import DirectRepeats
from bt5.rules.catalog.f2_near_perfect_repeats import (
    MIN_LENGTH_BP,
    SPACER_DECAY_BP,
    NearPerfectRepeats,
    risk_score,
)
from bt5.vector.kmers import ConstructKmerIndex
from conftest import construct, context, slot

discover()


def dna(n: int, seed: int = 3) -> str:
    rng = np.random.default_rng(seed)
    return "".join("ACGT"[i] for i in rng.integers(0, 4, n))


def mutate(unit: str, positions: tuple[int, ...]) -> str:
    """A copy of `unit` differing at exactly `positions`."""
    out = list(unit)
    for p in positions:
        out[p] = "ACGT"[("ACGT".index(out[p]) + 1) % 4]
    return "".join(out)


@pytest.fixture
def svc() -> Services:
    return Services(
        fold=None,
        kmer=ConstructKmerIndex,
        tables=None,  # type: ignore[arg-type]
        rng=np.random.default_rng(1),
    )


def planted(unit: str, copy: str, spacer_bp: int = 300) -> Construct:
    seq = dna(200) + unit + dna(spacer_bp, 5) + copy + dna(200, 7)
    return construct(seq[:400], seq[400:])


def hits(rule: NearPerfectRepeats, c: Construct, svc: Services, ctx=None):
    return rule.evaluate(c, ctx or context(), svc).breaches


class TestGroundTruth:
    """A planted repeat must be reported at the length and identity it has."""

    def test_reports_the_repeat_exactly_as_constructed(self, svc: Services) -> None:
        unit = dna(60, 9)
        found = hits(NearPerfectRepeats(), planted(unit, mutate(unit, (15, 35, 52))), svc)
        assert len(found) == 1
        assert found[0].detail["length"] == 60.0
        assert found[0].detail["mismatches"] == 3.0
        assert found[0].detail["identity"] == 0.95

    def test_identity_is_not_pinned_to_the_threshold(self, svc: Services) -> None:
        """The trap this rule was rewritten to avoid. Extending while merely
        above the floor bleeds into flanking sequence until enough mismatches
        accumulate to hit it, so EVERY finding comes back at ~90% and several
        bases too long -- an artefact of the stopping rule, not a measurement.
        """
        unit = dna(80, 13)
        one = hits(NearPerfectRepeats(), planted(unit, mutate(unit, (40,))), svc)[0]
        many = hits(NearPerfectRepeats(), planted(unit, mutate(unit, (10, 30, 50, 70))), svc)[0]
        assert one.detail["identity"] > many.detail["identity"]
        assert one.detail["identity"] > 0.98, "one mismatch in 80 bp is ~99%, not 90%"

    def test_one_physical_repeat_is_one_finding(self, svc: Services) -> None:
        """Adjacent seeds inside one repeat each extend to nearly the same
        extent, so a key on exact coordinates deduplicates nothing -- this
        reported the same 60 bp repeat three times, at 197, 198 and 200."""
        unit = dna(60, 9)
        assert len(hits(NearPerfectRepeats(), planted(unit, mutate(unit, (15, 35, 52))), svc)) == 1

    def test_a_pair_below_the_identity_floor_is_not_reported(self, svc: Services) -> None:
        unit = dna(60, 9)
        far = mutate(unit, tuple(range(2, 58, 4)))  # ~23% mismatched
        assert not hits(NearPerfectRepeats(), planted(unit, far), svc)

    def test_a_pair_below_the_length_floor_is_not_reported(self, svc: Services) -> None:
        unit = dna(24, 33)
        assert not hits(NearPerfectRepeats(), planted(unit, mutate(unit, (12,))), svc)

    def test_clean_sequence_is_clean(self, svc: Services) -> None:
        assert not hits(NearPerfectRepeats(), construct(dna(400, 21), dna(400, 22)), svc)


class TestBoundaryWithF1:
    """F1 owns exact repeats; F2 owns the ones exact matching cannot see."""

    def test_an_exact_pair_belongs_to_f1_alone(self, svc: Services) -> None:
        unit = dna(60, 9)
        c = planted(unit, unit)
        assert DirectRepeats().evaluate(c, context(), svc).breaches
        assert not hits(NearPerfectRepeats(), c, svc), "reporting it twice doubles the panel"

    def test_a_near_perfect_pair_is_invisible_to_f1(self, svc: Services) -> None:
        """The case that justifies the rule existing at all."""
        unit = dna(60, 9)
        c = planted(unit, mutate(unit, (10, 20, 30, 40, 50)))
        exact = DirectRepeats().evaluate(c, context(), svc).breaches
        assert not any(b.detail["length"] >= 40 for b in exact), (
            "an exact scan sees only the short fragments between mismatches"
        )
        assert hits(NearPerfectRepeats(), c, svc)


class TestRiskScore:
    def test_monotone_in_length(self) -> None:
        assert risk_score(80, 500) > risk_score(40, 500)

    def test_decays_with_spacer(self) -> None:
        assert risk_score(60, 100) > risk_score(60, 5000)

    def test_the_decay_length_is_the_brief_s_high_risk_boundary(self) -> None:
        """The brief puts high risk under 3 kb, so at 3 kb the factor is 1/e."""
        assert risk_score(MIN_LENGTH_BP, SPACER_DECAY_BP) == pytest.approx(np.exp(-1), rel=1e-6)

    def test_a_minimum_length_adjacent_pair_scores_one(self) -> None:
        assert risk_score(MIN_LENGTH_BP, 0) == pytest.approx(1.0)

    def test_the_score_travels_with_the_finding(self, svc: Services) -> None:
        unit = dna(60, 9)
        breach = hits(NearPerfectRepeats(), planted(unit, mutate(unit, (15, 35, 52))), svc)[0]
        assert breach.detail["risk_score"] > 0
        # `detail` carries the rounded value for display; `magnitude` is the
        # full float the weighted sum sees.
        assert breach.magnitude == pytest.approx(breach.detail["risk_score"], abs=1e-4)

    def test_the_unit_says_it_is_not_a_rate(self) -> None:
        """The literature supports the SHAPE of the relationship, not a
        calibrated frequency, and BT5 never prints a predicted rate."""
        assert "not a rate" in NearPerfectRepeats.unit


class TestGeometry:
    def test_a_near_perfect_repeat_spanning_the_origin_is_found(self, svc: Services) -> None:
        unit = dna(60, 11)
        copy = mutate(unit, (20, 45))
        seq = unit[30:] + dna(400) + copy + dna(400, 13) + unit[:30]
        c = Construct(
            seq,
            Topology.CIRCULAR,
            (Segment(Interval(0, len(seq)), SegmentKind.DESIGNABLE_CDS, "cds"),),
        )
        assert hits(NearPerfectRepeats(), c, svc), "a repeat across position 0 is still a repeat"

    def test_closer_copies_score_higher(self, svc: Services) -> None:
        unit = dna(60, 9)
        copy = mutate(unit, (15, 35, 52))
        near = hits(NearPerfectRepeats(), planted(unit, copy, spacer_bp=100), svc)[0]
        far = hits(NearPerfectRepeats(), planted(unit, copy, spacer_bp=2000), svc)[0]
        assert near.magnitude > far.magnitude


class TestExemptRegions:
    def test_a_pair_inside_whitelisted_repeats_is_not_reported(self, svc: Services) -> None:
        unit = dna(60, 23)
        copy = mutate(unit, (10, 30, 50))
        lead, mid = dna(100), dna(100, 29)
        seq = lead + unit + mid + copy + dna(100, 31)
        first = Interval(len(lead), len(lead) + len(unit))
        second_start = len(lead) + len(unit) + len(mid)
        c = Construct(
            seq,
            Topology.CIRCULAR,
            (
                Segment(first, SegmentKind.WHITELISTED_REPEAT, "LTR"),
                Segment(
                    Interval(second_start, second_start + len(copy)),
                    SegmentKind.WHITELISTED_REPEAT,
                    "LTR",
                ),
                Segment(Interval(0, len(lead)), SegmentKind.DESIGNABLE_CDS, "cds"),
            ),
        )
        assert not hits(NearPerfectRepeats(), c, svc)


class TestContract:
    def test_it_is_soft_so_a_preset_can_weight_it(self) -> None:
        """The packaging presets weight 2.F2 at 1.0, and `resolve()` refuses to
        weight anything that is not SOFT."""
        assert NearPerfectRepeats.enforcement is Enforcement.SOFT
        assert NearPerfectRepeats.default_weight > 0.0
        assert NearPerfectRepeats.weight_provenance.strip()

    def test_the_shipped_presets_now_bind_this_rule(self) -> None:
        from bt5.score.presets import LENTIVIRAL, resolve

        resolved = resolve(LENTIVIRAL)
        assert "f2_near_perfect_repeats" in resolved.weights
        assert "2.F2" not in resolved.unimplemented

    def test_it_is_weighted_below_the_exact_repeat_rule(self) -> None:
        """The threshold here is a convention, not a measured constant, and the
        detector has a stated blind spot. Level weighting would give a
        conventional cut the authority of a measured one."""
        assert NearPerfectRepeats.steering_weight < DirectRepeats.steering_weight

    def test_it_is_not_a_lattice_rule(self) -> None:
        assert NearPerfectRepeats().lattice_terms(None) is None

    def test_it_applies_in_every_context(self) -> None:
        for modality in Modality:
            assert NearPerfectRepeats().gate(slot(modality=modality))

    def test_absurd_parameters_are_refused(self) -> None:
        with pytest.raises(ValueError, match="outside"):
            NearPerfectRepeats(min_identity=0.1)
        with pytest.raises(ValueError, match="occur by chance"):
            NearPerfectRepeats(min_length=8)
        with pytest.raises(ValueError, match="match everywhere"):
            NearPerfectRepeats(seed=4)
        with pytest.raises(ValueError, match="exceeds min_length"):
            NearPerfectRepeats(min_length=20, seed=30)

    def test_it_is_registered_under_its_brief_row(self) -> None:
        assert get("f2_near_perfect_repeats") is NearPerfectRepeats
        assert NearPerfectRepeats.brief_ref == "2.F2"
