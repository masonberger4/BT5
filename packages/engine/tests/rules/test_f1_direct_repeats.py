"""F1: the repeat rule, and the recA- claim the report must not get wrong.

The scientific content here is Q4a. A recA- strain is standard for LVV and AAV
cloning and suppresses only the RecA-DEPENDENT pathway, which needs >~200-300 bp
of homology. The repeats codon choice creates are 15-100 bp, in the
RecA-INDEPENDENT regime the strain does not touch. A report that implies
otherwise is worse than no report, so several of these tests are about the
message text rather than the finding.
"""

from __future__ import annotations

import numpy as np
import pytest
from bt5.core.context import Modality
from bt5.core.registry import discover, get
from bt5.core.services import Services
from bt5.core.types import Construct, Interval, Segment, SegmentKind, Topology
from bt5.rules.catalog.f1_direct_repeats import DirectRepeats, _spacer, risk_band
from bt5.vector.kmers import ConstructKmerIndex
from conftest import construct, context, slot

discover()


def dna(n: int, seed: int = 3) -> str:
    """Non-repetitive filler. A repetitive pad would BE the finding."""
    rng = np.random.default_rng(seed)
    return "".join("ACGT"[i] for i in rng.integers(0, 4, n))


@pytest.fixture
def svc() -> Services:
    return Services(
        fold=None,
        kmer=ConstructKmerIndex,
        tables=None,  # type: ignore[arg-type]
        rng=np.random.default_rng(1),
    )


def planted(unit: str, spacer_bp: int, *, lead: int = 200) -> Construct:
    """Two copies of `unit`, `spacer_bp` apart, in a circular construct."""
    seq = dna(lead) + unit + dna(spacer_bp, 5) + unit + dna(200, 7)
    return construct(seq[: len(seq) // 2], seq[len(seq) // 2 :])


class TestRiskSurface:
    def test_a_short_distant_repeat_is_low(self) -> None:
        assert risk_band(15, 5000, tandem=False) == "low"

    def test_a_short_tandem_repeat_is_not_low(self) -> None:
        """Slipped-strand mispairing needs no loop, so the length floor does not
        apply to a tandem array."""
        assert risk_band(15, 0, tandem=True) == "moderate"

    def test_proximity_dominates_in_the_reca_independent_regime(self) -> None:
        assert risk_band(30, 10, tandem=False) == "high"
        assert risk_band(30, 5000, tandem=False) == "low"

    def test_a_substantial_repeat_is_never_low_however_far_apart(self) -> None:
        """189 bp of identity between two LTRs is real at any spacing."""
        assert risk_band(189, 3673, tandem=False) == "moderate"

    def test_the_bands_are_not_the_inverted_repeat_bands(self) -> None:
        """F3 exists separately for a reason: these bands are calibrated on the
        DELETION literature and would be a category error on a cruciform."""
        assert risk_band(40, 0, tandem=True) == "high"


class TestSpacer:
    def test_linear_gap(self) -> None:
        assert _spacer(Interval(0, 40), Interval(340, 380), 800, circular=False) == 300

    def test_circular_takes_the_short_way_round(self) -> None:
        """Measuring the long way round turns the most dangerous configuration
        into the safest-looking number in the report."""
        assert _spacer(Interval(20, 60), Interval(700, 740), 800, circular=True) == 80

    def test_touching_copies_are_zero(self) -> None:
        assert _spacer(Interval(0, 40), Interval(40, 80), 800, circular=False) == 0

    def test_the_breach_interval_takes_the_same_arc_the_spacer_measured(
        self, svc: Services
    ) -> None:
        """Two copies bound two arcs, and the finding must describe one of them.

        `_spacer` measures the short way round; the breach interval took the
        forward arc unconditionally, so on this fixture the message read
        "80 bp apart" while the interval covered the other 720 bp. That is not
        only cosmetic: F1 declares `LocalizationPolicy.PAIRED_SEGMENTS`, and
        `solver/repair.py` hands the interval straight to the repair window.
        """
        unit = dna(40, 11)
        seq = dna(20, 3) + unit + dna(640, 5) + unit + dna(60, 7)
        assert len(seq) == 800, "copies land at [20,60) and [700,740)"
        c = Construct(
            seq,
            Topology.CIRCULAR,
            (Segment(Interval(0, len(seq)), SegmentKind.DESIGNABLE_CDS, "cds"),),
        )
        found = DirectRepeats().evaluate(c, context(), svc).breaches
        planted_hit = [b for b in found if "80 bp apart" in b.message]
        assert planted_hit, "the wrap arc is the short one, so the spacer is 80 bp"
        iv = planted_hit[0].interval
        assert iv.end > c.length, "the short arc runs through the origin, so the interval wraps"
        assert iv.length == 2 * 40 + 80, "the interval spans both copies plus the gap it reports"


class TestDetection:
    def test_finds_a_planted_repeat_with_its_geometry(self, svc: Services) -> None:
        breaches = DirectRepeats().evaluate(planted(dna(40, 9), 300), context(), svc).breaches
        assert breaches
        found = breaches[0]
        assert found.detail["length"] == 40.0
        assert found.detail["spacer"] == 300.0

    def test_non_repetitive_sequence_is_clean(self, svc: Services) -> None:
        c = construct(dna(400, 21), dna(400, 22))
        assert not DirectRepeats().evaluate(c, context(), svc).breaches

    def test_a_tandem_array_is_flagged_as_tandem(self, svc: Services) -> None:
        unit = dna(30, 17)
        seq = dna(200) + unit * 4 + dna(200, 19)
        c = construct(seq[:250], seq[250:])
        breaches = DirectRepeats().evaluate(c, context(), svc).breaches
        assert any(b.detail["tandem"] == "yes" for b in breaches)
        assert any(b.detail["risk"] == "high" for b in breaches)

    def test_a_repeat_spanning_the_origin_is_found(self, svc: Services) -> None:
        """The case that justifies evaluating on the assembled circular
        construct rather than the CDS in isolation."""
        unit = dna(40, 11)
        seq = unit[20:] + dna(400) + unit + dna(400, 13) + unit[:20]
        c = Construct(
            seq,
            Topology.CIRCULAR,
            (Segment(Interval(0, len(seq)), SegmentKind.DESIGNABLE_CDS, "cds"),),
        )
        assert DirectRepeats().evaluate(c, context(), svc).breaches

    def test_a_hard_length_repeat_fails_the_evaluation(self, svc: Services) -> None:
        assert not DirectRepeats().evaluate(planted(dna(40, 9), 300), context(), svc).passes

    def test_a_short_repeat_is_reported_without_failing(self, svc: Services) -> None:
        """Info-band findings are worth reporting and are not defects."""
        result = DirectRepeats().evaluate(planted(dna(17, 31), 300), context(), svc)
        assert result.passes


class TestRecAClaim:
    """The report must never tell a user that recA- covers a short repeat."""

    def test_a_short_repeat_says_the_strain_does_not_cover_it(self, svc: Services) -> None:
        breach = DirectRepeats().evaluate(planted(dna(40, 9), 300), context(), svc).breaches[0]
        assert breach.detail["reca_strain_helps"] == "no"
        assert "does NOT cover this" in breach.message
        assert "RecA-independent" in breach.message

    def test_a_long_repeat_says_the_strain_does(self, svc: Services) -> None:
        breach = DirectRepeats().evaluate(planted(dna(250, 41), 300), context(), svc).breaches[0]
        assert breach.detail["reca_strain_helps"] == "yes"
        assert "does NOT cover" not in breach.message

    def test_the_threshold_is_the_reca_dependence_floor(self, svc: Services) -> None:
        short = DirectRepeats().evaluate(planted(dna(150, 43), 300), context(), svc).breaches[0]
        long = DirectRepeats().evaluate(planted(dna(250, 47), 300), context(), svc).breaches[0]
        assert short.detail["reca_strain_helps"] == "no"
        assert long.detail["reca_strain_helps"] == "yes"


class TestExemptRegions:
    """LTRs and ITRs are long perfect direct repeats BY CONSTRUCTION."""

    def whitelisted(self, unit: str, *, both: bool) -> Construct:
        lead, mid = dna(100), dna(100, 29)
        seq = lead + unit + mid + unit + dna(100, 31)
        first = Interval(len(lead), len(lead) + len(unit))
        second_start = len(lead) + len(unit) + len(mid)
        second = Interval(second_start, second_start + len(unit))
        exempt = [Segment(first, SegmentKind.WHITELISTED_REPEAT, "LTR")]
        if both:
            exempt.append(Segment(second, SegmentKind.WHITELISTED_REPEAT, "LTR"))
        covered = {s.interval for s in exempt}
        rest = [
            Segment(iv, SegmentKind.DESIGNABLE_CDS, "cds")
            for iv in (Interval(0, len(lead)),)
            if iv not in covered
        ]
        return Construct(seq, Topology.CIRCULAR, (*exempt, *rest))

    def test_a_pair_wholly_inside_exempt_regions_is_not_reported(self, svc: Services) -> None:
        c = self.whitelisted(dna(60, 23), both=True)
        assert not DirectRepeats().evaluate(c, context(), svc).breaches

    def test_a_pair_with_only_one_copy_exempt_is_still_reported(self, svc: Services) -> None:
        """A designed repeat that happens to match part of an ITR is still the
        design's problem."""
        c = self.whitelisted(dna(60, 23), both=False)
        assert DirectRepeats().evaluate(c, context(), svc).breaches


class TestFixability:
    def test_a_repeat_touching_the_cds_is_fixable(self, svc: Services) -> None:
        unit = dna(40, 9)
        seq = dna(200) + unit + dna(300, 5) + unit + dna(200, 7)
        c = construct(seq[:300], seq[300:])  # the first copy sits in the CDS
        breach = DirectRepeats().evaluate(c, context(), svc).breaches[0]
        assert breach.fixable_by_codon_choice

    def test_a_repeat_wholly_in_the_backbone_is_reported_but_unfixable(self, svc: Services) -> None:
        unit = dna(40, 9)
        seq = "ATGCTGTAA" + dna(200) + unit + dna(300, 5) + unit + dna(200, 7)
        c = construct("ATGCTGTAA", seq[9:])
        breaches = DirectRepeats().evaluate(c, context(), svc).breaches
        assert breaches, "a repeat in the user's own vector is still worth reporting"
        assert not any(b.fixable_by_codon_choice for b in breaches)


class TestContract:
    def test_it_is_hard_but_carries_no_objective_weight(self) -> None:
        assert DirectRepeats.enforcement.is_hard
        assert DirectRepeats.default_weight == 0.0

    def test_it_steers_the_dp_instead(self) -> None:
        """Repeats outrank GC in the only published model trained on real vendor
        outcomes, so the steering term is the highest in the catalog."""
        assert DirectRepeats.steering_weight >= 1.0

    def test_it_is_not_a_lattice_rule(self) -> None:
        """Whether a codon completes a repeat depends on the whole construct,
        not on a bounded suffix."""
        assert DirectRepeats().lattice_terms(None) is None

    def test_it_applies_in_every_context(self) -> None:
        """Every construct passes through a cloning host.

        `slot` is imported at module level, not here. There are four
        `conftest.py` files under packages/engine/tests and none is a package,
        so a deferred `from conftest import ...` resolves against whichever one
        pytest bound last -- this test passed alone and failed in the full suite
        against the VECTOR lane's conftest.
        """
        for modality in Modality:
            assert DirectRepeats().gate(slot(modality=modality))

    def test_absurd_parameters_are_refused(self) -> None:
        with pytest.raises(ValueError, match="distinguishable from chance"):
            DirectRepeats(min_len=4)
        with pytest.raises(ValueError, match="must not be below"):
            DirectRepeats(min_len=30, hard_len=20)

    def test_it_is_registered_under_its_brief_row(self) -> None:
        assert get("f1_direct_repeats") is DirectRepeats
        assert DirectRepeats.brief_ref == "2.F1"
