"""E8: multiplicity, which is the statement no other repeat rule can make.

Four rules in this catalog now look at repeated sequence, so nearly every test
here is about a boundary. F1 owns pairs on the plasmid at 15 bp and up, E5 owns
pairs in the ordered fragment at 12 bp and up, E7 owns tandem tracts, and this
owns "how many independent places could an oligo land?". A design where those
overlap is a design that reports one repeat four times.

`conftest` is imported at module level; see the note in test_f1_direct_repeats.
"""

from __future__ import annotations

import numpy as np
import pytest
from bt5.core.context import Modality
from bt5.core.registry import discover, get
from bt5.core.services import Services
from bt5.core.spec import Enforcement
from bt5.core.types import Construct, Interval, Segment, SegmentKind, Topology
from bt5.rules.catalog.e5_synthesis_repeats import SynthesisRepeats
from bt5.rules.catalog.e6_repeat_density import RepeatDensity
from bt5.rules.catalog.e7_short_tandem_repeats import ShortTandemRepeats
from bt5.rules.catalog.e8_kmer_uniqueness import (
    KMER_BP,
    MIN_MULTIPLICITY,
    KmerUniqueness,
    landing_sites,
)
from bt5.rules.catalog.f1_direct_repeats import DirectRepeats
from bt5.vector.kmers import ConstructKmerIndex
from conftest import construct, context, slot

discover()


def dna(n: int, seed: int = 3) -> str:
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


def hits(c: Construct, svc: Services, **kw: object):
    return KmerUniqueness(**kw).evaluate(c, context(), svc).breaches  # type: ignore[arg-type]


UNIT = dna(30, 51)


def at_sites(n: int) -> Construct:
    """One 30 bp unit placed at `n` well-separated sites in the CDS."""
    cds = "".join(UNIT + dna(400, 41 + i) for i in range(n))
    return construct(cds, dna(1500, 45))


class TestLandingSites:
    """Occurrences closer than k overlap: they are one place, not several."""

    def test_a_chain_of_overlapping_occurrences_is_one_site(self) -> None:
        """Linkage is to the previous OCCURRENCE, not the site that opened the
        cluster. Comparing against the cluster start reopens a site every k
        bases and turned (CAG)x20 into five sites rather than one stretch."""
        tandem = list(range(700, 750, 3))
        assert landing_sites(tandem, KMER_BP, 3000, circular=True) == (700,)

    def test_well_separated_occurrences_are_separate_sites(self) -> None:
        assert landing_sites([400, 830, 1260], KMER_BP, 3000, circular=True) == (400, 830, 1260)

    def test_occurrences_meeting_across_the_origin_are_one_site(self) -> None:
        assert landing_sites([2995, 2], KMER_BP, 3000, circular=True) == (2,)


class TestBoundariesWithTheRestOfTheFamily:
    @pytest.mark.parametrize("tract", ["CAG" * 20, "CAGGCT" * 25], ids=["(CAG)x20", "(CAGGCT)x25"])
    def test_a_tandem_array_is_e7s_alone(self, tract: str, svc: Services) -> None:
        """A tandem array is one contiguous stretch and one landing site, so it
        does not reach this rule's multiplicity floor -- which is what keeps the
        panel from carrying the same tract under two names."""
        c = construct(dna(700, 41) + tract + dna(700, 43), dna(1500, 45))
        assert ShortTandemRepeats().evaluate(c, context(), svc).breaches
        assert not hits(c, svc)

    def test_a_two_site_repeat_is_the_pair_rules_alone(self, svc: Services) -> None:
        c = at_sites(2)
        assert DirectRepeats().evaluate(c, context(), svc).breaches
        assert not hits(c, svc)

    def test_a_three_site_repeat_is_this_rules(self, svc: Services) -> None:
        """A pair rule can only ever say "these two are the same", three times
        over. "One sequence, three landing sites" is the statement that predicts
        the misassembly."""
        found = hits(at_sites(3), svc)
        assert len(found) == 1
        assert found[0].detail["multiplicity"] == 3.0

    def test_a_backbone_collision_is_visible_here_and_nowhere_else(self, svc: Services) -> None:
        """E5 and E6 see only the ordered fragment, so a 12-mer the insert shares
        with the user's own vector is invisible to them. It is exactly the
        collision that misassembles a Gibson junction."""
        shared = UNIT[:KMER_BP]
        c = construct(
            dna(400, 41) + shared + dna(300, 43),
            dna(400, 45) + shared + dna(400, 47) + shared + dna(300, 49),
        )
        assert not SynthesisRepeats().evaluate(c, context(), svc).breaches
        assert not RepeatDensity().evaluate(c, context(), svc).breaches
        assert hits(c, svc)


class TestOnePhysicalRepeatIsOneFinding:
    def test_a_thirty_bp_repeat_is_not_nineteen_findings(self, svc: Services) -> None:
        """A 30 bp region holds 19 distinct 12-mers, each at three sites. One
        breach per k-mer reported one physical repeat nineteen times -- the same
        failure E1 hit with homopolymers and E6 with overlapping windows."""
        assert len(hits(at_sites(3), svc)) == 1

    def test_the_finding_spans_the_repeated_region_not_one_kmer(self, svc: Services) -> None:
        breach = hits(at_sites(3), svc)[0]
        assert breach.detail["region_bp"] >= 30.0
        assert breach.interval.length >= 30

    def test_every_site_is_named(self, svc: Services) -> None:
        breach = hits(at_sites(4), svc)[0]
        assert breach.detail["multiplicity"] == 4.0
        assert "4 independent sites" in breach.message


class TestTheScalar:
    """Uniqueness by position -- deliberately a different definition from the
    breaches, because a tandem array IS non-unique even though E7 reports it."""

    def test_a_clean_construct_scores_zero(self, svc: Services) -> None:
        result = KmerUniqueness().evaluate(construct(dna(1500, 41), dna(1500, 43)), context(), svc)
        assert result.raw_score == 0.0
        assert result.passes

    def test_a_tandem_array_still_moves_the_scalar(self, svc: Services) -> None:
        """It raises no breach here, but it is not unique and the objective the
        DP steers on must say so."""
        c = construct(dna(700, 41) + "CAG" * 20 + dna(700, 43), dna(1500, 45))
        assert KmerUniqueness().evaluate(c, context(), svc).raw_score > 0.0
        assert not hits(c, svc)

    def test_the_scalar_rises_with_repetitiveness(self, svc: Services) -> None:
        def score(c: Construct) -> float:
            return KmerUniqueness().evaluate(c, context(), svc).raw_score

        assert score(at_sites(4)) > score(at_sites(2)) > 0.0

    def test_it_is_a_bounded_fraction(self, svc: Services) -> None:
        assert 0.0 <= KmerUniqueness().evaluate(at_sites(4), context(), svc).raw_score <= 1.0


class TestGeometry:
    def test_a_kmer_spanning_the_origin_is_counted(self, svc: Services) -> None:
        """The construct is circular, so a 12-mer across position 0 is a real
        12-mer and a real landing site."""
        unit = dna(20, 11)
        seq = unit[10:] + dna(400, 13) + unit + dna(400, 15) + unit + dna(400, 17) + unit[:10]
        c = Construct(
            seq,
            Topology.CIRCULAR,
            (Segment(Interval(0, len(seq)), SegmentKind.DESIGNABLE_CDS, "cds"),),
        )
        found = hits(c, svc)
        assert found
        assert max(float(b.detail["multiplicity"]) for b in found) >= 3


class TestExemptRegions:
    """LTRs make every k-mer in them recur, and no codon can change that."""

    def lentiviral(self, *, insert_collides: bool) -> Construct:
        ltr = dna(200, 81)
        insert = (ltr[:KMER_BP] if insert_collides else dna(KMER_BP, 83)) + dna(300, 85)
        seq = ltr + insert + dna(400, 87) + ltr
        second = len(ltr) + len(insert) + 400
        return Construct(
            seq,
            Topology.CIRCULAR,
            (
                Segment(Interval(0, len(ltr)), SegmentKind.WHITELISTED_REPEAT, "5' LTR"),
                Segment(
                    Interval(second, second + len(ltr)),
                    SegmentKind.WHITELISTED_REPEAT,
                    "3' LTR",
                ),
                Segment(
                    Interval(len(ltr), len(ltr) + len(insert)),
                    SegmentKind.DESIGNABLE_CDS,
                    "cds",
                ),
            ),
        )

    def test_the_ltrs_alone_raise_no_finding(self, svc: Services) -> None:
        assert not hits(self.lentiviral(insert_collides=False), svc)

    def test_but_an_insert_colliding_with_them_does(self, svc: Services) -> None:
        """The LTRs are counted when computing multiplicity even though they are
        not scored -- dropping them would hide exactly the collision the user
        can fix."""
        found = hits(self.lentiviral(insert_collides=True), svc)
        assert found
        assert any(b.fixable_by_codon_choice for b in found)

    def test_whitelisted_sequence_does_not_swamp_the_scalar(self, svc: Services) -> None:
        """Two 200 bp LTRs would otherwise put 400 non-unique positions into the
        score, dominating it with something no codon can move."""
        result = KmerUniqueness().evaluate(self.lentiviral(insert_collides=False), context(), svc)
        assert result.raw_score == 0.0


class TestContract:
    def test_it_is_soft_so_it_can_be_weighted(self) -> None:
        assert KmerUniqueness.enforcement is Enforcement.SOFT
        assert KmerUniqueness.default_weight > 0.0
        assert KmerUniqueness.weight_provenance.strip()

    def test_it_is_weighted_below_the_density_rule_it_overlaps(self) -> None:
        """A fragment repetitive at 9-mers is usually repetitive at 12-mers;
        weighting both at full strength counts one thing twice."""
        assert KmerUniqueness.default_weight < RepeatDensity.default_weight

    def test_it_does_not_stack_another_full_strength_repeat_steering_term(self) -> None:
        assert KmerUniqueness.steering_weight < DirectRepeats.steering_weight

    def test_it_is_not_a_lattice_rule(self) -> None:
        assert KmerUniqueness().lattice_terms(None) is None

    def test_it_applies_in_every_context(self) -> None:
        for modality in Modality:
            assert KmerUniqueness().gate(slot(modality=modality))

    def test_absurd_parameters_are_refused(self) -> None:
        with pytest.raises(ValueError, match="recurs by chance too often"):
            KmerUniqueness(k=10)
        with pytest.raises(ValueError, match="what unique MEANS"):
            KmerUniqueness(min_multiplicity=1)

    def test_the_multiplicity_floor_is_above_the_pairwise_case(self) -> None:
        """At 2 it would re-report every F1 and E5 pair under a third name."""
        assert MIN_MULTIPLICITY >= 3

    def test_it_is_registered_under_its_brief_row(self) -> None:
        assert get("e8_kmer_uniqueness") is KmerUniqueness
        assert KmerUniqueness.brief_ref == "2.E8"
