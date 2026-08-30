"""E7: the tract no pair-based rule can see.

The first test class is the one that justifies the rule existing. A pair scan
asks "are these two spans identical?", and a tandem array answers "yes, at every
offset" -- which the pair extension resolves by clamping to the period and then
dropping the match as too short. So a 90 bp array passes the entire repeat
family. That is a measurement, and it is pinned here so it stays one.

`conftest` is imported at module level; see the note in test_f1_direct_repeats.
"""

from __future__ import annotations

import numpy as np
import pytest
from bt5.core.context import Modality
from bt5.core.registry import discover, get
from bt5.core.services import Services
from bt5.core.spec import Enforcement
from bt5.core.types import Construct
from bt5.rules.catalog.e1_homopolymers import Homopolymers
from bt5.rules.catalog.e5_synthesis_repeats import SynthesisRepeats
from bt5.rules.catalog.e7_short_tandem_repeats import (
    HARD_TRACT_BP,
    WARN_TRACT_BP,
    ShortTandemRepeats,
    tandem_tracts,
)
from bt5.rules.catalog.f1_direct_repeats import DirectRepeats
from bt5.rules.catalog.f2_near_perfect_repeats import NearPerfectRepeats
from bt5.rules.vendors import VendorSelection
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


def with_tract(tract: str) -> Construct:
    return construct(dna(200, 41) + tract + dna(200, 43), dna(300, 45))


def hits(c: Construct, svc: Services, **kw: object):
    return ShortTandemRepeats(**kw).evaluate(c, context(), svc).breaches  # type: ignore[arg-type]


class TestTheGapThisRuleFills:
    """Every other repeat rule is blind to a tandem array."""

    @pytest.mark.parametrize(
        "tract",
        ["CAG" * 20, "AT" * 30, "CAGGCT" * 15],
        ids=["(CAG)x20", "(AT)x30", "(CAGGCT)x15"],
    )
    def test_the_pair_rules_see_nothing(self, tract: str, svc: Services) -> None:
        """Not a tuning miss -- structural. The pair extension clamps an
        overlapping match back to its period, so a 90 bp array of a 6 bp unit is
        reported as a 6 bp match and dropped under every min_len."""
        c = with_tract(tract)
        ctx = context()
        assert not DirectRepeats().evaluate(c, ctx, svc).breaches
        assert not NearPerfectRepeats().evaluate(c, ctx, svc).breaches
        assert not SynthesisRepeats().evaluate(c, ctx, svc).breaches

    @pytest.mark.parametrize(
        "tract",
        ["CAG" * 20, "AT" * 30, "CAGGCT" * 15],
        ids=["(CAG)x20", "(AT)x30", "(CAGGCT)x15"],
    )
    def test_and_this_rule_does(self, tract: str, svc: Services) -> None:
        assert hits(with_tract(tract), svc)


class TestTractGeometry:
    def test_reports_the_tract_at_its_minimal_period(self) -> None:
        assert tandem_tracts("CAGGCT" * 15, min_length=WARN_TRACT_BP) == [(0, 90, 6)]

    def test_an_array_is_not_reported_once_per_multiple_of_its_unit(self) -> None:
        """(CAG)x20 is period 3 and also period 6 and also period 9. Reporting
        each would be three findings that disagree about the unit."""
        assert tandem_tracts("CAG" * 20, min_length=WARN_TRACT_BP) == [(0, 60, 3)]

    def test_a_homopolymer_reduces_to_period_one(self) -> None:
        """Which is how it gets handed to E1 rather than reported twice."""
        assert tandem_tracts("A" * 30, min_length=WARN_TRACT_BP) == [(0, 30, 1)]

    def test_the_length_filter_cannot_drop_a_needed_container(self) -> None:
        """Containment requires the covering tract to span the covered one, so a
        container is never shorter -- which is what makes the performance filter
        safe rather than approximate."""
        for min_length in (0, WARN_TRACT_BP):
            assert tandem_tracts("AT" * 30, min_length=min_length)[0][2] == 2

    def test_the_finding_carries_the_unit_and_copy_count(self, svc: Services) -> None:
        """The unit is a ROTATION of the planted one and the tract is at least as
        long as planted, not exactly.

        Tracts are reported maximally, so flanking bases that happen to continue
        the period are part of the array -- here the filler ends in CT, which
        extends a (CAGGCT)x15 tract two bases to the left and rotates the unit
        read from its true start. Asserting the planted string would be asserting
        that the scan is NOT maximal, which is the same trap the ITR exemption
        hit when a 145 bp palindrome reported as a 147 bp stem.
        """
        breach = hits(with_tract("CAGGCT" * 15), svc)[0]
        unit = str(breach.detail["unit"])
        assert breach.detail["unit_bp"] == 6.0
        assert unit in "CAGGCT" * 2, f"{unit!r} is not a rotation of CAGGCT"
        assert breach.detail["tract_bp"] >= 90.0
        assert float(breach.detail["copies"]) >= 15.0


class TestBands:
    def test_a_tract_over_the_hard_limit_fails(self, svc: Services) -> None:
        c = with_tract("CAGGCT" * 25)  # 150 bp
        result = ShortTandemRepeats().evaluate(c, context(), svc)
        assert not result.passes
        assert result.breaches[0].detail["severity"] == "hard"

    def test_a_tract_under_it_is_reported_without_failing(self, svc: Services) -> None:
        result = ShortTandemRepeats().evaluate(with_tract("CAG" * 20), context(), svc)
        assert result.breaches
        assert result.passes

    def test_a_tract_under_the_warn_floor_is_not_reported(self, svc: Services) -> None:
        assert not hits(with_tract("CAG" * 4), svc)  # 12 bp

    def test_severity_rises_with_copy_number(self, svc: Services) -> None:
        short = hits(with_tract("CAG" * 9), svc)[0]
        long = hits(with_tract("CAG" * 30), svc)[0]
        assert long.magnitude > short.magnitude

    def test_the_hard_limit_is_above_the_warn_floor(self) -> None:
        assert HARD_TRACT_BP > WARN_TRACT_BP


class TestHomopolymersBelongToE1:
    """A homopolymer is a tandem array of period 1, and reporting it here at a
    20 bp threshold would sit next to E1's 9 nt finding and contradict it."""

    def test_a_long_homopolymer_is_e1s_alone(self, svc: Services) -> None:
        c = with_tract("A" * 30)
        assert Homopolymers().evaluate(c, context(), svc).breaches
        assert not hits(c, svc)

    def test_the_scan_still_finds_it_in_order_to_exclude_it(self) -> None:
        """Excluding period 1 from the SCAN would leave A x30 looking like a
        period-2 array and report it after all."""
        assert tandem_tracts("A" * 30, min_length=WARN_TRACT_BP)[0][2] == 1

    def test_refusing_to_scan_unit_one_is_refused(self) -> None:
        with pytest.raises(ValueError, match="e1_homopolymers"):
            ShortTandemRepeats(max_unit=1)


class TestScope:
    def test_a_tract_in_the_backbone_is_not_the_vendors_problem(self, svc: Services) -> None:
        c = construct(dna(400, 41), "CAGGCT" * 25 + dna(300, 45))
        assert not hits(c, svc)

    def test_findings_carry_parent_construct_coordinates(self, svc: Services) -> None:
        c = construct(dna(200, 41) + "CAGGCT" * 15 + dna(200, 43), dna(300, 45), cds_start=100)
        breach = hits(c, svc)[0]
        assert c.overlaps_editable(breach.interval)
        assert breach.interval.start >= 100

    def test_a_forced_tract_is_still_reported_as_fixable(self, svc: Services) -> None:
        """poly-Gln cannot be removed, but alternating CAA/CAG breaks the period,
        so the finding must not be routed to the advisor as a dead end."""
        breach = hits(with_tract("CAG" * 20), svc)[0]
        assert breach.fixable_by_codon_choice
        assert "alternate synonymous codons" in breach.message


class TestContract:
    def test_it_is_hard_repair_and_carries_no_objective_weight(self) -> None:
        assert ShortTandemRepeats.enforcement is Enforcement.HARD_REPAIR
        assert ShortTandemRepeats.default_weight == 0.0

    def test_it_steers_below_the_dispersed_repeat_rules(self) -> None:
        """A tract is usually forced by the protein, so the repair is to break
        the period rather than to steer away from the residues."""
        assert 0.0 < ShortTandemRepeats.steering_weight < SynthesisRepeats.steering_weight

    def test_it_is_not_a_lattice_rule(self) -> None:
        assert ShortTandemRepeats().lattice_terms(None) is None

    def test_it_applies_in_every_context(self) -> None:
        for modality in Modality:
            assert ShortTandemRepeats().gate(slot(modality=modality))

    def test_absurd_parameters_are_refused(self) -> None:
        with pytest.raises(ValueError, match="e1_homopolymers"):
            ShortTandemRepeats(max_unit=1)
        with pytest.raises(ValueError, match="two copies"):
            ShortTandemRepeats(max_unit=6, warn_tract=8)
        with pytest.raises(ValueError, match="must not be below"):
            ShortTandemRepeats(warn_tract=50, hard_tract=30)
        with pytest.raises(ValueError, match="unknown vendor"):
            ShortTandemRepeats(vendors=VendorSelection.of("acme"))

    def test_it_is_registered_under_its_brief_row(self) -> None:
        assert get("e7_short_tandem_repeats") is ShortTandemRepeats
        assert ShortTandemRepeats.brief_ref == "2.E7"
