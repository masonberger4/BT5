"""E9: whether the vendor will take the order at all.

The rule under test cannot be satisfied by changing the sequence, so the thing
worth testing hardest is not detection -- a length comparison is not subtle --
but the three properties that make the finding useful rather than merely true:
that it routes away from the solver, that it names somewhere the order COULD
go, and that it refuses to guess on the one bound the vendors leave ambiguous.
"""

from __future__ import annotations

import numpy as np
import pytest
from bt5.core.registry import discover, get
from bt5.core.services import Services
from bt5.core.spec import Enforcement
from bt5.core.types import Construct
from bt5.rules.catalog.e9_length_tiers import AMBIGUOUS, NARROWED, UNORDERABLE, LengthTiers
from bt5.rules.vendors import PROFILES, VendorSelection, orderable_keys
from conftest import construct, context


def dna(n: int, seed: int = 3) -> str:
    rng = np.random.default_rng(seed)
    return "".join("ACGT"[i] for i in rng.integers(0, 4, n))


@pytest.fixture
def svc() -> Services:
    from bt5.vector.kmers import ConstructKmerIndex

    return Services(
        fold=None,
        kmer=ConstructKmerIndex,
        tables=None,  # type: ignore[arg-type]
        rng=np.random.default_rng(42),
    )


def run(insert_bp: int, svc: Services, vendor: str = "twist_gene_fragment"):
    c: Construct = construct(dna(insert_bp), dna(300, 11))
    return LengthTiers(vendors=VendorSelection.of(vendor)).evaluate(c, context(), svc)


class TestTheFloor:
    """The bound nobody remembers, and the reason this rule exists."""

    def test_an_insert_under_300_bp_cannot_be_ordered_from_twist(self, svc: Services) -> None:
        ev = run(200, svc)
        assert not ev.passes
        assert ev.raw_score == UNORDERABLE
        assert "below the 300-5000 bp" in ev.breaches[0].message

    def test_and_it_names_the_one_configuration_that_would_take_it(self, svc: Services) -> None:
        """A finding no codon can act on is only useful if it says where to go."""
        breach = run(200, svc).breaches[0]
        assert breach.detail["alternatives"] == "idt_gblocks"
        assert "idt_gblocks" in breach.message

    def test_under_125_bp_nothing_accepts_it(self, svc: Services) -> None:
        breach = run(100, svc).breaches[0]
        assert breach.detail["alternatives"] == "none"
        assert "no configured vendor accepts this length" in breach.message

    def test_a_300_bp_insert_is_exactly_in_range(self, svc: Services) -> None:
        """The floor is inclusive; 300 is orderable and 299 is not."""
        assert run(300, svc).passes
        assert not run(299, svc).passes


class TestTheCeiling:
    def test_eblocks_stop_at_1500(self, svc: Services) -> None:
        ev = run(2000, svc, vendor="idt_eblocks")
        assert not ev.passes
        assert "above the 300-1500 bp" in ev.breaches[0].message

    def test_and_the_same_insert_is_fine_as_a_gene_fragment(self, svc: Services) -> None:
        assert run(2000, svc).passes
        assert "twist_gene_fragment" in str(
            run(2000, svc, vendor="idt_eblocks").breaches[0].detail["alternatives"]
        )


class TestItNeverGoesToTheSolver:
    def test_no_finding_is_fixable_by_codon_choice(self, svc: Services) -> None:
        """Every synonymous codon is three bases, so the mutation space contains
        no sequence of a different length. Sending this to the solver would
        exhaust that space and report infeasible on a design that is fine."""
        for bp in (100, 200, 6000):
            for b in run(bp, svc).breaches:
                assert b.fixable_by_codon_choice is False

    def test_the_rule_is_hard_check_and_never_weighted(self) -> None:
        discover()
        spec = get("e9_length_tiers")
        assert spec.enforcement is Enforcement.HARD_CHECK
        assert spec.default_weight == 0.0
        assert spec.steering_weight == 0.0, "no codon choice moves a length"


class TestTheAmbiguousBound:
    """The vendors publish a range without saying insert or insert-plus-adapters."""

    def test_just_under_the_floor_with_adapters_is_reported_as_unresolved(
        self, svc: Services
    ) -> None:
        # 280 bp of insert is under the 300 bp floor; 280 + 44 of adapter is over it.
        ev = run(280, svc, vendor="twist_gene_fragment_adapter_on")
        assert ev.raw_score == AMBIGUOUS
        assert ev.passes, "unresolved is a finding, not a proven rejection"
        assert "does not say whether it applies" in ev.breaches[0].message

    def test_just_under_the_ceiling_with_adapters_too(self, svc: Services) -> None:
        ev = run(4980, svc, vendor="twist_gene_fragment_adapter_on")
        assert ev.raw_score == AMBIGUOUS
        assert ev.breaches[0].detail["with_adapters_bp"] == 5024.0

    def test_an_adapter_free_order_is_never_ambiguous(self, svc: Services) -> None:
        """With no adapters the two readings are the same number."""
        assert run(280, svc).raw_score == UNORDERABLE
        assert run(4980, svc).passes

    def test_well_inside_the_range_adapters_change_nothing(self, svc: Services) -> None:
        assert run(900, svc, vendor="twist_gene_fragment_adapter_on").passes


class TestMultiVendorDisjunctive:
    """A selection of vendors is a menu, not a contract to satisfy all of them.

    #43 V3 decided multi-select is disjunctive: a fragment one selected vendor
    will build is not blocked because another refuses it. The rule still emits one
    breach per fragment -- it routes, it does not fail. Only when EVERY selected
    vendor refuses does the length become unorderable.
    """

    def _multi(self, insert_bp: int, svc: Services, *keys: str):
        c: Construct = construct(dna(insert_bp), dna(300, 11))
        return LengthTiers(vendors=VendorSelection.of(*keys)).evaluate(c, context(), svc)

    def test_over_gblocks_ceiling_but_under_twist_does_not_block(self, svc: Services) -> None:
        # 4002 bp: over gBlocks' 3000, inside Twist's 5000.
        ev = self._multi(4002, svc, "idt_gblocks", "twist_gene_fragment")
        assert ev.passes, "a fragment one selected vendor builds is not a rejection"
        assert len(ev.breaches) == 1, "one fragment, one routing finding -- no per-vendor loop"
        b = ev.breaches[0]
        assert b.magnitude == NARROWED
        assert b.detail["refusing"] == "idt_gblocks"
        assert b.detail["accepting"] == "twist_gene_fragment"
        assert "outside the range of idt_gblocks" in b.message
        assert "twist_gene_fragment, which you also selected, accepts it" in b.message

    def test_under_twist_floor_but_gblocks_takes_it_does_not_block(self, svc: Services) -> None:
        # 200 bp: under Twist's 300 floor, inside gBlocks' 125-3000. The floor case,
        # run in reverse -- gBlocks is the one that reaches lower.
        ev = self._multi(200, svc, "twist_gene_fragment", "idt_gblocks")
        assert ev.passes
        b = ev.breaches[0]
        assert b.magnitude == NARROWED
        assert b.detail["refusing"] == "twist_gene_fragment"
        assert b.detail["accepting"] == "idt_gblocks"

    def test_when_every_selected_vendor_refuses_it_is_unorderable(self, svc: Services) -> None:
        # 6000 bp: over both gBlocks' 3000 and Twist's 5000.
        ev = self._multi(6000, svc, "idt_gblocks", "twist_gene_fragment")
        assert not ev.passes
        b = ev.breaches[0]
        assert b.magnitude == UNORDERABLE
        assert b.detail["accepting"] == ""
        assert b.detail["refusing"] == "idt_gblocks, twist_gene_fragment"
        assert "outside the range of every selected configuration" in b.message
        assert "idt_gblocks, twist_gene_fragment" in b.message

    def test_a_length_both_accept_produces_no_breach(self, svc: Services) -> None:
        # 1500 bp: inside gBlocks (125-3000) and inside Twist (300-5000).
        ev = self._multi(1500, svc, "idt_gblocks", "twist_gene_fragment")
        assert ev.passes
        assert ev.breaches == ()


class TestTheRegistry:
    """The bug this guards against already happened once.

    Two vendor namespaces that validated their own keys independently let a run
    be spec'd for one vendor's limits and another's adapters, and each lookup
    succeeded. The old test here compared the two dicts' key sets; there is now
    one dict, so the same guarantee is asserted against the structure that
    replaced them -- an orderable profile carries EVERY vendor fact or none of
    them, and half a profile is refused at construction.
    """

    def test_an_orderable_profile_is_complete(self) -> None:
        for key in orderable_keys():
            p = PROFILES[key]
            assert p.adapters.vendor == key, "adapters must name their own configuration"
            assert p.length_bp is not None
            assert p.homopolymer_at is not None
            assert p.homopolymer_gc is not None
            assert p.global_gc is not None
            assert p.last_verified
            assert p.notes.strip()

    def test_a_half_specified_profile_cannot_be_constructed(self) -> None:
        """The structural replacement for comparing two dicts' key sets."""
        from bt5.rules.fragment import IDT_GBLOCKS
        from bt5.rules.vendors import VendorProfile

        with pytest.raises(ValueError, match="half-specified"):
            VendorProfile(
                key="idt_gblocks",
                vendor="IDT",
                product="gBlocks Gene Fragment",
                adapters=IDT_GBLOCKS,
                length_bp=(125, 3000),  # a length but no run limits and no GC band
            )

    def test_a_profile_cannot_carry_another_configurations_adapters(self) -> None:
        from bt5.rules.fragment import TWIST_ADAPTER_ON
        from bt5.rules.vendors import VendorProfile

        with pytest.raises(ValueError, match="one configuration, one name"):
            VendorProfile(
                key="idt_gblocks",
                vendor="IDT",
                product="gBlocks Gene Fragment",
                adapters=TWIST_ADAPTER_ON,
            )

    def test_none_is_not_an_orderable_configuration(self) -> None:
        with pytest.raises(ValueError, match="not orderable from anyone"):
            LengthTiers(vendors=VendorSelection.of("none"))

    def test_an_unknown_vendor_is_refused(self) -> None:
        with pytest.raises(ValueError, match="unknown vendor"):
            LengthTiers(vendors=VendorSelection.of("genscript_gentitan"))

    def test_every_range_is_non_empty_and_positive(self) -> None:
        for key in orderable_keys():
            lo, hi = PROFILES[key].length_bp
            assert 0 < lo < hi, key
