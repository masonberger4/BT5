"""E1: homopolymer limits, and the cases a linear per-codon check cannot see."""

from __future__ import annotations

import pytest
from bt5.core.context import Modality
from bt5.core.registry import discover, get
from bt5.core.services import Services
from bt5.core.spec import Enforcement
from bt5.core.types import Construct, Interval, Segment, SegmentKind, Topology
from bt5.rules.catalog.e1_homopolymers import Homopolymers, _maximal_runs
from bt5.rules.vendors import PROFILES
from conftest import construct, context, slot

discover()


def evaluate(rule: Homopolymers, c: Construct, svc: Services | None = None):
    return rule.evaluate(c, context(), svc)  # type: ignore[arg-type]


class TestMaximalRuns:
    def test_linear_runs(self) -> None:
        assert list(_maximal_runs("AAGGGT", circular=False)) == [
            (0, 2, "A"),
            (2, 3, "G"),
            (5, 1, "T"),
        ]

    def test_a_run_across_the_origin_is_one_run(self) -> None:
        """Two halves reported separately would let a 16 nt run past a 9 nt
        limit as two compliant runs of 8."""
        runs = {start: (length, base) for start, length, base in _maximal_runs("AAGGGAAA", True)}
        assert runs[5] == (5, "A"), "the wrapping run starts at 5 and is 5 long"
        assert 0 not in runs, "and its tail must not also be reported as its own run"

    def test_a_construct_of_one_base_is_one_run(self) -> None:
        assert list(_maximal_runs("AAAA", circular=True)) == [(0, 4, "A")]

    def test_the_same_sequence_linear_has_two_runs(self) -> None:
        runs = list(_maximal_runs("AAGGGAAA", circular=False))
        assert [r[0] for r in runs] == [0, 2, 5]

    def test_empty(self) -> None:
        assert list(_maximal_runs("", circular=True)) == []


class TestLatticeTerms:
    def test_forbids_the_first_run_that_is_too_long(self) -> None:
        terms = Homopolymers(max_at_run=9, max_gc_run=5).lattice_terms(None)
        assert terms.forbidden == ("A" * 10, "G" * 6)

    def test_lists_forward_motifs_only(self) -> None:
        """The solver closes the set under reverse complement, so poly-T and
        poly-C come free; listing them would report every run twice."""
        forbidden = set(Homopolymers().lattice_terms(None).forbidden)
        assert not any(m[0] in "TC" for m in forbidden)

    def test_a_lattice_rule_carries_no_steering_weight(self) -> None:
        """HARD_LATTICE is unreachable by construction; steering is a no-op."""
        assert Homopolymers.steering_weight == 0.0
        assert Homopolymers.enforcement is Enforcement.HARD_LATTICE


class TestEvaluate:
    def test_reports_one_breach_per_run_with_its_true_length(self) -> None:
        """A run of 12 contains three matches of A*10. Reporting three findings
        for one physical run makes the conflict panel unreadable."""
        c = construct("ATG" + "A" * 12 + "TAA")
        breaches = evaluate(Homopolymers(), c).breaches
        assert len(breaches) == 1
        assert breaches[0].interval == Interval(3, 15)
        assert "12 nt" in breaches[0].message

    def test_magnitude_grows_with_the_overrun(self) -> None:
        """One base over is a surcharge; twice the limit is a failed synthesis."""
        small = evaluate(Homopolymers(), construct("ATG" + "A" * 10 + "TAA")).breaches[0]
        large = evaluate(Homopolymers(), construct("ATG" + "A" * 20 + "TAA")).breaches[0]
        assert large.magnitude > small.magnitude

    def test_at_and_gc_have_different_limits(self) -> None:
        """The 10-vs-6 asymmetry is the clearest vendor evidence that G/C runs
        are chemically worse, so one shared limit would encode wrong physics."""
        rule = Homopolymers()
        assert not evaluate(rule, construct("ATG" + "A" * 9 + "TAA")).breaches
        assert evaluate(rule, construct("ATG" + "G" * 6 + "TAA")).breaches, (
            "six G is already over the limit that nine A is not"
        )

    def test_poly_t_is_caught_even_though_only_poly_a_is_listed(self) -> None:
        breaches = evaluate(Homopolymers(), construct("ATG" + "T" * 12 + "TAA")).breaches
        assert len(breaches) == 1
        assert breaches[0].detail["base_class"] == "A/T"

    def test_a_run_spanning_the_origin_is_caught(self) -> None:
        """The case that justifies evaluating on the assembled circular
        construct: a linear scan sees two compliant runs of 6."""
        seq = "T" * 6 + "ACGACGACG" + "T" * 6
        c = Construct(
            seq,
            Topology.CIRCULAR,
            (Segment(Interval(0, len(seq)), SegmentKind.DESIGNABLE_CDS, "cds"),),
        )
        breaches = evaluate(Homopolymers(), c).breaches
        assert len(breaches) == 1
        assert breaches[0].interval.wraps(len(seq))
        assert "12 nt" in breaches[0].message

    def test_the_same_sequence_linear_is_clean(self) -> None:
        seq = "T" * 6 + "ACGACGACG" + "T" * 6
        c = Construct(
            seq,
            Topology.LINEAR,
            (Segment(Interval(0, len(seq)), SegmentKind.DESIGNABLE_CDS, "cds"),),
        )
        assert not evaluate(Homopolymers(), c).breaches

    def test_a_run_across_the_cds_backbone_junction_is_caught(self) -> None:
        """Half in the CDS, half in the vector -- no per-codon check sees it,
        and it is partly recodable, so it is not the user's problem alone."""
        c = construct("ATGAAAAAA", "AAAAAACCC")
        breaches = evaluate(Homopolymers(), c).breaches
        assert len(breaches) == 1
        assert breaches[0].fixable_by_codon_choice, "the CDS half can be recoded"

    def test_a_run_wholly_in_the_backbone_is_reported_but_unfixable(self) -> None:
        c = construct("ATGCTGTAA", "GGGGGGGGCCC")
        breaches = evaluate(Homopolymers(), c).breaches
        assert breaches, "a run in the user's own vector is still worth reporting"
        assert not any(b.fixable_by_codon_choice for b in breaches), (
            "and no codon can shorten it, so the solver must not chase it"
        )

    def test_a_clean_construct_passes(self) -> None:
        assert evaluate(Homopolymers(), construct("ATGACGTACGTACGTTAA")).passes


class TestVendors:
    def test_the_vendor_changes_the_limits(self) -> None:
        """A 12 nt A-run is over IDT's 9 and under Twist's 13."""
        c = construct("ATG" + "A" * 12 + "TAA")
        assert evaluate(Homopolymers(vendor="idt_gblocks"), c).breaches
        assert not evaluate(Homopolymers(vendor="twist_gene_fragment"), c).breaches

    def test_explicit_limits_override_the_vendor(self) -> None:
        rule = Homopolymers(max_at_run=20, vendor="idt_gblocks")
        assert rule.max_at_run == 20
        assert rule.max_gc_run == PROFILES["idt_gblocks"].homopolymer_gc, (
            "the other limit still applies"
        )

    def test_an_unknown_vendor_is_refused(self) -> None:
        with pytest.raises(ValueError, match="unknown vendor"):
            Homopolymers(vendor="acme")

    def test_no_vendor_chosen_is_refused_rather_than_guessed(self) -> None:
        """`none` reaches E4-E7 legitimately -- they need only the synthesis
        scope. Here it is a request to check a sequence against the limits of a
        vendor nobody picked, and guessing which vendor to answer with is exactly
        the split default this registry exists to end."""
        with pytest.raises(ValueError, match="not orderable from anyone"):
            Homopolymers(vendor="none")

    def test_the_default_is_the_strictest_shipped_configuration(self) -> None:
        """The unified default must not LOOSEN a lattice bound.

        E1 is HARD_LATTICE, so its output when a run is permitted is silence --
        there is no finding to notice. That asymmetry is why the shared default
        resolves toward the tighter run limits even though the length range it
        brings with it is not the one most orders use.
        """
        default = Homopolymers()
        for p in PROFILES.values():
            if not p.is_orderable:
                continue
            assert default.max_at_run <= p.homopolymer_at, p.key
            assert default.max_gc_run <= p.homopolymer_gc, p.key

    def test_an_absurd_limit_is_refused(self) -> None:
        with pytest.raises(ValueError, match="forbid ordinary sequence"):
            Homopolymers(max_at_run=2)


def test_it_applies_in_every_modality() -> None:
    """An IVT mRNA template is still ordered as DNA, so there is no modality
    where a vendor's synthesis limits stop applying."""
    rule = Homopolymers()
    for modality in Modality:
        assert rule.gate(slot(modality=modality))


def test_it_is_registered_under_its_brief_row() -> None:
    assert get("e1_homopolymers") is Homopolymers
    assert Homopolymers.brief_ref == "2.E1"
