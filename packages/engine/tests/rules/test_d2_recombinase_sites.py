"""D2: recombination sites, and the two ways a both-strand scan misreports them."""

from __future__ import annotations

import pytest
from bt5.core.context import Modality
from bt5.core.registry import discover, get
from bt5.core.spec import Enforcement, RepairPolicy
from bt5.core.types import (
    Construct,
    Interval,
    Segment,
    SegmentKind,
    Topology,
    reverse_complement,
)
from bt5.rules.catalog.d2_recombinase_sites import (
    MAG_PARTIAL,
    MAG_SITE,
    RecombinaseSites,
)
from conftest import context, slot

discover()

#: Local helpers rather than conftest additions -- that file is shared with S3.

SPACER = "ATGTATGC"
LOXP = "ATAACTTCGTATA" + SPACER + "TATACGAAGTTAT"
FRT = "GAAGTTCCTATTC" + SPACER + "GTATAGGAACTTC"
ATTB1 = "ACAAGTTTGTACAAAAAAGCAGGCT"
ATTB2 = "ACCCAGCTTTCTTGTACAAAGTGGT"


def whole(seq: str, topology: Topology = Topology.LINEAR) -> Construct:
    return Construct(
        sequence=seq,
        topology=topology,
        segments=(Segment(Interval(0, len(seq)), SegmentKind.DESIGNABLE_CDS, "cds"),),
    )


def hits(rule: RecombinaseSites, c: Construct) -> list:
    return list(rule.evaluate(c, context(), None).breaches)  # type: ignore[arg-type]


def sites(rule: RecombinaseSites, c: Construct) -> list[str]:
    return [str(b.detail["site"]) for b in hits(rule, c)]


class TestDetection:
    @pytest.mark.parametrize(
        ("name", "site"),
        [
            ("loxP_family", LOXP),
            ("FRT", FRT),
            ("Gateway_attB1", ATTB1),
            ("Gateway_attB2", ATTB2),
        ],
    )
    def test_each_site_is_found_once(self, name: str, site: str) -> None:
        assert sites(RecombinaseSites(), whole("CCC" + site + "CCC")) == [name]

    def test_a_site_on_the_reverse_strand_is_found(self) -> None:
        assert sites(RecombinaseSites(), whole("CCC" + reverse_complement(ATTB1) + "CCC")) == [
            "Gateway_attB1"
        ]

    def test_a_clean_construct_passes(self) -> None:
        c = whole("ACGTACGTACGTACGTACGTACGT")
        assert hits(RecombinaseSites(), c) == []
        assert RecombinaseSites().evaluate(c, context(), None).passes  # type: ignore[arg-type]


class TestNoDoubleCounting:
    """loxP's arms are exact reverse complements of each other, so the whole 34 bp
    pattern matches its own reverse complement at the same coordinates. Keying the
    dedup on (span, strand) reports one physical site as two."""

    def test_a_palindromic_site_is_one_breach_not_two(self) -> None:
        assert sites(RecombinaseSites(), whole("CCC" + LOXP + "CCC")) == ["loxP_family"]

    def test_the_arms_of_a_complete_site_are_not_also_reported(self) -> None:
        """The two half-sites are that site again in smaller pieces."""
        assert "loxP_half_site" not in sites(RecombinaseSites(), whole("CCC" + LOXP + "CCC"))

    def test_a_lone_half_site_still_is(self) -> None:
        """An orphan arm is the finding worth having -- it does not arise by chance
        either, and it is what a partially-deleted cassette leaves behind."""
        breaches = hits(RecombinaseSites(), whole("CCC" + "ATAACTTCGTATA" + "CCCAAA"))
        assert [str(b.detail["site"]) for b in breaches] == ["loxP_half_site"]
        assert breaches[0].magnitude == MAG_PARTIAL

    def test_a_lone_half_site_does_not_fail_the_rule(self) -> None:
        """It cannot recombine on its own, so it is reported without blocking."""
        c = whole("CCC" + "ATAACTTCGTATA" + "CCCAAA")
        assert RecombinaseSites().evaluate(c, context(), None).passes  # type: ignore[arg-type]

    def test_partials_can_be_switched_off(self) -> None:
        c = whole("CCC" + "ATAACTTCGTATA" + "CCCAAA")
        assert hits(RecombinaseSites(report_partials=False), c) == []


class TestCircular:
    #: A loxP split across the origin: second arm first, first arm and spacer last.
    WRAPPED = "TATACGAAGTTAT" + "CCCAAACCCAAA" + "ATAACTTCGTATA" + SPACER

    def test_a_site_spanning_the_origin_is_found(self) -> None:
        c = whole(self.WRAPPED, Topology.CIRCULAR)
        breaches = hits(RecombinaseSites(), c)
        assert [str(b.detail["site"]) for b in breaches] == ["loxP_family"]
        assert breaches[0].interval.end > len(self.WRAPPED), "a wrapping interval"

    def test_the_same_sequence_linear_finds_only_the_two_arms(self) -> None:
        assert sites(RecombinaseSites(), whole(self.WRAPPED)) == [
            "loxP_half_site",
            "loxP_half_site",
        ]

    def test_a_wrapping_site_still_swallows_its_own_arms(self) -> None:
        """Containment has to be by residue. A site stored as [25, 59) on a 46 nt
        plasmid covers 25-45 and 0-12, and its first arm is stored as [0, 13) -- which
        `outer.start <= inner.start` rejects, reporting a lone half-site of the very
        site it belongs to."""
        assert "loxP_half_site" not in sites(
            RecombinaseSites(), whole(self.WRAPPED, Topology.CIRCULAR)
        )


class TestBxb1Collision:
    """brief.md:99 states the collision but gives no attB/attP sequence, so the rule
    surfaces it rather than pretending to detect a Bxb1 site."""

    def test_a_site_carrying_bsai_says_so(self) -> None:
        lox = "ATAACTTCGTATA" + "GGTCTCGC" + "TATACGAAGTTAT"
        breach = hits(RecombinaseSites(), whole("CCC" + lox + "CCC"))[0]
        assert breach.detail["bsai_collision"] == "yes"
        assert "GGTCTC" in breach.message

    def test_an_ordinary_site_does_not(self) -> None:
        breach = hits(RecombinaseSites(), whole("CCC" + LOXP + "CCC"))[0]
        assert breach.detail["bsai_collision"] == "no"

    def test_the_collision_is_declared_structurally(self) -> None:
        assert "d1_restriction_sites" in RecombinaseSites.conflicts_with


class TestCheckOnly:
    """brief.md:98 grades this H, check-only."""

    def test_it_is_hard_check(self) -> None:
        assert RecombinaseSites.enforcement is Enforcement.HARD_CHECK

    @pytest.mark.parametrize(
        "modality",
        [Modality.LENTIVIRAL, Modality.AAV, Modality.BACTERIAL_EXPRESSION, Modality.IVT_MRNA],
    )
    def test_it_applies_in_every_modality(self, modality: Modality) -> None:
        table = 11 if modality is Modality.BACTERIAL_EXPRESSION else 1
        s = slot(modality=modality, table=table) if table == 1 else slot(modality=modality)
        assert RecombinaseSites().gate(s)
        assert RecombinaseSites().enforcement_for(s) is Enforcement.HARD_CHECK

    def test_nothing_is_offered_to_the_solver(self) -> None:
        """A 25-48 bp site is present because someone put it there. Marking it fixable
        sends the solver after bases it cannot move until it reports infeasible on a
        design that was fine (core/spec.py:154-166)."""
        breaches = hits(RecombinaseSites(), whole("CCC" + LOXP + "CCC"))
        assert breaches
        assert all(not b.fixable_by_codon_choice for b in breaches)

    def test_it_is_not_a_lattice_rule(self) -> None:
        """The 8 nt spacer would expand to 4**8 = 65,536 patterns against a
        MAX_PATTERN_EXPANSION of 1,024, and a HARD_CHECK rule never reaches the
        automaton that would consume them."""
        assert RecombinaseSites().lattice_terms(context()) is None  # type: ignore[arg-type]

    def test_a_hard_rule_carries_no_weight(self) -> None:
        assert RecombinaseSites.default_weight == 0.0
        assert RecombinaseSites.steering_weight == 0.0

    def test_the_repair_policy_is_inert_and_declared(self) -> None:
        assert RecombinaseSites.repair is RepairPolicy.SINGLE_PASS

    def test_it_is_registered_under_its_brief_row(self) -> None:
        assert get("d2_recombinase_sites").brief_ref == "2.D2"


class TestMagnitudes:
    def test_a_full_site_outranks_a_half_site(self) -> None:
        assert MAG_SITE > MAG_PARTIAL

    def test_a_full_site_fails_the_rule(self) -> None:
        c = whole("CCC" + LOXP + "CCC")
        assert not RecombinaseSites().evaluate(c, context(), None).passes  # type: ignore[arg-type]
