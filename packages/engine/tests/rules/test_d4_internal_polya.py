"""D4: internal polyA, and the strand question that makes it worth having.

The expensive failure this rule guards is silent in both directions. An internal
polyA signal RAISES the expression number and CUTS functional titer 8-9x, so the
assay says the construct got better. And a scan that reads the wrong strand
returns clean on exactly the construct that has the problem.
"""

from __future__ import annotations

from bt5.core.context import Modality
from bt5.core.registry import discover, get
from bt5.core.spec import Enforcement
from bt5.core.types import reverse_complement
from bt5.rules.catalog.d4_internal_polya import (
    CANONICAL,
    InternalPolyA,
    _has_downstream_element,
)
from conftest import construct, context, slot

discover()

#: 24 nt of filler with no hexamer and no GU-rich element.
PAD = "CCAGCCAGCCAGCCAGCCAGCCAG"


def hits(rule: InternalPolyA, c, ctx):
    return rule.evaluate(c, ctx, None).breaches


def wrapped_hexamer(*, circular: bool):
    """Ends with AATA and starts with AA, so AATAAA exists only across the origin."""
    from bt5.core.types import Construct, Interval, Segment, SegmentKind, Topology

    seq = "AA" + PAD + "AATA"
    return Construct(
        seq,
        Topology.CIRCULAR if circular else Topology.LINEAR,
        (Segment(Interval(0, len(seq)), SegmentKind.DESIGNABLE_CDS, "cds"),),
    )


class TestStrand:
    """The regression guard for `strand_for` composing cassette orientation.

    `AATAAA` on the forward strand and `TTTATT` on the forward strand are the
    same physical signal read from opposite sides. Which one matters depends on
    which strand gets packaged, and that is the cassette's orientation composed
    with the slot's own preference -- not either one alone.
    """

    def forward_signal(self):
        return construct("ATG" + PAD + "AATAAA" + PAD + "TAA")

    def reverse_signal(self):
        """The forward strand carries TTTATT, so the SIGNAL is on the minus
        strand. Only a reverse-oriented cassette packages it."""
        return construct("ATG" + PAD + reverse_complement("AATAAA") + PAD + "TAA")

    def test_a_forward_cassette_finds_a_forward_signal(self) -> None:
        found = hits(InternalPolyA(), self.forward_signal(), context(cassette_orientation=1))
        assert any(b.detail["hexamer"] == "AATAAA" for b in found)

    def test_a_forward_cassette_ignores_the_reverse_one(self) -> None:
        found = hits(InternalPolyA(), self.reverse_signal(), context(cassette_orientation=1))
        assert not any(b.detail["hexamer"] == "AATAAA" for b in found), (
            "TTTATT on the packaged strand cannot polyadenylate anything"
        )

    def test_a_reverse_cassette_finds_the_reverse_signal(self) -> None:
        """The case that returned clean before `strand_for` composed the
        cassette orientation: the packaged genome is the reverse complement."""
        found = hits(InternalPolyA(), self.reverse_signal(), context(cassette_orientation=-1))
        assert any(b.detail["hexamer"] == "AATAAA" for b in found)

    def test_a_reverse_cassette_ignores_the_forward_one(self) -> None:
        found = hits(InternalPolyA(), self.forward_signal(), context(cassette_orientation=-1))
        assert not any(b.detail["hexamer"] == "AATAAA" for b in found)

    def test_orientation_alone_changes_the_answer(self) -> None:
        """Flip only the cassette orientation, touch nothing else."""
        c = self.reverse_signal()

        def canonical(orientation):
            return [
                b.detail["hexamer"]
                for b in hits(InternalPolyA(), c, context(cassette_orientation=orientation))
                if b.detail["canonical"] == "yes"
            ]

        assert canonical(1) == []
        assert canonical(-1) == ["AATAAA"]

    def test_the_reported_position_is_in_construct_coordinates(self) -> None:
        """A breach the user cannot find on their own map is not actionable,
        whichever strand found it."""
        c = self.reverse_signal()
        found = hits(InternalPolyA(), c, context(cassette_orientation=-1))
        breach = next(b for b in found if b.detail["hexamer"] == "AATAAA")
        assert c.sequence[breach.interval.start : breach.interval.start + 6] == reverse_complement(
            "AATAAA"
        )
        assert breach.interval.strand == -1, "and it says which strand it was found on"

    def test_a_slot_reading_the_antisense_composes_back_to_forward(self) -> None:
        """Both reversed is forward again: an antisense slot on a reverse
        cassette reads the cassette's own sense strand."""
        antisense = slot(strand_of_interest=-1)
        found = hits(
            InternalPolyA(), self.forward_signal(), context(antisense, cassette_orientation=-1)
        )
        assert any(b.detail["hexamer"] == "AATAAA" for b in found)


class TestDownstreamElement:
    """A hexamer alone is common; a hexamer with a downstream element cleaves."""

    def test_gt_repeat_escalates(self) -> None:
        with_dse = construct("ATG" + PAD + "AATAAA" + "CCCCCCCCCC" + "GTGTGTGT" + PAD + "TAA")
        breach = next(
            b for b in hits(InternalPolyA(), with_dse, context()) if b.detail["hexamer"] == "AATAAA"
        )
        assert breach.detail["downstream_element"] == "yes"
        assert breach.magnitude > 1.0

    def test_a_bare_hexamer_does_not_escalate(self) -> None:
        breach = next(
            b
            for b in hits(
                InternalPolyA(), construct("ATG" + PAD + "AATAAA" + PAD + "TAA"), context()
            )
            if b.detail["hexamer"] == "AATAAA"
        )
        assert breach.detail["downstream_element"] == "no"
        assert breach.magnitude == 1.0

    def test_the_window_is_bounded(self) -> None:
        """A GU-rich stretch 200 nt away is not this hexamer's downstream element."""
        far = construct("ATG" + PAD + "AATAAA" + "C" * 100 + "GTGTGTGT" + PAD + "TAA")
        breach = next(
            b for b in hits(InternalPolyA(), far, context()) if b.detail["hexamer"] == "AATAAA"
        )
        assert breach.detail["downstream_element"] == "no"

    def test_the_element_predicate(self) -> None:
        assert _has_downstream_element("AAGTGTAA")
        assert _has_downstream_element("AATGTGAA")
        assert _has_downstream_element("ACTTTTAC"), "a U-run of 4 in a 6 nt window"
        assert not _has_downstream_element("CCAGCCAGCCAG")


class TestSeverity:
    def test_canonical_hexamers_outrank_variants(self) -> None:
        canonical = hits(
            InternalPolyA(), construct("ATG" + PAD + "AATAAA" + PAD + "TAA"), context()
        )
        variant = hits(InternalPolyA(), construct("ATG" + PAD + "AATACA" + PAD + "TAA"), context())
        assert canonical[0].magnitude > variant[0].magnitude

    def test_variants_can_be_switched_off(self) -> None:
        c = construct("ATG" + PAD + "AATACA" + PAD + "TAA")
        assert hits(InternalPolyA(include_variants=True), c, context())
        assert not hits(InternalPolyA(include_variants=False), c, context())

    def test_a_variant_alone_still_passes(self) -> None:
        """Weak signals are reported, not treated as failures."""
        c = construct("ATG" + PAD + "AATACA" + PAD + "TAA")
        assert InternalPolyA().evaluate(c, context(), None).passes

    def test_a_canonical_hexamer_does_not(self) -> None:
        c = construct("ATG" + PAD + "AATAAA" + PAD + "TAA")
        assert not InternalPolyA().evaluate(c, context(), None).passes


class TestEnforcement:
    def test_hard_where_the_titer_cost_is_measured(self) -> None:
        rule = InternalPolyA()
        for modality in (Modality.LENTIVIRAL, Modality.AAV, Modality.GENOME_INTEGRATED):
            assert rule.enforcement_for(slot(modality=modality)).is_hard

    def test_soft_where_it_costs_a_little_expression(self) -> None:
        rule = InternalPolyA()
        assert rule.enforcement_for(slot(modality=Modality.PLASMID_TRANSIENT)) is Enforcement.SOFT

    def test_the_class_level_floor_is_soft_so_it_can_be_weighted(self) -> None:
        """The class attribute is the floor; `enforcement_for` is authoritative.
        A rule declared hard at class level could carry no objective weight at
        all, which would be wrong in the modalities where it is a preference."""
        assert InternalPolyA.enforcement is Enforcement.SOFT
        assert InternalPolyA.default_weight > 0.0
        assert InternalPolyA.weight_provenance.strip()

    def test_it_does_not_apply_to_bacterial_or_ivt(self) -> None:
        rule = InternalPolyA()
        assert not rule.gate(slot(modality=Modality.BACTERIAL_EXPRESSION))
        assert not rule.gate(slot(modality=Modality.IVT_MRNA))
        assert rule.gate(slot(modality=Modality.LENTIVIRAL))

    def test_it_is_not_a_lattice_rule(self) -> None:
        """`forbidden` is closed under reverse complement, which would also
        forbid TTTATT -- a sequence that cannot polyadenylate anything."""
        assert InternalPolyA().lattice_terms(None) is None


def test_a_hexamer_spanning_the_origin_is_found() -> None:
    """A signal assembled from the end and the start of a circular construct is
    invisible to any linear scan, and it is a real cleavage site: the plasmid
    does not know where the file's first character is.

    Found while writing these tests -- filler ending `...CCAGTAA` ahead of a CDS
    starting `ATG` spells the variant hexamer AGTAAA across position 0, which is
    why this is a test rather than a footnote.
    """
    c = wrapped_hexamer(circular=True)
    found = hits(InternalPolyA(), c, context())
    wrapping = [b for b in found if b.interval.start + 6 > c.length]
    assert wrapping, f"expected a hit across the origin, got {[b.interval for b in found]}"
    assert wrapping[0].detail["hexamer"] == "AATAAA"


def test_the_same_construct_linear_has_no_origin_hit() -> None:
    c = wrapped_hexamer(circular=False)
    assert not hits(InternalPolyA(), c, context())


def test_slot_role_travels_with_the_breach() -> None:
    """Two slots can disagree about which strand matters, so the conflict
    detector needs to know which context produced each finding."""
    c = construct("ATG" + PAD + "AATAAA" + PAD + "TAA")
    found = hits(InternalPolyA(), c, context(slot(role="target")))
    assert all(b.slot_role == "target" for b in found)


def test_it_is_registered_under_its_brief_row() -> None:
    assert get("d4_internal_polya") is InternalPolyA
    assert InternalPolyA.brief_ref == "2.D4"
    assert set(CANONICAL) == {"AATAAA", "ATTAAA"}
