"""Autodiscovery and the import-time contract."""

from __future__ import annotations

from bt5.core.registry import all_specs, discover, get
from bt5.core.spec import Enforcement, Evidence
from bt5.core.types import Interval


def test_discovery_finds_the_reference_rules() -> None:
    discover()
    ids = {s.id for s in all_specs()}
    assert {"d1_restriction_sites", "e2_gc_band"} <= ids


def test_registry_is_sorted_for_stable_output() -> None:
    discover()
    ids = [s.id for s in all_specs()]
    assert ids == sorted(ids)


def test_lookup_by_id() -> None:
    discover()
    assert get("e2_gc_band").enforcement is Enforcement.HARD_REPAIR


def test_gc_band_encodes_the_published_twist_bound_not_the_folklore_one() -> None:
    """The widely repeated 35-65% 50bp window figure has no vendor source; Twist's
    published High-Complexity trigger is 10-90%."""
    discover()
    spec = get("e2_gc_band")
    assert spec.evidence is Evidence.VENDOR_ASSERTED
    assert any("10%" in c.label or "90%" in c.label for c in spec.citations)


def test_restriction_rule_declares_forward_motifs_only() -> None:
    """The solver closes the pattern set under reverse complement, so a rule must
    not list both a motif and its distinct reverse complement or the site would
    be counted twice.

    Note every classic six-cutter is palindromic (revcomp == itself), so this has
    to be phrased as "no distinct revcomp pair", not "no revcomp present".
    """
    from bt5.core.types import reverse_complement

    discover()
    terms = get("d1_restriction_sites")().lattice_terms(None)
    assert terms.forbidden, "a HARD_LATTICE rule must declare motifs"

    forbidden = set(terms.forbidden)
    for motif in forbidden:
        rc = reverse_complement(motif)
        if rc != motif:
            assert rc not in forbidden, (
                f"{motif} and its reverse complement {rc} are both listed; the "
                f"solver already closes the set under revcomp"
            )


def test_a_backbone_site_is_reported_but_not_offered_to_the_solver() -> None:
    """`fixable_by_codon_choice` has to be COMPUTED from where the hit landed.

    Both reference rules scan the whole assembled construct, so they find things
    in the user's own backbone. Those are worth reporting and impossible to
    recode, and marking them fixable sends the solver after bases it may not
    touch until the mutation space is empty -- which surfaces as an
    infeasibility certificate for a design that was never infeasible.
    """
    from bt5.core.context import (
        BiosecurityVerdict,
        ContextSlot,
        DesignContext,
        HostId,
        Modality,
    )
    from bt5.core.types import Construct, Segment, SegmentKind, Topology

    discover()
    # EcoRI in the CDS, BamHI out in the backbone. Lengths kept a multiple of 3
    # so the CDS stays in frame.
    cds = "ATGAAAGAATTCAAACCCTAA"  # 21 nt, GAATTC at 6
    backbone = "TTTGGATCCTTT"  # GGATCC at 3 -> construct offset 24
    c = Construct(
        sequence=cds + backbone,
        topology=Topology.CIRCULAR,
        segments=(
            Segment(Interval(0, len(cds)), SegmentKind.DESIGNABLE_CDS, "cds"),
            Segment(Interval(len(cds), len(cds) + len(backbone)), SegmentKind.BACKBONE, "vector"),
        ),
    )
    ctx = DesignContext(
        slots=(ContextSlot("propagation", HostId.E_COLI_K12, Modality.PLASMID_TRANSIENT, 11),),
        cassette_orientation=1,
        seed=1,
        screen=BiosecurityVerdict("not_run"),
    )

    breaches = get("d1_restriction_sites")(("EcoRI", "BamHI")).evaluate(c, ctx, None).breaches
    by_enzyme = {b.detail["enzyme"]: b for b in breaches}
    assert {"EcoRI", "BamHI"} <= by_enzyme.keys(), "both sites must still be REPORTED"
    assert by_enzyme["EcoRI"].fixable_by_codon_choice, "inside the CDS: recode it"
    assert not by_enzyme["BamHI"].fixable_by_codon_choice, (
        "in the user's backbone: real, reported, and not the solver's to chase"
    )
