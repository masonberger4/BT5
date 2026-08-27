"""Autodiscovery and the import-time contract."""

from __future__ import annotations

from bt5.core.registry import all_specs, discover, get
from bt5.core.spec import Enforcement, Evidence


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
