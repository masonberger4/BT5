"""`design()` at its seams: the site, the flanks, the carried motif, the hashes.

These are the tests the walking skeleton left behind that are still about the
same seams, updated where the ranking increment changed what `design()` returns.
What the design REFUSES to claim moved to `test_increment.py`, which is where the
ranking, the panel, the baseline and the order file are exercised together --
this file stays about the parts of the pipeline the ranking did not touch.

The one subtlety worth reading twice is still the flank orientation. Tier A seeds
its automaton with the immutable backbone on either side of the insert so a
forbidden site formed half by the vector and half by the first codon is excluded
by construction. Those flanks are in CODING orientation, and on a reverse-strand
site a sign error there comes back silently clean.
"""

from __future__ import annotations

import dataclasses
import io
from typing import Any

import pytest
from bt5.codon.tables import FileTableProvider
from bt5.core.context import HostId, Modality
from bt5.core.result import VerificationError
from bt5.core.spec import Enforcement
from bt5.core.types import Interval, reverse_complement
from bt5.design import DesignError, design
from bt5.design.catalog import partition_forbidden, scored_objectives
from bt5.design.runner import _coding_flanks, _context
from bt5.rules.vendors import DEFAULT_SELECTION
from bt5.solver.catalog import build_rule_set, default_services
from bt5.vector import read_genbank
from bt5.vector.backbone import InsertionSite, VectorBackbone, insertion_site_from_interval
from bt5.verify import verify_construct

CODE = FileTableProvider().genetic_code(1)


def rules_and_partition(backbone: VectorBackbone, site: InsertionSite):  # noqa: ANN201
    """The solver rule set plus the design lane's usable/carried partition of its
    forbidden set -- the same objects design() builds internally."""
    ctx = _context(
        modality=Modality.LENTIVIRAL,
        hosts=[HostId.HEK293],
        table_id=1,
        cassette_orientation=site.strand,
    )
    rules = build_rule_set(ctx, default_services(seed=0), vendors=DEFAULT_SELECTION)
    usable, carried = partition_forbidden(rules.forbidden(), backbone, site)
    return rules, usable, carried


class TestEndToEnd:
    def test_every_candidate_is_verified_and_preserves_the_protein(
        self, backbone: VectorBackbone, fast: Any, protein: str
    ) -> None:
        """`optimize()` runs the independent validator on every candidate, so a
        candidate that reached the panel is one `verify_construct` accepted --
        the panel is not a place unverified designs can hide."""
        res = fast(backbone)
        assert res.result.candidates
        assert res.report.candidates == len(res.result.candidates)
        assert res.optimize_result.repair_outcome.clean
        for candidate in res.result.candidates:
            assert CODE.translate(candidate.cds)[:-1] == protein
        # One ORF: the design is single-cassette by construction.
        assert len(res.assembly.construct.translation_units) == 1

    def test_it_never_claims_completeness_or_a_baseline_it_was_not_given(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        res = fast(backbone)
        assert res.report.is_complete is False
        assert res.report.native_baseline_available is False
        assert res.result.native_baseline is None

    def test_the_wild_type_sentence_never_renders_without_a_wild_type(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """The native-sequence line is a claim about a real wild-type CDS; the
        always-on honesty disclaimer stays either way."""
        res = fast(backbone)
        assert "The native sequence is included as a candidate" not in res.rendered
        assert "never a predicted expression level" in res.rendered

    def test_the_genbank_round_trips(self, backbone: VectorBackbone, fast: Any) -> None:
        res = fast(backbone)
        # Through a handle, not the raw string: read_genbank path-tests a string
        # first, which OSErrors on text this long -- a vector-lane gap, not E1's.
        parsed = read_genbank(io.StringIO(res.genbank))
        assert parsed.length == res.assembly.construct.length


class TestSiteHandling:
    def test_forward_flanks_are_the_immediate_neighbours(self, backbone: VectorBackbone) -> None:
        site = insertion_site_from_interval(Interval(720, 780, 1), label="mcs", table_id=1)
        _rules, usable, _carried = rules_and_partition(backbone, site)
        k = max(len(m) for m in usable) - 1
        left, right = _coding_flanks(backbone, site, usable)
        assert left == backbone.sequence[720 - k : 720]
        assert right == backbone.sequence[780 : 780 + k]

    def test_reverse_strand_flanks_are_revcomp_of_the_opposite_side(
        self, backbone: VectorBackbone
    ) -> None:
        """The one place a sign error comes back silently clean: on a reverse
        site the coding-5' neighbour is the revcomp of the DOWNSTREAM backbone."""
        site = insertion_site_from_interval(Interval(720, 780, -1), label="mcs", table_id=1)
        forbidden = ("GCGGCCGC",)  # longest pattern 8 -> k = 7
        left, right = _coding_flanks(backbone, site, forbidden)
        assert left == reverse_complement(backbone.sequence[780:787])
        assert right == reverse_complement(backbone.sequence[713:720])

    def test_a_reverse_strand_design_preserves_the_protein(
        self, backbone: VectorBackbone, fast: Any, protein: str
    ) -> None:
        site = insertion_site_from_interval(Interval(720, 780, -1), label="mcs", table_id=1)
        res = fast(backbone, site=site)
        for candidate in res.result.candidates:
            assert CODE.translate(candidate.cds)[:-1] == protein
            # The insert on the plus strand is the reverse complement of the CDS.
            assert reverse_complement(candidate.cds) in candidate.construct.sequence

    def test_an_origin_spanning_site_verifies(self, backbone: VectorBackbone, fast: Any) -> None:
        """G2 geometry: an insert placed across position 0, with I9 still holding
        (design() would raise VerificationError if the backbone were touched)."""
        site = insertion_site_from_interval(Interval(1990, 2010, 1), label="wrap", table_id=1)
        res = fast(backbone, site=site)
        assert res.result.candidates


class TestCarriedMotif:
    def _site(self) -> InsertionSite:
        return insertion_site_from_interval(Interval(720, 780, 1), label="mcs", table_id=1)

    def test_the_backbone_xbai_is_carried_not_usable(self, backbone: VectorBackbone) -> None:
        _rules, usable, carried = rules_and_partition(backbone, self._site())
        assert "TCTAGA" in carried
        assert "TCTAGA" not in usable

    def test_unpartitioned_the_validator_would_refuse(
        self, backbone: VectorBackbone, fast: Any, protein: str
    ) -> None:
        """The negative: hand the carried motif to the validator and it refuses
        (I6), which is exactly why the partition exists."""
        res = fast(backbone)
        construct = res.assembly.construct
        _rules, usable, _carried = rules_and_partition(backbone, self._site())
        forbidden_with_carried = (*usable, "TCTAGA")
        with pytest.raises(VerificationError) as exc:
            verify_construct(
                construct, protein=protein, table_id=1, forbidden=forbidden_with_carried
            )
        assert exc.value.invariant == "I6"


class TestHonesty:
    def test_the_mandatory_degradations_are_present(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """Set equality moved to `test_no_degradation_arrives_unremarked`, which
        keeps the same protection while letting the CONTENT depend on the
        environment. What is pinned here is the subset that must hold on every
        machine: this run screened nothing, was given no wild-type CDS, and used
        a backbone that carries a forbidden motif no codon can remove."""
        res = fast(backbone)
        degradations = set(res.result.provenance.degradations)
        assert any("biosecurity screening: not_run" in d for d in degradations)
        assert any("no native baseline" in d for d in degradations)
        assert (
            "forbidden motif TCTAGA carried by the backbone, excluded from enforcement"
            in degradations
        )

    def test_the_skeleton_s_own_degradations_are_gone(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """Three sentences the walking skeleton emitted unconditionally. Each is
        now conditional on the thing genuinely being absent, and this run has a
        ranking and a null -- so none of them may appear."""
        degradations = " | ".join(fast(backbone).result.provenance.degradations)
        assert "ranking not computed" not in degradations
        assert "no null distribution and no percentiles" not in degradations

    def test_the_3_5_guard_keeps_hard_rules_out_of_the_scored_set(
        self, backbone: VectorBackbone
    ) -> None:
        site = insertion_site_from_interval(Interval(720, 780, 1), label="mcs", table_id=1)
        rules, _usable, _carried = rules_and_partition(backbone, site)
        slot = rules.ctx.active_slots[0]
        scored = scored_objectives(rules)
        for spec in scored:
            assert not spec.enforcement_for(slot).is_hard
        # d4 is SOFT by ClassVar but HARD_REPAIR here, so it must NOT be scored.
        assert "d4_internal_polya" not in {s.id for s in scored}
        repair_ids = {s.id for s in rules.repair_specs()}
        assert repair_ids >= {"d4_internal_polya", "e2_gc_band"}
        for spec in rules.repair_specs():
            assert spec.enforcement_for(slot) is Enforcement.HARD_REPAIR


class TestDeterminism:
    def test_same_seed_same_panel(self, backbone: VectorBackbone, fast: Any) -> None:
        """Every hash AND every rank. An unseeded null would leave the sequences
        identical and the ORDER of the panel irreproducible, which is the subtler
        half of the same failure -- the user reads the top-ranked one."""
        a = fast(backbone, seed=0)
        b = fast(backbone, seed=0)
        assert [c.design_hash for c in a.result.candidates] == [
            c.design_hash for c in b.result.candidates
        ]
        assert [c.scorecard.total for c in a.result.candidates] == [
            c.scorecard.total for c in b.result.candidates
        ]

    def test_a_different_backbone_gives_a_different_hash(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """design_hash salts on the backbone, so the same insert in two vectors
        does not collide -- the whole reason the CDS-only hash is not enough.

        Compared as SETS: the mutated base sits outside the insertion site, so
        the same CDSs are solved, but the rules evaluate against the whole
        construct and may rank them differently. The claim under test is about
        the hash, not the order.
        """
        mutated_seq = (
            backbone.sequence[:900]
            + ("T" if backbone.sequence[900] != "T" else "A")
            + backbone.sequence[901:]
        )
        other = dataclasses.replace(backbone, sequence=mutated_seq)
        a = fast(backbone)
        b = fast(other)
        assert {c.cds for c in a.result.candidates} == {c.cds for c in b.result.candidates}
        assert not {c.design_hash for c in a.result.candidates} & {
            c.design_hash for c in b.result.candidates
        }


def test_a_protein_without_the_initiator_is_refused(backbone: VectorBackbone) -> None:
    with pytest.raises(DesignError, match="initiator"):
        design(
            backbone=backbone,
            protein="KLVTAAF",
            table_id=1,
            modality=Modality.LENTIVIRAL,
            hosts=[HostId.HEK293],
        )
