"""The walking skeleton, end to end and at its seams.

`design()` is the first thing in BT5 that runs every lane against a real
construct, so these tests are as much about what it REFUSES to claim -- no
baseline, no complete report, no score -- as about the sequence it emits.
"""

from __future__ import annotations

import dataclasses
import io
from pathlib import Path

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
from bt5.structure.vienna import degradation_reason
from bt5.vector import read_genbank
from bt5.vector.backbone import InsertionSite, VectorBackbone, insertion_site_from_interval
from bt5.verify import verify_construct

ROOT = Path(__file__).resolve().parents[4]
MCS_PATH = ROOT / "tests" / "data" / "backbones" / "synthetic_mcs_ef1a.gb"

#: A 140-residue protein, initiator first. Long enough to fill the fragment.
PROTEIN = (
    "MKLVTAAFERSKSVQNYVVSTKDSPLYYLRKWVRSGYKFDCEEVGLREHQGPAATYTPTQAIWRLTLPSPLL"
    "NVDVWQNSCKSLQHTASWKKHRFGLFTLVISPLIRLGEVASLCGLCEHTATSEVKVCPIDCLQSPTSF"
)
CODE = FileTableProvider().genetic_code(1)


@pytest.fixture
def backbone() -> VectorBackbone:
    return read_genbank(MCS_PATH)


def run(bb: VectorBackbone, **kw: object) -> object:
    return design(
        backbone=bb,
        protein=PROTEIN,
        table_id=1,
        modality=Modality.LENTIVIRAL,
        hosts=[HostId.HEK293],
        **kw,  # type: ignore[arg-type]
    )


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
    def test_one_verified_candidate_protein_preserved(self, backbone: VectorBackbone) -> None:
        res = run(backbone)
        assert len(res.result.candidates) == 1  # type: ignore[attr-defined]
        assert res.report.candidates == 1  # type: ignore[attr-defined]
        assert res.optimize_result.repair_outcome.clean  # type: ignore[attr-defined]
        cds = res.result.candidates[0].cds  # type: ignore[attr-defined]
        assert CODE.translate(cds)[:-1] == PROTEIN
        # One ORF: the skeleton is single-cassette by construction.
        assert len(res.assembly.construct.translation_units) == 1  # type: ignore[attr-defined]

    def test_it_never_claims_completeness_or_a_baseline(self, backbone: VectorBackbone) -> None:
        res = run(backbone)
        assert res.report.is_complete is False  # type: ignore[attr-defined]
        assert res.report.native_baseline_available is False  # type: ignore[attr-defined]
        assert res.result.native_baseline is None  # type: ignore[attr-defined]

    def test_the_wild_type_sentence_never_renders(self, backbone: VectorBackbone) -> None:
        """The native-sequence line is a claim about a real wild-type CDS, which
        the skeleton has none of; the always-on honesty disclaimer stays."""
        res = run(backbone)
        assert "The native sequence is included as a candidate" not in res.rendered  # type: ignore[attr-defined]
        assert "never a predicted expression level" in res.rendered  # type: ignore[attr-defined]

    def test_every_objective_is_reported_unavailable(self, backbone: VectorBackbone) -> None:
        res = run(backbone)
        scored = res.result.candidates[0].scorecard  # type: ignore[attr-defined]
        assert scored.available == ()  # nothing scored
        assert scored.unavailable  # but the objectives are named, not dropped
        assert scored.total == 0.0

    def test_the_genbank_round_trips(self, backbone: VectorBackbone) -> None:
        res = run(backbone)
        # Through a handle, not the raw string: read_genbank path-tests a string
        # first, which OSErrors on text this long -- a vector-lane gap, not E1's.
        parsed = read_genbank(io.StringIO(res.genbank))  # type: ignore[attr-defined]
        assert parsed.length == res.assembly.construct.length  # type: ignore[attr-defined]


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

    def test_a_reverse_strand_design_preserves_the_protein(self, backbone: VectorBackbone) -> None:
        site = insertion_site_from_interval(Interval(720, 780, -1), label="mcs", table_id=1)
        res = run(backbone, site=site)
        cds = res.result.candidates[0].cds  # type: ignore[attr-defined]
        assert CODE.translate(cds)[:-1] == PROTEIN
        # The insert on the plus strand is the reverse complement of the CDS.
        assembled = res.assembly.construct  # type: ignore[attr-defined]
        assert reverse_complement(cds) in assembled.sequence

    def test_an_origin_spanning_site_verifies(self, backbone: VectorBackbone) -> None:
        """G2 geometry: an insert placed across position 0, with I9 still holding
        (design() would raise VerificationError if the backbone were touched)."""
        site = insertion_site_from_interval(Interval(1990, 2010, 1), label="wrap", table_id=1)
        res = run(backbone, site=site)
        assert len(res.result.candidates) == 1  # type: ignore[attr-defined]
        assert CODE.translate(res.result.candidates[0].cds)[:-1] == PROTEIN  # type: ignore[attr-defined]


class TestCarriedMotif:
    def _site(self) -> InsertionSite:
        return insertion_site_from_interval(Interval(720, 780, 1), label="mcs", table_id=1)

    def test_the_backbone_xbai_is_carried_not_usable(self, backbone: VectorBackbone) -> None:
        _rules, usable, carried = rules_and_partition(backbone, self._site())
        assert "TCTAGA" in carried
        assert "TCTAGA" not in usable

    def test_unpartitioned_the_validator_would_refuse(self, backbone: VectorBackbone) -> None:
        """The negative: hand the carried motif to the validator and it refuses
        (I6), which is exactly why the partition exists."""
        res = run(backbone)
        construct = res.assembly.construct  # type: ignore[attr-defined]
        _rules, usable, _carried = rules_and_partition(backbone, self._site())
        forbidden_with_carried = (*usable, "TCTAGA")
        with pytest.raises(VerificationError) as exc:
            verify_construct(
                construct, protein=PROTEIN, table_id=1, forbidden=forbidden_with_carried
            )
        assert exc.value.invariant == "I6"


class TestHonesty:
    def test_the_declared_degradations(self, backbone: VectorBackbone) -> None:
        """A set equality so a NEW silent degradation source fails the test. The
        fold degradation is conditioned on whether ViennaRNA is installed."""
        res = run(backbone)
        expected = {
            "ranking not computed: no null distribution and no percentiles",
            "protein-level biosecurity screening did not run",
            "single candidate only: no gallery",
            "forbidden motif TCTAGA carried by the backbone, excluded from enforcement",
        }
        fold = degradation_reason()
        if fold is not None:
            expected.add(fold)
        assert set(res.result.provenance.degradations) == expected  # type: ignore[attr-defined]

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
    def test_same_seed_same_hash(self, backbone: VectorBackbone) -> None:
        a = run(backbone, seed=0)
        b = run(backbone, seed=0)
        assert a.result.candidates[0].design_hash == b.result.candidates[0].design_hash  # type: ignore[attr-defined]

    def test_a_different_backbone_gives_a_different_hash(self, backbone: VectorBackbone) -> None:
        """design_hash salts on the backbone, so the same insert in two vectors
        does not collide -- the whole reason the CDS-only hash is not enough."""
        mutated_seq = (
            backbone.sequence[:900]
            + ("T" if backbone.sequence[900] != "T" else "A")
            + backbone.sequence[901:]
        )
        other = dataclasses.replace(backbone, sequence=mutated_seq)
        a = run(backbone)
        b = run(other)
        assert a.result.candidates[0].cds == b.result.candidates[0].cds  # same CDS ...
        assert (  # ... different backbone, so a different hash
            a.result.candidates[0].design_hash != b.result.candidates[0].design_hash  # type: ignore[attr-defined]
        )


def test_a_protein_without_the_initiator_is_refused(backbone: VectorBackbone) -> None:
    with pytest.raises(DesignError, match="initiator"):
        design(
            backbone=backbone,
            protein="KLVTAAF",
            table_id=1,
            modality=Modality.LENTIVIRAL,
            hosts=[HostId.HEK293],
        )
