"""D3: cryptic splice sites, and the fixed point CLAUDE.md 3.6 is about.

The load-bearing test in this file is
`TestFixedPoint::test_a_single_pass_ships_a_donor_the_fixed_point_removes`. Everything
else guards a property that, if it broke, would make that one pass for the wrong reason.
"""

from __future__ import annotations

import pytest
from bt5.codon.tables import FileTableProvider
from bt5.core.context import HostId, Modality
from bt5.core.registry import discover, get
from bt5.core.spec import Enforcement, LocalizationPolicy, RepairPolicy
from bt5.core.types import (
    Construct,
    Interval,
    Segment,
    SegmentKind,
    Topology,
    TranslationUnit,
    reverse_complement,
)
from bt5.rules.catalog.d3_splicing import (
    MAG_COARSE,
    MAG_FLAGGED,
    MAG_STRONG,
    V5_PEPTIDE,
    Splicing,
    _matches,
)
from bt5.solver.repair import RulePolicy, repair
from conftest import context, slot

discover()

#: Local helpers, NOT added to conftest.py: that file is shared with the other rules
#: session and an edit there is the one collision this split cannot absorb
#: (docs/buildout/README.md:62-72).


@pytest.fixture(scope="module")
def code() -> object:
    """NCBI table 1. Explicit and never defaulted (CLAUDE.md 3.1)."""
    return FileTableProvider().genetic_code(1)


def lenti() -> object:
    """A context whose slot makes this rule HARD_REPAIR (brief.md:223)."""
    return context(slot(modality=Modality.LENTIVIRAL))


def tagged(seq: str, protein: str, *, cds_start: int = 0) -> Construct:
    """A linear construct carrying a translation unit, so V5 can be found on the
    protein rather than on any one nucleotide encoding of it."""
    codons = tuple(Interval(cds_start + 3 * i, cds_start + 3 * i + 3) for i in range(len(protein)))
    return Construct(
        sequence=seq,
        topology=Topology.LINEAR,
        segments=(Segment(Interval(0, len(seq)), SegmentKind.DESIGNABLE_CDS, "cds"),),
        translation_units=(TranslationUnit(1, codons, protein, False, False),),
    )


def whole(seq: str, topology: Topology) -> Construct:
    return Construct(
        sequence=seq,
        topology=topology,
        segments=(Segment(Interval(0, len(seq)), SegmentKind.DESIGNABLE_CDS, "cds"),),
    )


def donors(rule: Splicing, c: Construct, ctx: object | None = None) -> list:
    ev = rule.evaluate(c, ctx or lenti(), None)  # type: ignore[arg-type]
    return [b for b in ev.breaches if b.detail.get("kind") == "donor"]


def strong(rule: Splicing, c: Construct, ctx: object | None = None) -> list:
    return [b for b in donors(rule, c, ctx) if b.magnitude >= MAG_STRONG]


class TestIupacMatch:
    def test_degenerate_codes(self) -> None:
        assert _matches("AGGTAAG", "ANGTRAG")
        assert _matches("ACGTGAG", "ANGTRAG")
        assert not _matches("CGGTAAG", "ANGTRAG"), "the -2 exonic A is required"
        assert not _matches("AGGTCAG", "ANGTRAG"), "R is A or G, not C"

    def test_length_mismatch_is_never_a_match(self) -> None:
        assert not _matches("AGGTAA", "ANGTRAG")


class TestDonorTiers:
    """brief.md:102 gives three shapes and brief.md:90 grades them by length."""

    def test_the_literal_blacklist_is_strong(self) -> None:
        for motif in ("GGTAAG", "GGTGAG"):
            hits = strong(Splicing(), whole("ATG" + motif + "CCCTTTAAA", Topology.LINEAR))
            assert len(hits) == 1, motif
            assert motif in hits[0].message

    def test_the_seven_mer_consensus_is_strong(self) -> None:
        """`AN|GT(A/G)AG` is 7 nt, and brief.md:90 puts >=7 nt in the hard tier."""
        hits = strong(Splicing(), whole("ATGAC" + "ACGTGAG" + "CCCTTTAAA", Topology.LINEAR))
        assert len(hits) == 1
        assert "consensus" in hits[0].message

    def test_the_coarse_five_mer_is_soft_only(self) -> None:
        """GTNNG is a 5-mer: brief.md:90 says <=5-mer is SOFT ONLY, so it must never
        reach the hard tier no matter how the rest of the rule changes."""
        c = whole("ATGCC" + "GTCCG" + "CCCTTTAAA", Topology.LINEAR)
        hits = donors(Splicing(report_coarse=True), c)
        assert [b.magnitude for b in hits] == [MAG_COARSE]

    def test_a_site_matching_two_shapes_is_one_breach(self) -> None:
        """GGTAAG matches the literal list AND the coarse screen; reporting both
        would double-count one site and skew the magnitude the solver reads."""
        assert len(donors(Splicing(), whole("ATG" + "GGTAAG" + "CCCTTTAAA", Topology.LINEAR))) == 1

    def test_the_coarse_screen_ships_off(self) -> None:
        """P(GTNNG) = 1/64 per position is ~78 hits on a 5 kb plasmid of random
        sequence. A screen that fires 78 times by chance reports noise, and
        brief.md:90 already grades a <=5-mer soft-only."""
        c = whole("ATGCC" + "GTCCG" + "CCCTTTAAA", Topology.LINEAR)
        assert donors(Splicing(), c) == []
        assert donors(Splicing(report_coarse=True), c) != []

    def test_a_clean_construct_passes_with_no_donors(self) -> None:
        assert donors(Splicing(), whole("ATGCCCTTTAAACCCTTTAAA", Topology.LINEAR)) == []


class TestPassesIsNotNotBreaches:
    def test_a_coarse_only_construct_still_passes(self) -> None:
        """solver/catalog.py:158-170: the solver chases only breaches from a rule that
        says it FAILED. A rule that failed on its own warn band sets repair chasing a
        threshold that was never crossed, and the search stagnates on a design the
        catalog accepts."""
        ev = Splicing(report_coarse=True).evaluate(
            whole("ATGCC" + "GTCCG" + "CCCTTTAAA", Topology.LINEAR), lenti(), None
        )  # type: ignore[arg-type]
        assert ev.passes
        assert any(b.magnitude == MAG_COARSE for b in ev.breaches)

    def test_a_strong_donor_fails(self) -> None:
        ev = Splicing().evaluate(
            whole("ATG" + "GGTAAG" + "CCCTTTAAA", Topology.LINEAR), lenti(), None
        )  # type: ignore[arg-type]
        assert not ev.passes


class TestStrand:
    """The property that makes this rule not a lattice rule. If it ever became
    reverse-complement symmetric, `lattice_terms() -> None` would be unjustified and
    the rule would be refusing designs over sites nothing transcribes."""

    SEQ = "ATG" + "GGTAAG" + "AAACCCTTTGGG" * 3

    def test_a_donor_is_found_on_the_forward_strand(self) -> None:
        assert len(strong(Splicing(), whole(self.SEQ, Topology.LINEAR))) == 1

    def test_the_reverse_complement_is_clean_when_read_forward(self) -> None:
        rc = whole(reverse_complement(self.SEQ), Topology.LINEAR)
        assert strong(Splicing(), rc) == [], "a donor's complement is not a donor"

    def test_the_same_sequence_is_dirty_when_the_cassette_is_reversed(self) -> None:
        """brief.md:244: the lentiviral genome and the target-cell mRNA are opposite
        strands in a reverse-oriented cassette. Without the orientation the splice
        analysis is exactly backwards -- these two tests are that claim."""
        rc = whole(reverse_complement(self.SEQ), Topology.LINEAR)
        ctx = context(slot(modality=Modality.LENTIVIRAL), cassette_orientation=-1)
        hits = strong(Splicing(), rc, ctx)
        assert len(hits) == 1
        assert hits[0].interval.strand == -1
        assert hits[0].interval.start == len(self.SEQ) - 1 - 9, "mapped back to forward coords"

    def test_it_reads_the_composed_orientation_not_the_raw_slot(self) -> None:
        """Two -1s compose to +1 (core/spec.py:289); a rule reading the slot field
        directly would get this backwards."""
        rc = whole(reverse_complement(self.SEQ), Topology.LINEAR)
        ctx = context(
            slot(modality=Modality.LENTIVIRAL, strand_of_interest=-1), cassette_orientation=-1
        )
        assert strong(Splicing(), rc, ctx) == []


class TestCircular:
    """CLAUDE.md 3.3: rules take a Construct, never a string. This is why."""

    SEQ = "AAG" + "CCCTTTAAACCC" * 3 + "GGT"  # ...GGT | AAG... across the origin

    def test_a_donor_spanning_the_origin_is_caught(self) -> None:
        hits = strong(Splicing(), whole(self.SEQ, Topology.CIRCULAR))
        assert len(hits) == 1
        iv = hits[0].interval
        assert iv.start == 37
        assert iv.end == 46
        assert iv.end > len(self.SEQ), "a wrapping interval, not a clamped one"

    def test_the_same_sequence_linear_is_clean(self) -> None:
        assert strong(Splicing(), whole(self.SEQ, Topology.LINEAR)) == []

    def test_a_linear_construct_reads_no_context_off_its_own_ends(self) -> None:
        """A GT in the first three bases has no exonic context; splicing it out of the
        far end of the string would invent a site from two ends that are not adjacent."""
        assert donors(Splicing(), whole("GTAAGCCCTTTAAACCC", Topology.LINEAR)) == []


class TestV5:
    """brief.md:105, the only evidence-A element of D3."""

    #: The standard encoding brief.md:105 names, carrying G|GTAAG at the Gly-Lys join.
    STANDARD = "GGTAAGCCTATCCCTAACCCTCTCCTCGGTCTCGATTCTACG"

    def test_the_standard_encoding_carries_a_donor_and_is_escalated(self) -> None:
        c = tagged("ATG" + self.STANDARD + "AAACCCTTT", V5_PEPTIDE, cds_start=3)
        hits = [b for b in donors(Splicing(), c) if b.detail.get("v5_tag") == "yes"]
        assert hits, "the 17/17 case must be found"
        assert all(b.magnitude == MAG_STRONG for b in hits)
        assert "17/17" in hits[0].message

    def test_the_tag_is_matched_on_the_protein_not_the_nucleotides(self) -> None:
        """The liability IS the standard encoding, so a nucleotide match would stop
        matching the moment repair recoded it -- the rule would report success exactly
        when it stopped working. A different encoding of the same peptide must still be
        recognised as V5."""
        recoded = "GGAAAACCAATCCCAAACCCACTCCTCGGACTCGATTCTACG"
        c = tagged("ATG" + recoded + "AAACCCTTT", V5_PEPTIDE, cds_start=3)
        rule = Splicing()
        assert rule._v5_spans(c), "the tag is still located after recoding"
        assert strong(rule, c) == [], "and this encoding carries no donor to escalate"

    def test_a_construct_with_no_translation_unit_has_no_v5_span(self) -> None:
        assert Splicing()._v5_spans(whole("ATGGGTAAGCCC", Topology.LINEAR)) == []


class TestAcceptors:
    """brief.md:103 pairs two context features. Either alone is not a finding."""

    #: Built to the geometry rather than by eye: the AG lands at 40, the branch point
    #: YTNAY at 16-20 (within 18-40 nt upstream) and a 10-nt pure-pyrimidine tract at
    #: 24-33 (within 5-40 nt upstream).
    BOTH = "AAACCCAAACCCAAAC" + "CTAAC" + "AAA" + "TCTCTCTCTC" + "CCAAAC" + "AG" + "CCCAAA"
    #: The same, with the branch point destroyed and the tract left intact.
    TRACT_ONLY = "AAACCCAAACCCAAAC" + "CAAAC" + "AAA" + "TCTCTCTCTC" + "CCAAAC" + "AG" + "CCCAAA"
    #: The same, with the tract destroyed and the branch point left intact.
    BRANCH_ONLY = "AAACCCAAACCCAAAC" + "CTAAC" + "AAA" + "AAAAAAAAAA" + "CCAAAC" + "AG" + "CCCAAA"

    def test_both_features_together_are_reported(self) -> None:
        hits = [
            b
            for b in Splicing(scan_acceptors=True)
            .evaluate(whole(self.BOTH, Topology.LINEAR), lenti(), None)
            .breaches
            if b.detail.get("kind") == "acceptor"
        ]  # type: ignore[arg-type]
        assert hits
        assert all(b.magnitude == MAG_FLAGGED for b in hits)

    def test_an_acceptor_is_never_hard(self) -> None:
        """brief.md:103 gives a flag threshold for acceptors and no hard one. The rule
        must not invent the cutoff the brief declined to state."""
        ev = Splicing(scan_acceptors=True).evaluate(
            whole(self.BOTH, Topology.LINEAR), lenti(), None
        )  # type: ignore[arg-type]
        assert all(
            b.magnitude < MAG_STRONG for b in ev.breaches if b.detail.get("kind") == "acceptor"
        )

    def test_a_bare_ag_is_not_a_finding(self) -> None:
        c = whole("AAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", Topology.LINEAR)
        ev = Splicing(scan_acceptors=True).evaluate(c, lenti(), None)  # type: ignore[arg-type]
        assert [b for b in ev.breaches if b.detail.get("kind") == "acceptor"] == []

    @pytest.mark.parametrize("seq_name", ["TRACT_ONLY", "BRANCH_ONLY"])
    def test_either_feature_alone_is_not_a_finding(self, seq_name: str) -> None:
        """BOTH, not either. A lone AG occurs every ~16 bases, and brief.md:103 pairs
        the two context features for exactly that reason -- reporting on one would bury
        every real finding under noise."""
        ev = Splicing(scan_acceptors=True).evaluate(
            whole(getattr(self, seq_name), Topology.LINEAR), lenti(), None
        )  # type: ignore[arg-type]
        assert [b for b in ev.breaches if b.detail.get("kind") == "acceptor"] == []

    def test_the_acceptor_scan_ships_off(self) -> None:
        """This rule's enforcement class is set by its donor half, so an acceptor
        breach emitted here would inherit HARD_REPAIR and set the solver chasing a
        site brief.md:103 only flags. Reachable, but not on by default."""
        ev = Splicing().evaluate(whole(self.BOTH, Topology.LINEAR), lenti(), None)  # type: ignore[arg-type]
        assert [b for b in ev.breaches if b.detail.get("kind") == "acceptor"] == []


class TestScoredPathIsHonest:
    def test_it_says_that_maxentscan_did_not_run(self) -> None:
        """No MaxEntScan model ships (brief.md:288, redistribution ambiguous). Reporting
        only the motif scan without saying so reads as a clean scan of everything the
        rule claims to cover."""
        ev = Splicing().evaluate(whole("ATGCCCTTTAAACCC", Topology.LINEAR), lenti(), None)  # type: ignore[arg-type]
        notices = [b for b in ev.breaches if "unavailable_reason" in b.detail]
        assert len(notices) == 1
        assert notices[0].magnitude == 0.0
        assert not notices[0].fixable_by_codon_choice, "no codon supplies a missing model"
        assert ev.passes, "an uncomputed objective is not a breach of the construct"

    def test_a_gated_off_slot_makes_no_claim_at_all(self) -> None:
        ev = Splicing().evaluate(
            whole("ATGGGTAAGCCCTTTAAA", Topology.LINEAR),
            context(slot(modality=Modality.IVT_MRNA)),
            None,  # type: ignore[arg-type]
        )
        assert ev.breaches == ()
        assert ev.n_evaluated == 0


class TestGating:
    @pytest.mark.parametrize(
        ("modality", "host", "table", "applies"),
        [
            (Modality.LENTIVIRAL, HostId.HEK293, 1, True),
            (Modality.GENOME_INTEGRATED, HostId.HEK293, 1, True),
            (Modality.AAV, HostId.HEK293, 1, True),
            (Modality.PLASMID_TRANSIENT, HostId.HEK293, 1, True),
            (Modality.IVT_MRNA, HostId.HEK293, 1, False),
            (Modality.BACTERIAL_EXPRESSION, HostId.E_COLI_K12, 11, False),
            (Modality.PLASMID_STABLE, HostId.S_CEREVISIAE, 1, False),
        ],
    )
    def test_eukaryotic_pol_ii_contexts_only(
        self, modality: Modality, host: HostId, table: int, applies: bool
    ) -> None:
        """brief.md:101 and brief.md:208."""
        assert Splicing().gate(slot(modality=modality, host=host, table=table)) is applies

    @pytest.mark.parametrize(
        ("modality", "enforcement"),
        [
            (Modality.LENTIVIRAL, Enforcement.HARD_REPAIR),
            (Modality.GENOME_INTEGRATED, Enforcement.HARD_REPAIR),
            (Modality.AAV, Enforcement.SOFT),
            (Modality.PLASMID_TRANSIENT, Enforcement.SOFT),
            (Modality.PLASMID_STABLE, Enforcement.SOFT),
        ],
    )
    def test_enforcement_follows_the_brief_row(
        self, modality: Modality, enforcement: Enforcement
    ) -> None:
        """brief.md:223: warn | warn | HARD (titer + safety) | warn | n/a |
        HARD (fusion transcripts)."""
        assert Splicing().enforcement_for(slot(modality=modality)) is enforcement


class TestSpecShape:
    def test_it_declares_the_fixed_point_policy(self) -> None:
        """CLAUDE.md 3.6. This rule is the first in the catalog to declare it."""
        assert Splicing.repair is RepairPolicy.FIXED_POINT

    def test_it_is_not_a_lattice_rule(self) -> None:
        """`forbidden` is closed under reverse complement by the solver, which would
        forbid CTTACC because GGTAAG is a donor -- refusing designs over a site that
        cannot fire. See d4_internal_polya.py:10-15 for the same argument."""
        assert Splicing().lattice_terms(lenti()) is None  # type: ignore[arg-type]

    def test_the_repair_window_reaches_the_whole_donor(self) -> None:
        """solver/catalog.py:241 hard-codes motif_len to 6 for every rule, so
        MOTIF_LEN_MINUS_1 would widen by 5 when a 9-mer donor needs 8. catalog.py:236
        reads `window` off the spec, and it must stay an int or catalog.py silently
        falls back to 50."""
        rule = Splicing()
        assert rule.localization is LocalizationPolicy.WINDOW_MINUS_1
        assert isinstance(rule.window, int)
        assert rule.window == 9

    def test_it_is_registered_under_its_brief_row(self) -> None:
        assert get("d3_splicing").brief_ref == "2.D3"

    def test_a_soft_floor_carries_its_weight_provenance(self) -> None:
        assert Splicing.enforcement is Enforcement.SOFT
        assert Splicing.weight_provenance.strip()


class TestFixedPoint:
    """CLAUDE.md 3.6, made executable.

    `MKGKGKPF` puts two GGTAAG donors inside one repair window. `repair()` takes the
    FIRST accepting candidate rather than the best (solver/repair.py:545-558), so the
    attempt on the first donor is spent on a candidate that clears the SECOND one --
    an improvement that leaves the target standing. `SINGLE_PASS` then retires the
    target unconditionally, cleared or not (repair.py:569-571), and the construct ships
    with a cryptic donor. `FIXED_POINT` leaves it eligible and comes back for it.
    """

    CDS = "ATGAAAGGTAAGGGTAAGCCCTTT"
    PROTEIN = "MKGKGKPF"

    def _assembler(self):
        def assemble(cds: str) -> Construct:
            return tagged(cds, self.PROTEIN)

        return assemble

    def _finder(self):
        """Faithful to `RuleSet.findings` (solver/catalog.py:181): a rule that PASSES
        contributes nothing. Built here rather than with the real RuleSet so the test
        isolates this rule instead of the whole catalog."""
        rule = Splicing(scan_acceptors=False)

        def find(c: Construct) -> tuple:
            ev = rule.evaluate(c, lenti(), None)  # type: ignore[arg-type]
            return () if ev.passes else tuple(ev.breaches)

        return find

    def _run(self, code: object, policy: RepairPolicy):
        return repair(
            self.CDS,
            self.PROTEIN,
            code,  # type: ignore[arg-type]
            assemble=self._assembler(),
            find_breaches=self._finder(),
            policies={
                "d3_splicing": RulePolicy(LocalizationPolicy.WINDOW_MINUS_1, policy, 9, 6, 0)
            },
            seed=0,
            max_iterations=200,
            raise_on_infeasible=False,
        )

    def test_the_fixture_starts_with_two_strong_donors_in_one_window(self, code: object) -> None:
        hits = strong(Splicing(), tagged(self.CDS, self.PROTEIN))
        assert [(b.interval.start, b.interval.end) for b in hits] == [(4, 13), (10, 19)]

    def test_a_single_pass_ships_a_donor_the_fixed_point_removes(self, code: object) -> None:
        """The exact failure CLAUDE.md 3.6 describes, as a differential."""
        single = self._run(code, RepairPolicy.SINGLE_PASS)
        fixed = self._run(code, RepairPolicy.FIXED_POINT)

        assert not single.clean, "single-pass must not reach a clean construct here"
        assert [b.magnitude for b in single.remaining] == [MAG_STRONG]
        assert single.remaining[0].interval.start == 4, "the target it retired uncleared"

        assert fixed.clean, "the fixed point must clear what one pass left behind"
        assert fixed.remaining == ()

    def test_both_policies_preserve_the_protein(self, code: object) -> None:
        """A repair that changed the protein would be a far worse bug than the donor,
        and `remaining == ()` alone would not catch it."""
        for policy in (RepairPolicy.SINGLE_PASS, RepairPolicy.FIXED_POINT):
            out = self._run(code, policy)
            assert code.translate(out.cds) == self.PROTEIN, policy  # type: ignore[attr-defined]

    def test_the_fixed_point_run_is_deterministic(self, code: object) -> None:
        """Seeded explicitly (CLAUDE.md 3.7); two runs must agree base for base."""
        assert self._run(code, RepairPolicy.FIXED_POINT).cds == (
            self._run(code, RepairPolicy.FIXED_POINT).cds
        )
