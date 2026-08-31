"""`optimize()` once the catalog is wired into it.

Three of these tests cover things that were previously unreachable rather than
merely untested: I9 never ran from `optimize()`, a backbone-carried motif had no
certificate, and "Tier B never ran" was byte-identical to "Tier B found nothing".
"""

from __future__ import annotations

import inspect

import pytest
from bt5.codon.tables import NcbiGeneticCode
from bt5.core.result import (
    InfeasibleConstraints,
    VerificationError,
)
from bt5.core.spec import Breach, LocalizationPolicy, RepairPolicy
from bt5.core.types import Construct, Interval
from bt5.rules.vendors import VendorSelection
from bt5.solver.catalog import build_rule_set, default_services, optimize_with
from bt5.solver.pipeline import optimize
from bt5.solver.repair import NO_RULES, no_rules, repair
from conftest import context, linear_cds, slot, with_backbone

PROTEIN = "M" + "KLIWQRSTVNDEYFPGHACM" * 2
LEFT = "GCTAGCACCATGGTGAGCAAGGGCGAGGAGCTGTTCACC"
RIGHT = "TAAGCGGCCGCTTAATTAAGCTTGCATGCCTGCAGGTCG"


def assembler(protein: str, left: str = LEFT, right: str = RIGHT):
    def build(cds: str) -> Construct:
        return with_backbone(cds, protein, left, right)

    return build


SHORT = "M" + "A" * 6
SHORT_CDS = "ATG" + "GCT" * 6 + "TAA"


def bare(cds: str) -> Construct:
    """No backbone: breach intervals are then CDS coordinates, which is what a
    repair-window test wants to talk about."""
    return linear_cds(cds, SHORT)


def no_breaches(_c: Construct) -> tuple[Breach, ...]:
    """A real finder that happens to find nothing. NOT the same as not running."""
    return ()


class TestNotRunIsNotClean:
    """Hazard 1. `optimize()` fabricated `RepairOutcome(cds, 0, True)` when no
    finder was supplied, which is byte-identical to a clean repair -- so
    forgetting the argument produced a construct that looked proven and had never
    been checked."""

    def test_find_breaches_has_no_default(self) -> None:
        param = inspect.signature(optimize).parameters["find_breaches"]
        assert param.default is inspect.Parameter.empty

    def test_the_opt_out_and_a_clean_repair_are_distinguishable(
        self, code: NcbiGeneticCode
    ) -> None:
        opted_out = optimize(
            PROTEIN, code, assemble=assembler(PROTEIN), find_breaches=NO_RULES, seed=1
        )
        checked = optimize(
            PROTEIN, code, assemble=assembler(PROTEIN), find_breaches=no_breaches, seed=1
        )
        # Identical in every field the old code could report...
        assert opted_out.repair_outcome.iterations == checked.repair_outcome.iterations == 0
        assert opted_out.repair_outcome.converged is checked.repair_outcome.converged is True
        assert opted_out.cds == checked.cds
        # ...and now distinguishable.
        assert opted_out.repair_outcome.ran is False
        assert checked.repair_outcome.ran is True

    def test_the_sentinel_is_a_named_value_not_a_lambda(self) -> None:
        assert NO_RULES is no_rules


class TestTheOracleIsArmed:
    def test_i9_catches_a_touched_backbone_through_optimize(self, code: NcbiGeneticCode) -> None:
        """The highest-value invariant in the oracle, and it never ran from
        `optimize()`: `original_backbone` was simply not forwarded."""
        reference = with_backbone("ATG" + "GCT" * 20 + "TAA", PROTEIN, LEFT, RIGHT)
        tampered = LEFT[:-1] + ("A" if LEFT[-1] != "A" else "T")

        with pytest.raises(VerificationError) as exc:
            optimize(
                PROTEIN,
                code,
                assemble=assembler(PROTEIN, left=tampered),
                find_breaches=NO_RULES,
                original_backbone=reference,
                seed=1,
            )
        assert exc.value.invariant == "I9"

    def test_without_a_reference_the_same_edit_ships(self, code: NcbiGeneticCode) -> None:
        """Not an endorsement -- the point is that arming it is what changed."""
        tampered = LEFT[:-1] + ("A" if LEFT[-1] != "A" else "T")
        res = optimize(
            PROTEIN,
            code,
            assemble=assembler(PROTEIN, left=tampered),
            find_breaches=NO_RULES,
            seed=1,
        )
        assert res.construct.sequence.startswith(tampered)


class TestTheImmutableRegionScreen:
    """A motif the user's own vector carries is a finding about the vector."""

    def test_a_backbone_motif_is_a_named_certificate_not_a_bare_i6_failure(
        self, code: NcbiGeneticCode
    ) -> None:
        with pytest.raises(InfeasibleConstraints) as exc:
            optimize(
                PROTEIN,
                code,
                assemble=assembler(PROTEIN),
                find_breaches=NO_RULES,
                forbidden=["GCGGCCGC"],  # NotI, present in RIGHT
                seed=1,
            )
        cert = exc.value.certificate
        assert cert.proof == "immutable_region"
        assert cert.minimal_conflicting_specs == ("GCGGCCGC",)
        assert cert.interval.length == 8

    def test_a_motif_the_solver_can_reach_is_not_screened_out(self, code: NcbiGeneticCode) -> None:
        """The screen must only fire on sequence no codon can change, or it
        would refuse every design Tier A was about to fix."""
        res = optimize(
            PROTEIN,
            code,
            assemble=assembler(PROTEIN),
            find_breaches=NO_RULES,
            forbidden=["GAATTC"],  # EcoRI, absent from both flanks
            seed=1,
        )
        assert "GAATTC" not in res.cds


class TestPoliciesReachTheLoop:
    def test_a_per_rule_localization_callable_is_consulted_per_breach(
        self, code: NcbiGeneticCode
    ) -> None:
        """One global `policy=` gave every rule the same repair window. The four
        HARD_REPAIR rules shipped today declare three different ones."""
        seen: list[str] = []

        def localization(spec_id: str) -> LocalizationPolicy:
            seen.append(spec_id)
            return LocalizationPolicy.PAIRED_SEGMENTS

        def find(c: Construct) -> tuple[Breach, ...]:
            return (
                Breach(
                    "stubborn",
                    Interval(3, 9),
                    1.0,
                    "never satisfied",
                    fixable_by_codon_choice=True,
                ),
            )

        with pytest.raises(InfeasibleConstraints):
            repair(
                SHORT_CDS,
                SHORT,
                code,
                assemble=bare,
                find_breaches=find,
                policy=localization,
                max_iterations=3,
            )
        assert seen == ["stubborn"] * 3

    def test_a_plain_policy_still_works(self, code: NcbiGeneticCode) -> None:
        """Non-catalog callers pass one enum and are unaffected."""

        def find(_c: Construct) -> tuple[Breach, ...]:
            return ()

        out = repair(
            SHORT_CDS,
            SHORT,
            code,
            assemble=bare,
            find_breaches=find,
            policy=LocalizationPolicy.WHOLE_SCOPE,
        )
        assert out.clean

    def test_the_worst_breach_is_chosen_in_normalised_units(self, code: NcbiGeneticCode) -> None:
        """Ranking by raw `magnitude` compares a nucleotide count against a GC
        fraction, so the rule reporting in the largest units monopolised the
        search and the other was never worked on."""
        targeted: list[str] = []

        def localization(spec_id: str) -> LocalizationPolicy:
            targeted.append(spec_id)
            return LocalizationPolicy.WINDOW_MINUS_1

        def find(_c: Construct) -> tuple[Breach, ...]:
            return (
                Breach("nt", Interval(3, 9), 4.0, "4 nt over", fixable_by_codon_choice=True),
                Breach("gc", Interval(9, 15), 0.05, "5% over", fixable_by_codon_choice=True),
            )

        with pytest.raises(InfeasibleConstraints):
            repair(
                SHORT_CDS,
                SHORT,
                code,
                assemble=bare,
                find_breaches=find,
                policy=localization,
                max_iterations=1,
            )
        # Both are at 1.0 of their own rule's scale, so the tie breaks on
        # position -- not on the eighty-fold difference in raw magnitude.
        assert targeted == ["nt"]

    def test_the_effective_policy_is_recorded(self, code: NcbiGeneticCode) -> None:
        res = optimize(
            PROTEIN,
            code,
            assemble=assembler(PROTEIN),
            find_breaches=no_breaches,
            repair_policy=RepairPolicy.SINGLE_PASS,
            seed=1,
        )
        assert res.repair_outcome.effective_repair_policy is RepairPolicy.SINGLE_PASS


class TestTheCatalogActuallyRuns:
    """The end-to-end claim: rules the repo shipped now decide what is emitted."""

    def test_a_catalog_driven_design_verifies_and_round_trips(self, code: NcbiGeneticCode) -> None:
        rules = build_rule_set(
            context(slot()),
            default_services(seed=5, autoload_fold=False),
            # D1 forbids NotI, which the RIGHT flank of this fixture carries;
            # that path has its own test above.
            include=lambda c: c.id != "d1_restriction_sites",
        )
        assert rules.forbidden()  # E1's homopolymer runs reach Tier A
        res = optimize_with(rules, PROTEIN, code, assemble=assembler(PROTEIN), seed=5)
        assert code.translate(res.cds).rstrip("*") == PROTEIN
        assert res.repair_outcome.ran is True
        lo, hi = rules.oracle_bounds().gc_bounds or (0.0, 1.0)
        gc = (res.cds.count("G") + res.cds.count("C")) / len(res.cds)
        assert lo <= gc <= hi

    def test_the_catalog_refuses_a_backbone_motif_it_forbids(self, code: NcbiGeneticCode) -> None:
        """With D1 in, the NotI site in this fixture's own flank is a finding --
        one a hand-written `forbidden` list would never have surfaced."""
        rules = build_rule_set(context(slot()), default_services(seed=5, autoload_fold=False))
        with pytest.raises(InfeasibleConstraints) as exc:
            optimize_with(rules, PROTEIN, code, assemble=assembler(PROTEIN), seed=5)
        assert exc.value.certificate.proof == "immutable_region"

    def test_the_selected_vendor_reaches_the_validator(self, code: NcbiGeneticCode) -> None:
        """E2 gates on the selection's band and I7 is handed the same numbers, so
        the two cannot enforce different contracts."""
        rules = build_rule_set(
            context(slot()),
            default_services(seed=5, autoload_fold=False),
            vendors=VendorSelection.of("twist_gene_fragment"),
            include=lambda c: c.id != "d1_restriction_sites",
        )
        e2 = next(s for s in rules.specs if s.id == "e2_gc_band")
        assert rules.oracle_bounds().gc_bounds == (e2.gc_min, e2.gc_max)  # type: ignore[attr-defined]
