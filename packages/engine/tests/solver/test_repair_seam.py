"""Seam behaviours of `repair()` that the fake assembler in test_repair.py cannot
show: cross-rule fairness, per-rule policy dispatch, the junction guard, the
candidate budget, and honest exits.

These use small hand-built fake finders rather than the real catalog -- each
isolates one mechanism -- so they are fast and deterministic. The real-catalog
proof lives in test_seam_against_the_catalog.py.
"""

from __future__ import annotations

import pytest
from bt5.codon.tables import FileTableProvider
from bt5.core.result import InfeasibleConstraints
from bt5.core.spec import Breach, LocalizationPolicy, RepairPolicy
from bt5.core.types import (
    Construct,
    Interval,
    Segment,
    SegmentKind,
    Topology,
    TranslationUnit,
)
from bt5.solver.repair import RulePolicy, repair


@pytest.fixture(scope="module")
def code() -> object:
    return FileTableProvider().genetic_code(11)


def all_cds_assembler(protein: str):
    """A construct that is one editable CDS from base 0, so every breach a finder
    plants is fixable and the only variable under test is which one is worked."""

    def assemble(cds: str) -> Construct:
        return Construct(
            sequence=cds,
            topology=Topology.LINEAR,
            segments=(Segment(Interval(0, len(cds)), SegmentKind.DESIGNABLE_CDS, "cds"),),
            translation_units=(
                TranslationUnit(
                    11, tuple(Interval(i, i + 3) for i in range(0, len(cds), 3)), protein, True
                ),
            ),
        )

    return assemble


class TestStarvationFreeSelection:
    """Worst-first-by-magnitude starves a small-magnitude rule. Breach count is
    the currency and rules are worked round-robin, so a GC deviation of 0.05 is
    targeted even beside a repeat of magnitude 4.0."""

    def test_a_tiny_magnitude_breach_is_targeted_beside_a_huge_one(self, code: object) -> None:
        # Codon 2 is Leu: CTG (GC 2/3) trips the e2 finder, CTA/TTA clear it.
        # Codons 10-11 are Trp (TGG, one synonym), so the f1 window is trivial and
        # the finder reports f1 UNCONDITIONALLY -- it can never be cleared.
        protein = "M" + "L" + "A" * 7 + "W" * 4 + "A" * 2  # 15 residues
        start = "ATG" + "CTG" + "GCC" * 7 + "TGG" * 4 + "GCC" * 2
        assert len(start) == 45

        def finder(c: Construct) -> tuple[Breach, ...]:
            out: list[Breach] = []
            gc = c.sequence[3:6].count("G") + c.sequence[3:6].count("C")
            if gc >= 2:  # e2, tiny magnitude, clearable by recoding codon 2
                out.append(
                    Breach("e2_gc_band", Interval(3, 6), 0.05, "gc", fixable_by_codon_choice=True)
                )
            # f1, huge magnitude, present no matter what codons are chosen.
            out.append(
                Breach(
                    "f1_direct_repeats", Interval(30, 36), 4.0, "rep", fixable_by_codon_choice=True
                )
            )
            return tuple(out)

        # f1 is unfixable in practice, so the pass must end infeasible -- but on
        # the tiny e2 breach having been cleared first, not last.
        with pytest.raises(InfeasibleConstraints) as exc:
            repair(
                start,
                protein,
                code,  # type: ignore[arg-type]
                assemble=all_cds_assembler(protein),
                find_breaches=finder,
                window=3,
                seed=0,
                max_iterations=200,
            )
        specs = exc.value.certificate.minimal_conflicting_specs
        assert "f1_direct_repeats" in specs
        assert "e2_gc_band" not in specs, (
            "the 0.05-magnitude e2 breach was targeted and cleared despite the "
            "4.0-magnitude f1 breach beside it -- worst-first would never reach it"
        )


class TestJunctionGuard:
    """A recoded block can create a forbidden motif that runs into the immutable
    suffix, and a motif of length L needs L-1 suffix bases to do it. The guard
    must scan by the longest PATTERN, never by a rule's motif_len."""

    def test_a_long_forbidden_motif_at_the_junction_is_never_emitted(self, code: object) -> None:
        # NotI, GCGGCCGC, is 8 nt and palindromic. The only codon that clears the
        # finder's breach is GCG for the last residue -- and GCG followed by the
        # immutable right flank "CGGCCGC" spells GCGGCCGC across the junction,
        # using one block base and seven flank bases. A guard keyed on motif_len=2
        # would scan only two flank bases and wave it through; the fix keys the
        # guard on the pattern length (8), scans seven, and rejects it. With no
        # other way to clear the breach, the pass ends infeasible rather than
        # shipping the site Tier A guaranteed away.
        protein = "MA"
        start = "ATGGCA"  # last codon GCA, not yet GCG
        assembler = all_cds_assembler(protein)

        def finder(c: Construct) -> tuple[Breach, ...]:
            if c.sequence[-3:] == "GCG":
                return ()
            return (Breach("r", Interval(3, 6), 1.0, "x", fixable_by_codon_choice=True),)

        # Sanity: GCG really is the escape the guard must block.
        assert finder(assembler("ATGGCG")) == ()

        with pytest.raises(InfeasibleConstraints):
            repair(
                start,
                protein,
                code,  # type: ignore[arg-type]
                assemble=assembler,
                find_breaches=finder,
                forbidden=["GCGGCCGC"],
                right_flank="CGGCCGC",
                window=3,
                motif_len=2,  # deliberately shorter than the pattern
                seed=0,
                max_iterations=20,
            )


class TestPolicyDispatch:
    """A per-rule RulePolicy overrides the scalar fallbacks. Window reach is the
    sharpest observable: the fixing codon sits outside the fallback window and
    inside the per-rule one."""

    def _finder_needing_a_far_codon(self):
        # The breach is at codon 1, but it only clears when codon 3 becomes GGG.
        def finder(c: Construct) -> tuple[Breach, ...]:
            if c.sequence[9:12] == "GGG":
                return ()
            return (Breach("r", Interval(3, 6), 1.0, "x", fixable_by_codon_choice=True),)

        return finder

    def test_the_fallback_window_cannot_reach_the_fixing_codon(self, code: object) -> None:
        protein = "MGGG"
        start = "ATGGGAGGAGGA"  # codon 3 is GGA, not GGG
        with pytest.raises(InfeasibleConstraints):
            repair(
                start,
                protein,
                code,  # type: ignore[arg-type]
                assemble=all_cds_assembler(protein),
                find_breaches=self._finder_needing_a_far_codon(),
                window=3,  # reaches codons 0-2 only
                seed=0,
                max_iterations=50,
            )

    def test_a_per_rule_window_reaches_it(self, code: object) -> None:
        protein = "MGGG"
        start = "ATGGGAGGAGGA"
        out = repair(
            start,
            protein,
            code,  # type: ignore[arg-type]
            assemble=all_cds_assembler(protein),
            find_breaches=self._finder_needing_a_far_codon(),
            window=3,  # same small fallback ...
            policies={  # ... overridden for this rule with a window that reaches codon 3
                "r": RulePolicy(
                    LocalizationPolicy.WINDOW_MINUS_1, RepairPolicy.FIXED_POINT, 9, 6, 0
                )
            },
            seed=0,
            max_iterations=50,
        )
        assert out.clean
        # No terminal stop in this fixture, so translate without stripping one.
        assert code.translate(out.cds) == protein  # type: ignore[attr-defined]


class TestCandidateBudget:
    """The candidate budget bounds find_breaches calls PER ITERATION on both
    branches -- so the cost is asserted in calls, never wall-clock, which would
    make the emitted sequence machine-dependent and break design_hash."""

    def test_per_iteration_calls_are_bounded_by_max_candidates(self, code: object) -> None:
        # Whole-CDS window over six Gly codons: 4**6 = 4096 variants. On the merge
        # base 4096 <= EXHAUSTIVE_LIMIT, so the exhaustive branch enumerated all
        # 4096 EVERY iteration; the budget routes anything over max_candidates to
        # guided random with exactly max_candidates trials.
        protein = "M" + "G" * 6
        start = "ATG" + "GGA" * 6
        calls = {"n": 0}

        def finder(c: Construct) -> tuple[Breach, ...]:
            calls["n"] += 1
            # Permanent, fixable, whole-CDS: the search runs its full budget every
            # iteration and never clears it.
            return (
                Breach("r", Interval(0, len(c.sequence)), 1.0, "x", fixable_by_codon_choice=True),
            )

        max_candidates, max_iterations = 256, 5
        with pytest.raises(InfeasibleConstraints):
            repair(
                start,
                protein,
                code,  # type: ignore[arg-type]
                assemble=all_cds_assembler(protein),
                find_breaches=finder,
                window=50,
                seed=0,
                max_candidates=max_candidates,
                max_iterations=max_iterations,
            )
        # 1 initial probe + at most max_candidates per iteration. The merge base
        # would have made ~4096 * 5 here.
        assert calls["n"] <= max_candidates * max_iterations + 1
