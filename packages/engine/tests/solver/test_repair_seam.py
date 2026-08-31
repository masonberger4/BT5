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
from bt5.core.spec import Breach
from bt5.core.types import (
    Construct,
    Interval,
    Segment,
    SegmentKind,
    Topology,
    TranslationUnit,
)
from bt5.solver.repair import repair


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
