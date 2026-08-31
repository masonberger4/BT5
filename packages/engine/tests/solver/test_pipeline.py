"""The optimize() pipeline: Tier-A GC steering, the honest 'not run' state, and
the I9 backbone-reference arming.

These test the WIRING the pipeline is responsible for -- that Tier A steers when
a band is given, that a skipped Tier B is legible rather than fabricated, and
that the backbone-untouched invariant is actually handed the reference. The
repair search itself is covered in test_repair.py and test_repair_seam.py.
"""

from __future__ import annotations

import pytest
from bt5.codon.tables import FileTableProvider
from bt5.core.result import VerificationError
from bt5.core.spec import Breach
from bt5.core.types import (
    Construct,
    Interval,
    Segment,
    SegmentKind,
    Topology,
    TranslationUnit,
)
from bt5.solver.lattice import cai_lattice_scorer, optimal_back_translate
from bt5.solver.pipeline import optimize
from bt5.verify import gc_fraction

GC_LO, GC_HI, WIN = 0.40, 0.60, 30
GC_RICH = "M" + "AGPRAAGGPPRRAGPRAAGGPPRR"
NORMAL = "M" + "KLIWQRSTVNDEYFPGHACMKLIW"


@pytest.fixture(scope="module")
def env() -> tuple:
    p = FileTableProvider()
    return p.genetic_code(11), p.usage("sharp_li_1987_ecoli_w")


def cds_only_assembler(protein: str):
    """One editable CDS from base 0; no backbone, no flanks."""

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


class TestGcSteeringIsWired:
    def test_steering_reaches_the_band_with_tier_b_off(self, env: tuple) -> None:
        """With `find_breaches=None` Tier B never runs, so an in-band result can
        only be Tier A steering -- which `optimize()` did not do before."""
        code, w = env
        res = optimize(
            GC_RICH,
            code,
            assemble=cds_only_assembler(GC_RICH),
            find_breaches=None,
            score=cai_lattice_scorer(w.w),
            gc_bounds=(GC_LO, GC_HI),
            gc_window=WIN,
            seed=7,
        )
        assert res.tier_b_ran is False
        assert GC_LO <= gc_fraction(res.cds) <= GC_HI
        assert code.translate(res.cds)[:-1] == GC_RICH


class TestTierBLegibility:
    def test_no_find_breaches_is_marked_not_run(self, env: tuple) -> None:
        code, _ = env
        res = optimize(NORMAL, code, assemble=cds_only_assembler(NORMAL), find_breaches=None)
        assert res.tier_b_ran is False
        # Distinguishable from a clean repair that ran and found nothing.
        assert res.repair_outcome.stop_reason == "not_run"

    def test_supplying_find_breaches_marks_tier_b_ran(self, env: tuple) -> None:
        code, _ = env

        def clean_finder(_c: Construct) -> tuple[Breach, ...]:
            return ()

        res = optimize(
            NORMAL, code, assemble=cds_only_assembler(NORMAL), find_breaches=clean_finder
        )
        assert res.tier_b_ran is True
        assert res.repair_outcome.stop_reason == "clean"


class TestI9IsArmed:
    """optimize() must hand the parsed backbone to verify_construct, or the
    highest-value invariant -- the backbone was not silently edited -- is never
    checked on the design path."""

    def _assembler_with_backbone(self, protein: str, backbone_seq: str):
        def assemble(cds: str) -> Construct:
            n = len(cds)
            return Construct(
                sequence=cds + backbone_seq,
                topology=Topology.LINEAR,
                segments=(
                    Segment(Interval(0, n), SegmentKind.DESIGNABLE_CDS, "cds"),
                    Segment(Interval(n, n + len(backbone_seq)), SegmentKind.BACKBONE, "bb"),
                ),
                translation_units=(
                    TranslationUnit(
                        11, tuple(Interval(i, i + 3) for i in range(0, n, 3)), protein, True
                    ),
                ),
            )

        return assemble

    def _reference(self, protein: str, code: object, backbone_seq: str) -> Construct:
        n = len(optimal_back_translate(protein, code))  # type: ignore[arg-type]
        return Construct(
            sequence="A" * n + backbone_seq,
            topology=Topology.LINEAR,
            segments=(
                Segment(Interval(0, n), SegmentKind.DESIGNABLE_CDS, "cds"),
                Segment(Interval(n, n + len(backbone_seq)), SegmentKind.BACKBONE, "bb"),
            ),
        )

    def test_a_mismatched_backbone_reference_is_caught(self, env: tuple) -> None:
        code, _ = env
        good = "ACGTACGTAC"
        assembler = self._assembler_with_backbone(NORMAL, good)
        reference = self._reference(NORMAL, code, "TTTTTTTTTT")  # differs from `good`
        with pytest.raises(VerificationError) as exc:
            optimize(
                NORMAL,
                code,
                assemble=assembler,
                find_breaches=None,
                original_backbone=reference,
            )
        assert exc.value.invariant == "I9"

    def test_a_matching_backbone_reference_passes(self, env: tuple) -> None:
        code, _ = env
        good = "ACGTACGTAC"
        assembler = self._assembler_with_backbone(NORMAL, good)
        reference = self._reference(NORMAL, code, good)
        res = optimize(
            NORMAL,
            code,
            assemble=assembler,
            find_breaches=None,
            original_backbone=reference,
        )
        assert res.construct.sequence.endswith(good)
