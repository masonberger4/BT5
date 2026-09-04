"""Tier B repair, GC steering, and the composed HARD_REPAIR guarantee."""

from __future__ import annotations

import pytest
from bt5.codon.tables import FileTableProvider
from bt5.core.result import InfeasibleConstraints, VerificationError
from bt5.core.spec import Breach, LocalizationPolicy
from bt5.core.types import (
    Construct,
    Interval,
    Segment,
    SegmentKind,
    Topology,
    TranslationUnit,
)
from bt5.solver.lattice import (
    achievable_gc_range,
    cai_lattice_scorer,
    gc_steering_scorer,
    optimal_back_translate,
    solve_with_gc_steering,
)
from bt5.solver.pipeline import optimize
from bt5.solver.repair import PER_TARGET_TOLERANCE, codon_span, localize, repair
from bt5.verify import gc_fraction

FORBID = ["GAATTC", "GGATCC", "GGTCTC"]
GC_LO, GC_HI, WIN = 0.40, 0.60, 30

GC_RICH = "M" + "AGPRAAGGPPRRAGPRAAGGPPRR"
AT_RICH = "M" + "KNIFKKNNIIFFKNIFKKNNIIFF"
NORMAL = "M" + "KLIWQRSTVNDEYFPGHACMKLIW"


@pytest.fixture(scope="module")
def env() -> tuple:
    p = FileTableProvider()
    return p.genetic_code(11), p.usage("sharp_li_1987_ecoli_w")


def make_assembler(protein: str):
    def assemble(cds: str) -> Construct:
        return Construct(
            sequence=cds,
            topology=Topology.LINEAR,
            segments=(Segment(Interval(0, len(cds)), SegmentKind.DESIGNABLE_CDS, "cds"),),
            translation_units=(
                TranslationUnit(
                    11,
                    tuple(Interval(i, i + 3) for i in range(0, len(cds), 3)),
                    protein,
                    True,
                ),
            ),
        )

    return assemble


def gc_breaches(c: Construct) -> tuple[Breach, ...]:
    out = []
    for start in range(0, max(1, c.length - WIN + 1), max(1, WIN // 5)):
        w = c.sequence[start : start + WIN]
        if len(w) < WIN:
            continue
        f = gc_fraction(w)
        if f < GC_LO:
            dev, side = GC_LO - f, "lower"
        elif f > GC_HI:
            dev, side = f - GC_HI, "upper"
        else:
            continue
        out.append(
            Breach(
                "e2_gc_band",
                Interval(start, start + WIN),
                dev,
                f"GC {f:.1%} at {start} ({side})",
                fixable_by_codon_choice=c.overlaps_editable(Interval(start, start + WIN)),
                detail={"binding_side": side},
            )
        )
    return tuple(out)


class TestGcSteering:
    """The steering term is a SIGNED Lagrangian, not a distance penalty."""

    def test_positive_multiplier_lowers_gc_and_negative_raises_it(self, env: tuple) -> None:
        code, _ = env
        down = optimal_back_translate(GC_RICH, code, score=gc_steering_scorer(1.0))
        up = optimal_back_translate(GC_RICH, code, score=gc_steering_scorer(-1.0))
        assert gc_fraction(down) < gc_fraction(up)

    def test_a_symmetric_penalty_would_be_wrong(self, env: tuple) -> None:
        """Regression for a real bug: penalising abs(gc - 0.5) per codon is
        symmetric, so on a protein whose codons all sit at 2/3 or 3/3 GC it pulls
        toward 2/3 and can RAISE global GC while appearing to steer it."""
        code, _ = env

        def symmetric(_i: int, codon: str, _p: str) -> float:
            gc = (codon.count("G") + codon.count("C")) / 3.0
            return 6.0 * abs(gc - 0.50)

        bad = optimal_back_translate(GC_RICH, code, score=symmetric)
        good = optimal_back_translate(GC_RICH, code, score=gc_steering_scorer(6.0))
        assert gc_fraction(good) < gc_fraction(bad)

    def test_achievable_range_is_bounded_by_amino_acid_composition(self, env: tuple) -> None:
        """Lys/Asn/Ile/Phe/Met codons carry at most one GC base each, so ~33% is
        a hard ceiling no matter how expensive GC is made."""
        code, _ = env
        lo, hi = achievable_gc_range(AT_RICH, code)
        assert hi == pytest.approx(1 / 3, abs=0.02)
        assert lo < hi

    def test_steering_reaches_the_band_when_it_is_reachable(self, env: tuple) -> None:
        code, w = env
        seq, _ = solve_with_gc_steering(
            GC_RICH,
            code,
            gc_bounds=(GC_LO, GC_HI),
            base_score=cai_lattice_scorer(w.w),
            forbidden=FORBID,
        )
        assert GC_LO <= gc_fraction(seq) <= GC_HI
        assert code.translate(seq)[:-1] == GC_RICH

    def test_unreachable_band_is_detected_rather_than_searched_for(self, env: tuple) -> None:
        code, w = env
        lo, hi = achievable_gc_range(AT_RICH, code, forbidden=FORBID)
        assert hi < GC_LO, "this fixture must have an unreachable band"
        seq, _ = solve_with_gc_steering(
            AT_RICH,
            code,
            gc_bounds=(GC_LO, GC_HI),
            base_score=cai_lattice_scorer(w.w),
            forbidden=FORBID,
        )
        # Returns the closest attainable encoding; the validator refuses later.
        assert gc_fraction(seq) == pytest.approx(hi, abs=0.02)

    def test_no_steering_applied_when_already_in_band(self, env: tuple) -> None:
        code, w = env
        seq, mult = solve_with_gc_steering(
            NORMAL,
            code,
            gc_bounds=(GC_LO, GC_HI),
            base_score=cai_lattice_scorer(w.w),
            forbidden=FORBID,
        )
        assert mult == 0.0
        assert GC_LO <= gc_fraction(seq) <= GC_HI


class TestLocalization:
    def test_window_policy_extends_by_window_minus_one(self) -> None:
        b = Breach("x", Interval(60, 90), 0.1, "", fixable_by_codon_choice=True)
        iv = localize(
            b,
            LocalizationPolicy.WINDOW_MINUS_1,
            window=30,
            motif_len=6,
            construct_length=300,
            circular=False,
        )
        assert iv.start == 60 - 29
        assert iv.end == 90 + 29

    def test_motif_policy_extends_by_motif_len_minus_one(self) -> None:
        b = Breach("x", Interval(60, 66), 1.0, "", fixable_by_codon_choice=True)
        iv = localize(
            b,
            LocalizationPolicy.MOTIF_LEN_MINUS_1,
            window=30,
            motif_len=6,
            construct_length=300,
            circular=False,
        )
        assert iv.start == 55
        assert iv.end == 71

    def test_whole_scope_covers_the_construct(self) -> None:
        b = Breach("x", Interval(10, 20), 1.0, "", fixable_by_codon_choice=True)
        iv = localize(
            b,
            LocalizationPolicy.WHOLE_SCOPE,
            window=30,
            motif_len=6,
            construct_length=300,
            circular=False,
        )
        assert iv == Interval(0, 300)

    def test_codon_span_returns_whole_codons_only(self) -> None:
        codon_map = tuple(Interval(i, i + 3) for i in range(0, 30, 3))
        first, last = codon_span(Interval(4, 11), codon_map, 30)
        assert (first, last) == (1, 4)  # bases 4-10 touch codons 1,2,3


class TestRepair:
    def test_repair_clears_gc_breaches(self, env: tuple) -> None:
        code, w = env
        assemble = make_assembler(GC_RICH)
        start = optimal_back_translate(
            GC_RICH, code, forbidden=FORBID, score=cai_lattice_scorer(w.w)
        )
        assert gc_breaches(assemble(start)), "fixture must start out of band"
        out = repair(
            start,
            GC_RICH,
            code,
            assemble=assemble,
            find_breaches=gc_breaches,
            forbidden=FORBID,
            window=WIN,
            seed=3,
        )
        assert out.clean
        assert not gc_breaches(assemble(out.cds))

    def test_repair_never_introduces_a_forbidden_motif(self, env: tuple) -> None:
        """THE critical interaction. Repair mutates codons to fix a windowed
        statistic, and a mutation can create a motif Tier A guaranteed away --
        including across a codon boundary. Tier B must never weaken Tier A."""
        code, w = env
        assemble = make_assembler(GC_RICH)
        start = optimal_back_translate(
            GC_RICH, code, forbidden=FORBID, score=cai_lattice_scorer(w.w)
        )
        out = repair(
            start,
            GC_RICH,
            code,
            assemble=assemble,
            find_breaches=gc_breaches,
            forbidden=FORBID,
            window=WIN,
            seed=3,
        )
        for motif in FORBID:
            assert motif not in out.cds

    def test_repair_preserves_the_protein(self, env: tuple) -> None:
        code, w = env
        assemble = make_assembler(GC_RICH)
        start = optimal_back_translate(
            GC_RICH, code, forbidden=FORBID, score=cai_lattice_scorer(w.w)
        )
        out = repair(
            start,
            GC_RICH,
            code,
            assemble=assemble,
            find_breaches=gc_breaches,
            forbidden=FORBID,
            window=WIN,
            seed=3,
        )
        assert code.translate(out.cds)[:-1] == GC_RICH

    def test_repair_is_deterministic_under_a_fixed_seed(self, env: tuple) -> None:
        code, w = env
        assemble = make_assembler(GC_RICH)
        start = optimal_back_translate(
            GC_RICH, code, forbidden=FORBID, score=cai_lattice_scorer(w.w)
        )

        def run() -> str:
            return repair(
                start,
                GC_RICH,
                code,
                assemble=assemble,
                find_breaches=gc_breaches,
                forbidden=FORBID,
                window=WIN,
                seed=11,
            ).cds

        assert run() == run()

    def test_a_clean_input_is_returned_untouched(self, env: tuple) -> None:
        code, w = env
        assemble = make_assembler(NORMAL)
        start = optimal_back_translate(
            NORMAL, code, forbidden=FORBID, score=cai_lattice_scorer(w.w)
        )
        out = repair(
            start,
            NORMAL,
            code,
            assemble=assemble,
            find_breaches=gc_breaches,
            forbidden=FORBID,
            window=WIN,
            seed=1,
        )
        assert out.cds == start
        assert out.iterations == 0

    def test_unfixable_breach_raises_with_a_certificate(self, env: tuple) -> None:
        """When the protein's own composition puts the band out of reach, repair
        cannot help and must say which rule conflicts rather than fail bare."""
        code, _ = env
        assemble = make_assembler(AT_RICH)
        start = optimal_back_translate(AT_RICH, code, forbidden=FORBID)
        with pytest.raises(InfeasibleConstraints) as exc:
            repair(
                start,
                AT_RICH,
                code,
                assemble=assemble,
                find_breaches=gc_breaches,
                forbidden=FORBID,
                window=WIN,
                seed=5,
                max_iterations=40,
            )
        assert "e2_gc_band" in exc.value.certificate.minimal_conflicting_specs


class TestPipelineComposition:
    """HARD_REPAIR is only real as the composition of all three tiers."""

    def test_optimize_produces_a_verified_in_band_construct(self, env: tuple) -> None:
        code, w = env
        res = optimize(
            GC_RICH,
            code,
            assemble=make_assembler(GC_RICH),
            find_breaches=gc_breaches,
            forbidden=FORBID,
            score=cai_lattice_scorer(w.w),
            gc_bounds=(GC_LO, GC_HI),
            gc_window=WIN,
            seed=7,
        )
        assert GC_LO <= gc_fraction(res.cds) <= GC_HI
        assert code.translate(res.cds)[:-1] == GC_RICH
        assert not any(m in res.cds for m in FORBID)

    def test_the_validator_refuses_rather_than_emitting_out_of_band(self, env: tuple) -> None:
        """The refusal IS the guarantee. Without it, 'hard' would mean 'we tried'."""
        code, _ = env
        with pytest.raises((VerificationError, InfeasibleConstraints)):
            optimize(
                AT_RICH,
                code,
                assemble=make_assembler(AT_RICH),
                find_breaches=gc_breaches,
                forbidden=FORBID,
                gc_bounds=(GC_LO, GC_HI),
                gc_window=WIN,
                seed=7,
            )

    def test_verification_is_on_by_default(self) -> None:
        import inspect

        from bt5.solver import pipeline

        assert inspect.signature(pipeline.optimize).parameters["_verify"].default is True


class TestAHopelessBreachDoesNotBurnTheBudget:
    """#111. A breach no codon choice can clear must stop being re-attacked
    while the sequence stands still -- without that give-up being dressed up as
    a proof.

    A failed iteration leaves `current` untouched, so re-targeting the same
    breach re-runs the same window over the same options with only a fresh RNG
    draw. `STAGNATION_TOLERANCE` alone bounds the GLOBAL run of failures, so one
    unclearable FIXED_POINT breach spent 100 iterations x `max_candidates`
    evaluations before giving up: 25,601 `find_breaches` calls, measured. In the
    design lane, where `find_breaches` is the wired catalog rather than this
    fixture's one-liner, that was ~260 s per solve and the candidate was then
    discarded anyway.
    """

    PROTEIN = "M" + "KLIWQRSTVNDEYFPGHACMKLIW"

    def _repair(self, code: object, **kw: object):
        """Repair against a breach that is reported forever, whatever the codons."""
        calls = {"n": 0}

        def find_breaches(construct: Construct) -> list[Breach]:
            calls["n"] += 1
            return [Breach("unfixable", Interval(30, 40), 1.0, "", fixable_by_codon_choice=True)]

        start = "".join(code.families()[a][0] for a in self.PROTEIN)  # type: ignore[attr-defined]
        out = repair(
            start,
            self.PROTEIN,
            code,  # type: ignore[arg-type]
            assemble=make_assembler(self.PROTEIN),
            find_breaches=find_breaches,
            seed=3,
            raise_on_infeasible=False,
            **kw,  # type: ignore[arg-type]
        )
        return out, calls["n"]

    def test_it_gives_up_within_the_per_target_budget(self, env: tuple) -> None:
        """The bound is PER TARGET, so it does not scale with the global
        tolerance. One breach, so at most `PER_TARGET_TOLERANCE` searches of it
        plus the selection pass that finds nothing left."""
        code, _ = env
        out, calls = self._repair(code)
        assert out.random_windows == PER_TARGET_TOLERANCE
        assert out.iterations <= PER_TARGET_TOLERANCE + 1
        # 1 initial survey + PER_TARGET_TOLERANCE searches of max_candidates each.
        assert calls == 1 + PER_TARGET_TOLERANCE * 256
        assert calls < 25_601, "the pre-#111 cost; a regression here is the bug returning"

    def test_giving_up_is_never_reported_as_convergence(self, env: tuple) -> None:
        """The whole risk of abandoning a target early. `exhausted_targets` sets
        `converged`, which claims the space was searched out and the breach
        cannot be cleared. Abandonment claims only that the search kept missing,
        so it must report `stagnation` and leave `converged` False -- otherwise
        this optimisation would manufacture an infeasibility proof."""
        code, _ = env
        out, _ = self._repair(code)
        assert out.stop_reason == "stagnation"
        assert out.converged is False
        assert out.clean is False
        assert out.remaining, "the breach is still there and must be carried"

    def test_it_still_repairs_what_it_can(self, env: tuple) -> None:
        """The guard must not cost a repair that was reachable. GC_RICH starts
        out of band and is clearable; abandonment must never fire on a search
        that is making progress, because any improvement clears the ledger."""
        code, w = env
        assemble = make_assembler(GC_RICH)
        start = optimal_back_translate(
            GC_RICH, code, forbidden=FORBID, score=cai_lattice_scorer(w.w)
        )
        assert gc_breaches(assemble(start)), "fixture must start out of band"
        out = repair(
            start,
            GC_RICH,
            code,
            assemble=assemble,
            find_breaches=gc_breaches,
            forbidden=FORBID,
            window=WIN,
            seed=3,
        )
        assert out.clean
        assert not gc_breaches(assemble(out.cds))
