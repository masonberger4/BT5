"""The conflict panel, which is a consequence of the data model rather than a feature.

Every Breach already carries an interval and a slot_role, so two rules demanding
opposite things over the same bases produce two Breaches over one Interval.
These tests are mostly about the cases where that is NOT enough: the origin, the
declared conflicts overlap can never find, and telling a user which of the two
rules they can actually do something about.
"""

from __future__ import annotations

from typing import ClassVar

from bt5.core.spec import (
    Breach,
    Citation,
    Direction,
    Enforcement,
    Evaluation,
    Evidence,
    LocalizationPolicy,
    RepairPolicy,
)
from bt5.core.types import Interval
from bt5.score import detect_conflicts, hard_versus_soft

N = 5000


def spec(
    spec_id: str,
    *,
    enforcement: Enforcement = Enforcement.SOFT,
    conflicts_with: tuple[str, ...] = (),
) -> type:
    class _Fake:
        id: ClassVar[str] = ""
        version: ClassVar[str] = "1.0.0"
        title: ClassVar[str] = ""
        enforcement: ClassVar[Enforcement] = Enforcement.SOFT
        evidence: ClassVar[Evidence] = Evidence.EVIDENCE_BACKED
        direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
        unit: ClassVar[str] = "au"
        citations: ClassVar[tuple[Citation, ...]] = (Citation("x", "https://example.org"),)
        last_verified: ClassVar[str] = "2026-08-28"
        weight_provenance: ClassVar[str] = "test"
        default_enabled: ClassVar[bool] = True
        default_weight: ClassVar[float] = 1.0
        steering_weight: ClassVar[float] = 0.0
        band: ClassVar[tuple[float, float] | None] = None
        localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.WHOLE_SCOPE
        repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
        cost_class: ClassVar[str] = "cheap"
        conflicts_with: ClassVar[tuple[str, ...]] = ()
        param_schema: ClassVar[dict[str, object]] = {}
        brief_ref: ClassVar[str] = "2.X1"
        engine_calibration: ClassVar[str | None] = None

    _Fake.id = spec_id
    _Fake.enforcement = enforcement
    _Fake.conflicts_with = conflicts_with
    _Fake.__name__ = spec_id
    return _Fake


def breach(
    spec_id: str,
    iv: Interval,
    *,
    magnitude: float = 1.0,
    fixable: bool = True,
    side: str = "",
) -> Breach:
    return Breach(
        spec_id=spec_id,
        interval=iv,
        magnitude=magnitude,
        message="",
        fixable_by_codon_choice=fixable,
        detail={"binding_side": side} if side else {},
    )


def ev(*breaches: Breach) -> Evaluation:
    spec_id = breaches[0].spec_id if breaches else "x"
    return Evaluation(spec_id=spec_id, passes=not breaches, raw_score=0.0, breaches=breaches)


SPECS = {
    "a_soft": spec("a_soft"),
    "b_soft": spec("b_soft"),
    "c_hard": spec("c_hard", enforcement=Enforcement.HARD_LATTICE),
}


class TestPositional:
    def test_two_rules_over_the_same_bases_is_a_conflict(self) -> None:
        conflicts = detect_conflicts(
            [ev(breach("a_soft", Interval(100, 130))), ev(breach("b_soft", Interval(120, 160)))],
            SPECS,
            length=N,
            circular=True,
        )
        assert len(conflicts) == 1
        assert conflicts[0].spec_ids == ("a_soft", "b_soft")
        assert conflicts[0].interval == Interval(100, 160)

    def test_one_rule_unhappy_alone_is_a_breach_not_a_conflict(self) -> None:
        conflicts = detect_conflicts(
            [ev(breach("a_soft", Interval(100, 130)), breach("a_soft", Interval(120, 160)))],
            SPECS,
            length=N,
            circular=True,
        )
        assert conflicts == ()

    def test_rules_that_never_meet_do_not_conflict(self) -> None:
        conflicts = detect_conflicts(
            [ev(breach("a_soft", Interval(100, 130))), ev(breach("b_soft", Interval(900, 960)))],
            SPECS,
            length=N,
            circular=True,
        )
        assert conflicts == ()

    def test_a_conflict_across_the_origin_is_found(self) -> None:
        """The linear predicate reports a plasmid as conflict-free precisely
        where it has one: a breach stored as [4980, 5020) and one at [10, 40)
        sit on the same bases."""
        conflicts = detect_conflicts(
            [ev(breach("a_soft", Interval(4980, 5020))), ev(breach("b_soft", Interval(10, 40)))],
            SPECS,
            length=N,
            circular=True,
        )
        assert len(conflicts) == 1
        assert conflicts[0].spec_ids == ("a_soft", "b_soft")

    def test_the_same_pair_on_a_linear_construct_does_not_conflict(self) -> None:
        conflicts = detect_conflicts(
            [ev(breach("a_soft", Interval(4980, 5020))), ev(breach("b_soft", Interval(10, 40)))],
            SPECS,
            length=N,
            circular=False,
        )
        assert conflicts == ()

    def test_the_reported_interval_stays_small_across_the_origin(self) -> None:
        """A conflict interval spanning most of the plasmid points at nothing.
        The merged span must be the short way round, in the canonical form."""
        conflicts = detect_conflicts(
            [ev(breach("a_soft", Interval(4980, 5020))), ev(breach("b_soft", Interval(10, 40)))],
            SPECS,
            length=N,
            circular=True,
        )
        iv = conflicts[0].interval
        assert iv.length <= 100, f"merged span {iv} is the long way round"
        assert iv.start < N <= iv.end, "and it is stored in the one wrapping representation"

    def test_a_chain_of_overlaps_is_one_conflict_naming_every_rule(self) -> None:
        """Splitting one contiguous problem into two conflicts would show the
        user two panels each naming half the rules involved."""
        conflicts = detect_conflicts(
            [
                ev(breach("a_soft", Interval(100, 140))),
                ev(breach("c_hard", Interval(130, 170))),
                ev(breach("b_soft", Interval(160, 200))),
            ],
            SPECS,
            length=N,
            circular=True,
        )
        assert len(conflicts) == 1
        assert conflicts[0].spec_ids == ("a_soft", "b_soft", "c_hard")


class TestKind:
    def test_opposing_bounds_are_an_opposing_gradient(self) -> None:
        conflicts = detect_conflicts(
            [
                ev(breach("a_soft", Interval(100, 150), side="upper")),
                ev(breach("b_soft", Interval(120, 160), side="lower")),
            ],
            SPECS,
            length=N,
            circular=True,
        )
        assert conflicts[0].kind == "opposing_gradient"

    def test_an_unfixable_member_makes_it_an_immutable_region(self) -> None:
        """It changes what the user can do, from 'relax one of these' to 'no
        codon anywhere resolves this', so it outranks the other labels."""
        conflicts = detect_conflicts(
            [
                ev(breach("a_soft", Interval(100, 150), side="upper")),
                ev(breach("b_soft", Interval(120, 160), side="lower", fixable=False)),
            ],
            SPECS,
            length=N,
            circular=True,
        )
        assert conflicts[0].kind == "immutable_region"

    def test_otherwise_mutually_exclusive(self) -> None:
        conflicts = detect_conflicts(
            [ev(breach("a_soft", Interval(100, 150))), ev(breach("b_soft", Interval(120, 160)))],
            SPECS,
            length=N,
            circular=True,
        )
        assert conflicts[0].kind == "mutually_exclusive"


class TestBinding:
    def test_a_hard_rule_binds_over_a_soft_one_regardless_of_magnitude(self) -> None:
        """Naming the soft rule would point the user at a slider that cannot
        help: the hard rule is not a preference."""
        conflicts = detect_conflicts(
            [
                ev(breach("a_soft", Interval(100, 150), magnitude=99.0)),
                ev(breach("c_hard", Interval(120, 160), magnitude=0.01)),
            ],
            SPECS,
            length=N,
            circular=True,
        )
        assert conflicts[0].binding_spec_id == "c_hard"

    def test_among_equals_the_largest_magnitude_binds(self) -> None:
        conflicts = detect_conflicts(
            [
                ev(breach("a_soft", Interval(100, 150), magnitude=0.2)),
                ev(breach("b_soft", Interval(120, 160), magnitude=0.9)),
            ],
            SPECS,
            length=N,
            circular=True,
        )
        assert conflicts[0].binding_spec_id == "b_soft"

    def test_hard_versus_soft_is_flagged_for_the_ui(self) -> None:
        conflicts = detect_conflicts(
            [
                ev(breach("a_soft", Interval(100, 150))),
                ev(breach("c_hard", Interval(120, 160))),
            ],
            SPECS,
            length=N,
            circular=True,
        )
        assert hard_versus_soft(conflicts[0], SPECS), (
            "not a trade-off the user can slide out of; showing it beside real "
            "trade-offs invites them to try"
        )

    def test_two_soft_rules_are_not_flagged_as_hard_versus_soft(self) -> None:
        conflicts = detect_conflicts(
            [ev(breach("a_soft", Interval(100, 150))), ev(breach("b_soft", Interval(120, 160)))],
            SPECS,
            length=N,
            circular=True,
        )
        assert not hard_versus_soft(conflicts[0], SPECS)


class TestDeclared:
    """Conflicts overlap can never find.

    NcoI CCATGG inside Kozak GCCACCATGG is the positional case. CpG depletion
    against a vendor GC floor is not: it is a disagreement about the whole
    sequence's composition, and waiting for two breaches to land in one window
    means never reporting it.
    """

    specs = {
        "cpg": spec("cpg", conflicts_with=("gc_floor",)),
        "gc_floor": spec("gc_floor", enforcement=Enforcement.HARD_REPAIR),
    }

    def test_declared_pair_is_reported_when_one_side_fires(self) -> None:
        conflicts = detect_conflicts(
            [ev(breach("cpg", Interval(100, 150)))], self.specs, length=N, circular=True
        )
        declared = [c for c in conflicts if c.kind == "declared"]
        assert len(declared) == 1
        assert declared[0].spec_ids == ("cpg", "gc_floor")

    def test_it_is_reported_once_not_once_per_direction(self) -> None:
        mutual = {
            "cpg": spec("cpg", conflicts_with=("gc_floor",)),
            "gc_floor": spec("gc_floor", conflicts_with=("cpg",)),
        }
        conflicts = detect_conflicts(
            [ev(breach("cpg", Interval(100, 150)))], mutual, length=N, circular=True
        )
        assert len([c for c in conflicts if c.kind == "declared"]) == 1

    def test_nothing_firing_reports_nothing(self) -> None:
        """Otherwise every construct carries the same declared conflicts and the
        panel is noise."""
        assert detect_conflicts([], self.specs, length=N, circular=True) == ()

    def test_a_conflict_with_an_absent_rule_is_not_reported(self) -> None:
        only_one = {"cpg": spec("cpg", conflicts_with=("gc_floor",))}
        conflicts = detect_conflicts(
            [ev(breach("cpg", Interval(100, 150)))], only_one, length=N, circular=True
        )
        assert conflicts == ()


def test_detection_does_not_read_the_registry() -> None:
    """The conflict panel is a reported artefact. It must not change because an
    unrelated import registered another rule, so `specs` is a parameter."""
    import inspect

    from bt5.score import conflicts as module

    source = inspect.getsource(module)
    assert "all_specs" not in source
    assert "from bt5.core.registry" not in source
