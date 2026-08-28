"""Design results, conflicts, and failure.

`native_baseline` being a FIELD of DesignResult -- not a UI afterthought -- is how
"the honest answer is often to use the native sequence" becomes structurally
impossible to forget. For homologous mammalian expression it frequently IS the
right answer, and no vendor tool will ever say so.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Literal

from bt5.core.spec import Breach
from bt5.core.types import Construct, Interval, Provenance


@dataclass(frozen=True, slots=True)
class ObjectiveScore:
    """A rule's score, normalised so sliders act on it linearly.

    `percentile` is the ONLY quantity the weighted sum operates on. Raw scores are
    incommensurable -- kcal/mol against CAI in [0,1] against integer motif counts --
    and summing them unnormalised leaves the sliders dead over most of their range
    while one term silently dominates.

    The null MUST be built on the assembled construct in the same context as the
    score it normalises, or "94th percentile" is measured against a distribution
    that never contained a backbone and the report line is simply false.
    """

    spec_id: str
    raw: float
    unit: str
    percentile: float  # [0, 1]; 1.0 == better than every null variant
    null_n: int
    null_mean: float
    null_sd: float
    null_kind: Literal["host_frequency", "uniform_synonymous"] = "host_frequency"
    windowed_fold_only: bool = True  # nulls never use whole-transcript MFE
    #: Non-empty when this objective could NOT be evaluated -- no folding engine
    #: installed, a required input absent. `raw` and `percentile` are then
    #: meaningless and must not enter the weighted sum.
    #:
    #: A field rather than an omission, because the two are not the same to a
    #: reader: a scorecard missing its highest-weight objective looks exactly
    #: like a scorecard where that objective was never configured, and the
    #: difference is whether the ranking means anything. Silently dropping the
    #: term is the failure the whole degradation vocabulary exists to prevent,
    #: and nothing else in the contract could express it.
    unavailable_reason: str = ""

    @property
    def available(self) -> bool:
        return not self.unavailable_reason

    @classmethod
    def unavailable(cls, spec_id: str, unit: str, reason: str) -> ObjectiveScore:
        """An objective that could not be evaluated, stated rather than dropped."""
        if not reason:
            raise ValueError("an unavailable objective must say why")
        return cls(
            spec_id=spec_id,
            raw=float("nan"),
            unit=unit,
            percentile=float("nan"),
            null_n=0,
            null_mean=float("nan"),
            null_sd=float("nan"),
            unavailable_reason=reason,
        )


@dataclass(frozen=True, slots=True)
class Relaxation:
    """A specific, costed way out of a conflict."""

    spec_id: str
    change: str  # "raise gc_max 0.60 -> 0.64"
    predicted_cost: Mapping[str, float] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class Conflict:
    """Two rules that cannot both be satisfied here.

    Surfaced, never silently resolved. Positional conflicts are DISCOVERED from
    overlapping breaches; structural ones (NcoI CCATGG inside Kozak GCCACCATGG,
    CpG depletion against a vendor GC floor) are DECLARED via Spec.conflicts_with,
    because overlap alone would miss them.
    """

    interval: Interval
    spec_ids: tuple[str, ...]
    kind: Literal["mutually_exclusive", "opposing_gradient", "immutable_region", "declared"]
    binding_spec_id: str
    relaxations: tuple[Relaxation, ...] = ()


@dataclass(frozen=True, slots=True)
class InfeasibilityCertificate:
    """A PROOF of infeasibility, not a guess.

    A zero-variant merged segment in the mutation space, or a dead automaton
    state, is a certificate. "No solution" alone is a useless product -- over
    constrained failure is the dominant real-world experience with constraint
    solvers, so the minimal conflicting set and ranked relaxations are the
    deliverable, not an extra.
    """

    interval: Interval
    protein_span: tuple[int, int]
    minimal_conflicting_specs: tuple[str, ...]
    proof: Literal["empty_mutation_space", "automaton_dead_state", "immutable_region"]
    relaxations: tuple[Relaxation, ...] = ()


class InfeasibleConstraints(Exception):
    """Raised instead of returning a sequence that violates a hard constraint."""

    def __init__(self, certificate: InfeasibilityCertificate) -> None:
        self.certificate = certificate
        super().__init__(
            f"infeasible over {certificate.interval}: "
            f"{', '.join(certificate.minimal_conflicting_specs)}"
        )


class VerificationError(Exception):
    """The independent validator refused to emit. Raised by verify_construct."""

    def __init__(self, invariant: str, detail: str) -> None:
        self.invariant = invariant
        self.detail = detail
        super().__init__(f"{invariant}: {detail}")


@dataclass(frozen=True, slots=True)
class ScoreCard:
    scores: tuple[ObjectiveScore, ...]
    hard_checks: tuple[Breach, ...] = ()
    total: float = 0.0

    @property
    def available(self) -> tuple[ObjectiveScore, ...]:
        """The scores a weighted sum may legitimately operate on."""
        return tuple(s for s in self.scores if s.available)

    @property
    def unavailable(self) -> tuple[ObjectiveScore, ...]:
        """Objectives that could not be evaluated. A report that does not show
        these is claiming a completeness it does not have."""
        return tuple(s for s in self.scores if not s.available)


@dataclass(frozen=True, slots=True)
class Candidate:
    """One design. `design_hash` goes on the tube label.

    Two runs producing two different sequences under one name is how a lab ends up
    with two tubes and an irreproducible result, so the content hash travels onto
    the report, the GenBank note and the order file.
    """

    label: str
    construct: Construct
    cds: str
    scorecard: ScoreCard
    design_hash: str
    codon_distance_to: Mapping[str, float] = field(default_factory=dict)


@dataclass(frozen=True)
class DesignResult:
    candidates: tuple[Candidate, ...]
    native_baseline: Candidate | None = None
    conflicts: tuple[Conflict, ...] = ()
    provenance: Provenance | None = None
