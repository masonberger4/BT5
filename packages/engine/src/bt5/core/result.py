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
