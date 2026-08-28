"""The rule protocol.

Every scientific rule in BT5 implements `Spec`. Three properties of this protocol
are load-bearing and each retires a specific failure mode:

- `Breach.interval` + `Breach.slot_role` make the conflict panel and the
  "which of the three contexts fired this" display FREE, because every rule is
  forced to produce them. Retrofitting attribution later would touch every rule.
- `lattice_terms()` decouples the solver from the rules completely. A rule can
  ship Tier-B-only on day 3 and upgrade itself to Tier-A on day 20 without the
  solver author touching a file.
- `Enforcement` makes it structurally impossible to "guarantee" a hard constraint
  with a penalty weight. `pipeline` asserts a HARD rule never reaches the
  weighted sum.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from enum import StrEnum
from typing import TYPE_CHECKING, ClassVar, Protocol, runtime_checkable

from bt5.core.types import Interval, Strand

if TYPE_CHECKING:
    from bt5.core.context import ContextSlot, DesignContext
    from bt5.core.services import Services
    from bt5.core.types import Construct


class Enforcement(StrEnum):
    """How a rule is enforced. NOT a severity -- a mechanism.

    HARD_LATTICE  guaranteed by construction in the Tier-A DP. Cannot be violated.
    HARD_REPAIR   not expressible in the lattice; enforced by Tier-B repair and
                  then PROVEN by the independent validator, which refuses to emit
                  on failure. Windowed GC content is the canonical member.
    HARD_CHECK    real, but not fixable by codon choice (an upstream AUG in the
                  user's own 5'UTR, a transmembrane-segment toxicity call, an AAV
                  size overflow). Reported and blocking, never chased by the solver.
    SOFT          a weighted objective. The ONLY class the weighted sum may see.
    REPORT_ONLY   surfaced to the user, never scored, never enforced.
    """

    HARD_LATTICE = "hard_lattice"
    HARD_REPAIR = "hard_repair"
    HARD_CHECK = "hard_check"
    SOFT = "soft"
    REPORT_ONLY = "report_only"

    @property
    def is_hard(self) -> bool:
        return self in (Enforcement.HARD_LATTICE, Enforcement.HARD_REPAIR, Enforcement.HARD_CHECK)

    @property
    def is_scored(self) -> bool:
        return self is Enforcement.SOFT


class Evidence(StrEnum):
    """Four levels, deliberately not three.

    Collapsing VENDOR_ASSERTED into FOLKLORE would ship every manufacturability
    rule disabled by default -- exactly backwards, since vendors enforce those at
    order time. Only FOLKLORE defaults off.
    """

    EVIDENCE_BACKED = "evidence_backed"
    CONTESTED = "contested"
    VENDOR_ASSERTED = "vendor_asserted"
    FOLKLORE = "folklore"


class Direction(StrEnum):
    """Which way is better, in the rule's own native units.

    BAND exists because a forced "lower is always better" scalar cannot represent
    CAI's 0.70-0.90 target, global GC 40-60%, or the two-sided 45-60% GC AT-window
    -- and a monotone weighted sum over a collapsed |deviation| drives CAI toward
    1.0, which is the precise failure the evidence refutes (max-CAI collapses to
    one codon per amino acid and produces perfect nucleotide repeats).
    """

    LOWER_IS_BETTER = "lower_is_better"
    HIGHER_IS_BETTER = "higher_is_better"
    BAND = "band"


class LocalizationPolicy(StrEnum):
    """How Tier-B widens a breach into a repair window.

    Declared as DATA and consumed generically, so ~45 rule authors do not each
    hand-roll DnaChisel's localisation heuristics and each get the
    two-overlapping-breaches suppression rule wrong independently.
    """

    MOTIF_LEN_MINUS_1 = "motif_len_minus_1"
    WINDOW_MINUS_1 = "window_minus_1"
    WHOLE_SCOPE = "whole_scope"
    PAIRED_SEGMENTS = "paired_segments"


class RepairPolicy(StrEnum):
    """SINGLE_PASS is unsafe for splice-site removal.

    Point-mutating one cryptic donor activates cryptic donors nearby (the A2UCOE
    case). A single-pass "remove motifs" step ships a construct whose donors were
    removed into NEW donors, and the validator passes it because the specific
    9-mer it was told to avoid is gone. FIXED_POINT re-scans until no new breach
    appears, with an iteration cap and a documented failure mode.
    """

    SINGLE_PASS = "single_pass"
    FIXED_POINT = "fixed_point"


@dataclass(frozen=True, slots=True)
class Citation:
    """A tuple of these, never a single URL.

    The rules that most need an honest evidence badge are the contested ones
    resting on two sources with OPPOSITE signs -- ZAP CpG, codon-pair bias,
    G-quadruplexes in mammalian cells. One string makes the badge dishonest
    precisely where honesty is the differentiator. `sign` records which way each
    source cuts.
    """

    label: str
    url: str
    year: int | None = None
    sign: str = "supports"  # "supports" | "refutes" | "qualifies"

    def __post_init__(self) -> None:
        if not self.url.startswith("https://"):
            raise ValueError(f"citation url must be https, got {self.url!r}")


@dataclass(frozen=True, slots=True)
class Breach:
    """One localized, attributable problem.

    The unit of everything downstream: the conflict panel, the infeasibility
    certificate, the per-window "which side is binding" display, and the report
    all consume ONLY Breaches.
    """

    spec_id: str
    interval: Interval
    magnitude: float  # rule-native; > 0 means worse
    message: str  # must name the exact offending substring
    #: Can the solver do anything about this by choosing different codons?
    #:
    #: REQUIRED, with no default, because both defaults are wrong. Defaulting
    #: True sends the solver after a uAUG in the user's own 5'UTR, an AAV size
    #: overflow or a toxic TM segment -- none of which any codon can move -- and
    #: it exhausts the mutation space and reports infeasible on a design that
    #: was fine. Defaulting False silently drops real, fixable motifs from the
    #: solver until the independent validator refuses to emit at the very end.
    #:
    #: A rule author knows this about their own finding and nobody downstream
    #: can recover it, so the contract makes them say it. False routes the
    #: breach to the advisor instead of the solver.
    fixable_by_codon_choice: bool
    slot_role: str | None = None  # "propagation" | "producer" | "target"
    detail: Mapping[str, float | str] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class Evaluation:
    """The result of running one rule against one construct."""

    spec_id: str
    passes: bool
    raw_score: float  # in the rule's NATIVE units, see `unit`/`direction` on the Spec
    breaches: tuple[Breach, ...] = ()
    windows: tuple[tuple[Interval, float], ...] = ()
    n_evaluated: int = 0  # denominator, for honest reporting
    binding_side: str | None = None  # "upper" | "lower" for BAND rules


@dataclass(frozen=True, slots=True)
class LatticeTerms:
    """A rule's OPT-IN to Tier-A exactness. Returning None means Tier-B only.

    `forbidden` patterns are closed under reverse complement by the SOLVER at
    automaton-construction time, so a rule lists only the forward motif and gets
    both strands guaranteed. This closure applies to `forbidden` ONLY -- directional
    scored models (MaxEntScan, the Salis promoter calculator, TIR/RBS, polyA
    hexamer + downstream element) are NOT reverse-complement symmetric and must
    read `packaged_strand` from the context instead.
    """

    forbidden: tuple[str, ...] = ()  # IUPAC; consumed by the Aho-Corasick automaton
    codon_weights: Mapping[str, float] | None = None
    codon_pair_weights: Mapping[tuple[str, str], float] | None = None
    positional: Callable[[int, str], float] | None = None  # h(codon_index, codon)


@runtime_checkable
class Spec(Protocol):
    """Every scientific rule implements this.

    Class-level metadata is what the UI renders from and what CI audits. A rule
    missing a citation or (if SOFT) a weight_provenance fails the contract test.
    """

    id: ClassVar[str]
    version: ClassVar[str]
    title: ClassVar[str]
    enforcement: ClassVar[Enforcement]
    evidence: ClassVar[Evidence]
    direction: ClassVar[Direction]
    unit: ClassVar[str]
    citations: ClassVar[tuple[Citation, ...]]
    last_verified: ClassVar[str]  # ISO date; vendor rules go stale within ~12 months
    weight_provenance: ClassVar[str]  # required non-empty for SOFT rules
    default_enabled: ClassVar[bool]  # FOLKLORE ships False
    #: Weight in the OBJECTIVE FUNCTION. Must be 0.0 for every hard rule -- the
    #: weighted sum only ever sees SOFT rules, and `pipeline` asserts it.
    default_weight: ClassVar[float]
    #: Strength of the Tier-A Lagrangian STEERING term. Distinct from
    #: default_weight and not part of the objective. A HARD_REPAIR rule such as
    #: windowed GC needs the DP nudged toward its band even though enforcement
    #: comes from repair plus the independent validator, never from a penalty.
    steering_weight: ClassVar[float]
    band: ClassVar[tuple[float, float] | None]
    localization: ClassVar[LocalizationPolicy]
    repair: ClassVar[RepairPolicy]
    cost_class: ClassVar[str]  # "cheap" | "moderate" | "expensive"; drives null sampling
    conflicts_with: ClassVar[tuple[str, ...]]  # DECLARED structural conflicts
    param_schema: ClassVar[Mapping[str, object]]  # JSON Schema; UI renders controls
    brief_ref: ClassVar[str]  # section id in docs/research/brief.md
    engine_calibration: ClassVar[str | None]  # e.g. "viennarna:rna_turner2004"

    def gate(self, slot: ContextSlot) -> bool:
        """Does this rule apply in this context slot?"""
        ...

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        """Enforcement can depend on delivery context.

        Internal polyA is HARD on the lentiviral sense strand but soft elsewhere;
        cryptic splice donors are HARD for lentiviral and genome-integrated but
        warn-only for plasmid; G-quadruplex is hard for bacterial CDS and soft for
        mammalian mRNA. One frozen level per rule is wrong in every job.
        """
        ...

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        """Evaluate against the ASSEMBLED construct. Must be pure."""
        ...

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """Opt in to Tier-A exactness, or return None for Tier-B repair only."""
        ...


def strand_for(ctx: DesignContext, slot: ContextSlot) -> Strand:
    """Resolve, in CONSTRUCT coordinates, the strand a directional model reads.

    For a reverse-oriented lentiviral cassette the packaged genome is the reverse
    complement of the cassette, so a rule that hard-codes the forward strand runs
    its polyA and splice-donor analysis exactly backwards and returns clean.
    Rules must call this rather than assuming.

    Two independent facts compose here, which is the whole reason this is a
    function and not an attribute read:

    - `ctx.cassette_orientation` -- which strand of the PLASMID the cassette
      reads on, decided by how the insert was cloned in.
    - `slot.strand_of_interest` -- which strand relative to the CASSETTE this
      slot's directional models care about. +1, the default, is the cassette's
      own sense strand; that is what a producer or target slot wants, because
      an internal polyA signal matters on the RNA that is actually made.

    Composing them is the only way a cassette-relative preference survives being
    cloned in backwards. `strand_of_interest` alone is meaningless without an
    orientation to measure it against, and this function previously returned it
    unchanged -- so `cassette_orientation` was a required field of the context
    that nothing in BT5 read, and the failure the docstring above warns about
    was live for every reverse-oriented cassette.

    On the common forward cassette the product is the identity, so nothing about
    a forward-oriented design changes.
    """
    composed = ctx.cassette_orientation * slot.strand_of_interest
    return 1 if composed == 1 else -1
