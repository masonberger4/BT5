"""E1 -- homopolymer run limits, guaranteed by the automaton.

The asymmetry is the interesting part and it is not arbitrary. IDT allows ~10
consecutive A/T but only ~6 consecutive G/C, and that 10-vs-6 gap is the
clearest vendor evidence that G/C runs are chemically worse rather than merely
repetitive: G-runs aggregate into quadruplexes on the solid support during
phosphoramidite synthesis. So this rule carries two limits, not one, and a
single "max homopolymer" parameter would encode the wrong physics.

HARD_LATTICE, and this is the case that shows what that buys. A homopolymer is
decidable from the last few bases, so the Tier-A automaton makes it UNREACHABLE
-- including runs created across a codon boundary, at the CDS/backbone junction,
and spanning the origin, none of which a per-codon check would see.

Only the A and G runs are listed in `forbidden`. The solver closes the pattern
set under reverse complement, so poly-T and poly-C come free; listing all four
would forbid each run twice and double-count every hit.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from typing import ClassVar

from bt5.core.context import ContextSlot, DesignContext
from bt5.core.registry import register
from bt5.core.services import Services
from bt5.core.spec import (
    Breach,
    Citation,
    Direction,
    Enforcement,
    Evaluation,
    Evidence,
    LatticeTerms,
    LocalizationPolicy,
    RepairPolicy,
)
from bt5.core.types import Construct, Interval
from bt5.rules.vendors import DEFAULT_VENDOR, PROFILES, orderable, orderable_keys

#: The limits the default configuration carries, for the schema's advertised
#: defaults. Read from the profile rather than restated, so a corrected vendor
#: number cannot leave the documented default disagreeing with the enforced one.
_DEFAULT = PROFILES[DEFAULT_VENDOR]


def _maximal_runs(seq: str, circular: bool) -> Iterator[tuple[int, int, str]]:
    """(start, length, base) for every maximal single-base run.

    On a circular construct a run spanning the origin is ONE run, not a head and
    a tail, so the sequence is first rotated to a position that begins a run.
    Reporting the two halves separately would let a 16 nt run past a 9 nt limit
    as two compliant runs of 8 -- the exact case a linear scan cannot see and
    the reason the automaton is worth having.
    """
    n = len(seq)
    if n == 0:
        return

    offset = 0
    if circular and seq[0] == seq[-1]:
        while offset < n and seq[offset] == seq[0]:
            offset += 1
        if offset == n:  # the whole construct is one base
            yield 0, n, seq[0]
            return

    rotated = seq[offset:] + seq[:offset] if offset else seq
    i = 0
    while i < n:
        j = i
        while j < n and rotated[j] == rotated[i]:
            j += 1
        yield (i + offset) % n, j - i, rotated[i]
        i = j


@register
class Homopolymers:
    id: ClassVar[str] = "e1_homopolymers"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "Homopolymer run limits (A/T and G/C separately)"
    enforcement: ClassVar[Enforcement] = Enforcement.HARD_LATTICE
    evidence: ClassVar[Evidence] = Evidence.VENDOR_ASSERTED
    direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
    unit: ClassVar[str] = "runs over limit"
    band: ClassVar[tuple[float, float] | None] = None
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "IDT gBlocks complexity: ~10 consecutive A/T but only ~6 consecutive G/C. "
            "The asymmetry is the vendor evidence that G/C runs are chemically worse, "
            "not merely repetitive",
            "https://www.idtdna.com/pages/products/genes-and-gene-fragments/double-stranded-dna-fragments/gblocks-gene-fragments",
            2026,
            sign="supports",
        ),
        Citation(
            "Twist gene synthesis FAQ: homopolymer limit moved from 14 to 30 bp between "
            "2023 and 2026 -- the reason these numbers carry a last_verified date",
            "https://www.twistbioscience.com/faq/gene-synthesis",
            2026,
            sign="qualifies",
        ),
        Citation(
            "G-quadruplex formation fouls solid-phase synthesis and stalls replication "
            "forks; E. coli mutation rates span 5.5e-5 to 2.7e-10 per cell per "
            "generation across G4 variants",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC10530614/",
            2023,
            sign="supports",
        ),
    )
    last_verified: ClassVar[str] = "2026-08-28"
    weight_provenance: ClassVar[str] = ""  # hard rule; never weighted
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.0
    steering_weight: ClassVar[float] = 0.0  # unreachable by construction; steering is a no-op
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.MOTIF_LEN_MINUS_1
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "cheap"
    conflicts_with: ClassVar[tuple[str, ...]] = ()
    brief_ref: ClassVar[str] = "2.E1"
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "max_at_run": {
                "type": "integer",
                "default": _DEFAULT.homopolymer_at,
                "minimum": 3,
            },
            "max_gc_run": {
                "type": "integer",
                "default": _DEFAULT.homopolymer_gc,
                "minimum": 3,
            },
            "vendor": {
                "type": "string",
                "default": DEFAULT_VENDOR,
                # Orderable only: "no vendor chosen" has no run limits to preset,
                # and inventing some would be answering an unasked question.
                "enum": list(orderable_keys()),
                "description": "Preset limits. Explicit max_at_run/max_gc_run override it.",
            },
        },
    }

    def __init__(
        self,
        max_at_run: int | None = None,
        max_gc_run: int | None = None,
        vendor: str = DEFAULT_VENDOR,
    ) -> None:
        p = orderable(vendor)
        assert p.homopolymer_at is not None  # both guaranteed by `orderable`
        assert p.homopolymer_gc is not None
        self.vendor = vendor
        self.max_at_run = p.homopolymer_at if max_at_run is None else max_at_run
        self.max_gc_run = p.homopolymer_gc if max_gc_run is None else max_gc_run
        if self.max_at_run < 3 or self.max_gc_run < 3:
            raise ValueError(
                f"run limits below 3 forbid ordinary sequence: "
                f"{self.max_at_run=} {self.max_gc_run=}"
            )

    def gate(self, slot: ContextSlot) -> bool:
        # Applies wherever DNA is synthesised, which is every modality BT5
        # serves -- an IVT mRNA template is still ordered as DNA.
        return True

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms:
        """The shortest run that is already too long, for one base of each pair.

        Forbidding `A * (limit + 1)` makes every longer run unreachable too: a
        run of 12 contains a run of 10. Listing T and C as well would be
        redundant -- the solver closes the set under reverse complement -- and
        would report every hit twice.
        """
        return LatticeTerms(forbidden=("A" * (self.max_at_run + 1), "G" * (self.max_gc_run + 1)))

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        """Report MAXIMAL runs, one breach each.

        Not `find_motifs` on the lattice terms, for two reasons. A run of 12 A
        contains three overlapping matches of `A * 10`, so that route reports
        one physical problem three times and the conflict panel shows three
        findings for one run. And it would put this rule and invariant I6 on the
        same code path -- the rule finds runs its own way, I6 checks the
        forbidden motifs, and the two agreeing means something.
        """
        breaches: list[Breach] = []
        for start, length, base in _maximal_runs(c.sequence, c.is_circular):
            kind = "A/T" if base in "AT" else "G/C"
            limit = self.max_at_run if base in "AT" else self.max_gc_run
            if length <= limit:
                continue
            iv = Interval(start, start + length)
            breaches.append(
                Breach(
                    spec_id=self.id,
                    interval=iv,
                    # Grows with the overrun: a run one base over the limit is a
                    # surcharge, one at twice the limit is a synthesis failure.
                    magnitude=float(length - limit),
                    message=(
                        f"{kind} homopolymer of {length} nt at {start}, "
                        f"over the {limit} nt {self.vendor} limit"
                    ),
                    # A run the user's own backbone already carries is real and
                    # worth reporting, and no codon can shorten it.
                    fixable_by_codon_choice=c.overlaps_editable(iv),
                    detail={
                        "base_class": kind,
                        "run_length": float(length),
                        "limit": float(limit),
                        "vendor": self.vendor,
                    },
                )
            )
        return Evaluation(
            spec_id=self.id,
            passes=not breaches,
            raw_score=float(len(breaches)),
            breaches=tuple(breaches),
            n_evaluated=c.length,
        )
