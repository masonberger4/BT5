"""F1 -- exact direct repeats in the assembled plasmid.

The rule BT5 exists to get right, and the one most often stated wrongly. Two
thresholds are routinely conflated:

  RecBCD MEPS, 23-27 bp   the floor for RecA-DEPENDENT recombination
  ~200 bp                 below which deletion is RecA-INDEPENDENT

A `recA-` strain -- Stbl3, NEB Stable, standard for LVV and AAV work -- suppresses
only the first. The repeats codon choice creates or removes are 15-100 bp,
squarely inside the second, where the strain gives nothing: a 28 bp pair still
recombined at 7.8e-7 to 3.1e-5 in FOUR different recA- strains. So short-repeat
avoidance is BT5's job and the report must never tell a user the strain covers it.

Risk is therefore a surface over (length, spacer), not a length cutoff. The
RecA-independent route is strongly proximity-sensitive -- inserting sequence
between two copies suppresses it -- and the worst configuration is a TANDEM
repeat, where the copies touch and slipped-strand mispairing needs no loop at
all. `reca_strain_helps` is reported per finding so the protocol advice and the
sequence finding cannot drift apart.

This rule owns its own thresholds and citations rather than importing the vector
lane's `repeat_risk`. Rules never import another lane; they receive `Services`.
That is not ceremony here -- the bands ARE the scientific claim, and they belong
in the file that carries the citations for them.
"""

from __future__ import annotations

from collections.abc import Mapping
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
from bt5.rules.exempt import both_arms_exempt

#: Report nothing shorter. A chance 15-mer in a 5 kb plasmid has expected
#: occurrence ~1e-5, so anything found at this length is designed, not luck.
INFO_BP = 15
#: The shortest repeat vendors and Gibson protocols both treat as a liability.
WARN_BP = 20
#: RecBCD MEPS is 23-27 bp. At or above this a repeat is a hard finding.
HARD_BP = 25
#: Above this, recombination becomes RecA-dependent -- the regime a recA- strain
#: actually suppresses.
RECA_DEPENDENT_BP = 200

#: Edge-to-edge separations at which the RecA-independent route stays efficient.
NEAR_SPACER_BP = 100
FAR_SPACER_BP = 1000

#: A repeat at least this long is never "far enough apart to ignore".
SUBSTANTIAL_BP = 100

#: A tandem array must not flood the report with thousands of pairs.
MAX_FINDINGS = 200


def _spacer(first: Interval, second: Interval, length: int, circular: bool) -> int:
    """Edge-to-edge separation, the short way round on a circular construct.

    Proximity is the variable that decides whether the RecA-independent route
    fires, so measuring it the long way round on a plasmid would turn the most
    dangerous configuration into the safest-looking number in the report.
    """
    gap = second.start - first.end
    if not circular:
        return max(0, gap)
    return max(0, min(gap, length - (second.end - first.start)))


def risk_band(length: int, spacer: int, *, tandem: bool) -> str:
    """Classify on the (length, spacer) surface. Three bands, deliberately coarse.

    The literature supports the ORDERING of these regimes and that proximity
    matters. It does not support a calibrated rate, and inventing one would be
    the kind of number this project refuses to print.
    """
    if tandem:
        # Before the length floor: slipped-strand mispairing needs no loop, so a
        # short tandem array is real even under the vendor threshold.
        return "high" if length >= WARN_BP else "moderate"
    if length < WARN_BP:
        return "low"
    if length >= RECA_DEPENDENT_BP:
        return "moderate" if spacer > FAR_SPACER_BP else "high"
    if spacer <= NEAR_SPACER_BP:
        return "high"
    if spacer <= FAR_SPACER_BP:
        return "moderate"
    return "moderate" if length >= SUBSTANTIAL_BP else "low"


@register
class DirectRepeats:
    id: ClassVar[str] = "f1_direct_repeats"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "Exact direct repeats, banded on length and spacer"
    enforcement: ClassVar[Enforcement] = Enforcement.HARD_REPAIR
    evidence: ClassVar[Evidence] = Evidence.EVIDENCE_BACKED
    direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
    unit: ClassVar[str] = "weighted repeat pairs"
    band: ClassVar[tuple[float, float] | None] = None
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "A 28 bp direct repeat pair still recombined at 7.8e-7 to 3.1e-5 in FOUR "
            "different recA- strains -- the strain is a mitigation for LONG repeats, "
            "not a fix for the ones codon choice controls",
            "https://www.genoscope.cns.fr/MGE/pubs/Oliveira_Mol_Biotechnol_2008.pdf",
            2008,
            sign="supports",
        ),
        Citation(
            "Below ~200 bp deletion is RecA-INDEPENDENT (slipped-strand / single-strand "
            "annealing) and unaffected by recA/recF/recJ/recO; RecBCD MEPS is 23-27 bp",
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5426353/",
            2017,
            sign="supports",
        ),
        Citation(
            "The RecA-independent route is strongly proximity-sensitive: inserting "
            "sequence between two copies suppresses it, which is why risk is a surface "
            "over (length, spacer) rather than a length cutoff",
            "https://www.pnas.org/doi/10.1073/pnas.111008398",
            2001,
            sign="supports",
        ),
        Citation(
            "Repetitive 9-mers per 100 bp and longest repetitive sequence are the two "
            "highest-importance features of the published synthesis success model "
            "(random forest, 1,076 real vendor outcomes, F1 0.928) -- repeats outrank GC",
            "https://pubs.acs.org/doi/10.1021/acssynbio.9b00460",
            2020,
            sign="supports",
        ),
    )
    last_verified: ClassVar[str] = "2026-08-28"
    weight_provenance: ClassVar[str] = ""  # hard rule; never weighted
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.0
    #: Not expressible in the automaton -- deciding whether a codon completes a
    #: repeat needs the whole construct, not the last few bases -- so the DP is
    #: steered and Tier B repairs. Weighted high because repeats outrank GC in
    #: the only published model trained on real vendor outcomes.
    steering_weight: ClassVar[float] = 1.0
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.PAIRED_SEGMENTS
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "moderate"
    conflicts_with: ClassVar[tuple[str, ...]] = ()
    brief_ref: ClassVar[str] = "2.F1"
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "min_len": {"type": "integer", "default": INFO_BP, "minimum": 8},
            "hard_len": {"type": "integer", "default": HARD_BP, "minimum": 8},
        },
    }

    def __init__(self, min_len: int = INFO_BP, hard_len: int = HARD_BP) -> None:
        if min_len < 8:
            raise ValueError(
                f"min_len {min_len} is below the length at which a repeat is "
                f"distinguishable from chance; an 8-mer occurs ~0.05 times per 1.5 kb"
            )
        if hard_len < min_len:
            raise ValueError(f"hard_len {hard_len} must not be below min_len {min_len}")
        self.min_len = min_len
        self.hard_len = hard_len

    def gate(self, slot: ContextSlot) -> bool:
        # Every construct BT5 designs passes through a cloning host, whatever it
        # is finally delivered as, so this applies in every context.
        return True

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """None. A repeat is not decidable from a bounded suffix: whether a codon
        completes one depends on the whole construct, including the backbone.
        Tier A is steered, Tier B repairs, and the independent validator refuses."""
        return None

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        index = svc.kmer.of(c, self.min_len)
        breaches: list[Breach] = []
        worst = 0

        for first, second in index.duplicates(self.min_len):
            if len(breaches) >= MAX_FINDINGS:
                break
            if both_arms_exempt(c, first, second):
                continue
            length = first.length
            spacer = _spacer(first, second, c.length, c.is_circular)
            tandem = spacer == 0
            risk = risk_band(length, spacer, tandem=tandem)
            worst = max(worst, length)
            reca_helps = length >= RECA_DEPENDENT_BP
            # The span from the first copy's start to the second's end: the
            # region a repair has to work within, and PAIRED_SEGMENTS is the
            # localisation policy that says so.
            span = Interval(first.start, max(first.end, second.end))

            breaches.append(
                Breach(
                    spec_id=self.id,
                    interval=span,
                    magnitude={"low": 0.3, "moderate": 1.0, "high": 2.0}[risk]
                    * (2.0 if length >= self.hard_len else 1.0),
                    message=(
                        f"{length} bp exact direct repeat at {first.start} and "
                        f"{second.start}"
                        + (", TANDEM (copies touch)" if tandem else f", {spacer} bp apart")
                        + f"; {risk} risk. "
                        + (
                            "Long enough that a recA- strain suppresses the RecA-dependent route."
                            if reca_helps
                            else "A recA- strain does NOT cover this: below ~200 bp "
                            "deletion is RecA-independent."
                        )
                    ),
                    fixable_by_codon_choice=c.overlaps_editable(first)
                    or c.overlaps_editable(second),
                    detail={
                        "length": float(length),
                        "spacer": float(spacer),
                        "risk": risk,
                        "tandem": "yes" if tandem else "no",
                        "reca_strain_helps": "yes" if reca_helps else "no",
                    },
                )
            )

        return Evaluation(
            spec_id=self.id,
            passes=worst < self.hard_len,
            raw_score=float(sum(b.magnitude for b in breaches)),
            breaches=tuple(breaches),
            n_evaluated=c.length,
        )
