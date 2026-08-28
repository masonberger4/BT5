"""B1 -- folding free energy of the Kudla window, and the only objective whose
published effect size earns a full weight.

The evidence, and why it outranks everything else in the catalog. Kudla 2009
built 154 synonymous GFP variants spanning a 250-fold expression range and found
that the folding free energy of a -4..+37 window around the start codon explained
**44% of the variance** (r = 0.66), and 59% in a second promoter system. In the
same dataset CAI gave r = 0.14, not significant, and whole-mRNA MFE r = 0.16,
also not significant. Nothing else BT5 computes comes close, which is why this
rule carries weight 1.0 while the codon-composition rules carry 0.2-0.3.

**The window spans the UTR/CDS junction, so it cannot be computed from the CDS
alone.** Four bases of it are the vector's 5' leader. That single fact is the
reason `Construct` exists, the reason gate G6 exists, and the reason this rule
has three ways to decline:

  * no folding engine installed;
  * the CDS sits too close to a linear construct's end for the window to fit;
  * **no annotated 5'UTR covering those four bases.**

The third is the one worth arguing about, because the tempting alternative is to
fold whatever backbone precedes the CDS. That backbone is often a promoter, which
is not transcribed and is not in the mRNA -- folding it and reporting the result
as the leader produces a number that looks entirely reasonable and is measuring a
molecule that does not exist. docs/PLAN.md is explicit: "with no annotated 5'UTR
the app degrades honestly and says the 5'-structure objective is unavailable -- it
must never silently fold the CDS alone and report a number as if the UTR were
there." An unavailable highest-weight objective is a worse-looking report and a
truer one.

`upstream=0` is refused at construction for the same reason: it would turn this
into a CDS-only fold that still called itself B1.

**Why `engine_calibration` is set here and is None on e5_synthesis_repeats.**
This field names the FOLD engine a rule's thresholds are measured on, and B1's
quantity really is a ViennaRNA kcal/mol. Declaring it means `check_engine_calibration`
refuses the run outright if some other engine is active -- the comparison would
succeed and mean nothing -- and, with no engine at all, returns this rule as
unrunnable so M3 reports `ObjectiveScore.unavailable` rather than dropping the
term. E5's Tm is not a folding energy, so declaring it there would have raised
against ViennaRNA on every run.

**Why `steering_weight` is 0.0 even though this is the objective that matters
most.** Tier A decides from a bounded suffix of the emitted sequence; a fold is
not that. Windowed GC gets a steering term because a running G+C count is cheap
and is *the same quantity* the rule scores. There is no such proxy for structure.
Early-codon GC or A-richness would steer on something B1 does not measure -- that
is B4's rule wearing B1's name -- so the honest arrangement is Tier C windowed
polish over the real window, and a steering term of zero.

**Gated to bacterial expression, on modality rather than host.** Kudla measured
Shine-Dalgarno-driven initiation in E. coli. Eukaryotic scanning initiation makes
cap-proximal structure the analogous term, which is B11's, with its own
position-dependent ladder; applying a -4..+37 bacterial window there is applying
the wrong ruler. And the gate reads `modality`, not `host`, because the E. coli
slot of a lentiviral job is a PROPAGATION slot -- the plasmid is maintained there,
never expressed -- so gating on "host is E. coli" would score a translation event
that does not happen.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import ClassVar

from bt5.core.context import ContextSlot, DesignContext, Modality
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
    strand_for,
)
from bt5.core.types import Construct, Interval

#: Kudla 2009's window, relative to the A of ATG. Duplicated from
#: bt5.structure.windows rather than imported: a rule never imports another
#: lane, it receives Services. The numbers are the paper's, not M6's.
KUDLA_UPSTREAM = 4
KUDLA_DOWNSTREAM = 37

#: GenBank spells the same feature both ways depending on the writer. The vector
#: lane's `utr_context` accepts the same pair; this is the GenBank vocabulary,
#: not that lane's invention, so agreeing on it is not a hidden dependency.
FIVE_PRIME_UTR_KINDS = ("5'UTR", "five_prime_UTR")

#: The engine and parameter set every threshold in this file is measured on.
ENGINE_CALIBRATION = "viennarna:rna_turner2004"


def leader_of(window: Interval, upstream: int) -> Interval:
    """The part of the window that is 5'UTR rather than coding sequence.

    On the minus strand the transcript's 5' end is at HIGHER construct
    coordinates, so the leader is the tail of the interval, not its head.
    Getting this backwards checks the wrong four bases for UTR annotation and
    then folds a window whose leader is really the far end of the CDS.
    """
    if window.strand == 1:
        return Interval(window.start, window.start + upstream, 1)
    return Interval(window.end - upstream, window.end, -1)


@register
class FivePrimeFolding:
    id: ClassVar[str] = "b1_five_prime"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "5' folding free energy over the Kudla -4..+37 window"
    enforcement: ClassVar[Enforcement] = Enforcement.SOFT
    evidence: ClassVar[Evidence] = Evidence.EVIDENCE_BACKED
    #: Less negative is less structure is better initiation. NOT a band: unlike
    #: CAI or GC there is no measured optimum in the middle, and Kudla's
    #: relationship is monotone across the whole 250-fold range they observed.
    direction: ClassVar[Direction] = Direction.HIGHER_IS_BETTER
    unit: ClassVar[str] = "kcal/mol"
    band: ClassVar[tuple[float, float] | None] = None
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "154 synonymous GFP variants over a 250-fold expression range: folding "
            "free energy of the -4..+37 window explained 44% of variance (r = 0.66), "
            "59% in a second promoter system, while CAI gave r = 0.14 (not "
            "significant) and whole-mRNA MFE r = 0.16 (not significant)",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC3902468/",
            2009,
            sign="supports",
        ),
        Citation(
            "'CAI has no value in predicting gene expression'; a deliberately "
            "high-CAI control expressed at ~15% of the best variant -- the negative "
            "result that makes a structure term rather than a codon term the "
            "highest-weight objective",
            "https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0007002",
            2009,
            sign="supports",
        ),
        Citation(
            "Boel 2016, 6,348 genes: codon content 3-5x more influential than "
            "structure on a soluble-protein readout under T7. Ships as a stated "
            "disagreement rather than a resolved one -- different promoter, "
            "different readout, and BT5 cannot adjudicate it",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC5054687/",
            2016,
            sign="contradicts",
        ),
        Citation(
            "All computable design features together explain 5-31% (mean ~14%) of "
            "protein-level variance over 244,000 designed sequences -- the ceiling "
            "this objective sits under, and the reason its output is a percentile "
            "against a null rather than a predicted expression level",
            "https://www.nature.com/articles/nbt.4238",
            2018,
            sign="qualifies",
        ),
    )
    last_verified: ClassVar[str] = "2026-08-28"
    weight_provenance: ClassVar[str] = (
        "1.0, the highest in the catalog, and the only weight in BT5 justified by a "
        "measured effect size rather than a ranking: r = 0.66 over 154 variants, 44% "
        "of expression variance, against CAI's non-significant r = 0.14 in the same "
        "dataset. It is weighted above the repeat family because those rules are "
        "backed by feature-importance orderings and vendor thresholds, which say "
        "WHICH factors matter but not HOW MUCH. Held at 1.0 rather than higher "
        "because Cambray's 244,000-sequence factorial puts the ceiling for all "
        "computable features together at ~14% of protein-level variance, and because "
        "Boel 2016 measured the opposite ordering on a different readout -- the "
        "disagreement is unresolved and a weight above every other objective "
        "combined would be asserting it settled."
    )
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 1.0
    #: Zero. See the module docstring: Tier A cannot evaluate a fold, and a cheap
    #: proxy for one would steer on a different quantity than this rule scores.
    steering_weight: ClassVar[float] = 0.0
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.WINDOW_MINUS_1
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "moderate"
    conflicts_with: ClassVar[tuple[str, ...]] = ()
    brief_ref: ClassVar[str] = "2.B1"
    engine_calibration: ClassVar[str | None] = ENGINE_CALIBRATION
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "upstream": {"type": "integer", "default": KUDLA_UPSTREAM, "minimum": 1},
            "downstream": {"type": "integer", "default": KUDLA_DOWNSTREAM, "minimum": 1},
        },
    }

    def __init__(self, upstream: int = KUDLA_UPSTREAM, downstream: int = KUDLA_DOWNSTREAM) -> None:
        if upstream < 1:
            raise ValueError(
                f"upstream={upstream} would fold the CDS alone and still call itself "
                f"B1. The measured window spans the UTR/CDS junction; a window that "
                f"does not is a different quantity with no published relationship "
                f"to expression"
            )
        if downstream < 1:
            raise ValueError(f"downstream={downstream} must reach into the CDS")
        self.upstream = upstream
        self.downstream = downstream

    def gate(self, slot: ContextSlot) -> bool:
        return slot.modality is Modality.BACTERIAL_EXPRESSION

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """None. A folding free energy is not decidable from a bounded suffix."""
        return None

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        slots = [s for s in ctx.active_slots if self.gate(s)]
        if not slots:
            return self._unavailable(
                c, "no bacterial expression slot in this context; B1 is E. coli only"
            )
        if svc.fold is None:
            return self._unavailable(
                c,
                "no folding engine installed, so the 5' folding objective cannot be "
                "computed; install ViennaRNA to enable it",
            )

        editable = sorted(c.editable)
        if not editable:
            return self._unavailable(c, "no designable CDS to take a 5' window from")

        measured: list[tuple[Interval, float]] = []
        for slot in slots:
            strand = strand_for(ctx, slot)
            cds = editable[0] if strand == 1 else editable[-1]
            window = self._window(Interval(cds.start, cds.end, strand), c)
            if window is None:
                return self._unavailable(
                    c,
                    f"the CDS starts too close to the end of this linear construct "
                    f"for the -{self.upstream}..+{self.downstream} window to fit",
                )
            leader = leader_of(window, self.upstream)
            if not self._annotated_leader(c, leader):
                return self._unavailable(
                    c,
                    f"no annotated 5'UTR covering the {self.upstream} bases upstream "
                    f"of the start codon. The measured window spans the UTR/CDS "
                    f"junction, and the sequence there may be promoter rather than "
                    f"leader -- folding it would report a molecule that is never "
                    f"transcribed. Annotate the 5'UTR in your map and re-run",
                    interval=leader,
                )
            measured.append((window, svc.fold.mfe_window(c.sequence, window).dg_kcal_mol))

        # The MOST NEGATIVE across gated slots: the most structured 5' end is the
        # binding one, and averaging would let a permissive slot hide it.
        worst = min(dg for _, dg in measured)
        return Evaluation(
            spec_id=self.id,
            passes=True,  # an objective, not a constraint: there is no published
            # cutoff for this window, only a monotone relationship
            raw_score=worst,
            windows=tuple(measured),
            n_evaluated=self.upstream + self.downstream,
        )

    def _window(self, cds: Interval, c: Construct) -> Interval | None:
        """`bt5.structure.windows.five_prime_window`, reimplemented here.

        Rules receive Services and never import another lane, and window geometry
        is not on the Services protocol. The arithmetic is small; what matters is
        that both copies return None rather than a clamped window when the leader
        does not exist, because a window silently shortened to what fits is a
        DIFFERENT window and not comparable to one that had its full leader.
        """
        if cds.length < self.downstream:
            return None
        span = self.upstream + self.downstream
        start = cds.start - self.upstream if cds.strand == 1 else cds.end - self.downstream
        if start < 0:
            if not c.is_circular:
                return None
            start += c.length
        if not c.is_circular and start + span > c.length:
            return None
        return Interval(start, start + span, cds.strand)

    def _annotated_leader(self, c: Construct, leader: Interval) -> bool:
        return any(
            f.interval.contains(leader, c.length, c.is_circular)
            for f in c.features
            if f.kind in FIVE_PRIME_UTR_KINDS
        )

    def _unavailable(
        self, c: Construct, reason: str, *, interval: Interval | None = None
    ) -> Evaluation:
        """NaN plus a breach carrying the reason.

        `Evaluation` has no field for "could not be computed" -- `ObjectiveScore.
        unavailable` is M3's type and is built downstream from this -- so the
        reason travels in the one string channel a rule has. NaN rather than 0.0
        because 0 kcal/mol is a real, meaningful value for this quantity: an
        unstructured 5' end, the best possible score. Reporting it for "we could
        not measure" would put the highest-weight objective at the top of the
        ranking precisely when it is unknown.

        B1 is the first rule to need this. If a second one does, a dedicated
        field on `Evaluation` is the right answer and is a MINOR contract change.
        """
        where = interval or (sorted(c.editable)[0] if c.editable else Interval(0, 1))
        return Evaluation(
            spec_id=self.id,
            passes=True,
            raw_score=float("nan"),
            breaches=(
                Breach(
                    spec_id=self.id,
                    interval=where,
                    magnitude=0.0,
                    message=f"5' folding objective unavailable: {reason}",
                    # Nothing about the CODING sequence causes this, and nothing
                    # about it fixes it -- the missing input is the vector's
                    # annotation or the engine install.
                    fixable_by_codon_choice=False,
                    detail={"unavailable_reason": reason, "calibration": ENGINE_CALIBRATION},
                ),
            ),
            n_evaluated=0,
        )
