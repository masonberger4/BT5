"""Percentiles against the null -- the only number BT5 lets out.

`bt5.score.null` has held `null_distribution`, `percentile_of` and `normalise`
since M3 landed, and nothing on the design path called any of them; the walking
skeleton shipped every objective `unavailable` with the reason "ranking not
computed". This module is that wiring, and the honesty that has to survive it.

What a percentile here claims, exactly: *this design's raw score for this rule
sits above this fraction of random synonymous variants of the same protein,
spliced into the same backbone, at the same insertion site, scored by the same
rule.* It is not a prediction that the protein will express, and it is not
comparable across proteins or across backbones. Everything below exists to keep
that sentence true.

Three decisions worth reading twice.

**The null is built on the ASSEMBLED CONSTRUCT.** `build` is threaded down to
`null_distribution` rather than scoring bare CDSs, so a junction-spanning or
origin-spanning finding counts in the null exactly as it counts in the design.
`null.py`'s own docstring says a null that never contained a backbone makes the
report line simply false; this is where that is honoured or lost.

**One null, shared by every candidate.** The null is a distribution over random
synonymous variants of the protein in this construct. It does not depend on
which candidate is being ranked -- anchoring it on candidate A and then scoring
candidate B against it is not an approximation, it is the same distribution. So
it is built once and every gallery member is ranked against it, which is also
what makes the gallery's percentiles comparable to each other.

**Null size comes from `cost_class`, and each objective carries its own `n`.**
`Spec.cost_class` is documented as the thing that "drives null sampling", and
`ObjectiveScore.null_n` is per-objective, so this needs no new contract. A cheap
rule gets the shipped `DEFAULT_NULL_N`; a moderate one gets fewer, because
`build_gallery` and the null share one 10 s budget (PLAN's G7) and 200 k-mer
index builds do not fit in it. The variants are drawn from ONE seeded stream, so
a moderate objective's null is a strict PREFIX of a cheap objective's -- the same
variants, fewer of them -- rather than a second, differently-shaped sample.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass

from bt5.codon.tables import CodonUsage, NcbiGeneticCode
from bt5.core.context import DesignContext
from bt5.core.result import ObjectiveScore, ScoreCard
from bt5.core.services import Services
from bt5.core.spec import Breach, Evaluation, Spec
from bt5.core.types import Construct
from bt5.score.null import DEFAULT_NULL_N, NullDistribution, normalise, null_distribution

#: Where a rule puts the reason its objective could not be computed. `Evaluation`
#: has no "could not be computed" field -- `ObjectiveScore.unavailable` is the
#: type that carries one, and it is built downstream from this -- so b1 and c1
#: both return NaN plus a breach whose `detail` holds the sentence.
UNAVAILABLE_DETAIL_KEY = "unavailable_reason"

#: Variants per objective, by `Spec.cost_class`. `cheap` takes the shipped size;
#: `moderate` is cut because a moderate rule (a k-mer index, a repeat scan) costs
#: one to two orders of magnitude more per evaluation and the whole run is
#: budgeted at 10 s. A `moderate` null of 40 still resolves a percentile to
#: 1.25%, which is finer than anything BT5 claims to distinguish between two
#: designs.
NULL_N_BY_COST: Mapping[str, int] = {
    "cheap": DEFAULT_NULL_N,
    "moderate": 40,
    "expensive": 0,
}

#: What a `cost_class` this module has no size for gets. Treated as expensive --
#: reported unavailable rather than run -- because the failure of guessing the
#: other way is a 10 s gate blown by a rule nobody sized.
UNKNOWN_COST_N = 0


def synonymous_map(code: NcbiGeneticCode) -> Mapping[str, tuple[str, ...]]:
    """codon -> the codons encoding the same amino acid under this table.

    Built from `code.families()`, which already excludes stop codons from every
    synonymous set (tables 27 and 28 make TGA both Trp and a stop). A stop codon
    is therefore absent from this map entirely, so `synonymous_variant` passes
    the terminal stop through untouched and every null variant terminates where
    the design does.
    """
    out: dict[str, tuple[str, ...]] = {}
    for codons in code.families().values():
        for codon in codons:
            out[codon] = codons
    return out


def null_size(spec: Spec, overrides: Mapping[str, int] | None = None) -> int:
    """How many variants this objective's null gets."""
    if overrides is not None and spec.id in overrides:
        return int(overrides[spec.id])
    return NULL_N_BY_COST.get(spec.cost_class, UNKNOWN_COST_N)


@dataclass(frozen=True, slots=True)
class Nulls:
    """One null per scorable objective, plus why the others have none.

    `unavailable` is a mapping rather than a dropped set: an objective with no
    null is one the ranking does not account for, and a scorecard that simply
    omits it looks exactly like one where the rule was never configured.
    """

    by_spec: Mapping[str, NullDistribution]
    unavailable: Mapping[str, str]
    seed: int


def unavailability(evaluation: Evaluation) -> str | None:
    """None when this raw score is a real measurement; otherwise why it is not.

    A rule that could not compute its objective returns NaN and says why in a
    breach -- `c1_cai._unavailable` explains the choice, and the argument is
    load-bearing here: 0.0 would read as a catastrophically rare-codon sequence
    and the band's midpoint as a design exactly on target, so NaN is the only
    honest value. Which means NaN must never reach `percentile_of`. Every
    comparison against NaN is False, so it would score `better = 0`, `ties = 0`
    and emerge as a confident 0.0 percentile -- "worse than every one of 200
    random variants" -- about a quantity nobody measured. This is the guard
    against that, and it is why the scorecard reports the objective unavailable
    with the rule's own sentence instead.
    """
    if not math.isnan(evaluation.raw_score):
        return None
    for breach in evaluation.breaches:
        reason = breach.detail.get(UNAVAILABLE_DETAIL_KEY)
        if isinstance(reason, str) and reason:
            return reason
    return "the rule returned no measurable score for this construct"


def build_nulls(
    anchor_cds: str,
    specs: Sequence[Spec],
    *,
    anchor: Mapping[str, Evaluation],
    code: NcbiGeneticCode,
    ctx: DesignContext,
    svc: Services,
    build: Callable[[str], Construct],
    seed: int,
    usage: CodonUsage | None,
    sizes: Mapping[str, int] | None = None,
) -> Nulls:
    """A null per objective, all drawn from one seeded stream.

    `anchor` is the first candidate's evaluations, and it decides which
    objectives get a null at all. An objective the rule could not compute on the
    real construct cannot be computed on 200 variants of it either, so building
    that null would spend the budget producing 200 NaNs -- and worse, would put a
    NullDistribution with a NaN mean on the report.

    `usage` decides the null's KIND, and the difference is reported rather than
    smoothed over. With a host codon-usage table the variants are sampled at
    host frequency, which is the null a reader means by "better than 94% of
    random synonymous versions of this gene". Without one they are sampled
    uniformly over synonyms, which is a different and flatter distribution --
    `NullDistribution.kind` carries which, and `ObjectiveScore.null_kind` puts it
    on the report.
    """
    synonyms = synonymous_map(code)
    weights = dict(usage.w) if usage is not None else None
    kind = "host_frequency" if weights is not None else "uniform_synonymous"

    by_spec: dict[str, NullDistribution] = {}
    unavailable: dict[str, str] = {}
    for spec in specs:
        evaluation = anchor.get(spec.id)
        reason = unavailability(evaluation) if evaluation is not None else "the rule did not run"
        if reason is not None:
            unavailable[spec.id] = reason
            continue
        n = null_size(spec, sizes)
        if n < 2:
            unavailable[spec.id] = (
                f"no null was built: cost class {spec.cost_class!r} is too expensive "
                f"to sample {DEFAULT_NULL_N} times inside the design budget, so this "
                f"objective is not ranked"
            )
            continue
        by_spec[spec.id] = null_distribution(
            anchor_cds,
            synonyms=synonyms,
            build=build,
            score=_raw_scorer(spec, ctx, svc),
            seed=seed,
            n=n,
            kind=kind,  # type: ignore[arg-type]
            # Every rule in the scored set reads windows or counts motifs; none
            # calls `FoldEngine.mfe`, which its own docstring reserves for report
            # time. Asserted rather than assumed -- see `test_the_null_never_folds
            # _whole_transcripts`.
            windowed_fold_only=True,
            weights=weights,
        )
    return Nulls(by_spec=by_spec, unavailable=unavailable, seed=seed)


def _raw_scorer(spec: Spec, ctx: DesignContext, svc: Services) -> Callable[[Construct], float]:
    """`spec.evaluate(...).raw_score`, in the rule's own native units.

    The raw score, never the breach count: `percentile_of` reads `direction` and
    `band` to decide which way is better, and a count collapses a BAND rule like
    c1's 0.70-0.90 CAI target into "fewer is better", which is the max-CAI
    failure this project is organised around.
    """

    def score(construct: Construct) -> float:
        return spec.evaluate(construct, ctx, svc).raw_score

    return score


def score_candidate(
    evaluations: Mapping[str, Evaluation],
    specs: Sequence[Spec],
    *,
    nulls: Nulls,
    hard_checks: Sequence[Breach] = (),
    extra_unavailable: Mapping[str, tuple[str, str]] | None = None,
) -> ScoreCard:
    """One candidate's scorecard: a percentile per objective, or a reason.

    Takes evaluations already computed rather than a construct to evaluate, so
    the caller can run `RuleSet.findings()` ONCE per candidate and spend that
    single catalog pass on both the HARD_CHECK findings and the scored
    objectives. Evaluating each rule twice per candidate is the difference
    between a design inside PLAN's 10 s bar and one outside it.

    Availability is decided per CANDIDATE and not once for the run: a rule can
    compute its objective on one design and not on another, and inheriting the
    anchor's answer would put a percentile on a number this candidate does not
    have.

    `extra_unavailable` carries objectives that never reached `specs` at all --
    a rule dropped by `check_engine_calibration` because no folding engine is
    installed is not in the rule set, and would otherwise vanish from the report
    rather than being named. Keyed by spec id, valued `(unit, reason)`.
    """
    scores: list[ObjectiveScore] = []
    for spec in specs:
        evaluation = evaluations.get(spec.id)
        reason = (
            unavailability(evaluation)
            if evaluation is not None
            else "the rule did not run against this candidate"
        )
        null = nulls.by_spec.get(spec.id)
        if reason is None and null is None:
            reason = nulls.unavailable.get(spec.id, "no null was built for this objective")
        if reason is not None or null is None or evaluation is None:
            scores.append(
                ObjectiveScore.unavailable(
                    spec.id, spec.unit, reason or "this objective was not ranked"
                )
            )
            continue
        scores.append(
            normalise(
                spec_id=spec.id,
                raw=evaluation.raw_score,
                unit=spec.unit,
                direction=spec.direction,
                null=null,
                band=spec.band,
            )
        )
    for spec_id, (unit, reason) in sorted((extra_unavailable or {}).items()):
        scores.append(ObjectiveScore.unavailable(spec_id, unit, reason))

    return ScoreCard(
        scores=tuple(scores),
        hard_checks=tuple(hard_checks),
        total=weighted_total(scores, specs),
    )


def weighted_total(scores: Sequence[ObjectiveScore], specs: Sequence[Spec]) -> float:
    """The weighted sum, over PERCENTILES and over available objectives only.

    Percentiles are the only commensurable quantity here -- kcal/mol against CAI
    in [0, 1] against integer motif counts is not a sum -- which is what
    `ObjectiveScore`'s docstring means by "the ONLY quantity the weighted sum
    operates on".

    Renormalised over what was actually evaluated, so a run missing its
    highest-weight objective does not silently score lower than one that has it.
    The number is a rank key across the candidates of ONE run and nothing more:
    two runs with different objectives available are not comparable, and no part
    of BT5 compares them.
    """
    weights = {spec.id: float(spec.default_weight) for spec in specs}
    usable = [(s, weights.get(s.spec_id, 0.0)) for s in scores if s.available]
    total_weight = sum(w for _s, w in usable)
    if total_weight <= 0.0:
        return 0.0
    return sum(s.percentile * w for s, w in usable) / total_weight
