"""The empirical null, and the percentile that is the entire UI promise.

BT5 never reports a predicted expression level. It reports where a design sits
against a distribution of random synonymous variants of the SAME protein, in the
SAME construct, scored the SAME way. That is a claim it can actually support, and
it is the only reason a number reaches the user at all.

Three properties make the percentile mean what it says, and each is enforced here
rather than left to the caller's diligence:

The null is built on the ASSEMBLED CONSTRUCT. `null_distribution` takes a
`build` callable rather than a bare sequence, so a variant is spliced into the
real backbone before it is scored. Scoring bare CDSs and calling the result a
percentile measures against a distribution that never contained a backbone, and
the report line is then simply false.

Direction is honoured. "Better" is not "larger". A rule whose target is a BAND --
CAI's 0.70-0.90, GC 40-60% -- is scored on distance to the band, because ranking
CAI as higher-is-better drives it to 1.0, and max-CAI collapsing to one codon per
amino acid is the documented failure this whole project is organised around.

The RNG is seeded explicitly, always. An unseeded null makes a percentile
irreproducible, which makes the number it produced unfalsifiable.
"""

from __future__ import annotations

import math
from bisect import bisect_right
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np

from bt5.core.result import ObjectiveScore
from bt5.core.spec import Direction
from bt5.core.types import Construct

NullKind = Literal["host_frequency", "uniform_synonymous"]

#: docs/PLAN.md sizes the null at 200-500 variants. 200 is the floor at which a
#: percentile has a resolution of half a percent, which is finer than anything
#: BT5 claims to distinguish.
DEFAULT_NULL_N = 200


@dataclass(frozen=True, slots=True)
class NullDistribution:
    """Scores of N random synonymous variants, and how they were produced."""

    values: tuple[float, ...]
    kind: NullKind
    #: Whether every score in here came from windowed folding only. Carried
    #: rather than defaulted: a null built with whole-transcript MFE is a
    #: different and far more expensive object, and one silently mixed with
    #: windowed scores is meaningless.
    windowed_fold_only: bool
    seed: int

    @property
    def n(self) -> int:
        return len(self.values)

    @property
    def mean(self) -> float:
        return sum(self.values) / self.n if self.values else 0.0

    @property
    def sd(self) -> float:
        if self.n < 2:
            return 0.0
        mu = self.mean
        return math.sqrt(sum((v - mu) ** 2 for v in self.values) / (self.n - 1))


def synonymous_variant(
    cds: str,
    synonyms: Mapping[str, Sequence[str]],
    rng: np.random.Generator,
    *,
    weights: Mapping[str, float] | None = None,
) -> str:
    """One random synonymous recoding, codon by codon.

    A codon the map offers no alternative for passes through untouched, so the
    variant encodes exactly the same protein under exactly the same table --
    which is what makes it a null for THIS design rather than a different one.
    Whether the terminal stop varies is the CALLER's choice, expressed by
    whether the map lists the other stops as its synonyms; both are defensible
    and the map is where that decision belongs.

    An invariant position consumes no randomness. That is deliberate rather than
    incidental: it means the stream depends only on positions that actually have
    a choice, so inserting a Met into a protein does not silently reshuffle every
    codon downstream of it and change the whole null.

    Building the weight table costs one pass over `synonyms`, so a caller drawing
    many variants from one table should use `null_distribution`, which hoists it
    out of the loop.
    """
    return _variant(cds, synonyms, rng, table=weight_table(synonyms, weights))


def weight_table(
    synonyms: Mapping[str, Sequence[str]],
    weights: Mapping[str, float] | None,
) -> dict[str, list[float] | None] | None:
    """Per-codon cumulative weights, normalised to end at 1.0.

    Hoisted out of the draw because it depends only on the table and the host,
    never on the variant: rebuilding it per codon per variant is what made the
    weighted path ~15x the cost of the uniform one, and that went unnoticed for
    as long as no host resolved to a table and the path never ran (#98).

    A codon whose options carry no weight at all maps to None, which the draw
    reads as "fall back to uniform over those options". That is the same
    behaviour as the previous `total <= 0` guard: a family absent from the host's
    table is unknown, not forbidden, and silently never emitting its codons would
    quietly remove them from the null's support.
    """
    if weights is None:
        return None
    table: dict[str, list[float] | None] = {}
    for codon, options in synonyms.items():
        if not options or len(options) == 1:
            continue
        running = 0.0
        cumulative: list[float] = []
        for option in options:
            running += max(0.0, weights.get(option, 0.0))
            cumulative.append(running)
        table[codon] = [c / running for c in cumulative] if running > 0 else None
    return table


def _variant(
    cds: str,
    synonyms: Mapping[str, Sequence[str]],
    rng: np.random.Generator,
    *,
    table: dict[str, list[float] | None] | None,
) -> str:
    """One variant against a prepared `weight_table`; None means uniform.

    The weighted draw is an inverse-CDF lookup on one `rng.random()` rather than
    `rng.choice(p=...)`. Same distribution, ~23x cheaper at these sizes, where
    `choice`'s argument validation and CDF construction dwarf a 2-6 element pick.
    `bisect_right` is capped at the last index so a cumulative that floating-point
    division left a hair under 1.0 cannot index past the end.
    """
    out: list[str] = []
    for i in range(0, len(cds), 3):
        codon = cds[i : i + 3]
        options = synonyms.get(codon)
        if not options or len(options) == 1:
            out.append(codon)
            continue
        cumulative = table.get(codon) if table is not None else None
        if cumulative is None:
            out.append(options[int(rng.integers(0, len(options)))])
            continue
        out.append(options[bisect_right(cumulative, float(rng.random()), hi=len(cumulative) - 1)])
    return "".join(out)


def null_distribution(
    cds: str,
    *,
    synonyms: Mapping[str, Sequence[str]],
    build: Callable[[str], Construct],
    score: Callable[[Construct], float],
    seed: int,
    n: int = DEFAULT_NULL_N,
    kind: NullKind = "host_frequency",
    windowed_fold_only: bool,
    weights: Mapping[str, float] | None = None,
) -> NullDistribution:
    """Score `n` random synonymous variants of `cds` in the assembled construct.

    `build` is what forces the null into the right context. It takes a candidate
    CDS and returns the construct that CDS produces -- backbone, junctions,
    origin and all -- so every null score is computed against the same
    surroundings as the design being ranked.

    `windowed_fold_only` has no default. A caller that has not thought about
    whether its scorer folds whole transcripts has not earned a percentile.
    """
    if n < 2:
        raise ValueError(f"a null of {n} variants cannot support a percentile")
    if kind == "host_frequency" and weights is None:
        raise ValueError(
            "null_kind='host_frequency' needs codon weights; pass weights=, or "
            "declare kind='uniform_synonymous' if uniform sampling is what you mean"
        )
    rng = np.random.default_rng(seed)
    picked = weights if kind == "host_frequency" else None
    # Built ONCE, not once per variant: `n` is 200 by default and the table is
    # the same for every draw, so hoisting it is the difference between one pass
    # over `synonyms` and 200 * len(cds)/3 of them.
    table = weight_table(synonyms, picked)
    values = tuple(score(build(_variant(cds, synonyms, rng, table=table))) for _ in range(n))
    return NullDistribution(
        values=values, kind=kind, windowed_fold_only=windowed_fold_only, seed=seed
    )


def band_deviation(raw: float, band: tuple[float, float]) -> float:
    """Distance outside a target band; 0.0 inside it.

    This is what makes BAND representable at all. Collapsing a two-sided target
    into "higher is better" is what drives CAI to 1.0 -- and one codon per amino
    acid is a perfect nucleotide repeat, which is the failure mode the evidence
    is clearest about.
    """
    low, high = band
    if low > high:
        raise ValueError(f"band is inverted: {band}")
    if raw < low:
        return low - raw
    if raw > high:
        return raw - high
    return 0.0


def percentile_of(
    raw: float,
    null: NullDistribution,
    direction: Direction,
    band: tuple[float, float] | None = None,
) -> float:
    """Fraction of the null this score beats, in [0, 1].

    Ties count half. A variant scoring exactly what the design scores is neither
    beaten nor beating, and awarding it wholly to one side turns a null that
    happens to contain the design itself into a 1.0 or a 0.0.
    """
    if direction is Direction.BAND:
        if band is None:
            raise ValueError("a BAND objective cannot be normalised without its band")
        mine = band_deviation(raw, band)
        theirs = [band_deviation(v, band) for v in null.values]
        better = sum(1 for v in theirs if mine < v)
        ties = sum(1 for v in theirs if mine == v)
    elif direction is Direction.HIGHER_IS_BETTER:
        better = sum(1 for v in null.values if raw > v)
        ties = sum(1 for v in null.values if raw == v)
    else:
        better = sum(1 for v in null.values if raw < v)
        ties = sum(1 for v in null.values if raw == v)
    return (better + 0.5 * ties) / null.n


def normalise(
    *,
    spec_id: str,
    raw: float,
    unit: str,
    direction: Direction,
    null: NullDistribution,
    band: tuple[float, float] | None = None,
) -> ObjectiveScore:
    """Turn a raw score plus its null into the only quantity sliders may sum.

    Raw scores are incommensurable -- kcal/mol against CAI in [0,1] against
    integer motif counts -- so a weighted sum over them leaves most sliders dead
    across most of their range while one term quietly dominates.
    """
    return ObjectiveScore(
        spec_id=spec_id,
        raw=raw,
        unit=unit,
        percentile=percentile_of(raw, null, direction, band),
        null_n=null.n,
        null_mean=null.mean,
        null_sd=null.sd,
        null_kind=null.kind,
        windowed_fold_only=null.windowed_fold_only,
    )
