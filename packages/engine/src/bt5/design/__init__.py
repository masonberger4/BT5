"""M11 design lane: one protein into one backbone, ranked and verified end to end.

`design()` is the single public entry point. It resolves an insertion site
(multiple cloning site first), sweeps the objective simplex for a panel of
candidates, repairs and PROVES each one with the independent validator, ranks
them by percentile against a random-synonymous null built on the assembled
construct, and returns a `SkeletonResult` carrying the panel, its QC report, the
rendered summary, the annotated GenBank and the vendor order CSV.

What it still refuses to pretend. No number here is a predicted expression
level, titer, yield or fold-improvement -- percentiles say where a design sits
against its own synonymous null and nothing more. An objective that could not be
evaluated is reported `unavailable` with its reason rather than dropped, because
a scorecard missing its highest-weight term looks exactly like one where that
term was never configured. `native_baseline` is populated only from a real
wild-type CDS the caller supplies, never from a back-translation. And
`QcReport.is_complete` is True only when there is genuinely nothing left to
state -- a biosecurity screen that did not run keeps it False, which is M8's to
change and not this lane's.
"""

from bt5.design.errors import DesignError
from bt5.design.gallery import DEFAULT_SWEEP_STEPS, SolveSpace, sweep_designs
from bt5.design.order import entries_for, order_csv, write_order_plate
from bt5.design.ranking import Nulls, build_nulls, score_candidate
from bt5.design.runner import (
    DEFAULT_GALLERY_SIZE,
    UNSCREENED,
    DesignInputs,
    SkeletonResult,
    design,
)

__all__ = [
    "DEFAULT_GALLERY_SIZE",
    "DEFAULT_SWEEP_STEPS",
    "UNSCREENED",
    "DesignError",
    "DesignInputs",
    "Nulls",
    "SkeletonResult",
    "SolveSpace",
    "build_nulls",
    "design",
    "entries_for",
    "order_csv",
    "score_candidate",
    "sweep_designs",
    "write_order_plate",
]
