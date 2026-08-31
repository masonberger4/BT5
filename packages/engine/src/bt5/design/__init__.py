"""M11 design lane: one protein into one backbone, verified end to end.

The walking skeleton. `design()` is the single public entry point; it resolves an
insertion site (multiple cloning site first), assembles the construct, runs the
HARD_REPAIR rules through the solver, proves the result with the independent
validator, and returns a `SkeletonResult` carrying the design, its QC report, the
rendered summary and the exported GenBank.

What it refuses to pretend: it scores nothing (objectives ship `unavailable`), it
carries no baseline (`native_baseline` is None), and `QcReport.is_complete` is
always False. Those are the ranking increment's, and stating them is the point.
"""

from bt5.design.errors import DesignError
from bt5.design.runner import DesignInputs, SkeletonResult, design

__all__ = ["DesignError", "DesignInputs", "SkeletonResult", "design"]
