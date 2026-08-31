"""Assemble the `Services` bundle every rule is handed, and say what degraded.

The rules take one `Services` object -- a folding engine (or None), a k-mer index
CLASS, a table provider, and an explicitly-seeded RNG. This module builds it and,
just as importantly, reports what could not be provided, because a folding engine
that is absent or installed-but-uncalibrated is not an error to swallow: it is a
degradation the report must name, or a scorecard silently drops its
highest-weight objective and looks complete.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import cast

import numpy as np

from bt5.codon.tables import FileTableProvider
from bt5.core.services import FoldEngine, Services
from bt5.structure.vienna import degradation_reason, load_fold_engine
from bt5.vector.kmers import ConstructKmerIndex


def build_services(*, seed: int, fold: FoldEngine | None = None) -> Services:
    """The one bundle passed to every rule.

    `fold` is loaded from ViennaRNA when the caller does not supply one; it comes
    back None when ViennaRNA is not installed, which is not a failure -- the rule
    that needs it reports its objective unavailable rather than fabricating a ΔG.
    The RNG is always explicitly seeded; the global RNG is banned in `src/`.
    """
    # `load_fold_engine` returns the concrete ViennaFold, the repo's FoldEngine
    # implementation; mypy does not accept its `version` property as the
    # protocol's ClassVar, so the cast states what is true at runtime.
    engine = cast("FoldEngine | None", fold if fold is not None else load_fold_engine())
    return Services(
        fold=engine,
        kmer=ConstructKmerIndex,
        # FileTableProvider returns NcbiGeneticCode/CodonUsage where the protocol
        # names GeneticCode/Mapping; the same # type: ignore every Services call
        # site in the repo carries, not a new gap.
        tables=FileTableProvider(),  # type: ignore[arg-type]
        rng=np.random.default_rng(seed),
    )


def engine_versions(fold: FoldEngine | None) -> Mapping[str, str]:
    """The engine versions to record in provenance. Empty when no engine loaded."""
    if fold is None:
        return {}
    return {fold.name: fold.version}


def fold_degradation() -> str | None:
    """What to record because folding is unavailable or uncalibrated, or None.

    Delegated to the structure lane so the exact wording and the calibration
    check live with the engine, not copied here where they would drift.
    """
    return degradation_reason()
