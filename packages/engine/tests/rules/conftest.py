"""Shared construction helpers for the rules lane.

Every rule is evaluated against a `Construct` and a `DesignContext`, never a
bare string -- that is the decision that makes junction-spanning and
origin-spanning hits impossible to miss. These helpers make building both cheap
so a rule test never has a reason to reach for a shortcut that would evaluate
the CDS in isolation.
"""

from __future__ import annotations

import numpy as np
import pytest
from bt5.core.context import (
    BiosecurityVerdict,
    ContextSlot,
    DesignContext,
    HostId,
    Modality,
)
from bt5.core.services import Services
from bt5.core.types import Construct, Interval, Segment, SegmentKind, Strand, Topology


def construct(
    cds: str,
    backbone: str = "",
    *,
    circular: bool = True,
    cds_start: int = 0,
) -> Construct:
    """A construct whose CDS is `cds`, optionally offset into a backbone.

    `cds_start` past the end of the sequence is how a CDS spanning the origin is
    expressed: the segment wraps, exactly as an insert cloned across position 0
    would be stored.
    """
    seq = cds + backbone if cds_start == 0 else backbone[:cds_start] + cds + backbone[cds_start:]
    end = cds_start + len(cds)
    segments = [Segment(Interval(cds_start, end), SegmentKind.DESIGNABLE_CDS, "cds")]
    if len(seq) > len(cds):
        lo, hi = (end, len(seq)) if cds_start == 0 else (0, cds_start)
        if hi > lo:
            segments.append(Segment(Interval(lo, hi), SegmentKind.BACKBONE, "vector"))
    return Construct(
        sequence=seq,
        topology=Topology.CIRCULAR if circular else Topology.LINEAR,
        segments=tuple(segments),
    )


def wrapping_construct(prefix: str, cds_head: str, cds_tail: str) -> Construct:
    """A construct whose CDS runs off the end and resumes at position 0.

    `prefix` is backbone, then `cds_tail` closes the sequence and `cds_head`
    opens it, so the CDS segment is stored as one interval with end > length.
    """
    seq = cds_head + prefix + cds_tail
    start = len(cds_head) + len(prefix)
    return Construct(
        sequence=seq,
        topology=Topology.CIRCULAR,
        segments=(
            Segment(
                Interval(start, start + len(cds_tail) + len(cds_head)),
                SegmentKind.DESIGNABLE_CDS,
                "cds",
            ),
            Segment(Interval(len(cds_head), start), SegmentKind.BACKBONE, "vector"),
        ),
    )


def slot(
    role: str = "producer",
    host: HostId = HostId.HEK293,
    modality: Modality = Modality.LENTIVIRAL,
    table: int = 1,
    strand_of_interest: Strand = 1,
) -> ContextSlot:
    return ContextSlot(role, host, modality, table, strand_of_interest)  # type: ignore[arg-type]


def context(*slots: ContextSlot, cassette_orientation: Strand = 1) -> DesignContext:
    return DesignContext(
        slots=slots or (slot(),),
        cassette_orientation=cassette_orientation,
        seed=42,
        screen=BiosecurityVerdict("not_run"),
    )


@pytest.fixture
def services() -> Services:
    """Services with no folding engine.

    None rather than a stub: a rule that needs folding must report its objective
    unavailable, and a stub returning a plausible number is the one thing the
    whole honesty apparatus exists to prevent.
    """
    from bt5.vector.kmers import ConstructKmerIndex

    class _Tables:
        def genetic_code(self, table_id: int) -> object: ...
        def usage(self, host: str) -> dict[str, float]:
            return {}

        def weights(self, host: str, kind: str) -> dict[str, float]:
            return {}

    return Services(
        fold=None,
        kmer=ConstructKmerIndex,  # type: ignore[arg-type]
        tables=_Tables(),  # type: ignore[arg-type]
        rng=np.random.default_rng(42),
    )
