"""Fixtures for the wired solver.

`packages/engine/tests/rules/conftest.py` has equivalents, but pytest does not
share a conftest across sibling directories, and copying its `services` fixture
would copy its `_Tables` double -- which returns `dict[str, float]` from
`usage()` where the real provider returns a `CodonUsage`. These use the real
providers, because the point of this lane is that the real ones get constructed.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from Bio import SeqIO
from bt5.codon.tables import FileTableProvider, NcbiGeneticCode
from bt5.core.context import (
    BiosecurityVerdict,
    ContextSlot,
    DesignContext,
    HostId,
    Modality,
)
from bt5.core.services import Services
from bt5.core.types import (
    Construct,
    Interval,
    Segment,
    SegmentKind,
    Topology,
    TranslationUnit,
)
from bt5.vector.kmers import ConstructKmerIndex

ROOT = Path(__file__).resolve().parents[4]
BACKBONE = ROOT / "tests" / "data" / "backbones" / "synthetic_lenti_ef1a.gb"


@pytest.fixture(scope="session")
def code() -> NcbiGeneticCode:
    return FileTableProvider().genetic_code(11)


@pytest.fixture
def services() -> Services:
    """Real providers, no folding engine.

    `fold=None` rather than a stub: every threshold in BT5 is a kcal/mol number,
    so a stub returning plausible energies is the one thing the honesty
    apparatus exists to prevent. B1 reports itself unavailable and that is the
    behaviour under test everywhere else.
    """
    return Services(
        fold=None,
        kmer=ConstructKmerIndex,
        tables=FileTableProvider(),  # type: ignore[arg-type]
        rng=np.random.default_rng(42),
    )


def slot(
    role: str = "propagation",
    host: HostId = HostId.E_COLI_K12,
    modality: Modality = Modality.BACTERIAL_EXPRESSION,
    table: int = 11,
) -> ContextSlot:
    return ContextSlot(role, host, modality, table)  # type: ignore[arg-type]


def context(*slots: ContextSlot) -> DesignContext:
    return DesignContext(
        slots=slots or (slot(),),
        cassette_orientation=1,
        seed=42,
        screen=BiosecurityVerdict("not_run"),
    )


@pytest.fixture
def ctx() -> DesignContext:
    return context()


@pytest.fixture
def lentiviral_ctx() -> DesignContext:
    """Propagation in E. coli, packaging in HEK293. Two slots, deliberately.

    D4 and D6 resolve their enforcement per slot, so a single-slot context
    cannot exercise the resolution at all.
    """
    return context(
        slot(),
        slot("producer", HostId.HEK293, Modality.LENTIVIRAL, 1),
    )


def linear_cds(cds: str, protein: str, table_id: int = 11) -> Construct:
    """A bare designable CDS, no backbone. The simplest thing a rule can see."""
    return Construct(
        sequence=cds,
        topology=Topology.LINEAR,
        segments=(Segment(Interval(0, len(cds)), SegmentKind.DESIGNABLE_CDS, "cds"),),
        translation_units=(
            TranslationUnit(
                table_id,
                tuple(Interval(i, i + 3) for i in range(0, len(cds), 3)),
                protein,
                True,
            ),
        ),
    )


def with_backbone(cds: str, protein: str, left: str = "", right: str = "") -> Construct:
    """A CDS spliced between two immutable flanks, circular."""
    seq = left + cds + right
    return Construct(
        sequence=seq,
        topology=Topology.CIRCULAR,
        segments=(
            Segment(Interval(0, len(left)), SegmentKind.BACKBONE, "left"),
            Segment(Interval(len(left), len(left) + len(cds)), SegmentKind.DESIGNABLE_CDS, "cds"),
            Segment(Interval(len(left) + len(cds), len(seq)), SegmentKind.BACKBONE, "right"),
        ),
        translation_units=(
            TranslationUnit(
                11,
                tuple(Interval(i, i + 3) for i in range(len(left), len(left) + len(cds), 3)),
                protein,
                True,
            ),
        ),
    )


@pytest.fixture(scope="session")
def lenti_record() -> object:
    return SeqIO.read(BACKBONE, "genbank")
