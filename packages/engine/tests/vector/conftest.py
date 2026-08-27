"""Shared fixtures for the vector lane.

The backbone fixture is a hand-specified GenBank file, NOT one produced by
bt5.vector's own writer. If the round-trip test read a file this lane had
written, it would prove only that the writer agrees with itself.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from Bio.Data import CodonTable
from bt5.vector import VectorBackbone, read_genbank

ROOT = Path(__file__).resolve().parents[4]
BACKBONE_PATH = ROOT / "tests" / "data" / "backbones" / "synthetic_lenti_ef1a.gb"


@pytest.fixture
def backbone_path() -> Path:
    return BACKBONE_PATH


@pytest.fixture
def backbone() -> VectorBackbone:
    return read_genbank(BACKBONE_PATH)


def make_cds(n_codons: int, *, seed: int = 11, table_id: int = 1) -> tuple[str, str]:
    """A valid ATG...stop CDS of `n_codons` codons, and the protein it encodes.

    Seeded explicitly: a global RNG here would make a failing example
    irreproducible, which is the whole point of the ban in CLAUDE.md.
    """
    table = CodonTable.unambiguous_dna_by_id[table_id]
    sense = sorted(table.forward_table)
    rng = np.random.default_rng(seed)
    picks = rng.integers(0, len(sense), n_codons - 2)
    cds = "ATG" + "".join(sense[i] for i in picks) + table.stop_codons[0]
    protein = "".join(table.forward_table[cds[i : i + 3]] for i in range(0, len(cds) - 3, 3))
    return cds, protein
