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


def translate(cds: str, *, table_id: int = 1) -> str:
    """Protein for a CDS, dropping a terminal stop."""
    table = CodonTable.unambiguous_dna_by_id[table_id]
    out = []
    for i in range(0, len(cds), 3):
        codon = cds[i : i + 3]
        if codon in table.stop_codons:
            break
        out.append(table.forward_table[codon])
    return "".join(out)


def resynonymise(cds: str, *, table_id: int = 1) -> str:
    """A DIFFERENT nucleotide sequence encoding the SAME protein.

    This is what re-optimisation actually does, and the case where annotation
    inside the CDS still describes exactly the residues it described before.
    """
    table = CodonTable.unambiguous_dna_by_id[table_id]
    synonyms: dict[str, list[str]] = {}
    for codon, aa in sorted(table.forward_table.items()):
        synonyms.setdefault(aa, []).append(codon)
    out = []
    for i in range(0, len(cds), 3):
        codon = cds[i : i + 3]
        if codon in table.stop_codons:
            out.append(codon)
            continue
        options = synonyms[table.forward_table[codon]]
        # pick a different codon where one exists, so the test is not vacuous
        out.append(next((c for c in options if c != codon), codon))
    return "".join(out)
