"""Shared fixtures for the score lane."""

from __future__ import annotations

import numpy as np
import pytest
from Bio.Data import CodonTable


@pytest.fixture(scope="session")
def synonyms() -> dict[str, list[str]]:
    """codon -> every codon encoding the same amino acid, table 11."""
    table = CodonTable.unambiguous_dna_by_id[11]
    by_aa: dict[str, list[str]] = {}
    for codon, aa in sorted(table.forward_table.items()):
        by_aa.setdefault(aa, []).append(codon)
    out = {c: by_aa[aa] for c, aa in table.forward_table.items()}
    for stop in table.stop_codons:
        out[stop] = list(table.stop_codons)
    return out


@pytest.fixture(scope="session")
def usage() -> dict[str, float]:
    """A deliberately skewed usage table, so weighting is observable."""
    table = CodonTable.unambiguous_dna_by_id[11]
    return {c: (10.0 if c.endswith("C") else 1.0) for c in table.forward_table}


def translate(cds: str, table_id: int = 11) -> str:
    table = CodonTable.unambiguous_dna_by_id[table_id]
    out = []
    for i in range(0, len(cds), 3):
        codon = cds[i : i + 3]
        if codon in table.stop_codons:
            break
        out.append(table.forward_table[codon])
    return "".join(out)


def make_cds(n_codons: int, *, seed: int = 3, table_id: int = 11) -> str:
    table = CodonTable.unambiguous_dna_by_id[table_id]
    sense = sorted(table.forward_table)
    rng = np.random.default_rng(seed)
    picks = rng.integers(0, len(sense), n_codons - 2)
    return "ATG" + "".join(sense[i] for i in picks) + table.stop_codons[0]
