"""The properties ARE the specification.

Note what is deliberately absent: any assertion that CAI must exceed a threshold,
or that an optimized sequence must beat a random back-translation on CAI. The
evidence establishes CAI has no predictive value for expression and must be
treated as a BAND, so a blocking CI gate rewarding max-CAI would train every
contributor to optimize the one metric the science refutes.
"""

from __future__ import annotations

import pytest
from Bio.Data import CodonTable
from bt5.core.result import VerificationError
from bt5.verify import expand_iupac, verify_solution
from hypothesis import given
from hypothesis import strategies as st
from strategies import STANDARD_AA, proteins, repetitive_proteins

BACTERIAL, STANDARD = 11, 1


def naive_back_translate(protein: str, table_id: int) -> str:
    """Deliberately dumb: first synonymous codon, every time.

    This is the ReferenceSolver in miniature -- obviously correct, never clever.
    It exists so the invariants have something to run against from day one, and
    it survives as a differential oracle.

    ONE rule beyond "pick the first codon", and it is a real product requirement:
    never emit a codon that is also a stop codon in the target table. NCBI tables
    27 and 28 (ciliate / karyorelict) make TGA both Trp AND a stop, so a naive
    lookup would encode Trp as TGA and terminate translation mid-protein. Every
    back-translator in BT5 must exclude stop codons from the synonymous set.
    """
    table = CodonTable.unambiguous_dna_by_id[table_id]
    stops = set(table.stop_codons)
    first: dict[str, str] = {}
    for codon, aa in sorted(table.forward_table.items()):
        if codon in stops:
            continue
        first.setdefault(aa, codon)
    missing = sorted(set(protein) - set(first))
    if missing:
        raise ValueError(f"no non-stop codon encodes {missing} in NCBI table {table_id}")
    return "".join(first[aa] for aa in protein)


ALL_TABLE_IDS = sorted(CodonTable.unambiguous_dna_by_id)


@pytest.mark.parametrize("table_id", ALL_TABLE_IDS)
def test_round_trip_holds_for_every_ncbi_table(table_id: int) -> None:
    """translate(back_translate(p, t), t) == p, for EVERY shipped table.

    NCBI table 12 reassigns CTG to Ser rather than Leu and table 4 makes TGA
    Trp rather than a stop. A tool that assumes table 1 emits a silently wrong
    protein that no assay catches for months.
    """
    table = CodonTable.unambiguous_dna_by_id[table_id]
    protein = "".join(sorted({aa for aa in table.forward_table.values() if aa in STANDARD_AA}))
    dna = naive_back_translate(protein, table_id)
    verify_solution(protein, dna, table_id=table_id)


@given(protein=proteins)
def test_naive_back_translation_round_trips(protein: str) -> None:
    verify_solution(protein, naive_back_translate(protein, BACTERIAL), table_id=BACTERIAL)


@given(protein=repetitive_proteins)
def test_repetitive_proteins_round_trip(protein: str) -> None:
    """Antibodies, (GGGGS)n linkers, His tags and tandem 2A peptides are the
    proteins people actually express, and the ones naive back-translation turns
    into perfect nucleotide repeats."""
    verify_solution(protein, naive_back_translate(protein, STANDARD), table_id=STANDARD)


@given(protein=proteins)
def test_output_length_is_three_times_protein(protein: str) -> None:
    assert len(naive_back_translate(protein, BACTERIAL)) == 3 * len(protein)


@given(protein=proteins)
def test_back_translation_is_deterministic(protein: str) -> None:
    a = naive_back_translate(protein, BACTERIAL)
    b = naive_back_translate(protein, BACTERIAL)
    assert a == b, "two runs of one protein must not produce two different tubes"


@given(protein=proteins, extra=st.integers(min_value=1, max_value=2))
def test_frame_violation_is_always_refused(protein: str, extra: int) -> None:
    """A length that is not a multiple of three must raise, never truncate.

    Biopython's Seq.translate() emits only a warning and silently truncates here,
    so a frame bug would pass a naive round-trip test.
    """
    dna = naive_back_translate(protein, BACTERIAL) + "A" * extra
    with pytest.raises(VerificationError) as exc:
        verify_solution(protein, dna, table_id=BACTERIAL)
    assert exc.value.invariant == "I2"


@given(protein=proteins)
def test_wrong_protein_is_always_caught(protein: str) -> None:
    """Swap the first codon for one encoding a DIFFERENT amino acid.

    Picking a codon at random is not enough -- TTC and TTT both encode Phe, so a
    naive "replace with TTT" is a synonymous change the oracle should NOT flag.
    """
    table = CodonTable.unambiguous_dna_by_id[BACTERIAL]
    stops = set(table.stop_codons)
    dna = naive_back_translate(protein, BACTERIAL)
    target = protein[0]
    replacement = next(
        c for c, aa in sorted(table.forward_table.items()) if aa != target and c not in stops
    )
    mutated = replacement + dna[3:]

    with pytest.raises(VerificationError) as exc:
        verify_solution(protein, mutated, table_id=BACTERIAL)
    assert exc.value.invariant == "I3"


def test_iupac_expansion() -> None:
    assert set(expand_iupac("RY")) == {"AC", "AT", "GC", "GT"}
    assert expand_iupac("ACGT") == ["ACGT"]
    with pytest.raises(ValueError, match="not an IUPAC code"):
        expand_iupac("AZ")
