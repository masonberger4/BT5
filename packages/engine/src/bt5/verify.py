"""The independent correctness oracle.

THIS FILE IS LOAD-BEARING. Read before editing; changes require the
`approved:oracle-change` label.

Why it lives in `src/` and not `tests/`: if verification lives in the test suite,
the agent that writes a broken optimizer also writes the test, and the bug ships.
Here it runs on EVERY optimize() call, so a wrong sequence raises in the user's
app, in every property example, in every golden test and in every benchmark row --
from ONE definition.

Why it may not import any lane module: three callers of one pure function is zero
independence. This file re-derives every invariant from Biopython and the standard
library only. A CI check enforces the import ban.

Invariants:
  I1  alphabet is ACGT only
  I2  frame: CDS length is a multiple of 3
  I3  round trip: translation equals the declared protein
  I4  initiator codon is valid for the declared table
  I5  stops: no interior in-frame stop; terminal stop matches declaration
  I6  forbidden motifs absent on the CIRCULAR construct, both strands, including
      junction- and origin-spanning hits, and including inside the backbone
  I7  GC band, global and windowed, with windows that wrap the origin
  I8  homopolymer / repeat ceilings across the whole construct
  I9  every backbone base is byte-identical to the input backbone   <-- highest value
  I10 cassette frame invariant across the assembled CDS
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import cast

from Bio.Data import CodonTable

from bt5.core.result import VerificationError
from bt5.core.types import (
    DNA_ALPHABET,
    Construct,
    Interval,
    SegmentKind,
    reverse_complement,
)

IUPAC_EXPANSION = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}


def expand_iupac(pattern: str) -> list[str]:
    """Expand an IUPAC pattern to concrete ACGT strings."""
    out = [""]
    for ch in pattern.upper():
        opts = IUPAC_EXPANSION.get(ch)
        if opts is None:
            raise ValueError(f"not an IUPAC code: {ch!r} in {pattern!r}")
        out = [prefix + o for prefix in out for o in opts]
    return out


def _codon_table(table_id: int) -> CodonTable.CodonTable:
    try:
        # Biopython ships no type stubs, so this lookup is Any without the cast.
        return cast("CodonTable.CodonTable", CodonTable.unambiguous_dna_by_id[table_id])
    except KeyError as exc:
        raise VerificationError("I3", f"unknown NCBI translation table {table_id}") from exc


def _translate(dna: str, table_id: int) -> str:
    """Translate WITHOUT Biopython's Seq.translate().

    Seq.translate() emits only a BiopythonWarning and SILENTLY TRUNCATES when the
    length is not a multiple of three, so a frame-length bug would sail through a
    naive round-trip test. We check length first and index the table directly.
    """
    if len(dna) % 3 != 0:
        raise VerificationError("I2", f"CDS length {len(dna)} is not a multiple of 3")
    table = _codon_table(table_id)
    aas: list[str] = []
    for i in range(0, len(dna), 3):
        codon = dna[i : i + 3]
        if codon in table.stop_codons:
            aas.append("*")
        else:
            try:
                aas.append(table.forward_table[codon])
            except KeyError as exc:
                raise VerificationError("I1", f"untranslatable codon {codon!r} at {i}") from exc
    return "".join(aas)


def _windows(seq: str, size: int, step: int, circular: bool) -> Iterable[tuple[int, str]]:
    n = len(seq)
    if n == 0:
        return
    if circular:
        doubled = seq + seq[: size - 1] if size > 1 else seq
        for start in range(0, n, step):
            yield start, doubled[start : start + size]
    else:
        for start in range(0, max(1, n - size + 1), step):
            yield start, seq[start : start + size]


def gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def longest_homopolymer(seq: str) -> tuple[str, int]:
    best_base, best = "", 0
    run_base, run = "", 0
    for ch in seq:
        if ch == run_base:
            run += 1
        else:
            run_base, run = ch, 1
        if run > best:
            best_base, best = run_base, run
    return best_base, best


def find_motifs(construct: Construct, patterns: Sequence[str]) -> list[tuple[str, int]]:
    """Locate forbidden motifs on BOTH strands of the assembled construct.

    Two bug classes dominate here and both are invisible to any per-codon test:
    motifs created ACROSS CODON BOUNDARIES, and motifs on the REVERSE STRAND. A
    circular construct adds a third: hits that span the origin. We search the
    doubled sequence and close the pattern set under reverse complement, so all
    three are structurally impossible to miss.
    """
    seq = construct.sequence
    n = len(seq)
    longest = max((len(p) for p in patterns), default=1)
    haystack = seq + seq[: longest - 1] if construct.is_circular else seq

    exempt: list[Interval] = list(construct.exempt)
    hits: list[tuple[str, int]] = []
    for pattern in patterns:
        for concrete in expand_iupac(pattern):
            for probe in {concrete, reverse_complement(concrete)}:
                start = haystack.find(probe)
                while start != -1:
                    pos = start % n if n else start
                    if not any(iv.start <= pos < iv.end for iv in exempt):
                        hits.append((pattern, pos))
                    start = haystack.find(probe, start + 1)
    return sorted(set(hits), key=lambda h: (h[1], h[0]))


def verify_construct(
    construct: Construct,
    *,
    protein: str,
    table_id: int,
    forbidden: Sequence[str] = (),
    gc_bounds: tuple[float, float] | None = None,
    gc_window: int = 50,
    max_homopolymer: int | None = None,
    max_repeat: int | None = None,
    original_backbone: Construct | None = None,
    expect_terminal_stop: bool = True,
) -> None:
    """Re-derive every invariant. Raise VerificationError on the first failure.

    This function REFUSES TO PASS a construct it cannot prove correct. Callers
    must treat an exception as "do not emit", never as advisory.
    """
    seq = construct.sequence

    # I1 -- alphabet
    bad = sorted(set(seq.upper()) - DNA_ALPHABET)
    if bad:
        raise VerificationError("I1", f"non-ACGT characters present: {bad}")

    # I2/I3/I4/I5 -- per translation unit
    for tu in construct.translation_units or ():
        cds = "".join(construct.slice(iv) for iv in tu.codon_map)
        if len(cds) % 3 != 0:  # I2
            raise VerificationError("I2", f"CDS length {len(cds)} is not a multiple of 3")

        observed = _translate(cds, tu.table_id)
        expected = tu.protein
        body = observed[:-1] if observed.endswith("*") else observed

        if observed.endswith("*") != tu.has_terminal_stop:  # I5
            raise VerificationError(
                "I5",
                f"terminal stop mismatch: sequence {'has' if observed.endswith('*') else 'lacks'} "
                f"one, declaration says {tu.has_terminal_stop}",
            )
        if "*" in body:  # I5
            raise VerificationError("I5", f"interior in-frame stop at codon {body.index('*')}")
        if body != expected:  # I3 -- the round trip
            for i, (a, b) in enumerate(zip(body, expected, strict=False)):
                if a != b:
                    raise VerificationError(
                        "I3",
                        f"translation differs from input protein at residue {i}: {a!r} != {b!r}",
                    )
            raise VerificationError(
                "I3", f"translation length {len(body)} != protein length {len(expected)}"
            )
        if cds and tu.starts_at_initiator:  # I4
            table = _codon_table(tu.table_id)
            if cds[:3] not in table.start_codons:
                raise VerificationError(
                    "I4",
                    f"{cds[:3]!r} is not a valid initiator for NCBI table {tu.table_id} "
                    f"(valid: {sorted(table.start_codons)})",
                )

    # I6 -- forbidden motifs, both strands, circular-aware, backbone included
    if forbidden:
        hits = find_motifs(construct, list(forbidden))
        if hits:
            pattern, pos = hits[0]
            raise VerificationError(
                "I6", f"forbidden motif {pattern!r} present at position {pos} ({len(hits)} total)"
            )

    # I7 -- GC band, global and windowed, windows wrap the origin
    if gc_bounds is not None:
        lo, hi = gc_bounds
        overall = gc_fraction(seq)
        if not (lo <= overall <= hi):
            raise VerificationError("I7", f"global GC {overall:.3f} outside [{lo}, {hi}]")
        step = max(1, gc_window // 5)
        for start, window in _windows(seq, gc_window, step, construct.is_circular):
            if len(window) < gc_window:
                continue
            frac = gc_fraction(window)
            if not (lo <= frac <= hi):
                raise VerificationError(
                    "I7", f"GC {frac:.3f} in {gc_window}bp window at {start} outside [{lo}, {hi}]"
                )

    # I8 -- homopolymer and repeat ceilings
    if max_homopolymer is not None:
        scan = seq + seq[:max_homopolymer] if construct.is_circular else seq
        base, run = longest_homopolymer(scan)
        if run > max_homopolymer:
            raise VerificationError(
                "I8", f"homopolymer run of {run} {base}s exceeds {max_homopolymer}"
            )
    if max_repeat is not None:
        seen: dict[str, int] = {}
        scan = seq + seq[:max_repeat] if construct.is_circular else seq
        for i in range(len(scan) - max_repeat + 1):
            kmer = scan[i : i + max_repeat]
            if kmer in seen:
                raise VerificationError(
                    "I8", f"repeat of length {max_repeat} at {seen[kmer]} and {i}: {kmer!r}"
                )
            seen[kmer] = i

    # I9 -- THE BACKBONE IS UNTOUCHED
    # The worst possible bug in a vector-context tool is silently editing the
    # user's own vector. This makes that a raised exception rather than a
    # shipped plasmid.
    if original_backbone is not None:
        for seg in construct.segments:
            if seg.kind is SegmentKind.DESIGNABLE_CDS:
                continue
            here = construct.slice(seg.interval)
            there = original_backbone.slice(seg.interval)
            if here != there:
                for i, (a, b) in enumerate(zip(here, there, strict=False)):
                    if a != b:
                        raise VerificationError(
                            "I9",
                            f"backbone segment {seg.label or seg.kind} was modified at offset {i} "
                            f"(construct position {seg.interval.start + i}): {b!r} -> {a!r}",
                        )
                raise VerificationError(
                    "I9", f"backbone segment {seg.label or seg.kind} changed length"
                )

    # I10 -- cassette frame invariant
    for tu in construct.translation_units or ():
        total = sum(iv.length for iv in tu.codon_map)
        if total % 3 != 0:
            raise VerificationError(
                "I10", f"assembled cassette spans {total} bases, not a multiple of 3"
            )


def verify_solution(
    protein: str,
    dna: str,
    *,
    table_id: int = 1,
    forbidden: Sequence[str] = (),
    require_initiator: bool = False,
) -> None:
    """Thin linear wrapper, kept so simple call sites and property tests stay short.

    `require_initiator` defaults to False because this wrapper is used for
    arbitrary peptide fragments (tags, linkers, 2A peptides) as well as complete
    ORFs. Pass True when verifying a full CDS.
    """
    from bt5.core.types import Segment, Topology, TranslationUnit

    # I2 must be checked BEFORE building the codon map: a length that is not a
    # multiple of 3 would otherwise produce an interval running past the end of
    # the sequence and surface as a ValueError from Construct.slice rather than
    # as the VerificationError the caller is required to handle.
    if len(dna) % 3 != 0:
        raise VerificationError("I2", f"CDS length {len(dna)} is not a multiple of 3")

    c = Construct(
        sequence=dna,
        topology=Topology.LINEAR,
        segments=(Segment(Interval(0, len(dna)), SegmentKind.DESIGNABLE_CDS, "cds"),),
        translation_units=(
            TranslationUnit(
                table_id=table_id,
                codon_map=tuple(Interval(i, i + 3) for i in range(0, len(dna), 3)),
                protein=protein,
                has_terminal_stop=len(dna) >= 3 and dna[-3:] in _codon_table(table_id).stop_codons,
                starts_at_initiator=require_initiator,
            ),
        ),
    )
    verify_construct(c, protein=protein, table_id=table_id, forbidden=forbidden)
