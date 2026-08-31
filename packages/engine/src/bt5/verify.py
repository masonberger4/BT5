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
  I7  GC band on each ORDERED (designable) span -- the DNA a vendor synthesises
  I8  homopolymer / repeat ceilings across the whole construct
  I9  every backbone base is byte-identical to the input backbone   <-- highest value
  I10 cassette frame invariant across the assembled CDS
"""

from __future__ import annotations

from collections.abc import Sequence
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
    max_homopolymer: int | None = None,
    max_repeat: int | None = None,
    original_backbone: Construct | None = None,
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

    # The CALLER's own declaration, checked against the construct's.
    #
    # These two arguments were required and then read by nothing: every check
    # below runs off `tu.protein` and `tu.table_id`, so the oracle validated the
    # assembler's claim about the construct rather than the caller's claim about
    # what it should be. Two live holes closed here. A construct with no
    # translation unit passed I3/I4/I5 VACUOUSLY -- armed and silent, the same
    # shape of bug I7 above had. And a caller asking for table 12 against a unit
    # declaring table 1 was verified under table 1 while believing CTG=Ser had
    # been checked, which is the highest-severity silent bug this project names.
    units = construct.translation_units or ()
    if protein and not units:
        raise VerificationError(
            "I3",
            f"a {len(protein)}-residue protein was declared but the construct carries "
            "no translation unit, so no round trip was verified",
        )
    if units and table_id not in {tu.table_id for tu in units}:
        declared = sorted({tu.table_id for tu in units})
        raise VerificationError(
            "I3",
            f"caller declared NCBI table {table_id} but no translation unit uses it "
            f"(units declare {declared}); table 12 reassigns CTG to Ser, so a table "
            "disagreement is a silently wrong protein",
        )

    # I2/I3/I4/I5 -- per translation unit
    for tu in units:
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

    # I7 -- GC band on each ORDERED (designable) span.
    #
    # SCOPE IS WHAT THIS INVARIANT USED TO GET WRONG. A GC band is a
    # MANUFACTURABILITY bound: it describes what a synthesis vendor's order-entry
    # checker accepts in one tube, and the tube holds the designable span. Nobody
    # synthesises the user's backbone. Measured across the whole construct, a
    # near-neutral backbone drags a GC-extreme insert into band, and the one
    # invariant whose job is refusing unbuildable DNA never fires on a vector
    # design: a 900 bp / 90% insert reads 0.629 against a 2 kb / 50% backbone and
    # passed, though both vendors deny the fragment.
    #
    # The 50 bp windowed band that used to live here is GONE, and its removal is a
    # MEASUREMENT rather than a preference. An 18-sequence ladder through two
    # vendors' order-entry checkers (docs/design/vendor-gc-calibration.md) put
    # sixteen probes carrying one extreme 50 bp window -- 4% to 96% local GC on an
    # otherwise neutral background -- through both, and all sixteen were accepted.
    # One probe accepted by both carries a 100 bp window at 23% GC, LOWER than a
    # probe refused at 21% GLOBAL. Only the overall GC of the ordered DNA
    # separates accepted from refused, so a windowed band refused DNA both vendors
    # manufacture, without comment. A window is never a breach.
    #
    # Two designable spans are two synthesis reactions, so each is held to the
    # band on its own and the construct passes only if every one of them does.
    # `slice` is wrap- and strand-aware, and both matter: a designable span may
    # itself cross the origin (stored as end > length), and a reverse-oriented
    # cassette's span comes back reverse-complemented -- harmless for GC, which is
    # invariant under it, but the bases must be the ORDERED ones either way.
    #
    # ASSUMPTION, written down so it is falsifiable: one designable span is one
    # ordered fragment. If the vendor lane ever splits a long span into several
    # tubes by length limit, a span-level average could hide an out-of-band tube
    # and this check would become weaker than the rule it backstops.
    #
    # Adapter sequence is NOT included: it is vendor data, and this module stays
    # independent of the lanes it validates. The extension point is a parameter
    # (`adapters: tuple[str, str]`), never an import.
    if gc_bounds is not None:
        lo, hi = gc_bounds
        spans = sorted(construct.editable)
        if not spans:
            raise VerificationError(
                "I7",
                "a GC band was requested but this construct declares no designable "
                "segment, so there is no ordered DNA to hold to it",
            )
        for iv in spans:
            frac = gc_fraction(construct.slice(iv))
            if not (lo <= frac <= hi):
                raise VerificationError(
                    "I7",
                    f"ordered GC {frac:.3f} over the {iv.length}bp designable span "
                    f"at {iv.start} outside [{lo}, {hi}]",
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
    # Required, never defaulted. CLAUDE.md 3.1: NCBI table 12 reassigns CTG to
    # Ser and table 4 makes TGA Trp, so a defaulted table is a silently wrong
    # protein that no assay catches for months. `TranslationUnit` has never had a
    # default; this wrapper quietly did.
    table_id: int,
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
