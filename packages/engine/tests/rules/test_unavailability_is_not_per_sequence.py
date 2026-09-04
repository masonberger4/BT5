"""No rule's unavailability may depend on the SEQUENCE.

THIS FILE IS A TRIPWIRE. It guards an assumption `design/ranking.py::build_nulls`
makes about the rules lane, and it lives here rather than in the design lane
because this is the side that would break it -- the same reasoning as
`test_fragment.py` and #64.

WHAT build_nulls ASSUMES. It decides whether an objective gets a null -- and
therefore whether ANY candidate in the panel can be scored on it -- from
`evaluations[picks[0]]` alone, the first sweep pick:

    evaluation = anchor.get(spec.id)
    reason = unavailability(evaluation) if ... else "the rule did not run"
    if reason is not None:
        unavailable[spec.id] = reason
        continue

That is sound only while every rule's availability is a property of the HOST and
CONTEXT rather than of the particular synonymous encoding. If a rule were
unavailable on `picks[0]` but computable on another candidate, that candidate
would report the objective `unavailable` even though its own raw score is fine.

`unavailability()` is non-None exactly when `raw_score` is NaN
(`design/ranking.py:129`), so "per-sequence unavailability" means precisely: a
rule that returns NaN for some synonymous encodings of a protein and a real
number for others. That is what this file measures.

WHY IT IS ONLY A TRIPWIRE. The under-reporting is unreachable today -- every
`unavailability()` reason in the catalog is host- or context-level -- and it
never leaks a false number in either direction. The failure mode is an honest
`unavailable` where a real score existed, which is the safe direction for BT5's
reporting posture. #100 asked that it not be lost when a per-sequence reason is
added; this is that.

WHEN THIS FILE FAILS, the named rule has gained a per-sequence unavailability
reason and `build_nulls` is now under-reporting it for every candidate that is
not `picks[0]`. Do NOT exempt the rule here. The fix is in `build_nulls`: decide
availability per candidate, or take the union across picks rather than the
anchor's answer.
"""

from __future__ import annotations

import math

import pytest
from Bio.Data import CodonTable
from bt5.core.registry import all_specs, discover
from bt5.core.services import Services
from bt5.core.spec import Spec
from conftest import construct, context, slot

discover()
SPECS = all_specs()

#: Long enough that windowed and structural rules have something to look at, and
#: rich enough in degenerate families that the two encodings below really differ.
PROTEIN = "MKLVTAAFERSKSVQNYVVSTKDSPLYYLRKWVRSGYKFDCEEVGLREHQGPAATYTPTQAIWRLTLPSPLL"

#: A neutral spacer so the CDS sits in an assembled construct rather than alone.
BACKBONE = "GGATCCAAGCTTGTCGACCTGCAGTTAACCGGTACCGAGCTCGAATTCACGCGTGGTACCTCTAGAGTCGAC"


def _encodings(table_id: int = 1) -> tuple[str, str]:
    """Two maximally different synonymous encodings of the same protein.

    First-listed codon per residue against last-listed. Same protein, same
    length, same frame -- so any difference in a rule's availability between them
    is a difference the SEQUENCE caused, which is the whole question.
    """
    table = CodonTable.unambiguous_dna_by_id[table_id]
    by_aa: dict[str, list[str]] = {}
    for codon, aa in sorted(table.forward_table.items()):
        by_aa.setdefault(aa, []).append(codon)
    first = "".join(by_aa[aa][0] for aa in PROTEIN)
    last = "".join(by_aa[aa][-1] for aa in PROTEIN)
    return first + table.stop_codons[0], last + table.stop_codons[0]


FIRST_CODONS, LAST_CODONS = _encodings()


def test_the_two_encodings_really_are_different_and_synonymous() -> None:
    """A precondition, not a formality. If the two encodings were equal the
    comparison below would hold vacuously for every rule, and this file would
    pass forever while measuring nothing."""
    assert FIRST_CODONS != LAST_CODONS
    assert len(FIRST_CODONS) == len(LAST_CODONS)
    table = CodonTable.unambiguous_dna_by_id[1]

    def translate(cds: str) -> str:
        return "".join(
            table.forward_table[cds[i : i + 3]]
            for i in range(0, len(cds) - 3, 3)
            if cds[i : i + 3] not in table.stop_codons
        )

    assert translate(FIRST_CODONS) == translate(LAST_CODONS) == PROTEIN
    differing = sum(
        FIRST_CODONS[i : i + 3] != LAST_CODONS[i : i + 3] for i in range(0, len(FIRST_CODONS), 3)
    )
    assert differing >= 30, f"only {differing} codons differ; the probe is too weak"


@pytest.mark.parametrize("spec", SPECS, ids=lambda s: s.id)
def test_availability_does_not_depend_on_the_encoding(spec: type[Spec], services: Services) -> None:
    ctx = context(slot())
    rule = spec()

    verdicts: dict[str, bool] = {}
    for label, cds in (("first-codon", FIRST_CODONS), ("last-codon", LAST_CODONS)):
        evaluation = rule.evaluate(construct(cds, BACKBONE), ctx, services)
        verdicts[label] = math.isnan(evaluation.raw_score)

    assert verdicts["first-codon"] == verdicts["last-codon"], (
        f"{spec.id}: availability differs between two synonymous encodings of the "
        f"same protein (NaN raw_score: {verdicts}). `build_nulls` decides this "
        f"objective's availability from picks[0] alone, so every other candidate "
        f"in the panel would now report it `unavailable` even where its own score "
        f"is real. See this file's docstring and #100 -- fix `build_nulls`, do not "
        f"exempt the rule here."
    )
