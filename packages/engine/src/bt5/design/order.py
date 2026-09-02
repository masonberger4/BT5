"""The vendor order file -- the last artefact between a design and real DNA.

`docs/PLAN.md` locks v1's output as "annotated construct + vendor order CSV", and
`bt5.score.order` has shipped `order_entries`, `write_csv` and `write_idt_plate`
since M3 with nothing on the design path calling them. This is that call.

Two choices are made here rather than left to the caller.

**The CDS is what gets ordered, not the construct.** The vendor synthesises the
fragment you clone in; the backbone is already on the bench. Ordering the whole
assembled plasmid would put the user's own vector -- kilobases of it -- on a
synthesis quote. `report.screening_burden` already sizes the colony-picking
burden on `len(candidate.cds)` for the same reason.

**Every candidate in the gallery is an order line.** A gallery whose members
cannot all be ordered is a gallery of one design and four pictures. Each line is
named `entry_name(construct_name, design_hash)`, so the five tubes that come
back are distinguishable -- which is the entire reason the design hash exists,
and why `write_csv` refuses two rows under one name.

The CSV is emitted unconditionally because it is stdlib `csv`. The IDT plate is
a real `.xlsx` and needs openpyxl from the `export` extra, so it stays a
function the caller invokes rather than a field on every result: a design run
must not fail because a spreadsheet library is absent.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import IO

from bt5.core.result import Candidate
from bt5.score.order import (
    DEFAULT_PLATE_NAME,
    OrderEntry,
    order_entries,
    write_csv,
    write_idt_plate,
)


def entries_for(candidates: Sequence[Candidate], *, construct_name: str) -> tuple[OrderEntry, ...]:
    """One order line per candidate, labelled with its design hash.

    The label carries the CANDIDATE's own label as well as the construct name,
    so a gallery reads as `pLV_native_baseline_a1b2c3d4e5f6` rather than five
    rows distinguishable only by a hash -- and `OrderEntry` still refuses two
    rows sharing a name, which is the check that actually holds.
    """
    return order_entries(
        (f"{construct_name}_{candidate.label}", candidate.design_hash, candidate.cds)
        for candidate in candidates
    )


def order_csv(candidates: Sequence[Candidate], *, construct_name: str) -> str:
    """The Name,Sequence CSV shipped alongside the GenBank."""
    return write_csv(entries_for(candidates, construct_name=construct_name))


def write_order_plate(
    candidates: Sequence[Candidate],
    target: str | Path | IO[bytes] | None = None,
    *,
    construct_name: str,
    plate_name: str = DEFAULT_PLATE_NAME,
) -> bytes:
    """The IDT eBlocks 96-well upload workbook, for a caller that wants one.

    Raises `OrderError` when openpyxl is absent, with the extra to install --
    `write_idt_plate`'s own refusal, passed through unchanged rather than
    downgraded to a CSV under an `.xlsx` name.
    """
    return write_idt_plate(
        entries_for(candidates, construct_name=construct_name),
        target,
        plate_name=plate_name,
    )
