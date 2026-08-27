"""Things BT5 wants the user to know about a design.

A note is not a log line. It is a typed statement with a location, what it bears
on, and what the user can do about it -- because the alternative, a tuple of free
text, forces every consumer (the GenBank COMMENT, the QC report, the UI) to
re-derive structure by matching on strings BT5 itself wrote.

The four kinds are deliberately few:

  unavailable   BT5 could not evaluate something. The 5' folding objective with
                no annotated 5'UTR is the canonical case. Reporting this is what
                stops the app from folding the CDS alone and presenting the
                number as if the UTR were there.
  assumption    BT5 proceeded on something the input did not state outright,
                such as a transcription start inferred from a promoter.
  liability     A real property of the sequence that carries risk. These are
                usually in the BACKBONE, which BT5 never edits, so the honest
                output is a located warning and a protocol suggestion.
  change        BT5 altered the map itself -- rotated the origin, dropped a
                feature whose coordinates no longer describe anything.

None of these ever carries a predicted expression level, titer or yield. BT5
reports what is there and what it could not do, never what will happen.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from bt5.core.types import Interval

NoteKind = Literal["unavailable", "assumption", "liability", "change"]


def format_span(interval: Interval, length: int) -> str:
    """1-based, GenBank-style, splitting an origin-spanning interval like join()."""
    if interval.end <= length:
        return f"{interval.start + 1}..{interval.end}"
    return f"{interval.start + 1}..{length},1..{interval.end - length}"


@dataclass(frozen=True, slots=True)
class DesignNote:
    """One statement about a design, addressed to the person holding the tube."""

    kind: NoteKind
    summary: str
    interval: Interval | None = None
    bears_on: str = ""
    action: str = ""

    def render(self, length: int | None = None) -> str:
        """One line of plain text, for the GenBank COMMENT and the QC report."""
        parts = [f"[{self.kind}]"]
        if self.interval is not None and length is not None:
            parts.append(f"{format_span(self.interval, length)}:")
        parts.append(self.summary.rstrip("."))
        if self.bears_on:
            parts.append(f"- bears on {self.bears_on}")
        if self.action:
            parts.append(f"- action: {self.action}")
        return " ".join(parts) + "."
