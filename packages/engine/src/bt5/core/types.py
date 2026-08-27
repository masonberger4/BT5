"""Geometry and the construct model.

The single most important decision in BT5 lives here: `Construct.editable` is the
CDS, and its complement is the backbone. Nothing outside `editable` may ever be
mutated, and invariant I9 in `bt5.verify` proves that byte-for-byte.

Coordinates are half-open [start, end), 0-based, in CONSTRUCT space. On a circular
construct an interval with end > length wraps the origin. There is exactly ONE
representation of a wrapping interval; GenBank join() is normalised into it on
import and re-split on export.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Literal

Strand = Literal[1, -1]

DNA_ALPHABET = frozenset("ACGT")

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    """Reverse complement. Defined here so every lane uses one implementation."""
    return seq.translate(_COMPLEMENT)[::-1]


class Topology(StrEnum):
    LINEAR = "linear"
    CIRCULAR = "circular"


class SegmentKind(StrEnum):
    """What a region of the construct is, for the optimizer and the scanners.

    DESIGNABLE_CDS is the only kind the optimizer may rewrite. The two EXEMPT
    kinds exist because some backbone features irreducibly violate rules that are
    otherwise correct:

    - ANNOTATED_INTRON: a deliberately placed chimeric/CMV intron is one of the
      most reliable mammalian expression levers, and a naive "remove all splice
      sites" pass would happily destroy it.
    - WHITELISTED_REPEAT: lentiviral LTRs are long perfect direct repeats and AAV
      ITRs are 145 bp palindromes. Both violate the repeat rules by construction.
      The answer is a strain/temperature protocol recommendation, not a redesign.
    """

    DESIGNABLE_CDS = "designable_cds"
    BACKBONE = "backbone"
    ANNOTATED_INTRON = "annotated_intron"
    WHITELISTED_REPEAT = "whitelisted_repeat"


@dataclass(frozen=True, slots=True, order=True)
class Interval:
    """Half-open [start, end), 0-based, construct coordinates.

    `end > construct.length` means the interval wraps the origin.
    """

    start: int
    end: int
    strand: Strand = 1

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ValueError(f"interval start must be >= 0, got {self.start}")
        if self.end <= self.start:
            raise ValueError(f"interval must be non-empty: [{self.start}, {self.end})")

    @property
    def length(self) -> int:
        return self.end - self.start

    def wraps(self, construct_length: int) -> bool:
        return self.end > construct_length

    def overlaps(self, other: Interval) -> bool:
        """Overlap in linear coordinates. Callers on circular constructs should
        normalise both intervals first (see Construct.normalise)."""
        return self.start < other.end and other.start < self.end

    def extended(self, by: int, construct_length: int, circular: bool) -> Interval:
        """Widen by `by` on both sides. This is how localisation policies turn a
        breach into a repair window.

        On a circular construct an extension that runs off the front wraps into
        the canonical representation (start inside the sequence, end past the
        end) rather than going negative -- a breach near the origin must still
        produce a usable repair window.
        """
        start = self.start - by
        end = self.end + by
        if circular:
            if start < 0:
                start += construct_length
                end += construct_length
        else:
            start = max(0, start)
            end = min(construct_length, end)
        return Interval(start, end, self.strand)


@dataclass(frozen=True, slots=True)
class Feature:
    """A GenBank feature. `qualifiers` is ordered and multi-valued so a
    round-trip preserves the source file byte-for-byte where it matters."""

    interval: Interval
    kind: str
    qualifiers: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    uid: str = ""


@dataclass(frozen=True, slots=True)
class Segment:
    """A typed region of the construct."""

    interval: Interval
    kind: SegmentKind
    label: str = ""

    @property
    def is_editable(self) -> bool:
        return self.kind is SegmentKind.DESIGNABLE_CDS

    @property
    def exempt_from_scanning(self) -> bool:
        """Annotated introns and whitelisted repeats are immutable AND not scanned."""
        return self.kind in (SegmentKind.ANNOTATED_INTRON, SegmentKind.WHITELISTED_REPEAT)


@dataclass(frozen=True, slots=True)
class TranslationUnit:
    """One ORF under one genetic code.

    `table_id` lives here and NOT on Segment. One ribosome reads one ORF under one
    code; a per-segment table is biologically incoherent within a reading frame and
    would pass verification while producing the wrong protein. NCBI table 12
    (Candida) reassigns CTG=Ser rather than Leu, so a wrong table is a silent
    protein-changing bug no assay catches for months.

    There is deliberately no default. A TranslationUnit cannot be constructed
    without an explicit table.
    """

    table_id: int
    codon_map: tuple[Interval, ...]
    protein: str
    has_terminal_stop: bool = True
    #: True for a complete ORF whose first codon must be a valid initiator for
    #: `table_id`. False for a mid-ORF fragment -- BT5 back-translates cassette
    #: elements (tags, linkers, 2A peptides) that legitimately do not begin with
    #: a start codon. I4 is only checked when this is True.
    starts_at_initiator: bool = True

    def __post_init__(self) -> None:
        if self.table_id <= 0:
            raise ValueError(f"table_id must be a positive NCBI table id, got {self.table_id}")


@dataclass(frozen=True, slots=True)
class Provenance:
    """Everything needed to reproduce a design. Persisted into the report, the
    GenBank note, and the order file."""

    app_version: str
    seed: int
    engine_versions: Mapping[str, str] = field(default_factory=dict)
    codon_table_name: str = ""
    constraint_set_hash: str = ""
    degradations: tuple[str, ...] = ()


@dataclass(frozen=True)
class Construct:
    """The ONLY thing a rule is ever evaluated against.

    There is no API anywhere in BT5 that evaluates a rule against a bare `str`.
    That single decision is what makes "evaluate on the assembled circular plasmid"
    the default rather than a feature someone has to remember to turn on --
    junction-spanning and origin-spanning hits are found because they cannot be
    missed, not because a rule author thought about them.
    """

    sequence: str
    topology: Topology
    segments: tuple[Segment, ...]
    translation_units: tuple[TranslationUnit, ...] = ()
    features: tuple[Feature, ...] = ()
    annotations: Mapping[str, str] = field(default_factory=dict)
    provenance: Provenance | None = None

    def __post_init__(self) -> None:
        bad = set(self.sequence.upper()) - DNA_ALPHABET
        if bad:
            raise ValueError(f"construct sequence contains non-ACGT characters: {sorted(bad)}")

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def is_circular(self) -> bool:
        return self.topology is Topology.CIRCULAR

    @property
    def editable(self) -> tuple[Interval, ...]:
        """The CDS. Its complement is the backbone, immutable by construction."""
        return tuple(s.interval for s in self.segments if s.is_editable)

    @property
    def exempt(self) -> tuple[Interval, ...]:
        """Regions immutable AND exempt from motif/repeat scanning."""
        return tuple(s.interval for s in self.segments if s.exempt_from_scanning)

    def slice(self, iv: Interval) -> str:
        """Wrap- and strand-aware extraction."""
        n = self.length
        if iv.end <= n:
            sub = self.sequence[iv.start : iv.end]
        elif self.is_circular:
            sub = self.sequence[iv.start :] + self.sequence[: iv.end - n]
        else:
            raise ValueError(f"interval {iv} runs past the end of a linear construct of length {n}")
        return reverse_complement(sub) if iv.strand == -1 else sub

    def tripled(self) -> tuple[str, int]:
        """Circular constructs linearised as seq*3 with offset L.

        Constraints are resolved only where they overlap [L, 2L), then
        `sequence[L:2L]` is emitted. This is how origin-spanning windows and
        motifs are handled without special-casing every rule.
        """
        if not self.is_circular:
            return self.sequence, 0
        return self.sequence * 3, self.length

    def is_editable(self, iv: Interval) -> bool:
        """True only if `iv` lies wholly inside a designable CDS segment."""
        return any(
            s.interval.start <= iv.start and iv.end <= s.interval.end
            for s in self.segments
            if s.is_editable
        )
