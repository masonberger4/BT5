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

    def _shifts(self, construct_length: int, circular: bool) -> tuple[int, ...]:
        """Offsets to try when comparing this interval with another.

        Either one may be the interval that wraps, so a circular comparison has
        to try the other a full turn in both directions.
        """
        if not circular:
            return (0,)
        return (0, construct_length, -construct_length)

    def overlaps(self, other: Interval, construct_length: int, circular: bool) -> bool:
        """Do the two intervals share at least one base?

        `construct_length` and `circular` are REQUIRED, not defaulted, because
        the linear answer is wrong exactly where it matters most. A feature
        stored as [4900, 5100) on a 5000 bp plasmid and a hit at [10, 40) sit on
        the same bases, and plain half-open comparison says they do not.

        There is no normalisation that would let the linear form work. Both of
        those intervals are ALREADY in BT5's one canonical representation, so
        the extra turn has to be tried at comparison time and cannot be
        pre-applied at construction time. An earlier version of this docstring
        told callers to normalise first via a `Construct.normalise` that was
        never written and could not have helped, so three lanes each hand-rolled
        the shift trial below instead -- and the one caller that trusted the
        docstring and called this method directly, the intron-overlap check in
        the vector lane, silently missed every origin-spanning intron.
        """
        return any(
            self.start < other.end + shift and other.start + shift < self.end
            for shift in self._shifts(construct_length, circular)
        )

    def contains(self, inner: Interval, construct_length: int, circular: bool) -> bool:
        """Does `inner` lie WHOLLY inside this interval? Wrap-aware, as above.

        Strand is deliberately ignored by both predicates: they answer "which
        bases", and a caller that also cares which strand compares `.strand`
        itself. Folding strand in here silently would make a reverse-strand hit
        inside a forward-strand feature look like no hit at all.
        """
        return any(
            self.start <= inner.start + shift and inner.end + shift <= self.end
            for shift in self._shifts(construct_length, circular)
        )

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
        """True only if `iv` lies wholly inside a designable CDS segment.

        Wrap-aware, and it has to be: an insert spanning the origin is stored as
        one segment with `end > length`, and the plain comparison this used to
        do called every codon PAST the origin backbone instead. It erred in the
        safe direction -- I9 was never at risk -- but it hands the solver a
        mutation space with a hole in it, which reads downstream as a design
        that cannot satisfy its own constraints rather than as a coordinate bug.
        The new answer is a strict superset of the old one, and identical on a
        linear construct.
        """
        return any(
            s.interval.contains(iv, self.length, self.is_circular)
            for s in self.segments
            if s.is_editable
        )

    def overlaps_editable(self, iv: Interval) -> bool:
        """Does `iv` touch any designable CDS base?

        This, not `is_editable`, is the honest answer to "can codon choice do
        anything about a finding here", and it is what belongs in
        `Breach.fixable_by_codon_choice`. A restriction site or a GC window that
        straddles the CDS/backbone junction is only partly inside the CDS, so
        `is_editable` says no -- yet recoding the codons on the CDS side really
        does destroy the site or move the window, and calling such a breach
        unfixable sends a repairable design to the advisor as a dead end.

        It lives on the contract because every rule needs it and each of the ~45
        of them hand-rolling the same wrap-aware scan is how they end up
        disagreeing about what the junction means.
        """
        return any(e.overlaps(iv, self.length, self.is_circular) for e in self.editable)
