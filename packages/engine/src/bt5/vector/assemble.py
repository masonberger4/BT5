"""Splice a designed CDS into a vector and produce the assembled `Construct`.

This is where the CDS/backbone boundary stops being a convention and becomes
structure. Every base the optimizer is allowed to touch ends up inside a
`DESIGNABLE_CDS` segment; everything else is `BACKBONE` or one of the two exempt
kinds. Invariant I9 then proves byte-identity of the complement, so "the optimizer
silently edited your vector" is a raised exception rather than a shipped plasmid.

`Assembly.reference` exists for exactly that check. It is built from the PARSED
INPUT vector, not from the assembled output, so it is an independent statement of
what the backbone was supposed to be.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace

from Bio.Data import CodonTable

from bt5.core.types import (
    DNA_ALPHABET,
    Construct,
    Feature,
    Interval,
    Segment,
    SegmentKind,
    Strand,
    TranslationUnit,
    reverse_complement,
)
from bt5.vector.backbone import (
    DEFAULT_EXEMPT_KINDS,
    DEFAULT_EXEMPT_LABELS,
    InsertionSite,
    UtrContext,
    VectorBackbone,
    VectorError,
    rotate_interval,
)
from bt5.vector.markers import is_recombination_site
from bt5.vector.notes import DesignNote
from bt5.vector.remap import IntervalRemapper

#: Filler for the CDS span of the I9 reference construct. Never compared -- I9
#: only slices non-CDS segments -- but `Construct` requires a valid alphabet, so
#: it cannot be N.
_REFERENCE_FILLER = "A"

#: INSDC feature keys whose coordinates describe RESIDUES, not bases. When the
#: redesigned CDS encodes the same protein these still describe exactly what they
#: did before, so they survive the redesign; a `primer_bind` or a restriction
#: site over the same span does not, because its bases changed underneath it.
#:
#: This is the difference between keeping and losing the annotation a user cares
#: about most. Real vector maps annotate cassette elements -- a signal peptide, a
#: Myc tag, a TM domain -- as sub-features INSIDE the ORF, and those elements are
#: exactly what BT5 back-translates as part of the whole CDS.
PROTEIN_LEVEL_KINDS = frozenset(
    {"cds", "sig_peptide", "mat_peptide", "transit_peptide", "propeptide"}
)


@dataclass(frozen=True)
class Assembly:
    """A designed CDS spliced into a vector, plus everything needed to check it."""

    construct: Construct
    #: Pass as `original_backbone=` to `verify_construct` to arm invariant I9.
    reference: Construct
    site: InsertionSite
    remapper: IntervalRemapper
    utr: UtrContext
    rotation: int
    compound_parts: Mapping[str, tuple[Interval, ...]]
    notes: tuple[DesignNote, ...] = ()

    @property
    def cds_interval(self) -> Interval:
        return self.remapper.insert_interval


def assemble(
    backbone: VectorBackbone,
    cds: str,
    *,
    protein: str,
    table_id: int,
    site: InsertionSite | None = None,
    label: str | None = None,
    has_terminal_stop: bool = True,
    exempt_kinds: Mapping[str, SegmentKind] = DEFAULT_EXEMPT_KINDS,
    exempt_labels: Sequence[str] = DEFAULT_EXEMPT_LABELS,
) -> Assembly:
    """Splice `cds` into `backbone` at `site`, returning the assembled construct.

    `table_id` is required and cross-validated against the vector's own
    `/transl_table`. A mismatch is an error, not a preference: NCBI table 12
    reassigns CTG to Ser and table 4 makes TGA a Trp, so quietly picking one
    produces a plasmid whose protein is wrong in a way no expression assay
    catches for months.
    """
    site = site if site is not None else backbone.find_insertion_site(label=label)
    _check_table(site, table_id)
    _check_cds(cds)
    _check_no_intron_in_cds(backbone, site, exempt_kinds, exempt_labels)

    rotation = 0
    if site.interval.end > backbone.length:
        if not backbone.is_circular:
            raise VectorError(f"insertion site {site.interval} runs past a linear vector")
        rotation = site.interval.start
        backbone = backbone.rotated(rotation)
        site = replace(
            site, interval=rotate_interval(site.interval, by=rotation, length=backbone.length)
        )

    utr = backbone.utr_context(site)
    remapper = IntervalRemapper(
        replaced=site.interval,
        new_insert_length=len(cds),
        old_length=backbone.length,
        circular=backbone.is_circular,
    )

    old = backbone.sequence
    # `cds` is always the CODING sequence, 5'->3' along the mRNA. For a reverse
    # oriented cassette the plus strand of the plasmid must carry its reverse
    # complement, and the codon map must run from high coordinates to low.
    # Inserting it verbatim would translate to a different protein entirely, and
    # I3 would be the only thing standing between that and a shipped plasmid.
    inserted = cds if site.strand == 1 else reverse_complement(cds)
    sequence = old[: site.interval.start] + inserted + old[site.interval.end :]
    if len(sequence) != remapper.new_length:  # pragma: no cover - arithmetic guard
        raise VectorError("assembled length disagrees with the remapper; refusing to emit")

    features, compound_parts, dropped = _remap_features(
        backbone,
        remapper,
        site,
        protein_preserved=_encodes_same_protein(backbone, site, table_id, protein),
    )
    features = (*features, _cds_feature(remapper.insert_interval, site, table_id))

    segments = _segments(
        backbone=backbone,
        remapper=remapper,
        features=features,
        exempt_kinds=exempt_kinds,
        exempt_labels=exempt_labels,
        cds_strand=site.strand,
    )

    codon_map = _codon_map(remapper.insert_interval)
    unit = TranslationUnit(
        table_id=table_id,
        codon_map=codon_map,
        protein=protein,
        has_terminal_stop=has_terminal_stop,
        starts_at_initiator=True,
    )

    notes = (*backbone.notes, *utr.notes, *dropped)
    construct = Construct(
        sequence=sequence,
        topology=backbone.topology,
        segments=segments,
        translation_units=(unit,),
        features=features,
        annotations=dict(backbone.annotations),
    )
    reference = Construct(
        # Built from the PARSED INPUT, so a splice that loses or duplicates a
        # backbone base shows up as an I9 failure rather than agreeing with
        # itself.
        sequence=old[: site.interval.start]
        + _REFERENCE_FILLER * len(inserted)
        + old[site.interval.end :],
        topology=backbone.topology,
        segments=segments,
    )
    return Assembly(
        construct=construct,
        reference=reference,
        site=site,
        remapper=remapper,
        utr=_remap_utr(utr, remapper),
        rotation=rotation,
        compound_parts=compound_parts,
        notes=notes,
    )


def _codon_map(insert: Interval) -> tuple[Interval, ...]:
    """One interval per codon, in TRANSLATION order.

    Forward: ascending construct coordinates. Reverse: descending, each interval
    carrying strand -1 so `Construct.slice` reverse-complements it. Ascending
    order on a reverse cassette yields the reverse-complement protein, which
    still round-trips through a naive per-codon check.
    """
    if insert.strand == 1:
        return tuple(Interval(i, i + 3, 1) for i in range(insert.start, insert.end, 3))
    return tuple(Interval(i - 3, i, -1) for i in range(insert.end, insert.start, -3))


# -- validation -----------------------------------------------------------


def _check_table(site: InsertionSite, table_id: int) -> None:
    detected = site.detected_table_id
    if detected is not None and detected != table_id:
        raise VectorError(
            f"the vector annotates /transl_table={detected} at {site.label!r} but the "
            f"design asked for table {table_id}; resolve this deliberately rather than "
            f"letting one of them win"
        )


def _encodes_same_protein(
    backbone: VectorBackbone, site: InsertionSite, table_id: int, protein: str
) -> bool:
    """Does the sequence being replaced already encode the protein being designed?

    True for the dominant case -- re-optimising a CDS in place -- and false when
    the user is swapping in a different protein, where an internal annotation's
    residue coordinates would describe nothing. Answered by translating, not
    assumed from lengths: two proteins of equal length are equally likely here.
    """
    old = backbone.slice(site.interval)
    if len(old) % 3 != 0:
        return False
    try:
        table = CodonTable.unambiguous_dna_by_id[table_id]
    except KeyError:
        return False
    residues: list[str] = []
    for i in range(0, len(old), 3):
        codon = old[i : i + 3]
        if codon in table.stop_codons:
            residues.append("*")
            continue
        aa = table.forward_table.get(codon)
        if aa is None:
            return False
        residues.append(aa)
    observed = "".join(residues)
    return observed.rstrip("*") == protein.rstrip("*") and "*" not in observed.rstrip("*")


def _is_in_frame_inside(feature: Feature, site: InsertionSite) -> bool:
    """A protein-level feature lying wholly inside the CDS, on a codon boundary."""
    if feature.kind.lower() not in PROTEIN_LEVEL_KINDS:
        return False
    iv = feature.interval
    cds = site.interval
    if not (cds.start <= iv.start and iv.end <= cds.end):
        return False
    offset = iv.start - cds.start if site.strand == 1 else cds.end - iv.end
    return offset % 3 == 0 and iv.length % 3 == 0


def _check_no_intron_in_cds(
    backbone: VectorBackbone,
    site: InsertionSite,
    exempt_kinds: Mapping[str, SegmentKind],
    exempt_labels: Sequence[str],
) -> None:
    """Refuse a vector that annotates an intron inside the CDS being replaced.

    There are no introns in a BT5 CDS. An intron placed to aid expression lives in
    the 5'UTR, which is backbone. Reaching this means the vector models something
    BT5 does not support, and dropping the feature silently -- which is what the
    remapper would otherwise do, since it overlaps the replaced span -- would
    delete an annotation the user put there on purpose.
    """
    covered = {p % backbone.length for p in range(site.interval.start, site.interval.end)}
    for feature in backbone.features:
        label = backbone.label_of(feature)
        if (
            _is_exempt(feature, label, exempt_kinds, exempt_labels)
            is not SegmentKind.ANNOTATED_INTRON
        ):
            continue
        span = range(feature.interval.start, feature.interval.end)
        if any(p % backbone.length in covered for p in span):
            raise VectorError(
                f"annotated intron {label!r} overlaps the designable CDS. There are no "
                f"introns inside a BT5 CDS: an intron placed to aid expression belongs in "
                f"the 5'UTR, which is backbone. Re-annotate it or move the insertion site."
            )


def _check_cds(cds: str) -> None:
    bad = sorted(set(cds.upper()) - DNA_ALPHABET)
    if bad:
        raise VectorError(f"designed CDS contains non-ACGT characters: {bad}")
    if len(cds) % 3 != 0:
        raise VectorError(f"designed CDS length {len(cds)} is not a multiple of 3")


# -- features -------------------------------------------------------------


def _remap_features(
    backbone: VectorBackbone,
    remapper: IntervalRemapper,
    site: InsertionSite,
    *,
    protein_preserved: bool,
) -> tuple[tuple[Feature, ...], dict[str, tuple[Interval, ...]], tuple[DesignNote, ...]]:
    """Move every backbone feature into assembled coordinates.

    The `source` feature is rewritten to span the new length; the old CDS feature
    is dropped because it is replaced.

    A feature overlapping the insert is dropped, with one exception. When the
    redesign encodes the SAME protein, an in-frame protein-level feature wholly
    inside the CDS still describes the same residues, so it is kept -- the codon
    count is unchanged, so its coordinates are unchanged too. Dropping it would
    silently delete the Myc tag and the signal peptide from a user's map on a
    routine re-optimisation.

    Everything else overlapping the insert is still dropped rather than clipped,
    because clipping asserts a boundary the source file never contained.
    """
    out: list[Feature] = []
    parts: dict[str, tuple[Interval, ...]] = {}
    dropped: list[DesignNote] = []

    for feature in backbone.features:
        if feature.interval == site.interval and feature.kind.upper() == "CDS":
            continue
        if feature.kind.lower() == "source":
            out.append(replace(feature, interval=Interval(0, remapper.new_length)))
            continue

        moved = remapper.interval(feature.interval)
        if moved is None:
            if protein_preserved and _is_in_frame_inside(feature, site):
                out.append(feature)
                continue
            label = backbone.label_of(feature)
            if is_recombination_site(f"{label} {' '.join(feature.qualifiers.get('note', ()))}"):
                dropped.append(
                    DesignNote(
                        kind="liability",
                        summary=(
                            f"{label!r} is a recombination site inside the CDS being "
                            f"redesigned; back-translation changes those bases and the "
                            f"site will no longer function"
                        ),
                        interval=feature.interval,
                        bears_on="downstream cloning",
                        action=(
                            "exclude this span from the design, or expect a subsequent "
                            "recombination reaction using it to fail"
                        ),
                    )
                )
                continue
            dropped.append(
                DesignNote(
                    kind="change",
                    summary=(
                        f"dropped feature {label!r} ({feature.kind}): it overlaps the "
                        f"redesigned CDS, so its coordinates no longer describe anything"
                    ),
                    bears_on="map fidelity",
                )
            )
            continue

        original_parts = backbone.compound_parts.get(feature.uid)
        if original_parts:
            moved_parts = [remapper.interval(p) for p in original_parts]
            if any(p is None for p in moved_parts):
                dropped.append(
                    DesignNote(
                        kind="change",
                        summary=(
                            f"dropped multi-part feature {backbone.label_of(feature)!r}: one "
                            f"of its parts overlaps the replaced CDS"
                        ),
                        bears_on="map fidelity",
                    )
                )
                continue
            parts[feature.uid] = tuple(p for p in moved_parts if p is not None)
        out.append(replace(feature, interval=moved))

    return tuple(out), parts, tuple(dropped)


def _cds_feature(interval: Interval, site: InsertionSite, table_id: int) -> Feature:
    return Feature(
        interval=interval,
        kind="CDS",
        qualifiers={
            "label": (site.label,),
            "codon_start": ("1",),
            # The table is printed into the exported map so the record carries the
            # code it was designed under, not the reader's assumption.
            "transl_table": (str(table_id),),
        },
        uid="cds",
    )


def _remap_utr(utr: UtrContext, remapper: IntervalRemapper) -> UtrContext:
    rbs = utr.ribosome_binding_site
    return replace(
        utr,
        five_prime=remapper.interval(utr.five_prime) if utr.five_prime else None,
        three_prime=remapper.interval(utr.three_prime) if utr.three_prime else None,
        # The RBS is upstream of the insert, so it survives; remapping it here
        # keeps every interval on the context in ONE frame. A caller that had to
        # remember which of them were still in vector coordinates would get it
        # wrong on the first reverse-oriented cassette.
        ribosome_binding_site=remapper.interval(rbs) if rbs else None,
    )


# -- segments -------------------------------------------------------------


def _is_exempt(
    feature: Feature,
    label: str,
    exempt_kinds: Mapping[str, SegmentKind],
    exempt_labels: Sequence[str],
) -> SegmentKind | None:
    kind = exempt_kinds.get(feature.kind.lower())
    if kind is not None:
        return kind
    lowered = label.lower()
    if any(token in lowered for token in exempt_labels):
        return SegmentKind.WHITELISTED_REPEAT
    return None


def _segments(
    *,
    backbone: VectorBackbone,
    remapper: IntervalRemapper,
    features: Sequence[Feature],
    exempt_kinds: Mapping[str, SegmentKind],
    exempt_labels: Sequence[str],
    cds_strand: Strand,
) -> tuple[Segment, ...]:
    """Type every base of the construct, then coalesce runs into segments.

    A per-base pass rather than interval algebra: features overlap freely in real
    plasmid maps, and getting precedence right (CDS beats exempt beats backbone)
    once here is far easier to check than doing it with interval arithmetic at
    every call site.
    """
    n = remapper.new_length
    cds = remapper.insert_interval
    kinds: list[SegmentKind] = [SegmentKind.BACKBONE] * n
    labels: list[str] = [""] * n

    for feature in features:
        label = backbone.label_of(feature)
        kind = _is_exempt(feature, label, exempt_kinds, exempt_labels)
        if kind is None:
            continue
        for p in _positions(feature.interval, n):
            kinds[p] = kind
            labels[p] = label

    for p in _positions(cds, n):
        kinds[p] = SegmentKind.DESIGNABLE_CDS
        labels[p] = "cds"

    out: list[Segment] = []
    start = 0
    for i in range(1, n + 1):
        if i == n or kinds[i] != kinds[start] or labels[i] != labels[start]:
            # The CDS segment carries the cassette's strand, so `Construct.editable`
            # slices back to the coding sequence rather than the plus strand.
            strand: Strand = cds_strand if kinds[start] is SegmentKind.DESIGNABLE_CDS else 1
            out.append(Segment(Interval(start, i, strand), kinds[start], labels[start]))
            start = i
    return tuple(out)


def _positions(iv: Interval, n: int) -> Sequence[int]:
    """Every construct position an interval covers, wrapping where it must."""
    return [p % n for p in range(iv.start, iv.end)]
