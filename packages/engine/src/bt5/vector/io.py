"""Vector file I/O: GenBank, SnapGene, FASTA in; GenBank out.

SnapGene `.dna` is READ-ONLY and always will be. Every open implementation of the
format is read-only -- Biopython handles only a few of the packet types -- so BT5
imports `.dna` and exports `.gb`. Promising a `.dna` round trip would be promising
something no open tool can deliver.

What "round trip" means here is worth stating precisely, because the honest claim
is narrower than the obvious one. Byte-identity with an ARBITRARY input file is
not achievable: Biopython's writer normalises LOCUS spacing, qualifier wrapping
and sequence case, so a file from SnapGene or Benchling comes back formatted
differently no matter what BT5 does. What is guaranteed, and tested, is:

  * semantic round trip -- sequence, topology, and every feature's location,
    key and qualifiers survive read -> write -> read unchanged; and
  * writer idempotence -- write(read(write(x))) is byte-identical to write(x),
    so BT5's own output is a fixed point and diffs between two BT5 exports show
    only real changes.
"""

from __future__ import annotations

import io as _io
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import IO, Any, cast

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from bt5.core.types import Construct, Feature, Interval, Topology
from bt5.vector.backbone import VectorBackbone, VectorError
from bt5.vector.locations import interval_to_location, location_to_interval, parts_to_location
from bt5.vector.notes import DesignNote

#: Record-level annotations carried through a round trip. Biopython preserves
#: all of these on read; anything left off this list is silently dropped on
#: export, which is how a user's plasmid map comes back with an empty ORGANISM
#: and a 1980 date. `date` and `data_file_division` populate the LOCUS line.
_SCALAR_ANNOTATIONS = (
    "molecule_type",
    "topology",
    "data_file_division",
    "date",
    "source",
    "organism",
    "definition",
    "comment",
)

#: `Construct.annotations` is a Mapping[str, str], so GenBank's list-valued
#: fields are joined on read and split again on write. "; " is GenBank's own
#: separator for all three. Stringifying the list instead makes the writer
#: non-idempotent, each pass re-quoting the last.
_LIST_ANNOTATIONS = ("keywords", "taxonomy", "accessions")
_LIST_SEPARATOR = "; "

#: Values GenBank uses to mean "absent"; carrying them forward would turn a
#: blank field into a literal dot on the next pass.
_EMPTY_VALUES = frozenset({"", "."})


def _seqfeature(location: Any, kind: str, qualifiers: dict[str, list[str]]) -> SeqFeature:
    # Biopython's constructor carries no annotations.
    return SeqFeature(location=location, type=kind, qualifiers=qualifiers)  # type: ignore[no-untyped-call]


def _annotations(record: SeqRecord) -> dict[str, Any]:
    """SeqRecord.annotations is typed `dict[str, str | int]` but genuinely holds
    lists at runtime (KEYWORDS, REFERENCES). One accessor keeps that discrepancy
    from leaking into the rest of the lane."""
    return cast("dict[str, Any]", record.annotations)


def _carry_annotations(source: Mapping[str, Any]) -> dict[str, str]:
    """Flatten a SeqRecord's record-level annotations into strings."""
    out: dict[str, str] = {}
    for key in _SCALAR_ANNOTATIONS:
        value = source.get(key)
        if value is not None and str(value) not in _EMPTY_VALUES:
            out[key] = str(value)
    for key in _LIST_ANNOTATIONS:
        value = source.get(key)
        parts = value if isinstance(value, list) else [value]
        joined = _LIST_SEPARATOR.join(
            str(v) for v in parts if v is not None and str(v) not in _EMPTY_VALUES
        )
        if joined:
            out[key] = joined
    return out


def _restore_annotations(record: SeqRecord, annotations: Mapping[str, str]) -> None:
    """The inverse of `_carry_annotations`, onto a record being written."""
    out = _annotations(record)
    for key in _SCALAR_ANNOTATIONS:
        if key == "definition":
            continue  # travels as SeqRecord.description
        if annotations.get(key):
            out[key] = annotations[key]
    for key in _LIST_ANNOTATIONS:
        if annotations.get(key):
            out[key] = annotations[key].split(_LIST_SEPARATOR)


def _read(handle: Any, fmt: str) -> SeqRecord:
    """Biopython's SeqIO.read is unannotated; this is the lane's only call site."""
    return cast("SeqRecord", SeqIO.read(handle, fmt))  # type: ignore[no-untyped-call]


def _topology_of(record: SeqRecord) -> Topology:
    """Topology is never guessed.

    A vector silently treated as linear loses every origin-spanning motif, and a
    linear fragment treated as circular invents a junction that does not exist.
    """
    value = str(_annotations(record).get("topology", "")).lower()
    if value == "circular":
        return Topology.CIRCULAR
    if value == "linear":
        return Topology.LINEAR
    raise VectorError(
        f"{record.id!r} does not declare a topology; GenBank records BT5 accepts "
        f"must say 'circular' or 'linear' on the LOCUS line"
    )


def backbone_from_record(record: SeqRecord, *, topology: Topology | None = None) -> VectorBackbone:
    """Convert a parsed SeqRecord into a VectorBackbone."""
    topo = topology if topology is not None else _topology_of(record)
    sequence = str(record.seq).upper()
    bad = sorted(set(sequence) - set("ACGT"))
    if bad:
        raise VectorError(
            f"{record.id!r} contains non-ACGT characters {bad}; BT5 needs an "
            f"unambiguous sequence to reason about motifs and folding"
        )

    length = len(sequence)
    circular = topo is Topology.CIRCULAR
    features: list[Feature] = []
    compound: dict[str, tuple[Interval, ...]] = {}
    notes: list[DesignNote] = []

    for index, feature in enumerate(record.features):
        uid = f"f{index:04d}"
        if feature.location is None:
            notes.append(
                DesignNote(
                    kind="change",
                    summary=f"dropped feature {index} ({feature.type}): it has no location",
                    bears_on="map fidelity",
                )
            )
            continue
        try:
            parsed = location_to_interval(feature.location, length=length, circular=circular)
        except ValueError as exc:
            notes.append(
                DesignNote(
                    kind="change",
                    summary=f"dropped feature {index} ({feature.type}): {exc}",
                    bears_on="map fidelity",
                )
            )
            continue
        qualifiers = {
            key: tuple(str(v) for v in values)
            for key, values in cast("Mapping[str, Iterable[Any]]", feature.qualifiers).items()
        }
        features.append(
            Feature(
                interval=parsed.interval, kind=str(feature.type), qualifiers=qualifiers, uid=uid
            )
        )
        if parsed.is_compound:
            compound[uid] = parsed.parts

    annotations = _carry_annotations(_annotations(record))
    annotations.setdefault("molecule_type", "DNA")
    annotations["topology"] = topo.value
    if record.description:
        annotations.setdefault("definition", record.description)

    return VectorBackbone(
        sequence=sequence,
        topology=topo,
        features=tuple(features),
        annotations=annotations,
        name=record.name or record.id or "vector",
        compound_parts=compound,
        notes=tuple(notes),
    )


def read_genbank(source: str | Path | IO[str]) -> VectorBackbone:
    """Parse a GenBank vector map."""
    handle = _open_text(source)
    record = _read(handle, "genbank")
    return backbone_from_record(record)


def read_snapgene(source: str | Path) -> VectorBackbone:
    """Import a SnapGene `.dna` map. Read-only by design; export is GenBank."""
    with Path(source).open("rb") as fh:
        record = _read(fh, "snapgene")
    return backbone_from_record(record)


def read_fasta(source: str | Path | IO[str], *, topology: Topology) -> VectorBackbone:
    """Parse a bare sequence. Topology is REQUIRED because FASTA cannot say."""
    handle = _open_text(source)
    record = _read(handle, "fasta")
    return backbone_from_record(record, topology=topology)


def read_vector(path: str | Path, *, topology: Topology | None = None) -> VectorBackbone:
    """Dispatch on file extension."""
    p = Path(path)
    suffix = p.suffix.lower()
    if suffix in (".gb", ".gbk", ".genbank", ".ape"):
        return read_genbank(p)
    if suffix == ".dna":
        return read_snapgene(p)
    if suffix in (".fa", ".fasta", ".fna", ".seq", ".txt"):
        if topology is None:
            raise VectorError(
                f"{p.name} is FASTA, which cannot record topology; pass "
                f"topology=Topology.CIRCULAR or Topology.LINEAR explicitly"
            )
        return read_fasta(p, topology=topology)
    raise VectorError(f"unrecognised vector file extension {suffix!r}")


def _open_text(source: str | Path | IO[str]) -> IO[str]:
    if isinstance(source, str | Path):
        path = Path(source)
        if path.exists():
            # Read eagerly rather than handing back an open handle: the caller
            # has no way to close it, and a leaked descriptor per parsed vector
            # is the kind of thing that only shows up under a batch run.
            return _io.StringIO(path.read_text(encoding="utf-8"))
        if isinstance(source, str):
            return _io.StringIO(source)
        raise VectorError(f"no such vector file: {path}")
    return source


# -- export ---------------------------------------------------------------


def _record_from(
    *,
    sequence: str,
    topology: Topology,
    features: Iterable[Feature],
    compound_parts: Mapping[str, tuple[Interval, ...]],
    annotations: Mapping[str, str],
    name: str,
) -> SeqRecord:
    length = len(sequence)
    accessions = annotations.get("accessions", "")
    record = SeqRecord(
        Seq(sequence),
        id=accessions.split(_LIST_SEPARATOR)[0] if accessions else ".",
        name=name[:16] or "vector",
        description=annotations.get("definition", ""),
    )
    _restore_annotations(record, annotations)
    out_annotations = _annotations(record)
    out_annotations["molecule_type"] = annotations.get("molecule_type", "DNA")
    out_annotations["topology"] = topology.value

    out: list[SeqFeature] = []
    for feature in features:
        parts = compound_parts.get(feature.uid)
        location = (
            parts_to_location(parts)
            if parts
            else interval_to_location(feature.interval, length=length)
        )
        out.append(
            _seqfeature(location, feature.kind, {k: list(v) for k, v in feature.qualifiers.items()})
        )
    record.features = out
    return record


def backbone_to_record(backbone: VectorBackbone) -> SeqRecord:
    return _record_from(
        sequence=backbone.sequence,
        topology=backbone.topology,
        features=backbone.features,
        compound_parts=backbone.compound_parts,
        annotations=backbone.annotations,
        name=backbone.name,
    )


def construct_to_record(
    construct: Construct,
    *,
    name: str = "bt5_design",
    compound_parts: Mapping[str, tuple[Interval, ...]] | None = None,
) -> SeqRecord:
    return _record_from(
        sequence=construct.sequence,
        topology=construct.topology,
        features=construct.features,
        compound_parts=compound_parts or {},
        annotations=construct.annotations,
        name=name,
    )


def write_genbank(record: SeqRecord, target: str | Path | IO[str] | None = None) -> str:
    """Write GenBank. Returns the text; also writes to `target` when given."""
    handle = _io.StringIO()
    SeqIO.write(record, handle, "genbank")
    text = handle.getvalue()
    if target is not None:
        if isinstance(target, str | Path):
            Path(target).write_text(text, encoding="utf-8")
        else:
            target.write(text)
    return text
