"""Vector file I/O and what "round trip" actually guarantees.

Byte-identity with an ARBITRARY input file is not achievable through Biopython's
writer, which normalises LOCUS spacing, qualifier wrapping and sequence case. The
two claims that ARE made, and tested here, are semantic round trip and writer
idempotence.
"""

from __future__ import annotations

import io
from pathlib import Path

import pytest
from bt5.core.types import Interval, Topology
from bt5.vector import (
    VectorError,
    backbone_to_record,
    read_fasta,
    read_genbank,
    read_vector,
    write_genbank,
)
from bt5.vector.backbone import VectorBackbone

MINIMAL_LINEAR = """LOCUS       tiny                      12 bp    DNA     linear   SYN 01-JAN-2026
FEATURES             Location/Qualifiers
     misc_feature    1..6
                     /label="left"
ORIGIN
        1 acgtacgtac gt
//
"""

NO_TOPOLOGY = MINIMAL_LINEAR.replace("DNA     linear  ", "DNA             ")


class TestParsing:
    def test_reads_sequence_and_topology(self, backbone: VectorBackbone) -> None:
        assert backbone.length == 2000
        assert backbone.topology is Topology.CIRCULAR
        assert set(backbone.sequence) <= set("ACGT"), "sequence is normalised to upper case"

    def test_reads_every_feature(self, backbone: VectorBackbone) -> None:
        kinds = [f.kind for f in backbone.features]
        assert kinds.count("LTR") == 2
        assert "5'UTR" in kinds
        assert "intron" in kinds
        assert kinds.count("CDS") == 1

    def test_origin_spanning_join_becomes_a_wrapping_interval(
        self, backbone: VectorBackbone
    ) -> None:
        ori = next(f for f in backbone.features if f.kind == "rep_origin")
        assert ori.interval == Interval(1800, 2080, 1)
        assert ori.interval.end > backbone.length, "the wrap is the whole point"
        assert ori.interval.length == 280

    def test_slicing_a_wrapping_feature_crosses_the_origin(self, backbone: VectorBackbone) -> None:
        ori = next(f for f in backbone.features if f.kind == "rep_origin")
        text = backbone.slice(ori.interval)
        assert len(text) == 280
        assert text == backbone.sequence[1800:] + backbone.sequence[:80]

    def test_topology_is_never_guessed(self) -> None:
        """A vector silently read as linear loses every origin-spanning motif."""
        with pytest.raises(VectorError, match="does not declare a topology"):
            read_genbank(io.StringIO(NO_TOPOLOGY))

    def test_ambiguous_bases_are_refused(self) -> None:
        record = MINIMAL_LINEAR.replace("acgtacgtac gt", "acgtacgtan gt")
        with pytest.raises(VectorError, match="non-ACGT"):
            read_genbank(io.StringIO(record))


class TestRoundTrip:
    def test_semantic_round_trip(self, backbone: VectorBackbone) -> None:
        text = write_genbank(backbone_to_record(backbone))
        again = read_genbank(io.StringIO(text))
        assert again.sequence == backbone.sequence
        assert again.topology == backbone.topology
        assert again.features == backbone.features

    def test_writer_is_idempotent(self, backbone: VectorBackbone) -> None:
        """BT5's own output is a fixed point, so a diff of two exports is real."""
        first = write_genbank(backbone_to_record(backbone))
        second = write_genbank(backbone_to_record(read_genbank(io.StringIO(first))))
        assert first == second

    def test_the_origin_spanning_join_is_written_back_as_a_join(
        self, backbone: VectorBackbone
    ) -> None:
        text = write_genbank(backbone_to_record(backbone))
        assert "join(1801..2000,1..80)" in text


class TestOtherFormats:
    def test_fasta_requires_an_explicit_topology(self, tmp_path: Path) -> None:
        path = tmp_path / "v.fasta"
        path.write_text(">v\nACGTACGTACGT\n")
        with pytest.raises(VectorError, match="cannot record topology"):
            read_vector(path)

    def test_fasta_with_a_topology_parses(self, tmp_path: Path) -> None:
        path = tmp_path / "v.fasta"
        path.write_text(">v\nACGTACGTACGT\n")
        backbone = read_fasta(path, topology=Topology.CIRCULAR)
        assert backbone.sequence == "ACGTACGTACGT"
        assert backbone.is_circular
        assert backbone.features == ()

    def test_read_vector_dispatches_on_extension(self, backbone_path: Path) -> None:
        assert read_vector(backbone_path).length == 2000

    def test_unknown_extension_is_refused(self, tmp_path: Path) -> None:
        path = tmp_path / "v.xyz"
        path.write_text("nonsense")
        with pytest.raises(VectorError, match="unrecognised vector file extension"):
            read_vector(path)
