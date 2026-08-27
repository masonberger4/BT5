"""Genetic codes and CAI."""

from __future__ import annotations

import pytest
from bt5.codon.tables import FileTableProvider, NcbiGeneticCode


@pytest.fixture(scope="module")
def provider() -> FileTableProvider:
    return FileTableProvider()


@pytest.fixture(scope="module")
def ecoli(provider: FileTableProvider) -> NcbiGeneticCode:
    return provider.genetic_code(11)


class TestGeneticCode:
    def test_synonymous_sets_exclude_stop_codons(self, provider: FileTableProvider) -> None:
        """NCBI tables 27/28 make TGA BOTH Trp and a stop. Offering TGA for Trp
        would terminate translation mid-protein."""
        for table_id in (27, 28):
            code = provider.genetic_code(table_id)
            assert "TGA" in code.stop_codons
            assert code.synonymous_codons("W") == ("TGG",), (
                f"table {table_id} must not offer TGA for Trp"
            )

    def test_every_table_offers_only_non_stop_codons(self, provider: FileTableProvider) -> None:
        from Bio.Data import CodonTable

        for table_id in sorted(CodonTable.unambiguous_dna_by_id):
            code = provider.genetic_code(table_id)
            stops = set(code.stop_codons)
            for aa in "ACDEFGHIKLMNPQRSTVWY":
                try:
                    codons = code.synonymous_codons(aa)
                except ValueError:
                    continue  # amino acid unencodable in this table
                assert not (set(codons) & stops), f"table {table_id}, {aa}: stop codon offered"

    def test_table_12_reassigns_ctg_to_serine(self, provider: FileTableProvider) -> None:
        """The canonical silent-wrong-protein bug: CTG is Leu in table 1 and Ser
        in table 12 (alternative yeast nuclear)."""
        assert provider.genetic_code(1).translate("CTG") == "L"
        assert provider.genetic_code(12).translate("CTG") == "S"

    def test_translate_rejects_a_bad_frame(self, ecoli: NcbiGeneticCode) -> None:
        with pytest.raises(ValueError, match="not a multiple of 3"):
            ecoli.translate("ATGA")

    def test_unknown_table_is_rejected(self, provider: FileTableProvider) -> None:
        with pytest.raises(ValueError, match="unknown NCBI translation table"):
            provider.genetic_code(999)


class TestCai:
    def test_preferred_codons_score_higher_than_rare_ones(
        self, provider: FileTableProvider, ecoli: NcbiGeneticCode
    ) -> None:
        u = provider.usage("sharp_li_1987_ecoli_w")
        assert u.cai("ATG" + "CTG" * 10, ecoli) > 0.99  # CTG: preferred Leu
        assert u.cai("ATG" + "CTA" * 10, ecoli) < 0.05  # CTA: rare Leu

    def test_known_rare_arginine_codons(self, provider: FileTableProvider) -> None:
        """AGA/AGG are the rare Arg codons that stall E. coli expression and
        motivate Rosetta/pRARE strains."""
        w = provider.usage("sharp_li_1987_ecoli_w").w
        assert w["CGT"] == 1.0
        assert w["AGA"] < 0.01
        assert w["AGG"] < 0.01

    def test_single_codon_families_are_excluded(
        self, provider: FileTableProvider, ecoli: NcbiGeneticCode
    ) -> None:
        """Met and Trp carry no information, so a poly-Met peptide has no
        informative codons and must degrade rather than raise."""
        u = provider.usage("sharp_li_1987_ecoli_w")
        assert u.cai("ATGATGATG", ecoli) == 0.0

    def test_reference_set_is_recorded(self, provider: FileTableProvider) -> None:
        """CAI is meaningless without a pinned reference set: swapping it raises
        every CAI value, which a benchmark would read as an improvement."""
        u = provider.usage("sharp_li_1987_ecoli_w")
        assert "Sharp" in u.reference_set
        assert u.citation_url.startswith("https://")

    def test_stop_codons_absent_from_the_index(self, provider: FileTableProvider) -> None:
        """CAI is a geometric mean over sense codons only."""
        w = provider.usage("sharp_li_1987_ecoli_w").w
        assert len(w) == 61
        for stop in ("TAA", "TAG", "TGA"):
            assert stop not in w

    def test_missing_host_is_a_clear_error(self, provider: FileTableProvider) -> None:
        with pytest.raises(FileNotFoundError, match="no codon usage table"):
            provider.usage("no_such_host")
