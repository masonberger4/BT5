"""Codon distance -- the space G4's gallery is selected in.

The choice of space is the finding, not a detail. Two designs can differ at 30%
of their BASES while every difference is third-position wobble, and a gallery
picked on base distance will present five of those as five options.
"""

from __future__ import annotations

import pytest
from bt5.score import codon_distance, distance_matrix, nucleotide_distance, pairwise_minimum


class TestCodonDistance:
    def test_identical_sequences_are_zero_apart(self) -> None:
        assert codon_distance("ATGAAACCC", "ATGAAACCC") == 0.0

    def test_every_codon_differing_is_one(self) -> None:
        assert codon_distance("ATGAAACCC", "TTTGGGTTT") == 1.0

    def test_it_counts_codons_not_bases(self) -> None:
        # One codon differs, by one base, out of three codons.
        assert codon_distance("ATGAAACCC", "ATGAAGCCC") == pytest.approx(1 / 3)

    def test_wobble_is_why_the_space_matters(self) -> None:
        """Six synonymous third-position changes: 22% of bases, 100% of codons.

        A gallery selected on nucleotide distance would rank this pair as barely
        different; in codon space it is maximally different, which is what the
        user is actually choosing between."""
        a = "CTTCTTCTTCTTCTTCTT"
        b = "CTGCTGCTGCTGCTGCTG"
        assert nucleotide_distance(a, b) == pytest.approx(6 / 18)
        assert codon_distance(a, b) == 1.0

    def test_a_length_mismatch_is_an_error_not_a_truncation(self) -> None:
        """Different lengths do not encode the same protein, so a distance
        between them is not the quantity anything here wants."""
        with pytest.raises(ValueError, match="different lengths"):
            codon_distance("ATGAAA", "ATGAAACCC")

    def test_a_partial_codon_is_refused(self) -> None:
        with pytest.raises(ValueError, match="whole number of codons"):
            codon_distance("ATGAA", "ATGAA")

    def test_empty_is_zero_not_an_error(self) -> None:
        assert codon_distance("", "") == 0.0


class TestPairwiseMinimum:
    def test_it_reports_the_closest_pair_not_the_average(self) -> None:
        """A gallery is only as diverse as its most similar pair: four good
        designs plus a near-duplicate is a gallery of four, and a mean hides it."""
        far_a, far_b = "ATGAAACCC", "TTTGGGTTT"
        near = "ATGAAACCG"  # one codon from far_a
        assert pairwise_minimum([far_a, far_b, near]) == pytest.approx(1 / 3)

    def test_a_single_design_is_trivially_diverse(self) -> None:
        assert pairwise_minimum(["ATGAAACCC"]) == 1.0
        assert pairwise_minimum([]) == 1.0


class TestDistanceMatrix:
    def test_it_is_symmetric_with_a_zero_diagonal(self) -> None:
        seqs = ["ATGAAACCC", "ATGAAGCCC", "TTTGGGTTT"]
        m = distance_matrix(seqs)
        assert all(m[i][i] == 0.0 for i in range(3))
        assert all(m[i][j] == m[j][i] for i in range(3) for j in range(3))

    def test_it_agrees_with_the_pairwise_function(self) -> None:
        seqs = ["ATGAAACCC", "ATGAAGCCC", "TTTGGGTTT"]
        m = distance_matrix(seqs)
        assert m[0][1] == pytest.approx(codon_distance(seqs[0], seqs[1]))
