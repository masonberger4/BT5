"""Exact repeats over the assembled construct, and the risk surface.

The cases that matter are the ones a length cutoff gets wrong: a short TANDEM
array (real, and below the vendor threshold), a long repeat whose copies are far
apart (real, and distance does not make it low), and the doubled-sequence
artefact that would otherwise report every k-mer in a circular plasmid.
"""

from __future__ import annotations

import numpy as np
import pytest
from bt5.core.types import Construct, Interval, Segment, SegmentKind, Topology
from bt5.vector import ConstructKmerIndex, repeat_risk


def dna(n: int, *, seed: int = 3) -> str:
    """Non-repetitive filler. A repetitive pad would BE the finding."""
    rng = np.random.default_rng(seed)
    return "".join("ACGT"[i] for i in rng.integers(0, 4, n))


def construct(seq: str, *, circular: bool = True, exempt: Interval | None = None) -> Construct:
    segments = [Segment(Interval(0, len(seq)), SegmentKind.BACKBONE)]
    if exempt is not None:
        segments = [
            Segment(Interval(0, exempt.start), SegmentKind.BACKBONE),
            Segment(exempt, SegmentKind.WHITELISTED_REPEAT, "ITR"),
            Segment(Interval(exempt.end, len(seq)), SegmentKind.BACKBONE),
        ]
    return Construct(
        sequence=seq,
        topology=Topology.CIRCULAR if circular else Topology.LINEAR,
        segments=tuple(segments),
    )


class TestDuplicates:
    def test_finds_a_planted_repeat(self) -> None:
        unit = dna(40, seed=9)
        seq = dna(200) + unit + dna(300, seed=5) + unit + dna(200, seed=7)
        pairs = ConstructKmerIndex.of(construct(seq), 20).repeat_pairs(20)
        assert any(p.length >= 40 for p in pairs)

    def test_a_circular_construct_does_not_report_every_kmer(self) -> None:
        """The doubling artefact: every position matches itself one length away."""
        seq = dna(2000)
        assert ConstructKmerIndex.of(construct(seq), 20).repeat_pairs(20) == []

    def test_finds_a_repeat_spanning_the_origin(self) -> None:
        unit = dna(40, seed=11)
        seq = unit[20:] + dna(400) + unit + dna(400, seed=13) + unit[:20]
        pairs = ConstructKmerIndex.of(construct(seq), 20).repeat_pairs(20)
        assert pairs, "a repeat straddling position 0 must still be found"

    def test_a_tandem_array_reports_its_period(self) -> None:
        unit = dna(30, seed=17)
        seq = dna(200) + unit * 4 + dna(200, seed=19)
        best = ConstructKmerIndex.of(construct(seq), 20).repeat_pairs(20)[0]
        assert best.tandem
        assert best.length == 30, "the period is interpretable; a smeared extension is not"

    def test_whitelisted_regions_can_be_excluded(self) -> None:
        """ITRs and LTRs are an accepted design feature, reported separately."""
        unit = dna(60, seed=23)
        seq = dna(100) + unit + dna(100, seed=29) + unit + dna(100, seed=31)
        c = construct(seq, exempt=Interval(100, 160))
        assert ConstructKmerIndex.of(c, 20).repeat_pairs(20)
        both = Interval(100, 160 + 100 + 60)
        assert not [
            p
            for p in ConstructKmerIndex.of(c, 20).repeat_pairs(20, exclude=(both,))
            if both.start <= p.first.start and p.second.end <= both.end
        ]

    def test_the_protocol_method_yields_interval_pairs(self) -> None:
        unit = dna(40, seed=37)
        seq = dna(200) + unit + dna(300, seed=41) + unit + dna(200, seed=43)
        first, second = next(iter(ConstructKmerIndex.of(construct(seq), 20).duplicates(20)))
        assert isinstance(first, Interval)
        assert isinstance(second, Interval)


class TestRiskSurface:
    def test_a_short_distant_repeat_is_low(self) -> None:
        assert repeat_risk(15, 5000) == "low"

    def test_a_short_tandem_repeat_is_not_low(self) -> None:
        """Slipped-strand mispairing needs no loop, so the length floor does not apply."""
        assert repeat_risk(15, 0, tandem=True) == "moderate"

    def test_proximity_dominates_in_the_reca_independent_regime(self) -> None:
        assert repeat_risk(30, 10) == "high"
        assert repeat_risk(30, 5000) == "low"

    def test_a_substantial_repeat_is_never_low_however_far_apart(self) -> None:
        """189 bp of identity between two LTRs is real at any spacing."""
        assert repeat_risk(189, 3673) == "moderate"

    @pytest.mark.parametrize(("length", "helps"), [(30, False), (150, False), (250, True)])
    def test_reca_only_helps_above_the_dependence_floor(self, length: int, helps: bool) -> None:
        unit = dna(length, seed=47)
        seq = dna(300) + unit + dna(300, seed=53) + unit + dna(300, seed=59)
        pair = next(
            p
            for p in ConstructKmerIndex.of(construct(seq), 20).repeat_pairs(20)
            if p.length >= length
        )
        assert pair.reca_strain_helps is helps


class TestInvertedRepeats:
    def test_finds_a_hairpin_stem(self) -> None:
        from bt5.core.types import reverse_complement

        stem = dna(25, seed=61)
        seq = dna(200) + stem + dna(30, seed=67) + reverse_complement(stem) + dna(200, seed=71)
        found = list(ConstructKmerIndex.of(construct(seq), 20).revcomp_pairs(20, 60))
        assert found
        assert any(b.strand == -1 for _, b in found)


class TestBiosecurity:
    def test_the_only_constructor_takes_a_construct(self) -> None:
        """No path accepts an external database; see the module docstring."""
        import inspect

        params = list(inspect.signature(ConstructKmerIndex.of).parameters)
        assert params == ["c", "k"]
