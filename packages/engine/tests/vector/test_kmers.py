"""Exact repeats over the assembled construct, and the risk surface.

The cases that matter are the ones a length cutoff gets wrong: a short TANDEM
array (real, and below the vendor threshold), a long repeat whose copies are far
apart (real, and distance does not make it low), and the doubled-sequence
artefact that would otherwise report every k-mer in a circular plasmid.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from bt5.core.types import (
    Construct,
    Interval,
    Segment,
    SegmentKind,
    Topology,
    reverse_complement,
)
from bt5.vector import ConstructKmerIndex, kmers, repeat_risk


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

    def test_two_matches_on_the_same_diagonal_are_both_found(self) -> None:
        """The memo holds one grown span per diagonal `b - a`, so two unrelated
        repeat pairs that happen to share an offset must not shadow each other.

        This is the one way the memo could lose a finding, so it is the one case
        worth pinning: both pairs are 40 bp with their copies 140 bp apart, in
        different parts of the construct.
        """
        first, second = dna(40, seed=61), dna(40, seed=67)
        seq = (
            first
            + dna(100, seed=71)
            + first
            + dna(220, seed=73)
            + second
            + dna(100, seed=79)
            + second
            + dna(100, seed=83)
        )
        pairs = ConstructKmerIndex.of(construct(seq), 20).repeat_pairs(20)
        starts = {p.first.start for p in pairs if p.length >= 40}
        assert starts == {0, 400}, f"expected both diagonal-140 pairs, got {sorted(starts)}"

    def test_a_kmer_at_three_sites_reports_every_adjacent_pair(self) -> None:
        unit = dna(30, seed=89)
        seq = unit + dna(90, seed=97) + unit + dna(90, seed=101) + unit
        pairs = ConstructKmerIndex.of(construct(seq, circular=False), 20).repeat_pairs(20)
        assert {p.first.start for p in pairs if p.length >= 30} == {0, 120}

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
    """A hairpin is not a direct repeat with a minus sign on it.

    A direct repeat is lost by deletion and codon choice fixes it. An inverted
    repeat extrudes a cruciform, stalls forks and is cleaved by SbcCD -- a
    different mechanism with a different answer, which is why it has its own
    type and stays out of `repeat_risk`.
    """

    def hairpin(self, arm: str, loop: int, *, lead: int = 200, tail: int = 200) -> str:
        return dna(lead) + arm + dna(loop, seed=67) + reverse_complement(arm) + dna(tail, seed=71)

    def test_finds_a_hairpin_stem(self) -> None:
        found = ConstructKmerIndex.of(
            construct(self.hairpin(dna(25, seed=61), 30)), 20
        ).inverted_repeats(20, 60)
        assert found
        assert all(p.second.strand == -1 for p in found)
        assert all(p.first.strand == 1 for p in found)

    def test_the_stem_is_reported_at_its_full_length(self) -> None:
        """Seeding at 20 and reporting the seed calls a 60 bp stem a 20 bp one."""
        arm = dna(60, seed=81)
        found = ConstructKmerIndex.of(construct(self.hairpin(arm, 10)), 20).inverted_repeats(20, 60)
        assert len(found) == 1, "one hairpin is one finding, not one per seed offset"
        assert found[0].stem >= 60
        assert found[0].loop == 10

    def test_one_hairpin_is_one_finding(self) -> None:
        """Without maximal extension this reported 26 near-identical pairs."""
        found = ConstructKmerIndex.of(
            construct(self.hairpin(dna(60, seed=81), 10)), 20
        ).inverted_repeats(20, 60)
        assert len(found) == 1

    def test_an_origin_spanning_hairpin_is_found(self) -> None:
        """The 3' arm lives past position n, so indexing only the first turn
        loses it -- silently, and on exactly the ITR layouts that matter."""
        arm = dna(30, seed=83)
        hairpin = arm + dna(10, seed=89) + reverse_complement(arm)
        cut = 35  # split the hairpin so it straddles position 0
        seq = hairpin[cut:] + dna(400, seed=97) + hairpin[:cut]
        found = ConstructKmerIndex.of(construct(seq), 20).inverted_repeats(20, 60)
        assert found, "an origin-spanning stem-loop is still a stem-loop"
        assert found[0].stem >= 20

    def test_the_same_hairpin_is_found_wherever_the_origin_sits(self) -> None:
        arm = dna(30, seed=83)
        hairpin = arm + dna(10, seed=89) + reverse_complement(arm)
        middle = dna(200) + hairpin + dna(200, seed=97)
        rotated = middle[240:] + middle[:240]
        a = ConstructKmerIndex.of(construct(middle), 20).inverted_repeats(20, 60)
        b = ConstructKmerIndex.of(construct(rotated), 20).inverted_repeats(20, 60)
        assert len(a) == len(b) == 1
        assert a[0].stem == b[0].stem
        assert a[0].loop == b[0].loop

    def test_a_perfect_palindrome_has_no_loop(self) -> None:
        arm = dna(40, seed=101)
        seq = dna(200) + arm + reverse_complement(arm) + dna(200, seed=103)
        found = ConstructKmerIndex.of(construct(seq), 20).inverted_repeats(20, 60)
        assert found[0].loop == 0
        assert found[0].perfect_palindrome

    def test_the_loop_ceiling_is_enforced(self) -> None:
        arm = dna(25, seed=107)
        assert (
            ConstructKmerIndex.of(construct(self.hairpin(arm, 200)), 20).inverted_repeats(20, 60)
            == []
        )

    def test_nested_stems_collapse_to_the_maximal_one(self) -> None:
        """A homopolymer arm pairs at every offset, so one physical stem seeds
        41 nested alignments. Reporting them all buries the finding in itself."""
        # Explicit non-pairing flanks: a stray A or T beside the arm gives the
        # homopolymer a second valid register, which is a real alternative stem
        # rather than a duplicate, and not what this test is about.
        seq = dna(149) + "G" + "A" * 40 + "T" * 40 + "G" + dna(149, seed=137)
        found = ConstructKmerIndex.of(construct(seq), 20).inverted_repeats(20, 60)
        assert len(found) == 1
        assert found[0].stem == 40
        assert found[0].perfect_palindrome

    def test_a_stem_is_found_from_its_outermost_seed_alone(self) -> None:
        """The invariant the inward-only extension rests on: every seed position
        is indexed, so the outermost one is always available and closing inward
        from it recovers the whole stem. Checked where the stem touches position
        0, which is where an outward pass would have been the only rescue."""
        arm = dna(50, seed=139)
        seq = arm + dna(10, seed=149) + reverse_complement(arm) + dna(300, seed=151)
        found = ConstructKmerIndex.of(construct(seq), 20).inverted_repeats(20, 60)
        assert len(found) == 1
        assert found[0].first.start == 0
        assert found[0].stem >= 50

    def test_an_exempt_region_can_be_excluded(self) -> None:
        """AAV ITRs are palindromic by construction; reporting them is noise."""
        arm = dna(60, seed=81)
        seq = self.hairpin(arm, 10)
        # An annotated ITR feature comfortably contains its palindrome; both
        # arms must fall inside for the exclusion to apply, as for direct repeats.
        itr = Interval(150, 380)
        index = ConstructKmerIndex.of(construct(seq), 20)
        assert index.inverted_repeats(20, 60)
        assert index.inverted_repeats(20, 60, exclude=[itr]) == []

    def test_a_direct_repeat_is_not_an_inverted_one(self) -> None:
        """The two scans must not leak into each other: only the direct one
        feeds a risk surface calibrated on the deletion literature."""
        unit = dna(40, seed=109)
        direct = dna(200) + unit + dna(300, seed=113) + unit + dna(200, seed=127)
        assert ConstructKmerIndex.of(construct(direct), 20).inverted_repeats(20, 400) == []

    def test_an_inverted_repeat_is_not_a_direct_one(self) -> None:
        arm = dna(40, seed=109)
        seq = dna(200) + arm + dna(300, seed=113) + reverse_complement(arm) + dna(200, seed=127)
        assert ConstructKmerIndex.of(construct(seq), 20).repeat_pairs(20) == []

    def test_the_arms_really_base_pair(self) -> None:
        """Checked against the sequence, not against the coordinates that
        produced it: a stem whose arms do not pair is not a stem."""
        arm = dna(45, seed=131)
        seq = self.hairpin(arm, 12)
        c = construct(seq)
        for pair in ConstructKmerIndex.of(c, 20).inverted_repeats(20, 60):
            assert c.slice(pair.first) == reverse_complement(
                seq[pair.second.start : pair.second.end]
            )
            assert pair.first.length == pair.second.length == pair.stem


class TestBiosecurity:
    def test_the_only_constructor_takes_a_construct(self) -> None:
        """No path accepts an external database; see the module docstring."""
        import inspect

        params = list(inspect.signature(ConstructKmerIndex.of).parameters)
        assert params == ["c", "k"]


class TestDirectRepeatCost:
    """A tandem array must not cost a second per scan.

    Every seed of one match lies on the same diagonal `b - a`, and both walks
    move the two positions together, so without a memo each of thousands of
    seeds re-walks the whole text. Measured on 6 kb of a 15 bp tandem unit:
    5,974 calls, 35.7M character comparisons, 1.9 seconds -- to return one pair.
    That is the workload BT5 exists for (antibodies, scFv/CAR, (GGGGS)n,
    duplicate 2A peptides), and three repeat rules each paid it, against a 10 s
    end-to-end budget.

    The bounds are generous on purpose: they exist to catch a return to
    re-walking every seed, not to police milliseconds on a shared runner.
    """

    def test_a_tandem_array_stays_affordable(self) -> None:
        import time

        seq = "GGTGGTGGTGGTAGC" * 400
        start = time.perf_counter()
        found = ConstructKmerIndex.of(construct(seq, circular=False), 12).repeat_pairs(12)
        elapsed = time.perf_counter() - start
        assert elapsed < 1.0, f"6 kb tandem array took {elapsed:.2f}s (was 1.9s per seed-walk)"
        assert found, "and it must still find the repeat it was always finding"

    def test_a_period_under_min_len_is_also_memoised(self) -> None:
        """(CAG)n is rejected on every seed because its period is under
        `min_len`, so the extent has to be recorded even when the pair is not --
        otherwise the rejected path re-walks the text thousands of times."""
        import time

        start = time.perf_counter()
        found = ConstructKmerIndex.of(construct("CAG" * 2000, circular=False), 12).repeat_pairs(12)
        elapsed = time.perf_counter() - start
        assert elapsed < 1.0, f"(CAG)x2000 took {elapsed:.2f}s"
        assert found == [], "a 3 bp period is below min_len and is e7's finding, not this one"

    def test_seeds_on_one_diagonal_are_grown_once(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """The invariant itself, counted rather than timed.

        A wall-clock bound loose enough not to flake on a shared runner is too
        loose to catch this case -- 3 kb of tandem inside a 9 kb plasmid ran in
        0.49s even while re-walking every seed. Counting the extensions is
        exact: one maximal match on one diagonal is one walk, however many
        thousands of seeds land inside it.
        """
        seq = dna(3000, seed=163) + "GGTGGTGGTGGTAGC" * 200 + dna(3000, seed=167)
        calls = 0
        original = ConstructKmerIndex._grow

        def counted(self: ConstructKmerIndex, a: int, b: int, min_len: int) -> tuple[int, int]:
            nonlocal calls
            calls += 1
            return original(self, a, b, min_len)

        monkeypatch.setattr(ConstructKmerIndex, "_grow", counted)
        ConstructKmerIndex.of(construct(seq), 15).repeat_pairs(15)
        assert calls < 20, (
            f"{calls} extensions for one tandem region: the diagonal memo is not holding, "
            f"and the cost is back to O(seeds x text)"
        )


class TestInvertedRepeatCost:
    """A plasmid is allowed to be AT-rich, and the scan has to survive it.

    Pairing every occurrence of every k-mer is quadratic, and it degrades
    exactly where a user most wants the answer: an 800 bp alternating-AT run
    produced 464,799 raw pairs and took 26 seconds. The bound is generous on
    purpose -- it exists to catch a return to quadratic, not to police
    milliseconds on a shared runner.
    """

    def test_a_pathological_at_run_stays_affordable(self) -> None:
        import time

        start = time.perf_counter()
        found = ConstructKmerIndex.of(construct("AT" * 2000), 20).inverted_repeats(20, 60)
        elapsed = time.perf_counter() - start
        assert elapsed < 10.0, (
            f"4 kb of alternating AT took {elapsed:.1f}s (was 27s when quadratic)"
        )
        assert len(found) <= 200, "the report must stay bounded however repetitive the input"

    def test_a_long_at_region_inside_a_normal_plasmid_is_cheap(self) -> None:
        seq = dna(3500, seed=151) + "AT" * 500 + dna(3500, seed=157)
        import time

        start = time.perf_counter()
        ConstructKmerIndex.of(construct(seq), 20).inverted_repeats(20, 60)
        assert time.perf_counter() - start < 5.0


class TestProtocolConformance:
    """`ConstructKmerIndex` must actually satisfy `KmerIndex`.

    It did not, until this test existed. `revcomp_pairs` was the name of the
    rich form and returned `list[InvertedRepeat]` where the frozen contract
    promises `Iterator[tuple[Interval, Interval]]`, so a rule reaching through
    `Services.kmer` would have been handed value objects where it expected
    tuples. Nothing caught it because nothing had yet consumed the protocol --
    a protocol with no conforming implementation and no consumer is a promise
    that has never once been checked.
    """

    def hairpin(self, arm: str, loop: int) -> str:
        return dna(200) + arm + dna(loop) + reverse_complement(arm) + dna(200)

    def test_the_static_conformance_assertion_is_where_mypy_can_see_it(self) -> None:
        """The type-level half of this guard lives in `kmers.py`, not here.

        mypy is configured over `packages/engine/src/bt5` only, so a
        `type[KmerIndex]` assertion written in a test file type-checks NOWHERE
        and passes whatever the implementation does -- which is the same shape
        of gap that let the mismatch exist. This test only checks the assertion
        is still present; mypy is what evaluates it.
        """
        source = Path(kmers.__file__).read_text()
        assert "_protocol_conformance: type[KmerIndex] = ConstructKmerIndex" in source

    def test_duplicates_yields_interval_pairs(self) -> None:
        seq = dna(100) + (unit := dna(40, seed=9)) + dna(100) + unit + dna(100)
        for first, second in ConstructKmerIndex.of(construct(seq), 20).duplicates(20):
            assert isinstance(first, Interval)
            assert isinstance(second, Interval)

    def test_revcomp_pairs_yields_interval_pairs_not_value_objects(self) -> None:
        arm = dna(30, seed=11)
        index = ConstructKmerIndex.of(construct(self.hairpin(arm, 10)), 20)
        pairs = list(index.revcomp_pairs(20, 60))
        assert pairs, "the fixture contains a hairpin"
        for first, second in pairs:
            assert isinstance(first, Interval)
            assert isinstance(second, Interval)

    def test_the_narrow_form_agrees_with_the_rich_one(self) -> None:
        """The adapter must not filter or reorder -- a rule and this lane's own
        report have to be looking at the same findings."""
        arm = dna(30, seed=11)
        index = ConstructKmerIndex.of(construct(self.hairpin(arm, 10)), 20)
        rich = index.inverted_repeats(20, 60)
        assert [(r.first, r.second) for r in rich] == list(index.revcomp_pairs(20, 60))

    def test_the_geometry_a_rule_needs_survives_the_narrowing(self) -> None:
        """Stem and loop are recoverable from the pair, which is why the
        contract can stay narrow without costing a caller anything."""
        arm = dna(30, seed=11)
        index = ConstructKmerIndex.of(construct(self.hairpin(arm, 10)), 20)
        repeat = index.inverted_repeats(20, 60)[0]
        first, second = next(iter(index.revcomp_pairs(20, 60)))
        assert first.length == repeat.stem
        assert second.start - first.end == repeat.loop
