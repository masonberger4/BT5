"""X1: what the protein decides before the solver is allowed an opinion.

The load-bearing test is `test_bounds_match_exhaustive_enumeration`. Everything
else in this file checks behaviour; that one checks the arithmetic is EXACT, by
enumerating every codon assignment of a short protein and comparing. If the
envelope were merely an estimate, the honesty contract -- outside the range is
PROVEN unreachable -- would be false, and every message built on it misleading.
"""

from __future__ import annotations

import itertools

import pytest
from bt5.cassette.envelope import DEFAULT_WINDOW_BP, envelope
from bt5.codon.tables import FileTableProvider

CODE = FileTableProvider().genetic_code(11)


def gc_percent(dna: str) -> float:
    return 100.0 * sum(b in "GC" for b in dna) / len(dna)


class TestExactness:
    def test_bounds_match_exhaustive_enumeration(self) -> None:
        """Brute force over every synonymous assignment of a short protein.

        Leu and Arg have six codons each and Gly four, so this is 6*6*4*2 = 288
        sequences -- small enough to enumerate and varied enough that an
        off-by-one in the partial-codon arithmetic would show.
        """
        protein = "LRGC"
        families = CODE.families()
        window = 3 * len(protein)
        env = envelope(protein, CODE, window=window)

        every = ["".join(combo) for combo in itertools.product(*(families[aa] for aa in protein))]
        assert len(every) == 288
        observed = [gc_percent(dna) for dna in every]

        assert env.lowest == pytest.approx(min(observed))
        assert env.highest == pytest.approx(max(observed))

    def test_partial_codon_windows_are_exact_too(self) -> None:
        """A window boundary mid-codon is where an estimate would go wrong."""
        protein = "WWWW"  # single codon TGG, so every window has one true answer
        env = envelope(protein, CODE, window=5)
        for w in env.windows:
            true = gc_percent(("TGG" * 4)[w.start : w.end])
            assert w.lowest == pytest.approx(true)
            assert w.highest == pytest.approx(true)


class TestWhatTheProteinForces:
    def test_the_three_reference_envelopes(self) -> None:
        """Cross-checked against an independent measurement of the same thing."""
        assert envelope("HHHHHH", CODE, window=18).lowest == pytest.approx(33.3, abs=0.1)
        assert envelope("HHHHHH", CODE, window=18).highest == pytest.approx(66.7, abs=0.1)

        linker = envelope("GGGGS" * 3, CODE, window=45)
        assert linker.lowest == pytest.approx(60.0, abs=0.1)
        assert linker.highest == pytest.approx(93.3, abs=0.1)

    def test_a_glycine_serine_linker_cannot_be_brought_below_60_percent(self) -> None:
        """The headline: the textbook repeat motif's hard constraint is GC."""
        assert envelope("GGGGS" * 3, CODE, window=45).lowest > 59.0

    def test_single_codon_residues_leave_no_choice_at_all(self) -> None:
        for aa, expected in (("W", 66.7), ("M", 33.3)):
            env = envelope(aa * 10, CODE, window=30)
            assert env.lowest == pytest.approx(expected, abs=0.1)
            assert env.highest == pytest.approx(expected, abs=0.1)
            assert env.windows[0].span == pytest.approx(0.0)

    def test_at_rich_residues_cap_the_ceiling(self) -> None:
        env = envelope("K" * 10, CODE, window=30)
        assert env.highest == pytest.approx(33.3, abs=0.1)


class TestTheHonestyContract:
    def test_a_target_outside_the_range_is_proven_unreachable(self) -> None:
        env = envelope("GGGGS" * 6, CODE)
        assert env.unreachable(20.0, 40.0), "a 60% floor cannot reach a 40% ceiling"

    def test_e2s_own_default_band_is_unsatisfiable_over_a_linker(self) -> None:
        """The case that motivates the feature.

        E2 ships a 40-60% band and is HARD_REPAIR, so the validator refuses to
        emit outside it. A construct carrying a GGGGS linker cannot satisfy it,
        and today that is discovered only at the very end of the pipeline.
        """
        flank = "AVLITKEQ" * 10
        protein = "M" + flank + "GGGGS" * 6 + flank
        bad = envelope(protein, CODE).unreachable(40.0, 60.0)
        assert bad, "the linker must be provably out of E2's band"
        worst = max(bad, key=lambda w: w.lowest)
        assert worst.lowest > 60.0
        assert "G" in "".join(c.residues for c in worst.holds_low)

    def test_an_empty_result_is_not_a_feasibility_guarantee(self) -> None:
        """Per-window exact does not compose into jointly achievable.

        Each window here admits the band on its own, so nothing is PROVEN
        unreachable -- which is all `unreachable()` claims. The name and the
        docstring carry that limit; this test pins that the API never implies
        more than it can support.
        """
        env = envelope("AVLITKEQ" * 20, CODE)
        assert env.unreachable(25.0, 65.0) == ()
        assert "PROVEN" in (env.unreachable.__doc__ or "")

    def test_a_wide_band_proves_nothing_even_over_a_hard_linker(self) -> None:
        """A 25-65% band tolerates the linker: its floor is 60-62%, inside it.

        Recorded because the obvious illustration is wrong -- the linker is a
        real constraint and NOT a violation of every band, and a message that
        said otherwise would be crying wolf.
        """
        assert envelope("GGGGS" * 6, CODE).unreachable(25.0, 65.0) == ()


class TestCulprits:
    def test_the_bound_names_residues_in_protein_coordinates(self) -> None:
        """Protein coordinates because every available action is a protein edit."""
        env = envelope("GGGGS" * 6, CODE)
        held = max(env.windows, key=lambda w: w.lowest).holds_low
        assert held
        for c in held:
            assert 0 <= c.start < c.end <= len(env.protein)
            assert c.residues == env.protein[c.start : c.end]

    def test_the_high_bound_names_the_residues_capping_it(self) -> None:
        env = envelope("K" * 10 + "GGGG", CODE, window=30)
        capped = min(env.windows, key=lambda w: w.highest)
        assert "K" in "".join(c.residues for c in capped.holds_high)


class TestForcedRepeats:
    def test_repeated_single_codon_runs_force_an_identical_stretch(self) -> None:
        assert envelope("WWAAAAWW", CODE, window=12).forced_repeat_bp == 6

    def test_a_gly_ser_linker_forces_no_repeat_at_all(self) -> None:
        """The finding, not a shortcoming: repeats are recodeable, GC is not.

        Every GGGGS copy can be encoded differently -- Gly has four codons and
        Ser six -- while the composition floor is immovable.
        """
        assert envelope("GGGGS" * 6, CODE).forced_repeat_bp == 0
        assert envelope("GGGGS" * 6, CODE).lowest > 59.0


class TestTheCodeIsRequired:
    def test_table_4_gives_tryptophan_a_choice_it_does_not_have_elsewhere(self) -> None:
        """Table 4 reads TGA as Trp rather than stop, so W stops being fixed.

        Under table 11 a poly-W stretch is pinned at 66.7% GC with no codon
        choice at all; under table 4 it can reach 33.3%. Same protein, different
        envelope -- which is why `code` is required and never defaulted.
        """
        table4 = FileTableProvider().genetic_code(4)
        assert envelope("WWWW", CODE, window=12).windows[0].span == pytest.approx(0.0)
        assert envelope("WWWW", table4, window=12).lowest == pytest.approx(33.3, abs=0.1)

    def test_table_2_raises_the_arginine_floor(self) -> None:
        """Table 2 makes AGA/AGG stops, removing Arg's two lowest-GC codons."""
        table2 = FileTableProvider().genetic_code(2)
        assert envelope("RRRR", CODE, window=12).lowest == pytest.approx(33.3, abs=0.1)
        assert envelope("RRRR", table2, window=12).lowest == pytest.approx(66.7, abs=0.1)

    def test_table_12_is_protein_critical_but_gc_neutral(self) -> None:
        """A negative result worth pinning.

        Table 12 moves CTG from Leu to Ser, which silently changes the PROTEIN
        and is the reason the table is never defaulted anywhere in BT5 -- but it
        does not move either family's GC bounds (Leu keeps a 0-GC and a 2-GC
        codon; Ser keeps a 1 and a 2), so the envelope is unchanged. The
        envelope is sensitive to some table differences and not others, and
        conflating the two would be a wrong argument for a right rule.
        """
        table12 = FileTableProvider().genetic_code(12)
        for aa in ("L", "S"):
            here, there = envelope(aa * 4, CODE, window=12), envelope(aa * 4, table12, window=12)
            assert here.lowest == pytest.approx(there.lowest)
            assert here.highest == pytest.approx(there.highest)

    def test_a_residue_with_no_codon_is_refused(self) -> None:
        with pytest.raises(ValueError, match="cannot be encoded under it"):
            envelope("MAGIC*", CODE)


class TestValidation:
    def test_an_empty_protein_is_refused(self) -> None:
        with pytest.raises(ValueError, match="no envelope"):
            envelope("", CODE)

    def test_a_non_positive_window_is_refused(self) -> None:
        with pytest.raises(ValueError, match="window must be positive"):
            envelope("MAAA", CODE, window=0)

    def test_a_protein_shorter_than_the_window_still_yields_one(self) -> None:
        env = envelope("MA", CODE, window=DEFAULT_WINDOW_BP)
        assert len(env.windows) == 1
        assert env.windows[0].end == 6
