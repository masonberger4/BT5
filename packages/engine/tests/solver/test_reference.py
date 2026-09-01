"""The reference solver, and the repeat-breaking requirement."""

from __future__ import annotations

import pytest
from bt5.codon.tables import FileTableProvider
from bt5.core.result import InfeasibleConstraints
from bt5.core.types import reverse_complement
from bt5.solver.reference import (
    IUPAC_CODES,
    MAX_PATTERN_EXPANSION,
    back_translate,
    cai_scorer,
    expand_forbidden,
    expand_iupac,
    longest_repeat,
    repeat_breaking_scorer,
)
from bt5.verify import verify_solution
from hypothesis import given
from hypothesis import strategies as st

FORBID = ["GAATTC", "GGATCC", "GGTCTC"]


@pytest.fixture(scope="module")
def env() -> tuple:
    p = FileTableProvider()
    return p.genetic_code(11), p.usage("sharp_li_1987_ecoli_w")


def test_forbidden_set_is_closed_under_reverse_complement(env: tuple) -> None:
    """Rules declare forward motifs only; the closure happens once, here."""
    assert expand_forbidden(["GGTCTC"]) == ("GAGACC", "GGTCTC")


class TestIupacExpansion:
    """`LatticeTerms.forbidden` is documented IUPAC; the engine has to mean it.

    Before this, a degenerate base reached `Automaton.__init__` and died there as
    a bare `KeyError` on `_BASE_INDEX` -- neither an expansion nor a refusal.
    """

    def test_a_degenerate_base_expands_to_its_acgt_set(self) -> None:
        assert expand_iupac("GGWCC") == ("GGACC", "GGTCC")

    def test_a_pure_acgt_pattern_expands_to_itself(self) -> None:
        """Every motif the catalog declares today is this case, and must be free."""
        assert expand_iupac("GGTCTC") == ("GGTCTC",)

    def test_expansion_precedes_the_reverse_complement_closure(self) -> None:
        """The order is the whole fix, and this is the case that proves it.

        `reverse_complement` is `str.translate`, which leaves an unmapped
        character ALONE. Complementing `RGATC` first yields `GATCR`, whose
        expansion is {GATCA, GATCG} -- a silently WRONG set, forbidding sequences
        nobody asked to forbid and permitting the two that were. Expanding first
        gives the true complements {GATCT, GATCC}.
        """
        got = expand_forbidden(["RGATC"])
        assert got == ("AGATC", "GATCC", "GATCT", "GGATC")
        assert "GATCA" not in got
        assert "GATCG" not in got

    def test_reverse_complement_alone_would_have_produced_the_wrong_set(self) -> None:
        """Pins the bug this replaced, so nobody reintroduces the old order."""
        assert reverse_complement("RGATC") == "GATCR"  # R passed through uncomplemented
        assert set(expand_iupac("GATCR")) != {"GATCT", "GATCC"}

    def test_a_character_that_is_not_iupac_is_refused_by_name(self) -> None:
        with pytest.raises(ValueError, match="not an IUPAC") as exc:
            expand_forbidden(["GGXCC"])
        assert "GGXCC" in str(exc.value)
        assert "'X'" in str(exc.value)

    def test_an_empty_pattern_is_refused(self) -> None:
        """An empty motif matches everywhere: it would forbid every codon."""
        with pytest.raises(ValueError, match="empty"):
            expand_forbidden([""])

    def test_an_unbounded_expansion_is_refused_rather_than_built(self) -> None:
        """A silent exponential reads as a hung solver. Refuse, and say why."""
        with pytest.raises(ValueError, match="over the cap") as exc:
            expand_forbidden(["N" * 8])
        assert "NNNNNNNN" in str(exc.value)
        assert "65536" in str(exc.value)
        assert str(MAX_PATTERN_EXPANSION) in str(exc.value)

    def test_a_pattern_exactly_at_the_cap_is_accepted(self) -> None:
        """The cap is a ceiling, not a strict bound -- 4**5 == 1024 must pass."""
        assert len(expand_iupac("N" * 5)) == MAX_PATTERN_EXPANSION

    def test_lowercase_is_accepted(self) -> None:
        assert expand_iupac("ggwcc") == ("GGACC", "GGTCC")

    def test_the_solver_and_the_oracle_agree_on_expansion(self) -> None:
        """The solver keeps its OWN table rather than importing the validator's.

        `tests/data_integrity/test_oracle_independence.py` keeps `verify.py` off
        every lane's code path; importing its expander here would defeat that
        pointing the other way, because one transposed row would then be
        invisible -- the design and its check would forbid the same wrong set.
        Two tables plus this test is the trade: divergence fails loudly here.
        """
        from bt5.verify import IUPAC_EXPANSION
        from bt5.verify import expand_iupac as oracle_expand_iupac

        assert IUPAC_CODES == IUPAC_EXPANSION
        vector = [
            "GGTCTC",  # BsaI, pure ACGT
            "GGWCC",  # AvaII
            "RGATC",
            "MAGGTRAGT",  # splice donor consensus
            "GCCRCCATGG",  # Kozak-like
            "BDHVN",  # every three- and four-fold code in one pattern
            "acgt",  # case
            "AYCGRTSWKM",
        ]
        for pattern in vector:
            assert list(expand_iupac(pattern)) == oracle_expand_iupac(pattern), pattern

    def test_a_forbidden_iupac_motif_is_avoided_on_both_strands(self, env: tuple) -> None:
        """End to end: the greedy solver honours a degenerate motif."""
        code, u = env
        protein = "MKLIWQRSTVNDEYFP"
        dna = back_translate(protein, code, forbidden=["GGWCC"], score=cai_scorer(u.w))
        assert code.translate(dna)[:-1] == protein
        for motif in ("GGACC", "GGTCC"):
            assert motif not in dna
            assert reverse_complement(motif) not in dna


#: Lengths cap at 5 so the widest pattern (NNNNN) is exactly MAX_PATTERN_EXPANSION,
#: keeping the property about expansion CORRECTNESS rather than about the cap.
iupac_patterns = st.text(alphabet=sorted(IUPAC_CODES), min_size=1, max_size=5)


@given(pattern=iupac_patterns)
def test_every_expansion_member_is_acgt_and_matches_position_by_position(
    pattern: str,
) -> None:
    members = expand_iupac(pattern)
    for member in members:
        assert set(member) <= set("ACGT")
        assert len(member) == len(pattern)
        for base, code in zip(member, pattern, strict=True):
            assert base in IUPAC_CODES[code]
    # Exactly the cartesian product: nothing dropped, nothing duplicated.
    expected = 1
    for code in pattern:
        expected *= len(IUPAC_CODES[code])
    assert len(set(members)) == len(members) == expected


@given(pattern=iupac_patterns)
def test_the_forbidden_set_is_closed_under_reverse_complement(pattern: str) -> None:
    """Closed over the EXPANDED set, which is the property rules rely on: a rule
    lists one forward motif and gets both strands of every sequence it denotes."""
    closure = expand_forbidden([pattern])
    assert set(expand_iupac(pattern)) <= set(closure)
    for motif in closure:
        assert set(motif) <= set("ACGT")
        assert reverse_complement(motif) in closure
    assert closure == tuple(sorted(closure))


def test_output_translates_back_and_avoids_forbidden_motifs(env: tuple) -> None:
    code, u = env
    protein = "MKLIWQRSTVNDEYFP"
    dna = back_translate(protein, code, forbidden=FORBID, score=cai_scorer(u.w))
    verify_solution(protein, dna, table_id=11, forbidden=FORBID, require_initiator=True)


def test_is_deterministic(env: tuple) -> None:
    """Two runs of one protein must never produce two different tubes."""
    code, u = env
    kw = {"forbidden": FORBID, "score": cai_scorer(u.w)}
    assert back_translate("MKLIWQ", code, **kw) == back_translate("MKLIWQ", code, **kw)


def test_junction_spanning_site_is_prevented(env: tuple) -> None:
    """A site formed half by the vector and half by the first codon must be
    caught during design, not discovered later by the validator."""
    code, _ = env
    left = "GGAAT"  # + TC.. would complete EcoRI GAATTC
    dna = back_translate("MK", code, forbidden=["GAATTC"], left_flank=left)
    assert "GAATTC" not in left + dna


def test_infeasible_constraints_raise_with_a_certificate(env: tuple) -> None:
    """ "No solution" is a useless product; the minimal conflicting set is the
    deliverable."""
    code, _ = env
    # Forbid every codon for Trp (only TGG) and it becomes unplaceable.
    with pytest.raises(InfeasibleConstraints) as exc:
        back_translate("MW", code, forbidden=["TGG"])
    cert = exc.value.certificate
    assert cert.proof == "empty_mutation_space"
    assert cert.protein_span == (1, 2)


class TestRepeatBreaking:
    """Max-CAI is actively harmful for the proteins people actually express."""

    PROTEIN = "MKLIWQ" + "GGGGS" * 3 + "HHHHHH"  # scFv-style linker + His tag

    def test_max_cai_produces_a_dangerous_repeat(self, env: tuple) -> None:
        code, u = env
        dna = back_translate(self.PROTEIN, code, forbidden=FORBID, score=cai_scorer(u.w))
        rep = longest_repeat(dna, min_len=8)
        assert rep is not None
        # Above the 20 bp vendor reject threshold AND the 23-27 bp RecBCD MEPS.
        assert len(rep[0]) >= 20, "this test documents the failure mode; it should be severe"

    def test_repeat_breaking_cuts_the_repeat_well_below_threshold(self, env: tuple) -> None:
        code, u = env
        dna = back_translate(
            self.PROTEIN, code, forbidden=FORBID, score=repeat_breaking_scorer(u.w)
        )
        verify_solution(self.PROTEIN, dna, table_id=11, forbidden=FORBID, require_initiator=True)
        rep = longest_repeat(dna, min_len=8)
        longest = len(rep[0]) if rep else 0
        assert longest < 15, f"longest repeat {longest} bp must stay well under the 20 bp floor"

    def test_repeat_breaking_beats_max_cai_on_repeat_length(self, env: tuple) -> None:
        code, u = env
        greedy = back_translate(self.PROTEIN, code, forbidden=FORBID, score=cai_scorer(u.w))
        broken = back_translate(
            self.PROTEIN, code, forbidden=FORBID, score=repeat_breaking_scorer(u.w)
        )
        g = longest_repeat(greedy, min_len=8)
        b = longest_repeat(broken, min_len=8)
        assert len(b[0]) < len(g[0]) if (g and b) else True

    def test_repeat_breaking_still_translates_correctly(self, env: tuple) -> None:
        """Breaking repeats must never change the protein."""
        code, u = env
        dna = back_translate(
            self.PROTEIN, code, forbidden=FORBID, score=repeat_breaking_scorer(u.w)
        )
        assert code.translate(dna)[:-1] == self.PROTEIN

    def test_repeat_breaking_costs_some_cai_and_that_is_correct(self, env: tuple) -> None:
        """Repeat-breaking OVERRIDES codon optimality rather than competing with
        it, so a CAI drop here is the intended trade, not a regression."""
        code, u = env
        broken = back_translate(
            self.PROTEIN, code, forbidden=FORBID, score=repeat_breaking_scorer(u.w)
        )
        cai = u.cai(broken, code)
        assert 0.6 < cai < 1.0, f"CAI {cai:.3f} should stay healthy while repeats are broken"


def test_longest_repeat_finds_nothing_in_a_unique_sequence() -> None:
    assert longest_repeat("ACGTACGGTTACCAGGATTCA", min_len=8) is None


class TestLongestRepeatMeasuresRatherThanBounds:
    """`longest_repeat` is the reference measurement the repeat rules rest on.

    It used to cap its search at `min(n // 2, 40)` and return the cap. Both
    halves under-reported on exactly the sequences the rules exist for, because
    one-codon-per-amino-acid back-translation of a repetitive protein produces a
    TANDEM array, where the copies overlap and neither bound holds.
    """

    #: What max-CAI does to a (G4S)3 linker: one codon per residue, three times.
    LINKER = "GGTGGTGGTGGTTCT" * 3

    def brute(self, seq: str, min_len: int = 8) -> tuple[str, int, int] | None:
        """Deliberately the slowest possible correct answer, for differencing."""
        for size in range(len(seq) - 1, min_len - 1, -1):
            seen: dict[str, int] = {}
            for i in range(len(seq) - size + 1):
                k = seq[i : i + size]
                if k in seen:
                    return (k, seen[k], i)
                seen[k] = i
        return None

    def test_overlapping_copies_are_found(self) -> None:
        """45 bp of period-15 tandem holds a 30 bp repeat; n // 2 says 22."""
        found = longest_repeat(self.LINKER)
        assert found is not None
        assert len(found[0]) == 30
        assert (found[1], found[2]) == (0, 15), "the copies overlap, which is the point"

    def test_a_long_tandem_is_not_truncated_at_forty(self) -> None:
        found = longest_repeat("GGTGGTGGTGGTTCT" * 8)
        assert found is not None
        assert len(found[0]) == 105

    def test_it_agrees_with_brute_force_on_the_linker(self) -> None:
        exact = self.brute(self.LINKER)
        found = longest_repeat(self.LINKER)
        assert exact is not None
        assert found is not None
        assert len(found[0]) == len(exact[0])

    @pytest.mark.parametrize("seed", [1, 2, 3, 5, 8, 13])
    def test_it_agrees_with_brute_force_on_random_sequence(self, seed: int) -> None:
        """Differential against the slow implementation, planted repeats and all.

        Seeded explicitly: an unseeded case that failed once could not be
        reproduced, which is the whole reason the global RNG is banned in src/.
        """
        import numpy as np

        rng = np.random.default_rng(seed)
        for _ in range(40):
            n = int(rng.integers(10, 160))
            seq = "".join("ACGT"[i] for i in rng.integers(0, 4, n))
            if rng.random() < 0.5:
                unit = "".join("ACGT"[i] for i in rng.integers(0, 4, int(rng.integers(8, 25))))
                gap = "".join("ACGT"[i] for i in rng.integers(0, 4, int(rng.integers(0, 30))))
                seq = seq + unit + gap + unit
            found, exact = longest_repeat(seq), self.brute(seq)
            assert (found is None) == (exact is None)
            if found is not None and exact is not None:
                assert len(found[0]) == len(exact[0]), f"disagreed on {seq!r}"

    def test_both_reported_copies_really_carry_the_repeat(self) -> None:
        """Checked against the sequence, not against the search that found it."""
        found = longest_repeat(self.LINKER)
        assert found is not None
        kmer, first, second = found
        assert self.LINKER[first : first + len(kmer)] == kmer
        assert self.LINKER[second : second + len(kmer)] == kmer
        assert first != second

    def test_min_len_is_still_a_floor(self) -> None:
        assert longest_repeat("ACGTACGGTTACCAGGATTCA", min_len=8) is None
        assert longest_repeat("AAAAAAAA", min_len=20) is None

    def test_a_sequence_too_short_to_repeat_returns_nothing(self) -> None:
        assert longest_repeat("ACGT", min_len=8) is None
        assert longest_repeat("", min_len=1) is None
