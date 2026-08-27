"""Tier-A exact DP, and the differential check against the reference solver.

The reference solver is greedy. This suite exists to prove the DP is worth having:
it finds optima greedy misses, and -- more importantly -- it finds SOLUTIONS
greedy misses, because a greedy choice can paint the sequence into a corner where
no codon for a later residue is placeable. Greedy then reports infeasibility that
is not real, which in a product means telling a user their design is impossible
when it is not.
"""

from __future__ import annotations

import pytest
from bt5.codon.tables import FileTableProvider
from bt5.core.result import InfeasibleConstraints
from bt5.solver.lattice import Automaton, cai_lattice_scorer, optimal_back_translate
from bt5.solver.reference import back_translate, cai_scorer
from bt5.verify import verify_solution

FORBID = ["GAATTC", "GGATCC", "GGTCTC", "AAGCTT"]


@pytest.fixture(scope="module")
def env() -> tuple:
    p = FileTableProvider()
    return p.genetic_code(11), p.usage("sharp_li_1987_ecoli_w")


class TestAutomaton:
    def test_detects_a_pattern(self) -> None:
        a = Automaton(("GAATTC",))
        assert a.consume(0, "AAGAATTCAA")[1] is True
        assert a.consume(0, "AAAAAAAAAA")[1] is False

    def test_detects_a_pattern_split_across_two_consume_calls(self) -> None:
        """This is the property that makes codon-boundary motifs findable: the
        state carries the partial match across the call boundary."""
        a = Automaton(("GAATTC",))
        state, hit = a.consume(0, "AAAGAA")
        assert hit is False
        _, hit2 = a.consume(state, "TTCAAA")
        assert hit2 is True, "a motif spanning two codons must still be caught"

    def test_suffix_matches_are_accepting(self) -> None:
        """Aho-Corasick failure links mean a state accepting a shorter pattern
        must accept even when reached via a longer one's trie path."""
        a = Automaton(("ATTC", "GAATTC"))
        assert a.consume(0, "GAATTC")[1] is True
        assert a.consume(0, "CATTC")[1] is True

    def test_empty_pattern_set_never_accepts(self) -> None:
        a = Automaton(())
        assert a.consume(0, "GAATTCGGATCC")[1] is False


class TestExactness:
    def test_output_is_valid(self, env: tuple) -> None:
        code, u = env
        protein = "MKLIWQRSTVNDEYFPGHACM"
        dna = optimal_back_translate(protein, code, forbidden=FORBID, score=cai_lattice_scorer(u.w))
        verify_solution(protein, dna, table_id=11, forbidden=FORBID, require_initiator=True)

    def test_dp_is_never_worse_than_greedy_on_cai(self, env: tuple) -> None:
        """CAI is a geometric mean, so summing -log(w) makes it a shortest path.
        The DP optimum must therefore dominate the greedy result."""
        code, u = env
        protein = "MKLIWQRSTVNDEYFPGHACMKLIWQ"
        dp = optimal_back_translate(protein, code, forbidden=FORBID, score=cai_lattice_scorer(u.w))
        greedy = back_translate(protein, code, forbidden=FORBID, score=cai_scorer(u.w))
        assert u.cai(dp, code) >= u.cai(greedy, code) - 1e-9

    def test_dp_solves_a_case_where_greedy_reports_false_infeasibility(self, env: tuple) -> None:
        """THE reason Tier A exists.

        Greedy commits to a locally-best codon and can reach a residue with no
        placeable codon left. It then raises InfeasibleConstraints -- telling the
        user the design is impossible when a valid design exists. The DP explores
        the whole lattice and finds it.

        This case was found by random search over tight constraint sets and is
        pinned here as a regression.
        """
        code, _ = env
        protein = "MDNECYA"
        motifs = ["CATC", "AATA", "GACA", "GGTC"]

        with pytest.raises(InfeasibleConstraints):
            back_translate(protein, code, forbidden=motifs, add_stop=False)

        dna = optimal_back_translate(protein, code, forbidden=motifs, add_stop=False)
        assert code.translate(dna) == protein
        for motif in motifs:
            assert motif not in dna

    def test_is_deterministic(self, env: tuple) -> None:
        code, u = env
        kw = {"forbidden": FORBID, "score": cai_lattice_scorer(u.w)}
        a = optimal_back_translate("MKLIWQRSTV", code, **kw)
        b = optimal_back_translate("MKLIWQRSTV", code, **kw)
        assert a == b, "two runs of one protein must not produce two different tubes"


class TestJunctions:
    def test_left_flank_site_is_excluded_by_construction(self, env: tuple) -> None:
        """A site formed half by the vector and half by the first codon."""
        code, _ = env
        left = "GGAAT"  # + TC.. completes EcoRI
        dna = optimal_back_translate("MK", code, forbidden=["GAATTC"], left_flank=left)
        assert "GAATTC" not in left + dna

    def test_right_flank_site_is_excluded_by_construction(self, env: tuple) -> None:
        """The last codon must not complete a motif with the backbone downstream."""
        code, _ = env
        right = "TTC"  # GAA + TTC completes EcoRI
        dna = optimal_back_translate(
            "MK", code, forbidden=["GAATTC"], right_flank=right, add_stop=False
        )
        assert "GAATTC" not in dna + right

    def test_a_forbidden_motif_inside_the_users_backbone_is_reported_not_ignored(
        self, env: tuple
    ) -> None:
        """That is a finding about the vector, not something codon choice can fix,
        so it raises rather than silently producing a construct that violates."""
        code, _ = env
        with pytest.raises(InfeasibleConstraints) as exc:
            optimal_back_translate("MK", code, forbidden=["GAATTC"], left_flank="AAGAATTCAA")
        assert exc.value.certificate.proof == "immutable_region"


class TestInfeasibility:
    def test_unplaceable_residue_raises_a_certificate(self, env: tuple) -> None:
        code, _ = env
        with pytest.raises(InfeasibleConstraints) as exc:
            optimal_back_translate("MW", code, forbidden=["TGG"], add_stop=False)
        cert = exc.value.certificate
        assert cert.proof == "automaton_dead_state"
        assert cert.protein_span == (1, 2)
