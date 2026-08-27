"""The reference solver, and the repeat-breaking requirement."""

from __future__ import annotations

import pytest
from bt5.codon.tables import FileTableProvider
from bt5.core.result import InfeasibleConstraints
from bt5.solver.reference import (
    back_translate,
    cai_scorer,
    expand_forbidden,
    longest_repeat,
    repeat_breaking_scorer,
)
from bt5.verify import verify_solution

FORBID = ["GAATTC", "GGATCC", "GGTCTC"]


@pytest.fixture(scope="module")
def env() -> tuple:
    p = FileTableProvider()
    return p.genetic_code(11), p.usage("sharp_li_1987_ecoli_w")


def test_forbidden_set_is_closed_under_reverse_complement(env: tuple) -> None:
    """Rules declare forward motifs only; the closure happens once, here."""
    assert expand_forbidden(["GGTCTC"]) == ("GAGACC", "GGTCTC")


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
