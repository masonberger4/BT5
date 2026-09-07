"""D1: forbidden restriction / Type IIS sites -- and what a HARD_LATTICE rule owes.

D1 is the catalog's reference pattern rule and was, until this file, the only
rule with no paired test. So this file is also documentation of a lattice rule's
contract, and it is organised around the three cases such a rule exists to make
unmissable and that a per-codon or single-strand check cannot see at all:

- a site formed ACROSS the CDS/backbone junction, which no scan of the designed
  CDS alone can find;
- a site SPANNING THE ORIGIN of a circular construct, which a linear scan splits
  into two clean halves;
- a non-palindromic site present only as its REVERSE COMPLEMENT, which the rule
  never scans for itself -- it lists forward motifs and the shared motif finder
  closes the set.

The last one is the point of CLAUDE.md 3.4. D1 lists `GGTCTC` and gets `GAGACC`
for free; a rule that scanned the reverse strand itself would report both and
double-count.
"""

from __future__ import annotations

from bt5.core.registry import discover
from bt5.core.spec import Enforcement
from bt5.core.types import (
    Construct,
    Interval,
    Segment,
    SegmentKind,
    Topology,
)
from bt5.rules.catalog.d1_restriction_sites import (
    SIX_CUTTERS,
    TYPE_IIS,
    RestrictionSites,
)
from conftest import construct, context, wrapping_construct

discover()

#: No six-cutter and no Type IIS site, forward or reverse complement.
CLEAN_CDS = "ATGAAACCCGTGACAAGTTAA"


def evaluate(rule: RestrictionSites, c: Construct):
    return rule.evaluate(c, context(), None)  # type: ignore[arg-type]


class TestConfiguration:
    def test_six_cutters_are_on_by_default_and_type_iis_is_not(self) -> None:
        """Golden Gate is not the assumed assembly method, so domestication of
        the Type IIS enzymes is opt-in rather than imposed on every design."""
        enabled = set(RestrictionSites().enzymes)
        assert enabled == set(SIX_CUTTERS)
        assert not enabled & set(TYPE_IIS)

    def test_type_iis_is_opt_in_by_name(self) -> None:
        rule = RestrictionSites(enzymes=("BsaI", "EcoRI"))
        assert rule.enzymes == {"BsaI": "GGTCTC", "EcoRI": "GAATTC"}

    def test_an_unknown_enzyme_name_is_dropped_rather_than_raising(self) -> None:
        """The name comes from `param_schema`, i.e. from a UI control, so an
        unknown one is a stale preset rather than a programming error."""
        assert RestrictionSites(enzymes=("EcoRI", "NotAnEnzyme")).enzymes == {"EcoRI": "GAATTC"}

    def test_no_arg_constructible(self) -> None:
        """`presets._unscored_enforcement` probes rules by calling `spec()`, and
        on TypeError falls back to the ClassVar floor. Pinned per rule (#82)."""
        assert RestrictionSites().enzymes


class TestLatticeContract:
    def test_declares_the_configured_sites_as_forbidden_motifs(self) -> None:
        terms = RestrictionSites(enzymes=("EcoRI", "BsaI")).lattice_terms(None)  # type: ignore[arg-type]
        assert terms is not None
        assert set(terms.forbidden) == {"GAATTC", "GGTCTC"}

    def test_lists_forward_motifs_only(self) -> None:
        """`GGTCTC` is listed; `GAGACC` is not. The solver closes the set under
        reverse complement at automaton-construction time, so listing both would
        double-count -- and D1's own docstring says it does not close it here."""
        forbidden = set(RestrictionSites(enzymes=("BsaI",)).lattice_terms(None).forbidden)  # type: ignore[arg-type,union-attr]
        assert forbidden == {"GGTCTC"}
        assert "GAGACC" not in forbidden

    def test_a_lattice_rule_carries_no_objective_or_steering_weight(self) -> None:
        """HARD_LATTICE is unreachable by construction: an objective weight would
        be the penalty-as-guarantee that CLAUDE.md 3.5 bans, and a steering term
        would nudge the DP toward something it already cannot reach."""
        assert RestrictionSites.enforcement is Enforcement.HARD_LATTICE
        assert RestrictionSites.default_weight == 0.0
        assert RestrictionSites.steering_weight == 0.0


class TestEvaluate:
    def test_a_clean_construct_passes_with_no_breaches(self) -> None:
        ev = evaluate(RestrictionSites(), construct(CLEAN_CDS))
        assert ev.passes
        assert ev.breaches == ()
        assert ev.raw_score == 0.0

    def test_a_site_in_the_cds_is_a_breach_the_solver_may_fix(self) -> None:
        c = construct("ATG" + "GAATTC" + "AAATAA")
        (breach,) = evaluate(RestrictionSites(), c).breaches
        assert breach.interval == Interval(3, 9)
        assert breach.fixable_by_codon_choice
        assert breach.detail["enzyme"] == "EcoRI"

    def test_the_message_names_the_enzyme_and_the_offending_motif(self) -> None:
        """`Breach.message` must name the exact offending substring: it is what
        the report prints and the only place a human can check the finding."""
        c = construct("ATG" + "GAATTC" + "AAATAA")
        (breach,) = evaluate(RestrictionSites(), c).breaches
        assert "EcoRI" in breach.message
        assert "GAATTC" in breach.message

    def test_a_site_the_backbone_already_carries_is_reported_but_not_chased(self) -> None:
        """Real, and worth a line in the report -- but no codon can remove it.
        Saying so is what stops the solver exhausting the mutation space over a
        base it is not allowed to touch."""
        c = construct(CLEAN_CDS, backbone="GGG" + "GAATTC" + "GGG")
        (breach,) = evaluate(RestrictionSites(), c).breaches
        assert not breach.fixable_by_codon_choice
        assert breach.detail["enzyme"] == "EcoRI"

    def test_a_site_formed_across_the_cds_backbone_junction_is_caught(self) -> None:
        """Neither half is a site. `GAAT` ends the CDS and `TC` opens the
        backbone, and only the ASSEMBLED construct shows the EcoRI site -- which
        is why a rule takes a Construct and never a bare string (CLAUDE.md 3.3).
        It is fixable, because the CDS half of it is editable."""
        cds = "ATGAAACCCGAAT"
        c = construct(cds, backbone="TCGGGCCC")
        (breach,) = evaluate(RestrictionSites(), c).breaches
        assert breach.interval == Interval(9, 15)
        assert breach.interval.start < len(cds) < breach.interval.end, "it straddles the junction"
        assert breach.fixable_by_codon_choice

    def test_a_site_spanning_the_origin_is_caught_and_reported_as_wrapping(self) -> None:
        """The case a linear scan cannot see: four bases close the sequence and
        two open it. Reported at its true start with an interval whose end runs
        past the construct length."""
        c = wrapping_construct(prefix="GGGCCCGGG", cds_head="TCAAAG", cds_tail="AAACGAAT")
        (breach,) = evaluate(RestrictionSites(), c).breaches
        assert breach.interval.start == 19
        assert breach.interval.wraps(c.length)
        assert breach.detail["enzyme"] == "EcoRI"

    def test_the_same_bases_laid_out_linearly_are_clean(self) -> None:
        """The control for the test above: circularity is what creates the site,
        so the same bases in a linear construct must report nothing. Without this
        the wrapping test would still pass if the finder simply matched twice."""
        seq = wrapping_construct(
            prefix="GGGCCCGGG", cds_head="TCAAAG", cds_tail="AAACGAAT"
        ).sequence
        linear = Construct(
            seq,
            Topology.LINEAR,
            (Segment(Interval(0, len(seq)), SegmentKind.DESIGNABLE_CDS, "cds"),),
        )
        assert evaluate(RestrictionSites(), linear).breaches == ()

    def test_a_non_palindromic_site_is_caught_on_the_reverse_strand(self) -> None:
        """`GAGACC` is the reverse complement of the listed `GGTCTC`. The rule
        never scans the reverse strand; the shared finder closes the set, and the
        breach is still attributed to the LISTED motif so the enzyme is named."""
        c = construct("ATG" + "GAGACC" + "AAATAA")
        (breach,) = evaluate(RestrictionSites(enzymes=("BsaI",)), c).breaches
        assert breach.interval == Interval(3, 9)
        assert breach.detail["enzyme"] == "BsaI"
        assert breach.detail["motif"] == "GGTCTC", "attributed to the forward motif that was listed"

    def test_the_palindromic_six_cutters_are_found_once_not_twice(self) -> None:
        """`GAATTC` is its own reverse complement. Closing the set must not turn
        one physical site into two findings."""
        c = construct("ATG" + "GAATTC" + "AAATAA")
        assert len(evaluate(RestrictionSites(), c).breaches) == 1

    def test_raw_score_is_the_number_of_sites(self) -> None:
        c = construct("ATG" + "GAATTC" + "AAA" + "GGATCC" + "TAA")
        ev = evaluate(RestrictionSites(), c)
        assert not ev.passes
        assert ev.raw_score == 2.0
        assert {b.detail["enzyme"] for b in ev.breaches} == {"EcoRI", "BamHI"}

    def test_a_disabled_enzyme_is_not_reported(self) -> None:
        """BsaI sites are legal unless the user opted into Type IIS domestication."""
        c = construct("ATG" + "GGTCTC" + "AAATAA")
        assert evaluate(RestrictionSites(), c).breaches == ()
        assert evaluate(RestrictionSites(enzymes=("BsaI",)), c).breaches
