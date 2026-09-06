"""D5: bacterial cryptic transcription -- the motif half, the geometry half, and
the seam between them.

The rule carries two enforcement mechanisms on purpose (see its module
docstring), so this file tests them as two things: patterns the Tier-A automaton
makes unreachable, and promoter GEOMETRY that only Tier-B repair can chase. The
seam is where the interesting failures live, so the boundary cases -- a spacer
one base outside the range, an AT-tract one base too far away -- are pinned
explicitly rather than left to a happy path.
"""

from __future__ import annotations

import re

import pytest
from bt5.core.context import Modality
from bt5.core.registry import discover
from bt5.core.spec import Breach, Enforcement, Evaluation, RepairPolicy
from bt5.core.types import Construct, reverse_complement
from bt5.rules.catalog.d5_cryptic_transcription import (
    AT_TRACT_GAP,
    SPACER_MAX,
    SPACER_MIN,
    CrypticTranscription,
    _within_one_mismatch,
)
from conftest import construct, context, slot, wrapping_construct

discover()

#: GC-rich filler: contains no -10-like or -35-like hexamer, so a test sequence
#: built from it fires only on what the test deliberately puts in.
PAD = "GGCCGG"
SPACER = "CGCGCGCGCGCGCGCGC"  # 17 bp, inside the 15-19 range


def run(rule: CrypticTranscription, c: Construct, ctx=None) -> Evaluation:
    return rule.evaluate(c, ctx or context(), None)  # type: ignore[arg-type]


def of_kind(ev: Evaluation, kind: str) -> list[Breach]:
    return [b for b in ev.breaches if b.detail.get("kind") == kind]


class TestVariantExpansion:
    def test_one_mismatch_of_a_hexamer_is_nineteen_strings(self) -> None:
        """The exact match plus six positions times three substitutions. Held as
        a precomputed tuple so the hot path is `str.find` rather than a
        per-position character comparison -- this runs once per repair
        candidate."""
        assert len(_within_one_mismatch("TTGACA", 1)) == 19

    def test_budget_zero_is_the_consensus_alone(self) -> None:
        assert _within_one_mismatch("TATAAT", 0) == ("TATAAT",)

    def test_budget_two_is_still_small_enough_to_enumerate(self) -> None:
        assert len(_within_one_mismatch("TTGACA", 2)) == 154


class TestContract:
    def test_the_classvar_is_lattice_and_the_probe_is_repair(self) -> None:
        """Not a contradiction -- two independent switches. `solver/catalog.py`
        reads the ClassVar when it collects forbidden patterns for the automaton
        and reads `enforcement_for` when it routes breaches to Tier-B. Returning
        HARD_LATTICE from the probe would leave the geometry breaches matching
        neither arm of the router, so they would be found and then dropped."""
        assert CrypticTranscription.enforcement is Enforcement.HARD_LATTICE
        assert CrypticTranscription().enforcement_for(slot()) is Enforcement.HARD_REPAIR

    def test_a_hard_rule_carries_no_objective_or_steering_weight(self) -> None:
        assert CrypticTranscription.default_weight == 0.0
        assert CrypticTranscription.steering_weight == 0.0

    def test_repair_is_fixed_point(self) -> None:
        """Mutating one element of a promoter creates -10-like hexamers nearby,
        and a new hexamer 5 bp downstream of an existing AT-tract is a NEW
        promoter a single pass would never re-scan for (CLAUDE.md 3.6)."""
        assert CrypticTranscription.repair is RepairPolicy.FIXED_POINT

    def test_no_arg_constructible_and_probeable(self) -> None:
        rule = CrypticTranscription()
        assert callable(rule.gate)
        assert callable(rule.enforcement_for)

    def test_declares_its_brief_row(self) -> None:
        assert CrypticTranscription.brief_ref == "2.D5"

    def test_declares_the_rules_it_opposes(self) -> None:
        """D5 removes AT-rich hexamers, raising local GC against f5's ceiling;
        d8's CpG fixes lower GC and walk back toward what D5 forbids."""
        assert set(CrypticTranscription.conflicts_with) == {"d8_cpg_depletion", "f5_at_window"}
        assert CrypticTranscription.id not in CrypticTranscription.conflicts_with

    @pytest.mark.parametrize("modality", list(Modality))
    def test_gates_on_every_modality_including_ivt_mrna(self, modality: Modality) -> None:
        """The trap. d3, d4 and d6 all gate OFF IVT_MRNA because they read the
        TRANSCRIPT. D5 reads the template plasmid during E. coli propagation,
        which an IVT mRNA construct still goes through -- brief.md:110 says
        "applies to ALL constructs regardless of final host"."""
        assert CrypticTranscription().gate(slot(modality=modality))


class TestLatticeTerms:
    def test_declares_the_three_pattern_shaped_parts(self) -> None:
        """The extended -10 is declared pre-expanded, so six strings cover three
        sub-rules."""
        forbidden = set(CrypticTranscription().lattice_terms(None).forbidden)
        assert forbidden == {
            "ATTATA",
            "TATACT",
            "TGATATAAT",
            "TGCTATAAT",
            "TGGTATAAT",
            "TGTTATAAT",
        }

    def test_every_declared_motif_is_pure_acgt(self) -> None:
        """The regression this rule was the first in the catalog to hit.

        `LatticeTerms.forbidden` is documented IUPAC and the SOLVER supports it
        (#73, closed by #77). But `design/catalog.py:77` still carries the interim
        guard from #71 and raises `DesignError` on any degenerate base before the
        solver ever sees the pattern, so a bare `TGNTATAAT` here failed 51
        end-to-end design tests. Declaring the expansion is the in-lane fix;
        this pins it."""
        forbidden = CrypticTranscription().lattice_terms(None).forbidden
        assert forbidden, "a HARD_LATTICE rule must declare motifs"
        for motif in forbidden:
            assert set(motif) <= set("ACGT"), f"{motif} is not pure ACGT"

    def test_does_not_list_tataat_because_the_closure_supplies_it(self) -> None:
        """`ATTATA` IS the reverse complement of `TATAAT`. Listing both would
        double-count, and CLAUDE.md 3.4 says list forward motifs only."""
        assert "TATAAT" not in CrypticTranscription().lattice_terms(None).forbidden


class TestMotifParts:
    def test_the_antisense_hexamer_is_a_breach(self) -> None:
        ev = run(CrypticTranscription(), construct("ATG" + PAD + "ATTATA" + PAD + "TAA"))
        assert [b.detail["motif"] for b in of_kind(ev, "motif")] == ["ATTATA"]
        assert not ev.passes

    def test_a_sense_strand_tataat_is_also_caught_by_the_closure(self) -> None:
        """The deliberate over-reach named in the module docstring: forbidding
        `ATTATA` necessarily forbids its reverse complement on the sense strand.
        Pinned so that removing it would fail a test rather than pass silently."""
        ev = run(CrypticTranscription(), construct("ATG" + PAD + "TATAAT" + PAD + "TAA"))
        assert of_kind(ev, "motif"), "TATAAT must be caught even though only ATTATA is listed"

    def test_the_sigma38_variant_is_a_breach(self) -> None:
        ev = run(CrypticTranscription(), construct("ATG" + PAD + "TATACT" + PAD + "TAA"))
        assert [b.detail["motif"] for b in of_kind(ev, "motif")] == ["TATACT"]

    def test_the_extended_minus_10_matches_every_n_expansion(self) -> None:
        """`TGnTATAAT` needs no -35, so each of the four expansions is its own
        promoter and each is declared in its own right. The message must name the
        substring actually present."""
        for base in "ACGT":
            motif = f"TG{base}TATAAT"
            ev = run(CrypticTranscription(), construct("ATG" + PAD + motif + PAD + "TAA"))
            found = [b for b in of_kind(ev, "motif") if b.detail["motif"] == motif]
            assert found, f"the extended -10 must match {motif}"
            assert found[0].detail["found"] == motif
            assert motif in found[0].message

    def test_a_motif_the_backbone_carries_is_reported_but_not_chased(self) -> None:
        c = construct("ATG" + PAD + "TAA", backbone=PAD + "TATACT" + PAD)
        breach = next(b for b in of_kind(run(CrypticTranscription(), c), "motif"))
        assert not breach.fixable_by_codon_choice

    def test_a_motif_spanning_the_origin_is_caught(self) -> None:
        """The case only the assembled circular construct shows: three bases
        close the sequence and three open it."""
        c = wrapping_construct(prefix=PAD, cds_head="ACT" + PAD, cds_tail=PAD + "TAT")
        ev = run(CrypticTranscription(), c)
        assert [b.detail["motif"] for b in of_kind(ev, "motif")] == ["TATACT"]


class TestSigma70Architecture:
    def _promoter(self, spacer_len: int, minus_35: str = "TTGACA") -> Construct:
        return construct("ATG" + PAD + minus_35 + "C" * spacer_len + "TATAAT" + PAD + "TAA")

    @pytest.mark.parametrize("spacer_len", [SPACER_MIN, 17, SPACER_MAX])
    def test_a_spacer_inside_the_range_is_a_promoter(self, spacer_len: int) -> None:
        ev = run(CrypticTranscription(), self._promoter(spacer_len))
        (breach,) = of_kind(ev, "sigma70_architecture")
        assert f"{spacer_len} bp spacer" in breach.message
        assert not ev.passes

    @pytest.mark.parametrize("spacer_len", [SPACER_MIN - 1, SPACER_MAX + 1])
    def test_a_spacer_outside_the_range_is_not(self, spacer_len: int) -> None:
        """The bounds are inclusive and exact. Geometry one base out is not a
        weaker promoter, it is a different thing, and reporting it would be
        noise the user cannot act on."""
        ev = run(CrypticTranscription(), self._promoter(spacer_len))
        assert of_kind(ev, "sigma70_architecture") == []

    def test_one_mismatch_in_the_minus_35_still_matches(self) -> None:
        """brief.md:111 says "within 1 mismatch", so an exact-only scan would
        miss most real cryptic promoters."""
        ev = run(CrypticTranscription(), self._promoter(17, minus_35="TTGACG"))
        assert of_kind(ev, "sigma70_architecture")

    def test_two_mismatches_in_the_minus_35_do_not(self) -> None:
        ev = run(CrypticTranscription(), self._promoter(17, minus_35="TTGCCG"))
        assert of_kind(ev, "sigma70_architecture") == []

    def test_the_breach_spans_both_elements_not_just_one(self) -> None:
        """The repair window has to be able to reach either element, so the
        breach must describe the whole architecture."""
        (breach,) = of_kind(run(CrypticTranscription(), self._promoter(17)), "sigma70_architecture")
        assert breach.interval.length == 6 + 17 + 6


class TestAtTract:
    def _tract(self, gap: int, tract: str = "TATTTAT") -> Construct:
        return construct("ATG" + PAD + tract + "C" * gap + "TATAAT" + PAD + "TAA")

    @pytest.mark.parametrize("tract", ["TATTTAT", "AATTT"])
    def test_a_tract_at_the_measured_offset_is_a_promoter(self, tract: str) -> None:
        ev = run(CrypticTranscription(), self._tract(AT_TRACT_GAP, tract))
        (breach,) = of_kind(ev, "at_tract")
        assert tract in breach.message
        assert "103/103" in breach.message

    @pytest.mark.parametrize("gap", [AT_TRACT_GAP - 1, AT_TRACT_GAP + 1])
    def test_a_tract_at_any_other_distance_is_not(self, gap: int) -> None:
        """ "Appropriately positioned" is the entire finding -- 103/103 held at
        this offset. A window instead of an exact offset would report the
        AT-richness of ordinary sequence."""
        assert of_kind(run(CrypticTranscription(), self._tract(gap)), "at_tract") == []

    @pytest.mark.parametrize(("tract", "gap_seq"), [("TATTTAT", "TGAC"), ("AATTT", "TGAC")])
    def test_the_papers_own_constructs_are_caught(self, tract: str, gap_seq: str) -> None:
        """The regression test for an off-by-one that shipped past a green suite.

        These are Warman & Grainger's Table 1 primers verbatim -- tract, then
        `TGAC`, then `TATAAT`. The 103/103 result was measured on exactly this
        geometry, so a rule citing it and not firing here is citing a result it
        cannot reproduce.

        `AT_TRACT_GAP` was 5, read straight from brief.md:111's "3' end sits
        exactly 5 bp upstream". That is a POSITION offset (-17 to -12), not a gap
        width: the paper's own primers put FOUR bases in between. At 5 the rule
        scored zero breaches on both constructs below, and fired instead on a
        one-base-shifted class nobody has tested."""
        assert len(gap_seq) == AT_TRACT_GAP, "the primer's gap is the rule's constant"
        c = construct("ATG" + PAD + tract + gap_seq + "TATAAT" + PAD + "TAA")
        assert of_kind(run(CrypticTranscription(), c), "at_tract")


class TestTheCalibrationAnchorIsNotATestVector:
    """The dengue-2 promoter is why the hazard is taken seriously. It is NOT
    something this rule detects, and that is pinned here so nobody later reads
    the citation as a recall claim."""

    def test_the_dengue_minus_35_is_three_mismatches_and_so_is_not_matched(self) -> None:
        """`TCAACG` vs `TTGACA` differs at three positions -- above the brief's
        1-mismatch budget, and above even the schema's maximum of 2."""
        assert sum(a != b for a, b in zip("TCAACG", "TTGACA", strict=True)) == 3
        for budget in (1, 2):
            assert "TCAACG" not in _within_one_mismatch("TTGACA", budget)

    def test_the_dengue_spacer_is_outside_the_range_the_brief_gives(self) -> None:
        """The published coordinates are -35 at nt 53 and -10 at nt 72, so the
        spacer is 13 bp -- not the 17 brief.md:111 states, and below SPACER_MIN
        either way. A second, independent reason the anchor cannot match."""
        assert 72 - (53 + 6) == 13
        assert not SPACER_MIN <= 13 <= SPACER_MAX


class TestScoredPathUnavailable:
    def test_every_run_says_the_promoter_calculator_did_not_run(self) -> None:
        """A rule that silently reports only what it could compute reads as a
        clean scan of everything it claims to cover."""
        clean = construct("ATG" + PAD + "TAA")
        (note,) = of_kind(run(CrypticTranscription(), clean), "unavailable")
        assert "Salis Promoter Calculator" in note.message
        assert "346 parameters" in note.message

    def test_the_note_is_weightless_and_never_chased(self) -> None:
        """Magnitude 0.0 and not fixable, so it lands on advisory rather than
        becoming a target the solver exhausts the mutation space against."""
        clean = construct("ATG" + PAD + "TAA")
        (note,) = of_kind(run(CrypticTranscription(), clean), "unavailable")
        assert note.magnitude == 0.0
        assert not note.fixable_by_codon_choice

    def test_a_clean_construct_still_passes(self) -> None:
        """The unavailability note must not, by itself, fail the rule -- that
        would make every design infeasible."""
        ev = run(CrypticTranscription(), construct("ATG" + PAD + "TAA"))
        assert ev.passes
        assert of_kind(ev, "motif") == []
        assert of_kind(ev, "sigma70_architecture") == []


class TestStrand:
    def test_the_geometry_reads_the_slots_strand_not_both(self) -> None:
        """A promoter reads one way, so this is a directional model and must
        follow `strand_for` (CLAUDE.md 3.4). A reverse-oriented cassette scans
        the other strand, and a promoter written on the forward strand is NOT
        one on the reverse -- the reverse-complement hazard is covered by the
        `ATTATA` motif instead, which is the point of it being a separate part."""
        c = construct("ATG" + PAD + "TTGACA" + SPACER + "TATAAT" + PAD + "TAA")
        forward = run(CrypticTranscription(), c, context(slot()))
        reverse = run(CrypticTranscription(), c, context(slot(), cassette_orientation=-1))
        assert of_kind(forward, "sigma70_architecture")
        assert of_kind(reverse, "sigma70_architecture") == []

    def test_minus_strand_positions_agree_between_interval_and_message(self) -> None:
        """A breach whose prose contradicts its own interval is not actionable.

        `emit` remaps the interval into construct coordinates, but the messages
        were formatted from raw reverse-strand scan indices -- so on a
        reverse-oriented cassette the -35 and -10 were printed at each other's
        positions."""
        promoter = "ATG" + PAD + "TTGACA" + SPACER + "TATAAT" + PAD + "TAA"
        c = construct(reverse_complement(promoter))
        ev = run(CrypticTranscription(), c, context(slot(), cassette_orientation=-1))
        (breach,) = of_kind(ev, "sigma70_architecture")
        reported = [int(tok) for tok in re.findall(r"hexamer at (\d+)", breach.message)]
        assert reported, "the message must name both element positions"
        for pos in reported:
            assert breach.interval.start <= pos < breach.interval.end, (
                f"position {pos} in the message falls outside {breach.interval}"
            )

    def test_the_breach_carries_the_slot_that_found_it(self) -> None:
        c = construct("ATG" + PAD + "TTGACA" + SPACER + "TATAAT" + PAD + "TAA")
        ev = run(CrypticTranscription(), c, context(slot(role="propagation")))
        (breach,) = of_kind(ev, "sigma70_architecture")
        assert breach.slot_role == "propagation"
