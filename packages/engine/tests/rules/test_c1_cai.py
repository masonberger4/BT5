"""C1: CAI as a band, the codons it must ignore, and the hosts it must refuse to score.

Two families of test carry most of the weight here, and neither is about arithmetic:

* **Unavailability.** `data/codon_usage/` ships one reference set against nine
  `HostId` values, so C1 has to say "unavailable" for seven of them. The failure
  this guards against is not a crash -- it is C1 quietly scoring a HEK293 design
  against E. coli weights and returning a number that looks entirely reasonable.
* **Never maximized.** The band's ceiling, the absent lattice term and the zero
  steering weight are one decision in three places. A test that only checked the
  floor would let all three drift toward a max-CAI optimizer.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
from bt5.codon.tables import FileTableProvider
from bt5.core.context import (
    BiosecurityVerdict,
    ContextSlot,
    DesignContext,
    HostId,
    Modality,
    SlotRole,
)
from bt5.core.registry import discover, get
from bt5.core.services import Services
from bt5.core.spec import Direction, Enforcement, Evaluation, Evidence
from bt5.core.types import Construct, Interval, Segment, SegmentKind, Strand, Topology
from bt5.rules.catalog.c1_cai import (
    BAND_HI,
    BAND_LO,
    CAI_BAND,
    CAI_REFERENCE_SET,
    CEILING_FRACTION_OF_HEADROOM,
    CodonAdaptationIndex,
    chance_cai,
)
from bt5.vector.kmers import ConstructKmerIndex
from conftest import construct, wrapping_construct

discover()

REFERENCE_SET = "sharp_li_1987_ecoli_w"

#: Sequences whose CAI against Sharp & Li's w-index is known, measured not guessed.
#: Each is ATG + 20 informative codons + TAA, so `n_evaluated` is 20 throughout and
#: a change in the exclusion logic shows up as a changed denominator, not only as a
#: changed score.
MAX_CAI = "ATG" + "GCTGAAGGTATCAAA" * 4 + "TAA"  # every codon its family's top: CAI 1.000
IN_BAND = "ATG" + "GCTGAGGGTATCAAA" * 4 + "TAA"  # CAI 0.763
RARE = "ATG" + "GCCGAGGGAATACTA" * 4 + "TAA"  # CAI 0.023


def _provider() -> FileTableProvider:
    return FileTableProvider()


def services(tables: object | None = None) -> Services:
    """Services carrying the REAL table provider unless a stub is passed.

    The shared `services` fixture in conftest stubs `weights()` to `{}`, which is a
    legitimate C1 input (it exercises the empty-reference-set path) but cannot
    exercise the arithmetic. Both are tested; neither substitutes for the other.
    """
    return Services(
        fold=None,
        kmer=ConstructKmerIndex,
        tables=tables or _provider(),  # type: ignore[arg-type]
        rng=np.random.default_rng(42),
    )


def slot(
    role: SlotRole = "producer",
    host: HostId = HostId.E_COLI_K12,
    modality: Modality = Modality.BACTERIAL_EXPRESSION,
    table: int = 11,
) -> ContextSlot:
    return ContextSlot(role, host, modality, table)


def context(*slots: ContextSlot, orientation: Strand = 1) -> DesignContext:
    return DesignContext(
        slots=slots or (slot(),),
        cassette_orientation=orientation,
        seed=42,
        screen=BiosecurityVerdict("not_run"),
    )


def evaluate(
    c: Construct,
    ctx: DesignContext | None = None,
    svc: Services | None = None,
    rule: CodonAdaptationIndex | None = None,
) -> Evaluation:
    return (rule or CodonAdaptationIndex()).evaluate(c, ctx or context(), svc or services())


def is_unavailable(ev: Evaluation) -> bool:
    """The shape B1 established: NaN, one breach, nothing the solver can chase."""
    return (
        math.isnan(ev.raw_score)
        and len(ev.breaches) == 1
        and ev.breaches[0].message.startswith("CAI objective unavailable:")
        and not ev.breaches[0].fixable_by_codon_choice
        and ev.n_evaluated == 0
        and ev.passes
    )


class TestSpecMetadata:
    def test_it_is_registered_under_its_brief_row(self) -> None:
        assert get("c1_cai") is CodonAdaptationIndex
        assert CodonAdaptationIndex.brief_ref == "2.C1"

    def test_it_is_the_soft_band_the_brief_section_header_asks_for(self) -> None:
        """brief.md:73 -- '2.C Codon composition (all S, soft bands, never
        maximized)'. Enforcement and Direction are that sentence, in code."""
        assert CodonAdaptationIndex.enforcement is Enforcement.SOFT
        assert CodonAdaptationIndex.direction is Direction.BAND

    def test_the_band_classvar_is_the_loosest_envelope_and_not_the_gate(self) -> None:
        """e2_gc_band's convention, and `solver/catalog.py:272` states it: "Read off
        the INSTANCE, never `Spec.band`." The bands are per host, so a single
        ClassVar can only be an envelope -- and it must be computed from `CAI_BAND`
        rather than transcribed, or adding a host silently makes it a lie."""
        assert CodonAdaptationIndex.band == (
            min(lo for lo, _ in CAI_BAND.values()),
            max(hi for _, hi in CAI_BAND.values()),
        )
        assert CodonAdaptationIndex.band == (0.0, 0.9553)
        assert "LOOSEST envelope" in CodonAdaptationIndex.weight_provenance

    def test_the_default_weight_matches_what_every_preset_already_assigns(self) -> None:
        """All three shipped presets weight 2.C1 at 0.2. A rule whose own default
        disagreed with the presets would make every preset an unexplained override."""
        assert CodonAdaptationIndex.default_weight == pytest.approx(0.2)

    def test_it_is_weighted_far_below_the_objective_with_a_measured_effect_size(self) -> None:
        """Kudla 2009 measured both on the same 154 variants: CAI r = 0.14 (not
        significant), 5' folding r = 0.66. The weights must reproduce that ordering."""
        b1 = get("b1_five_prime")
        assert CodonAdaptationIndex.default_weight < b1.default_weight

    def test_a_soft_rule_explains_its_weight_and_names_the_evidence(self) -> None:
        prov = CodonAdaptationIndex.weight_provenance
        assert "0.14" in prov, "the provenance must name the effect size it rests on"
        assert "0.66" in prov, "and the one it is being weighted against"

    def test_the_badge_is_contested_because_the_literature_is(self) -> None:
        """Kudla/Welch say CAI predicts nothing; Boel 2016 measured codon content as
        3-5x more influential than structure. BT5 cannot adjudicate that, and
        EVIDENCE_BACKED would badge the dispute settled."""
        assert CodonAdaptationIndex.evidence is Evidence.CONTESTED
        signs = {c.sign for c in CodonAdaptationIndex.citations}
        assert {"supports", "refutes"} <= signs, (
            "a CONTESTED badge is only honest if citations actually cut both ways"
        )

    def test_it_ships_enabled(self) -> None:
        """Only FOLKLORE defaults off, and all three presets weight this objective:
        shipping it disabled would silently drop a term every preset claims."""
        assert CodonAdaptationIndex.default_enabled

    def test_it_declares_the_repeat_rules_it_pulls_against(self) -> None:
        """Raising CAI collapses each family toward one codon, which is how a
        max-CAI sequence manufactures perfect repeats."""
        assert "e5_synthesis_repeats" in CodonAdaptationIndex.conflicts_with
        assert "f1_direct_repeats" in CodonAdaptationIndex.conflicts_with


class TestParamSchema:
    def test_neither_bound_advertises_a_default_it_cannot_honour(self) -> None:
        """e2_gc_band's reasoning, per host instead of per vendor: the gate is the
        slot host's own band, so a form materializing `default: 0.70` would show a
        HEK293 job E. coli's floor and then post it back as an explicit override."""
        props = CodonAdaptationIndex.param_schema["properties"]
        assert isinstance(props, dict)
        for name in ("cai_min", "cai_max"):
            assert "default" not in props[name], f"{name} must not advertise a default"
            assert "CAI_BAND" in props[name]["description"]


class TestNeverMaximized:
    """The three places the 'never maximized' decision has to hold at once."""

    def test_there_is_no_lattice_term(self) -> None:
        """A `codon_weights` term would be a monotone pull toward maximum CAI --
        the Tier-A DP maximizes what it is handed."""
        assert CodonAdaptationIndex().lattice_terms(None) is None

    def test_there_is_no_steering_term_either(self) -> None:
        """The other half of the same decision. A Lagrangian nudge on codon weights
        steers exactly as hard toward CAI 1.0 as a lattice term would."""
        assert CodonAdaptationIndex.steering_weight == 0.0

    def test_a_ceiling_at_one_is_refused(self) -> None:
        """brief.md:77 -- 'Never 1.0.' A band whose ceiling is 1.0 still looks
        two-sided while permitting the single-codon-per-amino-acid collapse."""
        with pytest.raises(ValueError, match="max-CAI collapse"):
            CodonAdaptationIndex(cai_max=1.0)

    def test_an_inverted_band_is_refused(self) -> None:
        with pytest.raises(ValueError, match="cai_min"):
            CodonAdaptationIndex(cai_min=0.9, cai_max=0.7)

    def test_the_ceiling_fires_on_a_perfectly_adapted_sequence(self) -> None:
        """The rule's whole point. Every codon is its family's top codon, CAI is
        1.0, and a rule that only had a floor would call this the ideal design."""
        ev = evaluate(construct(MAX_CAI))
        assert ev.raw_score == pytest.approx(1.0)
        assert ev.binding_side == "upper"
        assert not ev.passes
        assert "never a target to maximize" in ev.breaches[0].message

    def test_the_ceiling_breach_names_the_repeat_consequence(self) -> None:
        """A user told only 'CAI too high' has no reason to care. The measurable
        cost is the repeat structure it buys."""
        assert "repeats" in evaluate(construct(MAX_CAI)).breaches[0].message


class TestBand:
    def test_an_in_band_sequence_passes_with_no_breach(self) -> None:
        ev = evaluate(construct(IN_BAND))
        assert BAND_LO <= ev.raw_score <= BAND_HI
        assert ev.passes
        assert not ev.breaches

    def test_a_rare_codon_sequence_binds_on_the_floor(self) -> None:
        ev = evaluate(construct(RARE))
        assert ev.raw_score < BAND_LO
        assert ev.binding_side == "lower"
        assert "rare codons" in ev.breaches[0].message

    def test_magnitude_is_the_distance_outside_the_band(self) -> None:
        """Rule-native, and zero inside the band -- not |deviation from centre|,
        which would make a compliant design look like a small breach."""
        assert evaluate(construct(MAX_CAI)).breaches[0].magnitude == pytest.approx(0.1)
        assert evaluate(construct(RARE)).breaches[0].magnitude == pytest.approx(
            BAND_LO - 0.0231, abs=1e-3
        )

    def test_a_breach_is_fixable_by_codon_choice(self) -> None:
        """Codon choice is the ONLY thing that moves CAI, so a False here would
        route a trivially repairable finding to the advisor as a dead end."""
        assert evaluate(construct(MAX_CAI)).breaches[0].fixable_by_codon_choice

    def test_the_breach_carries_the_reference_set_it_was_scored_against(self) -> None:
        """CAI is meaningless without its reference set: the same CDS against
        another set is a different number, so one that travels without it is not
        reproducible."""
        detail = evaluate(construct(MAX_CAI)).breaches[0].detail
        assert detail["reference_set"] == REFERENCE_SET
        assert detail["host"] == HostId.E_COLI_K12
        assert detail["informative_codons"] == 20.0

    def test_a_widened_band_accepts_what_the_default_refuses(self) -> None:
        c = construct(RARE)
        assert not evaluate(c).passes
        assert evaluate(c, rule=CodonAdaptationIndex(cai_min=0.01, cai_max=0.99)).passes


class TestExcludedCodons:
    """Stops and single-codon families carry no codon-choice information."""

    def test_single_codon_families_change_neither_the_score_nor_the_denominator(self) -> None:
        """Four extra ATG/TGG codons are four more codons and zero more information."""
        plain = evaluate(construct(IN_BAND))
        padded = evaluate(construct("ATG" + "GCTGAGGGTATCAAA" * 4 + "ATGTGGATGTGG" + "TAA"))
        assert padded.raw_score == pytest.approx(plain.raw_score)
        assert padded.n_evaluated == plain.n_evaluated == 20

    def test_the_terminal_stop_is_excluded(self) -> None:
        """20 informative codons from a 22-codon ORF: ATG and TAA are both out."""
        assert evaluate(construct(IN_BAND)).n_evaluated == 20

    def test_single_codon_membership_is_read_from_the_table_not_hard_coded(self) -> None:
        """The case a hard-coded ('M', 'W') gets wrong. Under NCBI table 4 TGA is a
        second Trp codon, so Trp HAS a synonymous choice there and its codons carry
        information; under table 11 Trp is single-codon and TGG must drop out.
        Same codon, opposite treatment, decided by the table."""
        rule = CodonAdaptationIndex()
        provider = _provider()
        w = provider.weights(REFERENCE_SET, "cai")
        assert provider.genetic_code(4).synonymous_codons("W") == ("TGA", "TGG")
        assert rule._informative(("TGG",), provider.genetic_code(4), w) != []
        assert rule._informative(("TGG",), provider.genetic_code(11), w) == []

    def test_a_codon_with_no_weight_is_skipped_not_treated_as_zero(self) -> None:
        """ln(0) is -inf, which would collapse the whole geometric mean to 0.0 --
        a real, meaningful CAI value, reported for a missing weight."""
        rule = CodonAdaptationIndex()
        code = _provider().genetic_code(11)
        assert rule._informative(("GCT",), code, {"GCT": 0.0}) == []
        assert rule._informative(("GCT",), code, {}) == []


class TestUnavailable:
    """The three hosts BT5 still has no reference set for, and the honest report.

    Was seven. S6 shipped human, mouse and CHO sets (#90), so the mammalian hosts
    now compute -- see `TestMammalianHosts`. SF9, S_CEREVISIAE and P_PASTORIS were
    deferred deliberately: no shipped preset consumes them and their RefSeq
    coverage is materially messier, so they stay absent and stay honest.
    """

    @pytest.mark.parametrize(
        "host",
        [HostId.SF9, HostId.S_CEREVISIAE, HostId.P_PASTORIS],
    )
    def test_a_host_with_no_reference_set_reports_unavailable(self, host: HostId) -> None:
        """Scoring a CDS against another organism's weights because that table
        happened to be on disk is the failure this rule exists to refuse."""
        from bt5.core.context import LOCKED_TRANSLATION_TABLE

        ev = evaluate(
            construct(IN_BAND),
            context(slot("producer", host, Modality.LENTIVIRAL, LOCKED_TRANSLATION_TABLE[host])),
        )
        assert is_unavailable(ev)
        assert host in ev.breaches[0].message

    def test_unavailable_is_nan_and_never_a_plausible_number(self) -> None:
        """Every value in [0, 1] is a real CAI: 0.0 reads as a catastrophically
        rare-codon design and 0.8 reads as one exactly on target. Both would be
        affirmative false claims about a quantity nobody measured."""
        ev = evaluate(
            construct(IN_BAND), context(slot("producer", HostId.SF9, Modality.LENTIVIRAL, 1))
        )
        assert math.isnan(ev.raw_score)
        assert ev.binding_side is None

    def test_unavailable_does_not_read_as_a_breach_of_the_band(self) -> None:
        """`passes` stays True: the objective was not computed, which is a
        different statement from 'this construct is out of band'."""
        ev = evaluate(
            construct(IN_BAND),
            context(slot("producer", HostId.S_CEREVISIAE, Modality.LENTIVIRAL, 1)),
        )
        assert ev.passes
        assert ev.breaches[0].detail["unavailable_reason"]

    def test_every_mapped_host_has_a_table_actually_on_disk(self) -> None:
        """The guard on the map itself, and the one that matters when it grows.

        Adding a host here without its table in `data/codon_usage/` turns every
        honest 'unavailable' into a FileNotFoundError from inside a rule. Checked
        against the filesystem rather than a hard-coded list, so the guard keeps
        working as the map grows.
        """
        root = Path(__file__).resolve().parents[4]
        for host, key in CAI_REFERENCE_SET.items():
            path = root / "data" / "codon_usage" / f"{key}.json"
            assert path.exists(), f"{host} maps to {key!r}, which is not on disk"

    def test_the_deferred_hosts_are_absent_on_purpose(self) -> None:
        """Not an oversight: S6 deferred these three, and C1 says so rather than
        borrowing another organism's table."""
        for host in (HostId.SF9, HostId.S_CEREVISIAE, HostId.P_PASTORIS):
            assert host not in CAI_REFERENCE_SET


class TestMammalianHosts:
    """The four hosts S6's #90 took from `unavailable` to a real number.

    This is the wiring `docs/decisions/2026-09-02-s6-host-data-and-real-backbone.md`
    handed to S3: the data lane ships data, and `CAI_REFERENCE_SET` lives here.
    """

    @pytest.mark.parametrize(
        ("host", "expected_set"),
        [
            (HostId.HUMAN, "human_highly_expressed_refseq_w"),
            (HostId.HEK293, "human_highly_expressed_refseq_w"),
            (HostId.MOUSE, "mouse_highly_expressed_refseq_w"),
            (HostId.CHO, "cho_highly_expressed_refseq_w"),
        ],
    )
    def test_it_now_computes_against_the_right_reference_set(
        self, host: HostId, expected_set: str
    ) -> None:
        ev = evaluate(construct(IN_BAND), context(slot("producer", host, Modality.LENTIVIRAL, 1)))
        assert not is_unavailable(ev)
        assert not math.isnan(ev.raw_score)
        assert 0.0 < ev.raw_score <= 1.0
        assert CAI_REFERENCE_SET[host] == expected_set

    def test_hek293_shares_the_human_set_and_is_not_a_separate_table(self) -> None:
        """HEK293 is a Homo sapiens cell line, so this is the same same-species
        approximation BL21/K-12 already makes, one taxon down. Codon usage is a
        property of the organism's translational machinery, not the cell line.
        S6 shipped the human set FOR this mapping."""
        assert CAI_REFERENCE_SET[HostId.HEK293] == CAI_REFERENCE_SET[HostId.HUMAN]
        human = evaluate(
            construct(IN_BAND), context(slot("producer", HostId.HUMAN, Modality.LENTIVIRAL, 1))
        )
        hek = evaluate(
            construct(IN_BAND), context(slot("producer", HostId.HEK293, Modality.LENTIVIRAL, 1))
        )
        assert human.raw_score == pytest.approx(hek.raw_score)

    def test_a_mammalian_host_is_not_scored_against_e_coli(self) -> None:
        """The failure the unavailable path existed to prevent, now that a real
        alternative exists: a mammalian CAI must not be the E. coli number."""
        hek = evaluate(
            construct(IN_BAND), context(slot("producer", HostId.HEK293, Modality.LENTIVIRAL, 1))
        )
        ecoli = evaluate(
            construct(IN_BAND),
            context(slot("producer", HostId.E_COLI_K12, Modality.BACTERIAL_EXPRESSION, 11)),
        )
        assert hek.raw_score != pytest.approx(ecoli.raw_score)

    def test_the_mammalian_presets_no_longer_degrade_this_objective(self) -> None:
        """The point of the wiring. Both shipped mammalian presets key on HEK293,
        so before this they weighted 2.C1 at 0.2 and got `unavailable` back."""
        from bt5.score.presets import PRESETS, resolve

        for preset in PRESETS:
            resolved = resolve(preset)
            assert resolved.unimplemented == ()
            assert "c1_cai" in resolved.weights


class TestBandCalibration:
    """(0.70, 0.90) is E. coli's band and does not transfer unchanged.

    brief.md:77 offers those numbers as an EXAMPLE -- "e.g. 0.70-0.90, or +/-0.1
    of host median" -- and they were calibrated on a strong-bias organism. The
    mammalian w-tables are nearly flat (brief.md:206: "isochore GC, not
    selection"), so the same floor sits ~0.04 above chance instead of ~0.46 and
    stops discriminating. These tests re-derive every constant in `CAI_BAND` from
    the shipped tables so it cannot drift away from the data it came from.
    """

    def _chance(self, key: str, table_id: int) -> float:
        """The rule's own baseline, reached through its public helper.

        Deliberately NOT reimplemented here. The baseline the ceilings are scaled
        against has to use the same informative-family predicate as the metric being
        scored -- family size from the table first, positive weight second -- and a
        second copy in the test file is exactly how those two silently diverge, so
        that the number calibrating the ceiling stops describing the number the
        ceiling is applied to.
        """
        provider = FileTableProvider()
        return chance_cai(provider.weights(key, "cai"), provider.genetic_code(table_id))

    def _cds_near(self, target: float, key: str) -> tuple[str, float]:
        """A Leu-run CDS landing as near `target` as the family allows, plus its CAI.

        Built from the table rather than transcribed, so it cannot rot if the table
        is regenerated: a run of Leu codons split between the family's best and
        second-best, with the split solved for the target. Every codon is Leu, which
        has six codons under table 1, so all 300 are informative.
        """
        provider = FileTableProvider()
        w, code = provider.weights(key, "cai"), provider.genetic_code(1)
        best, alt = sorted(code.synonymous_codons("L"), key=lambda c: w[c])[:-3:-1]
        n = 300
        top, second = math.log(w[best]), math.log(w[alt])
        b = round(n * (math.log(target) - top) / (second - top))
        return "ATG" + alt * b + best * (n - b) + "TAA", math.exp((b * second + (n - b) * top) / n)

    def test_the_floor_barely_clears_chance_on_the_mammalian_tables(self) -> None:
        """The measurement the band change rests on. If this ever stops holding,
        the per-host band needs re-deriving rather than trusting."""
        ecoli = self._chance("sharp_li_1987_ecoli_w", 11)
        assert BAND_LO - ecoli > 0.4, "E. coli floor should sit far above chance"
        for key in (
            "human_highly_expressed_refseq_w",
            "mouse_highly_expressed_refseq_w",
            "cho_highly_expressed_refseq_w",
        ):
            assert BAND_LO - self._chance(key, 1) < 0.1, f"{key}: floor near chance"

    @pytest.mark.parametrize(
        ("host", "key"),
        [
            (HostId.HUMAN, "human_highly_expressed_refseq_w"),
            (HostId.HEK293, "human_highly_expressed_refseq_w"),
            (HostId.MOUSE, "mouse_highly_expressed_refseq_w"),
            (HostId.CHO, "cho_highly_expressed_refseq_w"),
        ],
    )
    def test_each_ceiling_is_that_hosts_share_of_its_own_headroom(
        self, host: HostId, key: str
    ) -> None:
        """Re-derives the constant rather than restating it: the ceiling is the
        same fraction of chance-to-1.0 that E. coli's 0.90 is of its own."""
        chance = self._chance(key, 1)
        expected = chance + CEILING_FRACTION_OF_HEADROOM * (1.0 - chance)
        assert CAI_BAND[host][1] == pytest.approx(expected, abs=5e-4)

    def test_the_ceiling_fraction_is_e_colis_published_0_90(self) -> None:
        ecoli = self._chance("sharp_li_1987_ecoli_w", 11)
        assert (
            pytest.approx((BAND_HI - ecoli) / (1.0 - ecoli), abs=5e-4)
            == CEILING_FRACTION_OF_HEADROOM
        )

    def test_e_colis_band_is_unchanged(self) -> None:
        """The published pair, untouched. A change that moved it would be a
        different rule, not a rescaling."""
        for host in (HostId.E_COLI_K12, HostId.E_COLI_BL21):
            assert CAI_BAND[host] == (BAND_LO, BAND_HI) == (0.70, 0.90)

    def test_the_rejected_rescaled_floor_is_re_derived_too(self) -> None:
        """The alternative the decision record turns down, pinned so it cannot drift
        while the accepted one stays fresh.

        Keeping E. coli's floor at the same share of headroom (0.6062) puts human's
        floor at ~0.864 -- above where a native human CDS sits -- so C1 would flag
        native sequence as "rare codons across the ORF" and hand the optimizer
        pressure to raise its CAI, which is what brief.md:13's Expi293F benchmark and
        brief.md:206 forbid. The shipped floor is inert instead, and the two must
        stay visibly different numbers.
        """
        ecoli = self._chance("sharp_li_1987_ecoli_w", 11)
        fraction = (BAND_LO - ecoli) / (1.0 - ecoli)
        human = self._chance("human_highly_expressed_refseq_w", 1)
        rescaled = human + fraction * (1.0 - human)
        assert rescaled == pytest.approx(0.864, abs=5e-4)
        assert CAI_BAND[HostId.HUMAN][0] == 0.0 < rescaled

    def test_the_two_host_maps_stay_in_step(self) -> None:
        """A host with a reference set but no calibrated band would otherwise be the
        next silent regression: added to one map, scored against whatever the other
        fell back to. `_band_for` refuses to guess, and this pins the pair so the
        refusal shows up here rather than as an unexplained `unavailable` in a run.
        """
        assert CAI_BAND.keys() == CAI_REFERENCE_SET.keys()

    def test_an_unmapped_host_gets_no_band_rather_than_e_colis(self) -> None:
        """The fallback that used to be there is the bug this map exists to fix."""
        assert CodonAdaptationIndex()._band_for(HostId.SF9) is None
        assert CodonAdaptationIndex()._band_for(HostId.S_CEREVISIAE) is None

    def test_the_mammalian_floor_is_inert_and_a_rare_cds_is_not_a_finding(self) -> None:
        """The behaviour the whole change exists for.

        A floor rescaled to keep E. coli's headroom would sit at ~0.864 -- above
        where a native human CDS sits -- so C1 would flag native sequence as
        "rare codons across the ORF" and hand the optimizer pressure to raise its
        CAI. brief.md:206 marks the CAI weight "very low" for CHO/HEK with default
        mode "Native or harmonize", and brief.md:13's Expi293F benchmark found
        optimization did not increase yields. So nothing here claims a low-CAI
        mammalian CDS is worse.
        """
        for host in (HostId.HEK293, HostId.MOUSE, HostId.CHO):
            assert CAI_BAND[host][0] == 0.0
            ev = evaluate(construct(RARE), context(slot("producer", host, Modality.LENTIVIRAL, 1)))
            assert ev.passes, f"{host}: an inert floor must not breach"
            assert ev.binding_side is None

    def test_the_mammalian_ceiling_still_bites(self) -> None:
        """Inert floor, live ceiling. Max-CAI collapse is MECHANICAL -- it drives
        each amino acid onto one codon and manufactures perfect direct repeats --
        so it transfers across organisms even where the floor does not."""
        provider = FileTableProvider()
        w, code = (
            provider.weights("human_highly_expressed_refseq_w", "cai"),
            provider.genetic_code(1),
        )
        # The human-optimal encoding of a Leu/Glu/Gly/Ile/Lys peptide, built from
        # the table so it cannot rot if the table is regenerated.
        best = {
            aa: max((c for c in code.synonymous_codons(aa) if c in w), key=lambda c: w[c])
            for aa in "LEGIK"
        }
        cds = "ATG" + "".join(best[a] for a in "LEGIK") * 4 + "TAA"
        ev = evaluate(
            construct(cds), context(slot("producer", HostId.HEK293, Modality.LENTIVIRAL, 1))
        )
        assert not ev.passes
        assert ev.binding_side == "upper"
        assert ev.breaches[0].detail["band_hi"] == CAI_BAND[HostId.HEK293][1]

    def test_the_widened_ceiling_is_what_actually_changed_for_hek293(self) -> None:
        """The one behavioural difference this change makes, in both directions.

        A max-CAI sequence breached before and breaches now, so it cannot tell the
        old band from the new one. The band moved for sequences BETWEEN the two
        ceilings: at CAI ~0.93 a HEK293 CDS breached E. coli's 0.90 and passes its
        own 0.9548, and just above 0.9548 it breaches again. Nothing else in the
        file pins where the ceiling sits.
        """
        hek = context(slot("producer", HostId.HEK293, Modality.LENTIVIRAL, 1))
        ceiling = CAI_BAND[HostId.HEK293][1]

        between, cai = self._cds_near(0.93, "human_highly_expressed_refseq_w")
        assert BAND_HI < cai < ceiling, f"fixture must sit between the ceilings: {cai}"
        assert evaluate(construct(between), hek).passes

        above, cai = self._cds_near(0.97, "human_highly_expressed_refseq_w")
        assert cai > ceiling, f"fixture must exceed the host ceiling: {cai}"
        ev = evaluate(construct(above), hek)
        assert not ev.passes
        assert ev.binding_side == "upper"

    def test_an_explicit_ceiling_equal_to_the_published_one_is_honoured(self) -> None:
        """A value-equality sentinel silently discarded this, and on the ceiling
        side it LOOSENED the limit: `cai_max=0.90` on a HEK293 job is a caller
        asking for the tighter published anti-max-CAI ceiling, and reading it as
        "unset" handed them 0.9548 instead -- permitting exactly the higher CAI the
        rule exists to refuse, through the rule's own parameter. The only way to get
        0.90 was to perturb the value.
        """
        between, cai = self._cds_near(0.93, "human_highly_expressed_refseq_w")
        c, hek = (
            construct(between),
            context(slot("producer", HostId.HEK293, Modality.LENTIVIRAL, 1)),
        )
        assert evaluate(c, hek).passes
        ev = evaluate(c, hek, rule=CodonAdaptationIndex(cai_max=BAND_HI))
        assert not ev.passes
        assert ev.binding_side == "upper"
        assert ev.breaches[0].detail["band_hi"] == BAND_HI

    def test_an_explicit_floor_equal_to_the_published_one_is_honoured(self) -> None:
        """The same defect on the side where it was merely a silent no-op. Asserted
        through `evaluate`, not `_band_for`: a band composed correctly and then not
        used would pass the implementation-level check."""
        hek = context(slot("producer", HostId.HEK293, Modality.LENTIVIRAL, 1))
        assert evaluate(construct(RARE), hek).passes
        ev = evaluate(construct(RARE), hek, rule=CodonAdaptationIndex(cai_min=BAND_LO))
        assert not ev.passes
        assert ev.binding_side == "lower"
        assert ev.breaches[0].detail["band_lo"] == BAND_LO

    def test_an_explicit_bound_overrides_only_its_own_side(self) -> None:
        """The user's own number wins, per side, so setting one does not discard
        the other -- e2_gc_band's discipline."""
        rule = CodonAdaptationIndex(cai_min=0.50)
        assert rule._band_for(HostId.HEK293) == (0.50, CAI_BAND[HostId.HEK293][1])
        assert rule._band_for(HostId.E_COLI_K12) == (0.50, BAND_HI)

    def test_a_band_that_only_inverts_once_composed_is_refused(self) -> None:
        """`__init__` cannot catch this -- the host is not known until evaluation --
        and a silently inverted band reports a max-CAI sequence as "rare codons"."""
        rule = CodonAdaptationIndex(cai_min=0.93)
        assert rule._band_for(HostId.HEK293) == (0.93, CAI_BAND[HostId.HEK293][1])
        with pytest.raises(ValueError, match="inverts"):
            rule._band_for(HostId.E_COLI_K12)


class TestGating:
    """Why the gate reads `role` -- the test that stops the wrong number shipping."""

    def test_an_ecoli_propagation_slot_does_not_supply_a_mammalian_job_with_a_cai(
        self,
    ) -> None:
        """THE case. A lentiviral job propagates its plasmid in E. coli and expresses
        the transgene in HEK293. E. coli is the one host BT5 has a w-table for, so a
        rule keyed on 'is any host E. coli' would find it, compute a confident CAI,
        and report it as the objective for a protein made in HEK293. The plasmid is
        only maintained in E. coli -- the transgene is never translated there."""
        ev = evaluate(
            construct(IN_BAND),
            context(
                slot("propagation", HostId.E_COLI_K12, Modality.LENTIVIRAL, 11),
                slot("producer", HostId.HEK293, Modality.LENTIVIRAL, 1),
            ),
        )
        # Before S6 shipped a human table this asserted `is_unavailable` -- the gate
        # was proved by ABSENCE, which stopped proving anything the moment HEK293
        # got a reference set. Now it is proved positively: the number exists, and
        # it was computed against the HUMAN set, never E. coli's.
        assert not is_unavailable(ev)
        assert ev.n_evaluated > 0

        # Which table answered, proved without needing a breach's detail: the
        # two-slot result must equal what HEK293 alone gives and differ from what
        # E. coli alone gives. The same CDS scored against a different reference
        # set is a different number -- that is why the set is part of CAI's
        # definition, and it is what makes this assertion bite.
        hek_only = evaluate(
            construct(IN_BAND),
            context(slot("producer", HostId.HEK293, Modality.LENTIVIRAL, 1)),
        )
        ecoli_only = evaluate(
            construct(IN_BAND),
            context(slot("producer", HostId.E_COLI_K12, Modality.BACTERIAL_EXPRESSION, 11)),
        )
        assert ev.raw_score == pytest.approx(hek_only.raw_score)
        assert ev.raw_score != pytest.approx(ecoli_only.raw_score)

    def test_a_propagation_only_context_has_nothing_to_score(self) -> None:
        ev = evaluate(
            construct(IN_BAND),
            context(slot("propagation", HostId.E_COLI_K12, Modality.LENTIVIRAL, 11)),
        )
        assert is_unavailable(ev)
        assert "no translating slot" in ev.breaches[0].message

    def test_producer_and_target_slots_are_both_scored(self) -> None:
        for role in ("producer", "target"):
            assert CodonAdaptationIndex().gate(slot(role))  # type: ignore[arg-type]
        assert not CodonAdaptationIndex().gate(slot("propagation"))

    def test_bacterial_expression_in_either_ecoli_strain_produces_a_number(self) -> None:
        """BL21 shares K-12's reference set deliberately -- a stated same-species
        approximation, not a silent one."""
        for host in (HostId.E_COLI_K12, HostId.E_COLI_BL21):
            ev = evaluate(construct(IN_BAND), context(slot("producer", host)))
            assert not math.isnan(ev.raw_score)
            assert ev.raw_score == pytest.approx(0.7632, abs=1e-3)


class TestServicesBoundary:
    """M4 reaches codon data only through the injected provider (CLAUDE.md 1)."""

    def test_the_rule_never_imports_the_codon_lane(self) -> None:
        """`Services` is what decouples the rules lane from M5. One convenience
        import of `CodonUsage.cai` would erase that and is invisible in behaviour."""
        root = Path(__file__).resolve().parents[4]
        source = (root / "packages/engine/src/bt5/rules/catalog/c1_cai.py").read_text()
        assert "bt5.codon" not in source
        assert "from bt5.core.services import" in source

    def test_it_computes_against_any_provider_honouring_the_protocol(self) -> None:
        """A stub with three codons and a real genetic code: exp(mean ln w) over
        two informative codons, computed from injected data alone."""

        class Stub:
            def genetic_code(self, table_id: int) -> object:
                return _provider().genetic_code(table_id)

            def usage(self, host: str) -> dict[str, float]:
                return {}

            def weights(self, host: str, kind: str) -> dict[str, float]:
                return {"GCT": 1.0, "GCC": 0.25, "ATG": 1.0}

        # ATG (excluded, single-codon) + GCT + GCC + TAA (excluded, stop).
        ev = evaluate(construct("ATGGCTGCCTAA"), svc=services(Stub()))
        assert ev.n_evaluated == 2
        assert ev.raw_score == pytest.approx(math.sqrt(1.0 * 0.25))

    def test_an_empty_reference_set_reports_unavailable_rather_than_zero(self) -> None:
        """conftest's shared `services` fixture stubs `weights()` to `{}`. A
        geometric mean over no weights is not 0.0; it is not defined."""

        class Empty:
            def genetic_code(self, table_id: int) -> object:
                return _provider().genetic_code(table_id)

            def usage(self, host: str) -> dict[str, float]:
                return {}

            def weights(self, host: str, kind: str) -> dict[str, float]:
                return {}

        ev = evaluate(construct(IN_BAND), svc=services(Empty()))
        assert is_unavailable(ev)
        assert "no relative adaptiveness weights" in ev.breaches[0].message

    def test_a_provider_that_raises_is_reported_not_propagated(self) -> None:
        """`FileTableProvider` raises FileNotFoundError for an absent table. A rule
        that let it escape would abort the whole design run over one objective."""

        class Raising:
            def genetic_code(self, table_id: int) -> object:
                return _provider().genetic_code(table_id)

            def usage(self, host: str) -> dict[str, float]:
                raise FileNotFoundError(f"no codon usage table for host {host!r}")

            def weights(self, host: str, kind: str) -> dict[str, float]:
                raise FileNotFoundError(f"no codon usage table for host {host!r}")

            def __repr__(self) -> str:
                return "Raising()"

        ev = evaluate(construct(IN_BAND), svc=services(Raising()))
        assert is_unavailable(ev)
        assert "could not be loaded" in ev.breaches[0].message


class TestConstructScope:
    """Evaluated on the assembled construct, never on a bare CDS string (CLAUDE.md 3.3)."""

    def test_a_cds_spanning_the_origin_is_read_as_one_orf(self) -> None:
        """The insert is stored as a single segment with end > length. A linear
        read would take the head and tail as two truncated, out-of-frame pieces."""
        ev = evaluate(wrapping_construct("ACGTACGTAC", "ATGGCTGAG", "GGTATCAAAGCTGAGGGT"))
        assert not math.isnan(ev.raw_score)
        assert ev.n_evaluated == 8, "9 codons across the origin, minus the ATG"

    def test_the_reverse_strand_is_read_when_the_cassette_is_cloned_in_backwards(
        self,
    ) -> None:
        """Directional and NOT revcomp-symmetric (CLAUDE.md 3.4): a reverse-cloned
        cassette's protein comes from the other strand, so the CAI must differ.
        Reading the plus strand regardless would score a protein nobody makes."""
        forward = evaluate(construct(IN_BAND), context(orientation=1))
        reverse = evaluate(construct(IN_BAND), context(orientation=-1))
        assert forward.raw_score != pytest.approx(reverse.raw_score)
        assert reverse.raw_score == pytest.approx(0.2545, abs=1e-3)

    def test_only_the_designable_cds_is_scored_not_the_backbone(self) -> None:
        """The vector is immutable, and its codon composition is not a fact about
        the protein being designed."""
        alone = evaluate(construct(IN_BAND))
        with_vector = evaluate(construct(IN_BAND, "GGGCCCGGGCCCGGGCCC" * 3))
        assert with_vector.raw_score == pytest.approx(alone.raw_score)

    def test_a_construct_with_no_designable_cds_reports_unavailable(self) -> None:
        bare = Construct(
            "ACGT" * 12,
            Topology.CIRCULAR,
            (Segment(Interval(0, 48), SegmentKind.BACKBONE, "vector"),),
        )
        ev = evaluate(bare)
        assert is_unavailable(ev)
        assert "no designable CDS" in ev.breaches[0].message

    def test_an_out_of_frame_cds_reports_unavailable_rather_than_a_number(self) -> None:
        """An out-of-frame CAI is a well-formed number about a protein that is not
        the one being made -- the most dangerous kind of wrong answer here."""
        ev = evaluate(construct(IN_BAND + "AA"))
        assert is_unavailable(ev)
        assert "whole number of codons" in ev.breaches[0].message

    def test_the_breach_spans_the_coding_scope(self) -> None:
        """WHOLE_SCOPE: no sub-interval is 'the' reason a geometric mean over the
        whole ORF left the band, so the finding is the CDS."""
        b = evaluate(construct(MAX_CAI)).breaches[0]
        assert b.interval == Interval(0, len(MAX_CAI))
        assert b.slot_role == "producer"


def test_it_applies_across_expression_modalities() -> None:
    """CAI is about translation, not delivery: a transgene is translated whether it
    arrived on a plasmid, a lentivirus or an AAV. The refusal for those modalities
    comes from the missing host table, not from the gate."""
    rule = CodonAdaptationIndex()
    for modality in Modality:
        assert rule.gate(slot("producer", HostId.HEK293, modality, 1))
