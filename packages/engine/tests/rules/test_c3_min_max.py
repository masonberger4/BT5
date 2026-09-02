"""C3: %MinMax, and the reason it reports unavailable on every host in this build.

The arithmetic is checked against Clarke & Clark 2008's own stated anchors --
+100 "only the most common codons", 0 "codon usage equal to the mean of all
possible codon choices", -100 "only the most rare codons" -- rather than against
a recorded snapshot, because a snapshot of a metric nobody can compute yet would
record whatever the first implementation happened to do.
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
)
from bt5.core.registry import discover, get
from bt5.core.services import GeneticCode, Services
from bt5.core.spec import Direction, Enforcement, Evaluation, Evidence
from bt5.core.types import Construct, Strand, Topology
from bt5.rules.catalog import c3_min_max
from bt5.rules.catalog.c3_min_max import (
    BAND_HI,
    BAND_LO,
    WINDOW_CODONS,
    MinMax,
    family_statistics,
    min_max_profile,
)
from bt5.vector.kmers import ConstructKmerIndex
from conftest import construct, wrapping_construct

discover()

#: The reference-set key the tests inject. Nothing ships under this name -- that
#: is the point: `MINMAX_REFERENCE_SET` is empty in the build, so every test that
#: wants the computing path has to say so explicitly.
FAKE_SET = "synthetic_frequencies_for_test"

#: E. coli's real NCBI table, and the one host+table pair that has ANY codon data.
ECOLI_TABLE = 11


# -- local helpers ------------------------------------------------------------
# Defined here, not in conftest: that file is shared with the other rules session
# and is read-only for both of us (docs/buildout/README.md).


def synthetic_frequencies(code: GeneticCode) -> dict[str, float]:
    """A complete codon -> usage frequency table with one clear favourite per family.

    Within every family the alphabetically first codon gets 0.7, the last 0.1 and
    the rest 0.2, so `most_common`/`rarest` below are predictable and a test can
    say "encode this entirely in the host's favourite codon" without hard-coding
    64 numbers that would silently rot if the table id changed.
    """
    freq: dict[str, float] = {}
    for a in "ACGT":
        for b in "ACGT":
            for c in "ACGT":
                codon = a + b + c
                if code.is_stop(codon):
                    freq[codon] = 0.0
                    continue
                aa = code.translate(codon)
                if not aa or aa == "*":
                    freq[codon] = 0.0
                    continue
                family = sorted(code.synonymous_codons(aa))
                if len(family) < 2:
                    freq[codon] = 1.0
                elif codon == family[0]:
                    freq[codon] = 0.7
                elif codon == family[-1]:
                    freq[codon] = 0.1
                else:
                    freq[codon] = 0.2
    return freq


def most_common(code: GeneticCode, aa: str) -> str:
    return sorted(code.synonymous_codons(aa))[0]


def rarest(code: GeneticCode, aa: str) -> str:
    return sorted(code.synonymous_codons(aa))[-1]


class FrequencyTables:
    """A `TableProvider` that honours the protocol's declared `usage` return type.

    `core/services.py:137` declares `usage(host) -> Mapping[str, float]`.
    `FileTableProvider.usage` returns a `CodonUsage` dataclass instead, so this
    stub is what a conforming provider looks like -- and the divergence is itself
    covered by `TestServicesBoundary`.
    """

    def __init__(self, freq: dict[str, float] | None = None, table_id: int = ECOLI_TABLE) -> None:
        self._code = FileTableProvider().genetic_code(table_id)
        self._freq = freq if freq is not None else synthetic_frequencies(self._code)

    def genetic_code(self, table_id: int) -> GeneticCode:
        return FileTableProvider().genetic_code(table_id)

    def usage(self, host: str) -> object:
        return self._freq

    def weights(self, host: str, kind: str) -> object:
        raise NotImplementedError(f"{kind} weights are not what C3 asks for")


def services(tables: object | None = None) -> Services:
    return Services(
        fold=None,
        kmer=ConstructKmerIndex,
        tables=tables or FrequencyTables(),  # type: ignore[arg-type]
        rng=np.random.default_rng(42),
    )


def slot(
    role: str = "producer",
    host: HostId = HostId.E_COLI_K12,
    modality: Modality = Modality.BACTERIAL_EXPRESSION,
    table: int = ECOLI_TABLE,
    strand_of_interest: Strand = 1,
) -> ContextSlot:
    return ContextSlot(role, host, modality, table, strand_of_interest)  # type: ignore[arg-type]


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
    rule: MinMax | None = None,
) -> Evaluation:
    return (rule or MinMax()).evaluate(c, ctx or context(), svc or services())


def is_unavailable(ev: Evaluation) -> bool:
    """The shape B1 established and C1 followed: NaN, one unfixable breach."""
    return (
        math.isnan(ev.raw_score)
        and len(ev.breaches) == 1
        and ev.breaches[0].message.startswith("%MinMax objective unavailable:")
        and not ev.breaches[0].fixable_by_codon_choice
        and ev.n_evaluated == 0
        and ev.passes
    )


@pytest.fixture
def wired(monkeypatch: pytest.MonkeyPatch) -> dict[HostId, str]:
    """Give every host a reference set, so the computing path can be exercised.

    Without this the rule is unavailable everywhere by design, and every band
    test below would pass for the wrong reason.
    """
    installed = dict.fromkeys(HostId, FAKE_SET)
    monkeypatch.setattr(c3_min_max, "MINMAX_REFERENCE_SET", installed)
    return installed


# -- the metric ---------------------------------------------------------------


class TestFamilyStatistics:
    def test_it_excludes_stops_and_single_codon_families(self) -> None:
        """Neither carries codon-choice information, so neither is in %MinMax."""
        code = FileTableProvider().genetic_code(ECOLI_TABLE)
        stats = family_statistics(code, synthetic_frequencies(code))
        for stop in ("TAA", "TGA", "TAG"):
            assert stop not in stats
        assert "ATG" not in stats  # Met: one codon under table 11
        assert "TGG" not in stats  # Trp: one codon under table 11
        assert "CTG" in stats  # Leu: six codons, informative

    def test_family_size_is_read_from_the_injected_table_not_hard_coded(self) -> None:
        """NCBI table 4 makes TGA a second Trp codon, so Trp IS informative there.

        Hard-coding ATG/TGG as the excluded families -- the obvious shortcut --
        would silently drop a real choice under table 4 (CLAUDE.md 3.1).
        """
        provider = FileTableProvider()
        table4 = provider.genetic_code(4)
        stats4 = family_statistics(table4, synthetic_frequencies(table4))
        assert table4.translate("TGA") == "W"
        assert "TGG" in stats4
        assert "TGA" in stats4

        table11 = provider.genetic_code(ECOLI_TABLE)
        stats11 = family_statistics(table11, synthetic_frequencies(table11))
        assert "TGG" not in stats11

    def test_a_row_is_the_codon_then_its_family(self) -> None:
        code = FileTableProvider().genetic_code(ECOLI_TABLE)
        freq = synthetic_frequencies(code)
        stats = family_statistics(code, freq)
        top = most_common(code, "L")
        used, mean, peak, trough = stats[top]
        assert used == freq[top] == 0.7
        assert peak == 0.7
        assert trough == 0.1
        assert trough < mean < peak

    def test_a_partly_covered_family_is_dropped_whole(self) -> None:
        """A partial family's mean is not the family's mean.

        Averaging over only the codons a reference set happens to list would ride
        silently into every window overlapping that amino acid.
        """
        code = FileTableProvider().genetic_code(ECOLI_TABLE)
        freq = synthetic_frequencies(code)
        del freq[most_common(code, "L")]
        stats = family_statistics(code, freq)
        assert not any(stats.get(codon) for codon in code.synonymous_codons("L"))
        assert stats.get(most_common(code, "R")) is not None


class TestProfile:
    """Clarke & Clark 2008's three defined anchors, and the window's behaviour."""

    @pytest.fixture
    def stats(self):  # noqa: ANN201
        code = FileTableProvider().genetic_code(ECOLI_TABLE)
        return code, family_statistics(code, synthetic_frequencies(code))

    def test_only_the_most_common_codons_is_plus_100(self, stats) -> None:
        code, table = stats
        profile = min_max_profile([most_common(code, "L")] * 40, table)
        assert profile[20] == pytest.approx(100.0)

    def test_only_the_rarest_codons_is_minus_100(self, stats) -> None:
        code, table = stats
        profile = min_max_profile([rarest(code, "L")] * 40, table)
        assert profile[20] == pytest.approx(-100.0)

    def test_a_rare_cluster_is_a_negative_peak(self, stats) -> None:
        """The paper's plotting convention, and what the band's sign depends on."""
        code, table = stats
        codons = [most_common(code, "L")] * 20 + [rarest(code, "L")] * 20
        profile = min_max_profile(codons, table, window=3)
        assert profile[3] > 0
        assert profile[36] < 0

    def test_the_window_shrinks_at_termini_rather_than_dropping_positions(self, stats) -> None:
        """brief.md:79. The profile is exactly as long as the CDS, ends included."""
        code, table = stats
        codons = [most_common(code, "L")] * 40
        assert len(min_max_profile(codons, table)) == 40
        assert min_max_profile(codons, table)[0] == pytest.approx(100.0)
        assert min_max_profile(codons, table)[-1] == pytest.approx(100.0)
        # Shorter than one window: still one value per codon, none dropped.
        assert len(min_max_profile(codons[:3], table, window=WINDOW_CODONS)) == 3

    def test_uninformative_codons_score_zero_not_undefined(self, stats) -> None:
        """Nothing in a poly-Trp window could have been chosen differently."""
        code, table = stats
        assert set(min_max_profile(["TGG"] * 20, table)) == {0.0}

    def test_the_two_equations_have_different_denominators(self, stats) -> None:
        """%Max divides by (X_max - X_avg), %Min by (X_avg - X_min).

        So an equal split of common and rare codons does NOT average to zero.
        This is Clarke & Clark's definition, not a rounding artefact -- a future
        reader "fixing" the asymmetry would be changing the metric.
        """
        code, table = stats
        codons = [most_common(code, "L")] * 20 + [rarest(code, "L")] * 20
        profile = min_max_profile(codons, table, window=3)
        assert abs(sum(profile) / len(profile)) > 1.0

    def test_a_window_wider_than_the_cds_flattens_the_profile(self, stats) -> None:
        code, table = stats
        codons = [most_common(code, "L")] * 10 + [rarest(code, "L")] * 10
        assert len({round(v, 9) for v in min_max_profile(codons, table, window=500)}) == 1

    def test_window_must_be_at_least_one_codon(self, stats) -> None:
        _, table = stats
        with pytest.raises(ValueError, match="at least 1 codon"):
            min_max_profile(["CTG"], table, window=0)


# -- the rule -----------------------------------------------------------------


class TestSpecMetadata:
    def test_the_band_edges_are_the_metrics_own_definition(self) -> None:
        """Neither edge is a threshold anybody picked, which is why they are safe.

        0 is Clarke & Clark's "codon usage equal to the mean of all possible codon
        choices"; -100 is the metric's definitional minimum and therefore cannot
        be breached. If either ever becomes a chosen number it needs its own
        evidence.
        """
        assert MinMax.band == (-100.0, 0.0)
        assert (BAND_LO, BAND_HI) == (-100.0, 0.0)

    def test_the_floor_can_never_bind(self) -> None:
        """A profile value below -100 is not representable, so the band is one-sided.

        Deliberate: brief.md:84 (C8) asks for native rare-codon clusters to be
        RETAINED, and Clarke & Clark's finding is that they cluster functionally.
        """
        assert MinMax().min_max_min == -100.0
        assert MinMax()._side(-100.0) is None

    def test_it_is_soft_and_never_steers(self) -> None:
        """steering_weight 0.0 is structural: see lattice_terms and C1's argument."""
        assert MinMax.enforcement is Enforcement.SOFT
        assert MinMax.direction is Direction.BAND
        assert MinMax.steering_weight == 0.0
        assert MinMax().lattice_terms(context()) is None

    def test_the_evidence_badge_matches_the_briefs_own_grade(self) -> None:
        """brief.md:79 grades C3 'B' = replicated but contested or single-lab."""
        assert MinMax.evidence is Evidence.CONTESTED

    def test_the_default_weight_matches_every_shipped_preset(self) -> None:
        from bt5.score.presets import PRESETS

        weighted = [
            entry.weight
            for preset in PRESETS
            for entry in preset.entries
            if entry.brief_ref == "2.C3"
        ]
        assert weighted == [0.3, 0.3, 0.3]
        assert MinMax.default_weight == 0.3

    def test_the_window_is_the_value_its_primary_source_measured_on(self) -> None:
        """Clarke & Clark used 18; PLAN.md:661 pins 18. Not the brief's 17-18 range."""
        assert WINDOW_CODONS == 18
        assert MinMax().window == 18
        assert MinMax.param_schema["properties"]["window"]["default"] == 18  # type: ignore[index]

    def test_it_rejects_an_inverted_or_out_of_range_band(self) -> None:
        with pytest.raises(ValueError, match="min_max_min"):
            MinMax(min_max_min=10.0, min_max_max=-10.0)
        with pytest.raises(ValueError, match="min_max_min"):
            MinMax(min_max_min=-500.0)
        with pytest.raises(ValueError, match="at least 1 codon"):
            MinMax(window=0)


class TestUnavailableInThisBuild:
    """The headline: no shipped reference set carries usage FREQUENCIES."""

    @pytest.mark.parametrize("host", list(HostId))
    def test_every_host_reports_unavailable(self, host: HostId) -> None:
        from bt5.core.context import LOCKED_TRANSLATION_TABLE

        code = FileTableProvider().genetic_code(LOCKED_TRANSLATION_TABLE[host])
        cds = most_common(code, "L") * 30
        ctx = context(slot(host=host, table=LOCKED_TRANSLATION_TABLE[host]))
        assert is_unavailable(evaluate(construct(cds), ctx))

    def test_e_coli_is_unavailable_too_and_that_is_the_finding(self) -> None:
        """C1 computes for the two E. coli hosts; C3 cannot, and the reason differs.

        C1 needs relative adaptiveness w, which is exactly what ships. %MinMax
        needs raw usage FREQUENCIES, and w normalises each family to its own peak
        and discards the peak, so the per-family scale %MinMax sums across
        families is not recoverable. The brief's premise that C3 would merely
        mirror C1's seven-host gap is what this test records as wrong.
        """
        ev = evaluate(construct("CTG" * 30), context(slot(host=HostId.E_COLI_K12)))
        assert is_unavailable(ev)
        assert "FREQUENCY" in ev.breaches[0].message
        assert "relative" in ev.breaches[0].message

    def test_the_reason_travels_where_m3_reads_it(self) -> None:
        ev = evaluate(construct("CTG" * 30))
        assert ev.breaches[0].detail["unavailable_reason"]
        assert ev.breaches[0].magnitude == 0.0

    def test_unavailable_is_not_a_breach_of_the_band(self) -> None:
        """`passes` stays True: nothing failed, an objective was not computed."""
        ev = evaluate(construct("CTG" * 30))
        assert ev.passes is True
        assert ev.binding_side is None

    def test_the_map_is_the_only_wiring_a_frequency_table_needs(self, wired: dict) -> None:
        """One entry in MINMAX_REFERENCE_SET makes the host live, and nothing else."""
        ev = evaluate(construct("CTG" * 30))
        assert not is_unavailable(ev)
        assert not math.isnan(ev.raw_score)


class TestBand:
    def test_all_favourite_codons_breaches_the_ceiling(self, wired: dict) -> None:
        code = FileTableProvider().genetic_code(ECOLI_TABLE)
        ev = evaluate(construct(most_common(code, "L") * 40))
        assert ev.raw_score == pytest.approx(100.0)
        assert ev.binding_side == "upper"
        assert not ev.passes
        assert ev.breaches[0].fixable_by_codon_choice

    def test_all_rare_codons_does_not_breach(self, wired: dict) -> None:
        """The floor is non-binding by construction -- the design decision.

        A -100 CDS is not something this rule objects to; C8 wants rare-codon
        structure retained, and nothing in the catalog's evidence penalises it.
        """
        code = FileTableProvider().genetic_code(ECOLI_TABLE)
        ev = evaluate(construct(rarest(code, "L") * 40))
        assert ev.raw_score == pytest.approx(-100.0)
        assert ev.binding_side is None
        assert ev.passes
        assert ev.breaches == ()

    def test_the_breach_message_never_predicts_expression(self, wired: dict) -> None:
        """CLAUDE.md: a CI gate bans prediction vocabulary. Say it here too."""
        code = FileTableProvider().genetic_code(ECOLI_TABLE)
        message = evaluate(construct(most_common(code, "L") * 40)).breaches[0].message
        lowered = message.lower()
        for banned in ("will increase", "predicted", "fold-improvement", "yield of"):
            assert banned not in lowered

    def test_the_profile_travels_for_the_report(self, wired: dict) -> None:
        """PLAN.md:661 asks for the profile to be printed, so it must be emitted."""
        code = FileTableProvider().genetic_code(ECOLI_TABLE)
        ev = evaluate(construct(most_common(code, "L") * 40))
        assert len(ev.windows) == 40
        assert ev.n_evaluated == 40
        assert all(iv.length == 3 for iv, _ in ev.windows)

    def test_the_worst_slot_binds_rather_than_the_average(self, wired: dict) -> None:
        """Averaging across slots would let a comfortable slot hide one at +100."""
        code = FileTableProvider().genetic_code(ECOLI_TABLE)
        mixed = context(
            slot(role="producer"),
            slot(role="target", host=HostId.E_COLI_BL21),
        )
        ev = evaluate(construct(most_common(code, "L") * 40), mixed)
        assert ev.binding_side == "upper"

    def test_detail_carries_the_reference_set(self, wired: dict) -> None:
        """The same CDS against another reference set is a different number."""
        code = FileTableProvider().genetic_code(ECOLI_TABLE)
        detail = evaluate(construct(most_common(code, "L") * 40)).breaches[0].detail
        assert detail["reference_set"] == FAKE_SET
        assert detail["window_codons"] == 18.0


class TestGating:
    def test_a_propagation_slot_is_not_translated(self) -> None:
        """C1's argument: a lentiviral job propagates in E. coli, expresses in HEK293.

        Keying on host would find the one host with codon data and report its
        number as the objective for a protein made somewhere else.
        """
        assert not MinMax().gate(slot(role="propagation"))
        assert MinMax().gate(slot(role="producer"))
        assert MinMax().gate(slot(role="target"))

    def test_propagation_only_context_is_unavailable_not_zero(self, wired: dict) -> None:
        ctx = context(slot(role="propagation"))
        ev = evaluate(construct("CTG" * 30), ctx)
        assert is_unavailable(ev)
        assert "propagation-only" in ev.breaches[0].message

    def test_it_is_scored_in_every_modality_every_preset_pins(self) -> None:
        """resolve() probes every admitted slot; a non-scored answer is a PresetError."""
        for modality in Modality:
            assert MinMax().enforcement_for(slot(modality=modality)) is Enforcement.SOFT


class TestConstructScope:
    def test_an_origin_spanning_cds_is_read_as_one_span(self, wired: dict) -> None:
        """Rules take a Construct precisely so this cannot be missed (CLAUDE.md 3.3)."""
        code = FileTableProvider().genetic_code(ECOLI_TABLE)
        top = most_common(code, "L")
        c = wrapping_construct("AAAAAA", top * 10, top * 10)
        ev = evaluate(c)
        assert not is_unavailable(ev)
        assert ev.raw_score == pytest.approx(100.0)

    def test_a_reverse_strand_cds_is_read_in_reading_order(self, wired: dict) -> None:
        """Not revcomp-symmetric: reading it forwards scores a protein nobody makes.

        Arg deliberately, not Leu: Leu's favourite codon here is CTA, whose reverse
        complement TAG is a STOP under table 11, so the reverse case would report
        "no informative codons" and the test would pass without ever exercising the
        reading-order path. AGA reverse-complements to TCT, a real Ser codon -- and
        the rarest of its family, so the two strands land at opposite ends of the
        scale rather than merely differing.
        """
        code = FileTableProvider().genetic_code(ECOLI_TABLE)
        forward = evaluate(construct(most_common(code, "R") * 20), context(slot()))
        reverse = evaluate(construct(most_common(code, "R") * 20), context(slot(), orientation=-1))
        assert not is_unavailable(forward)
        assert not is_unavailable(reverse)
        assert forward.raw_score == pytest.approx(100.0)
        assert reverse.raw_score == pytest.approx(-100.0)

    def test_an_out_of_frame_scope_is_unavailable_not_wrong(self, wired: dict) -> None:
        ev = evaluate(construct("CTGCTGCTGC"))
        assert is_unavailable(ev)
        assert "whole number of codons" in ev.breaches[0].message

    def test_no_designable_cds_is_unavailable(self, wired: dict) -> None:
        c = Construct(sequence="ACGT" * 10, topology=Topology.CIRCULAR, segments=())
        ev = evaluate(c)
        assert is_unavailable(ev)
        assert "no designable CDS" in ev.breaches[0].message


class TestServicesBoundary:
    def test_the_rule_never_imports_the_codon_lane(self) -> None:
        """`Services` is what decouples M4 from M5; an import edge would bypass it."""
        root = Path(__file__).resolve().parents[4]
        source = (root / "packages/engine/src/bt5/rules/catalog/c3_min_max.py").read_text()
        assert "bt5.codon" not in source
        assert "from bt5.core.services import" in source

    def test_a_provider_diverging_from_its_protocol_is_reported_not_crashed(
        self, wired: dict
    ) -> None:
        """`FileTableProvider.usage` returns CodonUsage, not the declared Mapping.

        A rule taking the protocol at its word must surface that as a stated
        unavailability, not an AttributeError from inside evaluate().
        """

        class Diverging(FrequencyTables):
            def usage(self, host: str):  # noqa: ANN201
                return FileTableProvider().usage("sharp_li_1987_ecoli_w")

        ev = evaluate(construct("CTG" * 30), svc=services(Diverging()))
        assert is_unavailable(ev)
        assert "codon to usage frequency" in ev.breaches[0].message

    def test_a_missing_table_is_reported_not_raised(self, wired: dict) -> None:
        class Missing(FrequencyTables):
            def usage(self, host: str):  # noqa: ANN201
                raise FileNotFoundError(f"no table for {host}")

        ev = evaluate(construct("CTG" * 30), svc=services(Missing()))
        assert is_unavailable(ev)
        assert "could not be loaded" in ev.breaches[0].message

    def test_an_empty_frequency_table_is_reported(self, wired: dict) -> None:
        ev = evaluate(construct("CTG" * 30), svc=services(FrequencyTables(freq={})))
        assert is_unavailable(ev)


# -- registration -------------------------------------------------------------


def test_it_is_registered_under_its_brief_row() -> None:
    assert get("c3_min_max") is MinMax
    assert MinMax.brief_ref == "2.C3"
    assert MinMax.id == Path(c3_min_max.__file__).stem


def test_no_shipped_reference_set_carries_frequencies() -> None:
    """The build-state assertion this whole rule is shaped around.

    When a frequency table lands this test is what tells the next session to
    revisit the rule's docstring and this file, rather than leaving an
    always-unavailable objective in place after it stopped being necessary.
    """
    assert c3_min_max.MINMAX_REFERENCE_SET == {}


def test_it_closes_the_last_unimplemented_objective() -> None:
    """2.C3 was the only brief_ref no rule claimed; shipping it empties the set."""
    from bt5.score.presets import PRESETS, resolve

    for preset in PRESETS:
        assert resolve(preset).unimplemented == ()
