"""The ranking increment, end to end: ranks, a panel, a baseline, an order file.

The walking skeleton's own docstring listed what it refused to do -- one
candidate, no gallery, every objective `unavailable`, `native_baseline` None, no
order CSV, no percentiles. This file is the inverse of that list, plus the
honesty each item has to survive.

The load-bearing tests here are the negative ones. `design()` producing five
sequences is easy; five sequences that are all the same design, ranked by a
percentile measured against a null that was never built, with a report that
calls itself complete, would pass a naive test of every positive claim above.
"""

from __future__ import annotations

import csv
import io
import math
from typing import Any

import pytest
from bt5.codon.tables import FileTableProvider
from bt5.core.context import HostId, Modality
from bt5.core.result import ObjectiveScore
from bt5.core.types import DNA_ALPHABET
from bt5.design import DesignError, design
from bt5.design.errors import DesignError as DesignErrorAlias
from bt5.design.gallery import DEFAULT_SWEEP_STEPS, SolveSpace, sweep_designs
from bt5.design.runner import DEFAULT_GALLERY_SIZE, UNSCREENED
from bt5.score.gallery import G4_MIN_PAIRWISE_DISTANCE, MAX_GALLERY, MIN_GALLERY
from bt5.score.null import DEFAULT_NULL_N
from bt5.score.steering import SWEEP_AXES, live_axes
from bt5.vector.backbone import VectorBackbone

CODE = FileTableProvider().genetic_code(1)


def _solve_space(backbone: VectorBackbone, protein: str) -> SolveSpace:
    """The same `SolveSpace` `design()` builds, for tests that need to sweep it
    directly rather than through the whole pipeline."""
    from bt5.design.catalog import partition_forbidden
    from bt5.design.runner import _coding_flanks, _context
    from bt5.design.sites import choose_site
    from bt5.solver.catalog import build_rule_set, default_services
    from bt5.vector import assemble

    services = default_services(seed=0)
    site = choose_site(backbone, table_id=1)
    ctx = _context(
        modality=Modality.LENTIVIRAL,
        hosts=[HostId.HEK293],
        table_id=1,
        cassette_orientation=site.strand,
    )
    rules = build_rule_set(ctx, services)
    usable, _carried = partition_forbidden(rules.forbidden(), backbone, site)
    left, right = _coding_flanks(backbone, site, usable)
    placeholder = "ATG" + "AAA" * (len(protein) - 1) + "TAA"
    return SolveSpace(
        protein=protein,
        code=CODE,
        assemble=lambda cds: (
            assemble(backbone, cds, protein=protein, table_id=1, site=site).construct
        ),
        forbidden=usable,
        seed=0,
        table_id=1,
        usage={},
        find_breaches=rules.breach_finder(),
        gc_bounds=rules.oracle_bounds().gc_bounds,
        left_flank=left,
        right_flank=right,
        policies=rules.policies(50),
        reference=assemble(backbone, placeholder, protein=protein, table_id=1, site=site).reference,
    )


#: Every degradation `design()` is allowed to emit, as the leading text of the
#: sentence. A NEW source of degradation that is not one of these fails
#: `test_no_degradation_arrives_unremarked` -- which is the property the
#: skeleton's set-equality test was protecting, kept while letting the CONTENT
#: vary with the environment (ViennaRNA present or not, host tables shipped or
#: not) instead of pinning a set that only holds on one machine.
KNOWN_DEGRADATIONS = (
    "ViennaRNA",
    "protein-level biosecurity screening:",
    "the ",  # the G4 shortfall and short-panel sentences
    "single candidate only:",
    "no codon usage table on file",
    "no native baseline:",
    "the native CDS was supplied",
    "objective ",
    "rule ",
    "forbidden motif ",
    "screening burden unavailable:",
)


def test_the_alias_is_the_same_error() -> None:
    """`bt5.design.DesignError` and `bt5.design.errors.DesignError` must not
    drift apart -- S5's CLI catches one of them."""
    assert DesignError is DesignErrorAlias


class TestTheGallery:
    def test_it_ships_a_panel_not_a_single_candidate(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        res = fast(backbone)
        assert len(res.result.candidates) >= 1
        assert res.report.candidates == len(res.result.candidates)
        if res.gallery is not None:
            assert len(res.result.candidates) == len(res.gallery.picks)
            assert res.gallery.swept >= res.gallery.distinct

    def test_every_candidate_encodes_the_protein_and_is_distinct(
        self, backbone: VectorBackbone, fast: Any, protein: str
    ) -> None:
        """Distinct is not decoration: `build_gallery` deduplicates before
        selecting, so two identical CDSs in the panel would mean the dedup was
        bypassed and a slot was spent proving two designs are the same one."""
        res = fast(backbone)
        seen = {candidate.cds for candidate in res.result.candidates}
        assert len(seen) == len(res.result.candidates)
        for candidate in res.result.candidates:
            assert CODE.translate(candidate.cds)[:-1] == protein

    def test_a_g4_shortfall_is_reported_and_never_lowered(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """G4's failure invalidates a PRODUCT decision, not a technical one, so
        the response is a degradation naming the distance actually reached --
        never a smaller threshold. The threshold itself is pinned here so a later
        change to it has to come through this test."""
        assert G4_MIN_PAIRWISE_DISTANCE == 0.15
        res = fast(backbone)
        if res.meets_g4:
            assert res.gallery.min_pairwise_distance >= G4_MIN_PAIRWISE_DISTANCE
        else:
            assert any(
                "does not meet gate G4" in d or "single candidate only" in d
                for d in res.result.provenance.degradations
            )

    def test_candidates_are_ranked_best_first(self, backbone: VectorBackbone, fast: Any) -> None:
        res = fast(backbone)
        totals = [candidate.scorecard.total for candidate in res.result.candidates]
        assert totals == sorted(totals, reverse=True)
        assert [c.label for c in res.result.candidates] == [
            f"design_{i + 1}" for i in range(len(totals))
        ]

    def test_the_exported_genbank_is_the_top_ranked_candidate(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """`assembly`, `optimize_result` and `genbank` all describe
        `candidates[0]`. A GenBank of a candidate other than the one the report
        ranks first is the worst possible failure here: it is silent, and the
        user orders it."""
        res = fast(backbone)
        winner = res.result.candidates[0]
        assert res.assembly.construct.sequence == winner.construct.sequence
        assert res.optimize_result.cds == winner.cds
        assert winner.design_hash in res.genbank

    def test_codon_distances_are_recorded_against_every_other_candidate(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        res = fast(backbone)
        labels = {c.label for c in res.result.candidates}
        for candidate in res.result.candidates:
            assert set(candidate.codon_distance_to) == labels - {candidate.label}
            for distance in candidate.codon_distance_to.values():
                assert 0.0 <= distance <= 1.0

    def test_the_shipped_defaults_are_the_ones_measured(self) -> None:
        """The conftest shrinks the sweep and the null so unit tests are cheap.
        This pins what the PRODUCT does, so that shrinking cannot drift into
        being the default."""
        assert DEFAULT_GALLERY_SIZE == 5
        assert MIN_GALLERY <= DEFAULT_GALLERY_SIZE <= MAX_GALLERY
        assert DEFAULT_SWEEP_STEPS == 1
        assert DEFAULT_NULL_N == 200


class TestPercentiles:
    def test_objectives_are_ranked_rather_than_reported_unavailable(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """The skeleton reported EVERY objective unavailable with the reason
        'ranking not computed'. That sentence must no longer be reachable."""
        res = fast(backbone)
        card = res.result.candidates[0].scorecard
        assert card.available, "no objective was ranked at all"
        for score in card.scores:
            assert "ranking not computed" not in score.unavailable_reason

    def test_every_percentile_is_a_fraction_of_a_real_null(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        res = fast(backbone)
        for score in res.result.candidates[0].scorecard.available:
            assert 0.0 <= score.percentile <= 1.0
            assert score.null_n >= 2
            assert not math.isnan(score.null_mean)
            assert score.null_kind in ("host_frequency", "uniform_synonymous")
            assert score.windowed_fold_only is True

    def test_an_unavailable_objective_keeps_its_reason(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """Named, never dropped. A scorecard missing its highest-weight objective
        looks exactly like one where that objective was never configured."""
        res = fast(backbone)
        for score in res.result.candidates[0].scorecard.unavailable:
            assert score.unavailable_reason
            assert math.isnan(score.percentile)

    def test_no_nan_reaches_a_reported_percentile(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """c1_cai returns NaN on this host (no CAI reference set for hek293 --
        the build ships only Sharp & Li's E. coli table). Every comparison
        against NaN is False, so an unguarded `percentile_of` would report it as
        a confident 0.0."""
        res = fast(backbone)
        for score in res.result.candidates[0].scorecard.available:
            assert not math.isnan(score.raw)
            assert not math.isnan(score.percentile)

    def test_the_null_is_built_on_the_assembled_construct(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """A null scored on bare CDSs measures against a distribution that never
        contained a backbone, and the report line is then simply false. The
        observable consequence: the null's spread reflects the construct the
        design sits in, so a null must exist for every ranked objective and carry
        the seed it was drawn with."""
        res = fast(backbone)
        assert res.nulls is not None
        for score in res.result.candidates[0].scorecard.available:
            null = res.nulls.by_spec[score.spec_id]
            assert null.n == score.null_n
            assert null.seed == res.inputs.seed

    def test_percentiles_are_comparable_across_the_panel(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """One null, shared. Ranking candidate B against a null anchored on B
        while candidate A is ranked against a null anchored on A would make the
        two percentiles incomparable -- and the gallery's whole job is to be
        compared."""
        res = fast(backbone)
        assert res.nulls is not None
        for candidate in res.result.candidates:
            for score in candidate.scorecard.available:
                assert score.null_n == res.nulls.by_spec[score.spec_id].n


class TestNativeBaseline:
    def test_it_stays_none_and_says_why_when_the_caller_supplies_nothing(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """`native_baseline` is a claim about a real wild-type CDS. BT5 will not
        manufacture one by back-translating, because 'do not optimize' compared
        against a sequence BT5 itself designed is not a comparison."""
        res = fast(backbone)
        assert res.result.native_baseline is None
        assert res.report.native_baseline_available is False
        assert "The native sequence is included as a candidate" not in res.rendered
        assert any("no native baseline" in d for d in res.result.provenance.degradations)

    def test_a_supplied_native_is_either_a_candidate_or_a_stated_refusal(
        self, backbone: VectorBackbone, fast: Any, protein: str, native_cds: str
    ) -> None:
        """Two legitimate outcomes and no third one. Either the native CDS is a
        first-class candidate a user can order, or the independent validator
        refused it and the report says which invariant -- BT5 does not offer a
        sequence it cannot prove, and it does not drop one silently either."""
        res = fast(backbone, native_cds=native_cds)
        baseline = res.result.native_baseline
        if baseline is None:
            assert any(
                "the native CDS was supplied" in d for d in res.result.provenance.degradations
            )
            assert res.report.native_baseline_available is False
            assert "The native sequence is included as a candidate" not in res.rendered
            return
        assert baseline.label == "native_baseline"
        assert baseline.cds == native_cds
        assert CODE.translate(baseline.cds)[:-1] == protein
        assert baseline.scorecard.scores, "the baseline is ranked like any candidate"
        assert res.report.native_baseline_available is True
        assert "The native sequence is included as a candidate" in res.rendered
        # "Do not optimize" is only a real option if the user can order the tube.
        assert any("native_baseline" in entry.name for entry in res.orders)
        assert baseline.design_hash in {entry.name.rsplit("_", 1)[-1] for entry in res.orders}

    def test_a_native_encoding_a_different_protein_is_refused(
        self, backbone: VectorBackbone, fast: Any, native_cds: str
    ) -> None:
        """And refused BEFORE the sweep runs. A caller error that costs 20 solves
        and a 200-variant null before it is reported is a caller error reported
        badly."""
        # Swap the second codon for one encoding a different residue.
        mutated = native_cds[:3] + "TGG" + native_cds[6:]
        with pytest.raises(DesignError, match="not a baseline"):
            fast(backbone, native_cds=mutated)

    def test_the_baseline_is_never_manufactured(self, backbone: VectorBackbone, fast: Any) -> None:
        """There is deliberately no parameter that asks BT5 to invent a native
        sequence, and no code path that back-translates one into the baseline
        slot. A manufactured baseline is a design wearing the word."""
        import inspect

        parameters = inspect.signature(design).parameters
        assert "native_cds" in parameters
        assert parameters["native_cds"].default is None
        assert not any("baseline" in name and name != "native_cds" for name in parameters)


class TestCompleteness:
    def test_is_complete_tracks_what_is_actually_missing(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """The skeleton hard-wired `is_complete` to False by always emitting
        three degradations. It is now derived, so the equivalence is the property
        worth pinning -- and it is what makes a True here mean something."""
        res = fast(backbone)
        expected = not res.report.unavailable and not res.report.degradations
        assert res.report.is_complete is expected

    def test_an_unscreened_run_can_never_be_complete(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """The screen is M8's to make real. Until it reports 'clear', the report
        must not claim completeness, and it must never print 'clear' for a screen
        that did not run."""
        res = fast(backbone)
        assert res.report.is_complete is False
        assert any(
            "biosecurity screening: not_run" in d for d in res.result.provenance.degradations
        )
        # The degradations are rendered, so the reader sees "not_run" and never
        # a "clear" this run did not earn.
        assert "biosecurity screening: not_run" in res.rendered
        assert "biosecurity screening: clear" not in res.rendered

    def test_a_clear_verdict_removes_only_the_screening_degradation(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """This lane RENDERS a verdict; it does not produce one. Handed a clear
        verdict it stops emitting the screening degradation and changes nothing
        else -- the other degradations are about other absences."""
        from bt5.core.context import BiosecurityVerdict

        unscreened = fast(backbone)
        screened = fast(backbone, screen=BiosecurityVerdict("clear", "test-db-1", ""))
        assert UNSCREENED.status == "not_run"
        removed = set(unscreened.result.provenance.degradations) - set(
            screened.result.provenance.degradations
        )
        assert all("biosecurity screening" in d for d in removed)
        assert removed

    def test_no_degradation_arrives_unremarked(self, backbone: VectorBackbone, fast: Any) -> None:
        """The skeleton pinned its degradations by set equality so a new silent
        one would fail. That set is environment-dependent now (ViennaRNA present
        or not, host tables shipped or not), so the same protection is expressed
        as a closed set of RECOGNISED sentences instead."""
        res = fast(backbone)
        for degradation in res.result.provenance.degradations:
            assert degradation.startswith(KNOWN_DEGRADATIONS), (
                f"unrecognised degradation: {degradation!r}. If this is a real new "
                f"source, add it to KNOWN_DEGRADATIONS deliberately."
            )


class TestOrderFile:
    def test_a_csv_is_emitted_with_one_row_per_candidate(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        res = fast(backbone)
        rows = list(csv.reader(io.StringIO(res.order_csv)))
        assert rows[0] == ["Name", "Sequence"]
        assert len(rows) - 1 == len(res.orders) == len(res.result.candidates)

    def test_every_ordered_sequence_is_the_cds_in_bare_acgt(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """A lowercase base or a stray IUPAC N reaches the synthesiser as an
        order, not as a question."""
        res = fast(backbone)
        ordered = {entry.sequence for entry in res.orders}
        assert ordered == {c.cds for c in res.result.candidates}
        for sequence in ordered:
            assert set(sequence) <= DNA_ALPHABET

    def test_every_tube_label_carries_its_design_hash(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """Two different sequences under one name is how a lab ends up with two
        tubes and an irreproducible result."""
        res = fast(backbone)
        names = [entry.name for entry in res.orders]
        assert len(set(names)) == len(names)
        for candidate in res.result.candidates:
            assert any(name.endswith(candidate.design_hash) for name in names)


class TestNoPredictionVocabulary:
    def test_the_report_states_what_a_percentile_is_not(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        res = fast(backbone)
        assert "never a predicted expression level" in res.rendered
        assert "percentile against a random-synonymous null; not a prediction" in res.rendered

    def test_no_score_field_carries_an_absolute_claim(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """`ObjectiveScore` has no field for a predicted level and must not grow
        one here; the percentile plus its null is the whole vocabulary."""
        res = fast(backbone)
        fields = set(ObjectiveScore.__dataclass_fields__)
        banned = {"expression", "titer", "yield", "fold_improvement", "predicted"}
        assert not (fields & banned)

        # A bare word scan over the rendering cannot express this property: the
        # DISCLAIMER contains "predicted expression level" and must, because
        # saying what BT5 refuses to report is the point. What must hold is that
        # the phrase appears ONLY there -- never on a line that reports a number.
        disclaimer = "BT5 reports ranks and percentiles"
        for line in res.rendered.splitlines():
            lowered = line.lower()
            if line.strip().startswith(disclaimer) or "explain 5-31%" in lowered:
                continue
            for phrase in ("predicted expression", "fold-improvement", "titer of"):
                assert phrase not in lowered, f"prediction vocabulary on a report line: {line!r}"


class TestSweepAxes:
    """The dead-axis claim, tested rather than asserted in a docstring."""

    def test_dropping_a_dead_axis_loses_no_design(
        self, backbone: VectorBackbone, protein: str
    ) -> None:
        """`live_axes` drops `codon_adaptation` when no host codon-usage table is
        on file, because its cost term is then identically zero. That is a claim
        about the reachable front, so it is measured: sweeping all four axes and
        sweeping the live ones must produce the SAME set of designs, and the
        dead axis must cost a full `optimize()` to prove it.
        """
        space = _solve_space(backbone, protein)
        assert live_axes({}) == ("repeat_avoidance", "gc_lean_at", "gc_lean_gc")
        assert live_axes({"CTG": 1.0}) == SWEEP_AXES

        everything, _ = sweep_designs(space, axes=SWEEP_AXES, steps=1, k=5)
        live, _ = sweep_designs(space, steps=1, k=5)
        assert {p.cds for p in everything.picks} == {p.cds for p in live.picks}
        assert live.swept == everything.swept - 1  # one solve saved, nothing lost

    def test_a_host_with_a_usage_table_keeps_the_adaptation_axis(self) -> None:
        """The axis is dropped because it is DEAD here, not because it is
        unwanted. When S6 ships a host table it comes back on its own."""
        assert "codon_adaptation" in live_axes({"CTG": 1.0, "TTA": 0.1})
