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
import re
from pathlib import Path
from typing import Any

import pytest
from bt5.codon.tables import FileTableProvider
from bt5.core.context import BiosecurityVerdict, HostId, Modality
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


def _scored_objectives(backbone: VectorBackbone) -> list[Any]:
    """The objectives `design()` scores for the `fast` fixture's context.

    Rebuilt rather than plumbed out of `SkeletonResult`, because a test that
    reads the runner's own list back cannot notice the runner passing the wrong
    one.
    """
    from bt5.design.catalog import scored_objectives
    from bt5.design.runner import _context
    from bt5.design.sites import choose_site
    from bt5.solver.catalog import build_rule_set, default_services

    site = choose_site(backbone, table_id=1)
    ctx = _context(
        modality=Modality.LENTIVIRAL,
        hosts=[HostId.HEK293],
        table_id=1,
        cassette_orientation=site.strand,
    )
    return list(scored_objectives(build_rule_set(ctx, default_services(seed=0))))


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


#: Every degradation `design()` is allowed to emit, as ONE keyed table:
#: pattern -> (a sentence it must match, a literal fragment its source must
#: contain).
#:
#: Keyed, and all three in one place, because the previous shape let a rewording
#: slip past every guard. The pattern lived in one tuple, its sample in a second,
#: and the source fragment in a third POSITIONAL tuple checked only for length.
#: Rewording the short-panel sentence updated the source; the pattern and its
#: sample drifted together, and the fragment anchored the one clause that had not
#: changed. Three guards, all green, one broken degradation.
#:
#: The chain below closes that: the pattern must match the sample, the fragment
#: must be a substring OF the sample, and the fragment must appear literally in
#: the module that emits it. A reworded sentence therefore breaks the fragment
#: check, and a fragment updated without the pattern breaks the substring check.
#: Each fragment must also cover more than a leading phrase -- the drift got
#: through on `"the sweep produced "`, a prefix that survived the rewording.
DEGRADATION_SOURCES: dict[str, tuple[str, str]] = {
    r"ViennaRNA .*": (
        "ViennaRNA is not installed, so folding objectives were not evaluated; no",
        "ViennaRNA is not installed, so folding objectives were not evaluated",
    ),
    r"protein-level biosecurity screening: \w+": (
        "protein-level biosecurity screening: not_run (did not run)",
        "protein-level biosecurity screening: ",
    ),
    r"the sweep produced \d+ distinct designs? from \d+ weight vectors? that solved, "
    r"below the \d+": (
        "the sweep produced 1 distinct design from 3 weight vectors that solved, "
        "below the 3 at which a panel offers a genuine choice.",
        " that solved, below the ",
    ),
    r"the \d+-candidate panel does not meet gate G4:": (
        "the 3-candidate panel does not meet gate G4: its minimum pairwise",
        "-candidate panel does not meet gate G4: its minimum ",
    ),
    r"single candidate only:": (
        "single candidate only: no weight vector in the sweep produced a design",
        "single candidate only: no weight vector in the sweep produced a design",
    ),
    r"no codon usage table on file": (
        "no codon usage table on file for host hek293; the null was sampled",
        "no codon usage table on file for host ",
    ),
    r"no native baseline:": (
        "no native baseline: the caller supplied no wild-type CDS, and BT5 will",
        "no native baseline: the caller supplied no wild-type CDS",
    ),
    r"the native CDS was supplied": (
        "the native CDS was supplied but is not offered as a candidate: the",
        "the native CDS was supplied but is not offered as a candidate",
    ),
    r"objective \S+ not ranked:": (
        "objective c1_cai not ranked: no CAI reference set for host hek293",
        " not ranked: ",
    ),
    r"rule \S+ not run:": (
        "rule b1_five_prime not run: its thresholds are calibrated against a folding "
        "engine that is not available",
        " not run: its thresholds are calibrated against a folding ",
    ),
    r"forbidden motif \S+ carried by the backbone": (
        "forbidden motif TCTAGA carried by the backbone, excluded from enforcement",
        " carried by the backbone, excluded from enforcement",
    ),
    r"screening burden unavailable:": (
        "screening burden unavailable: no published error-free length on file for x",
        "screening burden unavailable: no published error-free length on file ",
    ),
}

KNOWN_DEGRADATIONS = tuple(DEGRADATION_SOURCES)


def _assert_recognised(degradation: str) -> None:
    """A degradation a test PRODUCED must match a pattern the guard knows.

    This is the check that survives a rewording: the sentence comes from the
    source, not from a hand-written sample, so the two cannot drift together.
    Three of the four `_panel` outcomes never occur on the reference fixture, so
    without this their patterns are only ever checked against samples.
    """
    assert any(re.match(pattern, degradation) for pattern in KNOWN_DEGRADATIONS), (
        f"a real degradation matches no known pattern: {degradation!r}"
    )


#: The modules that emit a degradation sentence.
DEGRADATION_SOURCE_FILES = (
    "design/runner.py",
    "score/report.py",
    "structure/vienna.py",
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
        """Ranked on `comparable_totals`, NOT on `ScoreCard.total`.

        Asserting `total` descending would be asserting the invariant
        `comparable_totals` deliberately breaks: when two candidates have
        different available-objective sets their `total`s are not comparable, and
        a suite pinned to `total` would report that fix as a regression. The
        labels carry the order, so they are what to check.
        """
        from bt5.design.ranking import comparable_totals

        res = fast(backbone)
        candidates = res.result.candidates
        assert [c.label for c in candidates] == [f"design_{i + 1}" for i in range(len(candidates))]
        # The REAL objectives, which is what the runner passes. Handing
        # `comparable_totals` an empty spec list makes every weight 0.0 and every
        # key 0.0 -- an assertion that then holds for ANY ordering, including a
        # fully inverted one. That is how the first version of this test passed
        # while covering nothing.
        objectives = _scored_objectives(backbone)
        assert objectives, "no scored objective, so ranking cannot be tested"
        keys = comparable_totals({c.label: c.scorecard for c in candidates}, objectives)
        assert len(set(keys.values())) > 1, (
            "every candidate scored identically, so this run cannot tell a correct "
            "sort from an inverted one and the assertion below would be vacuous"
        )
        ordered = [keys[c.label] for c in candidates]
        assert ordered == sorted(ordered, reverse=True)

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
        three degradations. It is now DERIVED, so the property worth testing is
        that filling an absence actually removes its degradation and moves the
        report toward complete -- not that `is_complete` equals its own body,
        which is a tautology that cannot fail.
        """
        res = fast(backbone)
        supplied = fast(backbone, screen=BiosecurityVerdict("clear", "test-db-1", ""))
        removed = set(res.report.degradations) - set(supplied.report.degradations)
        assert removed, "filling an absence removed no degradation"
        assert len(supplied.report.degradations) < len(res.report.degradations)
        # Still False, and for reasons that are real rather than hard-wired: no
        # host codon-usage table ships for hek293, and no wild-type CDS was given.
        assert supplied.report.is_complete is False
        assert supplied.report.degradations

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
        # BOTH lists. `provenance.degradations` is what `design()` assembles;
        # `build_report` appends its own (the screening-burden line) to the
        # report's copy, so scanning only provenance leaves a whole source
        # uncovered -- in the test named for catching uncovered sources.
        seen = set(res.result.provenance.degradations) | set(res.report.degradations)
        assert set(res.report.degradations) >= set(res.result.provenance.degradations)
        for degradation in seen:
            assert any(re.match(pattern, degradation) for pattern in KNOWN_DEGRADATIONS), (
                f"unrecognised degradation: {degradation!r}. If this is a real new "
                f"source, add it to KNOWN_DEGRADATIONS deliberately."
            )

    def test_each_pattern_matches_its_own_sample(self) -> None:
        for pattern, (sample, _fragment) in DEGRADATION_SOURCES.items():
            assert re.match(pattern, sample), f"{pattern!r} does not match {sample!r}"

    def test_each_fragment_anchors_its_pattern_to_a_real_source_sentence(self) -> None:
        """The link that the previous shape was missing.

        `fragment in sample` ties the fragment to the pattern; `fragment in
        corpus` ties it to the code. Reword the sentence in the source and the
        second fails; update the fragment without the pattern and the first
        does.

        This is the SECONDARY net. It is a heuristic -- a fragment can only be
        as long as the longest literal run in an f-string, which for
        "objective {id} not ranked: " is thirteen characters -- so it cannot
        guarantee the fragment covers the clause someone rewords. The primary
        protection is that every sentence a test can actually PRODUCE is matched
        against these patterns at the point of production: the emitted set in
        `test_no_degradation_arrives_unremarked`, and the three panel sentences
        in `TestPanelDegradations` / `TestTheFallbackWhenNothingSolves`. That is
        what would have caught the rewording that got past the previous shape.
        """
        root = Path(__file__).resolve().parents[4] / "packages/engine/src/bt5"
        corpus = "".join((root / name).read_text() for name in DEGRADATION_SOURCE_FILES)
        for pattern, (sample, fragment) in DEGRADATION_SOURCES.items():
            assert len(fragment) >= 12, f"{pattern!r}: fragment {fragment!r} is too short"
            assert fragment in sample, f"{pattern!r}: {fragment!r} is not in its own sample"
            assert fragment in corpus, (
                f"{pattern!r}: no source in {DEGRADATION_SOURCE_FILES} emits "
                f"{fragment!r} any more -- the sentence was reworded and this "
                f"pattern is now stale."
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


class TestTheNullNeverFoldsWholeTranscripts:
    """`FoldEngine.mfe` is report-time only, and the null must honour that.

    Its own docstring says ~0.24 s at 1 kb and ~6.5 s at 3 kb, "so this must
    never run inside the interactive loop or the empirical null" — a 200-variant
    null calling it would be minutes. `build_nulls` hard-codes
    `windowed_fold_only=True`, and `null_distribution` refuses a default for that
    flag precisely because "a caller that has not thought about whether its
    scorer folds whole transcripts has not earned a percentile".

    **This is a forward-looking guard, and on this context it is vacuous today.**
    No scored objective touches the fold engine at all: `b1_five_prime` is the
    only rule that folds, it uses `mfe_window` rather than `mfe`, and it is
    gated OUT of a HEK293 producer slot. So the test proves the flag is safe
    rather than exercising a rule that could break it. That is worth having —
    it fails the moment a scored rule starts folding — but the vacuity is
    asserted below rather than left for a reader to discover, because a guard
    that is silently vacuous is indistinguishable from one that is passing.
    """

    def _counting_fold(self) -> tuple[object, dict[str, int]]:
        from bt5.solver.catalog import default_services

        inner = default_services(seed=0).fold
        if inner is None:
            pytest.fail(
                "ViennaRNA is not installed, so this test cannot prove anything. "
                "/bootstrap installs the [fold] extra."
            )
        calls: dict[str, int] = dict.fromkeys(("mfe", "mfe_window", "accessibility", "duplex"), 0)

        class Counting:
            """Wraps rather than replaces, so `check_engine_calibration` still
            sees the calibrated engine identity and no rule is dropped for the
            wrong reason — the point is to catch a rule that folds, not to hide
            one."""

            name = inner.name
            version = inner.version
            param_set = inner.param_set

            def mfe(self, seq: str) -> object:
                calls["mfe"] += 1
                raise AssertionError(
                    "a scored objective called FoldEngine.mfe, which is reserved for "
                    "report time. Inside the null that is ~0.24 s per variant at 1 kb; "
                    "at DEFAULT_NULL_N, minutes. Either the rule must use mfe_window, "
                    "or build_nulls must stop claiming windowed_fold_only."
                )

            def mfe_window(self, seq: str, iv: object) -> object:
                calls["mfe_window"] += 1
                return inner.mfe_window(seq, iv)

            def accessibility(self, seq: str, iv: object, u: int) -> object:
                calls["accessibility"] += 1
                return inner.accessibility(seq, iv, u)

            def duplex(self, a: str, b: str) -> object:
                calls["duplex"] += 1
                return inner.duplex(a, b)

        return Counting(), calls

    def test_no_scored_objective_calls_mfe(self, backbone: VectorBackbone, protein: str) -> None:
        from bt5.design.catalog import scored_objectives
        from bt5.design.runner import _context
        from bt5.design.sites import choose_site
        from bt5.solver.catalog import build_rule_set, default_services
        from bt5.solver.reference import back_translate
        from bt5.vector import assemble

        engine, calls = self._counting_fold()
        services = default_services(seed=0, fold=engine)  # type: ignore[arg-type]
        site = choose_site(backbone, table_id=1)
        ctx = _context(
            modality=Modality.LENTIVIRAL,
            hosts=[HostId.HEK293],
            table_id=1,
            cassette_orientation=site.strand,
        )
        rules = build_rule_set(ctx, services)
        objectives = scored_objectives(rules)
        assert objectives, "no scored objective to check"

        # Any valid CDS serves; a full design() run would add seconds and prove
        # nothing extra, since what is under test is what the RULES call.
        construct = assemble(
            backbone,
            back_translate(protein, CODE),
            protein=protein,
            table_id=1,
            site=site,
        ).construct
        for spec in objectives:
            spec.evaluate(construct, ctx, services)

        assert calls["mfe"] == 0  # the claim

        # And the vacuity, stated. If either of these changes the test above
        # starts doing real work and this assertion is what tells you.
        assert "b1_five_prime" not in {spec.id for spec in objectives}
        assert calls == dict.fromkeys(calls, 0), (
            f"a scored objective now touches the fold engine ({calls}); this test is "
            f"no longer vacuous, which is good — update the docstring."
        )


class _StubSpace:
    """The two attributes `_panel` and `sweep_designs` actually read.

    A real `SolveSpace` would solve, repair and verify for every weight vector,
    which is seconds per case and measures the solver rather than the branch
    under test. What these tests are about is which SENTENCE comes back for
    which shape of panel.
    """

    def __init__(self, sequences: list[str]) -> None:
        self._sequences = sequences
        self._calls = 0
        self.usage: dict[str, float] = {}

    def solve(self, weights: object) -> object:
        from bt5.solver.pipeline import OptimizeResult
        from bt5.solver.repair import RepairOutcome

        if not self._sequences:
            return None
        cds = self._sequences[min(self._calls, len(self._sequences) - 1)]
        self._calls += 1
        return OptimizeResult(
            cds=cds,
            construct=None,  # type: ignore[arg-type]
            repair_outcome=RepairOutcome(cds, 0, True, stop_reason="clean"),
        )


def _cds(*codons: str) -> str:
    """A CDS of ten codons, differing from its siblings at the given positions."""
    base = ["ATG"] + ["AAA"] * 9
    for i, codon in enumerate(codons):
        base[i + 1] = codon
    return "".join(base)


class TestPanelDegradations:
    """Which sentence comes back for which shape of panel.

    `Gallery.meets_g4` is False for two different reasons — too FEW candidates
    (`len(picks) < MIN_GALLERY`) or candidates too CLOSE
    (`min_pairwise_distance < 0.15`) — and `pairwise_minimum` returns 1.0 for a
    single sequence. Reporting the distance sentence for the count failure emits
    a literally false claim, in the one vocabulary this lane exists to keep
    honest.
    """

    def test_a_short_panel_makes_no_distance_claim(self) -> None:
        """The regression. Before the fix this said "its minimum pairwise codon
        distance is 100.0%, below the 15%" — of a one-candidate panel."""
        from bt5.design.runner import _panel

        gallery, picks, _solved, degradation = _panel(
            _StubSpace([_cds()]),  # type: ignore[arg-type]
            steps=1,
            k=5,
        )
        assert len(picks) == 1
        assert gallery is not None
        assert not gallery.meets_g4
        assert degradation is not None
        # The property is that no NUMBER is asserted as a distance. The sentence
        # may say a distance is not reported; it may not report one.
        assert "%" not in degradation
        assert f"below the {MIN_GALLERY}" in degradation
        _assert_recognised(degradation)

    def test_a_close_panel_reports_the_distance_it_reached(self) -> None:
        from bt5.design.runner import _panel

        # Three designs differing at one codon in ten -> 10% pairwise, under 15%.
        _gallery, picks, _solved, degradation = _panel(
            _StubSpace([_cds("GCA"), _cds("GCC"), _cds("GCG")]),  # type: ignore[arg-type]
            steps=1,
            k=5,
        )
        assert len(picks) == MIN_GALLERY
        assert degradation is not None
        assert "does not meet gate G4" in degradation
        assert "10.0%" in degradation
        _assert_recognised(degradation)

    def test_a_complete_short_panel_is_not_a_degradation(self) -> None:
        """`k` is a CEILING. A sweep that exhausted the front at three genuinely
        different designs has answered completely, and degrading it pinned
        `is_complete` to False on every run — the skeleton's hard-wired False
        wearing a different sentence."""
        from bt5.design.runner import _panel

        far = [
            _cds("GCA", "GCA", "GCA", "GCA"),
            _cds("TGC", "TGC", "TGC", "TGC"),
            _cds("GAT", "GAT", "GAT", "GAT"),
        ]
        gallery, picks, _solved, degradation = _panel(
            _StubSpace(far),  # type: ignore[arg-type]
            steps=1,
            k=5,
        )
        assert len(picks) == 3 < 5
        assert gallery is not None
        assert gallery.meets_g4
        assert degradation is None, f"a complete panel must not degrade: {degradation!r}"

    def test_the_shipped_lattice_can_reach_min_gallery(self) -> None:
        """The arithmetic behind the fix: with `k` as a ceiling the defaults must
        still be able to fill a real panel. `steps=1` gives one vector per axis,
        so 3 without a host usage table and 4 with one — both at or above
        `MIN_GALLERY`, and both below `DEFAULT_GALLERY_SIZE`, which is why `k`
        had to become a ceiling rather than a promise."""
        from bt5.score.gallery import simplex_weights

        for axes in (live_axes({}), SWEEP_AXES):
            vectors = len(simplex_weights(axes, DEFAULT_SWEEP_STEPS))
            assert vectors >= MIN_GALLERY
            assert vectors <= DEFAULT_GALLERY_SIZE


class TestTheFallbackWhenNothingSolves:
    """The branch `_panel` checks first and nothing covered.

    `_StubSpace` above cannot express it: it returns the same thing for a swept
    vector and for the unsteered `solve(None)`. This stub separates them, which
    is the whole distinction the fallback rests on.
    """

    class _NothingSweeps:
        def __init__(self, fallback: str | None) -> None:
            self._fallback = fallback
            self.usage: dict[str, float] = {}

        def solve(self, weights: object) -> object:
            from bt5.solver.pipeline import OptimizeResult
            from bt5.solver.repair import RepairOutcome

            if weights is not None or self._fallback is None:
                return None  # every swept vector infeasible or refused
            return OptimizeResult(
                cds=self._fallback,
                construct=None,  # type: ignore[arg-type]
                repair_outcome=RepairOutcome(self._fallback, 0, True, stop_reason="clean"),
            )

    def test_the_unsteered_solve_stands_alone_and_says_there_is_no_gallery(self) -> None:
        from bt5.design.runner import _panel

        gallery, picks, solved, degradation = _panel(
            self._NothingSweeps(_cds()),  # type: ignore[arg-type]
            steps=1,
            k=5,
        )
        assert gallery is None, "there is no panel, so there is no Gallery to report"
        assert picks == [_cds()]
        assert solved[_cds()] is not None
        assert degradation is not None
        assert degradation.startswith("single candidate only:")
        assert "%" not in degradation
        _assert_recognised(degradation)

    def test_nothing_at_all_raises_rather_than_returning_a_partial_result(self) -> None:
        """Refusing is the guarantee working. Returning an empty panel would let
        a caller index `candidates[0]` on nothing."""
        from bt5.design.runner import _panel

        with pytest.raises(DesignError, match="no candidate survived"):
            _panel(self._NothingSweeps(None), steps=1, k=5)  # type: ignore[arg-type]


class TestCandidateProvenance:
    """Where a caller finds a design's provenance, pinned in both directions (#99).

    #89 restructured the walking skeleton's single annotated result into a panel,
    and `annotate()` became a per-EXPORT step: it is called once, on the winner,
    and its output feeds the GenBank string only. Nothing attaches it back to a
    `Candidate`, so `candidate.construct` is the raw assembly for every candidate
    including the baseline. A review bot on `design/runner.py:567` read that as a
    provenance loss, and #99 asks the lane owner to decide whether it is one.

    Nothing tested it either way -- `.sequence` was compared and `.features` and
    `.annotations` never were -- so the decision had nothing to overturn and no
    way to notice a drift. This class is that missing pin. It does NOT decide
    #99; it states exactly what is true today so the choice is between two
    described states rather than two guesses.

    IF #99 DECIDES THE CURRENT SHAPE IS INTENDED, this class stays and the
    `Candidate` docstring in `core/result.py` should say annotation is an export
    step (a `core/` docstring, so via `/contract-change`).

    IF #99 DECIDES IT IS AN OVERSIGHT, `test_candidate_constructs_are_unannotated`
    is the test to invert -- and it names, in its own message, what the fix has to
    supply. Do not delete it; make it assert the annotated shape instead.
    """

    def test_design_hash_is_on_the_candidate_itself(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """The load-bearing half, and the reason the current shape may be fine.

        `Candidate.design_hash` is a first-class FIELD. A caller wanting
        provenance reads it directly and never has to parse a GenBank comment, so
        an unannotated `construct` costs nothing as long as this holds. If this
        test ever fails, the #99 question stops being a style call and becomes a
        real loss.
        """
        res = fast(backbone)
        for candidate in res.result.candidates:
            assert candidate.design_hash, f"{candidate.label} has no design_hash"
        hashes = {c.design_hash for c in res.result.candidates}
        assert len(hashes) == len(res.result.candidates), (
            "two candidates share a design_hash; the hash is the tube label and "
            "two tubes under one name is the failure core/result.py names"
        )

    def test_the_hash_reaches_the_report_and_the_genbank(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """What `core/result.py`'s Candidate docstring actually promises: "the
        content hash travels onto the report, the GenBank note and the order
        file". It says nothing about `Candidate.construct`, which is why the
        current shape is defensible -- the promise is about where the hash
        ARRIVES, and it arrives."""
        res = fast(backbone)
        winner = res.result.candidates[0]
        assert winner.design_hash in res.genbank
        assert winner.design_hash in res.rendered
        assert any(winner.design_hash in entry.name for entry in res.orders)

    def test_candidate_constructs_are_unannotated(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """THE #99 PIN. Every candidate carries the raw assembly, not the
        annotated one: no `bt5_origin` provenance stamps and no design comment.

        Inverting this is the fix if #99 decides annotation belongs on the
        candidate. What that fix must supply is exactly what this test looks for.
        """
        from bt5.vector.annotate import ORIGIN_QUALIFIER

        res = fast(backbone)
        panel = [*res.result.candidates]
        if res.result.native_baseline is not None:
            panel.append(res.result.native_baseline)
        for candidate in panel:
            construct = candidate.construct
            stamped = [f for f in construct.features if ORIGIN_QUALIFIER in f.qualifiers]
            assert not stamped, (
                f"{candidate.label}: construct now carries {len(stamped)} "
                f"{ORIGIN_QUALIFIER} stamps. If #99 was decided in favour of "
                f"annotating candidates, invert this test rather than deleting it."
            )
            assert "comment" not in (construct.annotations or {}), (
                f"{candidate.label}: construct now carries a provenance comment; "
                f"see #99 and invert this test."
            )

    def test_the_exported_construct_is_annotated_and_the_candidate_is_not(
        self, backbone: VectorBackbone, fast: Any
    ) -> None:
        """The two shapes side by side, which is what makes the #99 choice legible.

        Same sequence, different annotation: whatever is missing from
        `candidate.construct` is present in the export, so nothing is LOST -- it
        is only reachable from one object rather than two.
        """
        import io as _io

        from bt5.vector import read_genbank
        from bt5.vector.annotate import ORIGIN_QUALIFIER

        res = fast(backbone)
        winner = res.result.candidates[0]
        exported = read_genbank(_io.StringIO(res.genbank))

        assert exported.sequence == winner.construct.sequence, (
            "the export is a different sequence from the candidate it came from"
        )
        assert [f for f in exported.features if ORIGIN_QUALIFIER in f.qualifiers], (
            "the export carries no provenance stamps either; then the #99 question "
            "is not where annotation lives but whether it happens at all"
        )
        assert len(exported.features) >= len(winner.construct.features)
