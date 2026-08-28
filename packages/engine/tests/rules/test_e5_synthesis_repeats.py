"""E5: the repeats the VENDOR cares about, and the scope that makes it not-F1.

Two rules in this catalog scan for exact direct repeats. If E5 were F1 with
different constants the panel would report every repeat twice, so most of what
matters here is the boundary: E5 sees a linear ordered fragment plus vendor
adapters, F1 sees the assembled circular plasmid, and they are meant to
disagree at the origin, at the CDS/backbone junction, and on short GC-rich
repeats.

`slot` and friends are imported at module level, not inside a test. There are
four `conftest.py` files under packages/engine/tests and none is a package, so a
deferred `from conftest import ...` resolves against whichever one pytest bound
last.
"""

from __future__ import annotations

import numpy as np
import pytest
from bt5.core.context import Modality
from bt5.core.registry import discover, get
from bt5.core.services import Services
from bt5.core.spec import Enforcement
from bt5.core.types import Construct, Interval, Segment, SegmentKind, Topology
from bt5.rules.catalog.e5_synthesis_repeats import (
    ANNEAL_C,
    MIN_LENGTH_BP,
    SynthesisRepeats,
    duplex_tm,
    severity,
)
from bt5.rules.catalog.f1_direct_repeats import DirectRepeats
from bt5.rules.fragment import TWIST_FIVE_PRIME
from bt5.vector.kmers import ConstructKmerIndex
from conftest import construct, context, slot

discover()

#: Aperiodic, 14 bp, Tm 62.5 C. Under every published length limit and over the
#: annealing temperature: the case a length rule cannot see.
SHORT_HOT = "GCCGCGGCGCGCGT"
#: 22 bp, Tm 52.4 C. Over the 20 bp vendor limit and under the anneal: the case
#: the Tm criterion alone would miss.
LONG_COOL = "CGGGTAGCCAACTACTTAAGAC"


def dna(n: int, seed: int = 3) -> str:
    rng = np.random.default_rng(seed)
    return "".join("ACGT"[i] for i in rng.integers(0, 4, n))


@pytest.fixture
def svc() -> Services:
    return Services(
        fold=None,
        kmer=ConstructKmerIndex,
        tables=None,  # type: ignore[arg-type]
        rng=np.random.default_rng(1),
    )


def planted(unit: str, *, spacer_bp: int = 300) -> Construct:
    """Two copies of `unit` inside the CDS, with backbone after it."""
    cds = dna(200) + unit + dna(spacer_bp, 5) + unit + dna(200, 7)
    return construct(cds, dna(300, 11))


def hits(c: Construct, svc: Services, **kw: object):
    return SynthesisRepeats(**kw).evaluate(c, context(), svc).breaches  # type: ignore[arg-type]


class TestTheTmCriterion:
    """Length is what vendors publish; duplex stability is the mechanism."""

    def test_a_short_gc_rich_repeat_is_hard_and_invisible_to_length(self, svc: Services) -> None:
        c = planted(SHORT_HOT)
        found = hits(c, svc)
        assert found
        assert found[0].detail["severity"] == "hard"
        assert found[0].detail["duplex_tm_c"] >= ANNEAL_C
        assert "mis-prime" in found[0].message
        assert not DirectRepeats().evaluate(c, context(), svc).breaches, (
            "14 bp is under F1's floor -- if F1 caught it, E5's Tm branch would be "
            "reporting an already-reported finding rather than adding one"
        )

    def test_the_same_length_at_rich_repeat_is_only_a_warning(self, svc: Services) -> None:
        """The discrimination the whole criterion exists for: identical length,
        opposite verdict, because only one of them is a duplex at 60 C."""
        found = hits(planted("ATTATTAATTATAA"), svc)
        assert found
        assert found[0].detail["severity"] == "warn"
        assert found[0].detail["duplex_tm_c"] < ANNEAL_C

    def test_a_long_cool_repeat_is_hard_on_length_alone(self, svc: Services) -> None:
        """The converse: the Tm criterion ADDS to the published length rule, it
        does not replace it."""
        found = hits(planted(LONG_COOL), svc)
        assert found
        assert found[0].detail["duplex_tm_c"] < ANNEAL_C
        assert found[0].detail["severity"] == "hard"

    def test_either_criterion_alone_is_sufficient(self) -> None:
        assert severity(14, 65.0, hard_len=20, anneal_c=60.0) == "hard"
        assert severity(24, 40.0, hard_len=20, anneal_c=60.0) == "hard"
        assert severity(14, 40.0, hard_len=20, anneal_c=60.0) == "warn"
        assert severity(240, 40.0, hard_len=20, anneal_c=60.0) == "severe"

    def test_the_conditions_travel_with_every_tm(self, svc: Services) -> None:
        """A Tm without its salt and oligo concentration is not a number, for
        the same reason a dG without its energy parameters is not one."""
        found = hits(planted(SHORT_HOT), svc)
        assert "Na+ 50 mM" in str(found[0].detail["tm_conditions"])

    def test_the_pinned_tm_reproduces_the_documented_values(self) -> None:
        """The docstring's worked example. If the library's defaults move under
        us, this fails rather than silently reclassifying every repeat."""
        assert duplex_tm("GCGCGCGCGCGCGC") == pytest.approx(65.2, abs=0.2)
        assert duplex_tm("ACGTACGTACGTACGTACGT") == pytest.approx(53.1, abs=0.2)


class TestScopeIsTheFragment:
    """E5 scans what the vendor builds. F1 scans the assembled plasmid."""

    def test_a_repeat_shared_with_the_backbone_is_f1s_not_the_vendors(self, svc: Services) -> None:
        """The vendor never receives the backbone, so this cannot break their
        assembly PCR. It can still delete the plasmid, which is why F1 keeps it."""
        unit = dna(30, 9)
        c = construct(dna(200, 13) + unit, unit + dna(300, 17))
        assert not hits(c, svc)
        assert DirectRepeats().evaluate(c, context(), svc).breaches

    def test_a_repeat_spanning_the_origin_is_f1s_not_the_vendors(self, svc: Services) -> None:
        """The ordered fragment is linear; its two ends are never in one tube.

        The wrapping copy is split 10/10, so each linear piece is under E5's
        floor and only the circular scan can see the pair.

        Adapters are off here on purpose: this construct's 3' end plus the first
        two bases of the Twist adapter happen to reproduce a 12-mer from inside
        the insert, which is a real E5 finding and an irrelevant one to the
        question this test asks.
        """
        unit = dna(20, 11)
        seq = unit[10:] + dna(400) + unit + dna(400, 13) + unit[:10]
        c = Construct(
            seq,
            Topology.CIRCULAR,
            (Segment(Interval(0, len(seq)), SegmentKind.DESIGNABLE_CDS, "cds"),),
        )
        assert not hits(c, svc, vendor="none")
        assert DirectRepeats().evaluate(c, context(), svc).breaches

    def test_findings_carry_parent_construct_coordinates(self, svc: Services) -> None:
        """The fragment has its own coordinate space; a breach reported in it
        would point at the wrong bases in the report and the GenBank."""
        cds_start = 120
        unit = dna(30, 9)
        cds = dna(100) + unit + dna(200, 5) + unit + dna(100, 7)
        c = construct(cds, dna(300, 11), cds_start=cds_start)
        found = hits(c, svc)
        assert found
        iv = found[0].interval
        assert c.overlaps_editable(iv)
        assert iv.start >= cds_start, "a fragment offset was reported as a construct offset"

    def test_a_repeat_between_two_separate_fragments_is_not_a_synthesis_finding(
        self, svc: Services
    ) -> None:
        """Two designable spans are two tubes. The molecules never meet until
        after assembly, so neither vendor reaction can mis-prime on the other."""
        unit = dna(30, 9)
        seq = unit + dna(300, 5) + unit + dna(200, 7)
        c = Construct(
            seq,
            Topology.CIRCULAR,
            (
                Segment(Interval(0, len(unit)), SegmentKind.DESIGNABLE_CDS, "a"),
                Segment(
                    Interval(len(unit) + 300, len(unit) * 2 + 300),
                    SegmentKind.DESIGNABLE_CDS,
                    "b",
                ),
                Segment(Interval(len(unit), len(unit) + 300), SegmentKind.BACKBONE, "v"),
            ),
        )
        assert not hits(c, svc)
        assert DirectRepeats().evaluate(c, context(), svc).breaches


ADAPTER_ON = "twist_gene_fragment_adapter_on"


class TestAdapters:
    """Adapter sequence is synthesised with the fragment and is in no plasmid.

    Only on the ADAPTER-ON option, though. Twist states that adapter sequences
    are not added by default to Gene Fragments -- adapter-on and adapter-free
    are two choices made at checkout -- so every test here names the option it
    means rather than relying on the default.
    """

    def test_an_insert_reproducing_the_adapter_is_caught(self, svc: Services) -> None:
        c = construct(dna(150) + TWIST_FIVE_PRIME + dna(150, 5), dna(300, 11))
        found = hits(c, svc, vendor=ADAPTER_ON)
        assert found, "no whole-plasmid scan can see this: the adapter is not in the plasmid"
        assert any(b.detail["involves_adapter"] == "yes" for b in found)
        assert any("recode the insert side" in b.message for b in found)

    def test_the_default_order_carries_no_adapters(self, svc: Services) -> None:
        """The regression guard for the mixed-default bug.

        A plain Twist Gene Fragment is the ordered DNA and nothing else. When
        the default silently meant adapter-on, this construct produced a finding
        about a collision with sequence the user would never receive -- a false
        positive that costs real sequence freedom to repair.
        """
        c = construct(dna(150) + TWIST_FIVE_PRIME + dna(150, 5), dna(300, 11))
        assert not hits(c, svc)
        assert hits(c, svc, vendor=ADAPTER_ON), "and adapter-on must still find it"

    def test_without_adapters_there_is_nothing_to_collide_with(self, svc: Services) -> None:
        c = construct(dna(150) + TWIST_FIVE_PRIME + dna(150, 5), dna(300, 11))
        assert not hits(c, svc, vendor="none")

    def test_an_adapter_colliding_with_itself_is_not_the_designs_problem(
        self, svc: Services
    ) -> None:
        """A finding lying wholly inside vendor sequence has no construct
        coordinate and nothing BT5 can do about it."""
        assert not hits(construct(dna(400, 21), dna(300, 22)), svc, vendor=ADAPTER_ON)


class TestDetection:
    def test_reports_the_repeat_at_the_length_it_has(self, svc: Services) -> None:
        found = hits(planted(dna(30, 9)), svc)
        assert len(found) == 1
        assert found[0].detail["length"] == 30.0

    def test_a_pair_below_the_floor_is_not_reported(self, svc: Services) -> None:
        assert not hits(planted(dna(10, 33)), svc)

    def test_clean_sequence_is_clean(self, svc: Services) -> None:
        assert not hits(construct(dna(600, 21), dna(400, 22)), svc)

    def test_a_hard_finding_fails_the_evaluation(self, svc: Services) -> None:
        c = planted(dna(30, 9))
        assert not SynthesisRepeats().evaluate(c, context(), svc).passes

    def test_a_warning_does_not_fail_the_evaluation(self, svc: Services) -> None:
        c = planted("ATTATTAATTATAA")
        result = SynthesisRepeats().evaluate(c, context(), svc)
        assert result.breaches
        assert result.passes


class TestContract:
    def test_it_is_hard_repair_and_carries_no_objective_weight(self) -> None:
        assert SynthesisRepeats.enforcement is Enforcement.HARD_REPAIR
        assert SynthesisRepeats.default_weight == 0.0

    def test_it_steers_below_f1(self) -> None:
        """F1 steers on the whole plasmid and already suppresses most of what E5
        would; the increment is the short GC-rich repeat and the adapter."""
        assert 0.0 < SynthesisRepeats.steering_weight < DirectRepeats.steering_weight

    def test_it_declares_no_fold_engine_calibration(self) -> None:
        """`engine_calibration` names the FOLD engine a rule's thresholds are
        measured on. Putting the Tm parameter set there would make
        `check_engine_calibration` raise against ViennaRNA on every run."""
        assert SynthesisRepeats.engine_calibration is None

    def test_it_is_not_a_lattice_rule(self) -> None:
        assert SynthesisRepeats().lattice_terms(None) is None

    def test_it_applies_in_every_context(self) -> None:
        for modality in Modality:
            assert SynthesisRepeats().gate(slot(modality=modality))

    def test_absurd_parameters_are_refused(self) -> None:
        with pytest.raises(ValueError, match="60 C duplex"):
            SynthesisRepeats(min_len=MIN_LENGTH_BP - 1)
        with pytest.raises(ValueError, match="must not be below"):
            SynthesisRepeats(min_len=20, hard_len=15)
        with pytest.raises(ValueError, match="unknown vendor"):
            SynthesisRepeats(vendor="acme")
        with pytest.raises(ValueError, match="outside any real PCR"):
            SynthesisRepeats(anneal_c=5.0)

    def test_it_is_registered_under_its_brief_row(self) -> None:
        assert get("e5_synthesis_repeats") is SynthesisRepeats
        assert SynthesisRepeats.brief_ref == "2.E5"
