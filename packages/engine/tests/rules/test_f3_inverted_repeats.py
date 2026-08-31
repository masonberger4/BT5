"""F3: inverted repeats, which are not direct repeats with a minus sign.

The claim under test is mechanistic. A direct repeat is lost by deletion and
codon choice fixes it; an inverted repeat extrudes a cruciform, stalls the fork
and is cleaved by SbcCD, and the answer is usually a strain and a temperature.
That is why F3 has its own bands, its own citations and HARD_CHECK enforcement
rather than reusing F1's.
"""

from __future__ import annotations

import numpy as np
import pytest
from bt5.core.context import Modality
from bt5.core.registry import discover, get
from bt5.core.services import Services
from bt5.core.spec import Enforcement
from bt5.core.types import Construct, Interval, Segment, SegmentKind, Topology, reverse_complement
from bt5.rules.catalog.f3_inverted_repeats import AAV_STEM_BP, HARD_STEM_BP, InvertedRepeats, _loop
from bt5.vector.kmers import ConstructKmerIndex
from conftest import construct, context, slot

discover()


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


def hairpin(arm: str, loop: int, *, lead: int = 200) -> Construct:
    seq = dna(lead) + arm + dna(loop, 67) + reverse_complement(arm) + dna(200, 71)
    return construct(seq[: lead + 20], seq[lead + 20 :])


def hits(rule: InvertedRepeats, c: Construct, svc: Services, ctx=None):
    return rule.evaluate(c, ctx or context(), svc).breaches


class TestLoop:
    def test_a_plain_gap(self) -> None:
        assert _loop(Interval(100, 130), Interval(140, 170), 800, circular=True) == 10

    def test_abutting_arms_are_a_perfect_palindrome(self) -> None:
        assert _loop(Interval(100, 130), Interval(130, 160), 800, circular=True) == 0

    def test_a_loop_across_the_origin_is_measured_forward(self) -> None:
        """A negative or whole-plasmid loop would put the most extrudable
        hairpins in the least alarming band."""
        assert _loop(Interval(780, 800), Interval(10, 30), 800, circular=True) == 10


class TestDetection:
    def test_finds_a_hairpin_with_its_geometry(self, svc: Services) -> None:
        breaches = hits(InvertedRepeats(), hairpin(dna(30, 11), 10), svc)
        assert breaches
        found = breaches[0]
        assert found.detail["stem"] >= 30.0
        assert found.detail["loop"] == 10.0

    def test_non_palindromic_sequence_is_clean(self, svc: Services) -> None:
        c = construct(dna(400, 21), dna(400, 22))
        assert not hits(InvertedRepeats(), c, svc)

    def test_a_direct_repeat_is_not_reported_as_an_inverted_one(self, svc: Services) -> None:
        """The two are mechanically different and must not be conflated."""
        unit = dna(40, 9)
        seq = dna(200) + unit + dna(300, 5) + unit + dna(200, 7)
        c = construct(seq[:300], seq[300:])
        assert not hits(InvertedRepeats(), c, svc)

    def test_a_perfect_palindrome_is_named_as_one(self, svc: Services) -> None:
        """Arms abutting with no loop is the most extrudable configuration."""
        breach = hits(InvertedRepeats(), hairpin(dna(35, 13), 0), svc)[0]
        assert breach.detail["perfect_palindrome"] == "yes"
        assert "PERFECT PALINDROME" in breach.message

    def test_total_length_counts_both_arms_and_the_loop(self, svc: Services) -> None:
        breach = hits(InvertedRepeats(), hairpin(dna(30, 11), 10), svc)[0]
        stem, loop, total = (
            breach.detail["stem"],
            breach.detail["loop"],
            breach.detail["total"],
        )
        assert total == 2 * stem + loop


class TestBands:
    def test_a_long_stem_is_hard(self, svc: Services) -> None:
        breach = hits(InvertedRepeats(), hairpin(dna(35, 13), 10), svc)[0]
        assert breach.detail["severity"] in ("hard", "unbuildable")

    def test_the_total_palindrome_band_is_enforced_and_not_merely_reported(
        self, svc: Services
    ) -> None:
        """The hard band has two terms -- stem >= 30 OR arms+loop >= 60 -- but
        `passes` tested only the stem, so a pair this rule labels "hard" in its
        own breach message still reported passes=True. F3 is HARD_CHECK, so
        `passes` IS the enforcement surface: the second term was reported and
        never enforced.

        This fixture sits in exactly that gap: a 27 bp stem with a 10 bp loop is
        under the 30 bp stem threshold and over the 60 bp total.
        """
        c = hairpin(dna(25, 41), 10)
        ev = InvertedRepeats().evaluate(c, context(slot(modality=Modality.PLASMID_TRANSIENT)), svc)
        hard = [b for b in ev.breaches if b.detail["severity"] == "hard"]
        assert hard, "27 bp stem + 10 bp loop is 64 bp of palindrome, over the 60 bp band"
        assert not ev.passes, "a HARD_CHECK rule that calls a finding 'hard' must not pass it"

    def test_a_very_long_palindrome_says_do_not_build(self, svc: Services) -> None:
        """SbcCD cleaves long palindromes and destroys the replicon."""
        breach = hits(InvertedRepeats(), hairpin(dna(160, 19), 4), svc)[0]
        assert breach.detail["severity"] == "unbuildable"
        assert "Do not build" in breach.message
        assert "sbcC" in breach.message

    def test_severity_orders_by_stem_length(self, svc: Services) -> None:
        short = hits(InvertedRepeats(), hairpin(dna(16, 23), 10), svc)
        long = hits(InvertedRepeats(), hairpin(dna(160, 19), 4), svc)
        assert short[0].magnitude < long[0].magnitude

    def test_the_strain_protocol_travels_with_the_finding(self, svc: Services) -> None:
        """The answer to an unavoidable IR is a strain and a temperature, so the
        finding has to carry it -- a user told only 'palindrome' cannot act."""
        breach = hits(InvertedRepeats(), hairpin(dna(35, 13), 10), svc)[0]
        assert "30 C" in breach.message or "42 C" in breach.message


class TestAav:
    """AAV-GPseq mapped truncation hotspots exactly to inverted repeats."""

    def test_aav_uses_a_stricter_stem_threshold(self) -> None:
        assert AAV_STEM_BP < HARD_STEM_BP

    def test_a_stem_a_plasmid_tolerates_is_hard_for_aav(self, svc: Services) -> None:
        c = hairpin(dna(22, 29), 10)
        plasmid = context(slot(modality=Modality.PLASMID_TRANSIENT))
        aav = context(slot(modality=Modality.AAV, host=slot().host))

        assert InvertedRepeats().evaluate(c, plasmid, svc).passes
        assert not InvertedRepeats().evaluate(c, aav, svc).passes


class TestExemptRegions:
    """AAV ITRs are 145 bp palindromes BY DESIGN."""

    def itr_construct(self, arm: str, *, exempt_both: bool) -> Construct:
        lead = dna(200)
        loop = dna(10, 67)
        rc = reverse_complement(arm)
        seq = lead + arm + loop + rc + dna(200, 71)
        first = Interval(len(lead), len(lead) + len(arm))
        second_start = len(lead) + len(arm) + len(loop)
        second = Interval(second_start, second_start + len(rc))
        segments = [Segment(first, SegmentKind.WHITELISTED_REPEAT, "ITR")]
        if exempt_both:
            segments.append(Segment(second, SegmentKind.WHITELISTED_REPEAT, "ITR"))
        segments.append(Segment(Interval(0, len(lead)), SegmentKind.DESIGNABLE_CDS, "cds"))
        return Construct(seq, Topology.CIRCULAR, tuple(segments))

    def test_an_itr_palindrome_is_not_reported(self, svc: Services) -> None:
        c = self.itr_construct(dna(40, 31), exempt_both=True)
        assert not hits(InvertedRepeats(), c, svc)

    def test_a_half_exempt_hairpin_is_still_reported(self, svc: Services) -> None:
        c = self.itr_construct(dna(40, 31), exempt_both=False)
        assert hits(InvertedRepeats(), c, svc)


class TestContract:
    def test_it_is_hard_check_not_hard_repair(self) -> None:
        """The dominant cases -- ITRs, shRNA hairpins, an IR in the user's own
        backbone -- cannot be fixed by codon choice at all."""
        assert InvertedRepeats.enforcement is Enforcement.HARD_CHECK

    def test_a_hard_check_rule_does_not_steer_the_dp(self) -> None:
        """Steering toward something the solver must never chase is a category
        error, not a tuning choice."""
        assert InvertedRepeats.steering_weight == 0.0
        assert InvertedRepeats.default_weight == 0.0

    def test_fixability_needs_both_arms_in_the_cds(self, svc: Services) -> None:
        """Recoding one arm of a hairpin whose other arm is immutable backbone
        cannot break the pairing without changing the backbone."""
        arm = dna(30, 11)
        loop = dna(10, 67)
        seq = dna(200) + arm + loop + reverse_complement(arm) + dna(200, 71)
        cds_end = 200 + len(arm)  # the CDS stops between the two arms
        c = Construct(
            seq,
            Topology.CIRCULAR,
            (
                Segment(Interval(0, cds_end), SegmentKind.DESIGNABLE_CDS, "cds"),
                Segment(Interval(cds_end, len(seq)), SegmentKind.BACKBONE, "vector"),
            ),
        )
        breaches = hits(InvertedRepeats(), c, svc)
        assert breaches
        assert not any(b.fixable_by_codon_choice for b in breaches)

    def test_it_is_not_a_lattice_rule(self) -> None:
        assert InvertedRepeats().lattice_terms(None) is None

    def test_absurd_parameters_are_refused(self) -> None:
        with pytest.raises(ValueError, match="noise, not a finding"):
            InvertedRepeats(min_stem=4)
        with pytest.raises(ValueError, match="max_loop must be"):
            InvertedRepeats(max_loop=-1)

    def test_it_is_registered_under_its_brief_row(self) -> None:
        assert get("f3_inverted_repeats") is InvertedRepeats
        assert InvertedRepeats.brief_ref == "2.F3"
