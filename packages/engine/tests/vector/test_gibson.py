"""Junction uniqueness: can this construct be assembled unambiguously?

The findings that matter here are the ones a cloning tool cannot produce,
because they depend on which bases the optimizer chose. An arm is half backbone
and half insert, so whether an ambiguity is fixable at all is a question about
the CDS/backbone boundary -- which is this lane's whole subject.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest
from bt5.core.types import Construct, Interval, Topology, reverse_complement
from bt5.vector import VectorBackbone, assemble, insertion_site_from_interval
from bt5.vector.gibson import (
    MAX_ARM_BP,
    MIN_ARM_BP,
    MIN_ARM_TM_C,
    HomologyArm,
    TmConditions,
    build_arm,
    gc_fraction,
    insert_shared_repeats,
    junction_points,
    melting_temperature,
    occurrences,
    plan_junctions,
    shortest_usable_arm,
)
from bt5.vector.gibson import _shared_between as shared_between
from conftest import translate

GC_RICH = "GCTGACCTGGAACGTCTGCA"
AT_RICH = "ATATTTAAATATTTAAATATTTAAATATTTAAATATTTAAAT"


def padding(n: int, *, seed: int) -> str:
    """Random ACGT. Seeded: an unseeded pad is itself a repeat waiting to happen."""
    rng = np.random.default_rng(seed)
    return "".join("ACGT"[i] for i in rng.integers(0, 4, n))


def cds_from(protein: str, *, alternating: bool = True) -> str:
    """Encode a protein, either diversifying synonyms or collapsing to one codon.

    The collapsed form is the max-CAI failure mode: one codon per amino acid
    turns a repetitive protein into a perfect nucleotide repeat.
    """
    from Bio.Data import CodonTable

    table = CodonTable.unambiguous_dna_by_id[1]
    synonyms: dict[str, list[str]] = {}
    for codon, aa in sorted(table.forward_table.items()):
        synonyms.setdefault(aa, []).append(codon)
    out = ["ATG"]
    for i, aa in enumerate(protein[1:]):
        options = synonyms[aa]
        out.append(options[i % len(options)] if alternating else options[0])
    return "".join(out) + "TAA"


def plasmid(
    insert_region: str, *, before: str, after: str, circular: bool = True
) -> VectorBackbone:
    from bt5.core.types import Feature

    sequence = before + insert_region + after
    return VectorBackbone(
        sequence=sequence,
        topology=Topology.CIRCULAR if circular else Topology.LINEAR,
        features=(
            Feature(
                interval=Interval(len(before), len(before) + len(insert_region)),
                kind="CDS",
                qualifiers={"label": ("transgene",), "transl_table": ("1",)},
                uid="cds0",
            ),
        ),
        name="synthetic",
    )


def build(
    before: str,
    after: str,
    cds: str,
    placeholder: str | None = None,
    *,
    circular: bool = True,
):
    """Assemble `cds` into a backbone made of `before` + placeholder + `after`."""
    region = placeholder if placeholder is not None else cds
    bb = plasmid(region, before=before, after=after, circular=circular)
    site = insertion_site_from_interval(
        Interval(len(before), len(before) + len(region)), table_id=1
    )
    return assemble(bb, cds, protein=translate(cds), table_id=1, site=site)


class TestMeltingTemperature:
    def test_a_longer_overlap_is_warmer(self) -> None:
        seq = "ATGGCTAGCAAAGGAGAAGAACTTTTCACT"
        assert melting_temperature(seq[:15]) < melting_temperature(seq[:30])

    def test_gc_raises_the_melting_temperature(self) -> None:
        assert melting_temperature("GCGCGCGCGCGCGCGCGCGC") > melting_temperature(
            "ATATATATATATATATATAT"
        )

    def test_conditions_change_the_answer_and_are_reported(self) -> None:
        """A Tm without its conditions cannot be compared to anything."""
        seq = "ATGGCTAGCAAAGGAGAAGAACTTTTCACT"
        low = TmConditions(na_mm=10.0)
        assert melting_temperature(seq, low) < melting_temperature(seq)
        assert "Na+ 10 mM" in low.describe()
        assert "SantaLucia_Hicks_2004" in low.describe()

    def test_the_parameter_table_is_pinned_not_defaulted(self) -> None:
        """Biopython's default is the 1997 table; leaving it unset would let a
        library default move every number this module reports."""
        from Bio.SeqUtils import MeltingTemp  # noqa: N812

        seq = "ATGGCTAGCAAAGGAGAAGAACTTTTCACT"
        pinned = melting_temperature(seq)
        stale = MeltingTemp.Tm_NN(  # type: ignore[no-untyped-call]
            seq, nn_table=MeltingTemp.DNA_NN3, Na=50, Mg=0, dnac1=250, dnac2=250
        )
        assert pinned != pytest.approx(stale)

    def test_gc_fraction(self) -> None:
        assert gc_fraction("GCGC") == 1.0
        assert gc_fraction("ATAT") == 0.0
        assert gc_fraction("") == 0.0


class TestOccurrences:
    def construct(self, sequence: str, *, circular: bool = True) -> Construct:
        from bt5.core.types import Segment, SegmentKind

        return Construct(
            sequence=sequence,
            topology=Topology.CIRCULAR if circular else Topology.LINEAR,
            segments=(Segment(interval=Interval(0, len(sequence)), kind=SegmentKind.BACKBONE),),
        )

    def test_a_unique_sequence_occurs_once(self) -> None:
        c = self.construct(padding(500, seed=1) + GC_RICH + padding(500, seed=2))
        assert occurrences(c, GC_RICH, skip=Interval(500, 500 + len(GC_RICH))) == ()

    def test_a_duplicated_sequence_is_found(self) -> None:
        c = self.construct(
            padding(300, seed=1) + GC_RICH + padding(300, seed=2) + GC_RICH + padding(300, seed=3)
        )
        found = occurrences(c, GC_RICH, skip=Interval(300, 300 + len(GC_RICH)))
        assert [f.start for f in found] == [620]

    def test_a_reverse_complement_hit_counts(self) -> None:
        """What anneals in Gibson is a single-stranded overhang, so the
        reverse-complement copy is a real mis-annealing site."""
        c = self.construct(
            padding(300, seed=1)
            + GC_RICH
            + padding(300, seed=2)
            + reverse_complement(GC_RICH)
            + padding(300, seed=3)
        )
        found = occurrences(c, GC_RICH, skip=Interval(300, 300 + len(GC_RICH)))
        assert [(f.start, f.strand) for f in found] == [(620, -1)]

    def test_a_hit_spanning_the_origin_is_found(self) -> None:
        """The whole point of holding topology: a match across position 0."""
        half = len(GC_RICH) // 2
        c = self.construct(GC_RICH[half:] + padding(500, seed=4) + GC_RICH[:half])
        found = occurrences(c, GC_RICH)
        assert len(found) == 1
        assert found[0].end > c.length, "the match wraps"
        assert c.slice(found[0]) == GC_RICH

    def test_a_linear_construct_has_no_origin_junction(self) -> None:
        half = len(GC_RICH) // 2
        c = self.construct(GC_RICH[half:] + padding(500, seed=4) + GC_RICH[:half], circular=False)
        assert occurrences(c, GC_RICH) == ()


class TestArmGeometry:
    def test_the_arm_straddles_the_junction(self) -> None:
        """A split overlap: half the bases come from each fragment."""
        cds = cds_from("M" + "AKLEDGRT" * 12)
        assembly = build(padding(400, seed=5), padding(400, seed=6), cds)
        at = assembly.cds_interval.start
        arm = build_arm(assembly.construct, "5'", at, 20)
        assert arm is not None
        assert arm.interval.start < at < arm.interval.end
        assert arm.length == 20

    def test_the_arm_bases_match_the_construct(self) -> None:
        cds = cds_from("M" + "AKLEDGRT" * 12)
        assembly = build(padding(400, seed=5), padding(400, seed=6), cds)
        arm = build_arm(assembly.construct, "5'", assembly.cds_interval.start, 24)
        assert arm is not None
        assert arm.sequence == assembly.construct.slice(arm.interval)

    def test_an_arm_across_the_origin_wraps(self) -> None:
        cds = cds_from("M" + "AKLEDGRT" * 12)
        assembly = build("", padding(600, seed=7), cds)
        arm = build_arm(assembly.construct, "5'", 0, 20)
        assert arm is not None
        assert arm.interval.start > assembly.construct.length - 20
        assert len(arm.sequence) == 20

    def test_junction_points_are_the_cds_edges(self) -> None:
        cds = cds_from("M" + "AKLEDGRT" * 12)
        assembly = build(padding(400, seed=5), padding(400, seed=6), cds)
        names = [n for n, _ in junction_points(assembly)]
        coords = [c for _, c in junction_points(assembly)]
        assert names == ["5' backbone-insert", "3' insert-backbone"]
        assert coords == [assembly.cds_interval.start, assembly.cds_interval.end]


class TestArmSelection:
    def test_a_clean_junction_yields_a_short_unique_warm_arm(self) -> None:
        cds = cds_from("M" + "AKLEDGRTQWF" * 20)
        assembly = build(padding(600, seed=8), padding(600, seed=9), cds)
        arm = shortest_usable_arm(assembly.construct, "5'", assembly.cds_interval.start)
        assert arm is not None
        assert arm.usable
        assert MIN_ARM_BP <= arm.length <= MAX_ARM_BP

    def test_the_shortest_workable_arm_is_chosen(self) -> None:
        """Longer is not better: every extra base is another chance at a repeat."""
        cds = cds_from("M" + "AKLEDGRTQWF" * 20)
        assembly = build(padding(600, seed=8), padding(600, seed=9), cds)
        arm = shortest_usable_arm(assembly.construct, "5'", assembly.cds_interval.start)
        assert arm is not None
        shorter = build_arm(assembly.construct, "5'", assembly.cds_interval.start, arm.length - 1)
        assert shorter is None or not shorter.usable

    def test_an_at_rich_junction_needs_a_longer_overlap(self) -> None:
        """The temperature floor is what drives the chosen length up.

        Nearest-neighbour Tm rises with length whatever the GC, so an AT-rich
        junction is not unreachable -- it just costs bases. That is exactly why
        the vendor window has an upper bound as well as a lower one.
        """
        at_rich = build(
            padding(400, seed=10) + AT_RICH,
            padding(600, seed=11),
            cds_from("M" + "KNIKNIKNIK" * 20, alternating=False),
        )
        gc_rich = build(
            padding(400, seed=12) + GC_RICH * 2,
            padding(600, seed=13),
            cds_from("M" + "AGPRAWGPRA" * 20, alternating=False),
        )
        cold = shortest_usable_arm(at_rich.construct, "5'", at_rich.cds_interval.start)
        warm = shortest_usable_arm(gc_rich.construct, "5'", gc_rich.cds_interval.start)
        assert cold is not None
        assert warm is not None
        assert cold.gc < warm.gc
        assert cold.length > warm.length
        assert cold.warm_enough, "reachable, just not cheaply"

    def test_short_overlaps_at_an_at_rich_junction_are_rejected_as_cold(self) -> None:
        assembly = build(
            padding(400, seed=10) + AT_RICH,
            padding(600, seed=11),
            cds_from("M" + "KNIKNIKNIK" * 20, alternating=False),
        )
        short = build_arm(assembly.construct, "5'", assembly.cds_interval.start, MIN_ARM_BP)
        assert short is not None
        assert short.unique, "it is the temperature that disqualifies it, not ambiguity"
        assert short.tm_c < MIN_ARM_TM_C
        assert not short.usable

    def test_a_junction_too_cold_for_a_narrow_window_is_reported(self) -> None:
        """A user with a fixed 25 bp overlap protocol has a real constraint."""
        assembly = build(
            padding(400, seed=10) + AT_RICH,
            padding(600, seed=11),
            cds_from("M" + "KNIKNIKNIK" * 20, alternating=False),
        )
        plan = plan_junctions(assembly, max_bp=25)
        cold = [n for n in plan.notes if "below the" in n.summary]
        assert cold, "a junction that cannot reach the floor in the window must be said"
        assert not plan.usable
        assert all("AT-rich" in n.summary for n in cold)
        assert all(n.action for n in cold)

    def test_arms_at_one_junction_are_nested(self) -> None:
        """The property the fallback rests on.

        Every arm is centred on the same boundary, so the arm at length L is a
        substring of the arm at L+1. That makes uniqueness monotone: a duplicate
        of a longer arm contains the shorter arm too, so once an arm is unique
        every longer one is. If this stopped holding, "longest unique" would no
        longer be the warmest unambiguous arm.
        """
        cds = cds_from("M" + "AKLEDGRTQWF" * 25)
        assembly = build(padding(700, seed=70), padding(700, seed=71), cds)
        at = assembly.cds_interval.start
        previous: str | None = None
        for length in range(MIN_ARM_BP, MAX_ARM_BP + 1):
            arm = build_arm(assembly.construct, "j", at, length)
            assert arm is not None
            if previous is not None:
                assert previous in arm.sequence, f"arm at {length} does not contain {length - 1}"
            previous = arm.sequence

    def test_the_fallback_is_the_longest_unique_arm_not_the_shortest(self) -> None:
        """It is the warmest unambiguous arm the window has, so the temperature
        the note quotes is one the user could actually reach."""
        assembly = build(
            padding(400, seed=10) + AT_RICH,
            padding(600, seed=11),
            cds_from("M" + "KNIKNIKNIK" * 20, alternating=False),
        )
        at = assembly.cds_interval.start
        arm = shortest_usable_arm(assembly.construct, "j", at, max_bp=25)
        assert arm is not None
        assert arm.unique
        assert not arm.warm_enough
        assert arm.length == 25, "the shortest unique arm would quote a colder, unreachable Tm"
        longest = build_arm(assembly.construct, "j", at, 25)
        assert longest is not None
        assert arm.tm_c == pytest.approx(longest.tm_c)

    def test_a_capped_arm_is_shorter_than_the_window_allows(self) -> None:
        """On a linear construct an arm cannot grow past the ends, so the arm
        BT5 measured and the maximum it was asked for come apart -- the only
        case where quoting one in place of the other is visibly wrong."""
        assembly = self.capped()
        arm = shortest_usable_arm(
            assembly.construct, "j", assembly.cds_interval.start, max_bp=MAX_ARM_BP
        )
        assert arm is not None
        assert arm.length < MAX_ARM_BP, "the construct ends before a full-width arm fits"

    def capped(self):
        """A LINEAR construct whose 5' junction is too near the end for a full arm."""
        return build(
            AT_RICH[:12],
            padding(600, seed=72),
            cds_from("M" + "KNIKNIKNIK" * 20, alternating=False),
            circular=False,
        )

    def test_the_cold_note_quotes_the_capped_length_not_the_window(self) -> None:
        """The one case where the arm measured and the maximum asked for come
        apart -- and therefore the only one that can tell a correct report from
        a number that happens to coincide."""
        assembly = self.capped()
        plan = plan_junctions(assembly, max_bp=MAX_ARM_BP)
        arms = {a.junction: a for a in plan.arms}
        cold = [n for n in plan.notes if "below the" in n.summary]
        assert cold, "an AT-rich capped junction cannot reach the floor"
        checked = 0
        for note in cold:
            junction = next(j for j in arms if j in note.summary)
            arm = arms[junction]
            if arm.length == MAX_ARM_BP:
                continue  # indistinguishable here; the capped junction is the test
            assert f"at {arm.length} bp" in note.summary
            assert f"at {MAX_ARM_BP} bp" not in note.summary
            checked += 1
        assert checked, "no junction was actually capped, so nothing was proved"

    def test_the_cold_note_quotes_the_arm_it_chose(self) -> None:
        """Mixing the requested maximum with a different arm's measurement is a
        wrong number on the report, not a wording preference."""
        assembly = build(
            padding(400, seed=10) + AT_RICH,
            padding(600, seed=11),
            cds_from("M" + "KNIKNIKNIK" * 20, alternating=False),
        )
        plan = plan_junctions(assembly, max_bp=25)
        cold = [n for n in plan.notes if "below the" in n.summary]
        assert cold
        arms = {a.junction: a for a in plan.arms}
        for note in cold:
            junction = next(j for j in arms if j in note.summary)
            arm = arms[junction]
            assert f"at {arm.length} bp" in note.summary
            assert f"{arm.tm_c:.1f} C" in note.summary


class TestSharedArms:
    def arm(self, name: str, sequence: str, start: int) -> HomologyArm:
        return HomologyArm(
            junction=name,
            interval=Interval(start, start + len(sequence)),
            sequence=sequence,
            tm_c=55.0,
            gc=gc_fraction(sequence),
        )

    def test_two_identical_short_arms_are_reported(self) -> None:
        """The most ambiguous pair there is. A fixed 20 bp comparison window
        does not fit inside a 15 bp arm and misses this entirely."""
        seq = "ATGGCTAGCAAAGGA"
        shared = shared_between(None, [self.arm("5'", seq, 0), self.arm("3'", seq, 500)])  # type: ignore[arg-type]
        assert len(shared) == 1
        assert shared[0][:2] == ("5'", "3'")

    def test_a_partial_overlap_is_enough(self) -> None:
        """Two arms sharing their last 20 bases anneal to each other's partners
        just as happily as two identical arms."""
        tail = "GCTGACCTGGAACGTCTGCA"
        a = self.arm("5'", padding(15, seed=20) + tail, 0)
        b = self.arm("3'", tail + padding(15, seed=21), 500)
        assert len(shared_between(None, [a, b])) == 1  # type: ignore[arg-type]

    def test_a_reverse_complement_overlap_counts(self) -> None:
        """It lets the insert assemble backwards, not merely in the wrong place."""
        tail = "GCTGACCTGGAACGTCTGCA"
        a = self.arm("5'", padding(15, seed=20) + tail, 0)
        b = self.arm("3'", reverse_complement(tail) + padding(15, seed=21), 500)
        assert len(shared_between(None, [a, b])) == 1  # type: ignore[arg-type]

    def test_unrelated_arms_share_nothing(self) -> None:
        a = self.arm("5'", padding(30, seed=22), 0)
        b = self.arm("3'", padding(30, seed=23), 500)
        assert shared_between(None, [a, b]) == ()  # type: ignore[arg-type]


class TestInsertSharedRepeats:
    """The plan's rule stated directly: no exact repeat >=20 bp between the
    insert and any other part of the construct. A repeat gives a fragment a
    second place to anneal whether or not it sits at a junction."""

    #: 30 bp, non-periodic, in frame. Deliberately not a tandem like (GCC)n:
    #: a periodic motif is reported by its period, which would make this test
    #: about the repeat scan's tandem handling rather than about the boundary.
    SHARED = "GCTGACCTGGAACGTCTGCAGGTACATTGGC"[:30]

    def test_a_repeat_planted_in_the_backbone_and_the_insert_is_found(self) -> None:
        cds = "ATG" + self.SHARED + cds_from("M" + "KLEDGRTQWF" * 8)[3:]
        before = padding(400, seed=30) + self.SHARED + padding(200, seed=31)
        assembly = build(before, padding(500, seed=32), cds)
        pairs = insert_shared_repeats(assembly.construct, assembly.cds_interval)
        assert pairs, "a 30 bp insert/backbone identity must be reported"
        assert max(p.length for p in pairs) >= 30

    def test_the_crossing_repeat_is_reported_as_a_liability(self) -> None:
        cds = "ATG" + self.SHARED + cds_from("M" + "KLEDGRTQWF" * 8)[3:]
        before = padding(400, seed=30) + self.SHARED + padding(200, seed=31)
        plan = plan_junctions(build(before, padding(500, seed=32), cds))
        crossing = [n for n in plan.notes if "shared between the designed insert" in n.summary]
        assert crossing
        assert all(n.kind == "liability" for n in crossing)
        assert any("recA- strain does not suppress it" in n.summary for n in crossing), (
            "a 30 bp repeat is in the RecA-INDEPENDENT regime, and telling a user "
            "their strain covers it would be false"
        )

    def test_a_repeat_entirely_inside_the_insert_is_somebody_elses_finding(self) -> None:
        """Insert-internal repeats are a synthesis and stability problem, already
        reported by the repeat scan. Reporting them here would double-count."""
        motif = "GCC" * 10
        cds = "ATG" + motif + padding(60, seed=33).replace("T", "A") + motif + "TAA"
        cds = cds[: len(cds) // 3 * 3]
        assembly = build(padding(500, seed=34), padding(500, seed=35), cds)
        pairs = insert_shared_repeats(assembly.construct, assembly.cds_interval)
        assert pairs == ()

    def test_a_repeat_entirely_inside_the_backbone_is_not_this_rule(self) -> None:
        """Nothing codon choice does can fix it, and it is reported elsewhere."""
        motif = "GCTGACCTGGAACGTCTGCAGGTAC"
        cds = cds_from("M" + "KLEDGRTQWF" * 10)
        before = padding(300, seed=36) + motif + padding(300, seed=37) + motif
        assembly = build(before, padding(400, seed=38), cds)
        pairs = insert_shared_repeats(assembly.construct, assembly.cds_interval)
        assert pairs == ()

    def test_a_short_identity_is_below_the_vendor_threshold(self) -> None:
        motif = "GCTGACCTGGAA"  # 12 bp
        cds = "ATG" + motif + cds_from("M" + "KLEDGRTQWF" * 8)[3:]
        before = padding(400, seed=39) + motif + padding(200, seed=40)
        assembly = build(before, padding(500, seed=41), cds)
        assert insert_shared_repeats(assembly.construct, assembly.cds_interval) == ()


class TestPlanJunctions:
    def test_a_clean_construct_plans_two_usable_arms(self) -> None:
        cds = cds_from("M" + "AKLEDGRTQWF" * 25)
        assembly = build(padding(700, seed=50), padding(700, seed=51), cds)
        plan = plan_junctions(assembly)
        assert plan.usable
        assert [a.junction for a in plan.arms] == [
            "5' backbone-insert",
            "3' insert-backbone",
        ]
        assert not plan.notes

    def ambiguous_junction(self):
        """A construct whose entire 5' junction window occurs twice.

        The duplicate is half backbone and half insert, so no overlap length in
        the vendor range is unique and the ambiguity is partly the optimizer's
        to fix -- which is the case the whole module is built around.
        """
        cds = cds_from("M" + "AKLEDGRTQWF" * 25)
        before = padding(700, seed=52)
        window = before[-MAX_ARM_BP:] + cds[:MAX_ARM_BP]
        after = padding(300, seed=53) + window + padding(300, seed=54)
        return build(before, after, cds)

    def test_an_ambiguous_junction_is_reported(self) -> None:
        plan = plan_junctions(self.ambiguous_junction())
        ambiguous = [n for n in plan.notes if "is unique" in n.summary]
        assert ambiguous, "a duplicated junction window must be reported"
        assert not plan.usable

    def test_the_ambiguity_is_named_as_fixable_by_codon_choice(self) -> None:
        """The distinction that stops the solver chasing an unreachable constraint."""
        plan = plan_junctions(self.ambiguous_junction())
        ambiguous = [n for n in plan.notes if "is unique" in n.summary]
        assert any("different codons" in n.summary for n in ambiguous)
        assert all("re-run" in n.action for n in ambiguous)

    def test_a_junction_arm_always_reaches_the_designed_cds(self) -> None:
        """The geometry behind calling junction ambiguity fixable.

        An arm is centred on the boundary, so it always contains designed bases.
        A duplicate therefore has to match the insert half too, and changing
        codons there breaks the match. If this ever stopped holding -- an
        off-centre arm, or an insert shorter than half an overlap -- the
        "re-run" advice on the ambiguity note would become a false promise.
        """
        cds = cds_from("M" + "AKLEDGRTQWF" * 25)
        assembly = build(padding(700, seed=63), padding(700, seed=64), cds)
        construct = assembly.construct
        for _, at in junction_points(assembly):
            for length in (MIN_ARM_BP, 25, MAX_ARM_BP):
                arm = build_arm(construct, "j", at, length)
                assert arm is not None
                overlaps = any(
                    arm.interval.start < e.end and e.start < arm.interval.end
                    for e in construct.editable
                )
                assert overlaps, f"a {length} bp arm at {at} missed the CDS entirely"

    def test_the_ambiguity_note_does_not_promise_an_unreachable_fix(self) -> None:
        plan = plan_junctions(self.ambiguous_junction())
        for note in plan.notes:
            if "is unique" in note.summary:
                assert "spans the designed CDS" in note.summary

    def test_every_note_is_located_and_actionable(self) -> None:
        cds = cds_from("M" + "KNIKLEDGRT" * 25)
        at_rich = "ATATTTAAATATTTAAATATTTAAATATTTAAATATTTAAAT"
        assembly = build(padding(600, seed=54) + at_rich, padding(600, seed=55), cds)
        plan = plan_junctions(assembly)
        for note in plan.notes:
            assert note.bears_on == "assembly"
            assert note.action, "a liability with no action is just an alarm"
            assert note.kind in ("liability", "unavailable")

    def test_no_note_predicts_an_outcome(self) -> None:
        """BT5 reports what is present, never what will happen."""
        banned = ("predict", "titer", "yield", "fold-improvement", "expression level")
        cds = cds_from("M" + "KNIKLEDGRT" * 25)
        at_rich = "ATATTTAAATATTTAAATATTTAAATATTTAAATATTTAAAT"
        assembly = build(padding(600, seed=54) + at_rich, padding(600, seed=55), cds)
        for note in plan_junctions(assembly).notes:
            assert not any(word in note.summary.lower() for word in banned)

    def test_a_linear_construct_with_a_junction_at_the_very_end(self) -> None:
        """No arm can exist there, and saying so is better than reporting one."""

        cds = cds_from("M" + "AKLEDGRTQWF" * 20)
        assembly = build(padding(600, seed=56), padding(600, seed=57), cds)
        linear = replace(assembly.construct, topology=Topology.LINEAR)
        assert build_arm(linear, "5'", 2, 40) is None
