"""Whether the assembled construct can be built unambiguously by Gibson.

BT5 does not design an assembly. It has no opinion about which fragments you
order or where you cut, and it will not emit a protocol. What it does own is the
question codon choice actually controls: given this construct, is there a
homology arm at each junction that anneals in exactly one place?

Repeats and extreme GC are the documented cause of Gibson misassembly
(https://blog.addgene.org/plasmids-101-gibson-assembly), and NEB's stated
overlap window is 15-40 bp at a melting temperature above 48 C. Both halves of
that matter and they pull against each other: a short arm is more likely to
occur twice, a long one is more likely to contain a repeat, and an AT-rich
junction may not reach the temperature floor at any length in the window. So the
useful answer is not pass/fail at one arbitrary length. It is the SHORTEST arm
in the window that is both unique and warm enough -- and, when there is none,
which of the two it failed and whether codon choice can fix it.

That last part is the whole reason this lives in BT5 rather than in a cloning
tool. An arm is half backbone and half insert. When the ambiguity involves
insert bases the solver can break it by choosing different synonymous codons;
when it is backbone against backbone, nothing the optimizer does will help and
the honest output is a located warning. `Finding.fixable_by_codon_choice` is
that distinction, and it is what stops the solver chasing an unreachable
constraint into spurious infeasibility.

Melting temperatures here are nearest-neighbour (SantaLucia & Hicks 2004,
https://pubmed.ncbi.nlm.nih.gov/15139820/) and never travel without the salt and
strand conditions that produced them, for the same reason a fold energy never
travels without its parameter set: a Tm quoted without conditions cannot be
compared to anything, including a later run of BT5.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Literal

from Bio.SeqUtils import MeltingTemp  # noqa: N812

from bt5.core.types import Construct, Interval, reverse_complement
from bt5.vector.assemble import Assembly
from bt5.vector.kmers import VENDOR_UNIQUENESS_BP, ConstructKmerIndex, RepeatPair
from bt5.vector.notes import DesignNote

#: NEB's stated Gibson/NEBuilder overlap window and temperature floor.
#: vendor_asserted, last_verified 2026-08-27. Vendor numbers drift -- Twist moved
#: its homopolymer limit from 14 to 30 bp between 2023 and 2026 -- so this is a
#: parameter with a date on it, not a constant of nature.
MIN_ARM_BP = 15
MAX_ARM_BP = 40
MIN_ARM_TM_C = 48.0

JunctionName = Literal["5' backbone-insert", "3' insert-backbone"]


@dataclass(frozen=True, slots=True)
class TmConditions:
    """What a quoted melting temperature is a melting temperature AT.

    Defaults are the standard oligo-annealing conditions; they are the reported
    conditions, not a claim about the inside of an assembly reaction.
    """

    nn_table: str = "SantaLucia_Hicks_2004"
    na_mm: float = 50.0
    mg_mm: float = 0.0
    strand_nm: float = 250.0

    def describe(self) -> str:
        return (
            f"{self.nn_table}, Na+ {self.na_mm:g} mM, Mg2+ {self.mg_mm:g} mM, "
            f"strand {self.strand_nm:g} nM"
        )


DEFAULT_TM_CONDITIONS = TmConditions()


def melting_temperature(sequence: str, conditions: TmConditions = DEFAULT_TM_CONDITIONS) -> float:
    """Nearest-neighbour Tm. The 2004 table is pinned, never left to the default.

    Biopython's default is Allawi & SantaLucia 1997; leaving it unset would make
    every number here depend on a library default that can move under us.
    """
    # Biopython ships py.typed but leaves this constructor unannotated. The
    # ignore is confined to this one wrapper so no caller has to carry it.
    return float(
        MeltingTemp.Tm_NN(  # type: ignore[no-untyped-call]
            sequence,
            nn_table=MeltingTemp.DNA_NN4,
            Na=conditions.na_mm,
            Mg=conditions.mg_mm,
            dnac1=conditions.strand_nm,
            dnac2=conditions.strand_nm,
        )
    )


def gc_fraction(sequence: str) -> float:
    return (sequence.count("G") + sequence.count("C")) / len(sequence) if sequence else 0.0


@dataclass(frozen=True, slots=True)
class HomologyArm:
    """The overlap two fragments share at one junction, and where else it binds."""

    junction: str
    #: Centred on the boundary: half backbone, half insert, as a split overlap.
    interval: Interval
    sequence: str
    tm_c: float
    gc: float
    conditions: TmConditions = DEFAULT_TM_CONDITIONS
    #: Every OTHER place this exact sequence occurs, on either strand. A
    #: reverse-complement hit is a real mis-annealing site, not a curiosity:
    #: what anneals in Gibson is a single-stranded overhang.
    elsewhere: tuple[Interval, ...] = ()

    @property
    def length(self) -> int:
        return self.interval.length

    @property
    def unique(self) -> bool:
        return not self.elsewhere

    @property
    def warm_enough(self) -> bool:
        return self.tm_c >= MIN_ARM_TM_C

    @property
    def usable(self) -> bool:
        return self.unique and self.warm_enough


@dataclass(frozen=True, slots=True)
class JunctionPlan:
    """Per junction, the shortest usable arm -- or the best one and why it fails."""

    arms: tuple[HomologyArm, ...] = ()
    shared: tuple[tuple[str, str, Interval], ...] = ()
    insert_repeats: tuple[RepeatPair, ...] = ()
    notes: tuple[DesignNote, ...] = field(default_factory=tuple)

    @property
    def usable(self) -> bool:
        """True when every junction has an arm that is unique and warm enough."""
        return bool(self.arms) and all(a.usable for a in self.arms) and not self.shared


# -- geometry --------------------------------------------------------------


def junction_points(assembly: Assembly) -> tuple[tuple[str, int], ...]:
    """The two construct coordinates where designed sequence meets backbone."""
    cds = assembly.cds_interval
    return (("5' backbone-insert", cds.start), ("3' insert-backbone", cds.end))


def _arm_interval(construct: Construct, at: int, length: int) -> Interval | None:
    """An arm of `length` centred on the boundary at `at`.

    A split overlap is the standard layout, so half the arm comes from each
    fragment. On a circular construct an arm running off the front wraps into
    the one representation BT5 uses; on a linear one it simply does not exist.
    """
    n = construct.length
    start = at - length // 2
    if construct.is_circular:
        start %= n
    elif start < 0 or start + length > n:
        return None
    return Interval(start, start + length)


def occurrences(
    construct: Construct, needle: str, *, skip: Interval | None = None
) -> tuple[Interval, ...]:
    """Every place `needle` occurs in the construct, on either strand.

    Circular constructs are scanned on the doubled sequence so a match spanning
    the origin is found like any other; hits are folded back into construct
    coordinates. `skip` drops the arm's own position, which always matches.
    """
    n = construct.length
    text = construct.sequence * 2 if construct.is_circular else construct.sequence
    span = n if construct.is_circular else max(0, n - len(needle) + 1)
    rc = reverse_complement(needle)
    out: list[Interval] = []
    for i in range(span):
        window = text[i : i + len(needle)]
        if len(window) < len(needle):
            break
        strand = 1 if window == needle else (-1 if window == rc else 0)
        if strand == 0:
            continue
        if skip is not None and i == skip.start and strand == 1:
            continue
        out.append(Interval(i, i + len(needle), strand))  # type: ignore[arg-type]
    return tuple(out)


def build_arm(
    construct: Construct,
    junction: str,
    at: int,
    length: int,
    conditions: TmConditions = DEFAULT_TM_CONDITIONS,
) -> HomologyArm | None:
    """Measure the arm of a given length at a junction, or None if it cannot exist."""
    interval = _arm_interval(construct, at, length)
    if interval is None:
        return None
    sequence = construct.slice(interval)
    return HomologyArm(
        junction=junction,
        interval=interval,
        sequence=sequence,
        tm_c=melting_temperature(sequence, conditions),
        gc=gc_fraction(sequence),
        conditions=conditions,
        elsewhere=occurrences(construct, sequence, skip=interval),
    )


def shortest_usable_arm(
    construct: Construct,
    junction: str,
    at: int,
    *,
    min_bp: int = MIN_ARM_BP,
    max_bp: int = MAX_ARM_BP,
    conditions: TmConditions = DEFAULT_TM_CONDITIONS,
) -> HomologyArm | None:
    """The shortest arm in the window that is both unique and warm enough.

    Shortest rather than longest because every extra base is another chance to
    swallow a repeat.

    When nothing in the window satisfies both, the LONGEST arm that could be
    built comes back. That is not an arbitrary fallback: arms at one junction
    are NESTED, each centred on the same boundary, so the arm at length L is a
    substring of the arm at L+1. A duplicate of the longer arm therefore
    contains the shorter one too, which makes uniqueness monotone -- once an arm
    is unique, every longer arm is -- and Tm rises with length regardless. So
    the longest arm is simultaneously the warmest and, whenever any arm in the
    window is unique, a unique one. It is the best the window has to offer, and
    the temperature it reports is one the user could actually reach.

    `test_arms_at_one_junction_are_nested` is what guards that reasoning; if
    nesting ever stopped holding, this would silently start returning an
    ambiguous arm while a unique one existed.
    """
    longest: HomologyArm | None = None
    for length in range(min_bp, max_bp + 1):
        arm = build_arm(construct, junction, at, length, conditions)
        if arm is None:
            continue
        longest = arm
        if arm.usable:
            return arm
    return longest


def insert_shared_repeats(
    construct: Construct, insert: Interval, *, min_len: int = VENDOR_UNIQUENESS_BP
) -> tuple[RepeatPair, ...]:
    """Exact repeats of at least `min_len` shared between the insert and the rest.

    This is the uniqueness rule stated directly rather than inferred from the
    arms: a repeat anywhere between the insert and the backbone gives a fragment
    a second place to anneal, whether or not it happens to sit at a junction.
    Pairs with both copies inside the insert are somebody else's finding (they
    are a synthesis and stability problem, reported by the repeat scan); pairs
    with both copies in the backbone are not fixable by codon choice and are
    reported there too. What belongs here is the crossing pair.
    """
    index = ConstructKmerIndex.of(construct, min_len)
    out = [
        pair
        for pair in index.repeat_pairs(min_len, exclude=construct.exempt)
        if _inside(pair.first, insert) != _inside(pair.second, insert)
    ]
    return tuple(out)


def _inside(iv: Interval, outer: Interval) -> bool:
    return outer.start <= iv.start and iv.end <= outer.end


def plan_junctions(
    assembly: Assembly,
    *,
    min_bp: int = MIN_ARM_BP,
    max_bp: int = MAX_ARM_BP,
    conditions: TmConditions = DEFAULT_TM_CONDITIONS,
) -> JunctionPlan:
    """Check both junctions of an assembled construct, and say what is wrong."""
    construct = assembly.construct
    arms: list[HomologyArm] = []
    notes: list[DesignNote] = []

    for name, at in junction_points(assembly):
        arm = shortest_usable_arm(
            construct, name, at, min_bp=min_bp, max_bp=max_bp, conditions=conditions
        )
        if arm is None:
            notes.append(
                DesignNote(
                    kind="unavailable",
                    summary=(
                        f"the {name} junction is too close to the end of a linear "
                        f"construct for a {min_bp} bp overlap, so junction "
                        f"uniqueness could not be checked there"
                    ),
                    bears_on="assembly",
                )
            )
            continue
        arms.append(arm)
        notes.extend(_arm_notes(construct, arm, min_bp, max_bp))

    shared = _shared_between(construct, arms)
    for first, second, where in shared:
        notes.append(
            DesignNote(
                kind="liability",
                summary=(
                    f"the {first} and {second} arms share sequence, so a fragment "
                    f"can anneal at either junction and the insert can assemble "
                    f"in the wrong place or the wrong orientation"
                ),
                interval=where,
                bears_on="assembly",
                action="lengthen one overlap, or re-run so the solver diversifies the shared bases",
            )
        )

    repeats = insert_shared_repeats(construct, assembly.cds_interval, min_len=VENDOR_UNIQUENESS_BP)
    notes.extend(_repeat_notes(construct, repeats))

    return JunctionPlan(arms=tuple(arms), shared=shared, insert_repeats=repeats, notes=tuple(notes))


def _arm_notes(
    construct: Construct, arm: HomologyArm, min_bp: int, max_bp: int
) -> tuple[DesignNote, ...]:
    if arm.usable:
        return ()
    out: list[DesignNote] = []
    if not arm.unique:
        # An arm is CENTRED on the junction, so it always contains designed
        # bases. A duplicate has to match the insert half as well as the
        # backbone half, which means changing codons in that half breaks the
        # match. Junction ambiguity is therefore always the solver's to fix --
        # that is a consequence of the geometry, not an assumption, and it is
        # why this finding is worth a re-run rather than a shrug. (The solver
        # may still report infeasibility if the overlapping codons have no
        # synonyms at all, which is the designed answer to that.)
        out.append(
            DesignNote(
                kind="liability",
                summary=(
                    f"no overlap between {min_bp} and {max_bp} bp at the {arm.junction} "
                    f"junction is unique; at {arm.length} bp it still occurs "
                    f"{len(arm.elsewhere)} more time(s) in the construct, giving the "
                    f"fragment a second place to anneal. The overlap spans the "
                    f"designed CDS, so the solver can break it by choosing "
                    f"different codons"
                ),
                interval=arm.interval,
                bears_on="assembly",
                action="re-run: the ambiguity reaches into the designed CDS",
            )
        )
    if not arm.warm_enough:
        out.append(
            DesignNote(
                kind="liability",
                summary=(
                    f"the {arm.junction} junction is AT-rich: at {arm.length} bp the overlap "
                    f"melts at {arm.tm_c:.1f} C ({arm.gc:.0%} GC), below the "
                    f"{MIN_ARM_TM_C:.0f} C floor for a Gibson overlap "
                    f"({arm.conditions.describe()})"
                ),
                interval=arm.interval,
                bears_on="assembly",
                action=(
                    "use a longer overlap than the vendor window, or a different insertion point"
                ),
            )
        )
    return tuple(out)


def _shared_between(
    construct: Construct, arms: Sequence[HomologyArm]
) -> tuple[tuple[str, str, Interval], ...]:
    """Sequence common to two junctions' arms, which makes them interchangeable.

    Compared at the vendor uniqueness length rather than whole-arm, because a
    partial overlap is enough: two arms sharing their last 20 bases anneal to
    each other's partners just as happily as two identical arms.

    The window shrinks to fit the shorter arm. Comparing at a fixed 20 bp would
    make two IDENTICAL 15 bp arms -- the most ambiguous pair there is -- report
    as sharing nothing, because no 20 bp window fits inside either of them.
    """
    out: list[tuple[str, str, Interval]] = []
    for i, a in enumerate(arms):
        for b in arms[i + 1 :]:
            k = min(VENDOR_UNIQUENESS_BP, len(a.sequence), len(b.sequence))
            if k <= 0:
                continue
            windows: set[str] = set()
            for j in range(len(b.sequence) - k + 1):
                window = b.sequence[j : j + k]
                windows.add(window)
                windows.add(reverse_complement(window))
            for j in range(len(a.sequence) - k + 1):
                if a.sequence[j : j + k] in windows:
                    start = a.interval.start + j
                    out.append((a.junction, b.junction, Interval(start, start + k)))
                    break
    return tuple(out)


def _repeat_notes(construct: Construct, repeats: Sequence[RepeatPair]) -> tuple[DesignNote, ...]:
    return tuple(
        DesignNote(
            kind="liability",
            summary=(
                f"{pair.length} bp exact repeat shared between the designed insert and "
                f"the rest of the construct, {pair.spacer} bp apart; Gibson cannot tell "
                f"the two copies apart, and neither can the cell"
                + ("" if pair.reca_strain_helps else " -- and a recA- strain does not suppress it")
            ),
            interval=pair.first,
            bears_on="assembly",
            action="re-run: the insert half of the repeat is fixable by codon choice",
        )
        for pair in repeats
        if pair.risk != "low"
    )
