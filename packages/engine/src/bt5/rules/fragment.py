"""What the vendor actually builds.

Section 2.E of the brief carries its scope in its own header: manufacturability
is "evaluated on the synthesized fragment + adapters". That scope is not a
detail, it is the entire reason the E-series repeat rules are not duplicates of
the F-series ones, and getting it wrong would make BT5 report one 22 bp repeat
four times in one panel.

The two scopes differ in three ways that all matter:

**Extent.** The vendor synthesises the insert. They never see the backbone, so a
repeat shared between the insert and a promoter 4 kb away is a plasmid-stability
finding (F1's) and not a synthesis finding. Conversely a repeat wholly inside
the user's backbone is neither -- nobody is building that DNA.

**Topology.** The ordered fragment is a LINEAR molecule. The plasmid is
circular. A repeat spanning the origin is real for F1 and does not exist for the
vendor, because the two halves are never in the same tube. So the fragment is
built as a linear construct and the E-series rules legitimately disagree with
the F-series at position 0.

**Adapters.** Twist Gene Fragments ship with fixed flanking sequence, and it is
synthesised as part of the fragment. An insert whose 5' end happens to reproduce
the adapter is a real assembly failure that no whole-plasmid scan can see,
because the adapter is not in the plasmid at all.

Adapters are modelled as BACKBONE segments of the fragment construct. That is
what makes `Breach.fixable_by_codon_choice` come out right for free: recoding
can move the insert off a collision with an adapter, but nothing in BT5 can
change the adapter itself.
"""

from __future__ import annotations

from dataclasses import dataclass

from bt5.core.types import (
    Construct,
    Interval,
    Segment,
    SegmentKind,
    Topology,
)

#: Twist ADAPTER-ON Gene Fragment adapters, synthesised as part of the fragment.
#:
#: These belong to an OPTION, not to the product. Twist states plainly that
#: "adapter sequences are not added by default to Gene Fragments" -- adapter-on
#: and adapter-free are two choices made at checkout, so a plain Twist Gene
#: Fragment order carries none of this:
#: https://www.twistbioscience.com/faq/gene-synthesis/are-adapter-sequences-appended-ends-my-sequences
#:
#: VENDOR_ASSERTED and dated twice over: these are the adapters for orders placed
#: after 2021-11-02, so an older order carries a DIFFERENT pair. Exactly the kind
#: of constant that goes stale silently.
#: https://www.twistbioscience.com/faq/gene-synthesis/what-are-adapter-sequences-used-adapter-gene-fragments
TWIST_FIVE_PRIME = "CAATCCGCCCTCACTACAACCG"
TWIST_THREE_PRIME = "CTACTCTGGCGTCGATGAGGGA"


@dataclass(frozen=True, slots=True)
class Adapters:
    """Fixed sequence the vendor adds to every ordered fragment."""

    five: str = ""
    three: str = ""
    vendor: str = "none"

    @property
    def total(self) -> int:
        return len(self.five) + len(self.three)


NO_ADAPTERS = Adapters()

#: A plain Twist Gene Fragment order. No adapters -- but it still names the
#: vendor, because a finding has to say WHOSE fragment it is about and
#: `NO_ADAPTERS` would report the vendor as "none".
TWIST_GENE_FRAGMENT = Adapters(vendor="twist_gene_fragment")

#: The adapter-on option. Only this one carries the adapters.
TWIST_ADAPTER_ON = Adapters(TWIST_FIVE_PRIME, TWIST_THREE_PRIME, "twist_gene_fragment_adapter_on")

#: IDT ships gene fragments without adapters, so these carry only a vendor name.
IDT_EBLOCKS = Adapters(vendor="idt_eblocks")
IDT_GBLOCKS = Adapters(vendor="idt_gblocks")

# The registry that binds these to lengths, run limits and GC bands lives in
# `vendors.py`, which imports FROM here. Adapters are a molecular fact about the
# synthesised molecule and belong beside `fragments()`; which configurations
# exist and what else is true of them is a catalogue, and keeping the catalogue
# in a second module here is what produced the split default in the first place.


@dataclass(frozen=True, slots=True)
class Fragment:
    """One tube of synthetic DNA, and where it came from.

    `construct` is the fragment as a rule sees it: linear, with the ordered DNA
    designable and the adapters backbone. Rules evaluate against THAT, never
    against `sequence` as a bare string -- the same discipline the whole project
    rests on, applied one level down.
    """

    #: adapter5 + ordered DNA + adapter3.
    sequence: str
    #: Where the ordered DNA sits inside `sequence`. Equals [0, len) with no adapters.
    ordered: Interval
    #: Where the ordered DNA sits in the PARENT construct. May wrap the origin;
    #: the fragment itself never does, because it is a linear molecule.
    origin: Interval
    #: The fragment as a linear, adapter-annotated construct.
    construct: Construct
    #: Which adapter pair was spliced on, for the report.
    vendor: str

    def to_construct(self, iv: Interval) -> Interval | None:
        """Map a fragment-space interval back to parent-construct coordinates.

        Returns None for an interval lying WHOLLY inside a vendor adapter: that
        is the vendor's own sequence colliding with itself, which is their
        problem and not a finding about this design. An interval straddling the
        adapter/insert boundary IS the design's problem -- the insert end
        reproduces the adapter -- and is clamped to the ordered part so it
        carries a real parent coordinate.
        """
        lo = max(iv.start, self.ordered.start)
        hi = min(iv.end, self.ordered.end)
        if hi <= lo:
            return None
        delta = self.origin.start - self.ordered.start
        return Interval(lo + delta, hi + delta, iv.strand)

    def touches_adapter(self, iv: Interval) -> bool:
        return iv.start < self.ordered.start or iv.end > self.ordered.end


def fragments(c: Construct, adapters: Adapters = NO_ADAPTERS) -> list[Fragment]:
    """The tubes BT5 would order for this construct, one per designable span.

    One per span, not one concatenation of all of them. Two designable spans are
    two separate synthesis reactions, so a repeat shared BETWEEN them is not a
    synthesis liability for either -- the two molecules never meet until after
    assembly. (It is still a plasmid liability, and F1 reports it.)
    """
    out: list[Fragment] = []
    for iv in sorted(c.editable):
        ordered_seq = c.slice(iv)
        if not ordered_seq:
            continue
        seq = adapters.five + ordered_seq + adapters.three
        start = len(adapters.five)
        ordered = Interval(start, start + len(ordered_seq))

        segments = [Segment(ordered, SegmentKind.DESIGNABLE_CDS, "ordered")]
        if adapters.five:
            segments.append(Segment(Interval(0, start), SegmentKind.BACKBONE, "5' adapter"))
        if adapters.three:
            segments.append(
                Segment(
                    Interval(ordered.end, len(seq)),
                    SegmentKind.BACKBONE,
                    "3' adapter",
                )
            )
        out.append(
            Fragment(
                sequence=seq,
                ordered=ordered,
                origin=iv,
                construct=Construct(
                    sequence=seq,
                    topology=Topology.LINEAR,
                    segments=tuple(segments),
                ),
                vendor=adapters.vendor,
            )
        )
    return out
