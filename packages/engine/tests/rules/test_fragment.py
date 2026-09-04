"""One designable span is one ordered fragment -- the assumption I7 rests on.

THIS FILE IS A TRIPWIRE, NOT A UNIT TEST OF `fragments()`. It exists because the
oracle's I7 checks GC per designable SPAN, which is only equivalent to checking
it per ordered TUBE while `fragments()` emits exactly one tube per span.

If the vendor lane ever splits a long span into several ordered fragments -- by a
length limit, or because an insert exceeds a product's ceiling -- a span-level
average can hide an out-of-band tube. A 3 kb span at 0.50 overall could be two
1.5 kb tubes at 0.30 and 0.70, and I7 would pass it while both tubes sit outside
the band. The oracle would then be WEAKER than the rule it backstops, which
inverts the point of having an independent validator: E2 evaluates real
`Fragment` objects and would catch it; I7 would not.

AND THE ORACLE CANNOT NOTICE THAT CHANGE, BY CONSTRUCTION.
`tests/data_integrity/test_oracle_independence.py` forbids `verify.py` from
importing anything under `bt5.rules` -- by AST walk and again by raw substring on
the file text. So `verify.py` can never learn that `fragments()` started
splitting. There is no mechanism by which this breaks loudly, which is why the
assumption needs a test on the side that would actually make the change.

The assumption is already written down in `verify.py`'s I7 block as a falsifiable
claim rather than a definition. This is its other half, in the lane that owns
`fragments()`. See #64, filed from #62.

WHEN THIS FILE FAILS, DO NOT "FIX" IT BY UPDATING THE EXPECTED COUNT. A failure
here means splitting shipped and I7 is now weaker than E2. The fix is to pass the
split points into the oracle the way vendor facts already arrive -- as arguments
to `verify_construct`, never as an import -- so I7 checks tubes rather than spans.
"""

from __future__ import annotations

from bt5.core.types import Construct, Interval, Segment, SegmentKind, Topology
from bt5.rules.fragment import NO_ADAPTERS, fragments
from bt5.rules.vendors import PROFILES, orderable_keys


def largest_orderable_ceiling() -> int:
    """The longest single order any shipped profile accepts.

    Derived from the registry rather than hardcoded, so that a vendor raising or
    lowering a ceiling cannot quietly make the span below stop being oversized --
    which would leave this file passing while testing nothing.
    """
    ceilings = [
        PROFILES[key].length_bp[1]
        for key in orderable_keys()
        if PROFILES[key].length_bp is not None
    ]
    assert ceilings, "no orderable profile declares a length range; this probe is blind"
    return max(ceilings)


def construct_with_spans(*lengths: int) -> Construct:
    """A construct whose designable spans have the given lengths, separated by
    one backbone base each so they cannot merge into a single interval."""
    parts: list[str] = []
    segments: list[Segment] = []
    cursor = 0
    for i, length in enumerate(lengths):
        if i:
            parts.append("T")
            cursor += 1
        # Mid-GC so the span itself is unremarkable; this file is about COUNTS.
        body = ("ATGC" * (length // 4 + 1))[:length]
        parts.append(body)
        segments.append(
            Segment(Interval(cursor, cursor + length), SegmentKind.DESIGNABLE_CDS, f"cds{i}")
        )
        cursor += length
    return Construct(sequence="".join(parts), topology=Topology.LINEAR, segments=tuple(segments))


class TestOneSpanIsOneFragment:
    def test_a_span_longer_than_every_vendor_ceiling_is_still_one_fragment(self) -> None:
        """THE tripwire. This is the exact case a length-limit split would
        target, so it is the one that must be pinned."""
        ceiling = largest_orderable_ceiling()
        oversized = ceiling * 2
        c = construct_with_spans(oversized)

        frags = fragments(c, NO_ADAPTERS)

        assert len(frags) == len(c.editable) == 1, (
            f"a {oversized} bp span (every vendor ceiling is <= {ceiling}) produced "
            f"{len(frags)} fragments. If the vendor lane now splits by length, the "
            f"oracle's I7 checks spans while E2 checks tubes, and I7 is the weaker "
            f"of the two. See this file's docstring and #64 -- do not update this count."
        )
        assert len(frags[0].construct.sequence) == oversized

    def test_the_count_matches_the_span_count_for_several_spans(self) -> None:
        c = construct_with_spans(120, 240, 360)
        frags = fragments(c, NO_ADAPTERS)
        assert len(frags) == len(c.editable) == 3

    def test_each_fragment_carries_its_own_span_as_its_origin(self) -> None:
        """One-per-span is only meaningful if the mapping is the identity: a
        count that matched while the origins were scrambled would still break
        I7's span/tube equivalence."""
        c = construct_with_spans(120, 240, 360)
        frags = fragments(c, NO_ADAPTERS)
        assert [f.origin for f in frags] == sorted(c.editable)
        for f in frags:
            assert len(f.construct.sequence) == f.origin.end - f.origin.start

    def test_two_spans_are_two_molecules_not_one_concatenation(self) -> None:
        """The reason one-per-span is right in the first place: two designable
        spans are two synthesis reactions that never meet until after assembly,
        so a repeat shared BETWEEN them is not a synthesis liability for either.
        A concatenating implementation would pass the count test above only if it
        also collapsed the spans, so assert the sequences are separate."""
        c = construct_with_spans(60, 60)
        a, b = fragments(c, NO_ADAPTERS)
        assert len(a.construct.sequence) == 60
        assert len(b.construct.sequence) == 60
        assert a.origin != b.origin
