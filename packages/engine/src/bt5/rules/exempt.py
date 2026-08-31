"""Is this finding really the user's ITR again?

Shared by the repeat rules, and outside `catalog/` on purpose: autodiscovery
walks `bt5.rules.catalog` only, so a helper module there would be imported as a
rule and fail the contract test for having no citations.

`SegmentKind.WHITELISTED_REPEAT` exists because lentiviral LTRs are long perfect
direct repeats and AAV ITRs are 145 bp palindromes -- both violate the repeat
rules BY CONSTRUCTION, and the answer is a strain and a temperature, not a
redesign. Reporting them as defects buries every finding a user can act on
beneath the two they cannot.

The subtlety is that strict containment does not work here, and finding that out
is what this module is for. The repeat scans report MAXIMAL matches, so a stem
grows outward while the flanking bases happen to pair -- which they do by chance
roughly a quarter of the time per base. A 145 bp ITR palindrome therefore
reports as a 147 or 149 bp stem, spills a base or two past the annotated
feature, and a containment test says "not exempt" on the single most exempt
object in the file.

So coverage, not containment: a finding is the user's annotated feature when
almost all of it is. The threshold is deliberately high -- a genuinely new
repeat that merely abuts an ITR shares far less than 90% of its length with it.

This module also carries the wrap-aware PAIR GEOMETRY the repeat rules share
(`pair_span`). Same reason it holds the exemption helpers: both F1 and F3 ask
the same question of a pair of intervals on a circle, and answering it twice is
how the two answers drift apart.
"""

from __future__ import annotations

from collections.abc import Sequence

from bt5.core.types import Construct, Interval

#: How much of an arm must lie inside one exempt region before the finding is
#: that region rather than a design problem. High enough that abutting a feature
#: is not enough; loose enough to absorb a few bases of chance extension.
EXEMPT_COVERAGE = 0.9


def overlap_length(a: Interval, b: Interval, construct_length: int, circular: bool) -> int:
    """Bases shared by two intervals, wrap-aware. 0 when they are disjoint."""
    shifts = (0, construct_length, -construct_length) if circular else (0,)
    best = 0
    for shift in shifts:
        overlap = min(a.end, b.end + shift) - max(a.start, b.start + shift)
        best = max(best, overlap)
    return best


def coverage(
    iv: Interval, regions: Sequence[Interval], construct_length: int, circular: bool
) -> float:
    """The largest fraction of `iv` covered by any ONE region.

    One region rather than the union, because the question is "is this finding
    that annotated feature?" -- an arm straddling two different exempt features
    is a finding about the junction between them, which is real.
    """
    if iv.length <= 0 or not regions:
        return 0.0
    best = max(overlap_length(iv, r, construct_length, circular) for r in regions)
    return best / iv.length


def pair_span(first: Interval, second: Interval, construct_length: int, circular: bool) -> Interval:
    """The arc bounding a repeat pair -- the SAME arc its separation was measured on.

    Two copies on a circle bound two arcs, and the shorter one may be the one
    through the origin. F1's `_spacer` and F3's `_loop` both already measure the
    short way round; the breach interval did not, so a pair at 20 and 700 on an
    800 bp plasmid reported "80 bp apart" while its interval covered the other
    720 bp. Message and interval described opposite arcs.

    That is not cosmetic. Both rules declare `LocalizationPolicy.PAIRED_SEGMENTS`,
    and the solver hands this interval straight to the repair window -- so the
    wrong arc aims repair at the wrong 90% of the plasmid.

    Lives here rather than in either rule because the bug was a duplication bug:
    the same wrong line was written twice, and copying the fix twice would
    preserve the conditions for the two to drift apart again.
    """
    forward = Interval(first.start, max(first.end, second.end))
    if not circular or second.end > construct_length:
        # Linear, or a pair the k-mer index already emitted in wrapped form --
        # `forward` is then already the wrapping representation.
        return forward
    gap = second.start - first.end
    wrap = construct_length - (second.end - first.start)
    if 0 <= wrap < gap:
        # `end` past the construct length is the one representation `Interval`
        # defines for a wrapping interval (core/types.py:62).
        return Interval(second.start, first.end + construct_length)
    return forward


def both_arms_exempt(
    c: Construct,
    first: Interval,
    second: Interval,
    *,
    threshold: float = EXEMPT_COVERAGE,
) -> bool:
    """True when BOTH copies are substantially inside scan-exempt regions.

    Both, not either. A designed repeat that happens to match part of an ITR is
    still the design's problem, and exempting it because one end landed on a
    whitelisted feature would hide exactly the finding the user can fix.
    """
    regions = c.exempt
    if not regions:
        return False
    return all(
        coverage(iv, regions, c.length, c.is_circular) >= threshold for iv in (first, second)
    )
