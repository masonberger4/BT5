"""Findings that depend on the construct alone.

Everything here takes a `Construct` and nothing else, which is what lets the
same checks serve two callers that otherwise have nothing in common: a designed
construct on its way to an annotated export, and a plasmid a user handed over
with no CDS to design at all. Repeats destabilise a plasmid whether or not BT5
wrote any of it, so a survey path that re-implemented them would drift from the
design path exactly where it mattered least to notice.

Nothing here predicts an outcome. Each finding says what is present, where, and
what a user can do about it -- including, for the several cases where nothing
can be done by codon choice, saying that plainly.
"""

from __future__ import annotations

from collections.abc import Sequence

from bt5.core.types import Construct, Interval, SegmentKind
from bt5.vector.kmers import ConstructKmerIndex, RepeatPair
from bt5.vector.notes import DesignNote, format_span

#: Shortest exact repeat worth reporting; the vendor uniqueness floor.
MIN_REPEAT_BP = 20

#: A stem shorter than this is ubiquitous in any sequence and says nothing.
#: A cruciform needs a substantial one, and BT5 reports geometry rather than a
#: probability, so this is a reporting floor and not a risk threshold.
MIN_STEM_BP = 20

#: Beyond this the two arms are too far apart to be read as one hairpin.
MAX_LOOP_BP = 60


def identical_exempt_repeats(construct: Construct) -> tuple[DesignNote, ...]:
    """Report scan-exempt regions that are exact copies of each other.

    LTRs and ITRs are whitelisted precisely because they violate the repeat rules
    by construction, and no codon choice can change them -- so the deliverable is
    a strain and temperature protocol, not a redesign. Saying nothing would leave
    the user to discover it as a deletion after a maxiprep.
    """
    exempt = [s for s in construct.segments if s.kind is SegmentKind.WHITELISTED_REPEAT]
    out: list[DesignNote] = []
    for i, first in enumerate(exempt):
        for second in exempt[i + 1 :]:
            if construct.slice(first.interval) != construct.slice(second.interval):
                continue
            out.append(
                DesignNote(
                    kind="liability",
                    summary=(
                        f"{first.label or 'region'} and {second.label or 'region'} are "
                        f"identical over {first.interval.length} bp "
                        f"({format_span(second.interval, construct.length)})"
                    ),
                    interval=first.interval,
                    bears_on="plasmid stability",
                    action=(
                        "propagate in a recA- strain such as Stbl3 at 30 C; note that "
                        "this suppresses the RecA-dependent pathway only, which is the "
                        "one that acts on repeats this long"
                    ),
                )
            )
    return tuple(out)


def repeat_liabilities(
    construct: Construct, *, min_len: int = MIN_REPEAT_BP
) -> tuple[DesignNote, ...]:
    """Exact repeats anywhere in the assembled construct, clustered by region.

    Scoped to the whole construct on purpose: a repeat between the insert and the
    backbone is exactly as destabilising as one inside either, and it is only
    visible once they are assembled. Whitelisted ITR/LTR regions are excluded --
    they are an accepted design feature reported separately, not a finding.

    Overlapping pairs are CLUSTERED. A tandem array such as a Gal4 UAS or a TetO
    operator block yields one exact-repeat pair per offset -- fourteen of them on
    a real BiTE vector -- and fourteen near-identical notes bury the finding they
    are trying to deliver. One note per repetitive region, carrying the worst
    risk and the longest identity in it, is what a user can act on.
    """
    pairs = [
        p
        for p in ConstructKmerIndex.of(construct, min_len).repeat_pairs(
            min_len, exclude=construct.exempt
        )
        if p.risk != "low"
    ]
    out: list[DesignNote] = []
    for cluster in _cluster(pairs):
        worst = max(cluster, key=lambda p: (p.risk == "high", p.length))
        longest = max(p.length for p in cluster)
        lo = min(p.first.start for p in cluster)
        hi = max(p.first.end for p in cluster)
        tandem = any(p.tandem for p in cluster)
        if len(cluster) == 1:
            where = "in tandem" if worst.tandem else f"{worst.spacer} bp apart"
            summary = (
                f"{worst.length} bp exact repeat, {where} "
                f"(copies at {worst.first.start + 1} and {worst.second.start + 1})"
            )
        else:
            summary = (
                f"repetitive region of {hi - lo} bp containing {len(cluster)} exact "
                f"repeats up to {longest} bp" + (", some in tandem" if tandem else "")
            )
        out.append(
            DesignNote(
                kind="liability",
                summary=f"{summary}; {worst.risk} risk",
                interval=Interval(lo, hi),
                bears_on="plasmid stability and DNA synthesis",
                action=(
                    "a recA- strain such as Stbl3 suppresses this length class"
                    if worst.reca_strain_helps
                    else "a recA- strain does NOT help here: below about 200 bp, "
                    "deletion is RecA-independent"
                ),
            )
        )
    return tuple(out)


def _cluster(pairs: Sequence[RepeatPair]) -> list[list[RepeatPair]]:
    """Group pairs whose copies fall in the same stretch of sequence."""
    clusters: list[list[RepeatPair]] = []
    for pair in sorted(pairs, key=lambda p: p.first.start):
        # Cluster on the FIRST copy only. Two copies far apart are two regions;
        # spanning them would report the whole distance between them as
        # "repetitive", which on a lentiviral vector means the entire payload.
        span = pair.first
        for cluster in clusters:
            lo = min(p.first.start for p in cluster)
            hi = max(p.first.end for p in cluster)
            if span.start < hi and lo < span.end:
                cluster.append(pair)
                break
        else:
            clusters.append([pair])
    return clusters


def hairpin_liabilities(
    construct: Construct, *, min_stem: int = MIN_STEM_BP, max_loop: int = MAX_LOOP_BP
) -> tuple[DesignNote, ...]:
    """Inverted repeats long enough to extrude, reported with their geometry.

    Kept separate from the direct-repeat notes because the ACTION differs, which
    is the only reason a user cares which kind they have. A direct repeat is a
    deletion risk and codon choice fixes it. A stem-loop stalls replication
    forks, foul solid-phase synthesis, and is cleaved by SbcCD -- so the answer
    is a strain and a temperature, and telling someone to re-run the optimizer
    would waste their afternoon.

    Whitelisted regions are excluded: an AAV ITR IS a palindrome, by design, and
    reporting it as a liability every run trains the reader to skip the section.
    """
    stems = ConstructKmerIndex.of(construct, min_stem).inverted_repeats(
        min_stem, max_loop, exclude=construct.exempt
    )
    return tuple(
        DesignNote(
            kind="liability",
            summary=(
                f"{stem.stem} bp inverted repeat"
                + (
                    " with no loop -- a perfect palindrome"
                    if stem.perfect_palindrome
                    else f" separated by {stem.loop} bp"
                )
                + f" (second arm at {format_span(stem.second, construct.length)})"
            ),
            interval=stem.first,
            bears_on="plasmid stability and DNA synthesis",
            action=(
                "propagate at 30 C in a strain lacking sbcC; this is a "
                "hairpin, so re-running the optimizer will not remove it"
            ),
        )
        for stem in stems
    )
