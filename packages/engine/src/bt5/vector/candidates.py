"""Finding the insertion site on a vector somebody else annotated.

"Take the one CDS feature" does not survive contact with real maps. Across three
Addgene transfer vectors the transgene was never annotated as a CDS: what carried
the `CDS` key was a Myc tag, a TM domain, a four-codon Factor Xa site inside the
WPRE, a P2A peptide, and the bacterial markers. Annotation quality depends
entirely on whoever deposited the plasmid, so the detector has to work from the
sequence and treat annotation as evidence rather than instruction.

So: candidates are ORFs on both strands, scored by the context around them, with
every contributing signal recorded. The score is never the whole answer -- the
app shows a ranked list with reasons and the user confirms -- but it has to put
the transgene above the selection marker, which is exactly what a naive
longest-ORF rule fails to do. In an empty backbone the longest ORF in the
payload region runs from the cloning site straight through P2A into PuroR, so
without a marker penalty the tool would offer to codon-optimise the puromycin
resistance gene.

No signal here is a published rule and none is claimed as one; this is a
convenience for locating a site, and the user confirms it before anything is
designed.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal

from Bio.Data import CodonTable

from bt5.core.types import Feature, Interval, SegmentKind, Strand, reverse_complement
from bt5.vector.backbone import (
    DEFAULT_EXEMPT_KINDS,
    DEFAULT_EXEMPT_LABELS,
    InsertionSite,
    VectorBackbone,
)
from bt5.vector.markers import MIN_EXPRESSION_PROMOTER_BP, is_marker

CandidateKind = Literal["orf", "cloning_site"]

#: Scoring signals. Integers, not tuned weights: each is a yes/no piece of
#: evidence and the reasons travel with the score so a user can see why a
#: candidate ranked where it did.
SCORE_STARTS_AT_KOZAK = 3
SCORE_STARTS_AT_SIGNAL_PEPTIDE = 3
SCORE_INSIDE_PAYLOAD = 2
SCORE_AFTER_PROMOTER = 2
SCORE_BEFORE_POLYA = 2
SCORE_OVERLAPS_ANNOTATED_CDS = 2
SCORE_PER_100_AA = 1
SCORE_MAX_LENGTH_BONUS = 6
PENALTY_MARKER = -8
PENALTY_ORIGIN = -8

#: Below this a candidate is reported but not offered as a confident pick. An
#: empty backbone should land here: its best ORF is the selection cassette.
CONFIDENT_SCORE = 8

#: A multiple cloning site is tens to a couple of hundred bases. A larger
#: unannotated gap is not a cloning site, it is just unannotated sequence,
#: and offering it as somewhere to put a gene would be a guess dressed up.
MAX_CLONING_SITE_BP = 500

PROMOTER_SEARCH_BP = 2000
POLYA_SEARCH_BP = 5000


@dataclass(frozen=True, slots=True)
class SiteCandidate:
    """A possible place for the designed CDS, with its evidence."""

    interval: Interval
    kind: CandidateKind
    score: int
    reasons: tuple[str, ...]
    protein_length: int = 0
    label: str = ""

    @property
    def confident(self) -> bool:
        return self.score >= CONFIDENT_SCORE

    def as_site(self, *, table_id: int | None = None) -> InsertionSite:
        return InsertionSite(
            interval=self.interval,
            label=self.label or "candidate",
            source="explicit",
            detected_table_id=table_id,
        )


def _feature_text(backbone: VectorBackbone, feature: Feature) -> str:
    parts = [backbone.label_of(feature)]
    for key in ("note", "product", "gene"):
        parts.extend(feature.qualifiers.get(key, ()))
    return " ".join(parts)


def find_orfs(backbone: VectorBackbone, *, table_id: int, min_aa: int = 80) -> tuple[Interval, ...]:
    """Every ORF at least `min_aa` long, on both strands, origin-spanning included.

    Only the longest ORF ending at each stop is kept: a downstream in-frame ATG
    describes the same gene, and reporting both as separate candidates is noise.
    """
    table = CodonTable.unambiguous_dna_by_id[table_id]
    stops = set(table.stop_codons)
    n = backbone.length
    out: list[Interval] = []
    for strand in (1, -1):
        text = backbone.sequence if strand == 1 else reverse_complement(backbone.sequence)
        scan = text + text[: 2 * n] if backbone.is_circular else text
        best: dict[int, int] = {}  # stop position -> earliest start
        for start in range(n if backbone.is_circular else max(0, len(scan) - 2)):
            if scan[start : start + 3] != "ATG":
                continue
            i = start
            while i + 3 <= len(scan):
                if scan[i : i + 3] in stops:
                    if (i - start) // 3 >= min_aa and best.get(i, start + 1) > start:
                        best[i] = start
                    break
                i += 3
        for stop, start in best.items():
            end = stop + 3
            out.append(_to_genomic(start, end, strand, n))
    return tuple(out)


def _to_genomic(start: int, end: int, strand: Strand, length: int) -> Interval:
    """Map an ORF found on a scanning strand back to construct coordinates."""
    if strand == 1:
        return Interval(start % length, (start % length) + (end - start), 1)
    # The minus-strand scan runs along the reverse complement, so a hit at
    # [start, end) there covers [n - end, n - start) genomically -- and the
    # feature's 5' end is at the HIGH coordinate.
    lo = (length - end) % length
    return Interval(lo, lo + (end - start), -1)


def score_candidate(
    backbone: VectorBackbone, orf: Interval, *, table_id: int
) -> tuple[int, tuple[str, ...]]:
    """Score an ORF by the context around it, returning the reasons as well."""
    score = 0
    reasons: list[str] = []
    n = backbone.length
    covered = {p % n for p in range(orf.start, orf.end)}
    five_prime = orf.start if orf.strand == 1 else orf.end

    aa = orf.length // 3
    bonus = min(SCORE_MAX_LENGTH_BONUS, (aa // 100) * SCORE_PER_100_AA)
    score += bonus
    reasons.append(f"+{bonus} open reading frame of {aa} aa")

    for feature in backbone.features:
        text = _feature_text(backbone, feature)
        overlaps = any(
            p % n in covered for p in range(feature.interval.start, feature.interval.end)
        )
        kind = feature.kind.lower()

        if overlaps and kind == "cds" and is_marker(text):
            score += PENALTY_MARKER
            reasons.append(
                f"{PENALTY_MARKER} runs through {backbone.label_of(feature)!r}, a selection marker"
            )
        elif overlaps and kind == "cds":
            score += SCORE_OVERLAPS_ANNOTATED_CDS
            reasons.append(
                f"+{SCORE_OVERLAPS_ANNOTATED_CDS} contains the annotated CDS "
                f"{backbone.label_of(feature)!r}"
            )
        elif overlaps and kind == "rep_origin":
            score += PENALTY_ORIGIN
            reasons.append(f"{PENALTY_ORIGIN} runs through an origin of replication")

        if (
            kind == "regulatory"
            and "kozak" in text.casefold()
            and _touches(feature.interval, five_prime)
        ):
            score += SCORE_STARTS_AT_KOZAK
            reasons.append(f"+{SCORE_STARTS_AT_KOZAK} starts at an annotated Kozak sequence")
        if kind == "sig_peptide" and _starts_at(feature.interval, five_prime, orf.strand):
            score += SCORE_STARTS_AT_SIGNAL_PEPTIDE
            reasons.append(
                f"+{SCORE_STARTS_AT_SIGNAL_PEPTIDE} starts at the annotated signal peptide "
                f"{backbone.label_of(feature)!r}"
            )

    if _between_exempt_repeats(backbone, orf):
        score += SCORE_INSIDE_PAYLOAD
        reasons.append(f"+{SCORE_INSIDE_PAYLOAD} lies inside the ITR/LTR payload region")
    promoter = _upstream_promoter(backbone, orf)
    if promoter is not None:
        score += SCORE_AFTER_PROMOTER
        reasons.append(f"+{SCORE_AFTER_PROMOTER} sits downstream of {promoter!r}")
    polya = _downstream_polya(backbone, orf)
    if polya is not None:
        score += SCORE_BEFORE_POLYA
        reasons.append(f"+{SCORE_BEFORE_POLYA} sits upstream of {polya!r}")

    return score, tuple(reasons)


def _touches(iv: Interval, position: int) -> bool:
    """The feature brackets the start codon, within a codon either side."""
    return iv.start - 3 <= position <= iv.end + 3


def _starts_at(iv: Interval, position: int, strand: Strand) -> bool:
    edge = iv.start if strand == 1 else iv.end
    return abs(edge - position) <= 3


def _between_exempt_repeats(backbone: VectorBackbone, orf: Interval) -> bool:
    """AAV and lentiviral payloads sit between the ITR or LTR pair; markers do not."""
    repeats = [
        f.interval
        for f in backbone.features
        if DEFAULT_EXEMPT_KINDS.get(f.kind.lower()) is SegmentKind.WHITELISTED_REPEAT
        or any(t in backbone.label_of(f).casefold() for t in DEFAULT_EXEMPT_LABELS)
    ]
    if len(repeats) < 2:
        return False
    lo = min(r.end for r in repeats)
    hi = max(r.start for r in repeats)
    return lo <= orf.start and orf.end <= hi


def _upstream_promoter(backbone: VectorBackbone, orf: Interval) -> str | None:
    """The nearest promoter that is not driving a selection marker."""
    best: tuple[int, str] | None = None
    for feature in backbone.features:
        if feature.kind.lower() != "promoter":
            continue
        if is_marker(_feature_text(backbone, feature)):
            continue
        if feature.interval.strand != orf.strand:
            continue
        gap = (
            orf.start - feature.interval.end
            if orf.strand == 1
            else feature.interval.start - orf.end
        )
        if 0 <= gap <= PROMOTER_SEARCH_BP and (best is None or gap < best[0]):
            best = (gap, backbone.label_of(feature))
    return best[1] if best else None


def _downstream_polya(backbone: VectorBackbone, orf: Interval) -> str | None:
    best: tuple[int, str] | None = None
    for feature in backbone.features:
        if feature.kind.lower() not in ("polya_signal", "polya_site"):
            continue
        gap = (
            feature.interval.start - orf.end
            if orf.strand == 1
            else orf.start - feature.interval.end
        )
        if 0 <= gap <= POLYA_SEARCH_BP and (best is None or gap < best[0]):
            best = (gap, backbone.label_of(feature))
    return best[1] if best else None


def cloning_sites(backbone: VectorBackbone) -> tuple[SiteCandidate, ...]:
    """Unannotated gaps just downstream of a non-marker promoter.

    This is the answer for an empty backbone, where the honest output is not a
    transgene but "here is where one goes".
    """
    out: list[SiteCandidate] = []
    occupied = sorted(
        (f.interval.start, f.interval.end)
        for f in backbone.features
        if f.kind.lower() not in ("source", "primer_bind", "protein_bind", "regulatory")
    )
    for feature in backbone.features:
        if feature.kind.lower() != "promoter" or is_marker(_feature_text(backbone, feature)):
            continue
        if feature.interval.strand != 1:
            continue
        if feature.interval.length < MIN_EXPRESSION_PROMOTER_BP:
            # T7/T3/SP6 are ~19 bp and drive sequencing, not expression.
            continue
        start = feature.interval.end
        nxt = min((s for s, _ in occupied if s >= start), default=backbone.length)
        if not 20 <= nxt - start <= MAX_CLONING_SITE_BP:
            continue
        out.append(
            SiteCandidate(
                interval=Interval(start, nxt, 1),
                kind="cloning_site",
                score=0,
                reasons=(
                    f"unannotated {nxt - start} bp immediately downstream of "
                    f"{backbone.label_of(feature)!r}",
                ),
                label=f"cloning site after {backbone.label_of(feature)}",
            )
        )
    out.sort(key=lambda c: -c.interval.length)
    return tuple(out)


def suggest_insertion_sites(
    backbone: VectorBackbone, *, table_id: int, min_aa: int = 80, limit: int = 5
) -> tuple[SiteCandidate, ...]:
    """Ranked places the designed CDS could go, best first, with reasons."""
    scored: list[SiteCandidate] = []
    for orf in find_orfs(backbone, table_id=table_id, min_aa=min_aa):
        score, reasons = score_candidate(backbone, orf, table_id=table_id)
        scored.append(
            SiteCandidate(
                interval=orf,
                kind="orf",
                score=score,
                reasons=reasons,
                protein_length=orf.length // 3 - 1,
                label=f"ORF {orf.start + 1}..{orf.end}",
            )
        )
    scored.sort(key=lambda c: (-c.score, -c.protein_length, c.interval.start))
    top: Sequence[SiteCandidate] = scored[:limit]
    if not top or not top[0].confident:
        return (*top, *cloning_sites(backbone))
    return tuple(top)
