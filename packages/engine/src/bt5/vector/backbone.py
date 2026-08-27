"""The user's vector, parsed but not yet assembled with an insert.

The 5'UTR detection here feeds the highest-weight objective in BT5. The window
that explains 44-59% of expression variance in bacteria spans the UTR/CDS
junction, so it cannot be computed from the CDS alone. That makes "no annotated
5'UTR" a real state the app has to represent and degrade on, NOT a case to paper
over by folding the CDS by itself and reporting the number as though the UTR were
there. `UtrContext.five_prime is None` plus a degradation string is what that
honesty looks like in code.

Everything here is strand-aware. For a reverse-oriented cassette the 5' side of
the CDS is at HIGHER construct coordinates, and a detector that assumes otherwise
finds the 3'UTR and calls it the 5'UTR -- silently, on exactly the lentiviral
layouts this tool exists to serve.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from typing import Literal

from bt5.core.types import Feature, Interval, SegmentKind, Strand, Topology
from bt5.vector.markers import is_marker
from bt5.vector.notes import DesignNote

#: GenBank feature keys that mark an immutable, scan-exempt region. Introns are
#: exempt because a deliberately placed chimeric intron is one of the most
#: reliable mammalian expression levers and a naive splice-site sweep would
#: destroy it. LTRs and ITRs are exempt because they are long perfect repeats by
#: construction -- the answer is a strain and temperature protocol, not a redesign.
DEFAULT_EXEMPT_KINDS: Mapping[str, SegmentKind] = {
    "intron": SegmentKind.ANNOTATED_INTRON,
    "ltr": SegmentKind.WHITELISTED_REPEAT,
    "repeat_region": SegmentKind.WHITELISTED_REPEAT,
}

#: Label substrings that mark a whitelisted repeat even when the feature key does
#: not. AAV ITRs are routinely annotated as `misc_feature`.
DEFAULT_EXEMPT_LABELS: tuple[str, ...] = ("itr", "ltr", "inverted terminal repeat")

UtrSource = Literal["annotated_feature", "derived_from_promoter", "absent"]
SiteSource = Literal["annotated_cds", "explicit"]


class VectorError(ValueError):
    """The vector cannot be used as given."""


@dataclass(frozen=True, slots=True)
class InsertionSite:
    """The span the designed CDS replaces, in backbone coordinates."""

    interval: Interval
    label: str
    source: SiteSource
    #: From `/transl_table`, or None when the file did not say. NEVER defaulted
    #: to 1 here: NCBI table 12 reassigns CTG to Ser and table 4 makes TGA a Trp,
    #: so a guessed table is a silently wrong protein. The caller must supply one
    #: and `assemble()` cross-validates it against this.
    detected_table_id: int | None = None

    @property
    def strand(self) -> Strand:
        return self.interval.strand


@dataclass(frozen=True, slots=True)
class UtrContext:
    """Untranslated context around the insertion site.

    Both fields are honestly optional. A missing 5'UTR disables an objective; it
    does not get substituted with a guess.
    """

    five_prime: Interval | None = None
    three_prime: Interval | None = None
    five_prime_source: UtrSource = "absent"
    notes: tuple[DesignNote, ...] = ()

    @property
    def has_five_prime(self) -> bool:
        return self.five_prime is not None


def rotate_interval(iv: Interval, *, by: int, length: int) -> Interval:
    """Rotate an interval into a coordinate frame whose origin moved by `by`.

    Wrapping is preserved automatically: the start is taken modulo the length and
    the original span is re-added, so an interval that wrapped may stop wrapping
    and vice versa, which is exactly the intent.
    """
    start = (iv.start - by) % length
    return Interval(start, start + iv.length, iv.strand)


@dataclass(frozen=True)
class VectorBackbone:
    """A parsed vector: sequence, topology, features, and nothing designed yet."""

    sequence: str
    topology: Topology
    features: tuple[Feature, ...] = ()
    annotations: Mapping[str, str] = field(default_factory=dict)
    name: str = "vector"
    #: uid -> parts, for features whose GenBank location is genuinely
    #: discontiguous (not an origin wrap). Kept so export rebuilds them exactly
    #: instead of flattening a two-exon feature into one span.
    compound_parts: Mapping[str, tuple[Interval, ...]] = field(default_factory=dict)
    notes: tuple[DesignNote, ...] = ()

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def is_circular(self) -> bool:
        return self.topology is Topology.CIRCULAR

    def slice(self, iv: Interval) -> str:
        """Wrap- and strand-aware extraction, matching Construct.slice."""
        from bt5.core.types import reverse_complement

        n = self.length
        if iv.end <= n:
            sub = self.sequence[iv.start : iv.end]
        elif self.is_circular:
            sub = self.sequence[iv.start :] + self.sequence[: iv.end - n]
        else:
            raise VectorError(f"interval {iv} runs past the end of a linear vector of length {n}")
        return reverse_complement(sub) if iv.strand == -1 else sub

    def features_of(self, *kinds: str) -> tuple[Feature, ...]:
        wanted = {k.lower() for k in kinds}
        return tuple(f for f in self.features if f.kind.lower() in wanted)

    def label_of(self, feature: Feature) -> str:
        for key in ("label", "gene", "product", "note"):
            values = feature.qualifiers.get(key)
            if values:
                return values[0]
        return feature.kind

    # -- insertion site ----------------------------------------------------

    def find_insertion_site(self, *, label: str | None = None) -> InsertionSite:
        """Locate the CDS the design replaces.

        Ambiguity is an error, never a silent pick of the first match: choosing
        the wrong CDS in a two-cassette vector produces a plausible plasmid that
        expresses the wrong gene.
        """
        candidates = self.features_of("CDS")
        if label is not None:
            wanted = label.lower()
            candidates = tuple(f for f in candidates if wanted in self.label_of(f).lower())

        if not candidates:
            hint = f" matching {label!r}" if label else ""
            raise VectorError(
                f"no CDS feature{hint} in {self.name!r}; mark the insertion point "
                f"explicitly with insertion_site_from_interval()"
            )
        if len(candidates) > 1:
            labels = ", ".join(repr(self.label_of(f)) for f in candidates)
            raise VectorError(
                f"{len(candidates)} CDS features in {self.name!r} ({labels}); "
                f"pass label= to choose one"
            )

        cds = candidates[0]
        return InsertionSite(
            interval=cds.interval,
            label=self.label_of(cds),
            source="annotated_cds",
            detected_table_id=_transl_table(cds),
        )

    # -- untranslated context ---------------------------------------------

    def utr_context(self, site: InsertionSite, *, max_derived_utr: int = 500) -> UtrContext:
        """Find the 5' and 3' untranslated regions flanking the insertion site.

        Order of preference for the 5'UTR: an annotated `5'UTR` feature, then a
        span derived from the nearest upstream promoter, then nothing. The source
        travels with the answer so the report can say which it was.
        """
        notes: list[DesignNote] = []
        five, source = self._five_prime_utr(site, max_derived_utr, notes)
        three = self._three_prime_utr(site)

        if five is None:
            notes.append(
                DesignNote(
                    kind="unavailable",
                    summary=(
                        "no annotated 5'UTR and no upstream promoter, so the 5' "
                        "folding objective cannot be evaluated for this vector"
                    ),
                    bears_on="protein expression",
                    action=("annotate the 5'UTR or the promoter in your map and re-run"),
                )
            )
        elif source == "derived_from_promoter":
            notes.append(
                DesignNote(
                    kind="assumption",
                    summary=(
                        "5'UTR inferred from the upstream promoter, not annotated; "
                        "the transcription start is assumed, not measured"
                    ),
                    interval=five,
                    bears_on="protein expression",
                    action="annotate the real 5'UTR if you know it",
                )
            )
            if self._overlaps_intron(five):
                notes.append(
                    DesignNote(
                        kind="assumption",
                        summary=(
                            "the inferred 5'UTR contains an annotated intron, so the "
                            "mature 5'UTR after splicing is shorter than this span"
                        ),
                        interval=five,
                        bears_on="protein expression",
                    )
                )
        if three is None:
            notes.append(
                DesignNote(
                    kind="unavailable",
                    summary="no annotated 3'UTR or polyA signal downstream of the insert",
                    bears_on="protein expression",
                )
            )

        return UtrContext(
            five_prime=five,
            three_prime=three,
            five_prime_source=source,
            notes=tuple(notes),
        )

    def _five_prime_utr(
        self, site: InsertionSite, max_derived: int, notes: list[DesignNote]
    ) -> tuple[Interval | None, UtrSource]:
        annotated = self._nearest_upstream(site, ("5'UTR", "five_prime_UTR"))
        if annotated is not None:
            gap = self._gap_to_site(annotated.interval, site)
            if gap > 0:
                notes.append(
                    DesignNote(
                        kind="assumption",
                        summary=(
                            f"the annotated 5'UTR stops {gap} nt short of the start "
                            f"codon; that gap is treated as backbone, not as UTR"
                        ),
                        interval=annotated.interval,
                        bears_on="protein expression",
                    )
                )
            return annotated.interval, "annotated_feature"

        promoter = self._nearest_upstream(site, ("promoter",), skip_markers=True)
        if promoter is not None:
            derived = self._span_between(promoter.interval, site)
            if derived is not None and derived.length <= max_derived:
                return derived, "derived_from_promoter"
            if derived is not None:
                notes.append(
                    DesignNote(
                        kind="unavailable",
                        summary=(
                            f"the nearest promoter, {self.label_of(promoter)!r}, sits "
                            f"{derived.length} nt from the start codon, beyond "
                            f"max_derived_utr={max_derived}; it is probably driving a "
                            f"different transcript, so no 5'UTR was inferred"
                        ),
                        bears_on="protein expression",
                    )
                )
        return None, "absent"

    def _three_prime_utr(self, site: InsertionSite) -> Interval | None:
        annotated = self._nearest_downstream(site, ("3'UTR", "three_prime_UTR"))
        if annotated is not None:
            return annotated.interval
        polya = self._nearest_downstream(site, ("polyA_signal", "polyA_site", "regulatory"))
        if polya is None:
            return None
        return self._span_between(polya.interval, site, downstream=True)

    # -- geometry helpers --------------------------------------------------

    def _forward(self, site: InsertionSite) -> bool:
        return site.strand == 1

    def _upstream_edge(self, site: InsertionSite) -> int:
        """The construct coordinate of the 5' end of the insert."""
        return site.interval.start if self._forward(site) else site.interval.end

    def _downstream_edge(self, site: InsertionSite) -> int:
        return site.interval.end if self._forward(site) else site.interval.start

    def _oriented_gap(self, raw: int) -> int | None:
        """Resolve a signed gap into a distance, or None if it points the wrong way.

        On a circular vector the flanking feature may sit across the origin from
        the insert -- a 5'UTR abutting a CDS that straddles position 0 gives a raw
        gap of -length, which is a gap of zero. Linear arithmetic reports that as
        "not upstream" and loses the UTR on exactly the layouts this lane exists
        to handle. Anything more than half way round is treated as pointing the
        other way, so `nearest` still means nearest.
        """
        if not self.is_circular:
            return raw if raw >= 0 else None
        d = raw % self.length
        return d if d <= self.length // 2 else None

    def _distance_upstream(self, iv: Interval, site: InsertionSite) -> int | None:
        """How far `iv` sits upstream of the insert, or None if it is not upstream."""
        edge = self._upstream_edge(site)
        return self._oriented_gap(edge - iv.end if self._forward(site) else iv.start - edge)

    def _distance_downstream(self, iv: Interval, site: InsertionSite) -> int | None:
        edge = self._downstream_edge(site)
        return self._oriented_gap(iv.start - edge if self._forward(site) else edge - iv.end)

    def _gap_to_site(self, iv: Interval, site: InsertionSite) -> int:
        return self._distance_upstream(iv, site) or 0

    def _nearest_upstream(
        self, site: InsertionSite, kinds: Sequence[str], *, skip_markers: bool = False
    ) -> Feature | None:
        return self._nearest(site, kinds, self._distance_upstream, skip_markers=skip_markers)

    def _nearest_downstream(self, site: InsertionSite, kinds: Sequence[str]) -> Feature | None:
        return self._nearest(site, kinds, self._distance_downstream)

    def _nearest(
        self,
        site: InsertionSite,
        kinds: Sequence[str],
        measure: Callable[[Interval, InsertionSite], int | None],
        *,
        skip_markers: bool = False,
    ) -> Feature | None:
        wanted = {k.lower() for k in kinds}
        best: Feature | None = None
        best_d: int | None = None
        for f in self.features:
            if f.kind.lower() not in wanted or f.interval.strand != site.strand:
                continue
            if skip_markers and is_marker(
                f"{self.label_of(f)} {' '.join(f.qualifiers.get('note', ()))}"
            ):
                # An AmpR promoter is not this transcript's promoter.
                continue
            d = measure(f.interval, site)
            if d is None:
                continue
            if best_d is None or d < best_d:
                best, best_d = f, d
        return best

    def _span_between(
        self, other: Interval, site: InsertionSite, *, downstream: bool = False
    ) -> Interval | None:
        """The untranslated span between a flanking feature and the insert."""
        if downstream:
            lo, hi = (
                (self._downstream_edge(site), other.start)
                if self._forward(site)
                else (other.end, self._downstream_edge(site))
            )
        else:
            lo, hi = (
                (other.end, self._upstream_edge(site))
                if self._forward(site)
                else (self._upstream_edge(site), other.start)
            )
        if hi <= lo:
            if not self.is_circular:
                return None
            # The span runs across the origin; express it in the one wrapping
            # representation BT5 uses everywhere.
            hi += self.length
            if hi <= lo:
                return None
        return Interval(lo, hi, site.strand)

    def _overlaps_intron(self, iv: Interval) -> bool:
        return any(f.interval.overlaps(iv) for f in self.features_of("intron"))

    # -- transforms --------------------------------------------------------

    def rotated(self, by: int) -> VectorBackbone:
        """Move the origin `by` bases downstream. Circular vectors only.

        Used when the insertion site straddles the origin: rotating so the CDS is
        contiguous turns coordinate remapping back into simple arithmetic, at the
        cost of exporting a map whose origin sits somewhere new. That trade is
        recorded as a degradation rather than made silently.
        """
        if not self.is_circular:
            raise VectorError("only a circular vector can be rotated")
        n = self.length
        by %= n
        if by == 0:
            return self
        rotated_features = tuple(
            replace(f, interval=rotate_interval(f.interval, by=by, length=n)) for f in self.features
        )
        rotated_parts = {
            uid: tuple(rotate_interval(p, by=by, length=n) for p in parts)
            for uid, parts in self.compound_parts.items()
        }
        return replace(
            self,
            sequence=self.sequence[by:] + self.sequence[:by],
            features=rotated_features,
            compound_parts=rotated_parts,
            notes=(
                # Existing note intervals are in the OLD frame and must move too,
                # or a located warning ends up pointing at unrelated sequence.
                *(
                    replace(note, interval=rotate_interval(note.interval, by=by, length=n))
                    if note.interval is not None
                    else note
                    for note in self.notes
                ),
                DesignNote(
                    kind="change",
                    summary=(
                        f"the origin was rotated by {by} bp so the insert is contiguous; "
                        f"coordinates in this map differ from your input file"
                    ),
                    bears_on="map fidelity",
                ),
            ),
        )


def insertion_site_from_interval(
    interval: Interval, *, label: str = "insert", table_id: int | None = None
) -> InsertionSite:
    """Mark the insertion point by hand, for FASTA input or an unannotated map."""
    return InsertionSite(
        interval=interval, label=label, source="explicit", detected_table_id=table_id
    )


def _transl_table(feature: Feature) -> int | None:
    values = feature.qualifiers.get("transl_table")
    if not values:
        return None
    try:
        return int(values[0])
    except ValueError:
        return None
