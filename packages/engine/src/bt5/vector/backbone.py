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

It is also host-shaped in two directions rather than one. Cap-dependent scanning
makes the transcription start the 5' anchor in a eukaryote, so a promoter is the
right thing to derive a leader from. Bacteria do not scan: the ribosome is
recruited by the Shine-Dalgarno sequence, which is annotated as `RBS` (SnapGene)
or as `regulatory` with `/regulatory_class="ribosome_binding_site"` (INSDC). A
detector that knows only about promoters reports "the transcription start is
assumed" on a pET vector whose initiation element is annotated 30 nt from the
start codon -- true, but it buries the fact that the element that actually
matters is right there and known.
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

#: INSDC deprecated the standalone feature keys below in favour of `regulatory`
#: plus a `/regulatory_class`, but the deprecated keys are what SnapGene writes
#: and what every Addgene deposit therefore carries. Both spellings are matched,
#: and `regulatory` is matched BY CLASS: an untyped `regulatory` match would let
#: a downstream terminator be read as a polyA signal.
RBS_KINDS: tuple[str, ...] = ("RBS", "ribosome_binding_site")
RBS_CLASSES: tuple[str, ...] = ("ribosome_binding_site", "shine_dalgarno_sequence")
POLYA_KINDS: tuple[str, ...] = ("polyA_signal", "polyA_site")
POLYA_CLASSES: tuple[str, ...] = ("polya_signal_sequence", "polya_site")
PROMOTER_KINDS: tuple[str, ...] = ("promoter",)
PROMOTER_CLASSES: tuple[str, ...] = ("promoter",)

#: A Shine-Dalgarno sits roughly 5-13 nt from the start codon; the 16S contact
#: and the P-site cannot be much further apart than that. An annotated RBS well
#: outside that range is not a scoring question, it is evidence that the start
#: codon BT5 picked is not the one the ribosome uses -- which is this module's
#: business, since choosing the insertion site is what it does. Generous, so
#: that it fires on a real disagreement rather than on a spacing preference.
MAX_RBS_TO_START_BP = 30

UtrSource = Literal["annotated_feature", "derived_from_promoter", "derived_from_rbs", "absent"]
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

    Every field is honestly optional. A missing 5'UTR disables an objective; it
    does not get substituted with a guess.

    `ribosome_binding_site` is reported whenever one is annotated upstream, no
    matter which path produced `five_prime`. A bacterial initiation-rate model
    needs the Shine-Dalgarno position itself, not just a leader that happens to
    contain it, and re-deriving it downstream by re-scanning the features would
    mean two places that have to agree about strand and origin-wrap arithmetic.
    """

    five_prime: Interval | None = None
    three_prime: Interval | None = None
    five_prime_source: UtrSource = "absent"
    #: The annotated Shine-Dalgarno, in the same coordinate frame as the rest.
    ribosome_binding_site: Interval | None = None
    #: Bases between the 3' end of that RBS and the start codon, or None when no
    #: RBS is annotated. Reported, never scored -- the band is a rule, not a
    #: property of the vector.
    rbs_spacing: int | None = None
    notes: tuple[DesignNote, ...] = ()

    @property
    def has_five_prime(self) -> bool:
        return self.five_prime is not None

    @property
    def has_ribosome_binding_site(self) -> bool:
        return self.ribosome_binding_site is not None


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
        rbs = self._nearest_upstream(site, RBS_KINDS, classes=RBS_CLASSES)
        rbs_iv = rbs.interval if rbs is not None else None
        spacing = self._distance_upstream(rbs_iv, site) if rbs_iv is not None else None
        five, source = self._five_prime_utr(site, max_derived_utr, notes, rbs_iv)
        three = self._three_prime_utr(site)

        notes.extend(self._rbs_notes(rbs_iv, spacing, five, source))

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
            covers_rbs = (
                rbs_iv is not None
                and five is not None
                and self._contains(five, rbs_iv)
                # A Shine-Dalgarno too far away to reach this start codon is
                # not evidence about this start codon, however neatly it falls
                # inside the derived span.
                and spacing is not None
                and spacing <= MAX_RBS_TO_START_BP
            )
            notes.append(
                DesignNote(
                    kind="assumption",
                    summary=(
                        "5'UTR inferred from the upstream promoter, not annotated; "
                        "the transcription start is assumed, not measured"
                        + (
                            "; the annotated ribosome binding site lies inside "
                            "this span, so the initiation element is known even "
                            "though the leader's 5' end is not"
                            if covers_rbs
                            else ""
                        )
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
            ribosome_binding_site=rbs_iv,
            rbs_spacing=spacing,
            notes=tuple(notes),
        )

    def _five_prime_utr(
        self,
        site: InsertionSite,
        max_derived: int,
        notes: list[DesignNote],
        rbs: Interval | None,
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

        promoter = self._nearest_upstream(
            site, PROMOTER_KINDS, classes=PROMOTER_CLASSES, skip_markers=True
        )
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

        # No promoter, or one too far to be this transcript's. In a bacterium
        # that is not the end of the road: the ribosome is recruited by the
        # Shine-Dalgarno, not by the transcription start, so an annotated RBS
        # anchors a leader that is short but entirely real. Reporting "no 5'UTR"
        # while an annotated initiation element sits 30 nt from the start codon
        # throws away the better part of the one window that carries the
        # published signal.
        if rbs is not None:
            anchored = self._span_between(rbs, site, include_other=True)
            if anchored is not None and anchored.length <= max_derived:
                return anchored, "derived_from_rbs"
        return None, "absent"

    def _three_prime_utr(self, site: InsertionSite) -> Interval | None:
        annotated = self._nearest_downstream(site, ("3'UTR", "three_prime_UTR"))
        if annotated is not None:
            return annotated.interval
        polya = self._nearest_downstream(site, POLYA_KINDS, classes=POLYA_CLASSES)
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
        self,
        site: InsertionSite,
        kinds: Sequence[str],
        *,
        classes: Sequence[str] = (),
        skip_markers: bool = False,
    ) -> Feature | None:
        return self._nearest(
            site, kinds, self._distance_upstream, classes=classes, skip_markers=skip_markers
        )

    def _nearest_downstream(
        self, site: InsertionSite, kinds: Sequence[str], *, classes: Sequence[str] = ()
    ) -> Feature | None:
        return self._nearest(site, kinds, self._distance_downstream, classes=classes)

    def _nearest(
        self,
        site: InsertionSite,
        kinds: Sequence[str],
        measure: Callable[[Interval, InsertionSite], int | None],
        *,
        classes: Sequence[str] = (),
        skip_markers: bool = False,
    ) -> Feature | None:
        wanted = {k.lower() for k in kinds}
        wanted_classes = {c.lower() for c in classes}
        best: Feature | None = None
        best_d: int | None = None
        for f in self.features:
            if not _is_element(f, wanted, wanted_classes) or f.interval.strand != site.strand:
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
        self,
        other: Interval,
        site: InsertionSite,
        *,
        downstream: bool = False,
        include_other: bool = False,
    ) -> Interval | None:
        """The untranslated span between a flanking feature and the insert.

        `include_other` extends the span over the flanking feature itself. A
        promoter is transcribed from its 3' end, so the leader starts where the
        promoter stops; a Shine-Dalgarno is IN the leader and is the part of it
        that does the work, so anchoring on one has to take it in.
        """
        if downstream:
            far = other.end if include_other else other.start
            near = other.start if include_other else other.end
            lo, hi = (
                (self._downstream_edge(site), far)
                if self._forward(site)
                else (near, self._downstream_edge(site))
            )
        else:
            near = other.start if include_other else other.end
            far = other.end if include_other else other.start
            lo, hi = (
                (near, self._upstream_edge(site))
                if self._forward(site)
                else (self._upstream_edge(site), far)
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

    def _contains(self, outer: Interval, inner: Interval) -> bool:
        """Wrap-aware containment; both intervals live in this vector's frame.

        A leader derived across the origin is stored with `end > length`, while
        the RBS inside it is stored plainly, so plain comparison says no. Trying
        the inner interval shifted a full turn is what makes the two agree.
        """
        shifts = (0, self.length) if self.is_circular else (0,)
        return any(
            outer.start <= inner.start + shift and inner.end + shift <= outer.end
            for shift in shifts
        )

    def _rbs_notes(
        self,
        rbs: Interval | None,
        spacing: int | None,
        five: Interval | None,
        source: UtrSource,
    ) -> tuple[DesignNote, ...]:
        """What an annotated Shine-Dalgarno says about the site BT5 picked."""
        if rbs is None:
            return ()
        out: list[DesignNote] = []
        if source == "derived_from_rbs":
            out.append(
                DesignNote(
                    kind="assumption",
                    summary=(
                        "no annotated 5'UTR and no promoter within range, so the 5' "
                        "leader was anchored on the annotated ribosome binding site; "
                        "it covers the initiation element but stops there, not at "
                        "the real transcription start"
                    ),
                    interval=five,
                    bears_on="protein expression",
                    action="annotate the 5'UTR or the promoter for the full leader",
                )
            )
        elif five is not None and not self._contains(five, rbs):
            # The two 5' anchors disagree. Either the leader is not this
            # transcript's or the RBS is not -- and following the wrong one puts
            # the highest-weight objective on the wrong bases.
            out.append(
                DesignNote(
                    kind="liability",
                    summary=(
                        "an annotated ribosome binding site sits outside the 5'UTR "
                        "used here, so the two 5' annotations describe different "
                        "transcripts"
                    ),
                    interval=rbs,
                    bears_on="protein expression",
                    action="check that the insertion site is the ORF this RBS serves",
                )
            )
        if spacing is not None and spacing > MAX_RBS_TO_START_BP:
            out.append(
                DesignNote(
                    kind="liability",
                    summary=(
                        f"the annotated ribosome binding site stops {spacing} nt "
                        f"before the start codon, further than the {MAX_RBS_TO_START_BP} nt "
                        f"a Shine-Dalgarno reaches; the ribosome probably initiates "
                        f"at a different codon than the one BT5 is designing from"
                    ),
                    interval=rbs,
                    bears_on="protein expression",
                    action="check for an upstream in-frame start codon between the two",
                )
            )
        return tuple(out)

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


def _is_element(feature: Feature, kinds: set[str], classes: set[str]) -> bool:
    """Match a feature by its key, or by `/regulatory_class` under `regulatory`.

    `regulatory` is deliberately never matched on the key alone. It is a
    catch-all covering terminators, attenuators, operators and polyA signals
    alike, so an untyped match reads whichever of them happens to be nearest as
    the element being looked for.
    """
    kind = feature.kind.lower()
    if kind in kinds:
        return True
    if kind != "regulatory" or not classes:
        return False
    return any(v.lower() in classes for v in feature.qualifiers.get("regulatory_class", ()))


def _transl_table(feature: Feature) -> int | None:
    values = feature.qualifiers.get("transl_table")
    if not values:
        return None
    try:
        return int(values[0])
    except ValueError:
        return None
