"""D2 -- recombinase and recombination sites: loxP, FRT, Gateway att.

`brief.md:98` grades this **H, check-only**, and both halves of that matter.

**Why HARD_CHECK and not HARD_LATTICE.** These sites are 25-48 bp and, per
`brief.md:90`'s occurrence budget, a 25-mer occurs ~3e-9 times in 1.5 kb -- they never
arise by chance. A loxP in a construct is therefore something the user's backbone
*has*, deliberately, not something the optimiser accidentally wrote. Nothing about
codon choice removes a 34 bp site sitting in a vector, so chasing it would exhaust the
mutation space over bases the solver may not touch. HARD_CHECK is the enforcement class
for exactly that: *"real, but not fixable by codon choice ... reported and blocking,
never chased by the solver"* (`core/spec.py:38-41`).

**Why this rule scans both strands itself.** `CLAUDE.md` §3.4 says to list forward
motifs in `LatticeTerms.forbidden` and let the solver close the set. That mechanism is
unavailable here twice over. The loxP and FRT patterns carry an 8 nt spacer, and
`docs/decisions/2026-09-01-expand-forbidden-iupac.md` settled that N-spacers are
*"wildcards, not ambiguity, and enumerating them is the wrong mechanism"* -- 4^8 = 65,536
patterns against a `MAX_PATTERN_EXPANSION` of 1,024. And a HARD_CHECK rule never reaches
the automaton at all. `d6_non_b_dna` already scans both strands for the same reason (a
G4 is not a finite motif set either), and this rule copies its wrap-aware idiom rather
than inventing a second one.

**Bxb1 is deliberately not detected.** `brief.md:99` states that Bxb1 attB and attP both
contain `GGTCTC` (BsaI) and says to flag that collision -- but it does not give the attB
or attP sequences. Encoding a sequence from anywhere other than the cited source is the
defect `/verify-provenance` exists to catch, so the collision is declared through
`conflicts_with` and said out loud in the report, and site detection waits for a sourced
sequence. See `docs/decisions/2026-09-02-d2-recombinase-check-only.md`.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from typing import ClassVar

from bt5.core.context import ContextSlot, DesignContext
from bt5.core.registry import register
from bt5.core.services import Services
from bt5.core.spec import (
    Breach,
    Citation,
    Direction,
    Enforcement,
    Evaluation,
    Evidence,
    LatticeTerms,
    LocalizationPolicy,
    RepairPolicy,
)
from bt5.core.types import Construct, Interval, Strand, reverse_complement

#: `brief.md:98`'s own upper bound on a recombination site, and the overlap this scan
#: appends when linearising a circular construct so an origin-spanning site is found.
MAX_SITE = 48

#: `brief.md:99`, verbatim. The loxP inverted repeats are exact reverse complements of
#: each other; FRT's are not, which is why both arms are written out rather than one
#: being derived from the other.
SITE_PATTERNS: Mapping[str, tuple[str, str]] = {
    "loxP_family": (
        r"ATAACTTCGTATA[ACGT]{8}TATACGAAGTTAT",
        "loxP family (covers loxP, lox2272, lox5171, loxN, lox511)",
    ),
    "FRT": (r"GAAGTTCCTATTC[ACGT]{8}GTATAGGAACTTC", "FRT"),
    "Gateway_attB1": (r"ACAAGTTTGTACAAAAAAGCAGGCT", "Gateway attB1"),
    "Gateway_attB2": (r"ACCCAGCTTTCTTGTACAAAGTGGT", "Gateway attB2"),
}

#: Reported at a lower magnitude than a full site: a half-site is 13 bp, so
#: `brief.md:90` puts it at ~1e-5 expected occurrences per 1.5 kb -- still essentially
#: never by chance, but a half-site alone cannot recombine.
PARTIAL_PATTERNS: Mapping[str, tuple[str, str]] = {
    "loxP_half_site": (r"ATAACTTCGTATA|TATACGAAGTTAT", "loxP half-site"),
    "Gateway_core": (r"TTTGTACAAA[AG]", "Gateway shared core"),
}

#: `brief.md:99`: "Bxb1 attB and attP both contain GGTCTC (BsaI) -- flag this collision
#: explicitly when a user wants both a landing pad and BsaI-free Golden Gate."
BXB1_COLLISION_ENZYME = "BsaI"
BXB1_COLLISION_MOTIF = "GGTCTC"

MAG_SITE = 1.0
MAG_PARTIAL = 0.3


@register
class RecombinaseSites:
    id: ClassVar[str] = "d2_recombinase_sites"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "Recombinase and recombination sites (loxP, FRT, Gateway att)"
    #: `brief.md:98`: "H, check-only". Reported and blocking, never chased.
    enforcement: ClassVar[Enforcement] = Enforcement.HARD_CHECK
    evidence: ClassVar[Evidence] = Evidence.EVIDENCE_BACKED
    direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
    unit: ClassVar[str] = "sites"
    band: ClassVar[tuple[float, float] | None] = None
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "Expected-occurrence budget: a non-palindromic k-mer occurs ~2(L-k+1)/4^k "
            "times, so a 25-mer is ~3e-9 per 1.5 kb -- a recombination site is present "
            "because someone put it there, never by chance, which is what makes "
            "check-only the right enforcement rather than repair",
            "https://github.com/Edinburgh-Genome-Foundry/DnaChisel",
            2020,
            sign="supports",
        ),
    )
    last_verified: ClassVar[str] = "2026-09-02"
    weight_provenance: ClassVar[str] = ""  # hard rule; not weighted
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.0
    #: HARD_CHECK never reaches the Tier-A DP, so there is nothing to steer.
    steering_weight: ClassVar[float] = 0.0
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.MOTIF_LEN_MINUS_1
    #: Inert, and declared rather than defaulted so the choice is on the record: a
    #: HARD_CHECK rule is never repaired, so no repair can create new instances of what
    #: it removes and FIXED_POINT would describe a loop that does not run.
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "cheap"
    #: `brief.md:99`: Bxb1 attB and attP both carry BsaI's GGTCTC, so a user asking for
    #: a landing pad and BsaI-free Golden Gate has asked for two incompatible things.
    conflicts_with: ClassVar[tuple[str, ...]] = ("d1_restriction_sites",)
    brief_ref: ClassVar[str] = "2.D2"
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "report_partials": {
                "type": "boolean",
                "default": True,
                "description": (
                    "Also report lone half-sites and the Gateway shared core. A half-site "
                    "cannot recombine on its own but is still ~1e-5 per 1.5 kb by chance, "
                    "so its presence is informative."
                ),
            }
        },
    }

    def __init__(self, report_partials: bool = True) -> None:
        self.report_partials = report_partials
        self._compiled = {
            name: (re.compile(pattern), label)
            for name, (pattern, label) in {**SITE_PATTERNS, **PARTIAL_PATTERNS}.items()
        }

    def gate(self, slot: ContextSlot) -> bool:
        """Every context. A loxP is a loxP in any host: `brief.md:98` gives no gating,
        and an unnoticed recombination site is a hazard wherever the DNA ends up."""
        return True

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """None. See the module docstring: the 8 nt spacer would expand to 65,536
        patterns against a cap of 1,024, and a HARD_CHECK rule never reaches the
        automaton that would consume them anyway."""
        return None

    def _hits(self, c: Construct) -> list[tuple[str, str, Interval]]:
        """(name, label, interval) on both strands, wrap-aware.

        The both-strand scan is this rule's own because `forbidden` cannot carry these
        patterns; `d6_non_b_dna._hits` is the same idiom for the same reason.
        """
        n = c.length
        found: list[tuple[str, str, Interval]] = []
        strands: tuple[Strand, ...] = (1, -1)
        for strand in strands:
            seq = c.sequence if strand == 1 else reverse_complement(c.sequence)
            scan = seq + seq[:MAX_SITE] if c.is_circular else seq
            for name, (pattern, label) in self._compiled.items():
                if not self.report_partials and name in PARTIAL_PATTERNS:
                    continue
                for match in pattern.finditer(scan):
                    if match.start() >= n:
                        continue  # the wrap overlap, already reported at its real position
                    # `% n`, never `max(0, ...)`. A minus-strand site crossing the
                    # origin maps to a negative forward coordinate, and clamping it to
                    # zero would name an interval whose bases are not the site found.
                    lo = (match.start() if strand == 1 else n - match.end()) % n
                    found.append((name, label, Interval(lo, lo + len(match.group()), strand)))
        return found

    @staticmethod
    def _within(inner: Interval, outer: Interval, n: int) -> bool:
        """Containment by residue, so a wrapping site swallows its own arms.

        A plain `outer.start <= inner.start and inner.end <= outer.end` is wrong the
        moment the outer site crosses the origin: a loxP stored as [25, 59) on a 46 nt
        plasmid covers residues 25-45 and 0-12, and its first arm is stored as [0, 13),
        which fails that comparison and gets reported as a lone half-site of the very
        site it belongs to. Both spans are bounded by MAX_SITE, so the residue sets are
        small enough to build outright.
        """
        return {i % n for i in range(inner.start, inner.end)} <= {
            i % n for i in range(outer.start, outer.end)
        }

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        breaches: list[Breach] = []
        # Deduplicate on the SPAN and not on (span, strand). loxP's arms are exact
        # reverse complements of each other, so the whole 34 bp pattern matches its own
        # reverse complement at the same coordinates: keying on strand reports one
        # physical site as two and tells the user they have twice the problem they have.
        seen: set[tuple[str, int, int]] = set()
        full_spans = [iv for name, _, iv in self._hits(c) if name not in PARTIAL_PATTERNS]
        for name, label, iv in self._hits(c):
            key = (name, iv.start, iv.end)
            if key in seen:
                continue
            seen.add(key)
            partial = name in PARTIAL_PATTERNS
            # A lone half-site is the finding worth having; the two arms OF a complete
            # site are just that site, described again in smaller pieces.
            if partial and any(self._within(iv, f, c.length) for f in full_spans):
                continue
            carries_bsai = BXB1_COLLISION_MOTIF in c.slice(iv)
            breaches.append(
                Breach(
                    spec_id=self.id,
                    interval=iv,
                    magnitude=MAG_PARTIAL if partial else MAG_SITE,
                    message=(
                        f"{label} at {iv.start} on the "
                        f"{'+' if iv.strand == 1 else '-'} strand"
                        + (
                            "; a half-site cannot recombine alone, but it does not arise "
                            "by chance either"
                            if partial
                            else ""
                        )
                        + (
                            f"; this site also contains {BXB1_COLLISION_MOTIF} "
                            f"({BXB1_COLLISION_ENZYME}), so domesticating for Golden Gate "
                            "and keeping this site are incompatible"
                            if carries_bsai
                            else ""
                        )
                    ),
                    # HARD_CHECK's defining property. A 25-48 bp site is in the
                    # construct because someone put it there; no codon choice removes
                    # it, and telling the solver otherwise sends it after bases it
                    # cannot move until it reports infeasible on a design that is fine.
                    fixable_by_codon_choice=False,
                    detail={
                        "site": name,
                        "strand": float(iv.strand),
                        "partial": "yes" if partial else "no",
                        "bsai_collision": "yes" if carries_bsai else "no",
                    },
                )
            )
        return Evaluation(
            spec_id=self.id,
            passes=not any(b.magnitude >= MAG_SITE for b in breaches),
            raw_score=float(sum(b.magnitude for b in breaches)),
            breaches=tuple(breaches),
            n_evaluated=c.length,
        )
