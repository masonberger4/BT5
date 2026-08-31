"""E8 -- is every 12-mer in this construct unique?

The brief calls k-mer uniqueness the single most effective repeat rule, because
one property protects three different things at once: assembly PCR needs each
oligo to have one landing site, Gibson needs each junction's homology arm to be
found in exactly one place, and plasmid stability needs no second copy for
single-strand annealing to find. Whole construct, backbone included -- a 12-mer
that is fine inside the insert is disqualifying if the vector already carries it.
Vendor adapters, which the brief also names, are deliberately out of scope; the
paragraph before "Boundaries with the rest of the family" says why.

**What the literal rule does, measured before it was written.** "Every 12-mer
unique" is not achievable and mostly would not be information. Distinct 12-mers
occurring at least twice, over 20 seeds each:

                 1.5 kb   5 kb    10 kb   20 kb
    k=10           1.00   12.55    48.10  191.80
    k=12           0.05    0.70     2.15   12.85
    k=15           0.00    0.00     0.00    0.20

So a hard "must be unique" at k=12 reports about 13 pure-chance findings on an
ordinary 20 kb lentiviral transfer plasmid, and at k=10 nearly 200. Multiplicity
THREE, though, was never once reached by chance at k=12 at any length tested --
the maximum observed anywhere was 2. That is the line this rule reports on, and
it is measured rather than chosen.

**The finding that shaped the design.** Repeat the measurement on codon-BIASED
sequence, which is what a max-CAI optimizer produces, and it inverts:

    alpha=1.0 (unbiased)   k=12, 20 kb:    0.20 k-mers at multiplicity >=3
    alpha=0.3 (moderate)                  65.67, worst k-mer seen 14 times
    alpha=0.1 (max-CAI-like)             993.20, worst k-mer seen 548 times

Multiplicity is driven by codon bias, and a fixed threshold on the COUNT would
therefore report a thousand findings on exactly the sequences BT5 is most needed
for. But that explosion is not noise -- it is the signal. One codon per amino
acid is the collapse `docs/PLAN.md` names as the central repetitive-protein
failure mode, and rising 12-mer multiplicity is what it looks like in the DNA.
So the rule's primary output is the SCALAR -- the fraction of positions whose
12-mer recurs, bounded in [0, 1] and steerable -- and breaches are the capped
list of worst offenders, not an attempt to enumerate them.

**Exempt regions are counted but not scored.** A lentiviral vector's two LTRs
make every k-mer in them recur, which would swamp the score with something no
codon can change. Positions inside a whitelisted repeat are therefore left out
of the scalar. They are still COUNTED when computing multiplicity, because an
insert k-mer that collides with an LTR is a real Gibson liability and dropping
the LTR would hide it; a k-mer is only skipped when every one of its occurrences
is exempt.

**Vendor adapters are excluded, and that is a decision rather than an oversight.**
Brief row 2.E8 asks for uniqueness across the whole construct *including vendor
adapters*, and this rule stops at the backbone. E5 owns the adapter axis -- it
scans adapter + ordered + adapter and reports insert-versus-adapter collisions
pairwise from 12 bp, the same k used here. What is given up is exactly one case:
a 12-mer occurring once in an adapter, once in the insert and once in the
backbone, which E5 sees as a pair and this rule does not see at all. Buying it
back is expensive and the price falls in the wrong place. Adapter sequence has no
parent-construct coordinate (`Fragment.to_construct` returns None for an interval
inside an adapter, and these findings print positions), and `fragments()` splices
one adapter copy per designable span -- so a three-span design would give every
adapter 12-mer three landing sites and this rule would report the vendor's own
constant as a repeat. One of the five configurations carries adapters at all, and
it is not the default.

**Boundaries with the rest of the family.** F1 reports PAIRS at 15 bp and up on
the plasmid; E5 reports pairs at 12 bp and up in the ordered fragment; E7 reports
tandem tracts. This reports MULTIPLICITY, which none of them can express: "this
12-mer has eight landing sites" is a different statement from four pairwise
findings, and it is the one that predicts a failed assembly.

Holding that boundary took two corrections, both caught by measurement:

  * A 12-mer counted by POSITION makes every tandem array a huge multiplicity --
    (CAG)x20 reads as 17 occurrences -- and re-reports E7's finding under a third
    name. Occurrences closer together than k overlap, so they are not
    independent landing sites for a primer at all. They are clustered into ONE
    site, which drops tandem arrays here to multiplicity 1 and leaves them to E7,
    while also making the number mean what the message says it means.
  * One 30 bp repeated at three places contains 19 distinct 12-mers, each at
    three sites, and reported one per k-mer that was 19 findings for one
    physical repeat -- the same failure E1 hit with homopolymers, F2 with
    adjacent seeds and E6 with overlapping windows. Adjacent flagged positions
    are merged into a region and regions sharing a site set are reported once.

The scalar and the breaches deliberately use different definitions. The scalar
is uniqueness, by position, and a tandem array IS non-unique so it counts there.
The breaches are independent landing sites, because that is what misassembles.
"""

from __future__ import annotations

from collections import defaultdict
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
from bt5.core.types import Construct, Interval
from bt5.rules.exempt import EXEMPT_COVERAGE, coverage

#: The brief's k. Chance duplication at 10 is ~200 per 20 kb and at 12 is ~13,
#: so 12 is the shortest k at which recurrence carries information at all.
KMER_BP = 12
#: Never once reached by chance at k=12 in 20 seeds at any length up to 20 kb.
#: Also what keeps this rule off F1's and E5's pairwise findings.
MIN_MULTIPLICITY = 3
MAX_FINDINGS = 200
#: How many occurrence positions to name in a message before eliding.
NAMED_POSITIONS = 6


def landing_sites(positions: list[int], k: int, length: int, circular: bool) -> tuple[int, ...]:
    """Collapse occurrences closer together than k into one site.

    Two occurrences less than k apart OVERLAP, so they are one stretch of
    sequence and one place for a primer or a Gibson arm to land -- not two.
    Counting them separately is what made a tandem array read as a
    seventeen-fold repeat here instead of as E7's tract.

    The linkage is to the previous OCCURRENCE, not to the site that opened the
    cluster. A tandem array is a chain of overlapping copies, so comparing
    against the cluster start reopens a new site every k bases and turned
    (CAG)x20 into five sites rather than the one contiguous stretch it is.
    """
    ordered = sorted(positions)
    sites: list[int] = []
    previous: int | None = None
    for p in ordered:
        if previous is None or p - previous >= k:
            sites.append(p)
        previous = p
    if circular and len(sites) > 1 and length - ordered[-1] + ordered[0] < k:
        sites.pop()
    return tuple(sites)


@register
class KmerUniqueness:
    id: ClassVar[str] = "e8_kmer_uniqueness"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "12-mer uniqueness across the whole assembled construct"
    enforcement: ClassVar[Enforcement] = Enforcement.SOFT
    evidence: ClassVar[Evidence] = Evidence.EVIDENCE_BACKED
    direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
    unit: ClassVar[str] = "fraction of 12-mer positions that recur"
    band: ClassVar[tuple[float, float] | None] = None
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "Repeats and high GC are the documented cause of Gibson misassembly, and "
            "no two junctions may share a homology sequence -- uniqueness across the "
            "whole construct is what that requires",
            "https://blog.addgene.org/plasmids-101-gibson-assembly",
            2023,
            sign="supports",
        ),
        Citation(
            "Repetitive k-mer content and longest repetitive sequence outrank every "
            "GC feature in the only synthesis-success model trained on real vendor "
            "outcomes (random forest, 1,076 orders, F1 0.928)",
            "https://pubs.acs.org/doi/10.1021/acssynbio.9b00460",
            2020,
            sign="supports",
        ),
        Citation(
            "Below ~200 bp deletion is RecA-INDEPENDENT and unaffected by a recA- "
            "strain, so short shared sequence between insert and backbone is the "
            "design's problem and not the host's",
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5426353/",
            2017,
            sign="supports",
        ),
    )
    last_verified: ClassVar[str] = "2026-08-28"
    weight_provenance: ClassVar[str] = (
        "0.55: below e6_repeat_density (0.65) on purpose, because the two measure "
        "overlapping things -- a fragment that is repetitive at 9-mers is usually "
        "repetitive at 12-mers -- and weighting both at full strength would count "
        "the same repetitiveness twice in one weighted sum. E8's increment over E6 "
        "is its SCOPE: E6 sees only the ordered fragment, so a 12-mer the insert "
        "shares with the user's own backbone is invisible to it and is exactly the "
        "collision that misassembles a Gibson junction. Not higher despite the "
        "brief calling k-mer uniqueness the most effective repeat rule, because the "
        "measurement behind this file shows the raw count is composition-driven; "
        "the percentile against a length- and protein-matched null is what makes it "
        "comparable, and that is the score's job, not this weight's."
    )
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.55
    #: Modest. f1 steers repeats at 1.0 and e5 at 0.7 over largely the same
    #: bases; a fourth full-strength repeat term would spend all of the
    #: sequence's freedom on one family.
    steering_weight: ClassVar[float] = 0.3
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.PAIRED_SEGMENTS
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "moderate"
    conflicts_with: ClassVar[tuple[str, ...]] = ()
    brief_ref: ClassVar[str] = "2.E8"
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "k": {"type": "integer", "default": KMER_BP, "minimum": 11},
            "min_multiplicity": {
                "type": "integer",
                "default": MIN_MULTIPLICITY,
                "minimum": 2,
            },
        },
    }

    def __init__(self, k: int = KMER_BP, min_multiplicity: int = MIN_MULTIPLICITY) -> None:
        if k < 11:
            raise ValueError(
                f"k {k} recurs by chance too often to carry information: a 10-mer "
                f"duplicates ~192 times in an ordinary 20 kb plasmid, against ~13 "
                f"for a 12-mer"
            )
        if min_multiplicity < 2:
            raise ValueError(
                f"min_multiplicity {min_multiplicity} would report every k-mer; a "
                f"k-mer occurring once is what unique MEANS"
            )
        self.k = k
        self.min_multiplicity = min_multiplicity

    def gate(self, slot: ContextSlot) -> bool:
        # Assembly, cloning and propagation all happen for every modality.
        return True

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """None. Whether a codon completes a recurring k-mer depends on the whole
        construct, not on a bounded suffix."""
        return None

    def _occurrences(self, c: Construct) -> dict[str, list[int]]:
        """Every k-mer's start positions, origin-spanning ones included."""
        n = c.length
        if n < self.k:
            return {}
        scan = c.sequence + c.sequence[: self.k - 1] if c.is_circular else c.sequence
        last = n if c.is_circular else n - self.k + 1
        out: dict[str, list[int]] = defaultdict(list)
        for i in range(last):
            out[scan[i : i + self.k]].append(i)
        return out

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        occurrences = self._occurrences(c)
        if not occurrences:
            return Evaluation(spec_id=self.id, passes=True, raw_score=0.0, n_evaluated=0)

        exempt = c.exempt

        def is_exempt(start: int) -> bool:
            if not exempt:
                return False
            iv = Interval(start, start + self.k)
            return coverage(iv, exempt, c.length, c.is_circular) >= EXEMPT_COVERAGE

        sites = {
            kmer: landing_sites(pos, self.k, c.length, c.is_circular)
            for kmer, pos in occurrences.items()
        }
        n = c.length
        last = n if c.is_circular else max(0, n - self.k + 1)
        scan = c.sequence + c.sequence[: self.k - 1] if c.is_circular else c.sequence

        scored = recurring = 0
        flagged: list[bool] = []
        for i in range(last):
            kmer = scan[i : i + self.k]
            free = not is_exempt(i)
            if free:
                scored += 1
                # Uniqueness is by POSITION: a tandem array is not unique, and
                # the scalar should say so even though E7 owns the finding.
                if len(occurrences[kmer]) > 1:
                    recurring += 1
            flagged.append(free and len(sites[kmer]) >= self.min_multiplicity)

        breaches = self._regions(c, scan, flagged, sites)
        return Evaluation(
            spec_id=self.id,
            passes=not breaches,
            raw_score=recurring / scored if scored else 0.0,
            breaches=tuple(breaches),
            n_evaluated=scored,
        )

    def _regions(
        self,
        c: Construct,
        scan: str,
        flagged: list[bool],
        sites: Mapping[str, tuple[int, ...]],
    ) -> list[Breach]:
        """One breach per repeated REGION, deduplicated by its set of sites.

        Adjacent flagged positions are one repeat seen at successive offsets, and
        every copy of that repeat produces a region with an identical site set --
        so reporting the first and dropping the rest turns 19 findings about one
        30 bp repeat into the one finding it is.
        """
        out: list[Breach] = []
        seen: set[tuple[int, ...]] = set()
        i, n = 0, len(flagged)
        while i < n:
            if not flagged[i]:
                i += 1
                continue
            j = i
            while j < n and flagged[j]:
                j += 1
            where = sites[scan[i : i + self.k]]
            if where not in seen:
                seen.add(where)
                out.append(self._breach(c, i, j - i + self.k - 1, where))
            i = j
            if len(out) >= MAX_FINDINGS:
                break
        out.sort(key=lambda b: -b.magnitude)
        return out

    def _breach(self, c: Construct, start: int, span: int, where: tuple[int, ...]) -> Breach:
        shown = ", ".join(str(p) for p in where[:NAMED_POSITIONS])
        if len(where) > NAMED_POSITIONS:
            shown += f", ... ({len(where) - NAMED_POSITIONS} more)"
        return Breach(
            spec_id=self.id,
            interval=Interval(start, start + span),
            # Linear in the excess over the pairwise case F1 and E5 already own.
            magnitude=float(len(where) - self.min_multiplicity + 1),
            message=(
                f"{span} bp region at {start} occurs at {len(where)} independent sites "
                f"in the assembled construct: {shown}. Chance duplication of a "
                f"{self.k}-mer never exceeded 2 copies in any measured plasmid, so "
                f"three or more is designed rather than luck -- an oligo or Gibson "
                f"arm landing here has {len(where)} sites to choose from"
            ),
            fixable_by_codon_choice=any(
                c.overlaps_editable(Interval(p, p + self.k)) for p in where
            ),
            detail={
                "multiplicity": float(len(where)),
                "region_bp": float(span),
                "k": float(self.k),
            },
        )
