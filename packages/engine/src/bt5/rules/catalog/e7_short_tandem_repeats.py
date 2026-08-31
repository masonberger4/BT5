"""E7 -- short tandem repeats, the tract every other repeat rule cannot see.

This rule exists because of a measurement, not a taxonomy. Planted in a CDS and
run against the whole catalog as it stood:

    (CAG)x20    = 60 bp poly-Gln   F1=0  F2=0  E5=0  E1=0
    (AT)x30     = 60 bp            F1=0  F2=0  E5=0  E1=0
    (CAGGCT)x15 = 90 bp            F1=0  F2=0  E5=0  E1=0
    A x12       homopolymer        F1=0  F2=0  E5=0  E1=1

A 90 bp tandem array sails through every repeat rule in the catalog. The cause
is structural rather than a tuning miss: the pair scan seeds on k-mers and then
clamps an overlapping pair back to its PERIOD, so a 90 bp array of a 6 bp unit
is reported as a 6 bp match and dropped under every min_len. Pair rules ask
"are these two spans identical?", and a tandem array's answer is "yes, at every
offset", which is the same as no answer.

**Different mechanism, not just a different length.** A dispersed direct repeat
is lost by single-strand annealing between two copies. A tandem array is lost by
POLYMERASE SLIPPAGE -- the nascent strand melts and re-anneals one unit out of
register -- which needs no second site, no loop and no homology search. It
happens during synthesis of the fragment and again during replication in the
host, which is why the same tract appears in both a vendor's complexity screen
and the propagation literature.

**Unit length 2 to 6, and unit 1 is deliberately excluded.** A homopolymer is a
tandem array of period 1, so this scan finds every one of them. Reporting them
would be worse than useless: E1 bands them at 9 nt (A/T) and 5 nt (G/C) and
makes them UNREACHABLE in the Tier-A automaton, so a second finding here at a
20 bp threshold would sit next to E1's in the panel, contradict it, and look
authoritative. Every tract is reduced to its minimal period and dropped if that
period is 1.

**Scope is the ordered fragment**, as for E5 and E6 and for the reason in
`bt5.rules.fragment`. The propagation counterpart -- brief row 2.F6, stricter
still at 5 units for di- and trinucleotides -- is not yet implemented and would
scan the assembled plasmid instead.

The inverted-repeat clause of brief row E7 is F3's and is not duplicated here.
"""

from __future__ import annotations

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
from bt5.rules.fragment import Fragment, fragments
from bt5.rules.vendors import (
    DEFAULT_SELECTION,
    DEFAULT_VENDOR,
    VendorSelection,
    all_keys,
    require_selection,
)

#: Microsatellite unit lengths. 1 is excluded at report time, not here: the scan
#: needs period 1 in order to recognise a homopolymer and hand it to E1.
MIN_UNIT_BP = 2
MAX_UNIT_BP = 6
#: Report a tract at least this long.
WARN_TRACT_BP = 20
#: Above this the tract is a rejected order, not a surcharge.
HARD_TRACT_BP = 100
MAX_FINDINGS = 200


def tandem_tracts(
    seq: str, max_unit: int = MAX_UNIT_BP, min_length: int = 0
) -> list[tuple[int, int, int]]:
    """Maximal tandem tracts as (start, end, period), each at its MINIMAL period.

    Every array is periodic at multiples of its unit -- (CAG)x20 is period 3 and
    also period 6 -- so the same tract is found once per multiple. Reporting each
    would turn one tract into three findings that disagree about the unit. A
    tract is therefore dropped when a smaller-period tract covers it, which also
    reduces every homopolymer to period 1 so the caller can hand it to E1.

    `min_length` is a performance parameter and it matters more than it looks.
    Random sequence matches itself at a quarter of all offsets, so the raw scan
    finds roughly n/4 tiny tracts per period -- about 22,000 on a 15 kb construct
    -- and the containment pass is quadratic in what it keeps. Filtering first
    took the rule from 1,038 ms to single-digit ms at 15 kb.

    The filter is safe rather than approximate: containment requires the covering
    tract to span the covered one, so a container is never shorter than what it
    covers, and dropping everything below `min_length` cannot discard a container
    that a surviving tract still needs.
    """
    n = len(seq)
    found: list[tuple[int, int, int]] = []
    for period in range(1, max_unit + 1):
        i = 0
        while i + period < n:
            if seq[i] != seq[i + period]:
                i += 1
                continue
            j = i
            while j + period < n and seq[j] == seq[j + period]:
                j += 1
            if j + period - i >= min_length:
                found.append((i, j + period, period))
            i = j + 1

    found.sort(key=lambda t: (t[2], t[0]))
    kept: list[tuple[int, int, int]] = []
    for start, end, period in found:
        if any(smaller < period and s <= start and end <= e for s, e, smaller in kept):
            continue
        kept.append((start, end, period))
    return sorted(kept)


@register
class ShortTandemRepeats:
    id: ClassVar[str] = "e7_short_tandem_repeats"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "Short tandem repeat tracts (unit 2-6 bp) in the synthesized fragment"
    enforcement: ClassVar[Enforcement] = Enforcement.HARD_REPAIR
    evidence: ClassVar[Evidence] = Evidence.VENDOR_ASSERTED
    direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
    unit: ClassVar[str] = "weighted tandem tracts"
    band: ClassVar[tuple[float, float] | None] = None
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "Twist gene synthesis FAQ: repeated sequence is a named complexity "
            "trigger independent of length and GC, and short tandem arrays are the "
            "form that fails assembly by slippage rather than by mis-priming",
            "https://www.twistbioscience.com/faq/gene-synthesis",
            2026,
            sign="supports",
        ),
        Citation(
            "Below ~200 bp deletion is RecA-INDEPENDENT -- slipped-strand mispairing "
            "and single-strand annealing -- and is unaffected by recA/recF/recJ/recO, "
            "so a recA- strain does not cover a tandem tract either",
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5426353/",
            2017,
            sign="supports",
        ),
        Citation(
            "Longest repetitive sequence is one of the two highest-importance "
            "features of the only synthesis-success model trained on real vendor "
            "outcomes (random forest, 1,076 orders, F1 0.928)",
            "https://pubs.acs.org/doi/10.1021/acssynbio.9b00460",
            2020,
            sign="supports",
        ),
    )
    last_verified: ClassVar[str] = "2026-08-28"
    weight_provenance: ClassVar[str] = ""  # hard rule; never weighted
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.0
    #: Below f1's 1.0 and e5's 0.7. A tandem tract in a CDS is usually FORCED by
    #: the protein -- poly-Gln, poly-Gly, a (GGGGS) linker -- so the repair is to
    #: alternate synonymous codons and break the period, not to steer away from
    #: the residues. Steering hard here would push the DP against sequence it
    #: cannot legally change.
    steering_weight: ClassVar[float] = 0.5
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.WHOLE_SCOPE
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "cheap"
    conflicts_with: ClassVar[tuple[str, ...]] = ()
    brief_ref: ClassVar[str] = "2.E7"
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "max_unit": {"type": "integer", "default": MAX_UNIT_BP, "minimum": MIN_UNIT_BP},
            "warn_tract": {"type": "integer", "default": WARN_TRACT_BP},
            "hard_tract": {"type": "integer", "default": HARD_TRACT_BP},
            "vendors": {
                "type": "array",
                "minItems": 1,
                "uniqueItems": True,
                "items": {"type": "string", "enum": list(all_keys())},
                "default": [DEFAULT_VENDOR],
                "description": (
                    "Which vendor product the fragment is ordered as. Only the "
                    "adapter-on options carry adapters; a plain Gene Fragment "
                    "order is the ordered DNA and nothing else."
                ),
            },
        },
    }

    def __init__(
        self,
        max_unit: int = MAX_UNIT_BP,
        warn_tract: int = WARN_TRACT_BP,
        hard_tract: int = HARD_TRACT_BP,
        vendors: VendorSelection = DEFAULT_SELECTION,
    ) -> None:
        if max_unit < MIN_UNIT_BP:
            raise ValueError(
                f"max_unit {max_unit} leaves nothing to scan: unit 1 is a homopolymer "
                f"and belongs to e1_homopolymers, which bands it far more strictly"
            )
        if warn_tract < 2 * max_unit:
            raise ValueError(
                f"warn_tract {warn_tract} is under two copies of a {max_unit} bp unit, "
                f"so it would report sequence that is not a tandem array at all"
            )
        if hard_tract < warn_tract:
            raise ValueError(f"hard_tract {hard_tract} must not be below {warn_tract}")
        self.max_unit = max_unit
        self.warn_tract = warn_tract
        self.hard_tract = hard_tract
        self.vendors = require_selection(vendors)

    def gate(self, slot: ContextSlot) -> bool:
        return True

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None:
        """None, and this one is worth stating.

        A tandem tract IS decidable from a bounded suffix -- `hard_tract` bases
        of it -- so it looks like a lattice rule. It is not one: the automaton
        state would have to carry those bases to know how long the current tract
        already is, which is the same combinatorial explosion that keeps windowed
        GC out of Tier A. Steer, repair, and let the validator refuse.
        """
        return None

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        adapters = self.vendors.adapters
        breaches: list[Breach] = []
        worst = 0
        scanned = 0

        for frag in fragments(c, adapters):
            scanned += frag.ordered.length
            for start, end, period in tandem_tracts(frag.sequence, self.max_unit, self.warn_tract):
                if len(breaches) >= MAX_FINDINGS:
                    break
                length = end - start
                if period < MIN_UNIT_BP or length < self.warn_tract:
                    continue
                mapped = frag.to_construct(Interval(start, end))
                if mapped is None:
                    continue
                worst = max(worst, length)
                breaches.append(self._breach(frag, mapped, start, length, period))

        return Evaluation(
            spec_id=self.id,
            passes=worst <= self.hard_tract,
            raw_score=float(sum(b.magnitude for b in breaches)),
            breaches=tuple(breaches),
            n_evaluated=scanned,
        )

    def _breach(
        self, frag: Fragment, mapped: Interval, start: int, length: int, period: int
    ) -> Breach:
        unit = frag.sequence[start : start + period]
        copies = length / period
        hard = length > self.hard_tract
        return Breach(
            spec_id=self.id,
            interval=mapped,
            # Slippage rises with copy number, so the overrun is the magnitude
            # rather than a flat band, doubled once the tract is unbuildable.
            magnitude=(length - self.warn_tract + 1) / self.warn_tract * (2.0 if hard else 1.0),
            message=(
                f"{length} bp tandem repeat at {mapped.start}: {copies:.1f} copies of "
                f"{unit!r}. "
                + (
                    f"Over the {self.hard_tract} bp tract limit -- polymerase slippage "
                    f"during assembly makes this unbuildable as ordered"
                    if hard
                    else "Polymerase slippage deletes tracts like this during synthesis "
                    "and again during propagation"
                )
                + ". No pair-based repeat rule can see a tandem array, and a recA- "
                "strain does not cover it either. Where the protein forces the tract "
                "(poly-Gln, poly-Gly, a GGGGS linker), alternate synonymous codons to "
                "break the period rather than changing the residues."
            ),
            fixable_by_codon_choice=frag.construct.overlaps_editable(
                Interval(start, start + length)
            ),
            detail={
                "vendor": self.vendors.label,
                "tract_bp": float(length),
                "unit_bp": float(period),
                "unit": unit,
                "copies": round(copies, 1),
                "severity": "hard" if hard else "warn",
            },
        )
