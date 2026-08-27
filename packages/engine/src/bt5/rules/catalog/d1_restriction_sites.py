"""D1 -- forbidden restriction / Type IIS recognition sites.

REFERENCE RULE. Copy this file's shape when adding a pattern rule.

Note what this rule does NOT do:
  - it does not scan the reverse strand itself. It lists forward motifs only and
    the SOLVER closes the pattern set under reverse complement when building the
    automaton. Closing it here too would double-count.
  - it does not iterate or repair. Declaring HARD_LATTICE means the Tier-A DP
    makes the motif unreachable by construction, including across codon
    boundaries, at the CDS/backbone junction and spanning the origin.
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
from bt5.verify import find_motifs

#: Type IIS enzymes are DEFAULT-OFF: Golden Gate is not the assumed assembly
#: method here (Gibson is), so domestication is opt-in rather than imposed.
TYPE_IIS: dict[str, str] = {
    "BsaI": "GGTCTC",
    "BsmBI": "CGTCTC",
    "BbsI": "GAAGAC",
    "SapI": "GCTCTTC",
}

#: Common six-cutters, enabled by default because an accidental site in the
#: insert breaks downstream diagnostic digests.
SIX_CUTTERS: dict[str, str] = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "XhoI": "CTCGAG",
    "XbaI": "TCTAGA",
    "NotI": "GCGGCCGC",
}


@register
class RestrictionSites:
    id: ClassVar[str] = "d1_restriction_sites"
    version: ClassVar[str] = "1.0.0"
    title: ClassVar[str] = "Forbidden restriction / Type IIS sites"
    enforcement: ClassVar[Enforcement] = Enforcement.HARD_LATTICE
    evidence: ClassVar[Evidence] = Evidence.EVIDENCE_BACKED
    direction: ClassVar[Direction] = Direction.LOWER_IS_BETTER
    unit: ClassVar[str] = "sites"
    citations: ClassVar[tuple[Citation, ...]] = (
        Citation(
            "Expected-occurrence budget: a non-palindromic k-mer occurs ~2(L-k+1)/4^k "
            "times on both strands, so a 6-mer appears ~0.73x per 1.5 kb -- cheap to forbid",
            "https://github.com/Edinburgh-Genome-Foundry/DnaChisel",
            2020,
        ),
    )
    last_verified: ClassVar[str] = "2026-08-27"
    weight_provenance: ClassVar[str] = ""  # hard rule; not weighted
    default_enabled: ClassVar[bool] = True
    default_weight: ClassVar[float] = 0.0
    steering_weight: ClassVar[float] = 0.0  # guaranteed by construction; no steering needed
    band: ClassVar[tuple[float, float] | None] = None
    localization: ClassVar[LocalizationPolicy] = LocalizationPolicy.MOTIF_LEN_MINUS_1
    repair: ClassVar[RepairPolicy] = RepairPolicy.SINGLE_PASS
    cost_class: ClassVar[str] = "cheap"
    conflicts_with: ClassVar[tuple[str, ...]] = ("b8_kozak",)  # NcoI CCATGG lies inside GCCACCATGG
    brief_ref: ClassVar[str] = "2.D1"
    engine_calibration: ClassVar[str | None] = None
    param_schema: ClassVar[Mapping[str, object]] = {
        "type": "object",
        "properties": {
            "enzymes": {
                "type": "array",
                "items": {"type": "string"},
                "default": sorted(SIX_CUTTERS),
                "description": "Enzyme names whose recognition sites must be absent.",
            }
        },
    }

    def __init__(self, enzymes: tuple[str, ...] | None = None) -> None:
        table = {**SIX_CUTTERS, **TYPE_IIS}
        names = enzymes if enzymes is not None else tuple(SIX_CUTTERS)
        self.enzymes = {n: table[n] for n in names if n in table}

    def gate(self, slot: ContextSlot) -> bool:
        return True

    def enforcement_for(self, slot: ContextSlot) -> Enforcement:
        return self.enforcement

    def lattice_terms(self, ctx: DesignContext) -> LatticeTerms:
        return LatticeTerms(forbidden=tuple(self.enzymes.values()))

    def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation:
        hits = find_motifs(c, list(self.enzymes.values()))
        by_motif = {v: k for k, v in self.enzymes.items()}
        breaches = tuple(
            Breach(
                spec_id=self.id,
                interval=Interval(pos, pos + len(motif)),
                magnitude=1.0,
                message=f"{by_motif.get(motif, '?')} site {motif!r} at {pos}",
                detail={"enzyme": by_motif.get(motif, "?"), "motif": motif},
            )
            for motif, pos in hits
        )
        return Evaluation(
            spec_id=self.id,
            passes=not breaches,
            raw_score=float(len(breaches)),
            breaches=breaches,
            n_evaluated=c.length,
        )
