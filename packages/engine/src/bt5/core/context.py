"""Design context: the three simultaneous constraint sets.

The compound case that justifies the whole product -- a plasmid propagated in
E. coli, packaged into virus in a producer line, then transduced into a target
cell -- means THREE constraint sets apply to ONE sequence at the same time.

They are never collapsed. Each Spec is gated per slot, so every Breach knows
which context produced it, and two slots demanding opposite things produce two
Breaches over the same Interval. That is exactly the input the conflict detector
needs, which makes surfacing conflicts a CONSEQUENCE OF THE DATA MODEL rather
than a feature somebody has to remember to build.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Literal

from bt5.core.types import Strand

SlotRole = Literal["propagation", "producer", "target"]


class Modality(StrEnum):
    """How the construct is delivered. Drives enforcement_for()."""

    PLASMID_TRANSIENT = "plasmid_transient"
    PLASMID_STABLE = "plasmid_stable"
    LENTIVIRAL = "lentiviral"
    AAV = "aav"
    IVT_MRNA = "ivt_mrna"
    GENOME_INTEGRATED = "genome_integrated"
    BACTERIAL_EXPRESSION = "bacterial_expression"


class HostId(StrEnum):
    """Deeply curated hosts first; the rest are codon-table-level support."""

    E_COLI_K12 = "e_coli_k12"
    E_COLI_BL21 = "e_coli_bl21"
    HUMAN = "human"
    HEK293 = "hek293"
    CHO = "cho"
    S_CEREVISIAE = "s_cerevisiae"
    P_PASTORIS = "p_pastoris"
    SF9 = "sf9"
    MOUSE = "mouse"


#: NCBI genetic code table locked per host. A host cannot be used with a
#: mismatched table -- see ContextSlot.__post_init__. Table 11 is bacterial/plastid,
#: table 1 the standard code.
LOCKED_TRANSLATION_TABLE: Mapping[HostId, int] = {
    HostId.E_COLI_K12: 11,
    HostId.E_COLI_BL21: 11,
    HostId.HUMAN: 1,
    HostId.HEK293: 1,
    HostId.CHO: 1,
    HostId.S_CEREVISIAE: 1,
    HostId.P_PASTORIS: 1,
    HostId.SF9: 1,
    HostId.MOUSE: 1,
}


@dataclass(frozen=True, slots=True)
class ContextSlot:
    """One of the (up to three) simultaneous contexts."""

    role: SlotRole
    host: HostId
    modality: Modality
    table_id: int  # EXPLICIT. Never defaulted. See TranslationUnit for why.
    #: Which strand this slot's directional models read, RELATIVE TO THE
    #: CASSETTE -- not to the plasmid. +1 is the cassette's own sense strand,
    #: which is what a producer or target slot wants: an internal polyA signal
    #: matters on the RNA that actually gets made. Resolve it into construct
    #: coordinates with `spec.strand_for(ctx, slot)`, which composes it with
    #: `DesignContext.cassette_orientation`; reading this field directly is the
    #: bug that makes a reverse-cloned cassette's polyA analysis come back clean.
    strand_of_interest: Strand = 1
    enabled: bool = True

    def __post_init__(self) -> None:
        locked = LOCKED_TRANSLATION_TABLE.get(self.host)
        if locked is not None and self.table_id != locked:
            raise ValueError(
                f"host {self.host} is locked to NCBI translation table {locked}, "
                f"got {self.table_id}. A mismatched table silently changes the protein "
                f"(NCBI table 12 reassigns CTG=Ser rather than Leu)."
            )


@dataclass(frozen=True, slots=True)
class BiosecurityVerdict:
    """Result of screening the INPUT PROTEIN, before optimization.

    Protein-level screening is the one layer BT5's own output cannot defeat:
    the app's core function is producing a functionally identical sequence with
    maximally different nucleotides, which is the textbook method for evading
    nucleotide-homology screening.

    `status` is never "clear" when screening did not run. A NullScreen reports
    "not_run" so the report cannot imply a clean result that was never obtained.
    """

    status: Literal["not_run", "clear", "flag", "block"]
    database_version: str | None = None
    detail: str = ""

    @property
    def may_proceed(self) -> bool:
        return self.status != "block"


@dataclass(frozen=True)
class DesignContext:
    """Everything the optimizer needs that is not the sequence itself."""

    slots: tuple[ContextSlot, ...]
    cassette_orientation: Strand
    seed: int
    screen: BiosecurityVerdict
    strict_biosecurity: bool = True
    engine_versions: Mapping[str, str] = field(default_factory=dict)
    weights: Mapping[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.slots:
            raise ValueError("a DesignContext needs at least one ContextSlot")
        if len(self.slots) > 3:
            raise ValueError(f"at most three context slots, got {len(self.slots)}")
        roles = [s.role for s in self.slots]
        if len(set(roles)) != len(roles):
            raise ValueError(f"duplicate context slot roles: {roles}")

    @property
    def active_slots(self) -> tuple[ContextSlot, ...]:
        return tuple(s for s in self.slots if s.enabled)

    def slot(self, role: SlotRole) -> ContextSlot | None:
        return next((s for s in self.slots if s.role == role), None)
