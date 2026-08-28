"""Injected services.

Rules receive `Services`; they never import the codon, structure or vector lanes.
This is what lets the rules lane build on day 1 against the skeleton's providers
and become automatically correct on day 30 against the real ones, WITH AN EMPTY
DIFF in the rules lane.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar, Literal, Protocol

if TYPE_CHECKING:
    import numpy as np

    from bt5.core.types import Construct, Interval


@dataclass(frozen=True, slots=True)
class FoldEnergy:
    """A folding free energy that cannot travel without its provenance.

    Every literature threshold in BT5 -- the Boel -39 kcal/mol dual gate, the
    cap-proximal -30/-50/-60 ladder -- was calibrated against a specific engine
    and parameter set. Applying a ViennaRNA-calibrated number to another engine's
    output is the single most likely correctness bug in the folding feature, so
    the number carries its ruler everywhere, including into goldens.
    """

    dg_kcal_mol: float
    engine: str
    engine_version: str
    param_set: str
    temperature_c: float = 37.0
    dangles: int = 2
    #: Dot-bracket structure, or "" when the engine was not asked for one. An
    #: energy alone cannot answer a rule stated in terms of an individual
    #: hairpin's position -- the cap-proximal ladder, the -1 frameshift
    #: escalation and the AAV hairpin flag are all written that way -- so the
    #: structure has to be able to travel with the number rather than requiring
    #: a second fold to recover it.
    structure: str = ""
    #: Index in `structure` where the SECOND molecule begins, for an
    #: intermolecular duplex; None for an ordinary intramolecular fold.
    #: ViennaRNA's dimer structure carries no separator, so without this the two
    #: halves of a duplex are indistinguishable in the string -- and a duplex
    #: energy silently read as an intramolecular one is off by the energy of a
    #: covalent join that does not exist.
    duplex_split: int | None = None

    @property
    def calibration_key(self) -> str:
        return f"{self.engine}:{self.param_set}"

    @property
    def is_duplex(self) -> bool:
        return self.duplex_split is not None


class FoldEngine(Protocol):
    """RNA/DNA secondary structure. ViennaRNA is the bundled default."""

    name: ClassVar[str]
    version: ClassVar[str]
    param_set: ClassVar[str]

    def mfe(self, seq: str) -> FoldEnergy:
        """Whole-sequence minimum free energy. REPORT TIME ONLY.

        ViennaRNA is ~0.24 s at 1 kb and ~6.5 s at 3 kb, so this must never run
        inside the interactive loop or the empirical null.
        """
        ...

    def mfe_window(self, seq: str, iv: Interval) -> FoldEnergy:
        """Windowed fold. This is the interactive-loop and null-model primitive."""
        ...

    def accessibility(self, seq: str, iv: Interval, u: int) -> float | None:
        """Mean unpaired probability over a window, or None if unsupported."""
        ...

    def duplex(self, a: str, b: str) -> FoldEnergy:
        """Free energy of two molecules pairing with EACH OTHER.

        Distinct from folding their concatenation, which is not the same
        quantity: joining them end to end lets the model form a covalent bond
        that does not exist, and on a real Shine-Dalgarno / anti-SD pair that
        error is worth ~0.6 kcal/mol in the direction that makes a weak site
        look adequate.

        Required rather than optional because the contract has to be able to
        express the rules it was frozen to serve. The Salis TIR model -- one of
        the four families gate G0 validates the protocol against -- is defined
        in terms of dG_mRNA:rRNA and dG_standby, and neither `mfe` nor
        `mfe_window` can produce an intermolecular energy at any spelling.
        """
        ...


class KmerIndex(Protocol):
    """Repeat and uniqueness queries over the assembled construct.

    BIOSECURITY: constructed ONLY from a Construct. There is deliberately no
    constructor accepting an external sequence database, because pointing a
    homology-minimiser at an arbitrary target database turns BT5 into a
    general-purpose screening-evasion tool. A CI grep enforces this.
    """

    @classmethod
    def of(cls, c: Construct, k: int) -> KmerIndex: ...

    def duplicates(self, min_len: int) -> Iterator[tuple[Interval, Interval]]:
        """Direct repeat pairs at least `min_len` long."""
        ...

    def revcomp_pairs(self, min_stem: int, max_loop: int) -> Iterator[tuple[Interval, Interval]]:
        """Inverted repeats / potential hairpins."""
        ...


class GeneticCode(Protocol):
    table_id: int

    def translate(self, dna: str) -> str: ...
    def synonymous_codons(self, aa: str) -> tuple[str, ...]: ...
    def is_start(self, codon: str) -> bool: ...
    def is_stop(self, codon: str) -> bool: ...


class TableProvider(Protocol):
    """Genetic codes and host codon statistics."""

    def genetic_code(self, table_id: int) -> GeneticCode: ...
    def usage(self, host: str) -> Mapping[str, float]: ...
    def weights(
        self, host: str, kind: Literal["cai", "tai", "stai", "csc"]
    ) -> Mapping[str, float]: ...


@dataclass(frozen=True)
class Services:
    """Passed to every rule. Never imported around."""

    #: None when no folding engine is available. Widened deliberately: `mfe`
    #: returns a non-optional FoldEnergy, so an engine that cannot fold has only
    #: two ways to behave -- raise, which crashes the run, or invent a number,
    #: which is the one thing the whole honesty apparatus exists to prevent.
    #: A rule receiving None reports its objective unavailable; see
    #: `ObjectiveScore.unavailable`.
    fold: FoldEngine | None
    kmer: type[KmerIndex]
    tables: TableProvider
    rng: np.random.Generator  # always np.random.default_rng(seed); never global RNG
