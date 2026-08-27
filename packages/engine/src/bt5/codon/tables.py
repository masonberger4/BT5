"""Genetic codes, codon usage, and CAI.

CAI is implemented here as a DESCRIPTIVE STATISTIC and a soft band, never as a
maximization target. Kudla 2009 measured r = 0.14 (not significant) for CAI
against expression across 154 synonymous GFP variants, while 5' folding energy
explained 44% of the variance; Welch 2009 states outright that "CAI has no value
in predicting gene expression". Maximizing it collapses each amino acid to a
single codon, which produces perfect nucleotide repeats in exactly the repetitive
proteins people actually express.

Convention, pinned so our numbers are reproducible and comparable:
  - geometric mean of relative adaptiveness w over SENSE codons
  - ATG and TGG excluded (single-codon families carry no information)
  - stop codons excluded
  - pseudocount 0.5 for zero-count codons when building w from counts
Sharp & Li 1987: https://pubmed.ncbi.nlm.nih.gov/3547335/
"""

from __future__ import annotations

import json
import math
from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import Literal, cast

from Bio.Data import CodonTable


def _repo_root() -> Path:
    """Walk up to the directory holding pyproject.toml.

    A fixed `parents[n]` count is brittle -- it silently resolves to the wrong
    directory the moment the package moves, and the failure surfaces as a
    confusing FileNotFoundError rather than as a path bug.
    """
    for candidate in Path(__file__).resolve().parents:
        if (candidate / "pyproject.toml").exists():
            return candidate
    raise RuntimeError("could not locate the repository root (no pyproject.toml found)")


DATA_DIR = _repo_root() / "data"

#: Single-codon families contribute nothing to CAI and are excluded by convention.
SINGLE_CODON_AA = ("M", "W")


@dataclass(frozen=True)
class NcbiGeneticCode:
    """A genetic code, with stop codons excluded from every synonymous set.

    That exclusion is not a detail: NCBI tables 27 and 28 make TGA both Trp AND a
    stop, so a back-translator that picks TGA for Trp terminates translation
    mid-protein.
    """

    table_id: int

    @property
    def _table(self) -> CodonTable.CodonTable:
        return cast("CodonTable.CodonTable", CodonTable.unambiguous_dna_by_id[self.table_id])

    @property
    def stop_codons(self) -> tuple[str, ...]:
        return tuple(self._table.stop_codons)

    @property
    def start_codons(self) -> tuple[str, ...]:
        return tuple(self._table.start_codons)

    def translate(self, dna: str) -> str:
        if len(dna) % 3:
            raise ValueError(f"length {len(dna)} is not a multiple of 3")
        t = self._table
        out: list[str] = []
        for i in range(0, len(dna), 3):
            codon = dna[i : i + 3]
            out.append("*" if codon in t.stop_codons else t.forward_table[codon])
        return "".join(out)

    @cache  # noqa: B019
    def _synonymous(self) -> Mapping[str, tuple[str, ...]]:
        t = self._table
        stops = set(t.stop_codons)
        by_aa: dict[str, list[str]] = {}
        for codon, aa in sorted(t.forward_table.items()):
            if codon in stops:
                continue  # never offer a codon that also terminates translation
            by_aa.setdefault(aa, []).append(codon)
        return {aa: tuple(c) for aa, c in by_aa.items()}

    def families(self) -> Mapping[str, tuple[str, ...]]:
        """Every amino acid mapped to its non-stop synonymous codons."""
        return self._synonymous()

    def synonymous_codons(self, aa: str) -> tuple[str, ...]:
        try:
            return self._synonymous()[aa]
        except KeyError as exc:
            raise ValueError(
                f"amino acid {aa!r} has no non-stop codon in NCBI table {self.table_id}"
            ) from exc

    def is_start(self, codon: str) -> bool:
        return codon in self._table.start_codons

    def is_stop(self, codon: str) -> bool:
        return codon in self._table.stop_codons


@dataclass(frozen=True)
class CodonUsage:
    """Relative adaptiveness w per sense codon, plus its provenance.

    `reference_set` is load-bearing for reproducibility: CAI is meaningless
    without it, and swapping the reference set raises every CAI value, which a
    benchmark would happily read as an improvement.
    """

    host: str
    reference_set: str
    w: Mapping[str, float]
    citation_url: str = ""

    @classmethod
    def from_counts(
        cls,
        host: str,
        reference_set: str,
        counts: Mapping[str, int],
        code: NcbiGeneticCode,
        pseudocount: float = 0.5,
    ) -> CodonUsage:
        """Build w from raw codon counts over a highly-expressed reference set."""
        w: dict[str, float] = {}
        for codons in code.families().values():
            # Pseudocount before taking the family maximum, so a codon absent
            # from a small reference set gets a small non-zero w rather than
            # dropping out of the geometric mean entirely.
            adjusted = {c: counts.get(c, 0) + pseudocount for c in codons}
            peak = max(adjusted.values())
            for c, v in adjusted.items():
                w[c] = v / peak
        return cls(host=host, reference_set=reference_set, w=w)

    def cai(self, dna: str, code: NcbiGeneticCode) -> float:
        """Geometric mean of w over informative sense codons.

        Returns 0.0 for a sequence with no informative codons rather than raising,
        so a poly-Met or poly-Trp peptide degrades gracefully.
        """
        if len(dna) % 3:
            raise ValueError(f"length {len(dna)} is not a multiple of 3")
        logs: list[float] = []
        for i in range(0, len(dna), 3):
            codon = dna[i : i + 3]
            if code.is_stop(codon):
                continue
            aa = code.translate(codon)
            if aa in SINGLE_CODON_AA:
                continue
            weight = self.w.get(codon)
            if weight is None or weight <= 0.0:
                continue
            logs.append(math.log(weight))
        if not logs:
            return 0.0
        return math.exp(sum(logs) / len(logs))


class FileTableProvider:
    """Loads genetic codes from Biopython and usage tables from data/."""

    def __init__(self, data_dir: Path | None = None) -> None:
        self.data_dir = data_dir or DATA_DIR

    def genetic_code(self, table_id: int) -> NcbiGeneticCode:
        if table_id not in CodonTable.unambiguous_dna_by_id:
            raise ValueError(f"unknown NCBI translation table {table_id}")
        return NcbiGeneticCode(table_id)

    @cache  # noqa: B019
    def usage(self, host: str) -> CodonUsage:
        path = self.data_dir / "codon_usage" / f"{host}.json"
        if not path.exists():
            raise FileNotFoundError(f"no codon usage table for host {host!r} at {path}")
        raw = json.loads(path.read_text())
        prov = raw.get("_provenance", {})
        return CodonUsage(
            host=host,
            reference_set=prov.get("name", host),
            w=raw["w"],
            citation_url=prov.get("citation_url", ""),
        )

    def weights(self, host: str, kind: Literal["cai", "tai", "stai", "csc"]) -> Mapping[str, float]:
        if kind != "cai":
            raise NotImplementedError(f"{kind} weights are a later lane (M5 codon)")
        return self.usage(host).w
