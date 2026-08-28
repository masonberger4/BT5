"""ViennaRNA behind the `FoldEngine` protocol.

Every literature threshold BT5 uses was calibrated against ViennaRNA with a
specific parameter set, so the engine is pinned exactly and a version change is
a SCIENTIFIC change rather than a dependency bump -- it shifts every baseline and
makes old numbers incomparable to new ones. `CALIBRATED_VERSION` is what the
thresholds were measured against; `ViennaFold.version` is what is installed.
They are reported separately and never reconciled silently, because a guard that
hard-fails a comparison across differing versions can only exist if both numbers
survive to the point of comparison.

Three things measured here rather than assumed, on ViennaRNA 2.7.2:

  Temperature and dangles move the answer by 23-33%. 37 C gives -8.50 kcal/mol
  on a test sequence where 30 C gives -11.35, and dangles=0 gives -6.90 where
  dangles=2 gives -8.50. That is why `FoldEnergy` carries both, and why this
  engine refuses to produce a number without stamping them on it.

  ViennaRNA reads T as U silently, so DNA input does not error -- it just works.
  The conversion here is therefore not load-bearing for correctness today; it is
  explicit so that a future version which stops doing that is caught by a test
  rather than by a shifted baseline.

  Whole-sequence MFE is O(n^3) and slower than docs/PLAN.md estimates: 0.77 s at
  1 kb where the plan says ~0.24 s, and 8.9 s at 3 kb where it says ~6.5 s. So
  `mfe` is report-time only, and `mfe_window` is the primitive everything in a
  loop must use.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, ClassVar

from bt5.core.services import FoldEnergy
from bt5.core.types import Interval, reverse_complement

#: Matches `Spec.engine_calibration`'s documented form, "viennarna:rna_turner2004",
#: so `FoldEnergy.calibration_key` and a rule's declared calibration compare equal.
ENGINE_NAME = "viennarna"
PARAM_SET = "rna_turner2004"

#: The version BT5's thresholds were calibrated against. Pinned in pyproject as
#: viennarna==2.7.2. A bump carries `approved:algorithm-change` and regenerates
#: the baseline; see CLAUDE.md section 6.
CALIBRATED_VERSION = "2.7.2"

#: RNAplfold-style accessibility parameters. These change the number, so they are
#: named constants rather than call-site literals -- but note the `accessibility`
#: protocol returns a bare float, so unlike a FoldEnergy they cannot travel with
#: their result. Anything persisting an accessibility value must record them.
ACCESSIBILITY_WINDOW = 80
ACCESSIBILITY_MAX_BP_SPAN = 40


class FoldUnavailableError(RuntimeError):
    """ViennaRNA is not installed, so no folding number can be produced.

    Raised rather than returning a sentinel value. A fold objective that cannot
    be evaluated must degrade visibly -- see `load_fold_engine`, which returns
    None so the caller records a degradation instead of scoring a fake number.
    """


def _rna() -> Any:
    try:
        # ViennaRNA ships no py.typed marker. Confined to this one line, the
        # way the vector lane confines Biopython's untyped constructors, so no
        # caller and no config file has to carry it.
        import RNA  # type: ignore[import-untyped]
    except ImportError as exc:  # pragma: no cover - exercised only without the extra
        raise FoldUnavailableError(
            "ViennaRNA is not installed. It is declared as an optional extra: "
            'install with `uv pip install -e ".[fold]"`. Without it the 5\' folding '
            "objective is unavailable and must be reported as such, never estimated."
        ) from exc
    return RNA


def vienna_available() -> bool:
    """True when a folding number can actually be produced."""
    try:
        _rna()
    except FoldUnavailableError:
        return False
    return True


def installed_version() -> str:
    """The installed ViennaRNA version, or "" when it is absent."""
    try:
        return str(_rna().__version__)
    except FoldUnavailableError:
        return ""


def is_calibrated() -> bool:
    """True when the installed engine is the one the thresholds were measured on.

    Reported, not enforced. The benchmark comparability guard is what hard-fails
    across differing versions; this is how it knows.
    """
    return installed_version() == CALIBRATED_VERSION


@dataclass(frozen=True, slots=True)
class ViennaFold:
    """`FoldEngine` over ViennaRNA. Conditions are fields, not globals."""

    temperature_c: float = 37.0
    dangles: int = 2

    name: ClassVar[str] = ENGINE_NAME
    param_set: ClassVar[str] = PARAM_SET

    @property
    def version(self) -> str:
        return installed_version()

    # -- the protocol -----------------------------------------------------

    def mfe(self, seq: str) -> FoldEnergy:
        """Whole-sequence MFE. REPORT TIME ONLY -- O(n^3), ~8.9 s at 3 kb."""
        rna = _rna()
        _, dg = rna.fold_compound(_prepare(seq), self._md()).mfe()
        return self._energy(dg)

    def mfe_window(self, seq: str, iv: Interval) -> FoldEnergy:
        """Windowed fold: the interactive-loop and null-model primitive.

        Wrap-aware and strand-aware, because the 5' window of a reverse-oriented
        cassette sits at higher coordinates and may straddle the origin. The
        slice is taken here rather than by the caller so that every consumer gets
        the same answer for the same interval.
        """
        rna = _rna()
        _, dg = rna.fold_compound(_prepare(slice_of(seq, iv)), self._md()).mfe()
        return self._energy(dg)

    def accessibility(self, seq: str, iv: Interval, u: int) -> float | None:
        """Mean probability that a stretch of `u` bases is unpaired, over `iv`.

        Returns None when the window is shorter than the stretch being asked
        about, which is a question with no answer rather than a zero.
        """
        if u < 1:
            raise ValueError(f"stretch length must be positive, got {u}")
        rna = _rna()
        text = _prepare(slice_of(seq, iv))
        if len(text) < u:
            return None
        # 1-indexed: up[i][u] is P(the u bases ending at i are unpaired).
        up = rna.pfl_fold_up(text, u, ACCESSIBILITY_WINDOW, ACCESSIBILITY_MAX_BP_SPAN)
        values = [up[i][u] for i in range(u, len(text) + 1)]
        return sum(values) / len(values) if values else None

    # -- internals --------------------------------------------------------

    def _md(self) -> Any:
        md = _rna().md()
        md.temperature = self.temperature_c
        md.dangles = self.dangles
        return md

    def _energy(self, dg: float) -> FoldEnergy:
        return FoldEnergy(
            dg_kcal_mol=float(dg),
            engine=ENGINE_NAME,
            engine_version=installed_version(),
            param_set=PARAM_SET,
            temperature_c=self.temperature_c,
            dangles=self.dangles,
        )


def slice_of(seq: str, iv: Interval) -> str:
    """Wrap- and strand-aware extraction, matching `Construct.slice`.

    Duplicated deliberately: this engine takes a bare `str` because that is the
    protocol's shape, so it cannot call Construct.slice -- and two extractions
    that disagree about wrapping would put the fold window somewhere other than
    where the caller believes it is.
    """
    n = len(seq)
    # end > n means the interval wraps the origin -- the one representation BT5
    # uses everywhere. A bare str carries no topology, so wrapping is the only
    # reading available, and the caller is responsible for not asking for a
    # wrapped window on a linear molecule.
    wrapped = iv.end > n
    sub = seq[iv.start :] + seq[: iv.end - n] if wrapped else seq[iv.start : iv.end]
    return reverse_complement(sub) if iv.strand == -1 else sub


def _prepare(seq: str) -> str:
    """Upper-case and transcribe. ViennaRNA folds RNA; BT5 carries DNA."""
    text = seq.upper().replace("T", "U")
    bad = set(text) - {"A", "C", "G", "U"}
    if bad:
        raise ValueError(f"cannot fold non-ACGU characters: {sorted(bad)}")
    return text


def load_fold_engine(*, temperature_c: float = 37.0, dangles: int = 2) -> ViennaFold | None:
    """The engine, or None when folding is unavailable.

    None rather than a stub that returns plausible numbers. Every threshold in
    BT5 is calibrated in kcal/mol, so a stub's output would flow through the
    scorers, the null and the percentile unchallenged and come out the far end
    as a confident rank. The caller records `degradation_reason()` in
    `Provenance.degradations` and reports the objective as unavailable.
    """
    if not vienna_available():
        return None
    return ViennaFold(temperature_c=temperature_c, dangles=dangles)


def degradation_reason() -> str | None:
    """What to record in `Provenance.degradations`, or None when all is well."""
    if not vienna_available():
        return (
            "ViennaRNA is not installed, so folding objectives were not evaluated; "
            "no structure-derived score is included in this result"
        )
    if not is_calibrated():
        return (
            f"ViennaRNA {installed_version()} is installed but BT5's thresholds were "
            f"calibrated against {CALIBRATED_VERSION}; energies are comparable within "
            f"this run but not against baselines taken on another version"
        )
    return None
