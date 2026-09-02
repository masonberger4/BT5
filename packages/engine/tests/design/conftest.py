"""Shared fixtures for the design lane's tests.

The 140-residue protein and the synthetic MCS backbone are the same pair the
walking skeleton's tests used; they are lifted here because the ranking tests,
the order-file tests and the timing test all need them and none of them should
own the other's fixture.

Speed is a fixture concern in this lane and nowhere else. A full-default
`design()` sweeps one weight vector per live axis (3-4) and draws a 200-variant
null per cheap objective; the null is the expensive half, and it is what `fast`
shrinks. `FAST_SWEEP_STEPS` equals `DEFAULT_SWEEP_STEPS` because the shipped
lattice is already minimal -- the constant is kept so a future rise in the
default does not silently slow every unit test. Shrinking the null is a
test-runtime decision, never a claim about the product:
`test_the_shipped_defaults_are_the_ones_measured` pins the real values so this
cannot quietly become the default.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest
from bt5.core.context import HostId, Modality
from bt5.core.registry import all_specs, discover
from bt5.design import design
from bt5.vector import read_genbank
from bt5.vector.backbone import VectorBackbone

ROOT = Path(__file__).resolve().parents[4]
MCS_PATH = ROOT / "tests" / "data" / "backbones" / "synthetic_mcs_ef1a.gb"

#: A 140-residue protein, initiator first. Long enough to fill the fragment.
PROTEIN = (
    "MKLVTAAFERSKSVQNYVVSTKDSPLYYLRKWVRSGYKFDCEEVGLREHQGPAATYTPTQAIWRLTLPSPLL"
    "NVDVWQNSCKSLQHTASWKKHRFGLFTLVISPLIRLGEVASLCGLCEHTATSEVKVCPIDCLQSPTSF"
)

#: The 20 sense residues, in a fixed order so a seeded draw is reproducible.
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

#: One sweep step and a 12-variant null: the smallest values that still run the
#: sweep, the selection, the null and the percentile. NOT the shipped defaults.
FAST_SWEEP_STEPS = 1
FAST_NULL_N = 12

#: PLAN's G7 bar: 500 residues, end to end, in 10 s.
G7_PROTEIN_LENGTH = 500
G7_SECONDS = 10.0


def protein_of(length: int, *, seed: int = 0) -> str:
    """A pseudo-random protein of `length` residues, initiator first.

    Seeded explicitly with `np.random.default_rng` (CLAUDE.md 3.7); an unseeded
    draw here would make a timing failure irreproducible, which is the one thing
    a timing test cannot afford.

    Pseudo-random rather than a repeat of `PROTEIN`, because a tiled protein
    back-translates into a tandem array and would measure Tier B's repair of an
    artefact of the fixture rather than the cost of designing a 500-residue
    protein.
    """
    rng = np.random.default_rng(seed)
    return "M" + "".join(AMINO_ACIDS[i] for i in rng.integers(0, len(AMINO_ACIDS), length - 1))


def fast_null_sizes() -> dict[str, int]:
    """`FAST_NULL_N` for EVERY registered rule, derived rather than listed.

    This was a hand-written tuple of objective ids, and it rotted the moment the
    catalog grew: when b8, c3, d8 and e3 merged, each fell through to
    `NULL_N_BY_COST["cheap"]` and drew 200 variants instead of 12 -- in every
    test using the `fast` fixture. Since the memoised assembler builds `max(n)`
    constructs, that also took the assemblies from 12 to 200 per design. The
    lane's tests went from 75 s to over ten minutes, and nothing failed to say
    so; a slow fixture is invisible until someone times it.

    Deriving from `all_specs()` cannot drift: a rule added tomorrow is shrunk on
    the day it lands. Ids that are not scored in a given context are simply
    unused -- `null_size` only consults the override for objectives it is asked
    about -- so over-covering is free and under-covering is the bug.
    """
    discover()
    return {cls.id: FAST_NULL_N for cls in all_specs()}


@pytest.fixture
def backbone() -> VectorBackbone:
    return read_genbank(MCS_PATH)


@pytest.fixture
def protein_500() -> str:
    """The 500-residue protein gate G7 is stated at.

    A fixture rather than an import: `packages/engine/tests/design/` has no
    `__init__.py`, so `from conftest import protein_of` resolves through pytest's
    sys.path insertion and can pick up a different lane's conftest.
    """
    return protein_of(G7_PROTEIN_LENGTH)


@pytest.fixture
def protein() -> str:
    """The protein every design-lane test designs, as a fixture rather than a
    cross-file import: `packages/engine/tests/design/` is not a package, so
    `from .conftest import PROTEIN` does not resolve."""
    return PROTEIN


@pytest.fixture
def native_cds() -> str:
    """A synonymous CDS for PROTEIN, standing in for a real wild-type one.

    Built by the reference back-translator, so it arrives from OUTSIDE the
    optimizer -- which is what a native CDS is. It is the CALLER supplying it;
    `design()` still has no way to invent one.
    """
    from bt5.codon.tables import FileTableProvider as _Provider
    from bt5.solver.reference import back_translate

    return back_translate(PROTEIN, _Provider().genetic_code(1))


@pytest.fixture
def fast() -> Any:
    """`design()` with the sweep and the null shrunk to test size."""

    def run(bb: VectorBackbone, **kw: Any) -> Any:
        params: dict[str, Any] = {
            "protein": PROTEIN,
            "table_id": 1,
            "modality": Modality.LENTIVIRAL,
            "hosts": [HostId.HEK293],
            "sweep_steps": FAST_SWEEP_STEPS,
            "null_sizes": fast_null_sizes(),
        }
        params.update(kw)
        return design(backbone=bb, **params)

    return run
