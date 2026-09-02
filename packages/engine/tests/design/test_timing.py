"""Gate G7: a 500-residue design, end to end, in under 10 s.

Deliberately NOT marked `slow`. `scripts/gates.sh` and `ci.yml`'s engine job both
run `-m "not slow"`, so a `slow` mark here would make this the first use of a
marker whose only effect is that the test never runs -- a timing gate nothing
executes is worse than no timing gate.

Equally deliberately, this is a plain `pytest` assertion and not a benchmark.
`benchmarks/` does not exist in this repo despite `CLAUDE.md` 2 protecting two
files inside it; creating it is an owner decision under
`approved:algorithm-change` and is not this lane's to make.

A failure here is a FINDING about the design budget, not a flake. The budget is
spent in three places and each is tunable without touching the bar:
`DEFAULT_SWEEP_STEPS` (weight vectors, so full solves), `NULL_N_BY_COST`
(variants scored per objective) and `DEFAULT_GALLERY_SIZE` (candidates scored
and verified). If this goes red, one of those three is the answer -- never the
10 s.
"""

from __future__ import annotations

import time

from bt5.core.context import HostId, Modality
from bt5.design import design
from bt5.vector.backbone import VectorBackbone

#: PLAN's G7 bar: 500 residues, end to end, in 10 s. Duplicated from the conftest
#: rather than imported, so the number this test asserts is visible in the file
#: that asserts it.
G7_SECONDS = 10.0
G7_PROTEIN_LENGTH = 500


def test_a_500_residue_design_meets_the_g7_budget(
    backbone: VectorBackbone, protein_500: str
) -> None:
    assert len(protein_500) == G7_PROTEIN_LENGTH
    started = time.perf_counter()
    result = design(
        backbone=backbone,
        protein=protein_500,
        table_id=1,
        modality=Modality.LENTIVIRAL,
        hosts=[HostId.HEK293],
    )
    elapsed = time.perf_counter() - started

    # The design must be real, or the timing measures nothing.
    assert result.result.candidates
    assert result.genbank
    assert result.order_csv
    assert elapsed < G7_SECONDS, (
        f"a {G7_PROTEIN_LENGTH}-residue design took {elapsed:.1f} s against PLAN's "
        f"G7 bar of {G7_SECONDS:.0f} s. Tune DEFAULT_SWEEP_STEPS, NULL_N_BY_COST "
        f"or DEFAULT_GALLERY_SIZE -- never the bar."
    )
