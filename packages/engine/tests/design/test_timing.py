"""Gate G7: a 500-residue design, end to end, in under 10 s.

Deliberately NOT marked `slow`. `scripts/gates.sh` and `ci.yml`'s engine job both
run `-m "not slow"`, so a `slow` mark here would make this the first use of a
marker whose only effect is that the test never runs -- a timing gate nothing
executes is worse than no timing gate.

Equally deliberately, this is a plain `pytest` assertion and not a benchmark.
`benchmarks/` does not exist in this repo despite `CLAUDE.md` 2 protecting two
files inside it; creating it is an owner decision under
`approved:algorithm-change` and is not this lane's to make.

A failure here is a FINDING about the design budget, not a flake, and as of the
merge with main it is red for a reason OUTSIDE this lane. Measured at 500 aa on
the reference fixture: 23.5 s total, of which the sweep is 20.2 s and the null
3.1 s. The sweep is three full solves and each is dominated by Tier B repair.

`RuleSet.breach_finder()` (`solver/catalog.py`) is what repair calls once per
candidate, up to `max_candidates` per iteration. Its own docstring says it is
"SCOPED TO THE HARD_REPAIR RULES ... E8's k-mer index or B1's fold evaluated
here would be pure waste" -- but it returns `self.findings(c).repairable`, and
`findings()` walks every spec. Measured: 62.2 ms per candidate as implemented
against 5.1 ms for the 8 HARD_REPAIR specs it actually needs, a 12.2x waste,
with `f2_near_perfect_repeats` (19.7 ms) and `e8_kmer_uniqueness` (17.2 ms) paid
for and discarded every time. That is M1's lane, not this one.

The knobs this lane owns cannot close a 2.4x gap on their own:
`DEFAULT_SWEEP_STEPS` is already 1 (three vectors, the minimum for a panel),
dropping `DEFAULT_GALLERY_SIZE` below `MIN_GALLERY` stops being a gallery, and
`NULL_N_BY_COST` accounts for 3.1 s of the 23.5. Lowering `max_candidates` from
256 to 64 saves ~1.5 s per solve and buys repair quality with it.

So: never the bar, and not these knobs either. Reported rather than tuned.
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
