"""The solver seam, exercised against the REAL fixture and a REAL catalog rule.

Every test in `test_repair.py` drives `repair()` through a degenerate fake
assembler: one all-CDS segment at offset 0, no backbone, no flanks, one
homogeneous rule. That fake is blind to the whole class of defect this file
exists to catch -- a breach the solver must NOT chase, sitting in a backbone the
fake does not have.

So this file assembles a designed CDS into the synthetic lentiviral backbone
(the same fixture the vector lane round-trips) and runs the real `d4` internal
polyA rule over the assembled construct. That backbone carries polyA hexamers in
its two identical LTRs and its origin -- signals no codon choice can move -- and
the rule marks them `fixable_by_codon_choice=False`. Before the seam read that
flag, `target = max(breaches, ...)` picked the strongest of those backbone
signals, localized it to a window touching no editable codon, and aborted the
ENTIRE repair pass having done zero work, then raised a fabricated
`empty_mutation_space` certificate pointing 800 nt from the CDS.

The cross-lane import (vector + rules from a solver-lane test) is the one #58
sanctions for exactly this: the seam has no other consumer yet, and a fake
cannot prove the seam works on a real construct.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from bt5.codon.tables import FileTableProvider
from bt5.core.context import BiosecurityVerdict, ContextSlot, DesignContext, HostId, Modality
from bt5.core.services import Services
from bt5.core.spec import Breach
from bt5.core.types import Construct
from bt5.rules.catalog.d4_internal_polya import InternalPolyA
from bt5.solver.reference import back_translate
from bt5.solver.repair import repair
from bt5.vector import assemble, read_genbank
from bt5.vector.kmers import ConstructKmerIndex

ROOT = Path(__file__).resolve().parents[4]
BACKBONE_PATH = ROOT / "tests" / "data" / "backbones" / "synthetic_lenti_ef1a.gb"

#: The fixture's own transgene protein -- re-optimising it in place gives a
#: realistic assembled lentiviral construct with a real backbone on both sides.
PROTEIN = (
    "MLVTAAFARSKSVQNYVVSTKDSPLYYLRKWVRSGYKFDCEEVG"
    "LREHQGPAATYTPTQAIWRLTLPSPLLNVDVWQNSCKSLQHTASWKKHRFGLFTLVIS"
    "PLIRLGEVASLCGLCEHTATSEVKVCPIDCLQSPTSF"
)


@pytest.fixture(scope="module")
def code() -> object:
    return FileTableProvider().genetic_code(1)


def _lentiviral_context() -> DesignContext:
    """A producer slot on the lentiviral modality -- the one where an internal
    polyA is a measured functional-titer loss, so d4 escalates to HARD_REPAIR."""
    slot = ContextSlot(
        role="producer", host=HostId.HEK293, modality=Modality.LENTIVIRAL, table_id=1
    )
    return DesignContext(
        slots=(slot,), cassette_orientation=1, seed=0, screen=BiosecurityVerdict("not_run")
    )


class CountingD4Finder:
    """The real d4 rule as a BreachFinder, counting how often it is asked.

    UNFILTERED: it returns every d4 breach on the construct, including the
    unfixable ones in the LTRs and the origin. That is the point -- the seam
    must route those to `advisory` on its own, from the flag the rule sets.
    """

    def __init__(self, ctx: DesignContext, services: Services) -> None:
        self._rule = InternalPolyA()
        self._ctx = ctx
        self._services = services
        self.calls = 0

    def __call__(self, construct: Construct) -> tuple[Breach, ...]:
        self.calls += 1
        return self._rule.evaluate(construct, self._ctx, self._services).breaches


def _assembler(backbone: object, site: object):
    def assembler(cds: str) -> Construct:
        return assemble(backbone, cds, protein=PROTEIN, table_id=1, site=site).construct  # type: ignore[arg-type]

    return assembler


def test_an_unfixable_backbone_breach_no_longer_aborts_the_pass(code: object) -> None:
    """The regression. On the merge base this construct raised
    InfeasibleConstraints after ONE find_breaches call, the certificate pointing
    at an LTR base ~800 nt from the CDS. After E0.1 the pass reads
    `fixable_by_codon_choice`, carries the five backbone signals as advisories,
    and clears the four that lie in the CDS."""
    backbone = read_genbank(BACKBONE_PATH)
    site = backbone.find_insertion_site(label="transgene")
    services = Services(
        fold=None, kmer=ConstructKmerIndex, tables=FileTableProvider(), rng=np.random.default_rng(0)
    )
    finder = CountingD4Finder(_lentiviral_context(), services)
    assembler = _assembler(backbone, site)

    cds = back_translate(PROTEIN, code)  # type: ignore[arg-type]

    # The situation the seam must survive: real d4 breaches, some unfixable.
    breaches = finder(assembler(cds))
    finder.calls = 0
    unfixable = [b for b in breaches if not b.fixable_by_codon_choice]
    assert unfixable, "fixture must carry backbone polyA signals no codon can move"

    outcome = repair(
        cds,
        PROTEIN,
        code,  # type: ignore[arg-type]
        assemble=assembler,
        find_breaches=finder,
        forbidden=(),
        window=50,
        seed=0,
    )

    # It did NOT abort at the first breach: the zero-work abort called
    # find_breaches exactly once before raising.
    assert finder.calls > 1, "the pass must search, not abort on an unfixable backbone breach"
    assert outcome.clean, "the actionable (in-CDS) polyA signals are all removable"
    assert not outcome.remaining

    # The five unfixable backbone signals are carried, untouched, as advisories.
    assert len(outcome.advisory) == 5
    assert {b.spec_id for b in outcome.advisory} == {"d4_internal_polya"}
    assert all(not b.fixable_by_codon_choice for b in outcome.advisory)

    # And the protein is preserved, which is the whole non-negotiable.
    assert code.translate(outcome.cds)[:-1] == PROTEIN  # type: ignore[attr-defined]
