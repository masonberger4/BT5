"""Recorded instances of the frozen types, and the codec that revives them.

The amendment protocol calls for `test_backward_compat` re-parsing every
recorded fixture inside the amendment PR. This is what "re-parse" means for
BT5's contract: the fixture records the CONSTRUCTOR ARGUMENTS that built a value
under the frozen contract, and the test builds it again against today's code.

That is a sharper check than it sounds, because it fails on exactly the changes
that break real callers and passes on the ones that do not:

  renamed field        -> TypeError, unexpected keyword
  removed field        -> TypeError, unexpected keyword
  new REQUIRED field   -> TypeError, missing argument
  new defaulted field  -> constructs fine, which is the correct answer
  changed default      -> the recorded value still round-trips, and the
                          manifest catches the default itself

Fixtures are RECORDED, never hand-edited: `python tests/contract/regenerate.py`
rebuilds them from live objects. A hand-edited fixture is a fixture that agrees
with whatever broke it.
"""

from __future__ import annotations

import dataclasses
import enum
import importlib
import json
import math
from collections.abc import Mapping
from pathlib import Path
from typing import Any

FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures"


def _qualname(value: Any) -> str:
    cls = value if isinstance(value, type) else type(value)
    return f"{cls.__module__}.{cls.__qualname__}"


def encode(value: Any) -> Any:
    """A JSON-able record of `value`, keyed by constructor arguments."""
    if isinstance(value, enum.Enum):
        return {"__enum__": _qualname(value), "name": value.name}
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {
            "__type__": _qualname(value),
            "fields": {f.name: encode(getattr(value, f.name)) for f in dataclasses.fields(value)},
        }
    if isinstance(value, tuple):
        return {"__tuple__": [encode(v) for v in value]}
    if isinstance(value, Mapping):
        return {"__map__": {str(k): encode(v) for k, v in value.items()}}
    if isinstance(value, float) and not math.isfinite(value):
        # JSON has no NaN. ObjectiveScore.unavailable() fills raw and percentile
        # with NaN deliberately -- 0.0 would read as a real, terrible score --
        # so the codec has to carry it rather than round it away.
        return {"__float__": repr(value)}
    if isinstance(value, list):
        return [encode(v) for v in value]
    return value


def _resolve(path: str) -> Any:
    module, _, name = path.rpartition(".")
    obj = importlib.import_module(module)
    for part in name.split("."):
        obj = getattr(obj, part)
    return obj


def decode(record: Any) -> Any:
    """Rebuild a recorded value against TODAY's contract. May raise TypeError."""
    if isinstance(record, list):
        return [decode(v) for v in record]
    if not isinstance(record, dict):
        return record
    if "__enum__" in record:
        return getattr(_resolve(record["__enum__"]), record["name"])
    if "__tuple__" in record:
        return tuple(decode(v) for v in record["__tuple__"])
    if "__map__" in record:
        return {k: decode(v) for k, v in record["__map__"].items()}
    if "__float__" in record:
        return float(record["__float__"])
    if "__type__" in record:
        cls = _resolve(record["__type__"])
        return cls(**{k: decode(v) for k, v in record["fields"].items()})
    return record


def _construct() -> Any:
    """One assembled construct, exercising the geometry the contract is for.

    Deliberately circular, with a CDS that wraps the origin and a scan-exempt
    LTR: the cases where a linear reading of the coordinate model gives a
    plausible wrong answer are the ones worth freezing.
    """
    from bt5.core.types import (
        Construct,
        Feature,
        Interval,
        Provenance,
        Segment,
        SegmentKind,
        Topology,
        TranslationUnit,
    )

    seq = "ATGAAACCCGGGTTTACGTACGTACGTAAGCTTGGGCCCAAATTTCCCGGGTAA"
    return Construct(
        sequence=seq,
        topology=Topology.CIRCULAR,
        segments=(
            Segment(Interval(48, len(seq) + 12), SegmentKind.DESIGNABLE_CDS, "cds across origin"),
            Segment(Interval(12, 30), SegmentKind.WHITELISTED_REPEAT, "ltr"),
            Segment(Interval(30, 48), SegmentKind.BACKBONE, "vector"),
        ),
        translation_units=(
            TranslationUnit(
                table_id=11,
                codon_map=(Interval(48, 51), Interval(51, 54), Interval(54, 57)),
                protein="MK",
                has_terminal_stop=True,
                starts_at_initiator=True,
            ),
        ),
        features=(
            Feature(Interval(12, 30), "LTR", {"label": ("5'LTR",)}, uid="ltr5"),
            Feature(Interval(48, len(seq) + 12), "CDS", {"transl_table": ("11",)}, uid="cds"),
        ),
        annotations={"molecule_type": "ds-DNA", "topology": "circular"},
        provenance=Provenance(
            app_version="0.1.0",
            seed=42,
            engine_versions={"viennarna": "2.7.2"},
            codon_table_name="ncbi_11",
            constraint_set_hash="deadbeef",
            degradations=("no folding engine installed",),
        ),
    )


def specimens() -> dict[str, Any]:
    """One recorded value per frozen type that BT5 constructs.

    Protocols are absent on purpose: nothing constructs a `Spec` or a
    `FoldEngine`, so a fixture would record nothing. The manifest covers those,
    and covers them better -- an added protocol method is MAJOR there.
    """
    from bt5.core.context import (
        BiosecurityVerdict,
        ContextSlot,
        DesignContext,
        HostId,
        Modality,
    )
    from bt5.core.result import (
        Candidate,
        Conflict,
        DesignResult,
        InfeasibilityCertificate,
        ObjectiveScore,
        Relaxation,
        ScoreCard,
    )
    from bt5.core.services import FoldEnergy
    from bt5.core.spec import Breach, Citation, Evaluation, LatticeTerms
    from bt5.core.types import Interval

    construct = _construct()
    slot = ContextSlot("producer", HostId.HEK293, Modality.LENTIVIRAL, 1, -1)
    breach = Breach(
        spec_id="d1_restriction_sites",
        interval=Interval(50, 56),
        magnitude=1.0,
        message="HindIII site 'AAGCTT' at 50",
        fixable_by_codon_choice=True,
        slot_role="producer",
        detail={"enzyme": "HindIII", "motif": "AAGCTT"},
    )
    scored = ObjectiveScore(
        spec_id="c1_cai",
        raw=0.78,
        unit="au",
        percentile=0.82,
        null_n=200,
        null_mean=0.71,
        null_sd=0.04,
    )
    unavailable = ObjectiveScore.unavailable("b1_five_prime", "kcal/mol", "no folding engine")
    candidate = Candidate(
        label="candidate_1",
        construct=construct,
        cds="ATGAAATAA",
        scorecard=ScoreCard(scores=(scored, unavailable), hard_checks=(breach,), total=0.82),
        design_hash="a1b2c3d4e5f6",
        codon_distance_to={"native_baseline": 0.31},
    )

    return {
        "interval_wrapping": Interval(48, 66, 1),
        "interval_reverse": Interval(10, 40, -1),
        "construct_circular": construct,
        "context_slot_reverse": slot,
        "design_context": DesignContext(
            slots=(slot,),
            cassette_orientation=-1,
            seed=42,
            screen=BiosecurityVerdict("not_run", None, "screening did not run"),
            strict_biosecurity=True,
            engine_versions={"viennarna": "2.7.2"},
            weights={"c1_cai": 0.2},
        ),
        "citation": Citation(
            "Kudla 2009", "https://pmc.ncbi.nlm.nih.gov/articles/PMC3902468/", 2009, "supports"
        ),
        "breach": breach,
        "evaluation": Evaluation(
            spec_id="e2_gc_band",
            passes=False,
            raw_score=0.63,
            breaches=(breach,),
            windows=((Interval(0, 50), 0.63),),
            n_evaluated=54,
            binding_side="upper",
        ),
        "lattice_terms": LatticeTerms(forbidden=("GAATTC", "GGATCC")),
        "fold_energy": FoldEnergy(
            dg_kcal_mol=-12.3,
            engine="viennarna",
            engine_version="2.7.2",
            param_set="rna_turner2004",
            temperature_c=37.0,
            dangles=2,
            structure="((((....))))",
            duplex_split=None,
        ),
        "objective_score": scored,
        "objective_unavailable": unavailable,
        "relaxation": Relaxation("e2_gc_band", "raise gc_max 0.60 -> 0.64", {"c1_cai": -0.02}),
        "conflict": Conflict(
            interval=Interval(48, 58),
            spec_ids=("b8_kozak", "d1_restriction_sites"),
            kind="mutually_exclusive",
            binding_spec_id="d1_restriction_sites",
            relaxations=(Relaxation("b8_kozak", "accept an adequate Kozak", {}),),
        ),
        "infeasibility_certificate": InfeasibilityCertificate(
            interval=Interval(48, 60),
            protein_span=(0, 4),
            minimal_conflicting_specs=("d1_restriction_sites", "e2_gc_band"),
            proof="empty_mutation_space",
            relaxations=(),
        ),
        "candidate": candidate,
        "design_result": DesignResult(
            candidates=(candidate,),
            native_baseline=candidate,
            conflicts=(),
            provenance=construct.provenance,
        ),
    }


def record() -> None:
    """Write every specimen to `fixtures/`. Run via regenerate.py."""
    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    for name, value in specimens().items():
        payload = {"recorded_as": _qualname(value), "value": encode(value)}
        path = FIXTURE_DIR / f"{name}.json"
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def load_all() -> dict[str, dict[str, Any]]:
    return {path.stem: json.loads(path.read_text()) for path in sorted(FIXTURE_DIR.glob("*.json"))}
