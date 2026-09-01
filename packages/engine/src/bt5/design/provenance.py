"""Everything needed to reproduce a design, and the two content hashes.

`Provenance` has been a field of the result types since Wave 1, but nothing in
`src/` ever built one -- the design lane is its first home. Two hashes matter and
they answer different questions:

- `design_hash` goes on the tube label and the GenBank note. Two runs that
  produce two different sequences under one name is how a lab ends up with two
  tubes and an irreproducible result, so it hashes the emitted CDS -- and the
  BACKBONE, because the same CDS spliced into two different vectors is two
  different constructs and must not collide.
- `constraint_set_hash` records WHAT was enforced: every enabled rule's version
  and its per-slot enforcement, the usable forbidden set, the GC band, and the
  vendor selection. It is how a later run can tell whether it was designed under
  the same rules, independently of which sequence came out.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence

from bt5.core.services import FoldEngine
from bt5.core.spec import Spec
from bt5.core.types import Provenance
from bt5.rules.vendors import VendorSelection


def engine_versions(fold: FoldEngine | None) -> Mapping[str, str]:
    """The engine versions to record in provenance. Empty when no engine loaded."""
    if fold is None:
        return {}
    return {fold.name: fold.version}


#: The app version stamped into every design. A single source; the packaging lane
#: owns the real version string, and this is the value the engine records until
#: that is wired through.
APP_VERSION = "0.1.0"

#: Same width as `score.hashing.HASH_LENGTH`, deliberately: the constraint hash
#: sits beside the design hash and a reader should not have to wonder why one is
#: longer than the other.
_HASH_WIDTH = 12


def _slot_enforcement(rule: Spec, ctx: object) -> dict[str, str]:
    """Per active slot that the rule gates, the enforcement it resolves to.

    This is the thing `presets.resolve` cannot see -- it reads the static
    `enforcement` ClassVar, but a rule like d4 is SOFT there and HARD_REPAIR in a
    lentiviral slot -- so the constraint hash records the resolved value.
    """
    out: dict[str, str] = {}
    for slot in ctx.active_slots:  # type: ignore[attr-defined]
        if rule.gate(slot):
            out[slot.role] = rule.enforcement_for(slot).value
    return out


def constraint_set_hash(
    rules: Sequence[Spec],
    ctx: object,
    *,
    forbidden: Sequence[str],
    gc_bounds: tuple[float, float] | None,
    vendors: VendorSelection,
) -> str:
    """A stable digest of what this design enforces, over canonical JSON."""
    payload = {
        "rules": {
            rule.id: {
                "version": rule.version,
                "enforcement_by_slot": _slot_enforcement(rule, ctx),
            }
            for rule in rules
        },
        "forbidden": sorted(forbidden),
        "gc_bounds": list(gc_bounds) if gc_bounds is not None else None,
        "vendors": list(vendors.keys),
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:_HASH_WIDTH]


def design_hash_context(
    *,
    backbone_sequence: str,
    table_id: int,
    forbidden: Sequence[str],
    gc_bounds: tuple[float, float] | None,
    vendors: VendorSelection,
) -> Mapping[str, object]:
    """The context `design_hash(cds, context=...)` is salted with.

    Includes the BACKBONE sequence, without which two vectors carrying the same
    insert would produce the same hash -- the collision this argument exists to
    prevent, since `design_hash` otherwise sees only the CDS.
    """
    return {
        "backbone": backbone_sequence.upper(),
        "table_id": table_id,
        "forbidden": sorted(forbidden),
        "gc_bounds": list(gc_bounds) if gc_bounds is not None else None,
        "vendors": list(vendors.keys),
    }


def build_provenance(
    *,
    seed: int,
    table_id: int,
    fold: FoldEngine | None,
    constraint_hash: str,
    degradations: Sequence[str],
) -> Provenance:
    """The first `Provenance` constructed in `src/`."""
    return Provenance(
        app_version=APP_VERSION,
        seed=seed,
        engine_versions=engine_versions(fold),
        codon_table_name=f"ncbi_{table_id}",
        constraint_set_hash=constraint_hash,
        degradations=tuple(degradations),
    )
