# Code-reviewer notes for BT5

## Enforcement is per-slot, not per-class (score/presets.py)
`Spec.enforcement` (ClassVar) is only the FLOOR. The real routing is
`enforcement_for(slot) -> Enforcement`, which can escalate SOFT to HARD_REPAIR/
HARD_CHECK depending on `slot.modality` (e.g. `d4_internal_polya`: SOFT floor,
HARD_REPAIR on LENTIVIRAL/AAV/GENOME_INTEGRATED). Any guard in `score/` that
enforces CLAUDE.md §3.5 ("hard constraints never carry a weight") must ask
`enforcement_for` per slot admitted by the preset's `modality`, not read the
ClassVar alone — reading the ClassVar alone is exactly the bug fixed in #72
(LENTIVIRAL/AAV both shipped `WeightEntry("2.D4", ...)` past a ClassVar-only
guard). Watch for this pattern recurring anywhere a preset/objective assembly
step decides "is this rule scorable" — check it asks per-slot.

## `_slots_admitted_by`-style enumeration is deliberately over-broad
`packages/engine/src/bt5/score/presets.py::_slots_admitted_by` builds every
`(host, table_id) x SlotRole` pair for a fixed modality (from
`LOCKED_TRANSLATION_TABLE x SlotRole`), not just the "real" host for that
modality. This is intentional (documented in-file): a guard should err toward
refusing a weight rather than guessing a representative host. Confirmed
correct as of 2026-09 because every catalog rule's `gate`/`enforcement_for`
keys only on `slot.modality`, never `slot.host`/`slot.role` — re-check this
assumption if a future rule starts gating on host/role, since over-broad
enumeration could then cause spurious refusals.

## `spec()` instantiation inside a guard (try/except TypeError) is a known,
accepted gap, not a new one
`_unscored_enforcement` tries `spec()` to call `enforcement_for` per slot; if
the spec needs required constructor args it falls back to the ClassVar floor
only. As of 2026-09 every catalog rule's `__init__` has all-default args, so
this path is dead code today — but it is a silent hole if a future rule adds a
required constructor arg while also relying on `enforcement_for` to escalate.
Also getattr-based fallback for specs missing `gate`/`enforcement_for`
entirely (shouldn't happen for real catalog rules, which implement the `Spec`
Protocol, but silently permissive for anything that doesn't). Worth a
non-blocking mention if it recurs, not blocking on its own.

## Test inversions aren't automatically suppression
Flipping an assertion from "preset DOES weight X" to "preset does NOT weight
X" is legitimate when the old assertion was pinning the *bug* itself (i.e. the
production code changed to fix a real defect, and the test is updated to
match correct behavior) — especially when accompanied by a new invariant test
class that pins the general property going forward. Distinguish this from
suppression: suppression loosens a check to dodge a still-present failure;
here the underlying defect was actually fixed in prod code first.

## docs/decisions.md scientific-impact flagging
When a lane PR changes ranking/weight behavior in a shipped preset, a good
diff says so explicitly in `docs/decisions.md` under "Scientific impact" and
notes it is NOT "none", triggering §7b's owner-merge requirement rather than
self-merge. This is the correct pattern to look for — treat its absence on a
weight/ranking-affecting diff as a finding.
