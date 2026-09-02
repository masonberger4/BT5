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
`(host, table_id) x SlotRole` pair for a fixed modality, not just the "real"
host for that modality. Intentional (documented in-file): a guard should err
toward refusing a weight rather than guessing a representative host.
Re-check this if a future rule starts gating on `slot.host`/`slot.role`
instead of only `slot.modality`, since over-broad enumeration could then
cause spurious refusals.

## `spec()` instantiation inside a guard (try/except TypeError)
`_unscored_enforcement` tries `spec()` to call `enforcement_for` per slot; if
the spec needs required constructor args it falls back to the ClassVar floor
only. Dead code today (every catalog rule's `__init__` has all-default args)
but a silent hole if a future rule adds a required constructor arg while
relying on `enforcement_for` to escalate. Worth a non-blocking mention if it
recurs, not blocking on its own.

## Test inversions aren't automatically suppression
Flipping an assertion from "preset DOES weight X" to "does NOT weight X" is
legitimate when the old assertion was pinning the bug itself and prod code
changed first to fix it, especially with a new invariant test class pinning
the property going forward. Suppression loosens a check to dodge a
still-present failure; distinguish the two by whether prod code actually
changed.

## docs/decisions/ scientific-impact flagging
A lane PR that changes ranking/weight behavior in a shipped preset should say
so explicitly under "Scientific impact" (not "none"), triggering §7b's
owner-merge requirement. Treat its absence on a weight/ranking-affecting diff
as a finding.

## A thin adapter's "reportable errors" tuple is an allowlist, not a catch-all
`bt5/cli.py`'s `_REPORTABLE_ERRORS` (DesignError, InfeasibleConstraints,
VerificationError, VectorError, OrderError, OSError) misses at least two real
raise sites reachable from `bt5 design`: `FileTableProvider.genetic_code()`
raises bare `ValueError` for an unknown `--table-id`, and Biopython's
`SeqIO.read()` raises its own (non-`VectorError`) exceptions on a malformed
GenBank file passed to `--backbone`. Both escape as raw tracebacks instead of
a clean `bt5: ...` stderr line, though the process still exits nonzero so
scripted callers aren't broken — non-blocking UX gap, not a CLAUDE.md rule
violation, but worth naming every time a new CLI/adapter layer wraps engine
calls in a curated exception tuple: check the tuple against what the actual
callees (including one level of transitive dependency, e.g. the table
provider and the GenBank parser, not just `design()`'s own docstring) can
raise, not just what `design()`'s docstring documents.

## Reviewing prompt/doc diffs (docs/buildout/*.md session prompts)
This repo's docs culture is precise — file:line citations, quoted strings,
exported-symbol lists and SHAs have consistently checked out exactly against
source across multiple buildout reviews. Spend verification budget on the
gap classes that actually recur, not on re-confirming precision:
- diff each prompt's "never touch" list against the others and against
  `ls packages/engine/src/bt5/` (9 subdirs: cassette, codon, core, design,
  rules, score, solver, structure, vector), not just against the README table;
- check a shared fixture/support file one level up from a glob-partitioned
  test split (e.g. a lane-shared `conftest.py`) is declared read-only by
  every session whose glob doesn't cover it, not just some;
- diff a decision doc's "Rejected"/plural-attribution sentences ("S2 and S4
  do X") against each named prompt individually — a plural claim applied to
  only one file is the recurring failure shape;
- treat a commit-SHA-pinned "current main" claim as stale unless the doc says
  the SHA is point-in-time (living docs should; `docs/decisions/` ADRs are
  point-in-time by convention, don't flag those);
- no GitHub API access here — cross-check issue/PR numbers only against
  `git log --oneline` squash-merge `(#N)` trails, and say explicitly when a
  citation has no such trail rather than silently skipping it.
