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

## Reviewing prompt/doc diffs (docs/buildout/*.md session prompts)
These carry dense `file.py:123` / commit-SHA / issue-number citations meant
to be followed literally by unattended sessions. This repo's docs culture is
precise — on two separate buildout reviews (2026-09-01 initial, then a
follow-up after fixes) essentially every `file:line` citation, quoted
docstring/string, exported-symbol list, and skill model/effort front-matter
checked out exactly against source. Spend the verification budget instead on:

- **Diff "never touch" lists against each other and against `ls
  packages/engine/src/bt5/`** (currently 9 subdirs: cassette, codon, core,
  design, rules, score, solver, structure, vector), not just against the
  README's ownership table. A prior review caught one lane's file omitting
  two of nine from its enumeration that every sibling included — confirmed
  fixed in the 2026-09-01 buildout PR (S3 now lists `codon/`/`structure/`).
  Keep re-checking this on every future edit to these prompts; it is a
  copy-paste-prone spot.
- **Shared fixture/support files one level up from a glob-partitioned test
  directory** — e.g. `packages/engine/tests/rules/conftest.py`, shared by two
  rule-catalog sessions splitting `rules/catalog/{b,c}*.py` vs `{d,e,f}*.py`.
  Confirmed fixed 2026-09-01: both prompts and the README now declare it
  read-only for both sessions, with an explicit fallback (define locally, or
  open an issue for a truly-shared helper). Check every glob-partitioned
  split for this class of gap: a support file that matches no session's
  glob but that two sessions would plausibly both want to edit.
- **A decision doc's "Rejected" section is a commitment, not just a note —
  diff it against the actual prompt files.** The 2026-09-01
  `docs/decisions/2026-09-01-parallel-buildout-sessions.md` states "S2 and S4
  post a design note as their first PR comment instead [of plan mode]," but
  only S2's prompt (`s2-biosecurity-screen.md`) actually contains that
  instruction — S4's (`s4-rules-liabilities.md`, the one doing the
  cross-lane D3 splice-repair problem at xhigh effort, unattended) did not
  tell the session to do it. Confirmed fixed at a6c2250:
  `s4-rules-liabilities.md` now carries an equivalent gate before D3, and the
  README's "Not plan mode" paragraph names both sessions. Grep every session
  file for a phrase the decision doc attributes to it by name; do not assume a
  plural attribution ("S2 and S4 do X") was actually applied to both.
- **Commit SHAs describing "current" main state go stale within one PR's
  lifetime.** Check the cited SHA against `git log origin/main` directly. A
  *living* doc (a reusable prompt/README) is fine if it explicitly says the
  SHA is point-in-time and tells the reader to fetch fresh instead
  (this repo's buildout docs do this correctly); flag it only if a living doc
  states a stale fact as current without that qualifier. A `docs/decisions/`
  ADR-style file is point-in-time by convention — don't flag staleness there.
- **No GitHub API access in this sandbox.** Issue/PR number claims can only
  be cross-checked against local squash-merge commit messages (`git log
  --oneline`, real `(#N)` suffixes), not issue titles/bodies/state. Say so
  explicitly rather than silently skipping those citations — e.g. issues
  #45, #56, #69, #70, #74, #78, #80, #82, #83 in the 2026-09-01 buildout docs
  have no corresponding merged-commit trail to check against.
