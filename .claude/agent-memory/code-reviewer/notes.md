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

## Reviewing prompt/doc diffs (e.g. docs/buildout/*.md session prompts)
These carry dense `file.py:123` / commit-SHA / issue-number citations meant to
be followed literally by unattended sessions. On the 2026-09-01 buildout
review almost every citation checked out exactly (line ranges, exported
symbol lists, quoted docstrings, file sizes in KB, gate numbers from PLAN.md,
skill model/effort front-matter) — this repo's docs culture is precise, so
spend the verification budget on things that are unusually cheap to get wrong
and expensive to leave wrong:
- **Commit SHAs describing "current" state of `main` go stale within the same
  PR's lifetime** (main moved twice under one buildout PR in a few hours).
  Check the cited SHA against `git log origin/main` directly; don't assume a
  SHA written earlier in the same PR is still accurate later in review.
  Prefer/require "fetch fresh and branch" commands over a hardcoded SHA in any
  prompt meant to be reused.
- **Decision files under `docs/decisions/` are point-in-time ADR-style
  records** ("Evidence:" reflects state *at the time of the decision*) — do
  not flag an unqualified stale fact inside an already-merged decision file
  as a new defect; only flag it if a *living* doc (a prompt meant to be
  re-run, a README) states the same fact without a "checked as of / verify
  before trusting" qualifier.
- **Path-glob ownership splits ("session A owns `b*.py`, session B owns
  `d*.py`") reliably miss shared fixture files** that don't match anyone's
  glob, e.g. `packages/engine/tests/rules/conftest.py` here — neither S3 nor
  S4 explicitly owned it, and it's a plausible (if not certain) two-writer
  collision the mutex list didn't cover. Check for this class of gap
  specifically: shared `conftest.py`/`__init__.py`-adjacent support files one
  level up from a glob-partitioned test directory.
- **No GitHub API access in this sandbox** (github MCP tools are not actually
  wired even though instructions describe them) — issue/PR number claims can
  only be cross-checked against local squash-merge commit messages
  (`git log --oneline`, which carries the real `(#N)` suffix), not against
  issue titles/bodies/state. Say so explicitly rather than silently skipping
  those citations.
- **"Never touch" lists should be diffed against each other**, not just
  against the ownership table — in this review one lane's file omitted two
  of nine engine subdirectories from its explicit never-touch enumeration
  that every sibling lane's file included, a copy-paste-style asymmetry that
  is easy to miss reading files one at a time but obvious diffing them side
  by side.
