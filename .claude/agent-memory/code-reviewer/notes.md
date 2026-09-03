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
`bt5/cli.py`'s `_REPORTABLE_ERRORS` misses raise sites reachable from `bt5
design` (e.g. `FileTableProvider.genetic_code()`'s bare `ValueError`,
Biopython `SeqIO.read()` on a malformed `--backbone` file) — both escape as
raw tracebacks, though exit code stays nonzero. Non-blocking UX gap, not a
CLAUDE.md violation. Check the tuple against what callees one level of
transitive dependency down actually raise, not just the caller's docstring.

## Reviewing prompt/doc diffs (docs/buildout/*.md session prompts)
This repo's docs culture is precise — file:line citations, quoted strings,
exported-symbol lists and SHAs have consistently checked out exactly against
source. Recurring gap classes worth the verification budget:
- diff each prompt's "never touch" list against the others and against
  `ls packages/engine/src/bt5/` (9 subdirs), not just the README table;
- check a shared fixture one level up from a glob-partitioned test split is
  declared read-only by every session whose glob doesn't cover it;
- diff a decision doc's plural-attribution sentences ("S2 and S4 do X")
  against each named prompt individually — a plural claim applying to only
  one file is the recurring failure shape;
- cross-check issue/PR numbers only against `git log --oneline` squash-merge
  `(#N)` trails (no GitHub API here); say so when a citation has no trail.

## "read-only agent" is policy here, not mechanism — do not claim otherwise
`rule-auditor` holds `Agent` (to resolve `brief_ref` via `docs-miner`), which
also reaches `batch-editor` (holds `Edit`). Do NOT write this up as "the only
agent whose read-only-ness is not mechanical" — `code-reviewer`, `debugger`,
`docs-miner`, `gate-runner`, `security-reviewer` all hold `Bash`, which writes
via `sed -i` too. The `Agent` grant widens an existing surface, it doesn't
create a new class. `tools: Agent(docs-miner)` (parameterised) is dropped
silently by this CLI, leaving no `Agent` tool at all (fails closed, doesn't
scope). An agent holding `Bash` can reach the `claude` CLI on PATH regardless
of its tool list. Before calling a control "mechanical", name the mechanism
and check it holds. When a diff corrects a fact, grep the whole file for the
old value — a stale citation can survive in an example block below the fix.

## Rule-count/repair-count prose in `.claude/*.md` is a recurring drift point — verify, don't assume loss of a number is a defect
This repo has repeatedly patched `.claude/agents/`, `.claude/rules/`, `.claude/skills/`
prose that hardcoded "25 catalog rules" or "22 of 25 declare SINGLE_PASS" style counts
(574ea0e stopgap 15->25; #103 replaced a hardcoded count with an enumeration; a later
diff removed the "25" count and the "d3_splicing/b9_out_of_frame_atg/f5_at_window"
FIXED_POINT enumeration entirely, replacing with structural explanations). When review-
ing such a diff: (1) grep `git log --all -p -- <path>` for any historical number the
prose cites as "how long this bug existed" before calling it fabricated — the true prior
state can be one commit further back than a three-dot diff shows; (2) verify every
line:number and mechanism citation in the *replacement* text against live source (e.g.
`core/spec.py:231`, `solver/repair.py:417`, `solver/catalog.py`'s `policies()`) rather
than trusting the count it removed was accurate — sometimes the old enumeration was
itself already stale (e.g. a "four brief.md sections" list that omitted `2.C` even
though `c1_cai`/`c3_min_max` already cited it); (3) losing a concrete example (which
rules currently declare `FIXED_POINT`) in exchange for a rule that never goes stale is a
legitimate trade, not a diligence failure, as long as a canonical example survives
elsewhere (e.g. "splice donors are the canonical case" stays in `rule-add/SKILL.md`).
`.claude/agents/`, `.claude/rules/`, `.claude/skills/` are not in CLAUDE.md's lane table
or its §2 protected-path list — editing their prose needs no lane issue and no `approved:*` label.
