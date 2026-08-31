# Model + effort routing — setup notes

Configuration for AI coding sessions in this repo: which work runs on which model at
which effort, and why. Verified against Claude Code **2.1.251**.

No application, test or CI source was changed.

## The rule

**Capability failure raises the model. Diligence failure raises the effort.** Work that
fails because the model could not figure it out gets a stronger model. Work that fails
because files were skipped, a change went untested, or the agent bailed partway gets
higher effort. Conflating the two is how you end up paying opus to be thorough.

Volume work — code search, doc extraction, running gates, mechanical batch edits — runs
on a cheap model inside a subagent, so its output never enters the main window. That is
where the saving is.

## Routing table

| Agent | Model | Effort | Trigger | Failure mode this answers |
|---|---|---|---|---|
| `Explore` *(overrides built-in)* | haiku | low | Locating files, symbols, call sites | Neither — pure volume |
| `docs-miner` | sonnet | medium | Any answer living in `docs/` prose | CAPABILITY — resolving a `brief_ref` and spotting a superseded row |
| `gate-runner` | haiku | high | The full gate chain; why CI is red | DILIGENCE — running four of five gates and reporting green |
| `batch-editor` | sonnet | medium | One identical edit across ≥5 files | Volume, with a capability floor for staying ruff/mypy-clean |
| `code-reviewer` | sonnet | high | The branch diff before a PR | DILIGENCE — the rules are written down; missing one is inattention |
| `rule-auditor` | opus | xhigh | A rule Spec's provenance changed | CAPABILITY — does the cited source support the number? |
| `debugger` | opus | xhigh | A failure a first-pass fix already missed | CAPABILITY — a wrong model of the failure |
| `security-reviewer` | opus | xhigh | `/pre-pr`, on the biosecurity/CI surface | CAPABILITY — intent, not signature |
| `/architect` *(skill)* | opus | **ultracode** | Cross-lane, `core/`, MINOR-vs-MAJOR | CAPABILITY — fields and protocol methods classify oppositely |
| `/escalate` *(skill)* | opus | **ultracode** | A decision a normal attempt got wrong | CAPABILITY |

`code-reviewer` and `rule-auditor` carry `memory: project`
(`.claude/agent-memory/<name>/`) — committed and reviewable on purpose.

### What "ultracode" is

Read out of the 2.1.251 bundle: *"Current effort level: ultracode (xhigh + dynamic
workflow orchestration; this session only)"* and *"Enable ultracode for the session:
xhigh effort plus standing dynamic-workflow orchestration."*

So it is **`xhigh` effort plus a standing mandate to orchestrate `Workflow`s** — built on
`xhigh`, not `max`. It is session-scoped, selected by `/effort ultracode` or the boolean
`ultracode` settings key.

**Frontmatter cannot name it.** The agent schema validates `effort` as
`low|medium|high|xhigh|max` or an integer. So `/escalate` and `/architect` assemble it
from its two halves: `effort: xhigh` in frontmatter, and the orchestration mandate in the
body. Do not write `effort: ultracode` anywhere — it is not a valid value.

**Nothing uses `max`.** Opus 5's cost index is `{low 0.67, medium 0.76, high 1,
xhigh 1.6, max 1.7}`; `max` is ~6% above `xhigh`, and ultracode tops out at `xhigh`
because the extra capability comes from fanning out, not from thinking longer alone.

`ultracodeKeywordTrigger` is on by default: **typing "ultracode" in any prompt opts that
turn into the `Workflow` tool**, whether or not a skill was invoked.

## What was added

```
CLAUDE.md                     restructured, 177 -> 197 lines (statements kept, rationale moved)
.claude/settings.json         EXTENDED — no existing entry removed
.claude/agents/*.md           8 agents
.claude/skills/*/SKILL.md     8 skills
.claude/rules/*.md            4 path-scoped rules files
.claude/hooks/                6 scripts (5 hooks + statusline)
scripts/gates.sh              the gate chain; the only non-.claude file added
docs/decisions.md             durable session decisions
.gitignore                    local Claude state
SETUP-NOTES.md                this file
```

### What was kept from the existing config, and why

`.claude/settings.json` was deliberately authored and is **extended, not replaced**. All
55 `permissions.allow` entries and all 6 `permissions.deny` entries survive (verified by
diffing the key sets). Root `CLAUDE.md` keeps every rule's **statement**; only rationale
moved, because root CLAUDE.md is the only text re-read from disk after compaction.

### What was changed

- **`Write(.claude/**)` → `Edit(.claude/**)`.** File-path rules are consulted only for
  `Edit(path)` and `Read(path)`; a `Write(path)` rule is accepted, never consulted, and
  warns at startup.
- **Added the command forms that actually work:** `Bash(.venv/bin/pytest:*)` and friends.
  The existing `Bash(pytest:*)` rules are prefix matches against `/root/.local/bin`,
  whose interpreter has no numpy — they pre-approved the broken form.
- **Lane table globs corrected** to `packages/engine/src/bt5/…`; `git ls-files` returned
  0 for the previous `bt5/solver/`-style paths.
- **M7/M9/M10 marked *(planned, no files)*** — those directories do not exist.
- **`env.PYTHONHASHSEED=0`** — CI sets it, local did not, so an ordering bug reproduced
  on exactly one side.

## How to verify

```bash
bash .claude/verify-setup.sh
```

Checks frontmatter parses with both `name` and `description` (missing either silently
skips the file), the description token budget, every `paths:` glob matching real files,
every referenced command existing in a manifest, hook executability, and the hook probes
below.

### First session after merging this

**Agent and skill definitions load at session start.** They were created mid-session
here, so none of them were active while this was built: `Agent(subagent_type="docs-miner")`
returned *"Agent type 'docs-miner' not found"*, and a smoke-invoked `Explore` ran as the
**built-in** agent on Opus — which is precisely the cost the override exists to remove.

So one runtime check is still outstanding and must be done in a **fresh session**:

1. Run `/agents` (or `Agent(subagent_type="docs-miner")`) and confirm all eight resolve.
2. Smoke-invoke each and confirm the transcript shows the intended model — in particular
   that `Explore` now runs on haiku rather than inheriting the session model.
3. **If `effort` is rejected on a haiku agent**, drop the field from `Explore.md` and
   `gate-runner.md`; the model is already the cheap tier, and `scripts/gates.sh` is what
   actually stops `gate-runner` skipping a step. Everything else is unaffected.

Everything statically checkable — frontmatter parsing, `name`/`description` presence,
model and effort values, tool names, glob coverage, hook behaviour — is verified by
`.claude/verify-setup.sh` and passes.

Hook probes worth re-running by hand after any edit:

```bash
# must REWRITE
echo '{"tool_name":"Bash","tool_input":{"command":"pytest tests/invariants"}}' \
  | python3 .claude/hooks/compact_output.py
# must be SILENT (compound, explicit flag, or a full-output target)
echo '{"tool_name":"Bash","tool_input":{"command":"pytest -q && mypy"}}' \
  | python3 .claude/hooks/compact_output.py
# must ASK — this is the one whose failure mode is matching nothing
printf '{"tool_name":"Edit","tool_input":{"file_path":"%s/packages/engine/src/bt5/core/types.py"},"cwd":"%s"}' "$PWD" "$PWD" \
  | python3 .claude/hooks/protect_paths.py
```

## Tuning

Start conservative; promote only the agents that underperform. Signals:

1. **`code-reviewer`'s `memory: project`.** Switch to `local` the first time
   `git status --porcelain | grep agent-memory` shows a memory file inside a lane PR's
   diff, or the first conflict in it across two parallel branches — CLAUDE.md §5 names
   cross-PR conflicts as the most expensive failure here.
2. **`gate-runner` at haiku/high.** The cheapest bet, resting on the weakest assumption:
   that `effort` is honored on haiku. `scripts/gates.sh` is the insurance, but the agent
   still has to *report* correctly. Move to sonnet/medium after **one** instance of it
   reporting a gate green when `gates.sh` printed a non-zero line, or reporting exit 4 or
   10 as a test failure. A runner you have to double-check has negative value.
3. **The PostToolUse `ruff format` hook.** Drop it after three or more
   stale-`old_string` Edit failures in a week landing right after a successful Edit to
   the same file — reformatting moves lines the model is about to match. Replacement is
   free: one `.venv/bin/ruff format` inside `/pre-pr` reproduces the CI gate exactly, at
   one invocation per PR instead of one per edit.

## Caveats that will bite

- **Path-scoped rules drop out after compaction** until a matching file is read again.
  That is why every rule *statement* is in root CLAUDE.md and only rationale is in
  `.claude/rules/`. Do not move a statement out.
- **`/escalate` and `/architect` are turn-scoped.** The model and effort override applies
  to the rest of the current turn and is not saved; the next prompt reverts. For
  multi-turn work, re-invoke each turn.
- **Overriding `Explore` changes what loads.** Built-in `Explore` and `Plan` skip the
  CLAUDE.md hierarchy; a project agent named `Explore` may not. That is part of why
  CLAUDE.md is kept under 200 lines and why no `.claude/rules/` file is un-`paths:`ed.
- **Cloud and web sessions may ignore these settings.** `CLAUDE_EFFORT` and
  `CLAUDE_CODE_PROVIDER_MANAGED_BY_HOST=1` take precedence over `model`/`effortLevel`.
  These settings govern local CLI sessions.
- **`autoCompactWindow` is deliberately left at `auto`.** It takes a token count, and the
  CLI says `auto` is "tuned for your model and strongly recommended", while overriding
  "may result in high token usage". To pin it anyway:
  `"autoCompactWindow": "500k"` in `.claude/settings.json`.
- **`.claude/agent-memory/` is not gitignored.** `memory: project` is meant to be
  committed and reviewed; `/pre-pr` greps for it so it never lands unread. Trim it
  periodically — it is an instruction channel.

## Findings — reported, not fixed

All verified in this checkout. Each is application, test or CI source, so all are out of
scope for a configuration change. The first four are the ones with teeth.

1. **The biosecurity gate guards the wrong file.**
   `test_kmer_index_takes_no_external_database` reads only `core/services.py` and regexes
   the frozen `KmerIndex` Protocol. The only implementation is `ConstructKmerIndex.of` at
   `vector/kmers.py:158`. A `database=` parameter added there leaves the Protocol
   untouched, still conforms (a widened signature with a default does), needs no approval
   label (`check-approval-labels.sh` has no rule matching
   `packages/engine/src/bt5/vector/`), and evades `kmers.py:461`'s conformance assertion,
   which is under `if TYPE_CHECKING` — mypy-only, and mypy runs in no CI job.
2. **`HYPOTHESIS_PROFILE` is inert.** `tests/conftest.py` registers `ci`/`dev`/`nightly`
   then calls `settings.load_profile("dev")` unconditionally and never reads the variable.
   CI's `HYPOTHESIS_PROFILE: ci` on the `invariants` job has no effect: **every run, local
   and CI, is 50 examples, not 200.**
3. **Every `brief_ref` is unresolvable by literal search.** All 15 catalog values
   (`2.B1`, `2.E4`, …) return zero `grep -F` hits in `brief.md`. They are
   section-qualified: `### 2.E` plus a row `E4`, in one of two anchor shapes. The
   procedure is in `.claude/rules/rules-catalog.md`.
4. **`mypy` is in no CI job**, despite `strict = true` and the PR checklist asking for it.
5. **`test_d1_restriction_sites.py` does not exist** — the only catalog rule without a
   paired test, and `d1` calls itself the shape to copy.
6. **`brief.md:141` (E4) is struck through** and "corrected 2026-08-28"; a rule encoding
   the superseded threshold would pass all 11 contract assertions.
7. **`-p no:randomly` is a silent no-op** (`pytest-randomly` is not installed);
   **`-m "not slow"` deselects nothing** (the marker is applied zero times);
   **`fail_under = 85` has never been evaluated** (no `--cov` job);
   **`goldens-not-hand-edited` does not exist**.
8. **`verify.py:308` and `vector/survey.py:150` default `table_id: int = 1`**, against
   CLAUDE.md §3.1's "explicit and never defaulted". `solver/pipeline.py:52` does it right.
9. **The `changes` job's outputs are consumed by no job**, yet it is in
   `required-checks.needs`, where `skipped` counts as failure.
10. **The RNG grep is narrower than §3.7 implies** — it covers
    `seed|rand|randn|choice|randint` under `packages/engine/src/` only, so
    `np.random.random(`, `.shuffle(`, `.permutation(`, `.normal(` and `.uniform(` pass it.
