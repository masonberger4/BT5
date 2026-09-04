# Wave 2 — the scorecard the engine already knows how to build, and five rules it doesn't

Six prompts. **W0 runs alone first.** Then W1, W3, W4 and W5 run at the same time in four
separate Claude Code sessions; W2 cuts its branch only after W1 merges.

Each file in this directory is a complete prompt: open it, copy the whole thing, paste it
into a fresh session on this repo. This index holds what the six have to agree on. Read it
once before launching any of them; the session files assume you have.

## Why these six

`main` is green. Wave 1 (`docs/buildout/`) landed everything except S2, which the **owner
killed on purpose**: PR #87 closed unmerged because DNA-synthesis vendors screen orders
before synthesis, and that is where the risk gate belongs.

What is missing is again not machinery. It is the same shape of gap PR #89 closed for the
gallery: **`bt5/score/` exports capability that `bt5/design/` never calls.**

| What exists | Where | Who calls it |
|---|---|---|
| `presets.resolve()` → per-modality weights with provenance | `score/presets.py:206` | **nobody** — `ranking.py:317` uses `spec.default_weight` |
| `detect_conflicts` / `hard_versus_soft` | `score/conflicts.py:144,217` | **nobody** — `runner.py:630` hard-codes `conflicts=()` |
| host-frequency null | `ranking.py:168` | reachable, but `_host_usage` looks up the wrong filename (#98) |
| `Preset.strain_protocol` → `QcReport.strain_protocol`, rendered | `presets.py:330`, `report.py:298` | **nobody** — `runner.py:632` omits the argument |
| `RuleSet.gated_out` | `solver/catalog.py:130` | recorded, **never rendered** |
| `DesignContext.weights` (a frozen `core/` field) | `core/context.py:128` | **nobody reads it** |

Meanwhile the report tells every user, forever, that "protein-level biosecurity screening
did not run" — a degradation for a screen the project decided not to build, pinning
`QcReport.is_complete` to `False` by construction. That is the exact defect PR #89's own
decision record says it existed to remove. And `main` carries **no record at all** of the
decision to drop the screen; that record lived on PR #87's branch and died with it.

The rule catalog is also **half built**: `brief.md` §2 defines 50 numbered rules
(A1–A4 structural, B1–B11, C1–C10, D1–D8, E1–E12, F1–F9). Twenty-five have files.

Wave 2 closes the design path and adds the five highest-evidence missing rules.

**No new lane.** M9 `packages/server/` and M10 `apps/web/` stay out. A UI rendering a
scorecard with an empty conflict panel, unweighted objectives and a permanent false
degradation would be a UI built on a lie — which is `score/gallery.py`'s own argument
about G4.

## Ordering

```
  PRE-WAVE (alone, nothing else branched)
    W0  encoding sweep + the record that never landed
         |
         +------------+------------+------------+
  WAVE   W1 host      W3 D5        W4 B6/B7     W5 C7/C8
         reference    cryptic      initiation   composition
         sets (#98)   transcription
         |
         W2 design path  <- branches only after W1 MERGES
```

**Why the design lane is serialized.** W1 (#98) and W2 (preset weights) *each change the
winner*. #98 flips the null from uniform-synonymous to host-frequency, moving every
percentile; preset weights change the rank key directly. Landed together, a moved winner
cannot be attributed to either — and that attribution is exactly what the PR template's
"scientific impact" section and §7b's owner-merge gate exist to read.

W2 cuts from a fresh `main` after W1 lands. Cost: one PR cycle. What it buys: W2's diff
shows only the preset delta against an already-rebaselined `main`.

The three rules sessions have no such coupling and start immediately.

## Ownership matrix

Every writable path belongs to exactly one session. If your work needs a path you do not
own, **stop and open an issue** — do not reach into another lane.

| # | Session | Branch | Owns (write) | Label |
|---|---|---|---|---|
| W0 | [Encoding + records](w0-encoding-and-records.md) | `claude/w0-encoding-and-records` | **every text-I/O call site repo-wide**; `docs/decisions/2026-09-03-biosecurity-screen-dropped.md`; the guard step in `ci.yml`'s `python-quality`; `[tool.ruff.lint]` | `approved:contract-change`, `approved:oracle-change`, `approved:ci-change` |
| W1 | [Host reference sets](w1-host-reference-sets.md) | `claude/w1-host-reference-sets` | `bt5/codon/tables.py`, `bt5/rules/catalog/c1_cai.py`, `tests/codon/**`, `tests/rules/test_c1_cai.py`, **`design/runner.py:_host_usage` only** | none |
| W2 | [Design path truth](w2-design-path-truth.md) | `claude/w2-design-path-truth` | `bt5/design/**`, `bt5/score/**`, `tests/design/**`, `tests/score/**`, `tests/data_integrity/**` | `approved:oracle-change` |
| W3 | [D5 cryptic transcription](w3-d5-cryptic-transcription.md) | `claude/w3-d5-cryptic-transcription` | `rules/catalog/d5_*.py` + paired test; `tests/rules/test_d1_restriction_sites.py` (new) | none |
| W4 | [B6/B7 initiation](w4-b-initiation.md) | `claude/w4-b-initiation` | `rules/catalog/b6_*.py`, `b7_*.py` + paired tests | none |
| W5 | [C7/C8 composition](w5-c-composition.md) | `claude/w5-c-composition` | `rules/catalog/c7_*.py`, `c8_*.py` + paired tests | none |

### Two named exceptions

1. **W1 writes three files outside its lane** — `codon/tables.py` (M5), `c1_cai.py` (M4),
   and exactly one function inside W2's (`design/runner.py:_host_usage`). No other session
   touches any of the three. This is the device Wave 1 used for the shared rules conftest.
2. **`packages/engine/tests/rules/conftest.py` is read-only for W3, W4 and W5**, exactly as
   in Wave 1. It holds the helpers every rule test imports (`construct()`,
   `wrapping_construct()`, `slot()`, `context()`, the `services` fixture) and three rules
   sessions run concurrently. A helper only one session needs goes in that session's own
   test file; one that genuinely belongs to several becomes an issue, landed after they
   merge.

Rule *source* files collide with nothing — `core/registry.py` autodiscovers — but their
tests share a conftest, and that asymmetry is easy to miss.

### Branches are not pre-created, on purpose

Each session cuts its own branch from a freshly fetched `main`, with the command in its
prompt. Nobody makes the six up front.

A branch created before its session starts is born pointing at whatever `main` was that
day, and `main` moves. Wave 1's own first PR (#86) is the worked example: cut from
`628e130`, then #84 and #85 merged underneath it, costing a modify/delete conflict and
invalidating guidance in all six prompts. `git checkout -B` from a just-fetched
`origin/main` costs one command and cannot go stale.

**This matters more in Wave 2 than it did in Wave 1**, because W2 is *defined* as branching
after W1 merges.

### The global mutexes

1. **`pyproject.toml` — W0 only**, and only one line in `[tool.ruff.lint]`. Nobody else
   needs it: no session adds a dependency, and no session adds a package.
2. **`data/` — nobody.** W1 adds no file; the four reference sets already ship. W5's C7
   *wants* a CSC table and must ship reporting `unavailable` instead — see its prompt.
3. **`core/` — nobody.** Not one gap in this wave needs a field. That is a finding, not an
   omission: `DesignContext.weights`, `Conflict.relaxations`, `BiosecurityVerdict.status`
   and `ObjectiveScore.null_kind` were all built for exactly this and sit unused. If you
   believe you need a `core/` change, that is `/architect`, not a regeneration.
4. **`.github/` — W0 only**, one step inside an existing job.

## The un-draft queue

CI capacity is the binding constraint: 20 concurrent job slots, ~12 per Python PR, so
**at most 5 open non-draft PRs**. Drafts skip the expensive jobs and are free.

**Corrected arithmetic for this wave.** `ci.yml:11` triggers on `labeled`/`unlabeled` —
deliberately, so the approvals gate can clear when a label is applied. That means **each
label application re-runs the whole workflow**. W0 carries three labels and can burn ~48
slots over its life. That is the concrete reason W0 runs pre-wave with nothing else
non-draft.

- Open your PR as a **draft** and keep it there until you believe it is done.
- Before flipping to ready, list open non-draft PRs. **At 4 or more, stay in draft** and
  re-check after something merges. Four, not five, leaves a slot for a re-run.

| Order | PR | Why here |
|---|---|---|
| 1 | **W0** | Alone. Three labels, several full workflow runs. Nothing else ready until it merges. |
| 2 | **W1** | Owner merge, and must be alone among ranking changes. |
| 3–5 | **W3 / W4 / W5** | No labels, no protected paths. Ready in session-number order as slots free. |
| 6 | **W2** | Last by construction — it rebases onto W1 regardless, so being last costs nothing. |

## The risk that spans the wave: G7

End-to-end today is **7.37 s** at 500 aa against `G7_SECONDS = 10.0`
(`packages/engine/tests/design/test_timing.py:47`). Two independent forces push on that bar
this wave, and **no single session can see the total**:

- **W1 re-arms the `codon_adaptation` sweep axis.**
  `docs/decisions/2026-09-02-ranking-increment.md` measured that solve at **~2.4 s** when it
  merely rediscovered the unsteered design. With a real usage table it produces a
  *different* CDS and may trip repair at **~2.5 s** instead. 7.37 + 2.4 ≈ **9.8 s**. W1 will
  very likely be the first session to see that assertion fail.
- **W3–W5 add five rules**, and Tier-B repair evaluates up to `max_candidates` candidates
  per iteration with a **full catalog pass each**. Every new rule multiplies into the
  dominant cost term.

**The mitigation, and the wrong fixes are pre-refused.** Never raise `G7_SECONDS`. Never
mark the test `slow` — `.claude/rules/tests.md` records that `-m "not slow"` currently
deselects nothing, so marking it would make it the first use of a marker whose only effect
is that `gates.sh` and CI never run it. A timing gate nothing executes is worse than none.
The lever is the **number of solves**: `DEFAULT_SWEEP_STEPS` stays 1.

If four live axes plus five new rules still exceed the budget, that is a **finding about the
budget**, reported with the measured number — and it is the outcome PLAN itself anticipated.
`docs/PLAN.md:508` gives G7's fail consequence verbatim: *"re-allocate budget before rules
multiply cost."* This wave is that sentence coming true.

**Every session measures G7 before and after on its own branch and reports the delta in its
PR.** W2, landing last, carries the cumulative number.

## Inter-session contracts

What one session may rely on another not breaking.

- **`design()`'s signature stays frozen.** Keyword-only, `table_id` never defaulted. W2 may
  *add* fields to `SkeletonResult`; it may not remove or rename one. `bt5/cli.py` is built
  against it and no session in this wave owns the CLI.
- **`BiosecurityVerdict`'s shape stays frozen** (`core/context.py:97`). W2 changes what the
  *report* says about a `not_run` verdict. It does not touch the type, and it does not
  remove `design(screen=...)`.
- **Rule registration stays autodiscovery.** `core/registry.py` walks `bt5.rules.catalog`
  with `pkgutil`; `rules/__init__.py` and `rules/catalog/__init__.py` are both empty, so
  adding a rule edits **zero** shared files. No session may introduce a hand-maintained
  list — that would create the one collision this design avoids.
- **No new rule can break a preset, and this is checked.** The three shipped presets weight
  only `2.B1`, `2.C1`, `2.C3`, `2.D4`, `2.F2` — all already built — so no rule W3/W4/W5 adds
  can make `resolve()` raise. But **do not reuse a `brief_ref`**: `_index_by_ref`
  (`presets.py:120-127`) raises when two rules claim one, and that breaks every preset at
  import. W2 adds a `tests/data_integrity` assertion that all three presets resolve cleanly
  against the live catalog, so a violation fails a gate rather than a user's run.
- **A new SOFT rule is scored but unweighted, and that is expected.** Under W2's policy a
  rule no preset names reports its percentile at weight 0.0, flagged as outside the preset's
  claim. If your rule *should* count for a modality, **open an issue** naming the
  `brief_ref`, the proposed weight and the note. `score/presets.py` is W2's file.

## Shared rules

`CLAUDE.md` loads automatically and is not repeated here. What it does **not** say:

- **Bootstrap first.** A fresh checkout has no `.venv`; run `/bootstrap`. Every command uses
  `.venv/bin/…`. `gates.sh` exit **10** means no venv — BROKEN, not a code failure. Bare
  `pytest` exits 4 on a `conftest.py` import error and looks like a real failure.
- **`/pre-pr` cannot be self-invoked.** It is `disable-model-invocation: true`; the operator
  runs it, and a session must not replicate its steps by other means. After it runs and the
  final push lands, the attestation comment `/pre-pr <head-sha>` goes on the PR. An
  attestation names one commit and a later push makes it stale by design. Never attest a SHA
  that was not just reviewed.
- **Check main is green before you start, and do not trust a SHA written here.** At the time
  of writing main was `574ea0e`.
- **Decisions are one file per decision** under `docs/decisions/`, named
  `YYYY-MM-DD-slug.md`. Never append to a shared file. That directory replaced a single
  `docs/decisions.md` in #85 for exactly the reason this wave is split by write-ownership.
- **Spend context on judgment, not retrieval.** Route retrieval to `Explore` and
  `docs-miner`, gates to `gate-runner`, and keep the main thread for decisions only it can
  make.
- **Escalate on the right axis.** *Capability failure raises the model; diligence failure
  raises the effort.* If a first-pass fix already failed, that is `debugger`, not another
  attempt.

## Launch table

Sessions inherit `model: opus`, `effortLevel: high` from `.claude/settings.json`. Only the
**bold** cells are overrides.

| # | Permission mode | Model | Effort | Unattended? | Design note first? |
|---|---|---|---|---|---|
| W0 | `default` | opus | high | **no** — holds four protected paths | no |
| W1 | `acceptEdits` | opus | **xhigh** | yes | **yes** |
| W2 | `acceptEdits` | opus | **xhigh** | yes | **yes** |
| W3 | `acceptEdits` | opus | **xhigh** | yes | no |
| W4 | `acceptEdits` | opus | high | yes | no |
| W5 | `acceptEdits` | opus | high | yes | no |

**Not plan mode**, for the reason Wave 1's README gives: six plan-mode gates turn the
parallelism into a queue, and this repo already gates on `/pre-pr` → draft PR → owner merge.

**W1 and W2 carry the substitute** — each posts a **short design note as the first comment
on its own draft PR before implementing**, because each decides something the owner would
want to redirect while it is still cheap: W1 the G7 budget, W2 the semantics of
`is_complete`. That is the redirect opportunity plan mode would have given, at the point it
is still cheap, without stalling three other sessions behind one approval.

**Not `bypassPermissions` either.** `.claude/hooks/protect_paths.py` is what stops a session
editing `core/`, `verify.py`, `.github/`, `data/` or `pyproject.toml` without stopping to
think. W0 runs in `default` mode and needs a person within reach, because its work is
*defined* as editing protected paths.

## Merging

**No session in this wave self-merges.** W0 carries three `approved:*` labels — and §7b is
explicit that a label is sign-off on the change, not a licence to merge it unreviewed. W1
and W2 have non-"none" scientific impact, and W2 additionally changes what the app says it
did not do. Each rules PR changes emitted sequences.

Say in your PR why it qualifies for owner merge rather than self-merge, in the scientific
impact section.

## Out of scope

Stated so no session wanders into it:

- **M9 `packages/server/` and M10 `apps/web/`.** They need `pyproject.toml` (mypy `files`,
  hatch `packages`, pytest `testpaths`), and M10 additionally forces an always-run node job:
  `required-checks` fails on `skipped` (`ci.yml:315`), so a web job gated on the existing
  unconsumed `apps/web/**` filter would be `skipped` on every Python-only PR and deadlock
  the gate forever. Wave 3, after the scorecard is honest.
- **`core/result.Relaxation`** — costed ways out of a conflict. See W2's prompt for why it is
  deferred and what ships instead.
- **Issue #45's ε-constraint rearchitecture** (X2–X5, X7).
- **The oracle backlog:** #69, #70 (`verify.py`, `approved:oracle-change`), #52, #64.
- **`benchmarks/`** — the directory does not exist at all, despite `CLAUDE.md` §2 protecting
  `baseline.json` and `tolerances.yaml`. Creating it is an owner decision under
  `approved:algorithm-change`, and building a baseline from numbers about to move twice is
  wasted work. A Wave 3 decision, taken after the ranking stops moving.
- **#56** — no published IDT error-free length for gBlocks. Only the owner can close it, by
  finding the figure. Not a code task.
