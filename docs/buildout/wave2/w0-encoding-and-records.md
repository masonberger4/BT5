# W0 — text encoding, the guard, and the record that never landed

*Copy this whole file into a fresh Claude Code session on the BT5 repo.*

**Launch:** permission mode `default` · model opus · effort high (both the repo default) ·
**needs a person within reach** — your work is defined as editing protected paths.
**Do not run this in plan mode.** The gate is the draft PR.

**You run alone.** No other Wave 2 session is branched while this PR is open. Five sessions
are waiting on you, so keep the scope exactly as written and do not widen it.

---

You are fixing issue **#102**: every text file this repo reads or writes uses Python's
*locale-dependent* default encoding. On Linux CI that is UTF-8, so nothing surfaces. On a
Windows checkout it is cp1252, and the same command against the same commit behaves
differently.

## Read this first

Run `/bootstrap` — a fresh checkout has no `.venv`, and `gates.sh` exit **10** means BROKEN,
not a code failure. Then `CLAUDE.md`, then `docs/buildout/wave2/README.md`.

Your branch: **`claude/w0-encoding-and-records`**. Cut it yourself, from a **freshly
fetched** main:

```bash
git remote prune origin
git fetch -q origin main && git checkout -B claude/w0-encoding-and-records origin/main
```

`git fetch --prune origin main` does **not** clear a stale ref — pruning is bounded by the
refspec, so that form prunes only `origin/main` (CLAUDE.md §7a).

## Why this is sequenced first, alone

Three reasons, in descending order of force. You do not need to re-derive them, but you do
need to not undermine them.

1. **It is a mutex against every other session by file.** The un-encoded sites are
   distributed across `tests/contract/`, `packages/engine/tests/rules/`,
   `tests/data_integrity/`, `tests/vector/`, four files in engine `src/`, and
   `.github/scripts/`. W1 owns `codon/tables.py`; W2 owns `score/order.py` and
   `tests/{design,score}/`; W3–W5 own `tests/rules/`. Every session collides with you.
2. **Last is strictly worse than first.** New un-encoded sites are this repo's habit —
   `tests/data_integrity/test_no_expression_claims.py` reads with a bare `read_text()`
   *inside the honesty gate itself*. Running this last means re-sweeping a target five
   sessions just moved, and leaves the Windows failure in place for the whole wave.
3. **`score/order.py` writes the vendor order CSV.** A mojibake'd order file is a wrong
   tube — the failure class `Candidate.design_hash` exists to prevent.

## The reproduction, for your own confidence

On a Windows checkout of `main`, `/bootstrap` completed, gates run:

```
packages/engine/tests/score/test_presets.py:314: in test_every_weighted_ref_is_a_real_row_in_the_brief
    brief = (Path(__file__).resolve().parents[4] / "docs" / "research" / "brief.md").read_text()
E   UnicodeDecodeError: 'charmap' codec can't decode byte 0x9d in position 9402
```

`brief.md` is valid UTF-8 containing em-dashes and `×`. `PYTHONUTF8=1 .venv/bin/pytest`
on the same node id passes. That one test is the entire delta; the gate chain is otherwise
green on Windows apart from the G7 timing budget, which is machine-speed dependent and
unrelated (`CLAUDE.local.md` §3).

## What to build

### 1. `encoding="utf-8"` at every text-I/O call site

Find them first, and do not trust this prompt's counts:

```bash
git grep -n "encoding=" -- "*.py"
```

returns exactly **one** hit repo-wide today (`.claude/hooks/push_gate.py:82`), and it is a
Claude hook, not application code. Everything else — `Path.read_text()`,
`Path.write_text()`, bare `open()` on a text mode — inherits `locale.getencoding()`.

**This is `/cheap-pass` work.** `CLAUDE.md` routes ≥5 identical edits with the before/after
already decided to `batch-editor`. This is 30+ identical edits. Do not hand-pass it in the
main window — that is exactly the context-burn the delegation table exists to prevent.

Binary reads (`"rb"` / `"wb"`) take **no** `encoding=`. If the batch pass touches one,
that is a defect in the pass, not a file that needed it.

### 2. The guard that makes regression impossible — the real deliverable

A one-time sweep decays. What ships alongside it is a check that rejects the next bare
`read_text()`.

**Preferred:** ruff's unspecified-encoding rule in `[tool.ruff.lint] select` in
`pyproject.toml`. Check whether the pinned ruff has it out of preview before committing to
it — if it is preview-only, do not enable preview mode to get it.

**Fallback, and it is a good one:** a `grep` step inside the **existing** `python-quality`
job in `.github/workflows/ci.yml`, beside the "Ban global RNG state" step at
`ci.yml:76-89`. That step is the precedent for exactly this shape of guard, and copying it
means:

- **no new CI job**, so no `required-checks.needs` edit,
- **no new job slots** against a 20-slot cap in the wave whose binding constraint is CI
  capacity,
- and `approved:ci-change` shrinks from a new workflow to a few lines in one job.

Either way, prove the guard works: add a deliberately bare `read_text()` locally, watch the
check fail, remove it. Say in the PR that you did.

### 3. The decision record that never landed

Write `docs/decisions/2026-09-03-biosecurity-screen-dropped.md`.

`main` has **no record** that BT5 deliberately does not screen. PR #87 built a protein-level
biosecurity screen behind a `Screen` protocol and the owner closed it unmerged with:

> Closing unmerged. The owner determined BT5 doesn't need its own protein-level biosecurity
> screen — DNA synthesis vendors already screen orders before synthesis, which is where the
> actual risk gate belongs.

The original record was written on that branch and died with the closed PR. Reconstruct it
on `main`: what was built, why it was dropped, that `core/context.BiosecurityVerdict` is
untouched and stays frozen, and that `design(screen=...)` remains for a caller who *has*
screened elsewhere.

It lands here rather than in W2 so the record exists before five sessions read the repo and
wonder why the screen is absent. **W2 changes the behaviour; you write the history.** Do not
touch `design/runner.py` — that is W2's file and W2's call.

## What you must not do

- **Do not run `tests/contract/regenerate.py`.** Not even "just to check". It rewrites
  `manifest.json` and 17 fixtures **before it prints anything** and returns 0 on every path,
  so a speculative run destroys the local baseline and reports success. Every edit in this
  PR is on a *read* path, so `pytest tests/contract` must pass **without** regeneration. If
  the contract gate goes red, that is a signal you touched a write path — take it to the
  owner rather than regenerating green. `CLAUDE.md` §4: re-recording a contract fixture is
  not a fix.
- **Do not add `PYTHONUTF8=1` to `scripts/gates.sh` or to CI.** That is the workaround this
  PR exists to remove. `CLAUDE.local.md` §2 documents it as masking a real portability bug.
- **Do not widen scope.** No refactoring of the files you touch, no drive-by fixes. Five
  sessions are blocked behind this PR and every extra line is a review cycle they wait on.

## Files

Owned by you, and only you, for this PR:

- every `.py` with a bare `read_text` / `write_text` / text-mode `open`, wherever it lives
- `pyproject.toml` — `[tool.ruff.lint] select` only, if the ruff route works
- `.github/workflows/ci.yml` — the guard step inside `python-quality` only
- `docs/decisions/2026-09-03-biosecurity-screen-dropped.md` (new)

Do **not** touch `design/runner.py`'s logic, any rule file's logic, or `core/`.

## Labels and merge

This PR carries **three** protected-path labels:

- `approved:contract-change` — `tests/contract/**`
- `approved:oracle-change` — `tests/invariants/**`, `tests/data_integrity/**`
- `approved:ci-change` — `.github/**`

**You do not self-merge.** §7b is explicit: a label is sign-off on the *change*, not a
licence to merge it unreviewed. Say so in the PR.

Note that `ci.yml:11` triggers on `labeled`/`unlabeled`, so **each label application re-runs
the whole workflow**. Three labels is several full runs. That is why nothing else in the
wave is non-draft while you are open — apply the labels in one pass, not one at a time
across a day.

## Delegation

- **The sweep itself** → `/cheap-pass` (`batch-editor`). Non-negotiable; see above.
- **Finding the call sites** → `Explore`, or the `git grep` above run once and read.
- **Gates** → `gate-runner`.
- **PR #87's closing rationale** → `gh pr view 87` — it is in the comment thread, quoted
  above in full. You do not need to read the reverted diff.
- **Never bare-Read a source file over 20 KB** — `vector/backbone.py` (28.9 KB),
  `rules/vendors.py` (25.7), `vector/kmers.py` (21.2) all appear in the sweep. The
  `batch-editor` pass does not need to read them whole.

## Done means

- `bash scripts/gates.sh` → `ALL GATES PASSED`, exit 0.
- On Windows, `packages/engine/tests/score/test_presets.py::test_every_weighted_ref_is_a_real_row_in_the_brief`
  passes **without** `PYTHONUTF8=1`. If you are not on Windows, say so in the PR and name
  what you did check — do not claim a platform you did not run on.
- The new guard rejects a deliberately re-introduced bare `read_text()`, and you say in the
  PR that you watched it fail.
- `pytest tests/contract` passes **without** `regenerate.py` having been run.
- `git grep -n "encoding=" -- "*.py"` now returns a hit for every text-I/O site.
- The decision record exists and names PR #87 and the vendor-screening argument.
- Scientific impact section says **"none"** — byte-identical output on Linux — and explains
  that the change is about *which platforms produce that output*.
- Draft PR, `/pre-pr` from the operator, **owner merge**.
