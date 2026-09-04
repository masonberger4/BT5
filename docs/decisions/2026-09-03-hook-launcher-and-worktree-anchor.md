## 2026-09-03 — hooks launch through run_py.sh; protect_paths anchors on the target's git root

**Decided:** two independent defects made the CLAUDE.md §2 protected-path guard a silent
no-op, and both are fixed in this branch.

1. Every hook in `.claude/settings.json` was registered as `python3 <script>`. On Windows
   `python3` is the Microsoft Store stub — exit 49, empty stdout. Claude Code blocks a
   PreToolUse hook only on exit **2**, so 49 is a *non-blocking* error: one red line to the
   human, never entering the model's context, tool call proceeds. All five python hooks now
   launch through `.claude/hooks/run_py.sh`, which **runs** each candidate and requires
   `Python 3.` on stdout, then dispatches by a tier declared at the registration site —
   `--required` (exit 2, blocks), `--optional` (exit 0 + stderr), `--statusline` (banner on
   stdout), `--which` (shared resolver).
2. `protect_paths.py` was blind across worktree boundaries **even with a working
   interpreter**. Hook `cwd` is the *session's* cwd, and BT5 worktrees live inside the main
   repo at `.claude/worktrees/<slug>/`. `repo_root_for()` now walks up from the *target
   file* to the nearest `.git` (a directory in a checkout, a **file** in a linked worktree).

**Rejected:**

- **A `python3.exe` shim on PATH** (the cheap fix). Repairs one machine, is invisible to
  every other clone and to the committed configuration, and leaves the failure *class*
  intact — a guardrail that fails open and a verifier that passes vacuously. It also does
  nothing about defect 2.
- **Changing the hooks to `python`.** Absent from many Linux images; the web container has
  `python3` and often not `python`. No single name is correct on all four surfaces.
- **A mandatory per-hook attestation sentinel** (the strongest fail-closed proposal).
  Measured at 1877 ms per Bash tool call with hooks executed four times, and an ordinary
  `RuntimeError` inside `protect_paths.py` was misdiagnosed as "no interpreter" and blocked
  Edit/Write *on every surface*. Its own failure mode is worse than the hole it closes, and
  unlike that hole it has no detector. Bought the coverage with positive probes instead.
- **A `BT5_HOOKS_FAIL_OPEN` escape hatch.** A committed, documented, one-variable revert of
  a guardrail is CLAUDE.md §4 by the back door.
- **`BT5_PYTHON` falling back to the normal chain when it is wrong.** Silently running an
  interpreter other than the one an operator named is the same class of bug. It is
  authoritative; `--optional` keeps Bash alive so a stale value is never a brick.
- **Fail-closed on the `Bash` matcher.** Would brick a Python-less machine with no way to
  repair it from inside a session. Tiering is what makes fail-closed survivable.
- **Loosening `verify-setup.sh`'s `[ "$L" -lt 200 ]` to `-le`.** That is §4 loosening. The
  check was correctly reporting a real over-budget condition; CLAUDE.md was reflowed to 199
  instead, keeping the fact that `goldens-not-hand-edited` does not exist yet and dropping
  only the incidental `tests/goldens/` detail.

**Evidence:**

- `python3 --version` → exit **49**, message on stderr, stdout empty. `command -v python3`
  still returns **0** on it, so no existence check can detect this.
- Same payload, two interpreters: `python3 .claude/hooks/protect_paths.py` → exit 49, no
  output; `py .claude/hooks/protect_paths.py` → `{"permissionDecision": "ask", …}`. The
  Python was always correct; only the interpreter name was wrong.
- ~795 `"Python was not found"` records in `~/.claude/projects/C--Users-mason-BT5*`. Not one
  says a guardrail was off.
- Defect 2, measured before the fix with `cwd` = main repo and the file inside this
  worktree: `rc=0, 0 bytes` for `pyproject.toml`, `uv.lock`, `verify.py`, `core/types.py`
  and `.github/workflows/ci.yml`. Reverse direction (`cwd` = worktree, file in main) also
  `rc=0, 0 bytes`. After: all 16 combinations of {cwd main, worktree} × {file main,
  worktree} ask, and `README.md` stays silent from both.
- `verify-setup.sh` before: exit 1 with four FAILs, including a false `invalid JSON` about
  a well-formed `settings.json`, plus seven assertions whose pass condition was *empty
  stdout* — which is exactly what a dead interpreter produces. After: exit 0, 62 `ok`
  lines, and a negative control that proves the verifier can still fail.
- The new SessionStart probe was checked against the old hook: it reports `PARTIAL` with
  the pre-`repo_root_for()` version restored and `ok` with the fix. It is a detector that
  demonstrably detects.
- Windows trap worth recording: under Git Bash `$PWD` is `/c/Users/…`, which a native
  Windows Python resolves to a nonexistent `C:\c\Users\…`. Probes must use `pwd -W`. An
  early version of the test matrix reported a false MISS for exactly this reason.

**Still open, deliberately:** a *gutted* `protect_paths.py` exits 0 with no output under
`exec`, which is byte-identical to a clean pass. Covered by detection (SessionStart probe,
`verify-setup.sh`, and the CI shape check) rather than enforcement. The fail-closed
property itself rests on Claude Code treating exit 2 as blocking — read from the shipped
bundle, **not** executed end to end. That is a merge acceptance criterion, not an
assumption: set `BT5_PYTHON=/nonexistent/python`, attempt a `Write` to `pyproject.toml`,
confirm it is refused.

**Where:** branch `claude/hook-launcher-worktree-anchor`, in four commits — the launcher and
the five rewired command strings; the `repo_root_for()` worktree anchor; the verifier
repair, SessionStart detector and `/pre-pr` scheduling; and the `--optional` exit-2 clamp
from code review. Deliberately no SHAs here: this branch was rebased onto a moving `main`
mid-flight and an earlier version of this file cited commits that no longer existed.

The CI shape check (`.github/scripts/check-hook-commands.py` plus a step in `ci.yml`) is a
**separate, stacked** branch, `claude/ci-hook-command-guard`, because `.github/**` is a §2
protected path and carries `approved:ci-change` — which is owner sign-off, so bundling it
would have forced a label onto an otherwise unlabelled PR.
