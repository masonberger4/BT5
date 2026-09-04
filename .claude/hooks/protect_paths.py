#!/usr/bin/env python3
"""PreToolUse(Edit|Write): ask before touching a label-protected path.

THE FAILURE THIS EXISTS TO AVOID. `tool_input.file_path` is absolute -- the tool
schema requires it. The patterns in .github/scripts/check-approval-labels.sh are
`^`-anchored and repo-relative. An anchored relative regex never matches an absolute
path, so an unnormalized version of this hook matches nothing, exits 0, prints
nothing, and looks perfectly installed. Normalization is the whole point.

It covers two gaps that the repo's own tooling does not:
  * `Write` -- Claude Code consults file-path permission rules only for Edit(path)
    and Read(path). A Write(path) rule is accepted, never consulted, and warns at
    startup. A hook matching both Edit and Write is the only mechanism that covers it.
  * pyproject.toml / uv.lock -- check-approval-labels.sh says of them, verbatim,
    "protected by CLAUDE.md section 2 but no label is named for them, so they are
    deliberately NOT enforced here." CLAUDE.md section 5 calls a lockfile conflict
    across parallel PRs the single most expensive merge failure in this repo.

Decision is "ask", never "deny": a subagent cannot answer a prompt, and a block it
cannot satisfy gets routed around. Fails open on any unexpected input.
"""

from __future__ import annotations

import json
import os
import re
import sys

# (compiled pattern, required label) -- mirrors check-approval-labels.sh, plus the
# two entries that script deliberately leaves unenforced.
RULES: list[tuple[re.Pattern[str], str]] = [
    (re.compile(r"^packages/engine/src/bt5/verify\.py$"), "approved:oracle-change"),
    (re.compile(r"^tests/(invariants|data_integrity)/"), "approved:oracle-change"),
    (re.compile(r"^packages/engine/src/bt5/core/"), "approved:contract-change"),
    (re.compile(r"^tests/contract/"), "approved:contract-change"),
    (re.compile(r"^benchmarks/(baseline\.json|tolerances\.yaml)$"), "approved:algorithm-change"),
    (re.compile(r"^data/(genetic_codes|codon_usage)/"), "approved:data-change"),
    (re.compile(r"^\.github/"), "approved:ci-change"),
    (re.compile(r"^(pyproject\.toml|uv\.lock)$"), "NO LABEL EXISTS -- see CLAUDE.md section 5"),
    (re.compile(r"^\.claude/agent-memory/"), "reviewed agent memory"),
]

PATH_KEYS = ("file_path", "notebook_path", "path")


def repo_root_for(path: str, fallback: str) -> str:
    """Anchor on the git root of the TARGET FILE, not the session cwd.

    THE FAILURE THIS EXISTS TO AVOID. Hook `cwd` is the SESSION's cwd, not the project
    root, and BT5 worktrees live INSIDE the main repo at .claude/worktrees/<slug>/.
    A session launched in the main repo that edits a worktree file sends cwd=<main>, so
    relpath() yields ".claude/worktrees/<slug>/pyproject.toml" and no ^-anchored rule
    matches. The reverse -- a worktree session writing to the main checkout -- yields
    "../../../pyproject.toml", which relativize() maps to None. Both are silent passes.
    Both were reproduced: verify.py, core/types.py, pyproject.toml, uv.lock and
    .github/workflows/ci.yml all returned rc=0 with zero bytes.

    Walking up from the file to the nearest `.git` gives the right root in every case
    and costs no subprocess. `.git` is a directory in a normal checkout and a FILE in a
    linked worktree, so os.path.exists covers both. If the file is under no repo at all
    we fall back to `cwd`, i.e. exactly today's behaviour.
    """
    d = os.path.dirname(os.path.realpath(path))
    while True:
        if os.path.exists(os.path.join(d, ".git")):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            return os.path.realpath(fallback)
        d = parent


def relativize(raw: str, cwd: str) -> str | None:
    """Absolute tool path -> repo-relative POSIX path, or None if outside the repo."""
    if not raw:
        return None
    try:
        target = os.path.realpath(raw if os.path.isabs(raw) else os.path.join(cwd, raw))
        root = repo_root_for(target, cwd)
    except OSError:
        return None
    rel = os.path.relpath(target, root)
    if rel.startswith(os.pardir):
        return None
    return rel.replace(os.sep, "/")


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except Exception:
        return 0

    if payload.get("tool_name") not in ("Edit", "Write", "MultiEdit", "NotebookEdit"):
        return 0
    tool_input = payload.get("tool_input") or {}
    if not isinstance(tool_input, dict):
        return 0

    raw = next((tool_input[k] for k in PATH_KEYS if isinstance(tool_input.get(k), str)), None)
    if not raw:
        return 0

    rel = relativize(raw, payload.get("cwd") or os.getcwd())
    if rel is None:
        return 0

    for pattern, label in RULES:
        if pattern.search(rel):
            if label.startswith("NO LABEL"):
                reason = (
                    f"{rel} is protected by CLAUDE.md section 2, and no approved:* label "
                    f"covers it -- CI will NOT stop this. Section 5: every dependency is "
                    f"already declared; a lockfile conflict across parallel PRs is the most "
                    f"expensive merge failure in this repo. Open an issue instead."
                )
            elif label == "reviewed agent memory":
                reason = (
                    f"{rel} is committed agent memory -- an instruction channel that ships "
                    f"in the PR diff. Read the change before accepting it."
                )
            else:
                reason = (
                    f"{rel} is a protected path. The pull request will need the "
                    f"`{label}` label or the `approvals` job fails, blocking "
                    f"`required-checks`. See CLAUDE.md section 2."
                )
            json.dump(
                {
                    "hookSpecificOutput": {
                        "hookEventName": "PreToolUse",
                        "permissionDecision": "ask",
                        "permissionDecisionReason": reason,
                    }
                },
                sys.stdout,
            )
            return 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
