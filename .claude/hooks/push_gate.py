#!/usr/bin/env python3
"""PreToolUse: ask before pushing a tree /pre-pr has not gated.

CI capacity is this repo's binding constraint -- CLAUDE.md section 4: 20 concurrent job
slots, ~12 per Python PR, so at most 5 open non-draft pull requests. Pushing a red PR is
therefore the most expensive mistake available here, and it is a scheduling cost the
model cannot observe from inside a turn.

KEYED ON TREE STATE, NOT COMMIT RECENCY. "Has /pre-pr run since the last commit?"
inverts the correct order -- gate, commit, push -- so it would fire on every correct
workflow and pass only if you re-ran the chain after committing. A gate that fires on
the correct path gets disabled within a week. Instead /pre-pr records HEAD plus a hash
of the working tree, and this compares against exactly that.

Also matches the draft->ready flip (update_pull_request): ci.yml triggers on
`ready_for_review` and `synchronize`, and drafts skip the expensive jobs, so that call
is the one that actually spends the slots.

Decision is "ask", never "deny". Fails open.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys

MARKER = ".claude/.pre-pr-marker"
PR_TOOLS = {
    "mcp__github__create_pull_request",
    "mcp__github__update_pull_request",
    "mcp__github__update_pull_request_branch",
}


def _git(args: list[str], cwd: str) -> str | None:
    try:
        r = subprocess.run(
            ["git", *args],
            cwd=cwd,
            capture_output=True,
            text=True,
            encoding="utf-8",
            timeout=10,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return r.stdout if r.returncode == 0 else None


def tree_state(cwd: str) -> str | None:
    head = _git(["rev-parse", "HEAD"], cwd)
    status = _git(["status", "--porcelain"], cwd)
    diff = _git(["diff", "HEAD"], cwd)
    if head is None or status is None or diff is None:
        return None
    digest = hashlib.sha256((status + diff).encode("utf-8", "replace")).hexdigest()[:16]
    return f"{head.strip()}:{digest}"


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except Exception:
        return 0

    tool = payload.get("tool_name")
    tool_input = payload.get("tool_input") or {}
    if not isinstance(tool_input, dict):
        return 0

    if tool == "Bash":
        cmd = tool_input.get("command")
        if not isinstance(cmd, str) or "git push" not in cmd:
            return 0
    elif tool not in PR_TOOLS:
        return 0

    cwd = payload.get("cwd") or os.getcwd()
    current = tree_state(cwd)
    if current is None:
        return 0  # not a git repo, or git unavailable -- do not get in the way

    recorded = None
    marker_path = os.path.join(cwd, MARKER)
    try:
        with open(marker_path, encoding="utf-8") as fh:
            recorded = fh.read().strip()
    except OSError:
        pass

    if recorded == current:
        return 0

    why = "has never been run on this tree" if recorded is None else "last ran on a different tree"
    json.dump(
        {
            "hookSpecificOutput": {
                "hookEventName": "PreToolUse",
                "permissionDecision": "ask",
                "permissionDecisionReason": (
                    f"/pre-pr {why}. CI capacity is the binding constraint here "
                    f"(CLAUDE.md section 4: 20 slots, ~12 jobs per Python PR, at most 5 open "
                    f"non-draft PRs), so a red push is expensive. Run /pre-pr first, or "
                    f"approve to push anyway."
                ),
            }
        },
        sys.stdout,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
