#!/usr/bin/env python3
"""PreToolUse(Bash): shrink noisy output by appending NATIVE flags.

Why flags and not a pipe. Rewriting `pytest ...` into `pytest ... | filter` makes the
tool report the FILTER's exit status, so a failing suite comes back green and every
decision downstream is built on that. Native flags reduce output before it is produced
and leave the exit code untouched.

Conservative by construction:
  * only a single simple command is touched -- anything containing && || ; | > < &
    backticks or $( is left exactly as written, and so is a backgrounded run;
  * the command string is never re-tokenized, only appended to, so quoting survives
    (`pytest -k "not slow and vector"` must come through byte-for-byte);
  * an explicit user flag is never overridden -- `-v` and `-q` share one verbosity
    counter so appending `-q` would silently cancel a deliberate `-v`, and `--maxfail`
    is last-flag-wins;
  * `tests/contract` and `tests/data_integrity` are never capped with --maxfail: their
    full output is what classifies a change MINOR or MAJOR;
  * mypy, regenerate.py, check_amendment.py, gates.sh and git are never touched.

Fails open: any unexpected input exits 0 with no JSON and the command runs unchanged.
"""

from __future__ import annotations

import json
import re
import shlex
import sys

# Shell metacharacters that make this more than one simple command.
COMPOUND = re.compile(r"(\|\||&&|[;|<>&`]|\$\()")

# Never rewrite these, for the reasons in the module docstring.
NEVER = (
    "mypy",
    "regenerate.py",
    "check_amendment.py",
    "check-workflow-gate.py",
    "check-approval-labels.sh",
    "gates.sh",
)

PYTEST_HEADS = ("pytest", "python -m pytest", "python3 -m pytest", ".venv/bin/pytest")
FULL_OUTPUT_TARGETS = ("tests/contract", "tests/data_integrity")


def _head_matches(cmd: str, heads: tuple[str, ...]) -> bool:
    return any(cmd == h or cmd.startswith(h + " ") for h in heads)


def rewrite(cmd: str) -> str | None:
    """Return the rewritten command, or None to leave it alone."""
    stripped = cmd.strip()
    if not stripped or COMPOUND.search(stripped):
        return None
    if any(token in stripped for token in NEVER):
        return None

    try:
        argv = shlex.split(stripped)
    except ValueError:
        return None  # unbalanced quotes; do not touch it
    if not argv:
        return None

    extra: list[str] = []

    if _head_matches(stripped, PYTEST_HEADS):
        # All-or-nothing. -v, -q, --tb, --maxfail and -x all express a deliberate
        # decision about how much output the caller wants, and they interact:
        # -v and -q share one verbosity counter, so appending -q would silently
        # cancel a deliberate -v, and --maxfail is last-flag-wins, so appending
        # --maxfail=10 would override a deliberate --maxfail=1. If the caller has
        # said anything about output, stand down completely.
        if any(
            a in ("-v", "-vv", "-vvv", "--verbose", "-q", "--quiet", "-x")
            or a.startswith("--tb")
            or a.startswith("--maxfail")
            for a in argv
        ):
            return None
        # These suites' complete output IS the deliberation: the contract change
        # list is what classifies a change MINOR or MAJOR. Never truncate them.
        if any(t in stripped for t in FULL_OUTPUT_TARGETS):
            return None
        extra += ["--tb=short", "-q", "--maxfail=10"]

    elif _head_matches(stripped, ("ruff check", ".venv/bin/ruff check")):
        if not any(a.startswith("--output-format") for a in argv):
            extra += ["--output-format=concise"]

    elif stripped.startswith("uv pip install"):
        if not any(a in ("-q", "--quiet", "-v", "--verbose") for a in argv):
            extra += ["-q"]

    if not extra:
        return None
    return stripped + " " + " ".join(extra)


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except Exception:
        return 0

    if payload.get("tool_name") != "Bash":
        return 0
    tool_input = payload.get("tool_input") or {}
    if not isinstance(tool_input, dict):
        return 0
    if tool_input.get("run_in_background"):
        return 0

    command = tool_input.get("command")
    if not isinstance(command, str):
        return 0

    new_command = rewrite(command)
    if new_command is None or new_command == command:
        return 0

    updated = dict(tool_input)
    updated["command"] = new_command
    json.dump(
        {
            "hookSpecificOutput": {
                "hookEventName": "PreToolUse",
                "updatedInput": updated,
            }
        },
        sys.stdout,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
