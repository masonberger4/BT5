#!/usr/bin/env python3
"""Fail CI if a hook could silently not run.

Sibling of check-workflow-gate.py, and it makes the same argument in the opposite
direction. That script guards against a required check that never reports, which blocks
a PR forever with no error. This one guards against a hook whose interpreter is missing,
which does NOT block at all -- and that is worse, because it is silent.

History: every hook in .claude/settings.json was registered as `python3 <script>`. On
Windows `python3` is the Microsoft Store stub: exit 49, empty stdout. Claude Code blocks
a PreToolUse hook only on exit 2, so the tool call proceeded and nothing recorded that
protect_paths.py -- the sole mechanical enforcement of CLAUDE.md section 2 for the Write
tool and for pyproject.toml / uv.lock -- had not run. ~795 logged occurrences.

CI runners all carry a live python3, so this bug CANNOT reproduce here. Do not add an
executability probe; it would pass on every runner while the guardrail is dead on a
developer's machine, i.e. a third vacuous verifier. What CI can assert is SHAPE, which
is platform-independent. That is all this file does.

Note it parses settings.json rather than grepping the repo: compact_output.py contains
the literal string "python3 -m pytest" as DATA in PYTEST_HEADS.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
SETTINGS = ROOT / ".claude" / "settings.json"
HOOKDIR = ROOT / ".claude" / "hooks"
LAUNCHER = ".claude/hooks/run_py.sh"
BARE = {"python", "python3", "python3.11", "py", "node", "npx", "sh", "bash", "ruby", "perl"}

failures: list[str] = []


def fail(msg: str) -> None:
    failures.append(msg)


def commands(settings: dict[str, Any]) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    sl = settings.get("statusLine") or {}
    if isinstance(sl, dict) and isinstance(sl.get("command"), str):
        out.append(("statusLine", sl["command"]))
    for event, groups in (settings.get("hooks") or {}).items():
        for group in groups or []:
            for hook in (group or {}).get("hooks") or []:
                if hook.get("type") == "command" and isinstance(hook.get("command"), str):
                    out.append((f"{event}[{group.get('matcher', '*')}]", hook["command"]))
    return out


def main() -> int:
    settings = json.loads(SETTINGS.read_text(encoding="utf-8"))
    cmds = commands(settings)
    referenced: set[str] = set()

    for site, cmd in cmds:
        first = cmd.split()[0]
        if first in BARE and first != "bash":
            fail(
                f"{site}: command starts with the bare interpreter name {first!r}. "
                f"A hook whose interpreter is missing does not block -- it is SILENT. "
                f"Launch it as `bash {LAUNCHER} --required|--optional|--statusline <script>`."
            )
        for token in cmd.split():
            if token.startswith(".claude/hooks/"):
                referenced.add(token)
                if not (ROOT / token).exists():
                    fail(f"{site}: references {token}, which is not on disk.")
        if LAUNCHER in cmd and not re.search(r"\s--(required|optional|statusline|which)\s", cmd):
            fail(f"{site}: uses {LAUNCHER} with no tier flag.")

    on_disk = {
        f".claude/hooks/{p.name}"
        for p in HOOKDIR.iterdir()
        if p.suffix in (".py", ".sh") and p.name not in ("run_py.sh", "watch_open_prs.sh")
    }
    for orphan in sorted(on_disk - referenced):
        fail(
            f"{orphan} exists but no hook in settings.json runs it "
            f"-- somebody believes it is installed."
        )

    pp = ".claude/hooks/protect_paths.py"
    guarded = [c for _, c in cmds if pp in c and "--required" in c]
    if not guarded:
        fail(
            f"{pp} must be launched with --required. It is the only mechanical enforcement "
            f"of CLAUDE.md section 2 for the Write tool and for pyproject.toml / uv.lock, "
            f"and no CI job covers it -- claude-review grants no Edit/Write tool, so it "
            f"never runs on any runner."
        )

    for group in (settings.get("hooks") or {}).get("PreToolUse") or []:
        if pp in json.dumps(group) and "Write" not in (group.get("matcher") or ""):
            fail(
                "the protect_paths matcher no longer names Write. Write is the entire "
                "reason this hook exists: Claude Code accepts a Write(path) permission "
                "rule and never consults it."
            )

    src = (ROOT / pp).read_text(encoding="utf-8") if (ROOT / pp).exists() else ""
    for needed in (r"pyproject\.toml", r"uv\.lock", "def repo_root_for"):
        if needed not in src:
            fail(
                f"{pp} no longer contains {needed!r}. "
                f"check-approval-labels.sh says verbatim it does not enforce "
                f"pyproject.toml/uv.lock, and repo_root_for() is what makes the hook "
                f"see files inside .claude/worktrees/<slug>/ at all."
            )

    # Ask git, not the filesystem. On Windows core.filemode=false, so every one of these
    # files reports mode 0644 to stat() while the INDEX correctly records 100755 -- a
    # filesystem check produces eight false failures locally and cannot be run before
    # pushing. The index mode is also the thing that actually matters: it is what a POSIX
    # checkout materialises, and what verify-setup.sh's `[ -x ]` will judge there.
    want = sorted(referenced | {LAUNCHER})
    proc = subprocess.run(
        ["git", "ls-files", "-s", "--", *want], cwd=ROOT, capture_output=True, text=True
    )
    modes = {}
    for line in proc.stdout.splitlines():
        meta, _, name = line.partition("\t")
        if meta:
            modes[name.strip()] = meta.split()[0]
    for path in want:
        mode = modes.get(path)
        if mode is None:
            fail(f"{path} is referenced by a hook but is not tracked by git.")
        elif mode != "100755":
            fail(
                f"{path} is committed as {mode}, not 100755. verify-setup.sh asserts "
                f"[ -x ] on every file in .claude/hooks/, so this passes on Windows "
                f"(core.filemode=false hides it) and is red on every POSIX checkout. "
                f"Fix: git update-index --chmod=+x {path}"
            )

    for f in failures:
        print(f"FAIL {f}", file=sys.stderr)
    if failures:
        return 1
    print(f"ok  {len(cmds)} hook commands, all launched through {LAUNCHER}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
