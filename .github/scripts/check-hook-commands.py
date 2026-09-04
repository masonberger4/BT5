#!/usr/bin/env python3
"""Fail CI if a hook could silently not run, or has been quietly neutered.

Sibling of check-workflow-gate.py, and it makes the same argument in the opposite
direction. That script guards against a required check that never reports, which blocks
a PR forever with no error. This one guards against a hook that does not run, which does
NOT block at all -- and that is worse, because it is silent.

History: every hook in .claude/settings.json was registered as `python3 <script>`. On
Windows `python3` is the Microsoft Store stub: exit 49, empty stdout. Claude Code blocks
a PreToolUse hook only on exit 2, so the tool call proceeded and nothing recorded that
protect_paths.py -- the sole mechanical enforcement of CLAUDE.md section 2 for the Write
tool and for pyproject.toml / uv.lock -- had not run. ~795 logged occurrences.

TWO KINDS OF CHECK, AND WHY BOTH ARE HERE.

  SHAPE (static): CI runners all carry a live python3, so the python3-stub bug cannot
  reproduce here. Only shape can catch it, and shape is platform-independent.

  BEHAVIOUR (a probe): shape alone constrains SPELLING, not what the guard DECIDES. An
  earlier version of this file asserted only shape, and a security pass defeated it 11
  ways out of 13 -- a one-token matcher edit (`Write` -> `WriteX`), a `bash -c python3`
  wrapper, a trailing `|| true`, relocating the hook group from PreToolUse to
  PostToolUse, and gutting the hook body while leaving the asserted strings in a
  docstring, ALL reported ok. So this file also RUNS the guard and asserts its verdict.

  Note the distinction the earlier version got wrong: an EXECUTABILITY probe ("does
  python3 exist") would indeed be vacuous here, because it passes on every runner while
  the guardrail is dead on a developer's machine. A DECISION probe ("does the guard say
  ask for pyproject.toml") is not vacuous anywhere. .claude/hooks/session_start.sh
  already runs exactly that probe locally; this runs it in CI.

It parses settings.json rather than grepping the repo, because compact_output.py contains
the literal "python3 -m pytest" as DATA in PYTEST_HEADS.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
SETTINGS = ROOT / ".claude" / "settings.json"
HOOKDIR = ROOT / ".claude" / "hooks"
LAUNCHER = ".claude/hooks/run_py.sh"
PROTECT = ".claude/hooks/protect_paths.py"

# A hook command must be EXACTLY one of these two shapes and nothing else. This is a
# positive whitelist on purpose: the previous bare-name blacklist let `bash -c 'python3
# ...'`, an absolute interpreter path, and a trailing `|| true` through. Anything with a
# shell metacharacter, a redirection, a second command or an unexpected interpreter fails
# to match and is reported.
_DIR = r"(?:\$\{CLAUDE_PROJECT_DIR\}/|\$CLAUDE_PROJECT_DIR/)?"
SH_RE = re.compile(rf"^bash\s+{_DIR}(?P<script>\.claude/hooks/[A-Za-z0-9_.-]+\.sh)$")
PY_RE = re.compile(
    rf"^bash\s+{_DIR}\.claude/hooks/run_py\.sh"
    rf"\s+--(?P<tier>required|optional|statusline)"
    rf"\s+{_DIR}(?P<script>\.claude/hooks/[A-Za-z0-9_.-]+\.py)$"
)

failures: list[str] = []


def fail(msg: str) -> None:
    failures.append(msg)


def norm(token: str) -> str:
    """Strip a ${CLAUDE_PROJECT_DIR}/ prefix so a portable path compares equal."""
    return re.sub(r"^(?:\$\{CLAUDE_PROJECT_DIR\}/|\$CLAUDE_PROJECT_DIR/)", "", token)


def commands(settings: dict[str, Any]) -> list[tuple[str, str, str]]:
    """-> (event, matcher, command). Event matters: a PreToolUse guard moved to
    PostToolUse fires AFTER the write, where an "ask" verdict is meaningless."""
    out: list[tuple[str, str, str]] = []
    sl = settings.get("statusLine") or {}
    if isinstance(sl, dict) and isinstance(sl.get("command"), str):
        out.append(("statusLine", "", sl["command"]))
    for event, groups in (settings.get("hooks") or {}).items():
        for group in groups or []:
            matcher = (group or {}).get("matcher")
            for hook in (group or {}).get("hooks") or []:
                if hook.get("type") == "command" and isinstance(hook.get("command"), str):
                    # A matcher that is absent means "every tool", i.e. strictly broader
                    # coverage. Represent it as "*" rather than "" so it is never
                    # mistaken for a matcher that names nothing.
                    out.append((event, "*" if matcher is None else str(matcher), hook["command"]))
    return out


def check_shape(cmds: list[tuple[str, str, str]]) -> set[str]:
    referenced: set[str] = set()
    for event, matcher, cmd in cmds:
        site = f"{event}[{matcher}]"
        if not cmd.strip():
            fail(f"{site}: empty command string.")
            continue
        m = PY_RE.match(cmd.strip())
        if m:
            referenced.add(norm(m.group("script")))
            referenced.add(LAUNCHER)
            continue
        m = SH_RE.match(cmd.strip())
        if m:
            script = norm(m.group("script"))
            if script == LAUNCHER:
                fail(f"{site}: invokes {LAUNCHER} with no tier flag.")
            else:
                referenced.add(script)
            continue
        fail(
            f"{site}: command is not a recognised hook shape: {cmd!r}. "
            f"It must be exactly `bash .claude/hooks/<name>.sh` or "
            f"`bash {LAUNCHER} --required|--optional|--statusline .claude/hooks/<name>.py`. "
            f"Anything else -- a bare or absolute interpreter, `bash -c`, a redirection, "
            f"`|| true`, a second command -- can make the hook silently not run, or make "
            f"its exit status stop meaning what the tier says it means."
        )
    return referenced


def check_protect_paths(cmds: list[tuple[str, str, str]], settings: dict[str, Any]) -> None:
    sites = [(e, m, c) for e, m, c in cmds if PROTECT in c]
    if not sites:
        fail(f"{PROTECT} is not registered by any hook in settings.json.")
        return
    ok_site = False
    for event, matcher, cmd in sites:
        if event != "PreToolUse":
            fail(
                f"{PROTECT} is registered on {event}, not PreToolUse. A guard that runs "
                f"after the tool cannot stop it; an 'ask' verdict there is meaningless."
            )
            continue
        m = PY_RE.match(cmd.strip())
        if not m or m.group("tier") != "required":
            fail(
                f"PreToolUse[{matcher}]: {PROTECT} must be launched --required. It is the "
                f"only mechanical enforcement of CLAUDE.md section 2 for the Write tool "
                f"and for pyproject.toml / uv.lock, and no CI job covers it -- "
                f"claude-review grants no Edit/Write tool, so it never runs on a runner."
            )
            continue
        # Exact alternation membership, NOT a substring test: `WriteX`, `WriteFile` and
        # `NotebookWrite` all contain "Write" while protecting the Write tool from nothing.
        tokens = [t for t in re.split(r"[|\s]+", matcher) if t]
        if matcher != "*" and "Write" not in tokens:
            fail(
                f"PreToolUse[{matcher}]: the matcher does not name the Write tool as its "
                f"own alternation token. Write is the entire reason this hook exists: "
                f"Claude Code accepts a Write(path) permission rule and never consults it."
            )
            continue
        ok_site = True
    if not ok_site:
        fail(f"{PROTECT} has no PreToolUse registration that covers Write with --required.")

    # run_py.sh treats BT5_PYTHON as AUTHORITATIVE with no fallback, so a committed value
    # pointing at a shim that answers -V and does nothing else disables every hook here.
    env = settings.get("env") or {}
    if isinstance(env, dict) and "BT5_PYTHON" in env:
        fail(
            "settings.json sets env.BT5_PYTHON. run_py.sh treats it as authoritative with "
            "no fallback, so a committed value silences every python hook in the repo. "
            "It is a per-machine escape hatch for settings.local.json, never a committed one."
        )


def probe(payload: dict[str, Any], env_extra: dict[str, str] | None = None) -> tuple[int, str]:
    """Run the REAL guard through the REAL launcher and return (rc, stdout)."""
    env = dict(os.environ)
    env.update(env_extra or {})
    proc = subprocess.run(
        ["bash", LAUNCHER, "--required", PROTECT],
        cwd=ROOT,
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        env=env,
    )
    return proc.returncode, proc.stdout


def check_behaviour() -> None:
    """Shape cannot see a gutted hook, a `|| true`, or a downgraded run_py.sh. This can."""
    root = str(ROOT)
    target = str(ROOT / "pyproject.toml")

    rc, out = probe({"tool_name": "Write", "tool_input": {"file_path": target}, "cwd": root})
    if '"ask"' not in out:
        fail(
            f"the guard returned no 'ask' for a Write to pyproject.toml (rc={rc}, "
            f"stdout={out[:120]!r}). Shape is satisfied but the guard DECIDES NOTHING -- "
            f"check-approval-labels.sh does not enforce pyproject.toml/uv.lock either, so "
            f"nothing else is covering this."
        )

    # Off-root cwd: hook `cwd` is the SESSION's cwd, and worktrees live inside the repo at
    # .claude/worktrees/<slug>/. This is the one probe that fails on a pre-repo_root_for()
    # hook; a probe using the same value for cwd and file_path tests the configuration
    # that always worked.
    rc, out = probe(
        {"tool_name": "Write", "tool_input": {"file_path": target}, "cwd": str(ROOT.parent)}
    )
    if '"ask"' not in out:
        fail(
            f"the guard returned no 'ask' when the session cwd is not the repo root "
            f"(rc={rc}). A session started outside a worktree would write to protected "
            f"paths unprompted. See repo_root_for() in protect_paths.py."
        )

    rc, out = probe(
        {"tool_name": "Write", "tool_input": {"file_path": str(ROOT / "README.md")}, "cwd": root}
    )
    if out.strip():
        fail(f"the guard fired on README.md, a benign path (stdout={out[:120]!r}).")

    # Negative control. Proves --required still BLOCKS, i.e. that run_py.sh's exit 2 arm
    # has not been downgraded to exit 0 -- a one-character diff nothing else here sees.
    rc, _ = probe(
        {"tool_name": "Write", "tool_input": {"file_path": target}, "cwd": root},
        {"BT5_PYTHON": "/nonexistent/python"},
    )
    if rc != 2:
        fail(
            f"with a dead interpreter the --required tier exited {rc}, not 2. Claude Code "
            f"blocks a PreToolUse hook ONLY on exit 2, so the guard would fail OPEN and "
            f"silently -- exactly the original defect. See the --required arm of run_py.sh."
        )


def check_modes(referenced: set[str]) -> None:
    # Ask git, not the filesystem. On Windows core.filemode=false makes every hook report
    # 0644 to stat() while the index correctly holds 100755, so a filesystem check yields
    # false failures locally and cannot be run before pushing. The index mode is also the
    # one that matters: it is what a POSIX checkout materialises.
    want = sorted(referenced | {LAUNCHER})
    proc = subprocess.run(
        ["git", "ls-files", "-s", "--", *want], cwd=ROOT, capture_output=True, text=True
    )
    if proc.returncode != 0:
        fail(f"git ls-files failed ({proc.returncode}); cannot verify hook file modes.")
        return
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


def main() -> int:
    settings = json.loads(SETTINGS.read_text(encoding="utf-8"))
    cmds = commands(settings)

    referenced = check_shape(cmds)
    check_protect_paths(cmds, settings)

    on_disk = {
        f".claude/hooks/{p.name}"
        for p in HOOKDIR.iterdir()
        if p.suffix in (".py", ".sh") and p.name != "run_py.sh"
    }
    for orphan in sorted(on_disk - referenced):
        fail(
            f"{orphan} exists but no hook in settings.json runs it "
            f"-- somebody believes it is installed."
        )

    for path in sorted(referenced):
        if not (ROOT / path).exists():
            fail(f"a hook references {path}, which is not on disk.")

    check_modes(referenced)

    # Only probe if the wiring is sane; otherwise the probe's failure is a confusing
    # restatement of a shape failure already reported above.
    if not failures:
        check_behaviour()

    for f in failures:
        print(f"FAIL {f}", file=sys.stderr)
    if failures:
        return 1
    print(f"ok  {len(cmds)} hook commands; guard probed live and answers for pyproject.toml")
    return 0


if __name__ == "__main__":
    sys.exit(main())
