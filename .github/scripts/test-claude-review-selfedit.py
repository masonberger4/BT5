#!/usr/bin/env python3
"""Fixture test for claude-review-gate's verdict step in claude-review.yml.

That step is nine lines of branching shell deciding whether a pull request
passes the review gate, and as of 2026-09-04 it contains THE ONE PATH IN THIS
REPOSITORY'S CI THAT PASSES WITHOUT A VERDICT: a pull request editing
`claude-review.yml` cannot be reviewed at all, because the Claude GitHub App
refuses to mint a token for a workflow file differing from the default branch's
copy. That case used to be red forever -- unclearable by any push -- and now
passes if, and only if, `approved:ci-change` is on the pull request.

A branch that turns "no review happened" into a green check is exactly the kind
this repository has agreed not to ship untested. `test-attestation-matcher.py`
exists because its regex shipped two matching bugs "caught by reading, which
does not scale"; the review on #142 asked for the same treatment for the
`main-broken` predicate and got it. This is the third, and it guards the
loosest of the three.

THE SCRIPT UNDER TEST IS EXTRACTED FROM THE WORKFLOW, NOT COPIED HERE. A test
restating the branching would test the restatement -- which is how the
attestation matcher's second bug shipped. `git` and `gh` are stubbed so the
branches can be driven without a repository or a network; nothing else about
the script is altered.

Two properties are asserted structurally rather than behaviourally, because no
stub can observe them:

  * the label read must use the `pulls` endpoint. `issues/$PR/labels` returns
    the same labels but maps to the `issues` scope, which this workflow does not
    grant -- it would 403 into `set -e` and land straight back on the permanent
    red this branch removes.
  * the fault branch must still fail closed. `2026-09-03-ci-checks-that-can-go-
    green.md` rejected making this job neutral when it cannot produce a verdict;
    the self-edit branch is a narrow exception to that and must not become a
    general one.

Run by ci.yml's `python-quality` job, which is in `required-checks.needs`.
Requires `jq` (preinstalled on GitHub's ubuntu runners); exits 2 rather than
reporting a pass it did not perform.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

WORKFLOW = Path(".github/workflows/claude-review.yml")
STEP_NAME = "Gate on Claude verdict"

CLEAN = '{"blocking":false,"summary":"ok","blocking_findings":[],"provenance":"SUPPORTED"}'
UNSUPPORTED = '{"blocking":false,"summary":"bad","blocking_findings":[],"provenance":"UNSUPPORTED"}'
BLOCKING = '{"blocking":true,"summary":"x","blocking_findings":["y"],"provenance":"N/A"}'

# `git diff --quiet` exits 1 when the path differs, which is what the script
# reads as "this pull request edits claude-review.yml".
GIT_STUB = """#!/usr/bin/env bash
case "$1" in
  cat-file) exit 0 ;;
  diff) [ "${SELF_EDIT_FIXTURE:-no}" = "yes" ] && exit 1 || exit 0 ;;
esac
exit 0
"""

GH_STUB = """#!/usr/bin/env bash
[ "${GH_FAILS:-no}" = "yes" ] && { echo "gh: simulated API error" >&2; exit 1; }
echo "${LABELS_FIXTURE:-[]}"
"""

#: (name, env, expected exit). Every branch of the step, both directions.
CASES: list[tuple[str, dict[str, str], int]] = [
    # --- empty verdict: the three causes, which must NOT collapse together ---
    ("no token configured -> fail closed", {"OUT": "", "HAS_TOKEN": "false"}, 1),
    (
        "crash/timeout/max-turns -> fail closed",
        {"OUT": "", "HAS_TOKEN": "true", "SELF_EDIT_FIXTURE": "no"},
        1,
    ),
    # --- the self-edit branch, which is the reason this file exists ---
    (
        "self-edit WITH approved:ci-change -> pass",
        {
            "OUT": "",
            "HAS_TOKEN": "true",
            "SELF_EDIT_FIXTURE": "yes",
            "LABELS_FIXTURE": '["approved:ci-change"]',
        },
        0,
    ),
    (
        "self-edit, no labels at all -> fail",
        {
            "OUT": "",
            "HAS_TOKEN": "true",
            "SELF_EDIT_FIXTURE": "yes",
            "LABELS_FIXTURE": "[]",
        },
        1,
    ),
    (
        "self-edit, unrelated label -> fail",
        {
            "OUT": "",
            "HAS_TOKEN": "true",
            "SELF_EDIT_FIXTURE": "yes",
            "LABELS_FIXTURE": '["bug","approved:data-change"]',
        },
        1,
    ),
    # A label that CONTAINS the required one as a prefix must not satisfy it --
    # the same trap check-approval-labels.sh's `has_label` matches whole JSON
    # strings to avoid. A substring test would pass this and be wrong.
    (
        "self-edit, prefix-only near-miss label -> fail",
        {
            "OUT": "",
            "HAS_TOKEN": "true",
            "SELF_EDIT_FIXTURE": "yes",
            "LABELS_FIXTURE": '["approved:ci-changes"]',
        },
        1,
    ),
    (
        "self-edit, label read errors -> fail closed",
        {
            "OUT": "",
            "HAS_TOKEN": "true",
            "SELF_EDIT_FIXTURE": "yes",
            "GH_FAILS": "yes",
        },
        1,
    ),
    # --- a real verdict: untouched by the self-edit work, and proved so ---
    ("clean verdict -> pass", {"OUT": CLEAN, "HAS_TOKEN": "true"}, 0),
    ("UNSUPPORTED provenance -> fail", {"OUT": UNSUPPORTED, "HAS_TOKEN": "true"}, 1),
    ("blocking=true -> fail", {"OUT": BLOCKING, "HAS_TOKEN": "true"}, 1),
]


def extract_step() -> str:
    """The real `run:` body of the verdict step, straight out of the workflow."""
    import yaml

    doc = yaml.safe_load(WORKFLOW.read_text(encoding="utf-8"))
    steps = doc["jobs"]["claude-review-gate"]["steps"]
    for step in steps:
        if step.get("name") == STEP_NAME:
            return str(step["run"])
    raise SystemExit(
        f"::error::{WORKFLOW}: no step named {STEP_NAME!r}. Renaming it silently "
        f"turns this whole file into a no-op, so it is a failure, not a skip."
    )


def structural_checks(script: str) -> list[str]:
    problems = []
    if 'gh api "repos/$GH_REPO/pulls/$PR"' not in script:
        problems.append(
            "the label read must use the `pulls` endpoint. `issues/$PR/labels` "
            "returns the same labels but is mapped to the `issues` scope, which "
            "this workflow does not grant -- it would 403 into `set -e` and "
            "restore the permanent red the self-edit branch exists to remove."
        )
    if '"approved:ci-change"' not in script:
        problems.append(
            "the self-edit branch must match `approved:ci-change` as a COMPLETE "
            "JSON string (quotes included), or `approved:ci-changes` satisfies it."
        )
    if 'echo "::error::No verdict (crash, timeout, or --max-turns). Failing closed."' not in script:
        problems.append(
            "the fault branch must still fail closed. A crash, a timeout or a "
            "turn cap is a fault where a re-run can still deliver the verdict; "
            "2026-09-03-ci-checks-that-can-go-green.md rejected passing there, "
            "and the self-edit branch is a narrow exception, not a general one."
        )
    return problems


def main() -> int:
    if shutil.which("jq") is None:
        print("::error::jq not found; cannot exercise the verdict step.")
        return 2
    if not WORKFLOW.exists():
        print(f"::error::{WORKFLOW} not found (run from the repository root).")
        return 2

    script = extract_step()
    failures = 0

    for problem in structural_checks(script):
        print(f"::error::{problem}")
        failures += 1

    with tempfile.TemporaryDirectory() as tmp:
        stubs = Path(tmp) / "stubs"
        stubs.mkdir()
        for name, body in (("git", GIT_STUB), ("gh", GH_STUB)):
            path = stubs / name
            path.write_text(body, encoding="utf-8")
            path.chmod(0o755)
        script_path = Path(tmp) / "gate.sh"
        script_path.write_text(script, encoding="utf-8")

        for name, extra, expected in CASES:
            env = dict(os.environ)
            env["PATH"] = f"{stubs}{os.pathsep}{env['PATH']}"
            env.update({"BASE_SHA": "deadbeef", "GH_REPO": "o/r", "PR": "1", "GH_TOKEN": "x"})
            env.update(extra)
            proc = subprocess.run(
                ["bash", str(script_path)],
                env=env,
                capture_output=True,
                text=True,
                encoding="utf-8",
            )
            if proc.returncode == expected:
                print(f"  ok    {name}")
                continue
            failures += 1
            print(f"::error::{name}: exit {proc.returncode}, expected {expected}")
            for line in (proc.stdout + proc.stderr).splitlines():
                print(f"        {line}")

    if failures:
        print(f"::error::{failures} check(s) failed in the claude-review verdict step")
        return 1
    print(f"ok  {len(CASES)} verdict-step cases, 3 structural checks")
    return 0


if __name__ == "__main__":
    sys.exit(main())
