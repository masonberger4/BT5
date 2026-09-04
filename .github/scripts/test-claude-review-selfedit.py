#!/usr/bin/env python3
"""Fixture test for claude-review-gate's verdict step in claude-review.yml.

That step decides whether a pull request passes the review gate, and as of
2026-09-04 it contains THE ONE PATH IN THIS REPOSITORY'S CI THAT PASSES WITHOUT
A VERDICT: a pull request editing `claude-review.yml` cannot be reviewed at all,
because the Claude GitHub App refuses to mint a token for a workflow file
differing from the default branch's copy. That case used to be red forever --
unclearable by any push -- and now passes if, and only if, `approved:ci-change`
is on the pull request.

A branch that turns "no review happened" into a green check is exactly the kind
this repository has agreed not to ship untested. `test-attestation-matcher.py`
exists because its regex shipped two matching bugs "caught by reading, which
does not scale"; the review on #142 asked for the same treatment for the
`main-broken` predicate and got it. This is the third, and it guards the
loosest of the three.

THE SCRIPT UNDER TEST IS EXTRACTED FROM THE WORKFLOW, NOT COPIED HERE. A test
restating the branching would test the restatement -- which is how the
attestation matcher's second bug shipped.

GIT IS REAL HERE, NOT STUBBED, and that is the point of the file. The step
decides "did this pull request edit claude-review.yml" from the commit graph,
and the first version of it got that wrong in a way no stub could have shown:
it diffed the event payload's `base.sha` against `refs/pull/N/merge`, which
attributes anything `main` gained after the event TO THE PULL REQUEST. A
security audit reproduced it in a scratch repository; so these fixtures build
the same topologies with real commits and let real `git` answer. Only `gh` is
stubbed, because the label read is a network call with nothing to learn from.

Two properties are asserted structurally, because no fixture can observe them:

  * the diff must be anchored to `HEAD^1`/`HEAD^2`, never to `$BASE_SHA`.
  * the fault branch must still fail closed. `2026-09-03-ci-checks-that-can-go-
    green.md` rejected making this job neutral when it cannot produce a verdict;
    the self-edit branch is a narrow exception to that and must not become a
    general one.

Run by ci.yml's `python-quality` job, which is in `required-checks.needs`.
Requires `jq` and `git`; exits 2 rather than reporting a pass it did not perform.
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
TRACKED = ".github/workflows/claude-review.yml"

CLEAN = '{"blocking":false,"summary":"ok","blocking_findings":[],"provenance":"SUPPORTED"}'
UNSUPPORTED = '{"blocking":false,"summary":"bad","blocking_findings":[],"provenance":"UNSUPPORTED"}'
BLOCKING = '{"blocking":true,"summary":"x","blocking_findings":["y"],"provenance":"N/A"}'

GH_STUB = """#!/usr/bin/env bash
[ "${GH_FAILS:-no}" = "yes" ] && { echo "gh: simulated API error" >&2; exit 1; }
echo "${LABELS_FIXTURE:-[]}"
"""

GIT_ID = [
    "-c",
    "user.email=t@example.invalid",
    "-c",
    "user.name=t",
    "-c",
    "commit.gpgsign=false",
]


def git(repo: Path, *args: str) -> str:
    proc = subprocess.run(
        ["git", *GIT_ID, *args],
        cwd=repo,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    if proc.returncode != 0:
        raise SystemExit(f"::error::git {' '.join(args)} failed: {proc.stderr.strip()}")
    return proc.stdout.strip()


def write(repo: Path, rel: str, text: str) -> None:
    path = repo / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def build_repo(root: Path, *, main_edits: bool, pr_edits: bool, merge_ref: bool) -> Path:
    """A repository shaped like what `actions/checkout` leaves behind.

    On a `pull_request` event HEAD is `refs/pull/N/merge`: a merge commit whose
    FIRST parent is the base tip GitHub merged against and whose second is the
    pull request head. `merge_ref=False` builds a plain checkout instead, which
    has no second parent -- the case the step must refuse to judge.
    """
    # A multi-line body so `main` and the branch can BOTH edit this file without
    # colliding -- they touch opposite ends. That combination has to be
    # reachable: GitHub only publishes `refs/pull/N/merge` for a MERGEABLE pull
    # request, so a fixture whose two sides conflict could never occur in the
    # wild and would test nothing.
    original = [f"line {n}\n" for n in range(1, 11)]

    repo = root
    repo.mkdir(parents=True, exist_ok=True)
    git(repo, "init", "-q", "-b", "main")
    write(repo, TRACKED, "".join(original))
    write(repo, "other.txt", "one\n")
    git(repo, "add", "-A")
    git(repo, "commit", "-qm", "M1")
    base = git(repo, "rev-parse", "HEAD")

    git(repo, "checkout", "-q", "-b", "feature")
    write(repo, "other.txt", "two\n")
    if pr_edits:
        body = list(original)
        body[0] = "line 1, edited by the pull request\n"
        write(repo, TRACKED, "".join(body))
    git(repo, "add", "-A")
    git(repo, "commit", "-qm", "feature work")

    git(repo, "checkout", "-q", "main")
    git(repo, "reset", "-q", "--hard", base)
    if main_edits:
        # M2: `main` moves on AFTER the event fired, touching the same file at
        # the far end so the merge still succeeds.
        body = list(original)
        body[-1] = "line 10, edited on main, nothing to do with the PR\n"
        write(repo, TRACKED, "".join(body))
        git(repo, "add", "-A")
        git(repo, "commit", "-qm", "M2")

    if not merge_ref:
        git(repo, "checkout", "-q", "feature")
        return repo

    git(repo, "merge", "-q", "--no-ff", "-m", "merge", "feature")
    return repo


#: (name, repo spec or None, env, expected exit).
#: repo spec is (main_edits, pr_edits, merge_ref).
CASES: list[tuple[str, tuple[bool, bool, bool] | None, dict[str, str], int]] = [
    # --- empty verdict: the three causes, which must NOT collapse together ---
    (
        "no token configured -> fail closed",
        (False, False, True),
        {"OUT": "", "HAS_TOKEN": "false"},
        1,
    ),
    (
        "crash/timeout/max-turns -> fail closed",
        (False, False, True),
        {"OUT": "", "HAS_TOKEN": "true"},
        1,
    ),
    # --- THE REGRESSION THE AUDIT FOUND. `main` gained a commit touching this
    #     file after the event fired; the pull request never went near it. An
    #     anchor on the payload base reports SELF_EDIT=yes here and lets a
    #     FAULT-class empty verdict exit 0. It must fail closed.
    (
        "main edited the file, PR did not -> fail closed",
        (True, False, True),
        {"OUT": "", "HAS_TOKEN": "true", "LABELS_FIXTURE": '["approved:ci-change"]'},
        1,
    ),
    # --- the self-edit branch, which is the reason this file exists ---
    (
        "self-edit WITH approved:ci-change -> pass",
        (False, True, True),
        {"OUT": "", "HAS_TOKEN": "true", "LABELS_FIXTURE": '["approved:ci-change"]'},
        0,
    ),
    (
        "self-edit while main also moved -> pass",
        (True, True, True),
        {"OUT": "", "HAS_TOKEN": "true", "LABELS_FIXTURE": '["approved:ci-change"]'},
        0,
    ),
    (
        "self-edit, no labels at all -> fail",
        (False, True, True),
        {"OUT": "", "HAS_TOKEN": "true", "LABELS_FIXTURE": "[]"},
        1,
    ),
    (
        "self-edit, unrelated label -> fail",
        (False, True, True),
        {"OUT": "", "HAS_TOKEN": "true", "LABELS_FIXTURE": '["bug","approved:data-change"]'},
        1,
    ),
    # A label CONTAINING the required one as a prefix must not satisfy it -- the
    # trap check-approval-labels.sh's `has_label` matches whole JSON strings to
    # avoid. A substring test would pass this and be wrong.
    (
        "self-edit, prefix-only near-miss label -> fail",
        (False, True, True),
        {"OUT": "", "HAS_TOKEN": "true", "LABELS_FIXTURE": '["approved:ci-changes"]'},
        1,
    ),
    (
        "self-edit, label read errors -> fail closed",
        (False, True, True),
        {"OUT": "", "HAS_TOKEN": "true", "GH_FAILS": "yes"},
        1,
    ),
    # No second parent: not a merge ref, so the step cannot tell what the pull
    # request changed and must not guess.
    (
        "not a merge ref -> fail closed even though the file differs",
        (False, True, False),
        {"OUT": "", "HAS_TOKEN": "true", "LABELS_FIXTURE": '["approved:ci-change"]'},
        1,
    ),
    # --- a real verdict: untouched by the self-edit work, and proved so ---
    ("clean verdict -> pass", None, {"OUT": CLEAN, "HAS_TOKEN": "true"}, 0),
    ("UNSUPPORTED provenance -> fail", None, {"OUT": UNSUPPORTED, "HAS_TOKEN": "true"}, 1),
    ("blocking=true -> fail", None, {"OUT": BLOCKING, "HAS_TOKEN": "true"}, 1),
]


def extract_step() -> str:
    """The real `run:` body of the verdict step, straight out of the workflow."""
    import yaml

    doc = yaml.safe_load(WORKFLOW.read_text(encoding="utf-8"))
    for step in doc["jobs"]["claude-review-gate"]["steps"]:
        if step.get("name") == STEP_NAME:
            return str(step["run"])
    raise SystemExit(
        f"::error::{WORKFLOW}: no step named {STEP_NAME!r}. Renaming it silently "
        f"turns this whole file into a no-op, so it is a failure, not a skip."
    )


def structural_checks(script: str) -> list[str]:
    problems = []
    # Comments discuss $BASE_SHA at length -- explaining why it is NOT used is
    # the point of them. Strip whole-line comments before looking for it, so the
    # check reads the code and not the prose about the code.
    code = "\n".join(line for line in script.splitlines() if not line.lstrip().startswith("#"))
    if "BASE_SHA" in code:
        problems.append(
            "the self-edit diff must not be anchored to $BASE_SHA. That is the "
            "event payload's base, while HEAD is the merge ref -- so anything "
            "`main` gains after the event is attributed to this pull request, "
            "and a fault-class empty verdict takes the self-edit branch."
        )
    if 'git rev-parse --verify --quiet "HEAD^2"' not in script:
        problems.append(
            "the merge-ref guard must test HEAD^2. Every commit but the root has "
            "a first parent, so a HEAD^1 test verifies on a plain checkout and "
            "compares HEAD against its predecessor instead of against the base."
        )
    if 'git diff --quiet "HEAD^1" HEAD' not in script:
        problems.append(
            "the self-edit diff must run HEAD^1..HEAD -- on a merge ref that is "
            "exactly the base tip GitHub merged against."
        )
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
    for tool in ("jq", "git"):
        if shutil.which(tool) is None:
            print(f"::error::{tool} not found; cannot exercise the verdict step.")
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
        tmpdir = Path(tmp)
        stubs = tmpdir / "stubs"
        stubs.mkdir()
        gh = stubs / "gh"
        gh.write_text(GH_STUB, encoding="utf-8")
        gh.chmod(0o755)
        script_path = tmpdir / "gate.sh"
        script_path.write_text(script, encoding="utf-8")

        for i, (name, spec, extra, expected) in enumerate(CASES):
            cwd = tmpdir / f"repo{i}"
            if spec is None:
                cwd.mkdir()
            else:
                main_edits, pr_edits, merge_ref = spec
                build_repo(cwd, main_edits=main_edits, pr_edits=pr_edits, merge_ref=merge_ref)

            env = dict(os.environ)
            env["PATH"] = f"{stubs}{os.pathsep}{env['PATH']}"
            env.update({"GH_REPO": "o/r", "PR": "1", "GH_TOKEN": "x"})
            env.update(extra)
            proc = subprocess.run(
                ["bash", str(script_path)],
                cwd=cwd,
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
    print(f"ok  {len(CASES)} verdict-step cases, 6 structural checks")
    return 0


if __name__ == "__main__":
    sys.exit(main())
