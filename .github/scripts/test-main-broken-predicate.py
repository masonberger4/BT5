#!/usr/bin/env python3
"""Fixture test for `main-broken`'s superseded-run suppression in ci.yml.

`main-broken` is the safety net for a break that has already been pushed: a
required status check gates pull requests and can do nothing about `main`
itself. Its `if:` counts a `cancelled` job as broken, which is right in
principle -- a cancellation can hide a real break -- and wrong in practice for
the dominant cause, `cancel-in-progress` killing the still-running checks for
commit N when commit N+1 lands on top.

Six merges in ten minutes filed FIVE false issues that way -- #136, #137,
#138, #140, #141 -- while `main` was green throughout. Those five payloads are
the fixtures below. An alarm that cries wolf five times in ten minutes is worse
than no alarm, because the next real break arrives in the noise.

So the job now suppresses filing only when BOTH hold: no job reported
`failure`, and this commit is no longer the branch tip. The asymmetry is
deliberate and is what these fixtures pin:

    a real `failure`            -> files, whatever the tip is
    a cancellation AT the tip   -> files (nothing later will cover it)
    a cancellation, not the tip -> suppressed (superseded)

THE PROGRAM UNDER TEST IS EXTRACTED FROM THE WORKFLOW, NOT COPIED HERE: a test
that restates the predicate tests the restatement, and goes on passing while the
workflow's real copy drifts away from it. That is how the sibling attestation
matcher shipped two matching bugs before it was retired. The jq program is pulled
out of ci.yml and run through the real `jq`, and the surrounding shell guard is
checked structurally in the same file it lives in.

Run by ci.yml's `python-quality` job, which is in `required-checks.needs`.
Requires `jq` (preinstalled on GitHub's ubuntu runners); exits 2 if it is
absent rather than reporting a pass it did not perform.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

WORKFLOW = Path(".github/workflows/ci.yml")

JQ_CALL = re.compile(r"""jq\s+-e\s+'(?P<program>[^']+)'""")

#: The five real `toJSON(needs)` payloads, reduced to the field the predicate
#: reads. Every one of these filed an issue against a green `main`.
SUPERSEDED = {
    "#136 2037f29": {
        "changes": "success",
        "python-quality": "success",
        "mypy": "success",
        "invariants": "success",
        "contract": "success",
        "contract-freeze": "success",
        "python-tests": "cancelled",
    },
    "#137 3c9bf8b": {
        "changes": "success",
        "python-quality": "success",
        "mypy": "cancelled",
        "invariants": "success",
        "contract": "success",
        "contract-freeze": "success",
        "python-tests": "cancelled",
    },
    "#138 ad3f9f0": {
        "changes": "success",
        "python-quality": "success",
        "mypy": "success",
        "invariants": "success",
        "contract": "success",
        "contract-freeze": "success",
        "python-tests": "cancelled",
    },
    "#140 da01b19": {
        "changes": "success",
        "python-quality": "success",
        "mypy": "success",
        "invariants": "cancelled",
        "contract": "cancelled",
        "contract-freeze": "success",
        "python-tests": "cancelled",
    },
    "#141 b3ef0b2": {
        "changes": "success",
        "python-quality": "success",
        "mypy": "success",
        "invariants": "success",
        "contract": "success",
        "contract-freeze": "success",
        "python-tests": "cancelled",
    },
}

#: Must ALWAYS file, tip or not. A `failure` is a real break by definition.
MUST_FILE = {
    "a single failing test job": {"changes": "success", "python-tests": "failure"},
    "a failure among cancellations": {"mypy": "cancelled", "python-tests": "failure"},
    "everything failed": {"mypy": "failure", "python-tests": "failure"},
}


def fail(msg: str) -> None:
    print(f"FAIL: {msg}")
    sys.exit(1)


def extract_program() -> str:
    """The real jq predicate, out of the real workflow."""
    source = WORKFLOW.read_text(encoding="utf-8")
    matches = JQ_CALL.findall(source)
    if len(matches) != 1:
        fail(
            f"expected exactly one `jq -e '...'` call in {WORKFLOW}, found "
            f"{len(matches)}. If a second one was added, teach this test which "
            f"is the main-broken predicate rather than deleting the assertion."
        )
    return matches[0]


def any_failure(program: str, results: dict[str, str]) -> bool:
    """Run the extracted program exactly as the workflow does."""
    payload = json.dumps({job: {"result": r, "outputs": {}} for job, r in results.items()})
    proc = subprocess.run(
        ["jq", "-e", program],
        input=payload,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    if proc.returncode not in (0, 1):
        fail(f"jq errored on {payload}: {proc.stderr.strip()}")
    return proc.returncode == 0


def check_shell_guard() -> None:
    """The half jq cannot answer: that suppression needs BOTH conditions, and
    that the tip comparison uses the FULL sha.

    The full-sha point is not hypothetical. This job builds its issue title from
    `${SHA:0:9}`, so a nine-character comparison is the natural typo, and it
    would make two commits sharing a prefix indistinguishable.
    """
    source = WORKFLOW.read_text(encoding="utf-8")
    block = source[source.index("main-broken:") :]

    if '"$TIP" != "$SHA"' not in block:
        fail(
            'the tip comparison is no longer `"$TIP" != "$SHA"`. It must compare '
            "FULL shas -- the issue title uses ${SHA:0:9} and a truncated "
            "comparison would confuse two commits sharing a prefix."
        )
    if '-n "${TIP:-}"' not in block:
        fail(
            "the guard no longer requires TIP to be non-empty. A failed `gh api` "
            "call must fall through to FILING the issue: an alarm guarding main "
            "has to fail loud, not silent."
        )
    condition = block[: block.index("needs:")]
    for word in ("failure", "cancelled"):
        if f"'{word}'" not in condition:
            fail(
                f"main-broken's `if:` no longer mentions {word!r}. `cancelled` "
                f"must stay: a cancellation CAN hide a real break, and the "
                f"suppression is the narrow fix, not dropping the clause."
            )


def main() -> int:
    if not shutil.which("jq"):
        print("jq is not installed; cannot verify the predicate.")
        return 2
    if not WORKFLOW.exists():
        fail(f"{WORKFLOW} not found; run from the repository root.")

    program = extract_program()
    print(f"predicate under test: {program}")

    for name, results in SUPERSEDED.items():
        if any_failure(program, results):
            fail(
                f"{name}: predicate reports a failure, but every job in that real "
                f"payload was success or cancelled. It would file an issue for a "
                f"superseded run -- the exact false alarm this guard removes."
            )

    for name, results in MUST_FILE.items():
        if not any_failure(program, results):
            fail(
                f"{name}: predicate reports no failure on a payload containing "
                f"`failure`. A real break would be suppressed whenever the commit "
                f"is no longer the tip -- silently, which is the worst outcome."
            )

    check_shell_guard()

    print(
        f"ok  main-broken predicate: {len(SUPERSEDED)} superseded payloads suppressed, "
        f"{len(MUST_FILE)} real failures still filed, shell guard intact"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
