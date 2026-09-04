#!/usr/bin/env python3
"""Fixture test for the `/pre-pr` attestation matcher in pre-pr-attest.yml.

`pre-pr-attest` decides, for every pull request, whether `/pre-pr <head-sha>`
was actually claimed against the commit being merged. Once it is a required
context that single regex is the last thing standing between an unreviewed
commit and `main`, and it has already shipped two distinct matching bugs:

  1. A trailing character class of `([ \t]|$)` rejected every attestation that
     was not the entire comment body -- so a real review, correctly claimed,
     did not count if the author wrote a sentence around it.
  2. A leading `[ \t]*` allowed indented commands, which made THIS JOB'S OWN
     FAILURE MESSAGE a valid attestation: the help text prints
     `    /pre-pr <head-sha>` indented four spaces, so a maintainer pasting the
     log into a comment to ask about it silently turned the check green.

Both were found by reading. Neither would have survived this file.

THE PROGRAM UNDER TEST IS EXTRACTED FROM THE WORKFLOW, NOT COPIED HERE. A test
that restates the regex tests the restatement, which is exactly how (2) shipped.
`jq` uses Oniguruma, where `^` and `$` are line-anchored by default and Python's
`re` are not, so re-implementing the match in Python would test a translation
with different semantics on precisely the newline handling that matters. The
fixtures are therefore driven through the real `jq` with the real program, and
the two `author_association` allowlists are read from the real call sites.

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
import tempfile
from pathlib import Path

WORKFLOW = Path(".github/workflows/pre-pr-attest.yml")

# A real 40-character SHA, and a decoy sharing its first seven characters --
# the case that makes "the SHA must be the full 40" load-bearing rather than
# decorative.
HEAD = "a273fe1836e1ada8482bf00f489b5a4246c48eb9"
DECOY = "a273fe1000000000000000000000000000000000"

JQ_CALL = re.compile(
    r"""jq\s+-r\s+--arg\s+sha\s+"\$HEAD_SHA"\s+--arg\s+cmd\s+"\$1"\s+"""
    r"""--argjson\s+who\s+"\$2"\s+'(?P<program>.*?)'\s+/tmp/comments\.json""",
    re.DOTALL,
)
CALL_SITE = re.compile(
    r"""(?P<var>ATTESTED|BYPASSED)="\$\(claims\s+'(?P<cmd>[^']+)'\s+"""
    r"""'(?P<who>\[[^']*\])'\)\"""",
)


def extract() -> tuple[str, dict[str, tuple[str, list[str]]]]:
    """Pull the jq program and both call sites out of the workflow.

    Deliberately tight to the shipped invocation: if someone rewrites how
    `claims` is called, this fails loudly instead of silently testing a
    program the workflow no longer runs.
    """
    text = WORKFLOW.read_text(encoding="utf-8")
    match = JQ_CALL.search(text)
    if match is None:
        raise SystemExit(
            f"::error::could not find the claims() jq program in {WORKFLOW}. "
            f"If the invocation changed, update this test deliberately -- do "
            f"not delete it."
        )
    sites: dict[str, tuple[str, list[str]]] = {}
    for site in CALL_SITE.finditer(text):
        sites[site["var"]] = (site["cmd"], json.loads(site["who"]))
    missing = {"ATTESTED", "BYPASSED"} - sites.keys()
    if missing:
        raise SystemExit(f"::error::call sites not found in {WORKFLOW}: {sorted(missing)}")
    return match["program"], sites


def comment(body: str, association: str = "OWNER") -> dict[str, str]:
    return {"author_association": association, "body": body}


def run(program: str, comments: list[dict[str, str]], cmd: str, who: list[str]) -> int:
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False, encoding="utf-8") as fh:
        json.dump(comments, fh)
        path = fh.name
    try:
        result = subprocess.run(
            [
                "jq",
                "-r",
                "--arg",
                "sha",
                HEAD,
                "--arg",
                "cmd",
                cmd,
                "--argjson",
                "who",
                json.dumps(who),
                program,
                path,
            ],
            capture_output=True,
            text=True,
        )
    finally:
        Path(path).unlink(missing_ok=True)
    if result.returncode != 0:
        raise SystemExit(f"::error::jq failed: {result.stderr.strip()}")
    return int(result.stdout.strip())


# Each case is (name, comments, expected_count) against the ATTESTED call site.
# Expected counts, not booleans: the job branches on `!= 0`, but pinning the
# count catches a regex that starts matching twice as readily as one that stops
# matching at all.
ATTEST_CASES: list[tuple[str, list[dict[str, str]], int]] = [
    # -- what must count -------------------------------------------------
    ("a bare attestation on its own line", [comment(f"/pre-pr {HEAD}")], 1),
    (
        "an attestation with prose around it -- bug (1), which rejected these",
        [comment(f"Chain ran, all green.\n/pre-pr {HEAD}\nMerging shortly.")],
        1,
    ),
    ("the first line of the body", [comment(f"/pre-pr {HEAD}\nnotes follow")], 1),
    ("the last line, with no trailing newline", [comment(f"notes\n/pre-pr {HEAD}")], 1),
    ("a tab between command and sha", [comment(f"/pre-pr\t{HEAD}")], 1),
    ("several spaces between command and sha", [comment(f"/pre-pr   {HEAD}")], 1),
    ("trailing whitespace after the sha", [comment(f"/pre-pr {HEAD}  ")], 1),
    (
        "CRLF line endings, as a browser-submitted comment sends",
        [comment(f"a\r\n/pre-pr {HEAD}\r\nb")],
        1,
    ),
    ("a COLLABORATOR may attest", [comment(f"/pre-pr {HEAD}", "COLLABORATOR")], 1),
    ("a MEMBER may attest", [comment(f"/pre-pr {HEAD}", "MEMBER")], 1),
    (
        "two separate attestations are both counted",
        [comment(f"/pre-pr {HEAD}"), comment(f"re-attesting\n/pre-pr {HEAD}")],
        2,
    ),
    # -- what must NOT count ---------------------------------------------
    (
        "an INDENTED command -- bug (2): this job's own help text, pasted",
        [comment(f"the check says:\n    /pre-pr {HEAD}\nwhy is it red?")],
        0,
    ),
    (
        "a quote-reply, which is how the help text usually gets pasted",
        [comment(f"> /pre-pr {HEAD}")],
        0,
    ),
    ("a short sha", [comment(f"/pre-pr {HEAD[:7]}")], 0),
    ("a different sha sharing the first seven characters", [comment(f"/pre-pr {DECOY}")], 0),
    ("a sha with an extra character appended", [comment(f"/pre-pr {HEAD}0")], 0),
    ("the command with no sha at all", [comment("/pre-pr")], 0),
    ("the sha with no command", [comment(f"reviewed {HEAD}, looks fine")], 0),
    ("/pre-pr-bypass must not satisfy /pre-pr", [comment(f"/pre-pr-bypass {HEAD}")], 0),
    (
        "a CONTRIBUTOR -- every fork pull request author",
        [comment(f"/pre-pr {HEAD}", "CONTRIBUTOR")],
        0,
    ),
    ("a drive-by commenter", [comment(f"/pre-pr {HEAD}", "NONE")], 0),
    ("no comments at all", [], 0),
    # -- a residual, pinned so it cannot change silently -----------------
    # This one documents current behaviour rather than endorsing it: a fenced
    # code block does NOT protect a pasted attestation, because a fence line
    # still starts at column 0. Anchoring fixed the indented and quoted pastes,
    # not this one. If a future change makes fences safe, this expectation
    # should flip deliberately -- with a comment -- not drift.
    (
        "KNOWN RESIDUAL: a fenced code block does not protect a paste",
        [comment(f"```\n/pre-pr {HEAD}\n```")],
        1,
    ),
]

# Against the BYPASSED call site, whose allowlist is OWNER alone.
BYPASS_CASES: list[tuple[str, list[dict[str, str]], int]] = [
    ("the owner may waive", [comment(f"/pre-pr-bypass {HEAD}")], 1),
    ("a COLLABORATOR may not waive", [comment(f"/pre-pr-bypass {HEAD}", "COLLABORATOR")], 0),
    ("a MEMBER may not waive", [comment(f"/pre-pr-bypass {HEAD}", "MEMBER")], 0),
    ("a plain attestation is not a waiver", [comment(f"/pre-pr {HEAD}")], 0),
    ("an indented waiver does not count either", [comment(f"    /pre-pr-bypass {HEAD}")], 0),
]


def main() -> int:
    if shutil.which("jq") is None:
        print("::error::jq not found. This test drives the real matcher through jq; ")
        print("::error::without it the matcher is UNTESTED. Not reporting a pass.")
        return 2
    if not WORKFLOW.exists():
        print(f"::error::{WORKFLOW} not found -- run from the repository root")
        return 2

    program, sites = extract()
    failures: list[str] = []
    checked = 0

    for label, (cmd, who) in (("ATTESTED", sites["ATTESTED"]), ("BYPASSED", sites["BYPASSED"])):
        cases = ATTEST_CASES if label == "ATTESTED" else BYPASS_CASES
        for name, comments, expected in cases:
            checked += 1
            got = run(program, comments, cmd, who)
            if got != expected:
                failures.append(f"{label}: {name} -- expected {expected}, got {got}")

    for failure in failures:
        print(f"::error::{failure}")
    if failures:
        print(f"::error::attestation matcher: {len(failures)} of {checked} fixtures failed")
        return 1
    print(
        f"attestation matcher ok: {checked} fixtures, both call sites, program read from {WORKFLOW}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
