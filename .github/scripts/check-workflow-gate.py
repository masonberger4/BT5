#!/usr/bin/env python3
"""Guard the merge gate against the two ways it silently stops mattering.

docs/PLAN.md lists eight things that silently disable this system. Two of them
are properties of the workflow files themselves, and both fail SILENTLY -- no
error, no timeout, nothing red to click on:

  1. A `paths:` (or `paths-ignore:`) filter on a workflow that owns a required
     check. The workflow never triggers, so the check never reports, and the
     pull request shows "Expected -- Waiting for status to be reported" forever.

  2. A job that is not in the gate's `needs:`. It runs, it renders in the check
     list, a reviewer sees its red X -- and the merge box still goes green,
     because the gate only inspects `needs.*.result`. That is worse than not
     having the job: it looks like coverage and enforces nothing.

The required contexts are read from the RULESET SPEC rather than hard-coded, so
renaming the gate job without renaming the ruleset context (or the reverse) is
caught here instead of becoming a check that never reports.

This replaces a shell guard that scanned one hard-coded file with
`awk '/^on:/,/^jobs:/' | grep -qE '^\\s+paths:'`. That guard passed three
filters that deadlock identically: `paths-ignore:` (the form you reach for when
you want to skip CI on docs-only PRs, which is the usual reason anyone edits
this block at all), flow-style `pull_request: {paths: [...]}`, and a quoted
`"on":` key -- against which the awk range never opens, so ANY filter passed.
It also covered only ci.yml, so a second workflow got nothing.
"""

from __future__ import annotations

import json
import pathlib
import sys
from typing import Any

WORKFLOWS = pathlib.Path(".github/workflows")
RULESET = pathlib.Path(".github/rulesets/main-protection.json")

#: Jobs deliberately outside the merge gate. Empty on purpose: adding a name
#: here is a decision to let a job fail without blocking a merge, and it should
#: be visible in a diff rather than implied by an omission from `needs`.
NON_BLOCKING: frozenset[str] = frozenset()

PATH_FILTER_KEYS = ("paths", "paths-ignore")


def load_on(doc: dict[str, Any]) -> Any:
    """The `on:` block. PyYAML parses an unquoted `on:` key as the boolean True,
    which is why this cannot simply be `doc["on"]` -- and why a guard that reads
    the file as text has to cope with both spellings."""
    if "on" in doc:
        return doc["on"]
    return doc.get(True)


def required_contexts() -> set[str]:
    spec = json.loads(RULESET.read_text())
    return {
        check["context"]
        for rule in spec.get("rules", [])
        if rule.get("type") == "required_status_checks"
        for check in rule.get("parameters", {}).get("required_status_checks", [])
    }


def job_context(name: str, body: Any) -> str:
    """What GitHub calls this job in the checks list: its `name:` if set."""
    if isinstance(body, dict) and isinstance(body.get("name"), str):
        return body["name"]
    return name


def main() -> int:
    import yaml

    problems: list[str] = []
    contexts = required_contexts()
    if not contexts:
        problems.append(
            f"{RULESET}: no required_status_checks contexts found. Either the "
            f"ruleset stopped requiring a check, or this guard is reading the "
            f"wrong shape -- both mean nothing is gating a merge."
        )

    produced: dict[str, pathlib.Path] = {}
    for wf in sorted(WORKFLOWS.glob("*.y*ml")):
        doc = yaml.safe_load(wf.read_text()) or {}
        jobs = doc.get("jobs") or {}
        on = load_on(doc)

        if on is None:
            problems.append(f"{wf}: no `on:` block found; cannot verify its triggers")
        elif isinstance(on, dict):
            for event, cfg in on.items():
                if not isinstance(cfg, dict):
                    continue
                for key in PATH_FILTER_KEYS:
                    if key in cfg:
                        problems.append(
                            f"{wf}: on.{event}.{key} -- a workflow owning a required "
                            f"check must not filter on paths. The check never reports "
                            f"and the pull request blocks forever with no error. "
                            f"Filter inside the job instead."
                        )

        for name, body in jobs.items():
            produced[job_context(name, body)] = wf

        # Every job in a workflow that owns a gate must feed that gate.
        for name, body in jobs.items():
            if job_context(name, body) not in contexts:
                continue
            needs = (body or {}).get("needs") or []
            if isinstance(needs, str):
                needs = [needs]
            orphans = sorted(
                set(jobs) - set(needs) - {name} - NON_BLOCKING,
            )
            if orphans:
                problems.append(
                    f"{wf}: job(s) {', '.join(orphans)} are not in {name}.needs, so "
                    f"they can fail without blocking a merge. Add them to `needs:`, "
                    f"or to NON_BLOCKING in this script to say so deliberately."
                )

    for context in sorted(contexts):
        if context not in produced:
            problems.append(
                f"ruleset requires the check {context!r}, which no job in "
                f"{WORKFLOWS} produces. A required check that never reports blocks "
                f"every pull request permanently, with no error."
            )

    for problem in problems:
        print(f"::error::{problem}")
    if not problems:
        print(
            f"gate ok: contexts {sorted(contexts)} produced, "
            f"every job feeds the gate, no path filters"
        )
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
