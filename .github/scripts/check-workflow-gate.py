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
import re
import sys
from typing import Any

WORKFLOWS = pathlib.Path(".github/workflows")
RULESET = pathlib.Path(".github/rulesets/main-protection.json")

#: Jobs deliberately outside the merge gate. Adding a name here is a decision to
#: let a job fail without blocking a merge, and it belongs in a diff rather than
#: being implied by an omission from `needs`.
#:
#: `main-broken` runs ONLY on push-to-main and is skipped on every pull_request
#: event. The gate counts `skipped` as failure, so putting it in `needs` would
#: block every pull request permanently -- the exact deadlock this script exists
#: to catch, arrived at from the other direction.
#:
#: `rearm` (pre-pr-attest.yml) runs ONLY on issue_comment and produces no check
#: run of its own -- it re-runs the pull_request_target run that owns the check,
#: because an issue_comment run's check lands on main's tip and not on the
#: pull request's head. Same skipped-is-failure deadlock, and listing it here is
#: what keeps promoting `pre-pr-attest` to a required context the single ruleset
#: line that file's header promises it is.
NON_BLOCKING: frozenset[str] = frozenset({"main-broken", "rearm"})

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

    # context -> (workflow file, job body). The BODY is needed too: a workflow
    # can trigger on an event while the job producing the context is gated off
    # it, and those two failures need different messages.
    produced: dict[str, tuple[pathlib.Path, Any]] = {}
    triggers: dict[pathlib.Path, set[str]] = {}
    for wf in sorted(WORKFLOWS.glob("*.y*ml")):
        doc = yaml.safe_load(wf.read_text()) or {}
        jobs = doc.get("jobs") or {}
        on = load_on(doc)

        # `on:` is a mapping in every workflow here, but the schema also permits
        # a bare string or a list, and a guard that only understands one shape
        # silently passes the others -- the failure this file already documents
        # for `paths:`.
        if isinstance(on, dict):
            triggers[wf] = set(on)
        elif isinstance(on, str):
            triggers[wf] = {on}
        elif isinstance(on, list):
            triggers[wf] = set(on)
        else:
            triggers[wf] = set()

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
            produced[job_context(name, body)] = (wf, body)

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
            continue

        # THIRD silent killer, and the one that arrives without touching a file.
        # A merge queue evaluates required checks against the MERGE-GROUP ref,
        # not the pull request head. A required context whose workflow does not
        # trigger on `merge_group` therefore never reports there, and the queue
        # blocks forever with no error -- the same deadlock as a `paths:` filter,
        # reached through a repository SETTING that no diff will ever show.
        #
        # Required unconditionally rather than only when some other workflow
        # declares it. Declaring `merge_group:` while no queue exists is inert;
        # discovering the gap the day a queue is switched on is not.
        wf, body = produced[context]
        if "merge_group" not in triggers.get(wf, set()):
            problems.append(
                f"{wf}: produces the required check {context!r} but does not "
                f"trigger on `merge_group`. If a merge queue is ever enabled, "
                f"this check never reports on the merge-group ref and the queue "
                f"blocks with no error. Add `merge_group:` to its `on:` block "
                f"and make the job report there (see ci.yml's `approvals` job, "
                f"which exits 0 off a pull request for the same reason)."
            )
            continue

        # ...and the HALF OF THAT which is worse, because it fails green.
        # Triggering the workflow on `merge_group` is not enough: if the job
        # producing the context carries an `if:` that excludes the event, the job
        # lands `skipped` -- and a skipped check SATISFIES a required status
        # check. The queue then merges with this gate enforcing nothing, which is
        # the "looks like coverage, enforces nothing" failure this file opens by
        # rejecting. Blocking is loud; a silent pass is not.
        #
        # A textual check, not an evaluation -- GitHub expressions cannot be
        # evaluated here. The rule: a condition is suspect if it BRANCHES ON THE
        # EVENT (mentions `github.event_name` or `github.event.`) without
        # POSITIVELY admitting `merge_group`. A condition that does not
        # discriminate on the event cannot exclude one, so it passes.
        #
        # `if: always()` is why this is not simply "must contain merge_group".
        # ci.yml's `required-checks` carries exactly that, runs on every event
        # including a merge group, and a naive check flagged it -- which would
        # have blocked all of CI the moment this guard shipped. Caught by this
        # script's own negative test, not in review.
        #
        # And it demands `== 'merge_group'` rather than the mere substring,
        # because `!= 'merge_group'` CONTAINS the substring while meaning the
        # exact opposite: skipped on precisely the event this check exists to
        # guarantee coverage for. A substring test called that safe.
        #
        # Known limitation, stated rather than hidden: a positive form that is
        # not an equality -- say
        # `contains(fromJSON('["pull_request_target","merge_group"]'), github.event_name)`
        # -- is flagged even though it is correct. That is the safe direction to
        # be wrong in: the fix is one edit to the `if:`, whereas a false pass is
        # a merge gate that enforces nothing and never goes red.
        condition = str((body or {}).get("if") or "")
        branches_on_event = "github.event_name" in condition or "github.event." in condition
        admits_merge_group = re.search(r"==\s*['\"]merge_group['\"]", condition) is not None
        if branches_on_event and not admits_merge_group:
            problems.append(
                f"{wf}: triggers on `merge_group`, but the job producing the "
                f"required check {context!r} has an `if:` that does not positively "
                f"admit it (no `== 'merge_group'`; note a `!= 'merge_group'` names "
                f"the event while excluding it), so the job is skipped on that "
                f"event -- and a skipped check SATISFIES a required status check. "
                f"The queue would merge with this gate enforcing nothing. Add "
                f"`github.event_name == 'merge_group'` to the `if:`, or drop the "
                f"`if:` so the job always runs."
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
