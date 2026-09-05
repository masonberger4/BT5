---
name: ci-workflow-review
description: Patterns and gotchas seen reviewing .github/workflows/*.yml and check-workflow-gate.py in BT5
metadata:
  type: project
---

> **Note, 2026-09-05.** `pre-pr-attest.yml` and `claude-review.yml` were deleted
> (`docs/decisions/2026-09-05-remove-claude-usage-ci-checks.md`). Findings below that
> cite line numbers in those files are historical record, not live references.

Seen reviewing PR #104 (branch claude/github-ci-checks-2969e5, 2026-09-03), fixing the
`pre-pr-attest` re-arm no-op (issue_comment runs are default-branch runs, so a check
posted from one never lands on a PR head).

**Pattern worth recognizing again:** splitting one workflow into "reporter job on the
event whose check lands on the PR head" + "re-arm job on the event that can't report,
which instead calls `gh run rerun` on the reporter's run" is the correct fix for this
class of bug, not a red flag by itself. When reviewing a similar split, check:
- The `if:` conditions on the two jobs are mutually exclusive by `event_name` (never
  both true, never both false when a check is expected).
- `concurrency.group` includes `event_name` (or something equally distinguishing) —
  without it, an issue_comment run and the pull_request_target run it's meant to
  re-arm share one cancel-in-progress group and the comment cancels the very run
  it was posted to fix.
- The re-arm job selects the run to rerun by scoping on the CURRENT head_sha (fetched
  fresh via API, not from the event payload) and filters to `status == completed` —
  otherwise it can rerun a stale-head run or race against an in-flight rerun.
- A re-arm/no-report job's own failure (e.g. `gh run rerun` erroring) is invisible to
  the PR (it's a default-branch run) — that's by design for the check-gate, but it
  means operational failures of the rerun call are silently swallowed from the PR
  author's point of view. Flag as non-blocking, not blocking, given this repo's
  single-owner low-traffic operation.

**`check-workflow-gate.py`'s `NON_BLOCKING` set is a flat, workflow-unqualified set of
job names** (not `workflow.yml:jobname`). Adding a name there (e.g. `rearm`,
`main-broken`) exempts any job with that name in ANY workflow file from the
orphan-in-needs check, not just the one workflow that motivated the change. Pre-existing
design wart, not usually worth blocking on, but worth a one-line non-blocking note if a
generic name (like `rearm`) is added.

**Add up the numbers in a decision doc's Evidence section against the workflow comment it
paraphrases.** On PR #104 the workflow comment said six of the last fifteen runs died on
`--max-turns`, while the decision doc listed failures at "31, 31, 31, 32, 42" — five
entries for six failures, one `31` dropped. Real defect, corrected to
"31, 31, 31, 31, 32, 42".

Be careful *how* you argue it, though: my first write-up of this also claimed a mismatch
between "8 cited" and "fifteen runs", which was my own error — the turn-count list and the
run count measure different things, and a workflow comment listing only the runs that
*completed* is not inconsistent with a decision doc listing all *outcomes*. Check the
arithmetic of the objection before raising it. A miscounted objection to a miscount costs
more credibility than staying silent.

**A gate's own error message can be a valid input to the gate.** The single best find in
this area came from reading `pre-pr-attest.yml`'s matcher against the text the same job
prints on failure. The regex allowed `[ \t]*` before the command; the help text prints
`    /pre-pr <head-sha>` indented four spaces. So anyone counted as an attestor who pasted
the failing log into a comment to ask about it silently turned the check green with no
review having happened. Fixed by anchoring to the line start (PR #113).

Generalise it: **whenever a job both EMITS text and PARSES text, check whether its own
output satisfies its own parser.** Same question for log lines, bot comments, PR templates
and issue bodies. Nothing else in the review — not the original design, not a full pass by
this agent — caught it; it only appeared when the two halves of the file were read against
each other.

**Verify a claim about where a check run lands by looking at a real rollup**, not by
reasoning from the YAML. Two comments in `pre-pr-attest.yml` were confidently wrong about
this and survived multiple readings: "reports nothing at all" (every job emits a check run
named by its `name:`), and then the correction itself, which named one location when
`rearm` in fact reports in two — the PR head on `pull_request_target` where it lands
`skipped`, and the default-branch tip on `issue_comment` where it does the work.
`gh pr view <n> --json statusCheckRollup` settles it in one call.

**"The workflow triggers on the event" is not "the job runs on the event", and for a
required check the difference inverts the failure.** A workflow can declare
`on: {merge_group:}` while the job producing the required context is gated
`if: github.event_name == 'pull_request_target'`. The job is then `skipped` on
`merge_group` — and a skipped check **satisfies** a required status check. So the gate
passes while enforcing nothing, which is worse than the deadlock it was added to prevent,
because nothing goes red.

Raised on PR #114 as a non-blocking finding against the first cut of
`check-workflow-gate.py`'s `merge_group` guard, which only read the workflow-level `on:`
block. **Closed in `e87f2a4`** — the guard now inspects the producing job's `if:` too.
Check both halves when reviewing anything of this shape.

**And the trap inside that fix.** The obvious rule — "a gate job's `if:` must contain
`merge_group`" — flags `ci.yml`'s `required-checks`, whose `if: always()` runs on every
event including a merge group. Shipping that would have blocked all of CI immediately. The
guard therefore fires only when a condition *branches on the event*
(`github.event_name` / `github.event.`) without naming `merge_group`. It was the script's
own negative test that caught this, not a reading of it — when adding a heuristic guard,
write the negative test before trusting the positive result.

See also [[skill-frontmatter-governance]].
