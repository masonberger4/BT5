---
name: ci-workflow-review
description: Patterns and gotchas seen reviewing .github/workflows/*.yml and check-workflow-gate.py in BT5
metadata:
  type: project
---

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
- A re-arm/no-report job's own failure (e.g. `gh run rerun` erroring) is invisible to the
  PR (default-branch run) — flag as non-blocking, not blocking, given this repo's
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

Check the arithmetic of the objection itself before raising it — a prior draft of this
finding claimed a second mismatch ("8 cited" vs "fifteen runs") that was the reviewer's own
error, conflating two different counts. A miscounted objection to a miscount costs more
credibility than staying silent.

**A gate's own error message can be a valid input to the gate.** `pre-pr-attest.yml`'s
attestation regex allowed `[ \t]*` before the command, and the job's own failure text prints
`    /pre-pr <head-sha>` indented four spaces — so pasting the failing log into a comment
silently turned the check green with no review done. Fixed by anchoring to line start
(PR #113). Generalise: **whenever a job both EMITS text and PARSES text, check whether its
own output satisfies its own parser** (log lines, bot comments, PR templates, issue bodies).

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
`if: github.event_name == 'pull_request_target'` — the job then `skipped`s on `merge_group`,
and a skipped check **satisfies** a required status check, so the gate passes while
enforcing nothing (worse than the deadlock it was added to prevent).

Raised on PR #114 against the first cut of `check-workflow-gate.py`'s `merge_group` guard
(only read the workflow-level `on:`). **Closed in `e87f2a4`** — now inspects the producing
job's `if:` too. Check both halves when reviewing anything of this shape.

**And the trap inside that fix.** The obvious rule — "a gate job's `if:` must contain
`merge_group`" — flags `ci.yml`'s `required-checks`, whose `if: always()` runs on every
event including a merge group. Shipping that would have blocked all of CI immediately. The
guard therefore fires only when a condition *branches on the event*
(`github.event_name` / `github.event.`) without naming `merge_group`. It was the script's
own negative test that caught this, not a reading of it — when adding a heuristic guard,
write the negative test before trusting the positive result.

**Merge-ref self-edit/drift detection: anchor to `HEAD^1`/`HEAD^2`, never a payload SHA.**
PR #145 (`03a2008`→`c46d1b2`) fixed `git diff --quiet "$BASE_SHA" HEAD` (`base.sha` frozen
at event time) to `git diff --quiet "HEAD^1" HEAD`, guarded by
`git rev-parse --verify --quiet "HEAD^2"`. On the default `pull_request` checkout
(`refs/pull/N/merge`, **needs `fetch-depth: 0`** — a shallow clone lacks `HEAD^1`'s tree),
the merge ref's first parent is the base tip as of merge-ref computation, second is the PR
head; testing `HEAD^2` existence is how you detect "this is actually a merge ref" (plain
push/checkout has none). Diffing a frozen payload SHA against a freshly-recomputed merge
ref silently attributes anything `main` gained in between to the PR. Apply the same
scrutiny to any future "did this PR touch file X" check keyed off
`github.event.pull_request.base.sha`.

**`bypass_actors` on `main-protection.json` is a closed question.** Rejected twice: once by
`docs/research/github-setup.md:244,246,2828,2903` (the repo's own design doc, unread before
the first attempt), once in review (PR #145, `03a2008`, reverted in `c46d1b2`). Agents in
this repo authenticate as the repo owner (GitHub cannot tell an agent's token from the owner
at a keyboard), so any bypass on the admin role is held by every session. Treat a future
bypass-actor addition as blocking on sight.

See also [[skill-frontmatter-governance]].
