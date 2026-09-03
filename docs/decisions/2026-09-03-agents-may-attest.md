## 2026-09-03 — agents may run /pre-pr and post their own attestation

**Decided:** `/pre-pr` loses `disable-model-invocation: true`. An agent may run the whole
chain — gates, `code-reviewer`, the conditional security and provenance passes, the label
determination, the draft PR — and may post the `/pre-pr <sha>` comment at step 10, under
the condition step 10 already carries: **not** if a gate or a review came back blocking.

This supersedes the "Also unresolved" item in `2026-09-01-pre-pr-as-ci-gates.md:56`, which
left the toggle open and accepted `pre-pr-attest` red on every agent-opened PR as its cost.

**What the check now proves, and what it stops proving.** It stops meaning *"a human
vouched for this commit"*. It now means *"the agent's chain ran against this commit and did
not come back blocking"*. That is weaker, and it is the point — the goal is automated CI
checks, and an attestation only a human can post is not automatable. What survives is the
only property the workflow header ever claimed: **skipping becomes visible**
(`pre-pr-attest.yml:26-30`, which already conceded "Claude posts the comment"). A red
`pre-pr-attest` still separates "no chain ran" and "the chain came back blocking" from
"clean". `/pre-pr-bypass` stays OWNER-only and stays a human act, so the *waiver* path is
untouched — which is what keeps the two commands worth distinguishing.

**Rejected:**

- *Model-invocable chain, human-only attestation — gate step 10 on an explicit ask.* This
  closes the integrity concern completely and was the recommendation put to the owner. It
  lost on the owner's call: it leaves the owner commenting once per head forever, which is
  precisely the friction the 2026-09-01 note recorded, and is not "automated".
- *Leave the flag in place.* Keeps the 2026-09-01 item open and keeps every agent-opened
  PR red on a check no agent can clear.
- *Let the agent attest unconditionally.* Rejected by keeping step 10's existing
  condition. An attestation posted regardless of the chain's verdict proves only that the
  chain was invoked, which is not worth a check run.

**Consequence that had to be fixed in the same slice.** `rearm` originally skipped a run
whose `status` was not yet `completed`, on the reasoning that an in-flight run reads the
comment list itself. That is safe only while a human is driving, because a human who sees
the check still red just comments again. Unattended there is no second comment, and the
attesting party pushes and comments seconds apart — so the in-flight window is the normal
case, not the edge case. `rearm` now polls for a terminal state (bounded at 180s, inside
its own `timeout-minutes: 5`) and re-runs only if the run landed non-green.

**Evidence:**

- The open item this closes: `2026-09-01-pre-pr-as-ci-gates.md:56` — *"`/pre-pr` carries
  `disable-model-invocation: true`, so an agent cannot run it and therefore cannot honestly
  attest. Until that is decided, agent-opened PRs show `pre-pr-attest` red unless the owner
  runs `/pre-pr` or comments `/pre-pr-bypass <sha>`."*
- `code-reviewer` returned this change as its single BLOCKING finding, arguing the model
  could attest off an ambiguous instruction. Recorded because it is the real cost of this
  decision. Two of its supports do not hold: it cited "CLAUDE.md §3" for the
  `/pre-pr` vs `/pre-pr-bypass` split (§3 is the non-negotiable correctness rules; the
  split is argued at `pre-pr-attest.yml:32-34`), and it read the header's "Claude posts the
  comment" as conditioned on a human having typed the command first, which the header does
  not say.
- `security-reviewer` could not complete: three attempts, two on opus and one on sonnet,
  all terminated by an API safeguard on `[bio]`. The agent's own definition is written
  around "biosecurity posture", so the trigger is the agent, not the diff. A substitute
  permissions audit ran on a generic agent instead, and it earned its keep — see below.

**A real hole the substitute audit found, fixed in this slice.** The matcher allowed
`[ \t]*` before the command, and the attest job's own failure message prints
`    /pre-pr <head-sha>` indented four spaces. So a counted author who pasted the log into
a comment to ask about it silently turned the check green with no review having happened —
the gate defeating its own single stated property. The matcher is now anchored to the line
start. Verified both ways: the real OWNER attestation on #101 still matches through jq
against the live API, and the pasted-help-text case flips `True` -> `False`, while
quote-replies and `/pre-pr-bypass` stay correctly rejected.

Three residuals from the same audit are now documented in the workflow header rather than
fixed, because each needs a design decision rather than a patch: `author_association` is
resolved at read time (so adding a collaborator retroactively promotes their old
comments); only `issue_comment: [created]` is subscribed while bodies are re-read live (so
a counted author can edit an old comment into an attestation, costing auditability rather
than access); and `author_association` is not a permission check at all — `COLLABORATOR`
includes read-only and triage. The last is latent, not live: the repository has exactly one
collaborator, the owner, at `admin`. Two claims the audit could not settle from the files
were closed from live data: `claude[bot]` resolves to `author_association: NONE`, so bot
comments can never be counted, and a `pull_request_target` check run does attach to the PR
head, proven when re-running 33784923222 flipped #101's rollup entry.
- `rearm`'s six control-flow branches driven with `gh`/`jq`/`sleep` stubbed:
  in-flight→success (no re-run), in-flight→failure (re-run), in-flight→cancelled (re-run),
  deadline exceeded (warn, no re-run), no run for the head (notice), already green (no
  re-run). Plus both live-API paths against real head SHAs.

**Where:** PR #104, branch `claude/github-ci-checks-2969e5`, carrying `approved:ci-change`.
Owner merges (CLAUDE.md §7b: protected path).

**Left open:** `claude-review-gate` is still red by construction on any PR that edits
`claude-review.yml` — flagged at `2026-09-01-pre-pr-as-ci-gates.md:52` as *"a required
check red by construction cannot be satisfied by any push"*. Harmless while advisory, and
the blocker on ever promoting it. Not addressed here.
