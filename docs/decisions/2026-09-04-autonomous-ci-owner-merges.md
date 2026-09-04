## 2026-09-04 — every check is agent-clearable; the merge button is the owner's

**Decided:** three changes, from one owner instruction: *"Modify CI checks on GitHub so
that they can be fully done autonomously. I only want merge permission left to me, and I
can still delegate that to an agent if I so choose."*

1. **CLAUDE.md §7b loses the standing agent-merge permission.** It read *"An agent may
   squash-merge its own PR once CI is green"* with four carve-outs routed to the owner. It
   now reads: an agent takes a PR to green and stops, reporting readiness. The owner
   merges, and delegates per PR when they want to ("merge #N"). This is a **tightening**,
   and the only one here — everything else in this slice removes a human step.
2. **`main-protection.json` gains a bypass actor**: the Repository admin role,
   `bypass_mode: pull_request`. The owner can now merge past a red required check.
3. **`claude-review-gate` stops being red by construction** on a pull request that edits
   `claude-review.yml`. That branch now passes on a **checked** compensating control —
   `approved:ci-change` must be on the PR — instead of failing closed forever.

**Why 2 is the substance and not a loosening.** With `bypass_actors: []` the owner did not
in fact hold merge permission: they held it conditionally on GitHub Actions being healthy.
`pre-pr-attest.yml`'s own header spelled the trap out — a defect in that job fails every
open pull request *including the one that would repair it*, "no one, admin included, can
merge past it", and the only escape was to disable the ruleset out of band. That escape is
unreachable from a Claude session at all: the proxy blocks every repo-scoped GitHub REST
path (`scripts/apply-repo-settings.sh:5-7`), and the GitHub MCP server exposes no ruleset
tool. So the instruction "leave merge permission to me" could not be satisfied without
this change. `bypass_mode: pull_request` rather than `always` keeps direct pushes to
`main` blocked; the bypass reaches the merge box only.

**What was already autonomous, checked rather than assumed.** The audit found far less
broken than expected — three human-only steps, not a system of them. Against the live
check state of the five open pull requests (#131, #139, #142, #143, #144), every job in
`required-checks` is an ordinary code gate; `approvals` needs an `approved:*` label and
all five labels exist and are agent-appliable; `pre-pr-attest` was green on every open PR,
because agents comment as the owner and `2026-09-03-agents-may-attest.md` authorised
exactly that; and review threads are already being resolved by agents (#142, two threads,
both resolved). Nothing in that list needed changing, and nothing in it was changed.

**Rejected:**

- *A `/merge <head-sha>` comment command, OWNER-only, verifying the SHA is current and the
  required checks green before squashing.* Put to the owner as the recommended option
  against a plain owner-only policy; **the owner chose the plain one**. It would have made
  delegation work from a phone and put each delegation on the timeline, at the cost of a
  workflow holding `contents: write` on an `issue_comment` trigger — a merge button
  reachable by comment, in a repository that already declines to let Actions approve pull
  requests for the same class of reason.
- *Leaving `bypass_actors: []` and writing a runbook for the wedge instead.* The runbook's
  first step is unreachable from the environment the agents run in, so it is a runbook only
  the owner can execute, on a machine they may not be at, for a repository that is
  otherwise unblocked. Documenting an escape is not having one.
- *Making `claude-review-gate` neutral whenever it cannot produce a verdict.* Already
  rejected in `2026-09-03-ci-checks-that-can-go-green.md` — "a review that did not happen
  must not read as a pass" — and that rejection **stands**. The narrow branch taken here is
  a different shape: a crash, a timeout or a turn cap is a fault where a re-run can still
  deliver the verdict, so passing would lose a review that was owed; a self-edit is
  deterministic, and no push to that branch can ever produce one. The `else` branch is
  untouched and still fails closed.
- *Passing the self-edit case unconditionally.* The label is checked live, from the
  `pulls` endpoint, and its absence still fails. Read from `pulls` and not
  `issues/$PR/labels`: the latter maps to the `issues` scope, which the workflow does not
  grant, so it would 403 into `set -e` and land back on the permanent red being removed.
- *A `workflow_dispatch` workflow applying the ruleset with an admin PAT*, to close the
  last human step. It would be a workflow able to rewrite the rules that gate it, reachable
  by anyone who can edit `.github/` — the exposure `apply-repo-settings.sh` already reasons
  about when it sets `can_approve_pull_request_reviews=false`. The owner running a script
  occasionally is the cheaper side of that trade.

**Evidence, and what could not be verified from here.**

- `bypass_actors: []` blocking the owner is not inferred: rulesets grant admins no implicit
  override, and `pre-pr-attest.yml:18-22` states the consequence in the repository's own
  words.
- `actor_id: 5` for the built-in Repository admin role **could not be tested from this
  session** — `gh` is absent, and both `GET /repos/masonberger4/BT5` and
  `.../rulesets` return 403 through the proxy, authenticated or not. So
  `apply-repo-settings.sh` now **reads the stored bypass back** after the PUT and prints
  the resolved actor, warning loudly if the API stored none. A wrong id would not error —
  it would silently grant the bypass to a *different* role — which is the same failure mode
  the script already guards against for `integration_id`, handled the same way.
- The self-edit red is real and current: `claude-review-gate` on #143 failed twice
  (`HAS_TOKEN: true`, `OUT:` empty). That specific pair was a no-verdict fault rather than a
  self-edit, and is re-runnable; the self-edit case is the one recorded as unclearable in
  `2026-09-03-agents-may-attest.md`'s "Left open".

**Where:** branch `claude/autonomous-github-ci-boxw2y`, carrying `approved:ci-change`.
Owner merges (CLAUDE.md §2 protected path — and now §7b, which is what this PR changes).

**Left open:**

- **The ruleset spec is inert until applied.** `main-protection.json` is a spec, not the
  live config; the bypass does not exist until the owner runs
  `scripts/apply-repo-settings.sh` on an authenticated machine. Everything else in this
  slice takes effect on merge.
- `claude-review-gate` produced no verdict twice on #143, a README-only change. Advisory,
  and an agent can re-run it, so it blocks nothing — but the cause was not chased here.
- The three `author_association` residuals from `2026-09-03-agents-may-attest.md` are
  untouched and still open.
