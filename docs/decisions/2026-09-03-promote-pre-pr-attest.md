## 2026-09-03 — `pre-pr-attest` becomes a required context, with its two worst edges filed off

**Decided:** add `{ "context": "pre-pr-attest", "integration_id": 15368 }` to
`main-protection.json`, making it a second required context alongside `required-checks`.
Owner's call, taken with the risks below on the record rather than discovered later.

Four changes ship in the same diff because promotion is what makes each of them matter:

1. **`merge_group` coverage.** `ci.yml` declared `merge_group:`; this workflow did not. A
   merge queue evaluates required checks against the merge-group ref, so a required
   `pre-pr-attest` would never have reported there and the queue would block with no error
   — CLAUDE.md §9 reached through a repository *setting*, which no diff would ever show.
   The job now runs on `merge_group` and reports immediately, mirroring `ci.yml`'s
   `approvals` job. Not a hole: a pull request cannot enter the queue until the required
   checks pass on its own head, so everything in a merge group was already attested at the
   only SHA where an attestation means anything.
2. **`check-workflow-gate.py` now enforces that**, unconditionally, for every workflow
   owning a required context — and in *both* halves, which the first cut missed. Declaring
   `merge_group:` on the workflow is not enough: if the job producing the context carries
   an `if:` excluding the event it lands `skipped`, and a skipped check **satisfies** a
   required status check. That fails green, which is worse than the deadlock it replaces
   and is the "looks like coverage, enforces nothing" pattern this script opens by
   rejecting. `code-reviewer` caught the gap; the guard now checks the job condition too.

   The condition check is textual, since GitHub expressions cannot be evaluated here, and
   it fires only when a condition *branches on the event* without naming `merge_group`. A
   naive "must contain merge_group" version flagged `ci.yml`'s `required-checks`, which
   carries `if: always()` and does run in a queue — that would have blocked all of CI the
   moment the guard shipped. Caught by this script's own negative test, not in review.

   Verified in three directions: passes as shipped; rejects a workflow that drops the
   `merge_group:` trigger; rejects a workflow that keeps the trigger but narrows the
   producing job's `if:`.
3. **`rearm`'s two silent-green exits now fail.** "No run after 180s" and "still in flight
   after 180s" both returned 0, so `rearm` reported green over a pull request that stayed
   red — CLAUDE.md §9's signature, from a direction the guard does not model. Harmless
   while advisory; the thing itself once required. `rearm` gates nothing and is in
   `NON_BLOCKING`, so its red is pure signal and costs no merge.
4. **The failure help now tells the truth about who can attest and what to do.** It
   previously instructed fork contributors to post a comment that is silently discarded
   (`author_association` is `CONTRIBUTOR`, and `rearm`'s guard means no job even starts),
   and never mentioned the recoveries that work when commenting does not: close/reopen,
   draft/ready, or a manual re-run. `CLAUDE.md` §7b, which told agents there was exactly
   one required context, is corrected too — an agent checking only `required-checks` would
   conclude green and hit a merge refusal with no explanation attached to any job.

**Rejected:**

- *Wait out the stated week.* `pre-pr-attest.yml`'s own header required "a week of clean
  runs with the re-arm behaviour confirmed on a real push". The re-arm half is satisfied —
  `rearm` fired for real on #112 and #113 and issued genuine re-runs — but roughly zero of
  the seven days have elapsed, and the matcher and polling loop both changed today in
  `369429a`. Overridden deliberately by the owner; recorded here because the precondition
  was written down and is being knowingly skipped, which is worth more on the record than
  quietly deleting it.
- *Keep it advisory and surface skips another way* (e.g. a sticky comment when the head is
  unattested). Buys most of "skipping becomes visible" without putting a comment-driven,
  best-effort re-run on the critical path to every merge. Lost to the stated goal of
  unattended CI: a check nothing enforces is one an agent can ignore.
- *Report success unconditionally on `merge_group` without explaining why.* That is the
  "looks like coverage, enforces nothing" pattern `check-workflow-gate.py` exists to
  reject. The reasoning is written into the job instead.

**Accepted risks, none of them fixed here.** An adversarial audit raised these; they are
real and are being taken on knowingly:

- **The workflow can trap its own fix.** `pull_request_target` runs the base-branch copy,
  so a defect landing on `main` fails every open pull request including its repair, and
  `bypass_actors: []` leaves no merge-box override for anyone. `/pre-pr-bypass` cannot
  help — it is read *inside* the job, so it clears a missing attestation but not a job that
  died earlier. Escape is out of band: flip `main-protection` to
  `"enforcement": "disabled"`, merge, flip it back. Now written into the workflow header.
- **The full automated path has never completed end-to-end.** The re-run flipping a check
  green was proven manually on #101; `rearm` issuing re-runs was proven on #112 and #113;
  the two composed — push, comment, re-arm, green, unattended — have not been observed on a
  real pull request.
- **A required context is matched by NAME, and the name is not exclusive to us.** The
  ruleset pins `pre-pr-attest` to `integration_id: 15368`, which is the generic "GitHub
  Actions" app — not this workflow and not this job. A `pull_request`-triggered workflow
  runs the file *from the pull request's own branch*, so a PR that adds any workflow with
  a job named `pre-pr-attest` produces a second check run of that name, from that app, on
  that SHA, and can make it succeed trivially. `check-workflow-gate.py` cannot see it: it
  reads the workflows the repo ships, not ones a pull request introduces. **This is
  pre-existing and applies identically to `required-checks`** — promotion did not create
  it, it raised the payoff. Mitigated, not closed, by tightening the fork-PR workflow
  approval policy from GitHub's default `first_time_contributors` to
  `all_external_contributors` in `scripts/apply-repo-settings.sh`, so no outside
  contributor's workflow runs without a human approving it. Deserves its own issue.
- **`rearm` holds `actions: write` behind an `author_association` filter**, and that is
  not a permission check — `COLLABORATOR` includes read-only and triage collaborators.
  Inert while the repository has exactly one collaborator at `admin`, and bounded because
  the comment body reaches only a `contains()` expression, never a shell or a jq program.
  It becomes a live escalation the moment a second, lower-privileged collaborator is
  added. Recorded here because promotion put that job upstream of a blocking check.
- **`rearm`'s 180s deadline is wall-clock**, so it also absorbs queue time against
  CLAUDE.md §7's 20-slot ceiling. A saturated repository can exhaust it with nothing wrong.
- **`rearm` pins `HEAD_SHA` before its poll loop**, so a push landing mid-poll can have it
  re-run a superseded run.
- **`rearm` triggers on a bare `contains(body, '/pre-pr')`** with `cancel-in-progress`, so
  any later comment merely mentioning the string can cancel a live re-arm.
- **Fork pull requests now need a maintainer to attest after every push**, and an honest
  attestation means having run `scripts/gates.sh` — i.e. executed the contributor's code
  locally. `claude-review.yml` skips forks deliberately; this does not.
- **The matcher has no automated test.** The workflow used to claim one — "caught by the
  fixture test below" — and there is no fixture test, in that file or anywhere in the
  repository; the false claim is removed in this diff rather than left standing in the one
  file whose correctness now gates every merge. The regex has already shipped two distinct
  matching bugs: the trailing character class, and the leading-whitespace hole fixed in
  `369429a`. It is now the single point of decision for every merge, with no coverage at
  all. Writing that test is the highest-value follow-up on this file.

**Evidence:**

- Live ruleset `21678485` required only `required-checks` before this change, matching the
  checked-in spec.
- `rearm`'s two real executions: run 33802688221 (`re-running 33798210754 (failure) ... on
  #112`) and run 33803430526 (`re-running 33803025633 (failure) ... on #113`).
- Guard verified in both directions, per item 2 above.

**Where:** branch `claude/promote-pre-pr-attest`, carrying `approved:ci-change`. Owner
merges (CLAUDE.md §7b: protected path). **Not live until `scripts/apply-repo-settings.sh`
runs** — the JSON is a record, and that script is what PUTs it to GitHub.

**Left open:** the matcher test (highest value per line of anything on the risk list),
`rearm`'s `HEAD_SHA` refresh and anchored comment match, the fork-PR posture, and
`claude-review-gate` still being red by construction on any PR editing its own file.
