## 2026-09-04 — every check is agent-clearable; the merge button is the owner's

**Decided:** two changes, from one owner instruction: *"Modify CI checks on GitHub so that
they can be fully done autonomously. I only want merge permission left to me, and I can
still delegate that to an agent if I so choose."*

1. **CLAUDE.md §7b loses the standing agent-merge permission.** It read *"An agent may
   squash-merge its own PR once CI is green"* with four carve-outs routed to the owner. It
   now reads: an agent takes a PR to green and stops, reporting readiness; the owner
   merges, delegating per PR. This is a **tightening**, and the only one here.
2. **`claude-review-gate` stops being red by construction** on a pull request that edits
   `claude-review.yml`. That branch now passes on a **checked** compensating control —
   `approved:ci-change` must be on the PR — instead of failing closed forever.

**What was already autonomous, checked rather than assumed.** The audit found far less
broken than expected — and after the reversal below, only one thing actually needed fixing.
Against the live check state of the five open pull requests (#131, #139, #142, #143, #144),
every job in `required-checks` is an ordinary code gate; `approvals` needs an `approved:*`
label and all five labels exist and are agent-appliable; `pre-pr-attest` was green on every
open PR, because agents comment as the owner and `2026-09-03-agents-may-attest.md`
authorised exactly that; and review threads are already being resolved by agents (#142, two
threads). Nothing in that list needed changing, and nothing in it was changed.

---

## The bypass actor: proposed, chosen by the owner, then reverted

**This is the substance of the slice, and it is a record of getting it wrong.**

`main-protection.json` briefly gained the Repository admin role in `bypass_actors` with
`bypass_mode: pull_request`. **It is reverted; the file is byte-identical to `main`.**

**The argument for it.** With `bypass_actors: []` the owner does not hold merge permission
unconditionally — they hold it conditionally on Actions being healthy. `pre-pr-attest.yml`'s
header spells the trap out: a defect in that job fails every open pull request *including
the one that would repair it*, and the escape is out of band, unreachable from a Claude
session because the proxy blocks every repo-scoped REST path
(`scripts/apply-repo-settings.sh:5-7`).

**Why that argument fails, on two independent grounds.**

*First: the repository had already considered and rejected it, in the document this very
workflow cites as its design source.* Nobody checked before proposing it.

- `docs/research/github-setup.md:244` — *"**The bypass list is empty — including for you.**
  🔴 This is the single highest-value setting here. Rulesets do **not** implicitly exempt
  repo admins (the opposite of legacy branch protection). Red CI genuinely blocks the one
  human who can click Merge. AI agents produce plausible code that fails at the edges; CI
  being mandatory *for you* is the only thing standing between that and `main`."*
- `:246` — *"**Emergency escape hatch — do not add a bypass actor.** Bypass is per-*actor*,
  not per-*rule*…"*
- `:2828`, under **`### Never do`** — *"Add a bypass actor."* is the first item of six.
- `:2903`, §8.10, titled *"Everything blocks because CI itself broke"* — the exact scenario
  argued above — *"**Mitigation.** That is the intended cost, and it's the right one…
  **Never** solve it by adding yourself as a bypass actor."*

*Second, and specific to this repository: **agents authenticate as the owner.*** A
permissions audit verified this session's GitHub credential resolves to `masonberger4`, the
repository owner and therefore a holder of the Repository admin role. GitHub cannot
distinguish the owner at a keyboard from an agent carrying the owner's token. So a bypass on
that role is a bypass held by **every Claude session** — it converts a MECHANICAL guarantee
(nobody merges red) into a PROSE one (CLAUDE.md §7b asking agents not to), addressed to the
same actor that holds the credential. It therefore defeats the instruction that motivated
it: it does not leave merge permission to the owner, it hands merge-past-red-CI to every
agent.

The argument also weakened on its own terms. Change 1 makes the owner the only party who
merges anything, so requiring them at an authenticated machine for a rare un-wedge costs
almost nothing that owner-only merging does not already cost.

**Also true, and understated while it was in the tree:** bypass is per-actor, so at merge
time it applies to *every* rule in the ruleset, not only `required_status_checks` — including
`required_review_thread_resolution`, which `ci.yml:20-25` calls *"the one conversation gate a
zero-approval solo-owner repo actually relies on"* and which has no second holder.
(`allowed_merge_methods: ["squash"]` does have one: `allow_merge_commit=false` in
`apply-repo-settings.sh`.) The PR body, this record's first draft, §7b and the
`pre-pr-attest.yml` header all described the bypass as reaching only a red required check.

**What replaces it.** `scripts/apply-repo-settings.sh` now carries the wedge-recovery
runbook the research document prescribes — disable `main-protection` briefly, merge the
repair, re-enable, then read `/rulesets/{id}/history` — and its verification step now expects
**no** bypass actors and warns if it finds any. The `emergency-bypass` placeholder ruleset
named at `:246` was **not** created: it carries no rules and is permanently `disabled`, so it
does nothing the runbook does not, and a second ruleset in the list invites the misreading
that toggling *it* is the mechanism. The document's actual mechanism toggles
`main-protection`.

**Process note worth keeping.** The owner was asked to choose between "add a bypass actor"
and "leave it empty plus a runbook", and chose the former on a recommendation. Neither option
was the one the repository's own research prescribes, because `github-setup.md` had not been
read — `docs-miner` exists for exactly that and was not used until a code review forced it.
The question was re-put with the evidence and the credential fact, and the answer changed.
**Check `docs/research/` before proposing a change to a setting that document configures.**

---

## Rejected

- *A `/merge <head-sha>` comment command, OWNER-only, verifying the SHA is current and the
  required checks green before squashing.* Put to the owner as the recommended option; the
  owner chose plain owner-only merging. It would have made delegation work from a phone and
  put each delegation on the timeline, at the cost of a workflow holding `contents: write` on
  an `issue_comment` trigger — a merge button reachable by comment, in a repository that
  already declines to let Actions approve pull requests for the same class of reason.
- *Making `claude-review-gate` neutral whenever it cannot produce a verdict.* Already
  rejected in `2026-09-03-ci-checks-that-can-go-green.md` — "a review that did not happen
  must not read as a pass" — and that rejection **stands**. The narrow branch taken here is a
  different shape: a crash, a timeout or a turn cap is a fault where a re-run can still
  deliver the verdict, so passing would lose a review that was owed; a self-edit is
  deterministic, and no push to that branch can ever produce one. The fault branch is
  untouched and still fails closed.
- *Passing the self-edit case unconditionally.* The label is checked live, from the `pulls`
  endpoint, and its absence still fails. Read from `pulls` and not `issues/$PR/labels`: the
  latter maps to the `issues` scope, which the workflow does not grant, so it would 403 into
  `set -e` and land back on the permanent red being removed.
- *A `workflow_dispatch` workflow applying the ruleset with an admin PAT*, to close the last
  human step. A workflow able to rewrite the rules that gate it, reachable by anyone who can
  edit `.github/` — the exposure `apply-repo-settings.sh` already reasons about when it sets
  `can_approve_pull_request_reviews=false`.

## A real bug the audit found in the self-edit fix, and fixed here

The first cut detected the self-edit with `git diff --quiet "$BASE_SHA" HEAD`. `BASE_SHA` is
`base.sha` from the event payload, fixed when the event fired; `HEAD` is `refs/pull/N/merge`,
whose first parent is the base tip **as of the merge-ref computation**. Diffing one against
the other attributes anything `main` gained in between TO THE PULL REQUEST.

Reproduced on real commits — base M1 → M2 where M2 edits `claude-review.yml`, a branch off M1
touching only `other.txt`, HEAD the merge ref:

| predicate | verdict |
|---|---|
| `git diff $BASE_SHA HEAD` | **SELF_EDIT=yes** — wrong |
| `git diff $BASE_SHA...HEAD` (three-dot) | **SELF_EDIT=yes** — three-dot does not fix it |
| `git diff HEAD^1 HEAD` | SELF_EDIT=no — correct |

Reachable rather than theoretical: a review faults (six of fifteen runs once died on the turn
cap), `main` then merges a change to this file, someone clicks re-run failed jobs — and the
re-run replays the stale payload base against a freshly recomputed merge ref. With
`approved:ci-change` present, which any `.github/` pull request must carry anyway, a
**fault-class** empty verdict would take the self-edit branch and exit 0 — exactly the case
the fault branch exists to fail closed on.

The guard tests `HEAD^2`, not `HEAD^1`: every commit but the root has a first parent, so a
`HEAD^1` probe verifies on a plain checkout and would compare HEAD against its predecessor.
A second parent exists only on a merge ref, which is the sole condition under which `HEAD^1`
means "the base tip".

**Where:** branch `claude/autonomous-github-ci-boxw2y`, PR #145, carrying `approved:ci-change`.

**Left open:**

- `security-reviewer` **cannot be invoked at all**: its own system prompt trips a `[bio]`
  content filter on both `claude-opus-5` and `claude-sonnet-5`. Three independent sessions hit
  it on 2026-09-04, and `2026-09-03-agents-may-attest.md` recorded the same failure. One of
  `/pre-pr`'s four legs is therefore unrunnable, and a generic permissions audit was
  substituted again. Worth an issue: the substitution works but is not the documented chain.
- `claude-review-gate` produced no verdict twice on #143, a README-only change. Advisory, and
  re-runnable, so it blocks nothing — cause not chased.
- `.claude/verify-setup.sh`'s CLAUDE.md length check fails at 215 lines against a 200 budget.
  It already failed on `main` at 209; this slice adds 6. A trim pass is queued separately.
