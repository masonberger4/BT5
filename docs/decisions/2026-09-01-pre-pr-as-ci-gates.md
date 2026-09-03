## 2026-09-01 — The `/pre-pr` judgment layer becomes two CI checks, both advisory

**Decided:** `pre-pr-attest` (deterministic, no model) and `claude-review-gate` (adapted
from the design at `docs/research/github-setup.md:1163-1262`) ship as workflows, both
**advisory** — absent from `required-checks.needs` and from the ruleset, with the
promotion diff recorded in each file's header. `/pre-pr` gains step 10: comment
`/pre-pr <head-sha>` once the PR exists.

The problem: three mechanisms that looked load-bearing gated nothing.
`.claude/hooks/push_gate.py` emits `permissionDecision: "ask"` but `.claude/settings.json`
allow-lists every path it guards, so the prompt never surfaced — verified by pushing
`6dc4a0c` with no marker and no prompt. `agent-approval-check` and `claude-review-gate`
were designed in `docs/PLAN.md` and never built. §7b's carve-outs are prose no code
evaluates. A sweep found 255 merge-policy statements across the repo; none of the blocking
ones is enforced.

The cost, concretely: rule C1 merged as `628e130` with `/pre-pr` never having run, so
nobody checked whether its `(0.70, 0.90)` band and `default_weight = 0.2` are supported by
their citations. CI structurally cannot — `test_rule_contract.py` says so itself: it
"catches MISSING provenance; it cannot catch WRONG provenance."

**Rejected:**
- *Adding either job to `required-checks.needs`.* That gate counts `skipped` as failure
  (`ci.yml:295-320`), and both jobs skip on some events. `main-broken` is already excluded
  for exactly this reason. They become **parallel required contexts** when promoted, never
  inputs to the gate.
- *Using `.claude/.pre-pr-marker` as CI's evidence.* It is gitignored (`.gitignore:19`) and
  digests `git status` + `git diff HEAD`, which are empty and constant on a clean CI
  checkout. It cannot travel; a comment naming the head SHA can.
- *Exiting 0 when the review cannot run.* Considered for both the unconfigured-token case
  and the workflow-validation case. Rejected: once this is a required context, a green
  check on a PR that edits the reviewer reopens precisely the hole the upstream guard
  closes. A gate that never ran must not read as a pass.
- *A `NON_BLOCKING` entry in `check-workflow-gate.py`.* Verified unnecessary by running it:
  its orphan check only fires for a workflow containing a job whose context the ruleset
  requires.
- *Changing `Bash(gh pr comment:*)` to `Bash(gh pr comment *)`* on a review agent's
  suggestion. The `:*` form is the permission-rule syntax `.claude/settings.json` already
  uses throughout (`Bash(git push:*)`).

**Evidence:** three behaviours confirmed against live CI before merge — skip on a draft;
fail-closed with the token absent; fail-closed on workflow validation. The third was the
surprise: the Claude GitHub App refuses to mint a token for a workflow file not
byte-identical to the copy on the default branch, so the action self-skipped in 1.5s with
`outcome=success` and no `structured_output`. That is a supply-chain guard against a PR
editing its own reviewer, and it means **the first PR introducing a reviewer can never be
reviewed by it**. Fixture tests of the attestation matcher found a real bug: the trailing
character class `([ \t]|$)` rejected any attestation that was not the entire comment body.

**Known limitation, to resolve before promoting:** a PR touching `claude-review.yml`
always shows red on `claude-review-gate`. Harmless while advisory; a required check red by
construction cannot be satisfied by any push.

**Also unresolved:** `/pre-pr` carries `disable-model-invocation: true`, so an agent cannot
run it and therefore cannot honestly attest. Until that is decided, agent-opened PRs show
`pre-pr-attest` red unless the owner runs `/pre-pr` or comments `/pre-pr-bypass <sha>`.

> **Resolved 2026-09-03** — the flag is removed and an agent may attest, under step 10's
> existing "not if blocking" condition. See `2026-09-03-agents-may-attest.md` for what the
> check proves now and what it stopped proving.

**Where:** PR #84, merged as `81a63d8`.
