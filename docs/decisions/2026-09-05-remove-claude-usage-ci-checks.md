## 2026-09-05 — remove every CI check that spends Claude usage

**Decided:** delete `.github/workflows/claude-review.yml` and
`.github/workflows/pre-pr-attest.yml`, and drop `pre-pr-attest` from the required
contexts in `.github/rulesets/main-protection.json`. Owner's call, taken on cost
grounds. `required-checks` is now the ONLY required context.

Two checks spent Claude usage, in different currencies, and both are gone:

1. **`claude-review-gate`** (`claude-review.yml`) called
   `anthropics/claude-code-action@v1` on `claude-sonnet-5` with `--max-turns 100` and
   `timeout-minutes: 25`, on every non-draft push to a same-repo pull request. It was
   the direct spend: a full review session per push. It was **advisory** — absent from
   the ruleset and from `required-checks.needs` — so removing it blocks nothing and
   changes no merge outcome.
2. **`pre-pr-attest`** spent nothing on a runner; it is bash and `jq`. It spent Claude
   usage by *construction*: it was a required context that only went green when an
   agent ran `/pre-pr` — `gate-runner`, `code-reviewer`, and conditionally
   `rule-auditor` and `security-reviewer` — and posted `/pre-pr <head-sha>`. An
   attestation names one commit and goes stale on the next push, so the full subagent
   chain was owed again on every push. That is the larger of the two bills.

### What this costs, stated plainly

**Skipping `/pre-pr` is invisible again.** That was the entire point of the check:
`pre-pr-attest.yml`'s own header recorded that PR #79 merged as `628e130` with `/pre-pr`
never having run and nothing saying so. `.claude/.pre-pr-marker` cannot replace it — it
is gitignored (`.gitignore:19`) and digests `git status` + `git diff HEAD`, both empty
and constant on a clean CI checkout, which is why the comment-based check existed at all.
CLAUDE.md §8 is now discipline, not a gate.

**Nothing audits rule provenance.** `tests/data_integrity/test_rule_contract.py` catches
*missing* provenance; it cannot catch *wrong* provenance — it verifies citations exist,
are https and are labelled, but cannot read the cited paper and judge whether it supports
the threshold citing it. `claude-review-gate` was the only check that tried, and its
`provenance: UNSUPPORTED` verdict failed the build on its own. That gap is why rule C1
merged as `628e130` with its `(0.70, 0.90)` band and 0.2 weight never audited, and the
gap is open again. `/verify-provenance` remains as a manual rotation.

Both losses are accepted, not overlooked. The mitigation is the `/pre-pr` skill, which
survives — see below.

### The ordering hazard

`main-protection.json` is a **spec file, not live config**. The live ruleset still names
`pre-pr-attest` with `bypass_actors: []`. A required context whose workflow no longer
reports blocks every pull request permanently, with no error and no bypass — CLAUDE.md
§9's silent killer, and `pre-pr-attest.yml`'s header warned that the escape is out of
band. So:

1. Merge this PR **or** apply the ruleset first — either order is safe for *this* repo
   state, because `pull_request_target` runs the BASE branch copy: until this merges,
   the old workflow still reports on open PRs.
2. **The owner must then run `bash scripts/apply-repo-settings.sh`** from an
   authenticated machine. It cannot run from a Claude Code web session (the proxy blocks
   repo-scoped REST paths, and no MCP ruleset tool exists).
3. Until step 2 runs, every open pull request will sit on `pre-pr-attest` "Expected —
   Waiting for status to be reported" forever. `/pre-pr-bypass` does **not** help: the
   job that read it is gone.

This PR itself is gated by the base-branch copy of the check, so it needs
`/pre-pr <sha>` or `/pre-pr-bypass <sha>` from the owner to go green, or the ruleset
applied first.

### What was deliberately NOT removed

- **The `/pre-pr` skill.** The ask was to stop *CI* spending usage on every push, not to
  forbid running reviewers locally. The chain stays; only its step 10 (posting the
  attestation) is gone, since the check that read it is gone.
- **`main-broken`** — post-merge safety net, `push`-only, no Claude usage.
- **`required-checks`** and every job feeding it: ruff, mypy, invariants, contract,
  contract-freeze, python-tests, approvals.
- **The `all_external_contributors` fork-PR policy** in `apply-repo-settings.sh`. It
  guarded a check-run name collision against *any* required context; dropping
  `pre-pr-attest` narrowed the target but did not remove it.
- **`ci.yml`'s `merge_group:` trigger**, which `check-workflow-gate.py:177` requires
  unconditionally of any workflow producing a required context.

### Collateral edits, and why each was forced

- `.github/scripts/test-attestation-matcher.py` **deleted** — it read the real `jq`
  program out of `pre-pr-attest.yml`; with that file gone it fails rather than skips.
  Its `ci.yml` step went with it. Verified by running it before and after.
- `.github/scripts/check-workflow-gate.py` — `rearm` removed from `NON_BLOCKING`; that
  job no longer exists. `main-broken` stays.
- `.github/scripts/check-hook-commands.py` and `.claude/hooks/run_py.sh` — both
  justified a claim by naming `claude-review`. The claims survive; the justifications
  were reworded. Neither script ever opened a workflow file.
- `CLAUDE.md` §7b said the ruleset requires **two** contexts and warned that checking
  only `required-checks` was "the mistake this wording exists to prevent". Now one.

**Scientific impact: none.** No engine source, rule, threshold, genetic-code table or
output schema is touched. This changes only what CI checks, never what the app builds.
