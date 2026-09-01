---
name: pre-pr
description: Run the full local gate chain, then conditional review and security passes, determine the approval labels, and open a draft PR.
disable-model-invocation: true
allowed-tools: Bash, Read, Grep, Glob, Agent
---

# /pre-pr — the gate before the gate

## Verdict rules — read these first

- **Any gate exiting non-zero is RED.** There is no "probably unrelated".
- **`scripts/gates.sh` exit 10 is BROKEN, not FAIL** — there is no `.venv`. Run
  `/bootstrap`. Do not report a test failure.
- **pytest exits 2, 3, 4 and 5 are BROKEN**, not FAIL: collection error, internal error,
  usage error, and no-tests-collected. Never report any of them as PASS.
- **Never re-run a gate with looser flags to get green.** Never skip, `xfail`, loosen or
  delete a test (CLAUDE.md §4). If a property fails and also fails on the merge base, it
  is a pre-existing bug: record a fixture under `tests/data/regressions/`, open an issue,
  and say so in the PR.
- **mypy is a required CI job** (as of #63) and is in `required-checks.needs`. Running
  it here turns a merge-gate failure into a local one.
- **This is not CI parity on invariants.** Locally Hypothesis runs `dev` at 50
  examples; CI runs `ci` at 200. A property can pass here and fail there. Pytest's run
  header prints the profile and budget — read it rather than assuming.

## The chain, in order

**1. Gates.** Issue this as **one compound command**:

```bash
set -e; bash scripts/gates.sh
```

It must stay compound. The `PreToolUse` Bash hook appends truncating flags to *simple*
commands and stands down on compound ones — so a later session that "simplifies" this
into separate calls silently re-enables truncation on the one run that must be complete.

**2. Review.** Delegate the branch diff to `code-reviewer`.

**3. Security — only if warranted.** Run `security-reviewer` **iff**
`git diff --name-only origin/main...HEAD` touches any of
`packages/engine/src/bt5/vector/`, `packages/engine/src/bt5/core/services.py`,
`packages/engine/src/bt5/verify.py`, `packages/engine/src/bt5/cassette/`, `.github/`.
Conditional invocation is what stops three agents firing on one diff.

**4. Rule provenance — only if warranted.** Run `rule-auditor` **iff** the diff changes a
Spec's `citations`, `weight_provenance`, `enforcement`, `last_verified` or a threshold.

**5. Paired tests.** Every catalog rule needs one:

```bash
comm -23 \
  <(ls packages/engine/src/bt5/rules/catalog/*.py | xargs -n1 basename | sed 's/\.py$//' | grep -v '^__init__$' | sort) \
  <(ls packages/engine/tests/rules/test_*.py | xargs -n1 basename | sed 's/^test_//;s/\.py$//' | sort)
```

`d1_restriction_sites` is a known pre-existing gap. Anything else is blocking.

**6. Agent memory.** `git status --porcelain | grep agent-memory` — anything listed is an
unreviewed instruction channel about to be committed. Read the diff before proceeding.

**7. Labels.** Determine from `.github/scripts/check-approval-labels.sh`:

| path | label |
|---|---|
| `verify.py`, `tests/invariants/`, `tests/data_integrity/` | `approved:oracle-change` |
| `core/`, `tests/contract/` | `approved:contract-change` |
| `benchmarks/baseline.json`, `benchmarks/tolerances.yaml` | `approved:algorithm-change` |
| `data/genetic_codes/`, `data/codon_usage/` | `approved:data-change` |
| `.github/` | `approved:ci-change` |

State which labels the PR needs. The `approvals` job fails without them.

**8. Push marker.** Record `git rev-parse HEAD` plus a hash of `git status --porcelain`
and `git diff HEAD` to `.claude/.pre-pr-marker` so the push gate knows this ran against
this exact tree.

This marker is **local only** — it is gitignored (`.gitignore:19`) and its digest is empty
on a clean checkout, so it can never be CI's evidence. Step 10 is what makes this run
visible off this machine.

**9. PR.** Open it as a **draft** (CLAUDE.md §7 — drafts skip the expensive jobs, and CI
capacity is the binding constraint: 20 slots, ~12 per Python PR, so at most 5 open
non-draft PRs). Fill in `.github/pull_request_template.md`, including the **Scientific
impact** section: say what changed about the sequences the app produces, or "none".

**10. Attest, once the PR exists and the final push has landed.** Comment on the PR with
the full 40-character SHA of the head being reviewed:

```
/pre-pr <head-sha>
```

The `pre-pr-attest` check reads that comment. **Attest last** — an attestation names one
commit, and pushing again makes it stale on purpose, because a review of the previous tree
says nothing about this one. Do not attest a SHA you have not just reviewed; the whole
value of the check is that the claim is on the record and auditable.

If a gate or a review came back blocking and you are pushing anyway, do **not** attest —
say so in the PR instead and let the check stay red. Only the repo owner may waive it, with
`/pre-pr-bypass <head-sha>`.

## Report

```
GATES        <per-gate line from gates.sh>
REVIEW       <blocking count> blocking, <n> non-blocking
SECURITY     run | skipped (surface not touched)
PROVENANCE   run | skipped (no Spec provenance changed)
PAIRED TESTS ok | missing: <ids>
MEMORY       clean | <files>
LABELS       <labels needed>
ATTESTED     <head-sha> | withheld: <why>
VERDICT      READY (draft) | BLOCKED: <what>
```
