---
name: pre-pr
description: Run the full local gate chain, then conditional review and security passes, determine the approval labels, and open a draft PR.
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
- **This chain never merges, and neither does its caller.** It ends at a draft PR and an
  attestation. Taking the PR the rest of the way to green is yours too — labels, review
  threads, a re-run of a flaked job — but the merge is the owner's, delegated per PR
  (CLAUDE.md §7b). Report readiness; do not press the button.

## The chain, in order

**1. Gates.** Delegate the whole chain to `gate-runner` — this is the "Gates →
`gate-runner`" route CLAUDE.md names. It owns the Phase 0 environment probe and the exit
vocabulary above, and it runs on haiku, so the full ruff/mypy/pytest output never enters
this window. That is the entire point of the routing table in `SETUP-NOTES.md`: running
the chain inline spends the most expensive context in the repo on output nobody reads.

Take its report verbatim into the `GATES` line below. If it returns `ENV: MISSING` or
`ENV: INCOMPLETE`, that is **BROKEN**, not FAIL — run `/bootstrap`, then restart this
chain.

Also run the routing verifier, which `gates.sh` does not cover:

```bash
bash .claude/verify-setup.sh
```

It checks the things whose failure mode is silence — a hook that no longer runs, an
agent file with no `description`, a rules glob matching nothing. Its first check is a
negative control proving it can still fail; if that control does not report `rc=2`, the
verifier is blind and every `ok` after it is worthless. Treat non-zero as RED like any
other gate. It needs no `.venv` for the hook probes, but PyYAML (hence `/bootstrap`) for
the frontmatter sections — it says so rather than skipping silently.

Only if `gate-runner` is unavailable, run it here yourself:

```bash
bash scripts/gates.sh
```

A simple command is fine — `gates.sh` is on `compact_output.py`'s `NEVER` list
(`.claude/hooks/compact_output.py:35-42`), so the truncating-flags hook never rewrites
it, compound or not. What must not happen is running the gates *individually* to save
time: four of five green reported as green is the failure this step exists to prevent,
and `gates.sh` runs each one independently for exactly that reason.

**2. Review.** Delegate the branch diff to `code-reviewer`.

**3. Boundaries — only if warranted.** Run `boundary-reviewer` **iff**
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
