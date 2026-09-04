## 2026-09-03 — the two advisory gates were red for reasons no PR could fix

**Decided:** `pre-pr-attest` splits into two jobs. `pre-pr-attest` answers
`pull_request_target` only and is the sole reporter; a new `rearm` job answers
`issue_comment` and re-runs that run instead of reporting anything itself. The shared
`concurrency` group is keyed on the event as well. `claude-review-gate`'s `--max-turns`
goes 30 → 100. `main-broken` gains `mypy` in its `needs`.

**Rejected:**

- *Post a commit status on the head SHA from the `issue_comment` run.* Works, and is the
  usual pattern — but a commit status and a check run named `pre-pr-attest` are two
  different entries in the rollup, and the ruleset entry
  `{"context": "pre-pr-attest", "integration_id": 15368}` would then match ambiguously
  once promoted. Re-running keeps exactly one thing producing that name.
- *Create a check run on the head SHA via the Checks API from the comment run.* Same app
  (15368) and same name, so it would supersede — but it puts a second producer of the
  gate's name in the repo, and `check-workflow-gate.py` reasons about job names, not
  about check runs conjured by an API call.
- *Let `rearm` into `required-checks.needs`.* It is `skipped` on every pull request and
  the gate counts `skipped` as failure. That is the deadlock `check-workflow-gate.py`
  exists to catch, reached from the other side — hence `NON_BLOCKING` instead.
- *Make `claude-review-gate` neutral instead of failing closed when it cannot produce a
  verdict.* A review that did not happen must not read as a pass. The cap was the bug,
  not the fail-closed.
- *Bump `actions/checkout@v4` → v5 to clear the Node 20 deprecation warnings.* Warnings,
  not failures, and the v4 pin is a documented choice matching every job in `ci.yml`.
  Separate change.

**Evidence:**

- All 11 `pull_request_target` runs of `pre-pr-attest` on record are `failure`. The
  `issue_comment` re-runs are attributed to `main` and land their check on main's tip:
  run 33725505997 put a `pre-pr-attest` **failure** on `db066c1eb`, which is `main`.
  #92, #94, #97, #89 and #101 all merged with the check red.
- Re-running 33784923222 — the run the new `rearm` job would target for #101 — flipped
  that PR's `pre-pr-attest` from `FAILURE` to `SUCCESS` with the owner's attestation
  comment already in place and nothing else changed. That is the whole mechanism, tested
  on a real head.
- `claude-review-gate` turn counts: successes at 17, 27, 28; failures at 31, 31, 31, 31,
  32, 42 — six, matching the six of the last fifteen runs that died on the cap. Of those
  six, four hit `error_max_turns` at 31 and two had already produced a verdict — *"Claude reported a successful result
  after 42 turns, exceeding the configured maximum of 30"* (runs 33660939981,
  33696385723). Only 33663943738 failed on an actual finding.
- `main-broken.needs` listed `changes, python-quality, invariants, contract,
  contract-freeze, python-tests` — every required job except `mypy`.

**Where:** PR #104, branch `claude/github-ci-checks-2969e5`, carrying `approved:ci-change`.
Owner merges (CLAUDE.md §7b: protected path).

**Left open:** `bt5[screen]` declares `commec>=0.2`, which is not on PyPI, so any
whole-project resolve (`uv run`, `uv lock`) fails. CI never installs `[screen]` — it names
`[dev,export,fold]` explicitly — so no check is affected today. Worth an issue under §5
rather than a dependency edit.
