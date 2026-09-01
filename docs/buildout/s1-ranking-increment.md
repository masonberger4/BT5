# S1 — The ranking increment

*Copy this whole file into a fresh Claude Code session on the BT5 repo.*

**Launch:** permission mode `acceptEdits` · model opus · effort high (both are the
repo default from `.claude/settings.json`, so no override needed) · runs unattended.
**Do not run this in plan mode** — five other sessions are running in parallel and the
gate is the draft PR, not a plan approval.

---

You are closing the highest-value gap in BT5: the walking skeleton produces one
unranked candidate, and the product is ranks.

## Read this first

Run `/bootstrap` before anything else — a fresh checkout has no `.venv`, and
`gates.sh` exit **10** means BROKEN, not a code failure. Then read `CLAUDE.md`; it
loads automatically but the lane table and §3 correctness rules are load-bearing here.
Then read `docs/buildout/README.md` for the un-draft queue and the inter-session
contracts you are bound by.

Your branch: **`claude/s1-ranking-increment`**, cut from `main` (green at `628e130`).

## The situation

PR #71 landed `bt5/design/` and its docstring says exactly what it refuses to do:

> It ships one candidate, no gallery. It scores nothing — every objective is reported
> `unavailable` with a reason, never omitted, so the scorecard cannot look complete
> when it is not. It carries no baseline: `native_baseline` stays None …

Meanwhile `bt5/score/` already exports everything needed: `build_gallery`, `sweep`,
`greedy_max_min`, `simplex_weights`, `G4_MIN_PAIRWISE_DISTANCE`, `null_distribution`,
`percentile_of`, `normalise`, `synonymous_variant`, `order_entries`, `write_csv`,
`write_idt_plate`. `grep` will confirm nothing under `design/` calls any of them —
`design/runner.py:358` literally appends the note *"ranking not computed: no null
distribution and no percentiles"*.

**Your job is the wiring, and the honesty that has to survive it.**

## What to build

Target is `docs/PLAN.md:490-495`'s v1 bar. Use `docs-miner` to pull it verbatim rather
than reading `PLAN.md` inline (58 KB).

1. **Percentiles against the null.** Wire `null_distribution` / `percentile_of` into
   the design path so every scored objective reports a percentile against the
   random-synonymous null rather than `unavailable`. `DEFAULT_NULL_N` is the shipped
   size; PLAN's G3 gate wants the null under 2 s.
2. **A 5-candidate gallery.** Replace the single candidate with `build_gallery`.
   PLAN's G4 gate requires ≥15% pairwise codon distance — `G4_MIN_PAIRWISE_DISTANCE`
   and `greedy_max_min` already exist for exactly this. G4's failure invalidates a
   *product* decision, not a technical one; treat a diversity failure as a finding,
   never as a threshold to lower.
3. **`native_baseline` as a first-class candidate.** `core/result.py` already has the
   field and `docs/PLAN.md:132` calls it "a field, not a UI afterthought". It is a
   claim about a real wild-type CDS, so it is populated when the caller supplies one
   and stays `None` — with the reason stated — when they do not. Do not fabricate a
   baseline by back-translating.
4. **Confidence bands, and `QcReport.is_complete` meaning something.** Today it is
   always `False`. Make it true when the scorecard genuinely is complete, and keep
   every `ObjectiveScore.unavailable` carrying its reason rather than being dropped.
5. **The order CSV.** Wire `order_entries` + `write_csv` / `write_idt_plate` so the
   design emits the vendor order file alongside the GenBank.
   `docs/PLAN.md:198-200` locks the output as "annotated construct + vendor order CSV".
6. **A timing assertion.** PLAN's G7 bar is ≤10 s end to end at 500 aa. A plain
   `pytest` timing assertion, marked `slow` if it needs to be. **Do not create
   `benchmarks/`** — that directory does not exist and creating it is an owner
   decision under `approved:algorithm-change`.

## The line you must not cross

BT5 **never reports a predicted expression number** (`CLAUDE.md` §0,
`docs/PLAN.md:126-128`). A CI gate — `tests/data_integrity/test_no_expression_claims.py`
— bans prediction vocabulary from the schema, and it will fail you. Ranks,
percentiles, confidence bands. `native_baseline` — "don't optimize" — is a first-class
output, not a fallback.

This is the session where that pressure is highest, because ranking things is exactly
where a number wants to become a prediction. A percentile says *where this sequence
sits against its own synonymous null*. It does not say the protein will express.

## Files

**You own:** `packages/engine/src/bt5/design/**`, `packages/engine/src/bt5/score/**`,
`packages/engine/tests/design/**`, `packages/engine/tests/score/**`.

You get `score/` as well as `design/` because the two are inseparable for this change
and no other session writes score.

**Never touch:** `rules/`, `solver/`, `vector/`, `cassette/`, `codon/`, `structure/`,
`core/`, `verify.py`, `.github/`, `pyproject.toml`, `data/`, `tests/contract/`,
`tests/invariants/`, `tests/data_integrity/`.

If the wiring genuinely needs a new `core/` field, **stop**: open an issue naming the
type, and use `/architect` for the MINOR-vs-MAJOR call. `core/` is a global mutex —
another session may be waiting on it.

## Your contract with the other five

- **`design()`'s signature is frozen** (`design/runner.py:156-171`, keyword-only,
  `table_id` never defaulted). You may **add** fields to `SkeletonResult`; you may not
  remove or rename one. S5 is building a CLI against it right now.
- **Render whatever `BiosecurityVerdict` you are given** and never print "clear" for
  `not_run`. S2 is making the screen real in parallel; you consume it, you do not
  implement it.

## Delegation

- `docs-miner` — PLAN's v1 bar, the report contract, the banned-vocabulary rules.
  Never read `PLAN.md` or `brief.md` inline.
- `Explore` — every existing call site of `build_gallery`, `null_distribution`,
  `order_entries`, and every construction of `ObjectiveScore.unavailable`.
- `gate-runner` — `scripts/gates.sh` before each push.
- Main thread — the judgment: what a percentile is allowed to claim, when
  `is_complete` may be true, why `native_baseline` is an output and not a fallback.
  That is not delegable.
- `/architect` — only if this needs `core/`.

No file in your lane exceeds 20 KB, so bare `Read` is fine here. If you stray into
`vector/`, `backbone.py` (28.9 KB), `kmers.py` (21.2) and `assemble.py` (18.7) are over
the limit — use `offset`/`limit` or `Explore`.

## Done means

- `bash scripts/gates.sh` reaches `ALL GATES PASSED`. Exit 10 = run `/bootstrap`;
  pytest exits 2, 3, 4, 5 are BROKEN too, and 5 ("no tests collected") is never
  success.
- A design emits a 5-candidate gallery passing G4's ≥15% diversity bar, percentiles
  against the null, an annotated GenBank **and** an order CSV.
- `QcReport.is_complete` is `True` on a complete scorecard and `False` with reasons
  otherwise.
- `tests/data_integrity/test_no_expression_claims.py` passes — you added no prediction
  vocabulary.
- `/pre-pr` is clean, and the PR is **open as a draft**.
- You added a decision file at `docs/decisions/2026-XX-XX-<slug>.md`: what you decided,
  what you **rejected** and why, with evidence. One file per decision — never append to
  a shared one.

- **`/pre-pr` is run by the operator, not by you.** It is
  `disable-model-invocation: true`, so a session cannot self-invoke it and must not
  replicate its steps by other means. Ask for it when the branch is ready.
- **The attestation is posted last.** After `/pre-pr` and after the final push, comment
  the full 40-character head SHA on the PR:

  ```
  /pre-pr <head-sha>
  ```

  The advisory `pre-pr-attest` check reads that comment. An attestation names **one**
  commit, and pushing again makes it stale on purpose — a review of the previous tree
  says nothing about this one. Never attest a SHA that was not just reviewed; the whole
  value is that the claim is on the record. If a gate or review came back blocking and
  you are pushing anyway, do **not** attest — say so in the PR and let the check stay
  red. Only the owner may waive it, with `/pre-pr-bypass <head-sha>`.

**Do not self-merge.** Your scientific impact is non-"none" — you change what the app
produces, from one unranked candidate to a ranked gallery with an order file. Under
`CLAUDE.md` §7b that goes to the owner. Say in the PR that you know this.
