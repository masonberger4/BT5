# S2 — The biosecurity screen

*Copy this whole file into a fresh Claude Code session on the BT5 repo.*

**Launch:** permission mode `acceptEdits` · model opus · **effort xhigh** (override;
the repo default is high) · runs unattended.
**Do not run this in plan mode** — five other sessions are running in parallel and the
gate is the draft PR, not a plan approval.

---

You are making BT5's protein screen real. Today it reports `not_run` and the report
says so honestly; your job is to make it actually run, and to keep it just as honest
when it cannot.

## Read this first

Run `/bootstrap` before anything else — a fresh checkout has no `.venv`, and
`gates.sh` exit **10** means BROKEN, not a code failure. Then `CLAUDE.md`, then
`docs/buildout/README.md` for the un-draft queue and your inter-session contracts.

Your branch: **`claude/s2-biosecurity-screen`**. Cut it yourself, from a
**freshly fetched** main — do not reuse a branch someone made earlier, and do not
assume main is where it was when this prompt was written:

```bash
git remote prune origin
git fetch -q origin main && git checkout -B claude/s2-biosecurity-screen origin/main
```

`git fetch --prune origin main` does **not** clear a stale ref — pruning is bounded by
the refspec, so that form prunes only `origin/main` (CLAUDE.md §7a).

## The situation

`core/context.py:96-115` already defines the type, and its docstring states the threat
model better than any summary would:

> Protein-level screening is the one layer BT5's own output cannot defeat: the app's
> core function is producing a functionally identical sequence with maximally
> different nucleotides, which is the textbook method for evading nucleotide-homology
> screening.
>
> `status` is never "clear" when screening did not run. A `NullScreen` reports
> "not_run" so the report cannot imply a clean result that was never obtained.

`commec>=0.2` is **already declared** in the `screen` extra of `pyproject.toml`. You do
not add a dependency — `CLAUDE.md` §5 forbids it and S5 holds that file's mutex
anyway.

## What to build

In `bt5/cassette/` (lane M8, whose documented scope in `docs/PLAN.md:457` is the
protein validator, genetic-code selection, the tag/linker/2A/protease library, the
protein liability scan on the assembled fusion, and biosecurity screening):

1. **A screen behind a protocol**, the way `FoldEngine` is. The concrete commec-backed
   implementation and a `NullScreen` that reports `not_run`.
2. **Degradation that cannot lie.** `/bootstrap` installs `dev,fold,export` and
   **not** `screen`, so **the degraded path is the one CI actually exercises**. Getting
   that path right is most of this task. Mirror how the fold engine degrades — see
   `structure/vienna.py`'s `degradation_reason` and how `design/runner.py` consumes it.
3. **`database_version` recorded** whenever a real screen runs, so a verdict can be
   traced to what produced it.
4. **The `block` path honoured.** `BiosecurityVerdict.may_proceed` is
   `status != "block"`. A blocked design must not emit. `DesignContext` carries
   `strict_biosecurity: bool = True`.
5. **Protein validation** on the input, if scope allows — M8 owns it and the v1 flow is
   "protein → **validated** → screened → …".

## The failure this session exists to prevent

**A screen that reports "clear" when it never ran.** Every other bug here is
recoverable; that one hands a user a false assurance about a hazard. It is also
exactly the bug that passes every mechanical check — the types are right, the tests
are green, the field says "clear".

So: `not_run` is never upgraded by inference. Absent commec, absent a database, a
timeout, an exception — all of them are `not_run` or an error, never `clear`. Write
the test that would catch the opposite.

Related, from `CLAUDE.md` §9 and `.claude/rules/vector.md`: never add a "minimize
identity to a reference sequence" objective, and never let `KmerIndex` accept an
external database. Constraining the index to the assembled construct is the only thing
keeping BT5 from being a general-purpose screening-evasion tool. Your lane is the
*positive* screen — flagging hazards before an order goes out. Keep it that way.

## Files

**You own:** `packages/engine/src/bt5/cassette/**`,
`packages/engine/tests/cassette/**`.

**Never touch:** `design/`, `score/`, `rules/`, `solver/`, `vector/`, `codon/`,
`structure/`, `core/`, `verify.py`, `.github/`, `pyproject.toml`, `data/`,
`tests/contract/`, `tests/invariants/`, `tests/data_integrity/`.

`cassette/envelope.py` (the `FeasibleEnvelope` from issue #45's X1) is yours, but that
issue's remaining work — X2, X3, X4, X5, X7 — is **out of scope**. Leave it alone.

If you need a `core/` change, **stop**: open an issue naming the type and use
`/architect`. `core/` is a global mutex.

## Your contract with the other five

- **`BiosecurityVerdict`'s shape is frozen.** You change behaviour, not the type. S1 is
  wiring the report against `status`, `database_version` and `detail` right now.
- Post a **short design note as the first comment on your own draft PR** before you
  implement much — which protocol shape, where the degradation decision lives, what a
  timeout counts as. That is your review gate; it replaces plan mode without blocking
  the other five sessions.

## Delegation

- **`security-reviewer` fires on your diff by construction.** `/pre-pr` runs it iff the
  diff touches `vector/`, `core/services.py`, `verify.py`, `cassette/` or `.github/` —
  and `cassette/` is your whole lane. It is opus/xhigh and it judges *intent*, not
  signatures: "the change that passes every mechanical check and still weakens the
  posture." Write for that reader from the start.
- `docs-miner` — the brief's biosecurity section and PLAN's screening flow. Never read
  `brief.md` (63 KB) or `PLAN.md` (58 KB) inline.
- `Explore` — how `FoldEngine` degrades, and every construction of
  `BiosecurityVerdict` in the tree.
- `gate-runner` — gates before each push.
- **`/escalate`** rather than another attempt if the degradation design turns out
  wrong. This lane's defining risk is a capability failure, not a diligence one, and
  re-reading at higher effort does not fix those.

No file in your lane exceeds 20 KB.

## Done means

- `bash scripts/gates.sh` reaches `ALL GATES PASSED`. Exit 10 = run `/bootstrap`;
  pytest exits 2, 3, 4, 5 are BROKEN, and 5 is never success.
- With commec absent — the CI case — the verdict is `not_run` with a reason, and a
  test asserts it is **not** `clear`.
- With commec present, a real verdict carries a `database_version`.
- A `block` verdict refuses to emit.
- `security-reviewer` returns clean via `/pre-pr`, and the PR is **open as a draft**.
- You added a decision file at `docs/decisions/2026-XX-XX-<slug>.md`: decided,
  **rejected** and why, with evidence. One file per decision — never append to a
  shared one.

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

**Do not self-merge.** This changes what the app refuses to build, which under
`CLAUDE.md` §7b goes to the owner. Say in the PR that you know this.
