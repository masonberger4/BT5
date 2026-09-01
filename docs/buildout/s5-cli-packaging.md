# S5 — CLI and packaging (M7)

*Copy this whole file into a fresh Claude Code session on the BT5 repo.*

**Launch:** permission mode `default` · **model sonnet** (override; the repo default is
opus) · effort high · **needs a person within reach** — your job is editing
`pyproject.toml`, and that protected-path prompt needs someone to answer it.
**Do not run this in plan mode** — five other sessions are running in parallel and the
gate is the draft PR, not a plan approval.

---

You are giving BT5 a way to be used. Today a user has to write Python.

## Read this first

Run `/bootstrap` — a fresh checkout has no `.venv`, and `gates.sh` exit **10** means
BROKEN, not a code failure. Then `CLAUDE.md`, then `docs/buildout/README.md` for the
un-draft queue and your inter-session contracts.

Your branch: **`claude/s5-cli-packaging`**. Cut it yourself, from a
**freshly fetched** main — do not reuse a branch someone made earlier, and do not
assume main is where it was when this prompt was written:

```bash
git remote prune origin
git fetch -q origin main && git checkout -B claude/s5-cli-packaging origin/main
```

`git fetch --prune origin main` does **not** clear a stale ref — pruning is bounded by
the refspec, so that form prunes only `origin/main` (CLAUDE.md §7a).

## The situation

There is no `__main__.py`, no `[project.scripts]`, no CLI anywhere. The only way to run
BT5 end to end is:

```python
from bt5.design import design
from bt5.vector import read_genbank

result = design(
    backbone=read_genbank("plasmid.gb"),
    protein="MKTL...",
    table_id=11,
    modality=...,
    hosts=("human",),
    seed=42,
)
```

`docs/PLAN.md` never specifies a CLI — it goes straight to M9 (`packages/server/`,
FastAPI) and M10 (`apps/web/`), both **out of scope** for this buildout. The CLI is the
stand-in user surface, so you are designing it, not transcribing a spec.

## What to build

1. **`bt5 design`**, over `bt5.design.design()`. Take a backbone file and a protein,
   emit the annotated GenBank and the order CSV. `table_id` is **required and never
   defaulted** — `CLAUDE.md` §3.1, and the function already enforces it.
2. **`[project.scripts]`** in `pyproject.toml` wiring the entry point.
3. **`[tool.mypy] files`** extended. It currently reads
   `files = ["packages/engine/src/bt5"]`. Anything you add outside that path is
   **not typechecked**, and `mypy` is a required CI job as of #63 — so a new top-level
   package silently escapes the gate unless you extend this.
4. **The hatch wheel target**, if you add a package outside
   `packages/engine/src/bt5`. It currently reads `packages = ["packages/engine/src/bt5"]`.
5. **`packaging/`** — lane zero per `docs/PLAN.md:456`: the uv install flow, ViennaRNA,
   and the commec database flow. Document it; the Tauri sidecar and macOS notarization
   are later and out of scope.

**Do not add a dependency.** Every one is already declared — `server`, `screen`,
`export` and `fold` extras included. `CLAUDE.md` §5: a lockfile conflict across
parallel PRs is the single most expensive merge failure in this repo. If you think you
need one, open an issue instead.

## You hold the `pyproject.toml` mutex

You are the **only** session permitted to write `pyproject.toml`. Five others are
running against the same `main` and a conflict there is the expensive one.

`.claude/hooks/protect_paths.py` will prompt on every edit to it — the hook notes that
`check-approval-labels.sh` says of `pyproject.toml` and `uv.lock`, verbatim, *"protected
by CLAUDE.md section 2 but no label is named for them, so they are deliberately NOT
enforced here."* There is no `approved:*` label for this file. That means a human
judgement, which is why you run in `default` mode.

Keep your `pyproject.toml` diff **as small as it can possibly be**, and put it in its
own commit so it is trivially reviewable.

## Files

**You own:** `packages/engine/src/bt5/cli.py` (and a `__main__.py` if you want
`python -m bt5`), `packages/engine/tests/cli/**`, `packaging/**`, and
**`pyproject.toml`**.

**Never touch:** `design/`, `score/` (S1's), `cassette/` (S2's), `rules/` (S3's and
S4's), `solver/`, `vector/`, `codon/`, `structure/`, `core/`, `verify.py`, `.github/`,
`data/` (S6's), `tests/contract/`, `tests/invariants/`, `tests/data_integrity/`.

You are a **consumer** of every engine lane. If `design()` does not expose something
your CLI needs, that is an issue for S1, not an edit you make.

## Your contract with the other five

- **`design()`'s signature is frozen** for the duration
  (`design/runner.py:156-171` — keyword-only, `table_id` never defaulted). Build
  against it.
- **S1 is enriching `SkeletonResult` right now** — adding a gallery, percentiles,
  `native_baseline` and an order CSV. It has agreed to **add** fields, never remove or
  rename one. So build your CLI against today's fields and expect new ones to appear;
  do not block on S1, and do not design around fields that do not exist yet.
- **Extending `[tool.mypy] files` changes what CI typechecks for every other lane.**
  That is cross-lane even though only you write the file. Run **`/architect`** on that
  decision before making it.

## Delegation

This is the least judgment-dense session in the buildout — which is why it runs
sonnet. Lean on delegation and keep the main thread for the two real decisions (the
CLI's shape, and the mypy/hatch change).

- `Explore` — `design()`'s full surface, `SkeletonResult`'s fields, and how
  `read_genbank` / `write_genbank` are called today.
- `gate-runner` — gates before each push. `mypy` is a required CI job; run it locally
  rather than discovering it in CI.
- **`/architect`** — the `[tool.mypy] files` and hatch `packages` change, once.
- **`/cheap-pass`** — if argument plumbing turns into ≥5 identical edits with the
  before/after already decided.
- `docs-miner` — `docs/PLAN.md:456` on what M7 packaging owns. Never read `PLAN.md`
  (58 KB) inline.

**Over the 20 KB bare-Read limit**, if you stray: `vector/backbone.py` (28.9 KB),
`vector/kmers.py` (21.2), `vector/assemble.py` (18.7). Use `offset`/`limit` or
`Explore`.

## Done means

- `bash scripts/gates.sh` reaches `ALL GATES PASSED`, **including `mypy`** on whatever
  you added. Exit 10 = `/bootstrap`; pytest exits 2, 3, 4, 5 are BROKEN, and 5 is never
  success.
- `bt5 design` runs end to end from a shell on the synthetic fixture at
  `tests/data/backbones/synthetic_mcs_ef1a.gb` and writes a GenBank.
- The `pyproject.toml` diff is minimal and in its own commit.
- Nothing you added escapes `mypy`.
- The PR is **open as a draft**.
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

**Do not self-merge.** `pyproject.toml` has no `approved:*` label precisely because it
is a human call, and a new user-facing surface is a product decision. Goes to the owner.
