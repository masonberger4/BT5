# BT5 — contract for AI coding sessions

Design: `docs/PLAN.md`. Science: `docs/research/brief.md`. Reach both via `docs-miner`,
never inline; rationale that isn't here lives in `.claude/rules/` (per-lane) and
`docs/decisions/` (settled calls, incl. rejected ones). Python 3.11, one package
in src-layout at `packages/engine/src/bt5`. **Every command uses `.venv/bin/…`, never the
bare name.** No `.venv`? Run `/bootstrap`.

## 0. What this project is

BT5 back-translates proteins and codon-optimizes DNA **in the context of the assembled
construct**, then exports an annotated GenBank and a vendor order file. Computable design
features explain ~14% of protein-level variance, so BT5 **never reports a predicted
expression number** — only ranks, percentiles against a random-synonymous null, and bands.
`native_baseline` ("don't optimize") is a first-class output; a CI gate bans prediction
vocabulary from the schema.

## 1. Stay in your lane

One lane per PR; you own a **declared set of paths** and edit no other. Under `bt5/`: M1
`solver/`, M2 `vector/`, M3 `score/`, M4 `rules/catalog/`, M5 `codon/`, M6 `structure/`,
M8 `cassette/`, M11 `design/`. M7 is `packaging/` (repo root, not under `bt5/`) plus
`bt5/cli.py`. M9/M10 (`packages/server/`, `apps/web/`) are planned. A lane is usually one
directory but not always — the live matrix is `docs/buildout/wave2/README.md`, where W2 owns
`design/` **and** `score/`; take your path set from there, not from the count. Rules never
touch the solver or the oracle; nothing but M1 touches the solver. Cross-lane changes need an
issue first; blocked on another lane, code against the protocol in `bt5/core/` and a recorded
fixture — never reach into their directory.

## 2. Protected paths

Each needs its `approved:*` label: `verify.py` and `tests/{invariants,data_integrity}/**` →
`oracle-change`; `core/**` and `tests/contract/**` → `contract-change`; `benchmarks/`
baselines → `algorithm-change`; `data/{genetic_codes,codon_usage}/**` → `data-change`;
`.github/**` → `ci-change`. The edit hook says so too, but read-only reviewers never fire it.
`pyproject.toml`/`uv.lock`: **never add a dependency**; on a lockfile conflict rebase and
regenerate — never hand-merge it.

`core/` is frozen, and `tests/contract/surface.py:22-25` is the classifier. **MINOR** = a new
type, a defaulted field, an enum member, a field GAINING a default. **MAJOR** = a removal,
rename, changed annotation/signature, a changed default *value*, a field LOSING its default,
a new protocol method — needs an RFC and a deprecation shim. `/contract-change`; classify FIRST.

`.claude/agent-memory/**` carries no label but the edit hook asks anyway: it is an instruction
channel that ships inside the PR diff, so read any change to it before you accept it.

## 3. Rules that are not negotiable

1. **The genetic code table is explicit, never defaulted** — a wrong table is a silently
   wrong protein no assay catches for months.
2. **Never emit a codon that is also a stop codon in the target table.**
3. **Never evaluate a rule against a bare string.** Rules take a `Construct`; that is what
   makes junction-, origin-spanning and reverse-strand hits impossible to miss.
4. **Never scan the reverse strand yourself for motif rules.** List forward motifs in
   `LatticeTerms.forbidden` and let the solver close the set. Directional scored models are
   NOT revcomp-symmetric — they must read `slot.strand_of_interest`.
5. **Hard constraints are never a penalty weight.** Use `HARD_LATTICE`, `HARD_REPAIR` or
   `HARD_CHECK`, each with `default_weight` 0.0. `steering_weight` is a declared channel the
   solver does not read yet — today's only Tier-A steering is `solve_with_gc_steering`.
6. **Splice-site removal must use `RepairPolicy.FIXED_POINT`.** A single pass ships a
   construct whose donors were removed *into* new donors, and the validator passes it.
7. **Seed every RNG explicitly** with `np.random.default_rng(seed)`. Global `np.random.*`
   and stdlib `random` are banned in engine source and tests; CI greps for both.
8. **ViennaRNA is pinned.** A bump is a scientific change: `approved:algorithm-change` plus a
   baseline regeneration. Never put a ΔG in a byte-exact snapshot.
9. **Suppression is not a fix.** Never skip, disable, `xfail` or loosen a test or weaken a
   Hypothesis property; `--snapshot-update` and re-recording a contract fixture are no better.
   A property failing on your PR **and** on the merge base is a pre-existing bug: record a
   fixture under `tests/data/regressions/`, open an issue, say so in the PR — **owner** merges.
10. **Never give a workflow that owns a required check a `paths:` filter**, and never add a
   CI job without adding it to `required-checks.needs` — or, if it is skipped on
   `pull_request` (as `main-broken` is), to `NON_BLOCKING` in `check-workflow-gate.py`, since
   a `skipped` dependency counts as failure. A required check that never reports blocks the
   PR forever, with no error anywhere.
11. **Never add a "minimize identity to a reference sequence" objective**, or let `KmerIndex`
   accept an external database — constraining the index to the assembled construct is all
   that keeps BT5 from being a general-purpose screening-evasion tool.
12. **Never report a predicted expression level, titer, yield or fold-improvement.**

## 4. Branching, PRs and merging

Branch from `main`. Squash only. Fill in the PR template, including "scientific impact" —
what changed about the sequences the app produces, not just the code. Open as a **draft**
until done — it cannot be merged by mistake. A draft does **not** save CI: every job runs on
drafts too (`claude-review-gate`, the one job that skipped them, went on 2026-09-05), and CI
capacity binds (~5 non-draft PRs max). Before pushing run `/pre-pr`. After a merge:
`git remote prune origin` (NOT `--prune origin main`, which prunes only `origin/main`), then
`git fetch origin main && git checkout -B <branch> origin/main` — `remote prune` fetches
nothing, so without the fetch the new branch starts at the **pre-merge** commit.

**An agent may squash-merge its own PR once CI is green** — green means the ONE required
context on the PR's CURRENT head: `required-checks`. `pre-pr-attest` was a second required
context until 2026-09-05, when it and `claude-review-gate` went for spending Claude usage on
every push (`docs/decisions/2026-09-05-remove-claude-usage-ci-checks.md`). Nothing now proves
`/pre-pr` ran against the commit being merged, so the line above is discipline, not a gate,
and skipping it is invisible again. Green is necessary, not sufficient — four things go
to the owner instead: an `approved:*` label on a §2 path (the label signs off the change, not
the merge); a recorded pre-existing bug (§3.9); a non-"none" scientific impact or any change
to what the app REFUSES to build; an unresolved review thread. Say in the PR why it qualified.

## 5. Context management

Volume work goes to a subagent so its output never enters this window. Capability failure
raises the model; diligence failure raises the effort. `docs/{research,design}/**` and
`PLAN.md` → `docs-miner` (read `decisions/` and `buildout/` directly) · code →
`Explore` · gates → `gate-runner` · ≥5 identical edits → `/cheap-pass` · a missed fix →
`debugger` · before a PR → `/pre-pr` · cross-lane, `core/`, or a MINOR/MAJOR call →
`/architect` · stuck → `/escalate`. **Never bare-Read a source file over 20 KB** — the big
ones are in `rules/`, `design/`, `solver/`, `vector/`; use `offset`/`limit`, or delegate.
**Append settled decisions** to `docs/decisions/` each work slice as one new
`YYYY-MM-DD-slug.md`: what survives compaction is what lives on disk, and one shared file
collides across sessions.

**Compaction must preserve:** the task and its acceptance criteria; lane, branch and
`approved:*` label; every touched file path; the failing gate, node ids and **exact error
text**; decisions made **and rejected**, with reasons; any MINOR/MAJOR call and whether
`regenerate.py` ran; rule ids touched. Discard raw output and `docs/` already summarised.
On resume: `.venv/bin/python`, `bash scripts/gates.sh`, `git diff --name-only`.
