# BT5 — contract for AI coding sessions

Read this before touching anything. `docs/PLAN.md` has the design and
`docs/research/brief.md` the science — reach both via `docs-miner`, never inline.

## 0. What this project is

BT5 back-translates proteins and codon-optimizes DNA **in the context of the
assembled construct**, balancing protein expression, DNA synthesizability, viral titer
and plasmid stability, then exports an annotated GenBank and a vendor order file.

**What it refuses to claim.** All computable design features together explain only
5–31% (mean ~14%) of protein-level variance, and nine benchmarked commercial
optimizers were a coin flip against native sequence. So BT5 **never reports a
predicted expression number**. It reports ranks, percentiles against a
random-synonymous null, and confidence bands. `native_baseline` — "don't optimize"
— is a first-class output. A CI gate bans prediction vocabulary from the schema.

## Stack

Python 3.11 only. One package: `packages/engine/src/bt5` (src-layout), root
`pyproject.toml`. `apps/web/`, `packages/server/`, `packaging/` and `benchmarks/`
appear in the plan but **do not exist on disk**. **Every command uses `.venv/bin/…`,
never the bare name** — bare names resolve to an interpreter with no numpy, so bare
`pytest` exits 4 on a `conftest.py` import error, which looks like a code failure and
is not one. Fresh checkouts have no `.venv`: run `/bootstrap`.

## 1. Module ownership — stay in your lane

| Lane | You own | You must NOT edit |
|---|---|---|
| M1 solver | `packages/engine/src/bt5/solver/` | any rule file |
| M2 vector | `packages/engine/src/bt5/vector/` | the solver |
| M3 score | `packages/engine/src/bt5/score/` | rules or solver |
| M4 rules | `packages/engine/src/bt5/rules/catalog/` | the solver, the oracle |
| M5 codon | `packages/engine/src/bt5/codon/` | rules |
| M6 structure | `packages/engine/src/bt5/structure/` | rules |
| M7 packaging | `packaging/` *(planned, no files)* | engine source |
| M8 cassette | `packages/engine/src/bt5/cassette/` | the solver |
| M9 server | `packages/server/` *(planned, no files)* | `apps/web/` |
| M10 web | `apps/web/` *(planned, no files)* | the engine |

Cross-lane changes need an issue first. Blocked on another lane? Code against the
protocol in `bt5/core/` and a recorded fixture — never reach into their directory.

## 2. Files you must NEVER modify without the matching label

- `packages/engine/src/bt5/verify.py` — the oracle (`approved:oracle-change`)
- `packages/engine/src/bt5/core/**` — the frozen contract (`approved:contract-change`)
- `tests/contract/**` — the RECORD of that contract (`approved:contract-change`)
- `tests/invariants/**`, `tests/data_integrity/**` (`approved:oracle-change`)
- `benchmarks/baseline.json`, `benchmarks/tolerances.yaml` (`approved:algorithm-change`)
- `data/genetic_codes/**`, `data/codon_usage/**` (`approved:data-change`)
- `.github/**` (`approved:ci-change`)
- `pyproject.toml` / `uv.lock` — every dependency is already declared

## 2a. `core/` is frozen

`contract-freeze` classifies your branch against `main` by asking one question:
**who breaks?**

- **MINOR** — a new type, a new **defaulted** field, a new enum member, a field that
  gains a default. Regenerate, commit the manifest and fixtures, done.
- **MAJOR** — a removed or renamed anything, a changed annotation or default, a field
  that **loses** its default, a changed signature, a new protocol method. Needs an RFC,
  a deprecation shim, the two-window rule, and `pytest tests/contract` passing
  **without regenerating the fixtures**.

Use `/contract-change`: classification must precede regeneration, because
`regenerate.py` writes first and returns 0 on every path. Detail in
`.claude/rules/contract-core.md` and `docs/rfcs/README.md`.

## 3. Correctness rules that are not negotiable

1. **The genetic code table is explicit and never defaulted.** NCBI table 12
   reassigns CTG to Ser rather than Leu; table 4 makes TGA Trp. A wrong table is a
   silently wrong protein no assay catches for months.
2. **Never emit a codon that is also a stop codon in the target table.** Tables 27
   and 28 make TGA both Trp and a stop.
3. **Never evaluate a rule against a bare string.** Rules take a `Construct` — that is
   what makes junction-, origin-spanning and reverse-strand hits impossible to miss.
4. **Never scan the reverse strand yourself for motif rules.** List forward motifs in
   `LatticeTerms.forbidden`; the solver closes the set under reverse complement.
   Directional scored models are NOT revcomp-symmetric and must read
   `slot.strand_of_interest` — hard-coding strand 1 makes a reverse-oriented lentiviral
   cassette's polyA and splice analysis exactly backwards.
5. **Hard constraints are never enforced by a penalty weight.** Use `HARD_LATTICE`
   (guaranteed by the automaton), `HARD_REPAIR` (repair plus the independent validator,
   which refuses to emit) or `HARD_CHECK` (real but unfixable by codon choice).
   `default_weight` must be 0.0 for all three; `steering_weight` nudges the DP.
6. **Splice-site removal must use `RepairPolicy.FIXED_POINT`.** A single pass ships a
   construct whose donors were removed *into* new donors, and the validator passes it
   because the specific 9-mer is gone.
7. **Seed every RNG explicitly.** `np.random.default_rng(seed)`. Global
   `np.random.seed` / `random.seed` are banned in `src/`. The CI grep is narrower than
   this rule — do not rely on it to catch you.

## 4. Suppression is not a fix

Never skip, disable, `xfail` or loosen a test to get green. Never weaken a Hypothesis
property. If a property fails on your PR and reproduces on the merge base, it is a
pre-existing bug: file it as a fixture under `tests/data/regressions/` plus an issue, and
say so in the PR. The owner merges.

`--snapshot-update` is not a fix either, and neither is re-recording a contract
fixture. *(`goldens-not-hand-edited` is named in the plan but does not exist yet;
`tests/goldens/` holds only `.gitkeep`.)*

## 5. Dependencies and lockfiles

Every dependency is declared in `pyproject.toml`. **Do not add one.** If you genuinely
need something new, open an issue — a lockfile conflict across parallel PRs is the
single most expensive merge failure in this repo. If you hit a conflict: rebase and
regenerate, never hand-merge.

## 6. ViennaRNA is pinned deliberately

A version bump is a **scientific change**, not a dependency bump: it carries
`approved:algorithm-change` and regenerates the baseline. Never put a ΔG in a
byte-exact snapshot. Why: `.claude/rules/tests.md`.

## 7. Branching and PRs

Branch from `main`. One lane per PR. Squash merge only. Fill in the PR template,
including the "scientific impact" section — say what changed about the sequences the app
produces, not just the code.

Open your PR as a **draft** until you believe it is done; drafts skip the expensive CI
jobs, and CI capacity is the binding constraint (20 concurrent job slots, ~12 per Python
PR, so at most 5 open non-draft PRs at a time).

## 7a. Merged branches leave a stale ref behind

`git fetch --prune origin main` does **not** fix it: pruning is bounded by the refspec.
Prune the whole remote instead:

```bash
git remote prune origin            # or: git fetch --prune origin   (no refspec)
git fetch -q origin main
git checkout -B <your-branch> origin/main
```

## 8. Before you push

Run `/pre-pr`, or `bash scripts/gates.sh` directly — ruff, ruff format, mypy,
invariants, data_integrity, contract and engine tests, each independently. One validated
push beats three speculative ones.

Exit **10** = no usable venv (BROKEN, run `/bootstrap`). pytest exits **2, 3, 4, 5** are
BROKEN too — 5 is "no tests collected" and is never success. **mypy is local-only**: no
CI job runs it, so this is the only thing between a strict-mode error and `main`.

## 9. Never

- Add a `paths:` filter to `on:` in a workflow that owns a required check — a
  required check that never reports blocks the PR forever with no error.
- Add a CI job without adding it to `required-checks.needs`.
- Add a "minimize identity to a reference sequence" objective, or let `KmerIndex`
  accept an external database. BT5's output is already the textbook method for evading
  nucleotide-homology screening; constraining the index to the assembled construct is
  what keeps it from being a general-purpose evasion tool.
- Report a predicted expression level, titer, yield or fold-improvement.

## Delegation

Volume work goes to a subagent so its output never enters this window. Capability
failure raises the model; diligence failure raises the effort.

- **Anything in `docs/`** → `docs-miner`. Never read `brief.md` (63 KB), `PLAN.md`
  (58 KB), `github-setup.md` (121 KB) or `design-review-verdicts.json` (107 KB) inline.
- **Finding code** → `Explore`. **Running the gates** → `gate-runner`.
- **One identical edit across ≥5 files** → `/cheap-pass` (forks `batch-editor`).
- **Before a PR** → `/pre-pr`, which chains `gate-runner` → `code-reviewer`, adding
  `rule-auditor` and `security-reviewer` only when the diff warrants them.
- **A failure a first-pass fix already missed** → `debugger`. Not for first-pass bugs.
- **Cross-lane, `core/`, or a MINOR/MAJOR call** → `/architect`; **stuck on a decision**
  → `/escalate`. Both run opus at ultracode and orchestrate rather than answer alone.
- **Never bare-Read a source file over 20 KB** (`vector/backbone.py` 28.9 KB,
  `rules/vendors.py` 25.7 KB, `vector/kmers.py` 21.2 KB, `catalog/e4_gc_extent.py`
  19.1 KB, `vector/assemble.py` 18.7 KB) — use `offset`/`limit`, or delegate.
- **Append settled decisions to `docs/decisions.md`** at the end of each work slice:
  what survives compaction is what lives on disk.

## Compact instructions

A compaction summary must preserve:

1. The current task and its acceptance criteria.
2. The lane, the branch, and any `approved:*` label the change will need.
3. Every touched file path.
4. The failing gate, the failing node ids, and the **exact error text**, not a paraphrase.
5. Decisions made **and decisions rejected**, each with its reason.
6. Any MINOR/MAJOR classification reached, and whether `regenerate.py` has run.
7. Rule ids added or changed.

Discard raw command output, directory listings and `docs/` quotations already summarised.
On resume: check `.venv/bin/python`, run `bash scripts/gates.sh`, then `git diff --name-only`.
