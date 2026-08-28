# BT5 — contract for AI coding sessions

Read this before touching anything. Read `docs/PLAN.md` for the full design and
`docs/research/brief.md` for the science behind any rule you implement.

## 0. What this project is

BT5 back-translates proteins and codon-optimizes DNA **in the context of the
assembled construct**, balancing protein expression, DNA synthesizability, viral
titer and plasmid stability, then exports an annotated GenBank and a vendor order
file.

**What it refuses to claim.** All computable design features together explain only
5–31% (mean ~14%) of protein-level variance, and nine benchmarked commercial
optimizers were a coin flip against native sequence. So BT5 **never reports a
predicted expression number**. It reports ranks, percentiles against a
random-synonymous null, and confidence bands. `native_baseline` — "don't optimize"
— is a first-class output. A CI gate bans prediction vocabulary from the schema.

## 1. Module ownership — stay in your lane

| Lane | You own | You must NOT edit |
|---|---|---|
| M1 solver | `bt5/solver/` | any rule file |
| M2 vector | `bt5/vector/` | the solver |
| M3 score | `bt5/score/` | rules or solver |
| M4 rules | `bt5/rules/catalog/` | the solver, the oracle |
| M5 codon | `bt5/codon/` | rules |
| M6 structure | `bt5/structure/` | rules |
| M7 packaging | `packaging/` | engine source |
| M8 cassette | `bt5/cassette/` | the solver |
| M9 server | `packages/server/` | `apps/web/` |
| M10 web | `apps/web/` | the engine |

Cross-lane changes need an issue first. If you are blocked on another lane, code
against the protocol in `bt5/core/` and a recorded fixture — never reach into
their directory.

## 2. Files you must NEVER modify without the matching label

- `packages/engine/src/bt5/verify.py` — the oracle (`approved:oracle-change`)
- `packages/engine/src/bt5/core/**` — the frozen contract (`approved:contract-change`)
- `tests/invariants/**`, `tests/data_integrity/**` (`approved:oracle-change`)
- `benchmarks/baseline.json`, `benchmarks/tolerances.yaml` (`approved:algorithm-change`)
- `data/genetic_codes/**`, `data/codon_usage/**` (`approved:data-change`)
- `.github/**` (`approved:ci-change`)
- `pyproject.toml` / `uv.lock` — every dependency is already declared

## 3. Correctness rules that are not negotiable

1. **The genetic code table is explicit and never defaulted.** NCBI table 12
   reassigns CTG to Ser rather than Leu; table 4 makes TGA Trp. A wrong table is a
   silently wrong protein no assay catches for months.
2. **Never emit a codon that is also a stop codon in the target table.** Tables 27
   and 28 make TGA both Trp and a stop.
3. **Never evaluate a rule against a bare string.** Rules take a `Construct`. That
   is what makes junction-spanning, origin-spanning and reverse-strand hits
   impossible to miss rather than something you remembered to handle.
4. **Never scan the reverse strand yourself for motif rules.** List forward motifs
   in `LatticeTerms.forbidden`; the solver closes the set under reverse complement.
   Directional scored models (MaxEntScan, Salis TIR, promoter calculator, polyA
   downstream element) are NOT revcomp-symmetric — they must read
   `slot.strand_of_interest`. Hard-coding strand 1 makes a reverse-oriented
   lentiviral cassette's polyA and splice analysis exactly backwards.
5. **Hard constraints are never enforced by a penalty weight.** Use
   `HARD_LATTICE` (guaranteed by the automaton), `HARD_REPAIR` (repair plus the
   independent validator, which refuses to emit), or `HARD_CHECK` (real but not
   fixable by codon choice). `default_weight` must be 0.0 for all three; use
   `steering_weight` if the DP needs nudging.
6. **Splice-site removal must use `RepairPolicy.FIXED_POINT`.** Point-mutating one
   cryptic donor activates cryptic donors nearby; a single pass ships a construct
   whose donors were removed into new donors, and the validator passes it because
   the specific 9-mer is gone.
7. **Seed every RNG explicitly.** `np.random.default_rng(seed)`. Global
   `np.random.seed` / `random.seed` are banned in `src/` and CI greps for them.

## 4. Suppression is not a fix

Never skip, disable, `xfail` or loosen a test to get green. Never weaken a
Hypothesis property. If a property fails on your PR and reproduces on the merge
base, it is a pre-existing bug: file it as a fixture under
`tests/data/regressions/` plus an issue, and say so in the PR. The owner merges.

`--snapshot-update` is not a fix either. `goldens-not-hand-edited` regenerates
every snapshot from scratch and diffs, so a hand-edited golden fails.

## 5. Dependencies and lockfiles

Every dependency is declared in `pyproject.toml` in PR #0. **Do not add one.** If
you genuinely need something new, open an issue — a lockfile conflict across
parallel PRs is the single most expensive merge failure in this repo, and a
hand-resolved lockfile satisfies neither `--locked` nor reality. If you hit a
conflict: rebase and regenerate, never hand-merge.

## 6. ViennaRNA is pinned deliberately

Energy parameters determine every ΔG, so an unannounced version change silently
shifts every benchmark baseline and makes old results incomparable. A version bump
is a **scientific change**, not a dependency bump: it carries
`approved:algorithm-change` and regenerates the baseline. Never put a ΔG in a
byte-exact snapshot.

## 7. Branching and PRs

Branch from `main`. One lane per PR. Squash merge only. Fill in the PR template,
including the "scientific impact" section — say what changed about the sequences
the app produces, not just the code.

Open your PR as a **draft** until you believe it is done; drafts skip the
expensive CI jobs, and CI capacity is the binding constraint (20 concurrent job
slots, ~12 per Python PR, so at most 5 open non-draft PRs at a time).

## 7a. Merged branches leave a stale ref behind

Head branches are auto-deleted on merge, so after every merge the local
remote-tracking ref for your branch points at a commit that no longer exists on
the remote. Anything comparing HEAD against it then reports phantom unpushed
commits. Always fetch with `--prune`:

```bash
git fetch --prune origin main
git checkout -B <your-branch> origin/main
```

`git config remote.origin.prune true` makes every fetch do it, but the setting
is per-clone and does not survive a fresh checkout, so the flag is the reliable
form.

## 8. Before you push

Run these locally. One validated push beats three speculative ones.

```bash
. .venv/bin/activate
ruff check . && ruff format --check .
mypy
pytest tests packages/engine/tests -q
```

## 9. Never

- Add a `paths:` filter to `on:` in a workflow that owns a required check — a
  required check that never reports blocks the PR forever with no error.
- Add a CI job without adding it to `required-checks.needs`.
- Add a "minimize identity to a reference sequence" objective, or let
  `KmerIndex` accept an external database. BT5's output is already the textbook
  method for evading nucleotide-homology screening; constraining the index to the
  assembled construct is what keeps it from being a general-purpose evasion tool.
- Report a predicted expression level, titer, yield or fold-improvement.
