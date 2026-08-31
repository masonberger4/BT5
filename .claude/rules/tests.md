---
paths:
  - "tests/**"
  - "packages/engine/tests/**"
---

# Tests

Root `CLAUDE.md` §4 is the rule: suppression is not a fix. This is the rationale, plus
the list of gates in this repo that look load-bearing and are not.

## Why suppression is banned rather than discouraged

A skipped test does not record that the behaviour is broken — it records nothing, and it
looks the same as a test that passes. If a Hypothesis property fails on your branch and
also reproduces on the merge base, it is a pre-existing bug: record it as a fixture under
`tests/data/regressions/`, open an issue, say so in the PR, and let the owner merge. Do
not weaken the property to get green.

`--snapshot-update` is not a fix either, for the same reason a regenerated contract
fixture is not: re-recording a value does not make the old caller work.

## Gates that do not do what they appear to

Do not build a workflow, an agent instruction or a claim of CI parity on any of these:

- **`HYPOTHESIS_PROFILE` is inert.** `tests/conftest.py` registers `ci` (200 examples),
  `dev` (50) and `nightly` (2000), then calls `settings.load_profile("dev")`
  unconditionally and never reads the environment variable. CI sets
  `HYPOTHESIS_PROFILE: ci` on the `invariants` job and it has no effect: **every run,
  local and CI, is 50 examples.** The one upside is that the two sides agree.
- **`-p no:randomly` is a silent no-op** — `pytest-randomly` is not in the `dev` extra,
  and `-p no:<uninstalled>` exits 0 without warning.
- **`-m "not slow"` deselects nothing** — the `slow` marker is registered in
  `pyproject.toml` and applied zero times in the tree.
- **`fail_under = 85` has never been evaluated** — `pytest-cov` is installed but no job
  passes `--cov`.
- **`goldens-not-hand-edited` does not exist.** `tests/goldens/` holds only `.gitkeep`
  and syrupy is declared with no snapshot fixtures anywhere.
- **The RNG grep is narrower than §3.7 implies.** CI greps `packages/engine/src/` only,
  and its alternation covers `seed|rand|randn|choice|randint` — so `np.random.random(`,
  `.shuffle(`, `.permutation(`, `.normal(` and `.uniform(` pass it. Seed every RNG
  explicitly with `np.random.default_rng(seed)` because it is correct, not because the
  grep would catch you.
- **mypy is in no CI job.** It is `strict` over `packages/engine/src/bt5` and is the only
  thing that checks `kmers.py:461`'s `KmerIndex` conformance assertion. Running it
  locally is not optional.

## Environment

`pytest` must be `.venv/bin/pytest`. The bare name resolves to `/root/.local/bin`, whose
interpreter has no numpy: it exits **4** on a `conftest.py` import error, which looks
nothing like an environment problem. Exit codes 2, 3, 4 and 5 all mean BROKEN, not
FAILED — 5 in particular ("no tests collected") must never be read as success.

CI sets `PYTHONHASHSEED=0`; the project settings set it locally too, so a set- or
dict-ordering bug reproduces on both sides rather than one.

`pyproject.toml` promotes `Bio.BiopythonWarning` to an error on purpose: Biopython's
`Seq.translate()` silently truncates a non-multiple-of-3 sequence, so a frame-length bug
would pass a naive round-trip test.

## Layout

Root `tests/` holds the cross-cutting gates — `contract/`, `data_integrity/`,
`invariants/`. Per-lane unit tests live in `packages/engine/tests/<lane>/`, mirroring the
source tree. Every catalog rule needs a paired `packages/engine/tests/rules/test_<id>.py`;
`d1_restriction_sites` is the one pre-existing gap.
