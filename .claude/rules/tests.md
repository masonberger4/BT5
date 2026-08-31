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

- **`HYPOTHESIS_PROFILE` is honored, as of #63.** `tests/conftest.py` reads it and
  raises on an unregistered name; `pytest_report_header` prints the profile and
  `max_examples` so a silent budget regression is visible on the first run. Locally you
  get `dev` (50 examples); CI's `invariants` job sets `ci` (200). So **local is a
  quarter of CI's budget** — a property that passes locally has had far less search.
- **`-p no:randomly` is a silent no-op** — `pytest-randomly` is not in the `dev` extra,
  and `-p no:<uninstalled>` exits 0 without warning.
- **`-m "not slow"` deselects nothing** — the `slow` marker is registered in
  `pyproject.toml` and applied zero times in the tree.
- **`fail_under = 85` has never been evaluated** — `pytest-cov` is installed but no job
  passes `--cov`.
- **`goldens-not-hand-edited` does not exist.** `tests/goldens/` holds only `.gitkeep`
  and syrupy is declared with no snapshot fixtures anywhere.
- **The RNG grep is broad, as of #63.** It matches any `np.random.*` attribute call
  except the explicit-generator constructors (`default_rng`, `Generator`, `SeedSequence`,
  `PCG64`), and separately bans importing stdlib `random` at all — across engine source
  **and** the test tree, because an unseeded draw in a strategy is just as
  irreproducible. Seed with `np.random.default_rng(seed)`.
- **mypy is a required CI job, as of #63**, in `required-checks.needs`. It installs the
  `fold` extra, unlike `python-quality`, because `structure/vienna.py`'s
  `# type: ignore[import-untyped]` does not suppress `import-not-found` when ViennaRNA is
  absent. Run it locally too — it is the only check on `kmers.py:461`'s conformance
  assertion, which lives under `if TYPE_CHECKING`.

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
