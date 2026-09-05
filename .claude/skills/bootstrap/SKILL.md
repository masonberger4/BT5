---
name: bootstrap
description: Create the local .venv and install the dev, fold and export extras.
disable-model-invocation: true
allowed-tools: Bash
---

# /bootstrap

```bash
uv venv --python 3.11 .venv
uv pip install --python .venv/bin/python -e ".[dev,fold,export]"
```

**All three extras, always.** `[fold]` pins `viennarna==2.7.2`; without it, local `mypy`
reports `Unused "type: ignore"` on the ViennaRNA shims — errors that are correct only
because the dependency is missing. Installing a subset produces a type-check result you
cannot trust.

**This takes minutes.** Never run it from a hook, and never run it unprompted on a user's
behalf mid-task: `uv pip install` is a state change, and CLAUDE.md §2 makes dependency
events a human decision.

## Verify

```bash
.venv/bin/python -c "import numpy, Bio, hypothesis, bt5; print('ok')"
```

After this, **every** command uses the venv binaries — `.venv/bin/pytest`,
`.venv/bin/ruff`, `.venv/bin/mypy`, `.venv/bin/python`. Never the bare names: those
resolve to `/root/.local/bin`, whose interpreter has no numpy, so bare `pytest` exits 4
on a `conftest.py` import error and bare `mypy` emits phantom errors.

`uv.lock` is not checked in; `uv pip install -e` resolves from `pyproject.toml`. Do not
create a lockfile (CLAUDE.md §2).
