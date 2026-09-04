# Packaging (M7)

Lane zero per `docs/PLAN.md:456`: the `uv` install flow, ViennaRNA, and the
`commec` biosecurity database flow. This buildout covers those three. The
Tauri desktop sidecar and macOS notarization are a later packaging change on
top of the same `uv`-installed core (`docs/PLAN.md:177-178`) and are out of
scope here.

## Installing BT5

BT5 is one src-layout package, `bt5`, built with hatchling. There is no
published package yet, so today's install is always from a checkout:

```bash
uv venv --python 3.11 .venv
uv pip install --python .venv/bin/python -e .
.venv/bin/bt5 --help
```

That gives you the `bt5` command and the runtime dependencies every design
needs (`numpy`, `biopython`, `pyahocorasick`, `pydantic`, `pyyaml`) — no
extras. Add extras for what you actually plan to do:

| Extra | Installs | Needed for |
|---|---|---|
| `fold` | `viennarna==2.7.2` | Folding-derived objectives (5' cap structure, hairpins). Without it, BT5 still runs and says so: see [ViennaRNA](#viennarna) below. |
| `export` | `openpyxl` | `bt5 design --out-order` when the target format is a vendor plate workbook (`bt5.score.order.write_idt_plate`); the plain `Name,Sequence` CSV the CLI writes today needs neither. |
| `screen` | *(nothing — currently commented out)* | Protein-level biosecurity screening. **Not installable today:** `commec` is not published on PyPI, so the extra was commented out in `pyproject.toml` under [#106](https://github.com/masonberger4/BT5/issues/106) — it made every whole-project `uv run`/`sync`/`lock` fail in the repo root. Nothing imports `commec` yet, so nothing is lost. See [commec database](#commec-database) below. |
| `server` | `fastapi`, `uvicorn`, `sse-starlette` | M9 (`packages/server/`) — not built yet; installing this extra today buys you nothing. |

Combine extras as needed, e.g. `uv pip install -e ".[fold,export]"`. Contributors
install `.[dev,fold,export]` — see root `README.md` and `/bootstrap` — because
`dev` pulls in the test toolchain (`pytest`, `hypothesis`, `ruff`, `mypy`, …)
that an end user of the `bt5` command does not need.

`uv` verified against this project: `uv 0.8.17`, `Python 3.11.15`
(`docs/PLAN.md:230`). `ViennaRNA` ships prebuilt wheels for the platforms `uv`
resolves against — no source build, no known platform gaps, per the same line.

## ViennaRNA

Pinned exactly (`viennarna==2.7.2`), not as a floor: the folding objectives'
literature thresholds — the Boël −39 kcal/mol dual gate, the cap-proximal
−30/−50/−60 ladder — were calibrated against ViennaRNA's energy parameters, so
BT5 applies them directly rather than transferring across engines
(`docs/PLAN.md:186-191`). A version bump is a scientific change, not a
dependency bump (`CLAUDE.md` section 6): it needs `approved:algorithm-change`
and a rebuilt benchmark baseline, never a routine Renovate merge.

Without the `fold` extra installed, `bt5 design` still completes — it does not
refuse to run — but every folding-derived objective reports `unavailable` with
a stated reason, and the run's `Degradations` section says so explicitly. See
`bt5.structure.vienna.degradation_reason()`. There is no code path in which a
missing folding engine is silently treated as "no folding liability found".

## commec database

`commec` screens the input **protein**, before optimization, for biosecurity
concerns (`core.context.BiosecurityVerdict`). It ships in two tiers
(`docs/PLAN.md:664`):

- **Biorisk-only** — HMM search against a curated hazard list. Under 1 GB,
  runs on a laptop, and is the tier BT5 ships **on by default** once the
  `screen` extra is installed.
- **Full protein/nucleotide similarity** — comparison against a much larger
  reference corpus. Roughly 600 GB. Explicit opt-in only; never fetched or run
  implicitly by anything in this repository.

**Today's actual state, stated plainly so this section cannot be read as more
than it is:** protein-level screening is not wired into `design()` yet — every
run reports `BiosecurityVerdict("not_run", ...)`, and `bt5 design`'s output
never prints "clear" for that status (`core/context.py`,
`design/runner.py:152`). Wiring the screen into the design path is M2/S2's
lane (`bt5/cassette/`), not this one. What follows is the install-time half of
that flow — getting the database present and ready for the day the design
path calls into it — not a claim that screening runs today.

Installing `commec` (the `screen` extra) does not, by itself, fetch either
database tier — `commec` is a scanner, and a scanner without a database scans
nothing. Fetch the biorisk-only database following
[commec-databases](https://github.com/ibbis-bio/commec-databases/)'s own
instructions before relying on any screening result; do not attempt the 600 GB
full tier unless you specifically need cross-organism similarity search and
have the disk and bandwidth for it. `docs/PLAN.md` does not prescribe a BT5-side
flag or environment variable for choosing between the two tiers — because
nothing in this buildout invokes `commec` yet, there is nothing here for such a
flag to control. Recording that gap rather than inventing a flag with no
consumer: see `docs/decisions/` for this session's decision log.

## Out of scope here

- **Tauri sidecar / desktop packaging.** The process boundary to a future
  desktop wrapper is a hard OpenAPI contract over M9's server
  (`docs/PLAN.md:177-178`), which does not exist yet either. A `uv`-installed
  core is what that wrapper would sidecar; this document is that core's story.
- **macOS notarization.** Downstream of the Tauri packaging step above.
