# Buildout sessions — closing BT5 to PLAN's v1 bar

Six prompts, written to be run **at the same time** in six separate Claude Code
sessions. Each file in this directory is a complete prompt: open it, copy the whole
thing, paste it into a fresh session on this repo.

This index holds what the six have to agree on. Read it once before launching any of
them; the session files assume you have.

## Why these six

`main` is green. The engine is ~17k lines with zero stubs and one deliberate
`NotImplementedError` (`codon/tables.py:201`, tAI/stAI/CSC, deferred). What is missing
is not machinery — it is **wiring and surface**.

PR #71 landed the walking skeleton (`bt5/design/`, lane M11) and its own docstring
lists what it refuses to do: one candidate and no gallery, every objective shipped
`unavailable`, `native_baseline` None, `BiosecurityVerdict("not_run")`, no order CSV,
no percentiles. Meanwhile `bt5/score/` already exports `build_gallery`,
`null_distribution`, `percentile_of`, `order_entries` and `write_idt_plate`. Nothing
calls them from the design path.

`docs/PLAN.md:490-495` defines the target:

> protein → validated → **screened** → CDS planned → spliced into the circular
> backbone → mutation space over the CDS only → Tier-A DP → Tier-B repair →
> independent verify_construct → **normalized scorecard → 5-candidate gallery →
> annotated GenBank + order CSV**, under 10 s

Six sessions close that, plus a CLI so a user does not have to write Python.

## Ownership matrix

Every writable path belongs to exactly one session. If your work needs a path you do
not own, **stop and open an issue** — do not reach into another lane.

| # | Session | Branch | Owns (write) | Mutex |
|---|---|---|---|---|
| S1 | [Ranking increment](s1-ranking-increment.md) | `claude/s1-ranking-increment` | `bt5/design/**`, `bt5/score/**`, `packages/engine/tests/{design,score}/**` | — |
| S2 | [Biosecurity screen](s2-biosecurity-screen.md) | `claude/s2-biosecurity-screen` | `bt5/cassette/**`, `packages/engine/tests/cassette/**` | — |
| S3 | [Rules: translation](s3-rules-translation.md) | `claude/s3-rules-translation` | `rules/catalog/{b,c}*.py` + paired tests | — |
| S4 | [Rules: liabilities](s4-rules-liabilities.md) | `claude/s4-rules-liabilities` | `rules/catalog/{d,e,f}*.py`, `rules/vendors.py`, `rules/_provenance.json` + paired tests | — |
| S5 | [CLI + packaging](s5-cli-packaging.md) | `claude/s5-cli-packaging` | `bt5/cli.py`, `packaging/**`, `packages/engine/tests/cli/**` | **`pyproject.toml`** |
| S6 | [Host data + backbone](s6-host-data.md) | `claude/s6-host-data` | `data/**`, `tests/data/backbones/**` | **`data/`** |

### The three global mutexes

1. **`pyproject.toml` — S5 only.** Every dependency is already declared (`server`,
   `screen`, `export` extras included), so nobody needs a new one. But
   `[project.scripts]`, `[tool.mypy] files` and the hatch `packages` list live there.
   `CLAUDE.md` §5 calls a lockfile conflict across parallel PRs the single most
   expensive merge failure in this repo.
2. **`data/` — S6 only.** Carries `approved:data-change` and goes to the owner.
3. **`core/` — nobody, by default.** `tests/contract/regenerate.py` rewrites
   `manifest.json` and 17 fixtures; two sessions regenerating in parallel conflict
   irreconcilably. Claim it by opening an issue naming the type, wait for it, then use
   `/contract-change` — **classification before regeneration**, because
   `regenerate.py` writes first and returns 0 on every path.

## Launch table

Sessions launched in this repo inherit `model: opus`, `effortLevel: high` from
`.claude/settings.json`. Only the **bold** cells are overrides.

| # | Permission mode | Model | Effort | Unattended? |
|---|---|---|---|---|
| S1 | `acceptEdits` | opus | high | yes |
| S2 | `acceptEdits` | opus | **xhigh** | yes |
| S3 | `acceptEdits` | opus | high | yes |
| S4 | `acceptEdits` | opus | **xhigh** | yes |
| S5 | `default` | **sonnet** | high | **no** |
| S6 | `default` | opus | high | **no** |

**Not plan mode.** Plan mode ends by blocking on approval, so six plan-mode sessions
means six gates one person must clear before any code exists — it turns the
parallelism into a queue. It is also redundant: this repo already gates on
`/pre-pr` → draft PR → owner merge, and `CLAUDE.md` §7b governs what an agent may
merge alone.

**Not `bypassPermissions` either.** `.claude/hooks/protect_paths.py` is what stops a
session editing `core/`, `verify.py`, `.github/`, `data/` or `pyproject.toml` without
stopping to think. Its decision is deliberately "ask, never deny". Bypassing
permissions disarms the mechanism keeping six sessions off each other's protected
paths.

S5 and S6 run in `default` mode and need a person within reach, because their work is
*defined* as editing protected paths — they are the two mutex holders. The other four
run unattended.

## The un-draft queue

CI capacity is the binding constraint: 20 concurrent job slots, ~12 per Python PR, so
**at most 5 open non-draft PRs**. Drafts skip the expensive jobs and are free.

- Open your PR as a **draft** and keep it there until you believe it is done.
- Before flipping to ready, list open non-draft PRs. **At 4 or more, stay in draft**
  and re-check after something merges. Four, not five, leaves a slot for a re-run.
- Two ready at once: the **lower session number goes first**.

## Inter-session contracts

What one session may rely on another not breaking.

- **`design()`'s signature is frozen** (`design/runner.py:156-171` — keyword-only,
  `table_id` never defaulted). S1 may *add* fields to `SkeletonResult`; it may not
  remove or rename one. S5's CLI is built against it.
- **`BiosecurityVerdict`'s shape is frozen** (`core/context.py:96-115`). S2 changes
  behaviour, not the type. S1 renders whatever status it gets and never prints
  "clear" for `not_run`.
- **Rule registration stays autodiscovery.** `core/registry.py` walks
  `bt5.rules.catalog` with `pkgutil`; `rules/__init__.py` and
  `rules/catalog/__init__.py` are both empty, so adding a rule edits **zero** shared
  files. Neither S3 nor S4 may introduce a hand-maintained list — that would create
  the one collision this design avoids.
- **S6 ships data only.** Wiring a new reference set into `CAI_REFERENCE_SET` is S3's
  edit to `c1_cai.py`, taken up only if S6 merged first; otherwise a follow-up issue.
  A rule with no reference set reports `unavailable`, per the pattern
  `docs/decisions.md` records for C1.

## Shared rules

`CLAUDE.md` loads automatically and is not repeated here. What it does **not** say:

- **Bootstrap first.** A fresh checkout has no `.venv`; run `/bootstrap`. Every
  command uses `.venv/bin/…`. `gates.sh` exit **10** means no venv — that is BROKEN,
  not a code failure. Bare `pytest` exits 4 on a `conftest.py` import error and looks
  like a real failure.
- **`main` is green at `628e130`.** Issues #80 and #83 claim otherwise and are **false
  alarms**: `ci.yml`'s `main-broken` job fires on `cancelled` as well as `failure`,
  and `concurrency.cancel-in-progress` cancels the prior run on rapid merges. Every
  job in both issues reads `cancelled`, none `failure`. Do not chase them and do not
  "fix" main. (Fixing the job itself needs `approved:ci-change` and is out of scope.)
- **`docs/decisions.md` is newest-first**, so every session inserts at the same place
  and conflicts are expected. Resolve by keeping **both** entries in date order. Never
  drop another session's entry.
- **Spend context on judgment, not retrieval.** These are long sessions; what ends one
  early is a window full of file dumps and gate output. Route retrieval to `Explore`
  and `docs-miner`, gates to `gate-runner`, and keep the main thread for decisions
  only it can make.
- **Escalate on the right axis.** *Capability failure raises the model; diligence
  failure raises the effort.* Re-reading a file at higher effort does not fix a
  capability failure, and a bigger model does not fix a skipped checklist. If a
  first-pass fix already failed, that is `debugger`, not another attempt.

## Routing

The fleet every session inherits:

| Route | Model / effort | Use it for |
|---|---|---|
| `Explore` | haiku / low | Locating call sites and symbols. Never for judging correctness. |
| `docs-miner` | sonnet / medium | Any `docs/` claim. Verbatim quotes with file:line; flags superseded rows. |
| `gate-runner` | haiku / high | `scripts/gates.sh` and CI diagnosis. |
| `code-reviewer` | sonnet / high | The complete branch diff, via `/pre-pr`. Has project memory. |
| `batch-editor` via `/cheap-pass` | sonnet / medium | ≥5 identical edits, before/after already decided. |
| `rule-auditor` | opus / xhigh | Does the cited source actually support the number? |
| `security-reviewer` | opus / xhigh | Biosecurity posture and the CI trust boundary. |
| `debugger` | opus / xhigh | A failure a first-pass fix already missed. |
| `/architect`, `/escalate` | opus / xhigh + ultracode | Cross-lane, `core/`, MINOR-vs-MAJOR, or a decision already got wrong. |

**`/pre-pr` fires reviewers conditionally, by path.** `security-reviewer` runs iff the
diff touches `vector/`, `core/services.py`, `verify.py`, `cassette/` or `.github/`.
`rule-auditor` runs iff a Spec's `citations`, `weight_provenance`, `enforcement`,
`last_verified` or a threshold changed. S2 trips the first by construction; S3 and S4
trip the second on every PR. Budget for that opus/xhigh pass rather than being
surprised by it.

**Never bare-Read a source file over 20 KB.** Each session file lists the ones its own
lane will actually meet.

## Out of scope

Stated so no session wanders into it:

- **M9 `packages/server/` and M10 `apps/web/`.** Note `ci.yml:53-54` already carries an
  `apps/web/**` paths filter whose output no job consumes, so a web lane needs a new CI
  job under `approved:ci-change`.
- **Issue #45's ε-constraint rearchitecture** (X2, X3, X4, X5, X7).
- **The oracle backlog:** #69, #70 (`verify.py`, `approved:oracle-change`), #82, #52, #64.
- **`benchmarks/`** — the directory does not exist at all, despite `CLAUDE.md` §2
  protecting `baseline.json` and `tolerances.yaml`. Creating it is an owner decision
  under `approved:algorithm-change`. S1 carries a plain timing assertion instead.
