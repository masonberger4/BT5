## 2026-09-01 — the remaining buildout runs as six parallel sessions, split by write-ownership

**Decided:** `docs/buildout/` holds six copy-pasteable session prompts (S1 ranking
increment, S2 biosecurity screen, S3 rules 2.B/2.C, S4 rules 2.D/2.E/2.F, S5 CLI +
packaging, S6 host data) plus a README carrying the ownership matrix, the un-draft
queue, the inter-session contracts and the routing table. Scope is the engine to
PLAN's v1 bar (`docs/PLAN.md:490-495`) plus a CLI.

The split axis is **write-ownership, not topic**: every writable path belongs to
exactly one session, and each prompt names both what it owns and what it must never
touch. Three global mutexes are assigned by name — `pyproject.toml` to S5, `data/` to
S6, and `core/` to nobody (claimed by opening an issue first, because
`tests/contract/regenerate.py` rewrites `manifest.json` and 17 fixtures and two
parallel regenerations cannot be merged).

Each prompt also pre-applies `CLAUDE.md` §Delegation to its own lane's sub-tasks
rather than pointing at the paragraph, and states its permission mode, model and
effort. S2 and S4 override effort to xhigh; S5 overrides model down to sonnet.

**Rejected:**
- *Plan mode for the sessions.* Plan mode ends by blocking on approval, so six
  simultaneous plan-mode sessions become six gates one person must clear before any
  code exists — it converts the parallelism into a queue. The repo already gates on
  `/pre-pr` → draft PR → owner merge. S2 and S4 post a design note as their first PR
  comment instead, which gives the same chance to redirect without blocking.
- *`bypassPermissions`.* `.claude/hooks/protect_paths.py` is what keeps six sessions
  off each other's protected paths; its decision is deliberately "ask, never deny".
  Bypassing disarms exactly the mechanism the parallel design depends on.
- *Splitting by topic (one session per rule family, one per module).* Topic boundaries
  cut across files — a "vendor rules" session and a "repeat rules" session would both
  write `rules/vendors.py`. Ownership boundaries do not.
- *A wider four-session split.* Fewer sessions means larger PRs against a CI gate that
  admits at most 5 non-draft PRs anyway; the constraint is slots, not sessions, and
  drafts are free.
- *Two sequential waves.* Roughly doubles wall-clock to remove conflicts that the
  ownership matrix already removes.
- *Including M9 `packages/server/` and M10 `apps/web/`.* Both need `pyproject.toml`
  (S5's mutex) and a web lane additionally needs a new CI job under
  `approved:ci-change` — `ci.yml:53-54` already carries an `apps/web/**` paths filter
  whose output no job consumes.
- *Creating `benchmarks/` for PLAN's G7 timing bar.* The directory does not exist at
  all despite `CLAUDE.md` §2 protecting `baseline.json` and `tolerances.yaml`;
  creating it is an owner decision under `approved:algorithm-change`. S1 carries a
  plain `pytest` timing assertion instead.

**Evidence:**
- The gap is wiring, not machinery: `bt5/score/__init__.py:16-50` exports
  `build_gallery`, `null_distribution`, `percentile_of`, `order_entries` and
  `write_idt_plate`; nothing under `bt5/design/` calls any of them, and
  `design/runner.py:358` appends the note "ranking not computed: no null distribution
  and no percentiles".
- Rules are collision-free in parallel: `core/registry.py` autodiscovers via
  `pkgutil.walk_packages` and `@register`; `rules/__init__.py` and
  `rules/catalog/__init__.py` are both 0 bytes, so adding a rule edits no shared file.
- CI capacity is the binding constraint: `CLAUDE.md` §7 gives 20 concurrent job slots
  at ~12 per Python PR, so at most 5 open non-draft PRs. The queue targets 4, leaving
  a slot for a re-run.
- Issues #80 and #83 ("main is red") are **false alarms**, and the prompts say so:
  `ci.yml`'s `main-broken` job fires on `contains(needs.*.result, 'cancelled')` as
  well as `'failure'`, while `concurrency.cancel-in-progress: true` cancels the prior
  run on rapid merges to `main`. Every job in both issues reads `cancelled`, none
  `failure`; `main` is green at `628e130` (workflow run conclusion `success`).
- `data/genetic_codes/` does not exist on disk despite `CLAUDE.md` §2 protecting it,
  and `data/codon_usage/` holds one file. S6 is told to report this, not to invent
  files that satisfy the path.

**Where:** branch `claude/project-buildout-prompts-s2ocnr`.
