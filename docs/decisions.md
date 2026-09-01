# Decisions

Settled decisions from working sessions, appended at the end of each work slice.
What survives compaction is what lives on disk — a decision that exists only in a
conversation is gone at the next compaction.

**Scope, against `docs/rfcs/`.** RFCs record *contract amendments*: they are
load-bearing in CI, `check_amendment.py` reads the manifest they correspond to, and a
MAJOR change is unmergeable without one. This file records *session decisions that no
gate enforces* — what was tried, what was rejected, and why. If a decision changes
`bt5/core/`, it belongs in an RFC, and this file just points at it.

**Format.** Newest first. One entry per decision, not per session.

```
## YYYY-MM-DD — one-line summary
**Decided:** what will happen.
**Rejected:** the alternatives, each with the reason it lost.
**Evidence:** file:line, a command's output, or a measurement.
**Where:** PR / branch / commit, if there is one.
```

---

## 2026-08-31 — Model and effort routing, and eight findings it surfaced

**Decided:** Route work by failure mode — capability failure raises the model, diligence
failure raises the effort. Eight subagents (`Explore`, `docs-miner`, `gate-runner`,
`batch-editor`, `code-reviewer`, `rule-auditor`, `debugger`, `security-reviewer`), eight
skills, four path-scoped rules files, six hooks. Session default `opus` / `high`.
`/escalate` and `/architect` run opus at ultracode and orchestrate rather than answer
alone. Full rationale in `SETUP-NOTES.md`.

**Rejected:**
- *A `ci-triage` agent.* Its input — the diff — is already in the main window, and a job
  log is 2–10k tokens against a subagent's fixed overhead. Folded into `gate-runner`, so
  local and remote failure share one exit-code vocabulary.
- *A line-number map for `brief.md`.* Verified there are two anchor shapes, not three, so
  a resolution procedure covers every case — and line numbers rot on the next `docs:`
  commit with no gate to catch it.
- *`max` effort anywhere.* Ultracode is `xhigh` + orchestration by the platform's own
  definition, and Opus 5's cost index puts `max` only ~6% above `xhigh`. The capability
  comes from fanning out, not from more thinking in one context.
- *Pinning `autoCompactWindow`.* It takes a token count, and the CLI's own text says
  `auto` is "tuned for your model and strongly recommended"; overriding "may result in
  high token usage". Left at `auto`; the statusline surfaces context pressure instead.
- *An `architect` subagent.* `CLAUDE_CODE_MAX_SUBAGENT_SPAWN_DEPTH=1`, so a subagent
  cannot orchestrate. Moved to `/architect`, which runs in the main loop.

**Evidence:** see the findings list in `SETUP-NOTES.md`; each carries a `file:line`.

**Superseded in part:** #63 landed on `main` while this branch was open and fixed three of
the ten findings — `HYPOTHESIS_PROFILE` is now honored, the RNG grep is broadened to any
`np.random.*` plus stdlib `random` across source and tests, and `mypy --strict` is a
required CI job. This branch is merged with it and every affected statement in
`CLAUDE.md`, `.claude/rules/tests.md`, `.claude/rules/vector.md`, the `gate-runner` and
`security-reviewer` agents, `/pre-pr`, `scripts/gates.sh` and the SessionStart hook was
corrected. The lesson worth keeping: **config that asserts repo facts goes stale like
code**, and nothing gates it — `.claude/verify-setup.sh` checks structure, not truth.

**Where:** branch `claude/model-effort-routing-ovsb5r`.

---

## 2026-09-01 — The §3.5 weighting guard asks per slot, and 2.D4 leaves two presets (#72)

**Decided:** `score/presets.py` `resolve()` guards on `enforcement_for(slot)` for every
slot the preset's modality admits, not on the `enforcement` ClassVar. The ClassVar is a
FLOOR — `d4_internal_polya` declares `SOFT` (right for the plasmid case) while its
`enforcement_for` returns `HARD_REPAIR` on every packaged modality — so the old guard
asked whether a rule is hard *everywhere* instead of whether it is hard *here*, and
passed. `LENTIVIRAL` and `AAV` both shipped `WeightEntry("2.D4", 1.0)` through it. Those
entries are removed; that is the actual fix, and the resolver change is what stops it
being reintroduced.

**No signature change.** `resolve()` still takes `specs`, not a `DesignContext`.
`Preset.modality` is the pin, and every `enforcement_for` in the catalog keys on
`slot.modality` alone. Requiring a `DesignContext` would also be circular:
`ResolvedPreset.weights` is an *input* to `DesignContext.weights`, so the context does
not exist yet at the moment `resolve()` runs.

**Rejected:**
- *Naming one representative host per modality to build the probe slot.* `ContextSlot`
  cannot be built from a modality alone — `role` and `host` have no defaults and
  `__post_init__` locks `table_id` to the host — so a single probe means guessing a
  host. It answers correctly today (no catalog rule keys `gate` or `enforcement_for` on
  `host` or `role`) and goes silently wrong the first time one does, which is the same
  shape as the bug being fixed: asking a narrower question than the one that decides the
  answer. `_slots_admitted_by` enumerates `LOCKED_TRANSLATION_TABLE` × `SlotRole`
  instead, so nothing is guessed and nothing is defaulted.
- *Guarding on `is_hard`.* `is_scored` is `is SOFT`, so `is_hard` would have started
  admitting `REPORT_ONLY` into the weighted sum — a weakening smuggled in as a fix.
  The guard refuses anything not scored, as before.
- *Keeping the ClassVar guard and documenting the per-slot check as the consumer's job*
  (option 2 in #72). `bt5/design/catalog.py` already does that, but a guard that has to
  be re-implemented by every consumer is not a guarantee. Left untouched — it is PR #71's
  file.
- *Deleting `_POLYA_NOTE`'s science with its entry.* The 8–9× functional titer loss is
  why d4 is HARD in these modalities; it is now a comment above the presets explaining
  why 2.D4 is deliberately absent, plus the reason the old guard let it through.

**Also corrected:** `LENTIVIRAL.rationale` claimed weight went to "internal polyA on the
packaged strand, cryptic splice donors". Neither was weighted after this change, and the
splice claim was never true in any form — there is no splice-donor rule in the catalog at
all (15 rules, none for splicing), so the prose asserted an objective that does not
exist. It now names internal polyA only, and says it carries no weight *because* it is
hard.

**Scientific impact — not "none".** d4 leaves the weighted sum for lentiviral and AAV
designs, so candidate ranking moves for both. That is the correct direction: the
constraint is enforced by Tier-B repair plus the independent validator, which refuses to
emit, and weighting it as well both double-counts it and tells the user a guarantee is a
trade-off. Owner merges under §7b.

**Where:** branch `claude/hard-rule-weighting-guard-nitkbe`; lane M3, `score/` only.
