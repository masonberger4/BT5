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
