# Code-reviewer memory — BT5

- [CI workflow review checklist](ci-workflow-review.md) — reporter/re-arm split, `NON_BLOCKING` scoping, and the big one: a gate whose own error message satisfies its own parser
- [Skill frontmatter governance](skill-frontmatter-governance.md) — `/pre-pr`'s `disable-model-invocation` removal was decided 2026-09-03; check `docs/decisions/` before flagging a frontmatter change
- [Text I/O encoding pattern](text-io-encoding.md) — bbffdb8/#102's `encoding="utf-8"` sweep and its AST-scan guard; what to re-check if the guard file is touched again, and the four-protected-path-prefix label reminder
- [General review notes](notes.md) — `score/presets.py` enforcement-per-slot facts (incl. `gate is None` widens, doesn't weaken, the probe), duplicate-guard-across-protected-path pattern, doc/prompt-review checklists, `.claude/*.md` count-drift handling
- [score/null.py PR #112 techniques](score-null-pr112.md) — empirical RNG/chi-square verification, caller-chain "unreachable code" tracing, PR-description-vs-committed-test gaps
