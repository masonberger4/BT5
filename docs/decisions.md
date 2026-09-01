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

## 2026-09-01 — `expand_forbidden` supports IUPAC, with a solver-local expansion table

**Decided:** `LatticeTerms.forbidden` is documented IUPAC (`core/spec.py:196`) but
`expand_forbidden` assumed ACGT, so a degenerate base died as a bare `KeyError` on
`lattice._BASE_INDEX`. Option 1 from #73 — make the engine match the contract — over
option 2 (narrow the docstring to ACGT), because option 2 edits `core/spec.py` and needs
`approved:contract-change` plus an owner merge for what is an M1 bug.

Three parts:

1. **Expansion precedes the reverse-complement closure**, and the order is load-bearing.
   `core.types.reverse_complement` is `str.translate`, which leaves an unmapped character
   ALONE, so revcomp-first does not fail on a degenerate pattern — it half-complements.
   `RGATC` comes back `GATCR`, expanding to {GATCA, GATCG} when the true complements are
   {GATCT, GATCC}: the wrong sequences forbidden, the right ones permitted, silently.
   Pinned by `test_expansion_precedes_the_reverse_complement_closure`.
2. **`MAX_PATTERN_EXPANSION = 1024`, per pattern**, computed from the code widths before
   any string is built. Five fully-ambiguous positions or ten two-fold ones; every
   consensus motif a rule would plausibly declare (`MAGGTRAGT` 4, `GCCRCCATGG` 2, `GGWCC`
   2) sits orders of magnitude below. It excludes the N-spacer interrupted palindromes
   (BstXI `CCANNNNNNTGG` 4,096; XcmI `CCANNNNNNNNNTGG` 262,144) — those are wildcards,
   not ambiguity, and enumerating them is the wrong mechanism. The cost is paid twice
   downstream: an Aho-Corasick state per pattern prefix, then a 64-codon transition row
   per state in `lattice._codon_transitions`, so an all-N 8-mer reads as a hung solver.
3. **A `ValueError` naming the pattern** on a non-IUPAC character, an empty pattern, and
   an over-cap expansion. Eliminating the bare `KeyError` was the point of #73 either way.

**Rejected:**
- *Importing `verify.IUPAC_EXPANSION` / `verify.expand_iupac` into the solver.*
  `tests/data_integrity/test_oracle_independence.py` exists to keep the oracle off every
  lane's code path; sharing the expander defeats it pointing the other way. One
  transposed row would then be invisible, because the design and its check would forbid
  the same wrong set — the differential stops differentiating. The solver keeps its own
  `IUPAC_CODES`, and `test_the_solver_and_the_oracle_agree_on_expansion` asserts the two
  agree on a shared vector of eight patterns, so divergence fails loudly instead of
  silently. Two tables plus one test beats one table and no check.
- *A cap on the TOTAL expanded set rather than per pattern.* The blowup is exponential
  per pattern and merely linear across them; a rule listing fifty pure-ACGT six-cutters
  is fine today and must stay fine. A per-pattern cap is also the only one whose error
  message can name the offender.
- *Option 2 (refuse IUPAC, narrow the docstring).* Protected path, owner merge, and it
  removes a capability the contract advertises rather than supplying it.

**Scientific impact: none.** The only two rules populating `forbidden` today —
`d1_restriction_sites` (six-cutters) and `e1_homopolymers` (A/G runs) — emit pure ACGT,
which expands to itself. `expand_forbidden(["GGTCTC"]) == ("GAGACC", "GGTCTC")` still
holds. No shipped sequence changes.

**Follow-up for the design lane:** #71's `build_catalog` guard asserting lattice terms
are ACGT-only and raising `DesignError` is obsolete once this lands, and can be dropped.

**Evidence:** `packages/engine/src/bt5/solver/reference.py:35` (`expand_iupac`,
`expand_forbidden`), consumers unchanged at `reference.py:85`, `lattice.py:145`,
`repair.py:449`; `bash scripts/gates.sh` all eight gates exit 0;
`HYPOTHESIS_PROFILE=ci` on `tests/invariants` and the two new properties.

**Where:** branch `claude/iupac-forbidden-expansion-orzxz4`, closes #73.

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
