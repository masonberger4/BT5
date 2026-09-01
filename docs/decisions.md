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

**Evidence:** `packages/engine/src/bt5/solver/reference.py:44` (`IUPAC_CODES`), `:80`
(`MAX_PATTERN_EXPANSION`), `:83` (`expand_iupac`), `:126` (`expand_forbidden`); the three
consumers call it unchanged at `reference.py:185`, `lattice.py:145`, `repair.py:449`, and
are unaffected because expansion preserves pattern LENGTH, which is all
`_creates_forbidden`'s window bound and `repair.py`'s `guard_len` read. `bash
scripts/gates.sh` all eight gates exit 0; `HYPOTHESIS_PROFILE=ci` on `tests/invariants`
and on the two new properties.

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

---

## 2026-09-01 — C1 (CAI): a soft band that reports unavailable for seven of nine hosts

**Decided:** `c1_cai` ships as `Enforcement.SOFT` / `Direction.BAND`, band `(0.70, 0.90)`
per brief.md:77, `default_weight` 0.2 (matching what all three presets already assign
2.C1), `steering_weight` 0.0, `lattice_terms()` → `None`. 2.C1 now resolves in every
preset; only 2.C3 is left in `ResolvedPreset.unimplemented`.

The band's **ceiling is the operative half**. Kudla 2009 measured CAI at r = 0.14 (n.s.)
against 5' folding's r = 0.66 on the same 154 variants, and Welch 2009's deliberately
high-CAI control expressed at ~15% of the best variant. `Direction`'s own docstring names
this rule as the reason `BAND` exists.

`Evidence.CONTESTED`, not `EVIDENCE_BACKED`: the brief grades C1 "A (that it's weak)" —
grade A for the claim that CAI is a *weak* predictor — while Boel 2016 measured codon
content as 3–5x *more* influential than structure on a different readout. BT5 cannot
adjudicate that, and the citations carry both signs.

**Unavailability is the feature, not a caveat.** `data/codon_usage/` holds exactly one
reference set (`sharp_li_1987_ecoli_w.json`) against nine `HostId` values, so C1 computes
for the two E. coli hosts and reports the objective unavailable for HUMAN, HEK293, CHO,
S_CEREVISIAE, P_PASTORIS, SF9 and MOUSE — i.e. unavailable under `lentiviral_hek293` and
`aav_hek293`, a number only under `ecoli_expression`. NaN plus a reason-carrying breach,
copying `b1_five_prime._unavailable`.

**Rejected:**
- *A `codon_weights` lattice term, or a non-zero `steering_weight`.* Both are monotone
  pulls toward maximum CAI — the Tier-A DP maximizes what it is handed — which is exactly
  what "never maximized" (brief.md:73) forbids. A band is not expressible in the automaton
  either: it decides from a bounded codon suffix, and CAI is a geometric mean over the ORF.
- *Adding a mammalian codon usage table so the HEK293 path lights up.* `data/codon_usage/**`
  is protected (`approved:data-change`), and a human highly-expressed reference set is an
  evidence-bearing decision with its own provenance burden. Filed as its own issue.
- *Falling back to the E. coli table for a mammalian host.* It is the one table on disk, so
  the fallback would always succeed and always be wrong — a plausible-looking number
  measuring nothing. This is the failure the unavailable path exists to prevent.
- *Gating on `modality` (B1's choice) or on "host is E. coli".* Gating on host is the
  dangerous one: a lentiviral job propagates in E. coli and expresses in HEK293, so a
  host-keyed rule would find the one host with a w-table and report its CAI as the
  objective for a protein made elsewhere. The gate is `slot.role != "propagation"`.
- *Hard-coding ATG/TGG as the excluded single-codon families.* NCBI table 4 makes TGA a
  second Trp codon, so Trp carries information there. Family size is read from the injected
  table: `len(code.synonymous_codons(aa)) > 1`.
- *Importing `CodonUsage.cai` from M5.* `exp(mean ln w)` is recomputed in the rule so the
  lane keeps no import edge into the codon lane; `Services` is what decouples M4 from M5.
- *Editing `packages/engine/src/bt5/score/`.* The predicted ripple did not occur — the
  `unimplemented` assertions in `test_presets.py` use synthetic `fake_spec` lists, not the
  live registry. Nothing under `score/` was touched.

**Note:** `TableProvider.weights()` keys on the REFERENCE SET, not on `HostId` — every
caller in the tree passes `"sharp_li_1987_ecoli_w"`. C1 therefore carries an explicit
`CAI_REFERENCE_SET` map; a host absent from it has no reference set and reports
unavailable rather than raising `FileNotFoundError`. BL21 shares K-12's entry as a stated
same-species approximation, and the reference set travels in `detail["reference_set"]`.

**Where:** branch `claude/rule-c1-cai-enxfre`.
