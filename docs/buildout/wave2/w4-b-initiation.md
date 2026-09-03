# W4 — B6 and B7: translation initiation, and the folklore next to it

*Copy this whole file into a fresh Claude Code session on the BT5 repo.*

**Launch:** permission mode `acceptEdits` · model opus · effort high (both the repo default,
so no override needed) · runs unattended.
**Do not run this in plan mode** — three other sessions are running in parallel and the gate
is the draft PR.

---

You are building the two bacterial initiation rules the brief rates highest and the catalog
does not have. Four of eleven 2.B rules exist (B1, B2, B8, B9); you are adding two more.

## Read this first

Run `/bootstrap` — a fresh checkout has no `.venv`, and `gates.sh` exit **10** means BROKEN,
not a code failure. Then `CLAUDE.md`, then `.claude/rules/rules-catalog.md` (which governs
your lane specifically), then `docs/buildout/wave2/README.md`.

Your branch: **`claude/w4-b-initiation`**. Cut it yourself, from a **freshly fetched** main:

```bash
git remote prune origin
git fetch -q origin main && git checkout -B claude/w4-b-initiation origin/main
```

`git fetch --prune origin main` does **not** clear a stale ref — pruning is bounded by the
refspec, so that form prunes only `origin/main` (CLAUDE.md §7a).

## What to build

Pull both rows verbatim with `docs-miner` before you type a threshold — **never read
`brief.md` (63 KB) inline, and never paraphrase a number.**

### B6 — RBS/TIR model (Salis) · `brief.md:66` · evidence A

`r ∝ exp(−β·ΔG_tot)`, median ~2.3-fold error. This is one of the best-cited bacterial models
in the whole brief.

### B7 — internal initiation · `brief.md:67` · evidence B

Score every internal `ATG`/`GTG`/`TTG` with the TIR model; flag any internal TIR greater
than **10%** of the intended start's.

B7 consumes B6, so build B6 first and B7 against it.

## The trap, and it is the whole reason this session is separate

**B7 must NOT penalize internal Shine-Dalgarno as a "pausing" signal.**

The brief says so in B7's own row, and `brief.md:28` lists internal-SD pausing among the
**overturned** claims — it was an artifact. `docs/PLAN.md:113` echoes it: "Several standard
rules are folklore. Internal Shine-Dalgarno 'pausing' was an artifact…"

B7 is a **TIR-based flag**, never a motif ban. A session that ships a "remove internal SD
sites" rule has reintroduced folklore this project explicitly rejects, and it would sail
through CI because nothing mechanical distinguishes a well-cited motif rule from a
cargo-culted one. That distinction is `rule-auditor`'s job and yours.

Say in the rule's docstring, in one sentence, what B7 is *not*. The next session to read it
should not have to rediscover this.

## Degradation: no folding engine, no number

B6 needs ViennaRNA (the `fold` extra) for `ΔG_tot`. Where the engine is absent, the
objective reports **`unavailable` with a reason** — it does not guess, and it does not
return 0.0.

`b1_five_prime` and `b2_structure_windows` are the precedents; read how they word it, and
note `bt5.structure.vienna.degradation_reason()` already exists for exactly this. A rule
whose thresholds are calibrated against a specific folding engine also carries
`engine_calibration`, and `core/registry.check_engine_calibration` **raises** rather than
silently skipping when the active engine differs — because a skipped rule is a missing
constraint nobody sees. Decide whether B6 needs that field and say why in the docstring.

**ΔG never goes in a byte-exact snapshot** (`CLAUDE.md` §6, `.claude/rules/tests.md`).
ViennaRNA is pinned deliberately and a version bump is a *scientific* change carrying
`approved:algorithm-change`. Assert on bands and orderings, not on digits.

## If the full model is too large for one PR

Ship **B6's ΔG decomposition** and **B7's internal-start enumeration**, with the ranking
term reported `unavailable`, and say so plainly in the PR and the docstring.

A half-model that reports a number is worse than one that reports its absence — the same
argument `c1_cai._unavailable` makes for returning NaN rather than 0.0, and the same shape
`d3_splicing` already ships for its absent MaxEntScan model. Then file the issue for the
remainder.

What you must **not** do is invent a simplified TIR formula and cite Salis for it.

## Enforcement class — decide it deliberately

`CLAUDE.md` §3.5: hard constraints are **never** enforced by a penalty weight, and
`default_weight` must be 0.0 for `HARD_LATTICE`, `HARD_REPAIR` and `HARD_CHECK`.

The brief marks B6 as `S` and B7 as `H/S`. A scored TIR objective is SOFT; an internal start
whose TIR exceeds the intended start's by a wide margin may be repair-shaped. Both are
**bacterial** — use `gate` / `enforcement_for` so they do not fire in a mammalian slot, the
way `b8_kozak` gates per host and `d6_non_b_dna` escalates on bacterial.

Note that `enforcement_for(slot)` and the class-level `enforcement` ClassVar are read at
different points, and disagreeing carelessly is how two shipped presets came to weight a
rule that was hard in the modality they were pinned to (#72, fixed in #76). If your rule
escalates per slot, its ClassVar is the **floor**.

## Wave contracts you must honour

- **Declare `conflicts_with`.** W2 is lighting up the conflict panel this wave; a rule that
  opposes an existing one and does not declare it is invisible there. B7 plausibly opposes
  `b8_kozak` (which already declares `b9_out_of_frame_atg`), and B6's ΔG term plausibly
  opposes `b1_five_prime`. Check, and declare what you find.
- **Do not reuse a `brief_ref`.** `_index_by_ref` (`score/presets.py:120-127`) raises when
  two rules claim one, and that breaks every preset at import.
- **No preset weights your rule, and that is fine.** The three shipped presets weight only
  `2.B1`, `2.C1`, `2.C3`, `2.D4`, `2.F2`. Under W2's policy a new SOFT rule reports its
  percentile at weight 0.0, flagged as outside the preset's claim. If B6 *should* count for
  `ecoli_expression` — and note that `BACTERIAL`'s rationale already argues B1 carries the
  highest weight in any BT5 preset precisely because its published effect size justifies it
  — then **open an issue** naming the `brief_ref`, the proposed weight and the note.
  `score/presets.py` is W2's file, not yours.
- **`packages/engine/tests/rules/conftest.py` is READ-ONLY.** W3 and W5 are editing tests in
  the same directory. A helper only you need goes in your own test file.
- **Do not add a data file.** `data/` is a mutex nobody holds this wave.

## G7 — the shared budget

`packages/engine/tests/design/test_timing.py:47` sets `G7_SECONDS = 10.0` at 500 aa, and
`main` was at **7.37 s** before this wave. Tier-B repair evaluates up to `max_candidates`
candidates per iteration with a **full catalog pass each**, so every rule you add multiplies
into the dominant cost term — and B6 is a *folding* rule, which is the expensive kind. W1 is
re-arming a sweep axis worth ~2.4 s in parallel.

`FoldEngine.mfe` is reserved for report time by its own docstring (~0.24 s at 1 kb, ~6.5 s
at 3 kb) and **must never run inside the empirical null**; `windowed_fold_only=True` is
hard-coded in `design/ranking.py` for that reason, and a test proves it by running every
scored objective against a fold engine whose `mfe` raises. **Your rule must survive that
test** — if B6 calls `mfe`, it will fail, and the fix is to use windows, not to weaken the
test.

Measure G7 before and after on your branch and **put the delta in your PR.**

**Never raise `G7_SECONDS`. Never mark the test `slow`** — `-m "not slow"` currently
deselects nothing, so the marker's only effect would be that `gates.sh` and CI never run the
gate. If the budget is exceeded, that is a **finding**: `docs/PLAN.md:508` gives G7's fail
consequence verbatim as *"re-allocate budget before rules multiply cost."*

## How to add a rule here

**Use `/rule-add`.** It scaffolds the rule with its paired test, its provenance and a
resolvable `brief_ref`, and it knows the assertions in
`tests/data_integrity/test_rule_contract.py` your rule must satisfy.

Every rule ships with citations, an evidence badge, `last_verified`, and — if SOFT — a
`weight_provenance` explaining its default. `weight_provenance` for a SOFT rule is not a
formality: it is the sentence `rule-auditor` will hold to the citation.

## Files

- `packages/engine/src/bt5/rules/catalog/b6_*.py`, `b7_*.py` (new)
- `packages/engine/tests/rules/test_b6_*.py`, `test_b7_*.py` (new)
- `docs/decisions/2026-09-03-*.md`

Nothing else. Not `presets.py`, not `conftest.py`, not `data/`, not `structure/`.

## Delegation

- **Every threshold and every effect size** → `docs-miner`, verbatim, with `brief.md:LINE`.
- **How `b1_five_prime` words its unavailability and its windowing** → `Explore`.
- **Gates** → `gate-runner`.
- **`/pre-pr` will fire `rule-auditor`** (opus/xhigh) by construction. Its question is "does
  the cited source actually support the number?" — and for B7 it will also ask whether you
  reintroduced the internal-SD folklore. Have both answers in the docstring.

## Done means

- `bash scripts/gates.sh` → `ALL GATES PASSED`, exit 0.
- B6 and B7 each ship with their paired test, citations, and an honest `unavailable` path
  when no folding engine is present.
- B7's docstring says in one sentence what it is **not** — a motif ban on internal
  Shine-Dalgarno — and cites `brief.md:28`.
- No ΔG appears in a byte-exact snapshot.
- `conflicts_with` declared and justified.
- Your G7 delta is in the PR.
- A decision record naming the enforcement class for each rule and why, plus what you
  rejected — including, if it applies, the decision to ship a partial model with its
  ranking term `unavailable`.
- Scientific impact: **non-"none"** — this changes emitted sequences. **Owner merge; you do
  not self-merge.**
