# W3 — D5: bacterial cryptic transcription (and the one rule with no test)

*Copy this whole file into a fresh Claude Code session on the BT5 repo.*

**Launch:** permission mode `acceptEdits` · model opus · effort **xhigh** (override) · runs
unattended.
**Do not run this in plan mode** — three other sessions are running in parallel and the gate
is the draft PR.

---

You are building **D5**, the highest-value unbuilt rule in the catalog and the only one with
a disaster-grade calibration anchor: a dengue-2 cryptic promoter that made a clone
**uncloneable**.

Twenty-five of `brief.md` §2's fifty rules have files. D5 is the one the evidence most
clearly demands, and it carries a sub-rule **no competitor implements**.

## Read this first

Run `/bootstrap` — a fresh checkout has no `.venv`, and `gates.sh` exit **10** means BROKEN,
not a code failure. Then `CLAUDE.md`, then `.claude/rules/rules-catalog.md` (which governs
your lane specifically), then `docs/buildout/wave2/README.md`.

Your branch: **`claude/w3-d5-cryptic-transcription`**. Cut it yourself, from a **freshly
fetched** main:

```bash
git remote prune origin
git fetch -q origin main && git checkout -B claude/w3-d5-cryptic-transcription origin/main
```

`git fetch --prune origin main` does **not** clear a stale ref — pruning is bounded by the
refspec, so that form prunes only `origin/main` (CLAUDE.md §7a).

## What D5 is

`brief.md:110-112`. Pull the rows verbatim with `docs-miner` before you type a single
threshold — **never read `brief.md` (63 KB) inline, and never paraphrase a number.** At
least one row elsewhere in the brief is struck through and superseded (E4's extent
thresholds, corrected 2026-08-28 as below the chance floor), so "the brief says" is not good
enough; you need the current row with its line number.

Five parts, four classical and one novel:

- **(a)** hexamer within 1 mismatch of `TTGACA` + a 15–19 bp spacer + a hexamer within 1
  mismatch of `TATAAT`
- **(b)** extended −10, `TGnTATAAT`
- **(c)** the AT-tract rule — the brief reports *103/103 randomized appropriately positioned
  AT-tracts were active*
- **(d)** σ38 variant `TATACT`
- **antisense promoter (novel):** forbid `ATTATA` — the reverse complement of `TATAAT` — on
  the sense strand. Only **4.76%** of 484,741 natural *E. coli* CDSs contain it, yet
  **77.28%** of clean CDSs can silently acquire it by synonymous substitution.

That last asymmetry is the whole argument for the rule: nature avoids it, and naive codon
optimization walks straight into it.

## The scope call, already made

**Build the motif/spacer machinery and the antisense rule now. Do NOT reimplement the Salis
Promoter Calculator.**

The brief names it — 346 parameters, validated on 17,396 promoters — as the *quantitative*
scorer for D5. It is a model, not a threshold, and reimplementing it inside this PR would be
a second, unvalidated model wearing the citation of a real one.

Instead, **report that scoring path `unavailable` with a reason**, exactly as
`d3_splicing` already does for its absent MaxEntScan model (`d3_splicing.py` — read how it
words the unavailability, and match it). Then **open an issue** naming the Promoter
Calculator, what it would score, and what it would replace.

A half-model that reports a number is worse than one that reports its absence. That is the
same argument `c1_cai._unavailable` makes for returning NaN rather than 0.0.

## The trap: do not scan the reverse strand yourself

`CLAUDE.md` §3.4 is not advice here, it is the thing you will most plausibly get wrong.

The antisense rule is a **forward motif on the sense strand**. `ATTATA` is *already* the
revcomp of `TATAAT`; the rule is "do not put this hexamer on the coding strand", not "scan
the other strand for promoters". So:

- List `ATTATA` (and the classical motifs) in `LatticeTerms.forbidden` and **let the solver
  close the revcomp set.**
- Do not write your own reverse-strand scan. Getting this backwards makes the rule both
  wrong *and* a duplicate of the closure the solver already performs.
- A directional scored model — if you build one for the spacer geometry — **must** read
  `slot.strand_of_interest`, because scored directional models are not revcomp-symmetric.

`docs/decisions/2026-09-01-expand-forbidden-iupac.md` records how IUPAC motifs are expanded
*before* the revcomp closure. If any of your hexamer-with-one-mismatch sets is expressed as
IUPAC, that ordering is load-bearing — read it.

## Enforcement class — decide it deliberately

`CLAUDE.md` §3.5: hard constraints are **never** enforced by a penalty weight, and
`default_weight` must be 0.0 for `HARD_LATTICE`, `HARD_REPAIR` and `HARD_CHECK`.

The five parts do not all want the same class, and saying so in the rule's docstring is part
of the deliverable:

- A bare forbidden hexamer with no context is lattice-shaped (`d1_restriction_sites` and
  `e1_homopolymers` are the precedents).
- A **spacer geometry** — two hexamers at 15–19 bp — is not something the lattice automaton
  can express as a forbidden string, so it is repair- or check-shaped.
- Anything you cannot fix by codon choice (a match falling in backbone-carried sequence) is
  `HARD_CHECK` by definition.

Look at how `d2_recombinase_sites` argued itself into `HARD_CHECK`
(`docs/decisions/2026-09-02-d2-recombinase-check-only.md`) before you choose.

Note also that D5's hazard is **bacterial**. Use `gate` / `enforcement_for` the way
`d6_non_b_dna` does for its bacterial escalation, rather than firing everywhere.

## Also in this PR: the one rule with no test

`packages/engine/src/bt5/rules/catalog/d1_restriction_sites.py` is the **only rule in the
catalog with no paired test file.** It is d-series and it is one file, so it is yours.

Write `packages/engine/tests/rules/test_d1_restriction_sites.py`. It is a `HARD_LATTICE`
rule and one of PR #0's hand-written reference rules, so the test is also documentation of
what a lattice rule's contract looks like. Cover the junction-spanning and origin-spanning
cases — `conftest.py`'s `wrapping_construct()` exists for exactly that.

## How to add a rule here

**Use `/rule-add`.** It scaffolds the rule with its paired test, its provenance and a
resolvable `brief_ref`, and it knows the assertions in
`tests/data_integrity/test_rule_contract.py` that your rule must satisfy. Hand-writing a
rule file means rediscovering those one CI cycle at a time.

Every rule ships with citations, an evidence badge, `last_verified`, and — if SOFT — a
`weight_provenance` explaining its default.

## Wave contracts you must honour

- **Declare `conflicts_with`.** W2 is lighting up the conflict panel this wave. A rule that
  opposes an existing one and does not declare it is invisible there. D5 forbidding
  AT-rich hexamers plausibly opposes `f5_at_window` and `e2_gc_band`; check, and declare
  what you find.
- **Do not reuse a `brief_ref`.** `_index_by_ref` (`score/presets.py:120-127`) raises when
  two rules claim one, and that breaks every preset at import. W2 is adding a
  `data_integrity` assertion that all three presets resolve cleanly, so a violation fails a
  gate rather than a user's run — but it will fail *your* PR.
- **No preset weights your rule, and that is fine.** The three shipped presets weight only
  `2.B1`, `2.C1`, `2.C3`, `2.D4`, `2.F2`. Under W2's policy a new SOFT rule reports its
  percentile at weight 0.0, flagged as outside the preset's claim. If D5 *should* count for
  `ecoli_expression`, **open an issue** naming the `brief_ref`, the proposed weight and the
  note — `score/presets.py` is W2's file, not yours.
- **`packages/engine/tests/rules/conftest.py` is READ-ONLY.** W4 and W5 are editing tests in
  the same directory. A helper only you need goes in your own test file.
- **Do not add a data file.** `data/` is a mutex nobody holds this wave.
- **Report your G7 delta.** See below.

## G7 — the shared budget

`packages/engine/tests/design/test_timing.py:47` sets `G7_SECONDS = 10.0` at 500 aa, and
`main` was at **7.37 s** before this wave. Tier-B repair evaluates up to `max_candidates`
candidates per iteration with a **full catalog pass each**, so every rule you add multiplies
into the dominant cost term — and W1 is re-arming a sweep axis worth ~2.4 s in parallel.

Measure G7 before and after on your branch and **put the delta in your PR.** No session can
see the total; W2 lands last and carries it.

**Never raise `G7_SECONDS`. Never mark the test `slow`** — `.claude/rules/tests.md` records
that `-m "not slow"` currently deselects nothing, so the marker's only effect would be that
`gates.sh` and CI never run the gate. If the budget is exceeded, that is a **finding**:
`docs/PLAN.md:508` gives G7's fail consequence verbatim as *"re-allocate budget before rules
multiply cost."*

## Files

- `packages/engine/src/bt5/rules/catalog/d5_*.py` (new — name it
  `<brief_id>_<slug>.py` so `ls` answers "did we implement D5?")
- `packages/engine/tests/rules/test_d5_*.py` (new)
- `packages/engine/tests/rules/test_d1_restriction_sites.py` (new)
- `docs/decisions/2026-09-03-*.md`

Nothing else. Not `presets.py`, not `conftest.py`, not `data/`.

## Delegation

- **Every threshold** → `docs-miner`, verbatim, with `brief.md:LINE`. Never inline.
- **How an existing rule words an unavailable scoring path** → `Explore` on
  `d3_splicing.py`.
- **Gates** → `gate-runner`.
- **`/pre-pr` will fire `rule-auditor`** (opus/xhigh) by construction — a Spec's citations
  and thresholds changed. Budget for it. Its question is "does the cited source actually
  support the number?", so have the answer ready in `weight_provenance` and the docstring.

## Done means

- `bash scripts/gates.sh` → `ALL GATES PASSED`, exit 0.
- D5 ships with its paired test, its citations, and its Promoter-Calculator path reported
  `unavailable` with a reason — plus the issue filed for it.
- `test_d1_restriction_sites.py` exists and covers junction- and origin-spanning cases.
- `conflicts_with` declared and justified.
- Your G7 delta is in the PR.
- A decision record naming the enforcement class you chose for each of D5's five parts and
  **why**, plus what you rejected.
- Scientific impact: **non-"none"** — this changes emitted sequences. **Owner merge; you do
  not self-merge.**
