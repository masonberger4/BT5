# W5 — C7 and C8: codon composition, and the first rule that can read the native

*Copy this whole file into a fresh Claude Code session on the BT5 repo.*

**Launch:** permission mode `acceptEdits` · model opus · effort high (both the repo default,
so no override needed) · runs unattended.
**Do not run this in plan mode** — three other sessions are running in parallel and the gate
is the draft PR.

---

Two of ten 2.C rules exist (C1, C3). You are adding two more — and one of them is the first
rule in the catalog that can consume the native CDS PR #89 wired into `design()`.

## Read this first

Run `/bootstrap` — a fresh checkout has no `.venv`, and `gates.sh` exit **10** means BROKEN,
not a code failure. Then `CLAUDE.md`, then `.claude/rules/rules-catalog.md` (which governs
your lane specifically), then `docs/buildout/wave2/README.md`.

Your branch: **`claude/w5-c-composition`**. Cut it yourself, from a **freshly fetched** main:

```bash
git remote prune origin
git fetch -q origin main && git checkout -B claude/w5-c-composition origin/main
```

`git fetch --prune origin main` does **not** clear a stale ref — pruning is bounded by the
refspec, so that form prunes only `origin/main` (CLAUDE.md §7a).

## What to build, in this order

Pull both rows verbatim with `docs-miner` before you type a threshold — **never read
`brief.md` (63 KB) inline, and never paraphrase a number.**

### 1. C8 — rare-codon cluster preservation · `brief.md:84` · evidence B · **start here**

Retain **≥80%** of native rare-codon clusters when a native CDS is supplied (the brief cites
DeepCodon at 80–90%; ICOR/generic under 50%).

**This is newly buildable and nothing else in the catalog is.** `design(native_cds=...)`
landed in PR #89 as a keyword-only parameter, and C8 is the first rule that can read it. It
needs only host usage — which ships, four reference sets — plus the native CDS. No new data.

**The unavailability path is the common case and must read correctly.** With no native CDS
supplied, C8 reports `unavailable` with a reason that says **"not applicable here"**, not
"could not be evaluated" and certainly not a failure. Most runs will have no native. Read
how `c1_cai._unavailable` and `d3_splicing` word their reasons and match the register.

Note the design lane's own posture, from `docs/decisions/2026-09-02-ranking-increment.md`:
there is deliberately **no way to ask BT5 to invent a native** — "'do not optimize' compared
against a sequence BT5 itself designed is not a comparison." Your rule inherits that. Never
back-translate one to have something to compare against.

### 2. C7 — CSC, codon stability coefficients · `brief.md:83` · evidence A

Per-codon Pearson r between mRNA half-life and codon occurrence, as a **separate** objective
from elongation terms. The brief is explicit about the separation; do not fold it into a
CAI-shaped term.

**The constraint that shapes this rule: the CSC tables are data, and `data/` is a mutex
nobody holds this wave.**

So build C7 against `Services` with the table **absent**, so it ships reporting
`unavailable` with a reason, and **open the `approved:data-change` issue** naming the source
you would use, what organisms it covers, and how it would be validated.

**Do not inline a table into the rule file to dodge the mutex.** That is the failure mode
issue #78 was filed to prevent and `docs/decisions/2026-09-01-c1-cai-soft-band.md` argues
against at length: a reference set *is* the scientific work, and a table smuggled into a
rule file has no provenance, no `_provenance` block, no build script, and no reviewer who
can audit it down to the accession. `data/codon_usage/build_reference_set.py` is what a
defensible table looks like.

C1 is your model here in every respect: it computes for the hosts whose reference set ships
and reports `unavailable` for the rest, and that is *correct behaviour, not a bug* — it is
why the unavailable path exists.

## The two rules next door that you must NOT build

Named explicitly so no one wanders in and so `rule-auditor` is not surprised:

- **C9 (CPB/CPS)** — evidence **C**, contested sign. `brief.md:29` calls codon-pair-bias
  deoptimization "largely a CpG/UpA artifact" with a null fitness result cited, and the
  brief's own C9 row says **near-zero weight by default**, exposed only for deliberate viral
  deoptimization and always paired with explicit CpG/UpA reporting. If it is ever built, it
  **ships disabled**. Not this PR.
- **C6 (ENc/Nc)** — the brief says "Descriptive only — **never an objective**." A rule that
  scores it would contradict its own citation.

`CLAUDE.md`'s statement that rules resting on folklore ship disabled is not decoration; C9
is the clearest live example in the whole brief.

## Enforcement class — decide it deliberately

`CLAUDE.md` §3.5: hard constraints are **never** enforced by a penalty weight, and
`default_weight` must be 0.0 for `HARD_LATTICE`, `HARD_REPAIR` and `HARD_CHECK`.

Both C7 and C8 are composition objectives — SOFT is almost certainly right for both, which
means both need a **`weight_provenance` that explains the default**, and that sentence is
what `rule-auditor` will hold to the citation. "80% because the brief says 80%" is not
provenance; *why that number, measured how, on what* is.

C8's ≥80% is a **retention target**, not a hard floor: a protein whose native has no rare
codon clusters cannot fail it, and one whose clusters fall in immutable positions must not
be refused for it. Make sure the metric degrades sensibly at the edges and pin those edges
with tests.

## Wave contracts you must honour

- **Declare `conflicts_with`.** W2 is lighting up the conflict panel this wave; a rule that
  opposes an existing one and does not declare it is invisible there. C8's "keep the rare
  codons the native had" is in direct tension with C1's CAI band and C3's %MinMax — check
  what `c1_cai.py:416` already declares and mirror it where it is symmetric.
- **Do not reuse a `brief_ref`.** `_index_by_ref` (`score/presets.py:120-127`) raises when
  two rules claim one, and that breaks every preset at import. W2 is adding a
  `data_integrity` assertion that all three presets resolve cleanly, so a violation fails a
  gate — but it will fail *your* PR first.
- **No preset weights your rule, and that is fine.** The three shipped presets weight only
  `2.B1`, `2.C1`, `2.C3`, `2.D4`, `2.F2`. Under W2's policy a new SOFT rule reports its
  percentile at weight 0.0, flagged as outside the preset's claim. Note that all three
  presets already carry `_NATIVE_NOTE` on their C1 and C3 entries — so if C8 *should* count,
  **open an issue** naming the `brief_ref`, the proposed weight and the note.
  `score/presets.py` is W2's file, not yours.
- **`packages/engine/tests/rules/conftest.py` is READ-ONLY.** W3 and W4 are editing tests in
  the same directory. A helper only you need goes in your own test file. If you need a
  native-CDS fixture, define it in your own file — it is very likely a helper the whole lane
  will want later, and *that* is an issue to open after both merge, not a conftest edit now.
- **Do not add a data file.** `data/` is a mutex nobody holds this wave. This is the
  constraint C7 is designed around, not an obstacle to route past.

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

## How to add a rule here

**Use `/rule-add`.** It scaffolds the rule with its paired test, its provenance and a
resolvable `brief_ref`, and it knows the assertions in
`tests/data_integrity/test_rule_contract.py` your rule must satisfy.

Every rule ships with citations, an evidence badge, `last_verified`, and — if SOFT — a
`weight_provenance` explaining its default.

**Seed every RNG explicitly** with `np.random.default_rng(seed)` if either rule samples.
Global `np.random.*` and any stdlib `random` import are banned in engine source and the test
tree, and CI greps for both (`ci.yml:76-89`).

## Files

- `packages/engine/src/bt5/rules/catalog/c7_*.py`, `c8_*.py` (new)
- `packages/engine/tests/rules/test_c7_*.py`, `test_c8_*.py` (new)
- `docs/decisions/2026-09-03-*.md`

Nothing else. Not `presets.py`, not `conftest.py`, not `data/`, not `c1_cai.py` — **W1 owns
`c1_cai.py` this wave** and is editing it.

## Delegation

- **Every threshold and every effect size** → `docs-miner`, verbatim, with `brief.md:LINE`.
- **How C1 reaches its reference set through `Services`, and how it words unavailability** →
  `Explore` on `c1_cai.py` (read-only for you).
- **Gates** → `gate-runner`.
- **`/pre-pr` will fire `rule-auditor`** (opus/xhigh) by construction. It caught a real
  resolver defect in S6's reference sets *before* that PR opened, which is why S6 ran it
  early. Consider doing the same with C8's ≥80% retention metric — it is the number most
  likely to be defensible in the brief and wrong in the implementation.

## Done means

- `bash scripts/gates.sh` → `ALL GATES PASSED`, exit 0.
- C8 computes against a supplied native CDS and reports **"not applicable"** — in that
  register — when none is supplied.
- C7 ships reporting `unavailable` with a named reason, with **no table inlined anywhere**,
  and the `approved:data-change` issue is filed naming the source.
- `conflicts_with` declared and justified against C1 and C3.
- Your G7 delta is in the PR.
- A decision record naming the enforcement class for each rule and why, what you rejected,
  and an explicit line saying C9 and C6 were **not** built and why.
- Scientific impact: **non-"none"** for C8 (it changes emitted sequences); C7 alone would be
  "none" while it reports `unavailable`, so state the two separately rather than averaging
  them. **Owner merge; you do not self-merge.**
