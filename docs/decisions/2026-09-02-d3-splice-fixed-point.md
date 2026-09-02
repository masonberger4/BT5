## 2026-09-02 — D3 splicing: the fixed point is the solver's, and HARD_REPAIR's oracle promise is unfunded

*This file doubles as the D3 design note required by `docs/buildout/s4-rules-liabilities.md:130-138`
— the review gate that replaces plan mode for this session. It was written and committed
**before** `d3_splicing.py` was implemented, and is its own commit so it can be read
against an empty diff.*

**Decided:** `d3_splicing` ships as a directional, per-slot detector declaring
`enforcement = SOFT` (floor) with `enforcement_for` escalating to `HARD_REPAIR` on the
two modalities the brief grades hard, `repair = RepairPolicy.FIXED_POINT`, and
`lattice_terms() -> None`. Five sub-decisions, each with the alternative it beat:

### 1. The fixed point is over the solver's re-targeting loop, not a loop inside the rule

The rule is a **pure, idempotent detector**. It does not iterate and it does not repair.
`solver/repair.py` owns the loop: `RuleSet.breach_finder()` (`solver/catalog.py:192-207`)
wraps the whole catalog into one `BreachFinder`, the solver mutates codons itself, and
re-calls the finder on every candidate — up to `max_candidates=256` per iteration. There
is **no per-rule repair callback** anywhere in the seam.

So "iterate to a fixed point" (`brief.md:104`) is discharged by one ClassVar. The fixed
point is over *the set of breaches `evaluate()` emits on the assembled construct*, and it
is reached when the rule stops emitting them.

**Rejected:** *a `while` loop inside `evaluate()` that removes donors until none remain.*
This is the obvious approach and it is wrong three times over: `evaluate()` is required to
be pure (`core/spec.py:253`), it receives a `Construct` and has no mutation space or codon
table to work with, and the solver would still call it 256 times per iteration — so the
loop would run inside every candidate evaluation. `CLAUDE.md` §3.6 exists to make the
policy declarative precisely so no rule hand-rolls this.

### 2. Convergence and non-convergence are the solver's, and non-convergence is silent

Detection is not mine either. `converged` is true **only** when `_select_target` returns
`None` — `stop_reason == "exhausted_targets"` (`solver/repair.py:601`). The two failure
exits are `"stagnation"` (`STAGNATION_TOLERANCE = 100` consecutive non-improving
iterations) and `"iterations"` (`MAX_ITERATIONS = 1000`).

The behaviour on non-convergence is the part worth flagging: `solver/pipeline.py:122`
passes `raise_on_infeasible=False`, so `RepairNotConverged` is **never raised in
production**. Repair returns its best attempt with `remaining` non-empty and
`converged=False`, and the pipeline relies on `verify_construct` to be the refusal.

### 3. The load-bearing gap: `verify_construct` cannot see a splice donor

`Enforcement.HARD_REPAIR` promises, verbatim at `core/spec.py:36-38`, *"enforced by
Tier-B repair and then PROVEN by the independent validator, which refuses to emit on
failure."* **For this rule that promise is unfunded**, and it is not something I can fund
from the rules lane:

- `RuleSet.forbidden()` short-circuits at `solver/catalog.py:149-150` —
  `if spec.enforcement is not Enforcement.HARD_LATTICE: continue`. A `HARD_REPAIR` rule's
  `lattice_terms().forbidden` therefore **never reaches invariant I6**.
- `OracleBounds` (`solver/catalog.py:70-93`) carries exactly one field, `gc_bounds`, and
  `_gc_bounds()` hard-matches `spec.id != "e2_gc_band"`.
- `verify.py` may not import `bt5.rules` at all (AST-enforced, `verify.py:12-14`).

Net: if repair stops on stagnation, nothing refuses, and a construct ships with a cryptic
donor `converged=False` already recorded. That is a real hole in what the app refuses to
build — a `CLAUDE.md` §7b owner call — and it is an **M1/oracle-lane** change
(`solver/catalog.py` or `verify.py` + `approved:oracle-change`), not an edit this session
may make. **Filed as a follow-up rather than worked around.** It is the same shape as the
already-recorded `max_repeat` gap at `solver/catalog.py:86-90`.

**Rejected:** *declaring `enforcement = HARD_LATTICE` and returning the donor motifs from
`lattice_terms()` so I6 does prove them absent.* It buys the oracle backstop and loses
more than it buys. `LatticeTerms.forbidden` is closed under reverse complement by the
solver (`core/spec.py:188-193`), and splicing is directional: forbidding `CTTACC` because
`GGTAAG` is a donor refuses designs over a site that cannot fire on the transcribed
strand — exactly the argument `d4_internal_polya.py:10-15` already makes for `AATAAA` vs
`TTTATT`. It also cannot express a score threshold, and it makes `FIXED_POINT` moot,
since `HARD_LATTICE` rules never enter `repair_specs()`. Choosing it would satisfy the
oracle by contradicting `CLAUDE.md` §3.6, which is not a trade this lane gets to make.

**Rejected:** *setting `enforcement = HARD_LATTICE` on the ClassVar while returning
`HARD_REPAIR` from `enforcement_for`.* This genuinely would get both — `forbidden()` reads
the ClassVar, `repair_specs()` reads the method — but it is a trick that works by the two
call sites disagreeing, and it still forbids the reverse complements. A seam bug is not a
design.

### 4. MaxEntScan cannot be shipped, so the scored path reports `unavailable`

`brief.md:288` says MaxEntScan redistribution is *"ambiguous"* and recommends training a
max-entropy model on GENCODE instead; `brief.md:295` lists SpliceAI (CC BY-NC) among
*"licensing landmines to avoid."* No score matrices exist in this repo, and `data/` is
S6's mutex.

So the 9-mer donor / 23-mer acceptor scored path is **model-injected**: with no model it
reports the objective unavailable via the `c1_cai.py:490-525` pattern — `raw_score` NaN,
`passes=True`, `fixable_by_codon_choice=False`, and the reason in
`detail["unavailable_reason"]`. NaN rather than 0.0 because every value in the score range
is a real, meaningful number of bits, and 0.0 would read as an affirmative "no site here"
that nobody measured. `passes=True` because nothing about the construct failed.

The **enforceable core needs no licensed data**: the literal donor motifs from
`brief.md:102` (`GGTAAG`, `GGTGAG`, and the `AN|GT(A/G)AG` consensus, i.e. `ANGTRAG`).
Those are what repair actually chases, and they are what the fixed point iterates over.

**Rejected:** *vendoring MaxEntScan's matrices.* Licensing, and it is `data/`.
**Rejected:** *a plausible hand-rolled stand-in for the score.* This is the one thing the
whole honesty apparatus exists to prevent — the `services` fixture ships `fold=None`
rather than a stub for exactly this reason (`tests/rules/conftest.py:93-100`).

### 5. Thresholds, gating and the V5 case — all verbatim, none superseded

`brief.md:101-105`, confirmed clean: **no `~~strikethrough~~`, no "corrected",
"superseded" or "provisional" anywhere in lines 88-132.** (The two such markers in the file
are `brief.md:141` for E4 and `brief.md:158` for F4 — neither in the D section.)

- Donor: MaxEntScan 9-mer (3 exonic + 6 intronic), **flag > 3 bits, hard-constrain > 6-8**.
- Acceptor: MaxEntScan 23-mer (20 intronic + 3 exonic), **flag > 3 bits**. No hard cutoff
  is given for the acceptor — the `> 6-8` band is on the donor line only, and the rule
  does not invent one.
- Acceptor context: an AG preceded within **5-40 nt** by a **≥10-nt window ≥80%
  pyrimidine** and a branch-point-like **`YTNAY` 18-40 nt** upstream.
- `GTNNG` is a 5-mer, and `brief.md:90`'s occurrence budget says **≤5-mer → soft only**.
  It is reported, never hard, and never placed in `forbidden`.
- Gating (`brief.md:101` header, `brief.md:223` matrix): off entirely for
  `BACTERIAL_EXPRESSION`, `IVT_MRNA`, and the yeast hosts (heterologous CDS);
  `HARD_REPAIR` for `LENTIVIRAL` (*titer + safety*) and `GENOME_INTEGRATED` (*fusion
  transcripts*); `SOFT` (*warn*) for `PLASMID_TRANSIENT`, `PLASMID_STABLE`, `AAV`.
- **V5 (H, evidence A)** — `brief.md:105`: the standard V5 encoding contains `G|GTAAG` and
  spliced in **17/17** genes tested; **13/17** randomly chosen genes showed aberrant
  splicing from vector/tag context. The scan detects the V5 **peptide**
  (`GKPIPNPLLGLDST`, `brief.md:190`) in the translated CDS rather than any one nucleotide
  encoding — the whole point is that the *standard encoding* is the liability, so matching
  a nucleotide string would miss every re-encoding and match nothing after repair. A donor
  overlapping those codons is escalated.

**Honesty note carried into the rule's docstring**, `brief.md:321`: *"No published dataset
cleanly quantifies the titer cost of a cryptic splice donor inside a therapeutic ORF —
every quantified case is confounded. Say so."*

### The test that fails under a single pass

`_breach_key` is `(spec_id, interval.start, interval.end)` (`solver/repair.py:343-345`),
and `SINGLE_PASS` retires that key after one attempt (`repair.py:569-571`). So the §3.6
failure is reproducible exactly, not approximately: a repair that turns a **strong** donor
into a **weak** donor *within the same reported window* keeps the breach count at 1 while
the magnitude-sum strictly falls, which `_accepts` (`repair.py:270-312`) accepts — and the
new donor lands on the **same interval**, hence the **same key**, which `SINGLE_PASS` has
already retired. The construct ships with a donor the repair itself created.

The test drives `repair()` directly with a hand-built `policies={...: RulePolicy(...)}`
(the convention at `tests/solver/test_repair_seam.py:180-190`), once with `SINGLE_PASS`
and once with `FIXED_POINT`, and asserts the first leaves a donor and the second clears
it. Graded magnitudes are therefore load-bearing, not decoration.

**Evidence:** `brief.md:90, 101-105, 190, 208, 223, 288, 295, 321`; `core/spec.py:36-38,
106-117, 188-193, 253, 261-290`; `core/context.py:83`; `solver/repair.py:270-312, 343-345,
405-428, 489-580, 601, 624-628`; `solver/catalog.py:70-93, 138-154, 149-150, 192-207,
218-244, 246-255`; `solver/pipeline.py:122`; `verify.py:12-14, 165-169`;
`d4_internal_polya.py:10-15`; `c1_cai.py:490-525`. Baseline before any edit:
`bash scripts/gates.sh` → `ALL GATES PASSED`, 1240 passed.

**Also found, not fixed (stale doc citations, `.claude/` is not this lane):**
`.claude/rules/rules-catalog.md:63` and `.claude/skills/rule-add/SKILL.md:19` both cite
`solver/repair.py:174` as `FIXED_POINT`'s default; line 174 is prose inside `localize()`'s
docstring and the real default is `repair.py:417`. `.claude/rules/rules-catalog.md:57`
cites `core/spec.py:259` for `strand_of_interest`; that line is blank — the field is
`core/context.py:83` and the accessor is `strand_for()` at `core/spec.py:261-290`.
`LatticeTerms`' own docstring (`core/spec.py:191`) says to read `packaged_strand` from the
context; no such attribute exists.

---

## Addendum, same day — both optional scans ship OFF, and why that is not a perf hack

Found while implementing, after the note above was committed. Recorded here rather than
in a second file because it is the same decision: what `d3_splicing` actually enforces.

**Symptom.** With `report_coarse` and `scan_acceptors` defaulting on, the first
end-to-end design test (`tests/design/test_skeleton.py::TestEndToEnd`) went from passing
to not finishing in 500 s. `evaluate()` cost 8.7 ms on a 5 kb circular construct and
emitted **116 breaches, of which 2 were strong**. It runs once per candidate, up to 256
per repair iteration, so 114 breaches of noise were being re-derived tens of thousands of
times and were also dominating `_aggregate`'s breach count — the cross-rule currency the
solver steers on.

**Decided, and the scientific reason comes first in both cases:**

1. **`report_coarse=False`.** `GTNNG` is 5 nt with two free positions, so
   P = 1/4³ = 1/64 per position — **~78 hits on a 5 kb plasmid of random sequence**, and
   the measurement produced exactly 78. A screen that fires 78 times by chance on a clean
   construct reports noise, not findings. `brief.md:90` already grades a ≤5-mer soft-only;
   this goes one step further and makes it opt-in. This is the same failure mode as the
   E4 row `brief.md:141` struck through on 2026-08-28 — *a threshold below the chance
   floor* — arrived at independently, which is why it is called out rather than quietly
   defaulted.
2. **`scan_acceptors=False`.** The stronger reason is not cost. `brief.md:103` **flags**
   acceptors and deliberately states no hard cutoff, but `enforcement_for` is per-slot,
   not per-breach — so an acceptor breach emitted by this rule inherits the donor half's
   `HARD_REPAIR` on a lentiviral slot and the solver starts *chasing* sites the brief
   never authorised constraining. Emitting them by default would let a rule enforce a
   constraint its own evidence row declines to state.

**Rejected:**
- *Keeping both on and accepting the cost.* Not a cost question. On (2) it changes what
  the app refuses to build; on (1) it buries the two real findings under 78 chance hits.
- *`fixable_by_codon_choice=False` on acceptor breaches to route them to `advisory`.*
  This is the tempting one-line fix and it is a lie: a codon change genuinely can remove
  a cryptic acceptor inside a CDS, and `Breach`'s own docstring (`core/spec.py:154-166`)
  says the field exists because only the rule author knows the answer and nobody
  downstream can recover it. Mis-stating it to get a routing side effect is exactly the
  abuse the field was given no default to prevent.
- *A `while` loop or an internal cache in `evaluate()`.* `evaluate()` must stay pure
  (`core/spec.py:253`); it is called on candidate constructs that are then discarded.
- *Splitting the acceptor scan into a `REPORT_ONLY` rule now.* That is the right home —
  one enforcement class per rule is the contract — but it is a second rule against the
  same `brief_ref` row and a shape no catalog rule currently has. **Filed as a follow-up
  for the owner**, not decided unilaterally at the end of a session.

**Also fixed:** the scan moved from `range(len(text))` to `str.find`. Same results,
0.54 ms instead of 8.7 ms at the shipped defaults (16x), which matters because
`cost_class = "cheap"` is a promise this rule makes to the repair loop.

**Evidence:** measured on a seeded 5 kb random construct, `np.random.default_rng(0)`:
8.7 ms / 116 breaches / 2 strong before, 0.54 ms / 3 breaches / 2 strong after.
`tests/design` 16 passed in 9.18 s, having previously exceeded 500 s.

**Where:** branch `claude/s4-rules-liabilities`, session S4 of the six-way buildout.

---

## Addendum, same day — `ANGTRAG` is seven characters but not a 7-mer

`/verify-provenance` returned SUPPORTED for D3 with one **UNSUPPORTED** threshold, and it
was right.

The consensus `AN|GT(A/G)AG` was tiered `MAG_STRONG` on the reasoning that it is 7 nt and
`brief.md:90` says *"≥7 nt → hard (cheap)"*. That reads the budget as a rule about
**length**. It is a rule about **expected occurrences** — the line states them explicitly
(5-mer ≈ 2.9, 6-mer ≈ 0.73, 7-mer ≈ 0.18 per 1.5 kb) — and two of the seven positions are
degenerate:

```
A  N  G  T  R  A  G      N is 4-fold, R is 2-fold
P = (1/4)^5 x 1 x (1/2) = 1/2048      a true 7-mer = 1/16384
effective k = log4(2048) = 5.5        8x commoner than a 7-mer
```

Measured on random 50% GC DNA: **2.85 occurrences per 1.5 kb**, which is `brief.md:90`'s
own figure for a **5-mer (2.9)** — and that line puts a ≤5-mer at **SOFT ONLY**. The
literal blacklist measured 0.70 against the brief's 0.73 for a 6-mer, so it sits exactly
where the budget expects.

**Decided:** a consensus-only match reports at `MAG_FLAGGED`. `MAG_STRONG` — the tier that
drives `passes=False` and Tier-B repair — is reserved for the literal 6-mers `GGTAAG` and
`GGTGAG`, which `brief.md:90` permits to be hard *"only if strongly evidenced"* and which
the V5 17/17 evidence-A result supplies, and for any donor inside a V5 tag.

**What this changes in behaviour:** D3 now hard-constrains roughly a quarter as many sites.
The consensus is still reported at every hit, so nothing is hidden — it stops *refusing*
constructs over a pattern that arises ~3 times per 1.5 kb by chance. Pinned by
`TestDonorTiers::test_the_degenerate_consensus_is_not_a_seven_mer`, which asserts the
arithmetic so a later "it's 7 nt, promote it" cannot land unchallenged.

**Rejected:** *keeping `MAG_STRONG` because the brief lists the consensus alongside the
literal motifs.* Listing is not grading; `brief.md:90` is the grading rule and it is
explicitly about occurrence rate.

**Evidence:** `brief.md:90, 102`; measured over 200 seeded 1.5 kb sequences,
`np.random.default_rng(0)`. `bash scripts/gates.sh` ALL GATES PASSED, 1445 passed.
