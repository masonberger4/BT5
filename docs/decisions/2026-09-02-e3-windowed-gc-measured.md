## 2026-09-02 — E3: no windowed GC threshold gates, including IDT's own

**Decided:** `e3_windowed_gc` ships `Enforcement.SOFT` and `Evidence.CONTESTED`. It
reports three windowed readings — IDT's 100 bp / 30% floor, GenScript's 100 bp 25-65%
band, and Twist's 50 bp 10%/90% trigger — and **enforces none of them**.

This is the one decision in this session where the brief and the repo's own measurements
disagree, so the reasoning is set out in full.

### What `brief.md:140` says, graded A

> Windowed GC 50 bp: hard-fail any window <10% or >90% (Twist High-Complexity trigger);
> warn <25% or >75%. Windowed GC 100 bp: warn outside 25-65% (GenScript GenTitan — the
> only vendor publishing a windowed rule)

### What `docs/design/vendor-gc-calibration.md` measured

18 probes, two vendors, 2026-08-28. Eight carried **one 50 bp window swept 4% → 96% GC**
on a byte-identical 50% background. All eight were green at IDT and Standard at Twist —
**16 accepting verdicts out of 16**, over a range far wider than the 10%/90% trigger. The
document's conclusion: *"A single extreme 50 bp window is irrelevant to both vendors."*

So the 50 bp hard-fail is a threshold no vendor was observed to apply. Encoding it at
`HARD_REPAIR` would make BT5 refuse to emit constructs both vendors manufacture without
comment. That is the same defect `docs/design/repeats.md` §2 records for Twist's published
*repeat* trigger, and the same shape as the E4 row `brief.md:141` struck through for being
below the chance floor. **Decided: reported, never enforced.**

### The second step, which is easy to miss

The obvious repair is to keep the one windowed rule a vendor *does* state — IDT's
remediation text, *"a window of 100 bases ... Redesign this region to have a GC content
greater than 30%"* — and gate on that instead. The calibration even recommends it: *"Keep
a windowed floor, at IDT's geometry."*

**That is still wrong, and the calibration contains the disproof three paragraphs
earlier.** Explaining why the LOC probes passed, it says: *"a 50 bp window at 4% GC,
diluted by 50 bp of 50% background, is ~27% across 100 bp — barely under the 30% floor,
**scoring a few points and staying green**."*

Those probes **tripped IDT's own windowed floor and were accepted anyway**. IDT's floor is
a contributor to a complexity *score*; what denies an order is that score reaching 24, and
for GC the score is driven by **global** GC (`score = 1.40 × GC% − 83.8`, denial at 77%).
A `HARD_REPAIR` rule built on the 30% floor would refuse constructs IDT scores green — the
same mistake as the 50 bp trigger, one step further in.

This was caught by building the rule as `HARD_REPAIR` first and noticing that a fixture
reproducing the calibration's own LOC probe (`CLEAN + "AT"*25 + CLEAN`) came back
`passes=False`. The rule was refusing a probe the document records as accepted.

**Decided: the whole rule is SOFT.** Global GC is what gates; `e2_gc_band` owns that.
`default_weight = 0.4` carries the whole of this rule's influence, and
`steering_weight = 0.1` sits below `e2_gc_band`'s 0.5 and `f5_at_window`'s 0.25 so three
GC steering terms do not count one preference three times.

### Rejected

- **Encoding `brief.md:140` as written.** Refuses 16 of 16 measured-acceptable probes.
- **`HARD_REPAIR` on IDT's 100 bp / 30% floor.** See above — the calibration's own LOC
  results are the counterexample.
- **Dropping the refuted numbers entirely.** They stay, at the `note` tier, with the
  refutation in the message and a `Citation(sign="refutes")` on the spec. Deleting them
  would leave the *absence* of a rule with no recorded reason, and the next session would
  re-add it from `brief.md` — which is exactly how a struck-through threshold comes back.
- **Editing `e2_gc_band`.** E2 carries `brief_ref = "2.E2/2.E3"` and the calibration says
  E2 as shipped is *"wrong in every dimension at once"* — but that rework is **X7**, a
  change to what the app refuses to build, and an owner call. This rule adds nothing hard
  and takes nothing away, so it cannot make E2 worse. **Follow-up for the owner:** X7, and
  narrowing E2's `brief_ref` to `2.E2` once this rule owns `2.E3`.
- **`Evidence.EVIDENCE_BACKED`.** `brief.md:140` grades the row A, but that grade attaches
  to published numbers the measurement contradicts. CONTESTED is the honest badge.

### Shared geometry

`window_gc` and `merge_regions` live in `f5_at_window` and are imported here. Both rules
turn a per-window array into regions, and the wrap case in `merge_regions` — a stretch
crossing the origin arriving as two runs — is subtle enough that two copies would drift.
It is the bug F5 shipped with for one commit.

**Evidence:** `brief.md:140, 141`; `docs/design/vendor-gc-calibration.md` — "Method",
"Three findings that change what BT5 should encode", "The combined result", "What this
says about E2 as shipped — X7, settled"; `docs/design/repeats.md` §2.
`pytest packages/engine/tests/rules/test_e3_windowed_gc.py` 28 passed;
`bash scripts/gates.sh` ALL GATES PASSED, 1410 passed.

**Where:** branch `claude/s4-rules-liabilities`, session S4 of the six-way buildout.
