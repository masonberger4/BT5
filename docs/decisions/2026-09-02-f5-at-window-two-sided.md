## 2026-09-02 — F5 is a propagation rule, so the vendor GC calibration does not govern it

**Decided:** `f5_at_window` ships `HARD_REPAIR`, `Direction.BAND` with
`band = (0.45, 0.60)` and a hard-fail outside `(0.35, 0.65)` GC per 100 nt window,
reporting the binding side per region.

### The distinction that decides everything else

`docs/design/vendor-gc-calibration.md` measured 18 probes at two vendors and concluded
that **no vendor gates on a windowed GC bound** — sixteen of sixteen probes carrying one
50 bp window swept from 4% to 96% GC were accepted by both. It is tempting to read that
as retiring F5 too.

It does not. That calibration measured **what a vendor refuses to synthesise**. F5 is
about **what survives propagation in E. coli**: `brief.md:151` scopes section 2.F to
*"every construct that passes through a cloning host"*, and the evidence is a toxicity
measurement — *"Toxic horizontally-acquired E. coli genes are 63-68% AT vs a non-toxic
55% control"* (`brief.md:159`, grade **A**). Twist will happily manufacture a fragment
that then kills the strain carrying it. Two different questions about the same number, and
the module docstring says so explicitly so the next reader does not collapse them.

### The conflict is declared, not resolved

`brief.md:159` puts it in bold: *"This directly conflicts with vendor GC ceilings."*
Suppressing cryptic AT-rich promoters pushes GC up; the vendor floor and IDT's fitted
denial above ~77% global GC push it down. `e2_gc_band` already carries the same Nature
citation with `sign="qualifies"` and names `f5_at_window` in its `conflicts_with`; this
rule now names it back, so the conflict panel sees a declared pair rather than two rules
quietly fighting inside the solver.

**Rejected:** *picking a winner and widening one of the two bands to remove the overlap.*
Neither side is wrong; the brief asks for a two-sided band precisely because both
mechanisms are real. Resolving it in a rule file would bury an owner-level trade-off in a
module constant.

### `binding_side` is a requirement

`brief.md:159` ends *"and show which side is binding per window."* A `|deviation|`
scalar cannot say whether a window needs GC raised or lowered, and with two different
mechanisms behind the two sides that is the only actionable part of the finding. Every
breach carries `detail["binding_side"]`, and `Evaluation.binding_side` — a field that
exists on the type for exactly this — reports the worst window's side. Pinned by a test
that puts an AT-rich region and a GC-rich region on one construct and asserts both sides
appear.

### One breach per region, and the wrap bug that found

A 5 kb plasmid at 30% GC has ~5,000 offending windows. Emitting one breach each would
make this rule's breach **count** — the currency `_aggregate` steers on
(`solver/repair.py:253-267`) — swamp every other rule in the catalog for what is one
contiguous problem. So offending windows are merged into regions and each region is
reported once, at its worst window.

Building those regions from a `np.diff` edge scan over the mask is correct on a line and
wrong on a circle: **an offending stretch that crosses the origin arrives as two runs**,
one ending at the last window start and one beginning at the first. The rule reported it
as two regions and named two "worst" windows for one problem. Caught by
`TestRegionMerging::test_one_contiguous_bad_stretch_is_one_breach`, which is circular by
default; fixed by joining the first and last runs when the mask is set at both ends and
the construct is circular.

### Vectorised, because HARD_REPAIR is expensive

`breach_finder` calls this once per candidate, up to 256 per repair iteration. A Python
sweep over 5,000 window starts cost **4.65 ms**; a cumulative sum over the tripled text
plus a masked edge scan costs **1.42 ms**, and the region loop now runs over regions
rather than windows. That matters because an earlier version of `d3_splicing` in this
same branch took `tests/design` from 9 seconds to not finishing inside 500.

**Rejected:** *leaving it in pure Python and finding out in CI.* The failure mode is a
timeout, not a red assertion, which is the expensive kind to diagnose.

**Evidence:** `brief.md:151, 159`; `docs/design/vendor-gc-calibration.md` "The combined
result"; `e2_gc_band.py:128-134, 151`; `solver/repair.py:253-267`.
`pytest packages/engine/tests/rules/test_f5_at_window.py` 31 passed;
`bash scripts/gates.sh` ALL GATES PASSED, 1382 passed.

**Where:** branch `claude/s4-rules-liabilities`, session S4 of the six-way buildout.

---

## Addendum, same day — the audit found two real defects here

`/verify-provenance` over the five new rules returned **UNSUPPORTED** for this one. Both
findings were checked and both were right.

### 1. `SINGLE_PASS` was unsafe — this rule needs `FIXED_POINT`

The original justification was that *"a windowed statistic cannot recreate itself the way
a splice donor can, because the change is monotone in the quantity being measured."* That
holds for a **one-sided** bound — `e3_windowed_gc` is a floor and is genuinely safe under
`SINGLE_PASS` — and it is **false for a two-sided band**.

Consecutive 100 nt windows overlap by up to 99 nt, so raising GC to clear a window below
`hard_lo` raises it in every overlapping window too, and can push a neighbour that was
already near `hard_hi` straight past it. Measured on a fixture with a 0% GC window beside a
64% one: **lifting the low side created 79 new windows breaching the high side.**

That is exactly `CLAUDE.md` §3.6's test — *"if your repair can create new instances of what
it removes"* — so `FIXED_POINT` is mandatory, not a preference. This rule is now the second
in the catalog to declare it, alongside `d3_splicing`. Cost: `tests/design` 9.2 s → 14.7 s,
since a `FIXED_POINT` breach is never retired. Pinned by
`TestRepairPolicy::test_clearing_the_low_side_can_breach_the_high_side`, which asserts the
overshoot rather than describing it.

### 2. `EVIDENCE_BACKED` let the unevidenced half of the band wear the evidenced half's grade

`brief.md:159` grades the row **A**, and the first version took that at face value. But the
grade covers the measurement, and **the measurement is one-sided**: toxic
horizontally-acquired *E. coli* genes run 63-68% AT against a 55% AT control. That supports
a **floor** on GC. It says nothing whatever about a GC **ceiling**.

The 65% upper bound comes from `brief.md:159`'s own resolution instruction — *"resolve as a
two-sided band ... hard-fail outside 35-65%"* — not from the cited paper. And BT5's own
18-probe ladder cuts against a hard ceiling from the other side: 80% global GC was
**accepted by Twist as Standard** and produced no windowed finding at IDT.

Badging the whole band `EVIDENCE_BACKED` is the same trap `e3_windowed_gc` is `CONTESTED`
for, and grading the two rules differently on the same structure was an inconsistency the
audit caught. **Now `CONTESTED`**, with a third `Citation(sign="refutes")` naming the
calibration so the ceiling's weaker footing is on the record rather than in a comment.

**The band itself is unchanged.** The brief instructs it and F5 is a propagation rule, not
a synthesis one; what changed is the badge, the repair discipline, and the honesty of the
citation set. Enforcement staying `HARD_REPAIR` on an instruction rather than a measurement
is worth an owner's eye.

**Evidence:** `brief.md:159`; `docs/design/vendor-gc-calibration.md` "The combined result"
rows for `GLB_gc80`; `CLAUDE.md` §3.6. `bash scripts/gates.sh` ALL GATES PASSED, 1445
passed.

---

## Addendum, same day — the same unfunded oracle promise D3 has, undisclosed here until now

Found by an independent code-review pass, not by me: `d3_splicing`'s decision doc names a
real gap — `HARD_REPAIR`'s contract (`core/spec.py:36-38`) promises the construct is
*"PROVEN by the independent validator, which refuses to emit on failure"*, but that promise
is unfunded for a rule the oracle does not know how to check. F5 is the **other** new
`HARD_REPAIR` rule in this branch and has the identical gap, and I said nothing about it
here. That is an inconsistency in how honestly the two decision docs represent the same
class of risk, not two different facts.

**Confirmed independently**, not taken on the reviewer's word:
`solver/catalog.py:269-296`'s `_gc_bounds()` hard-matches `spec.id != "e2_gc_band": continue`
— it reads **only** `e2_gc_band`'s resolved band into `OracleBounds.gc_bounds`, which is the
one number `verify.py`'s I7 invariant checks. F5's windowed 35-65% band never reaches it.

**Scope, so this is not read as bigger than it is.** This is a **pre-existing, repo-wide
pattern**, not a defect this branch introduced: `e5_synthesis_repeats` and
`f1_direct_repeats` have the identical gap, and `OracleBounds`'s own docstring
(`solver/catalog.py:70-93`) already discloses it for `max_repeat`. Nothing about F5's own
logic is wrong. What was missing is the same sentence D3 already has: if Tier-B repair
stops on stagnation or the iteration cap — production always passes
`raise_on_infeasible=False`, per `d3_splicing`'s own decision doc — a construct violating
F5's hard band can still ship, silently, because nothing downstream re-derives the check.

**Not fixed here, for the same reason D3's isn't**: closing it means feeding a second
rule's band into `OracleBounds`, which is `solver/catalog.py` — the M1/oracle lane, not
this one, and it needs `approved:oracle-change`. **Filed as the same follow-up D3's doc
already raises**, now naming both rules rather than one.

**Evidence:** `core/spec.py:36-38`; `solver/catalog.py:70-93, 269-296`; `verify.py`'s I7;
`docs/decisions/2026-09-02-d3-splice-fixed-point.md` §3.
