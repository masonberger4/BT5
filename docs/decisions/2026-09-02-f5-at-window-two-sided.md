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
