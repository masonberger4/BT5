# 2026-09-02 — The GC lean must aim inside the band, not at its edge

**Lane:** M3 score (S1). **Branch:** `claude/s1-ranking-increment`.
Follows `2026-09-02-ranking-increment.md`; written after `main` merged S3/S4/S5/S6
into this branch and the design path stopped terminating.

## What happened

`main` brought nine new catalog rules. One of them, **`f5_at_window`** (S4), is a
two-sided GC band — 45-60% GC per 100 nt, hard-fail outside 35-65% — enforced
**HARD_REPAIR**.

The gallery's diversity axis was a flat per-codon preference for GC or for AT.
The Tier-A DP minimises a SUM, so a flat preference is not a nudge: it drives
every codon to the extreme it can reach. `oracle_bounds()` hands the solver e2's
band, currently 0.28-0.77, so Tier A emitted a ~70% GC CDS — outside f5's window
along the sequence's whole length. Tier B was then asked to repair a sequence
that was out of band everywhere, and did not converge.

Measured on the reference fixture at 140 aa, machine idle:

| sweep axis | before | after |
|---|---|---|
| unsteered | 0.6 s | 0.6 s |
| `repeat_avoidance` | 0.6 s | 0.6 s |
| `gc_lean_gc` | **>70 s, aborted** | 0.1 s |
| `gc_lean_at` | **>70 s, aborted** | 0.5 s |

CI's `python-tests` ran for **1 h 53 m** on the merged head before it was
cancelled. The lane's own suite went from 75 s to not finishing.

## Decided: the leans aim at targets inside the band

`lean_targets(gc_bounds)` returns the 35% and 65% points of the declared band,
and the cost is `|gc(codon) - target|`. Three things follow, and the first is the
one that matters:

**The term changes shape, not just size.** `|gc - target|` is bounded and
non-monotonic, so the DP settles AT the target instead of running away from it.
Shrinking a flat preference would only have moved the hang, not removed it.

**The targets stay clear of the band's edges,** which is where the tighter
windowed rules live. On e2's 0.28-0.77 they land at 45.2% and 59.9% GC — inside
f5's preferred 45-60% window. `test_the_shipped_band_aims_inside_f5s_window`
pins that concretely, and `test_a_lean_aims_inside_the_band_never_at_its_edge`
pins the general property.

**The axis still separates the designs,** which was the whole reason for it:
measured 40.2% GC for the AT lean against 60.0% for the GC lean, with the
unsteered solve at 44.7%. A codon's GC is quantised to {0, ⅓, ⅔, 1} and both
targets sit between ⅓ and ⅔, so those are the two values the leans choose
between — which is also why the first version of the paired test proved nothing
and had to be rewritten over a codon set spanning both.

## Rejected

**Dropping the GC-lean axes.** The fastest fix, and it would take the panel back
to a CAI-versus-repeats mixture — which measurement already showed clusters,
because both are dominated by the same relative-adaptiveness table. That fails
G4 for a reason unrelated to the constraint set.

**Reading f5's constants from `rules/catalog/f5_at_window.py`.** It would aim the
targets exactly, and it couples M3 to one rule's internals. Deriving from the
band the design already passes in costs nothing and cannot go stale when a rule
retunes. The repo has a scar from the duplicate-vendor-namespace bug that says
what a second copy of someone else's constants is worth.

**Lowering `G4_MIN_PAIRWISE_DISTANCE` or the G7 bar.** Neither was ever on the
table; both are recorded here only because a hang makes them tempting.

## Still red, and not this lane's to fix

`test_a_500_residue_design_meets_the_g7_budget` fails at **23.5 s against PLAN's
10 s bar**. The breakdown is sweep 20.2 s, null 3.1 s, everything else 0.2 s, and
the sweep is three full solves dominated by Tier B.

`RuleSet.breach_finder()` (`solver/catalog.py`) is called once per repair
candidate, up to `max_candidates` per iteration. **Its own docstring says it is
"SCOPED TO THE HARD_REPAIR RULES … E8's k-mer index or B1's fold evaluated here
would be pure waste"** — but it returns `self.findings(c).repairable`, and
`findings()` walks every spec. Measured per candidate:

- as implemented, all 23 specs: **62.2 ms**
- the 8 HARD_REPAIR specs it actually consumes: **5.1 ms**
- waste: **12.2x**, with `f2_near_perfect_repeats` (19.7 ms) and
  `e8_kmer_uniqueness` (17.2 ms) evaluated and discarded every time

A SOFT rule cannot return a repairable breach, so scoping the walk to
`repair_specs()` is behaviour-preserving. That is M1's file and needs an issue
under this repo's cross-lane rule, so it is reported rather than fixed here.

This lane's own knobs cannot close a 2.4x gap: `DEFAULT_SWEEP_STEPS` is already
1 (three vectors, the minimum that makes a panel), `DEFAULT_GALLERY_SIZE` below
`MIN_GALLERY` stops being a gallery, and the null is 3.1 s of the 23.5.
`max_candidates` 256 → 64 saves ~1.5 s per solve and pays for it in repair
quality. So the test stays red and says why.
