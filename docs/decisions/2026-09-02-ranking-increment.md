# 2026-09-02 — The ranking increment: percentiles, a panel, a baseline, an order file

**Lane:** M11 design + M3 score (S1). **Branch:** `claude/s1-ranking-increment`.
**Scientific impact:** non-"none" — the app now emits a ranked panel with an order
file where it emitted one unranked candidate. Under `CLAUDE.md` §7b that goes to
the owner, not to a self-merge.

PR #71 landed `bt5/design/` and its docstring listed what it refused to do: one
candidate, no gallery, every objective `unavailable`, `native_baseline` None, no
percentiles, no order CSV. `bt5/score/` already exported `build_gallery`,
`null_distribution`, `percentile_of`, `normalise`, `order_entries` and
`write_csv`, and nothing on the design path called any of them. This is that
wiring. What follows is what I decided, and what I rejected.

## Decided

### The null is sized by `cost_class`, and each objective carries its own `n`

`Spec.cost_class` is documented as the thing that "drives null sampling", and
`ObjectiveScore.null_n` is already per-objective, so this needed no new contract
field. `NULL_N_BY_COST` gives `cheap` the shipped `DEFAULT_NULL_N` (200) and
`moderate` 40; an unknown cost class gets none, and the objective is reported
unavailable rather than silently costing the budget.

Every objective's null is drawn from ONE seeded stream, which makes a moderate
objective's 40 variants a strict **prefix** of a cheap objective's 200 — the same
variants, fewer of them. That is a property of `null_distribution` re-seeding
`default_rng(seed)` per call with the same anchor CDS and synonym map, and it is
what lets a single memoised assembler serve every objective: `max(n)` assemblies
instead of `n × objectives`.

### One null, shared by every candidate

The null is a distribution over random synonymous variants of the protein in this
construct. It does not depend on which candidate is being ranked — anchoring it
on candidate A and ranking candidate B against it is not an approximation, it is
the same distribution. Building it once is also what makes the panel's
percentiles comparable **to each other**, which is the gallery's whole job.

### GC lean is a sweep axis, because CAI-versus-repeats alone is not a sweep

`build_gallery` needs a `solve(weights)` whose output moves as the weights move.
The solver ships `cai_scorer` and `repeat_breaking_scorer`, and both are dominated
by the same relative-adaptiveness table: at every mixture the argmin is very
nearly the same codon, so the panel clusters and G4 fails for a reason that has
nothing to do with the constraint set.

`score/steering.py` adds `gc_lean_at` / `gc_lean_gc` as two ends of one axis.
Third positions are where synonymous choice lives, so leaning GC the two ways
moves the CODON choice — the space `codon_distance` measures and G4 gates on —
while `solve_with_gc_steering` keeps both ends inside E2's band. The lean chooses
*where in the band* a design sits; the band is still enforced by Tier B plus the
independent validator.

`REPEAT_STEERING_PENALTY` is 4.0, not `repeat_breaking_scorer`'s 100.0. At 100
the repeat term swamps every other term at any non-zero weight, most of the
simplex collapses onto one design, and the sweep stops sweeping. 4.0 still
outranks the whole [0, 1] range of a codon-preference difference. This is a
steering weight in `CLAUDE.md` §3.5's sense and enforces nothing: E5 and F1 are
HARD_REPAIR regardless of what weight the sweep picked.

### NaN never becomes a percentile

`c1_cai` and `b1_five_prime` return NaN plus a breach carrying the reason when
they cannot compute their objective — `c1_cai._unavailable` argues the case, and
it is right: 0.0 would read as a catastrophically rare-codon sequence and the
band's midpoint as a design exactly on target.

Which means NaN must never reach `percentile_of`. Every comparison against NaN is
False, so it would count `better = 0`, `ties = 0` and emerge as a **confident
0.0** — "worse than every one of 200 random variants" — about a quantity nobody
measured. Nothing raises; the report just becomes false. `ranking.unavailability`
is the guard, applied twice: once before a null is built (an objective the rule
cannot compute on the real construct cannot be computed on 200 variants of it
either) and once per candidate, because a rule can compute its objective on one
design and not on another.

### `native_baseline` is populated only from a real wild-type CDS

`design(native_cds=...)` is a new defaulted keyword-only parameter. There is
deliberately no way to ask BT5 to invent one: "do not optimize" compared against
a sequence BT5 itself designed is not a comparison. A supplied native CDS goes
through the same `verify_construct` as every design, and a native the validator
refuses is reported as refused with the invariant named — never dropped, never
back-translated into something that would pass. It is on the order file, because
"do not optimize" is only a real option if the user can order the tube.

### `is_complete` is derived, so a `True` means something

The skeleton emitted three degradations unconditionally ("ranking not computed",
"single candidate only", "screening did not run"), which pinned `is_complete` to
False by construction. Each is now conditional on the thing genuinely being
absent. `design()` also gained a `screen: BiosecurityVerdict` parameter defaulting
to `not_run`, so this lane RENDERS whatever verdict it is handed and never
upgrades one — making the screen real is S2's, per the buildout's inter-session
contract.

### A G4 shortfall is a finding

If the panel's minimum pairwise codon distance lands under 15%, the report says
so and names the distance actually reached. G4's failure invalidates a **product**
decision, not a technical one: if the sweep cannot produce genuinely distinct
designs then the gallery is not a feature and a UI built on it is a lie. The
threshold is pinned by a test so a later change to it has to come through review.

## Rejected

**Lowering `G4_MIN_PAIRWISE_DISTANCE` to whatever the sweep reaches.** The
obvious way to make the panel "pass". `score/gallery.py`'s own docstring says a
G4 failure means re-planning around ε-constraint enumeration *before* a UI is
built on it. Reported, never lowered.

**Sweeping at `build_gallery`'s shipped `steps=8`.** Four axes at 8 steps is 165
lattice points, so 165 full `optimize()` calls — Tier A, Tier B and a validator
run each. That is not a 10 s design. `DEFAULT_SWEEP_STEPS = 3` (20 vectors) is
what fits, and the sampling was never what makes a gallery diverse — Das & Dennis
(1997) is the reason the SELECTION is greedy max-min on sequence distance and not
the weights.

**Scoring the null on bare CDSs.** Much cheaper, and it would make "94th
percentile" a measurement against a distribution that never contained a backbone.
`null_distribution` takes a `build` callable precisely so this cannot be done by
accident, and threading the real assembler through is the whole point.

**Evaluating the catalog twice per candidate** — once for HARD_CHECK findings and
once for raw scores. `RuleSet.findings()` returns both from one pass;
`score_candidate` therefore takes evaluations rather than a construct.

**Borrowing the E. coli codon table for a mammalian host's null.** No host usage
table ships for `hek293` (the build carries only Sharp & Li's E. coli w-index), so
the null degrades to uniform-synonymous sampling and SAYS so —
`NullDistribution.kind` carries which, and it reaches the report through
`ObjectiveScore.null_kind`. Reusing E. coli weights would produce a
plausible-looking number measuring nothing, which is the same failure
`c1_cai` already refuses by returning NaN. Wiring a real reference set is S6's
data and S3's edit to `c1_cai.py`, per the buildout README.

**Reaching into `CAI_REFERENCE_SET` in `rules/catalog/c1_cai.py`** to make the
null's weights match c1's. Out of lane, and the buildout README names that edit
as S3's.

**Marking the G7 timing test `slow`.** `.claude/rules/tests.md` records that
`-m "not slow"` currently deselects nothing because the marker is applied zero
times. Marking this test would make it the first use of a marker whose only
effect is that `gates.sh` and CI's engine job never run it. A timing gate nothing
executes is worse than no timing gate, so it runs in the suite and a failure is a
finding about the budget.

**Creating `benchmarks/`.** `CLAUDE.md` §2 protects two files in a directory that
does not exist; creating it is an owner decision under
`approved:algorithm-change`. A plain `pytest` assertion carries G7 instead.

**Adding a field to `core/`.** Nothing here needed one: `ObjectiveScore` already
carries `null_n`, `null_kind`, `windowed_fold_only` and `unavailable_reason`;
`Candidate` already carries `codon_distance_to`; `DesignResult` already carries
`native_baseline`. `core/` is a global mutex in this buildout and another session
may be waiting on it.

## Contracts honoured

- `design()`'s signature is frozen. Every new parameter is keyword-only with a
  default, so every call that worked before still works. `SkeletonResult` gained
  `gallery`, `orders`, `order_csv` and `nulls`; nothing was removed or renamed.
  S5's CLI is built against this type.
- `BiosecurityVerdict`'s shape is untouched. This lane renders the status it is
  given and never prints "clear" for `not_run`.
- No file under `rules/`, `solver/`, `vector/`, `cassette/`, `codon/`,
  `structure/`, `core/`, `.github/`, `data/`, `tests/contract/`,
  `tests/invariants/`, `tests/data_integrity/` or `pyproject.toml` was modified.
