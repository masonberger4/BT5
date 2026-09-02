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

### A panel smaller than the one requested is not a degradation

**Superseded once, then corrected — see "Corrected after review" below.** The
first version of this made a short panel emit a degradation. That was wrong: with
the shipped lattice it fired on every run and pinned `is_complete` to False by
construction.

What stands: `gallery_size` is a **ceiling**, not a promise. PLAN asks for 3-8
genuinely different candidates, so a sweep that exhausted the front at three has
answered completely and degrades nothing. `Gallery.swept` and `Gallery.distinct`
carry how many vectors solved and how many distinct designs they reached, for a
caller that wants to show "3 is all there was" rather than "we chose 3 from many".
Only a panel below `MIN_GALLERY` degrades.

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

**Sweeping densely at all.** This is the decision the measurements changed, and
the first draft of this file had it wrong.

`DEFAULT_SWEEP_STEPS` began at 3 (20 weight vectors) on the reasoning that
`build_gallery`'s shipped `steps=8` is 165 solves and 20 "looked affordable".
Measured at 500 aa on the reference backbone that was **36.2 s against a 10 s
bar** — and the profile said something more useful than "too slow":

| weight vectors | distinct designs | time |
|---|---|---|
| 3 | 3 | 5.2 s |
| 4 | 3 | 7.5 s |
| 6 | 3 | 10.4 s |
| 20 | 3 | 33.6 s |

Density bought **nothing**. Seventeen of the twenty solves rediscovered designs
the first three had already found, and the minimum pairwise codon distance was
0.399 at every setting — two and a half times G4's bar. The front these axes
reach on this protein has three points on it.

Cost is Tier B, not the DP: a solve whose CDS is already clean costs ~0.09 s, one
that needs repair ~2.5 s, because `repair()` evaluates up to `max_candidates`
candidates per iteration and each is a full catalog pass. So the only lever that
matters is the NUMBER OF SOLVES.

`DEFAULT_SWEEP_STEPS = 1` — the one-hot corners, each axis pushed alone. Those
are the extremes of the trade-off, which is what `greedy_max_min` seeds itself
from anyway; mixtures interpolate between designs the corners already reach. Das
& Dennis (1997) is the reason the SELECTION was never the sampling's job. End to
end at 500 aa: **7.37 s**, panel of 3, G4 met at 0.399. It stays a parameter of
`design()` for a protein whose front is richer.

**Sweeping a provably dead axis.** With no host codon-usage table on file,
`codon_adaptation`'s cost term is identically zero, so its corner of the simplex
solves to the unsteered design and `sweep` pays a full `optimize()` — ~2.4 s — to
rediscover it. `live_axes` drops it. That this loses no design is an empirical
claim, so it is a test and not a docstring: sweeping all four axes and sweeping
the live ones produce the same design set, one solve apart. The axis comes back
on its own when S6 ships a host table.

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


## Corrected after review

`/pre-pr`'s review pass found two blocking defects. Both were mine, both were in
the honesty vocabulary this session exists to protect, and both are recorded here
rather than quietly fixed.

### The G4 degradation reported a distance for a failure that was not about distance

`Gallery.meets_g4` is False for **two** reasons — fewer than `MIN_GALLERY`
candidates, *or* a minimum pairwise distance under 15% — and `pairwise_minimum`
returns **1.0** for a single sequence. `_panel` treated both as the distance
failure, so a one-candidate panel emitted:

> the 1-candidate panel does not meet gate G4: its minimum pairwise codon
> distance is **100.0%, below the 15%** …

A literally false sentence, in the one vocabulary this lane is for. A two-candidate
panel is the same bug wearing a plausible number ("41.8%, below the 15%"), which is
worse because nothing looks wrong. `_panel` now branches on the two conditions
separately and the short-panel sentence makes **no distance claim at all**;
`SkeletonResult.meets_g4`'s docstring carried the same conflation and is fixed.

### The short-panel degradation replaced one permanent False with another

`DEFAULT_GALLERY_SIZE` is 5 and the shipped lattice yields at most 4 weight
vectors (3 without a host usage table), so `len(picks) < k` on **every** run — the
degradation added to report a short panel fired unconditionally, and
`QcReport.is_complete` was therefore False by construction. The skeleton's
hard-wired `False` had been replaced by a different permanent `False`, which is
precisely what this increment claimed to remove.

The fix is not to lower `DEFAULT_GALLERY_SIZE`. `k` is a **ceiling**
(`greedy_max_min` returns every point when there are at most `k`), PLAN asks for
3-8 candidates, and capping at 3 would refuse a richer protein a fuller panel. A
sweep that exhausted the front at three genuinely different designs has answered
**completely**, so it no longer degrades; `Gallery.swept` and `Gallery.distinct`
carry the counts for a caller that wants to show them. Only a panel below
`MIN_GALLERY` degrades now.

**`is_complete` still does not reach True with the data that ships**, and that is
the honest answer rather than a remaining defect: every run carries at least the
biosecurity screen (S2's) and the missing codon usage table (S6's), and a
non-E. coli host also carries c1's missing CAI reference set (S3's).

The usage-table gap is worse than "S6 has not shipped the data yet", and the
re-review caught me understating it. `_host_usage` reads
`data/codon_usage/{host}.json`, but the one file that ships is named by
REFERENCE SET (`sharp_li_1987_ecoli_w.json`), not by host id — so all nine
`HostId` values miss, `kind="host_frequency"` is unreachable in the shipped
build, `SolveSpace.usage` is always empty, and `live_axes` therefore always drops
`codon_adaptation`. The "provably dead axis" is dead by a lookup-key mismatch as
much as by absent data. `c1_cai` gets this right by looking up through
`CAI_REFERENCE_SET`, which is why CAI works for E. coli while the null never
does. Every sentence the report emits about it is true; the mapping is what is
missing, and wiring it is S3's edit to `c1_cai.py` plus S6's data — not this
lane's. What changed is that
the `False` is now *derived from five named absences* instead of produced by a
sentence that could never be absent — and a test proves the mechanism by supplying
a clear verdict and watching exactly that degradation disappear.

### Also corrected

- **A rank key that compared incomparable totals.** `ScoreCard.total` renormalises
  over the objectives *that candidate* could evaluate, so a candidate measured on
  strictly less could outrank one measured on more — and the winner is what gets
  exported to GenBank and put on the order file. `comparable_totals` now computes
  the sort key over the objectives available to **every** candidate. Unreachable on
  today's catalog; the ranking is what the user acts on, and "unreachable today" is
  not a property a report can carry.
- **NaN in the null's own values was unguarded.** `unavailability` guarded the
  candidate's raw score only. A rule that computes on the anchor and fails on a
  variant would put NaNs in the null, which `percentile_of` counts as neither
  better nor a tie — silently deflating the percentile while `null_mean` renders as
  `nan`. The null is now discarded with a stated reason.
- **A comment cited a test that did not exist.** `build_nulls` hard-codes
  `windowed_fold_only=True` and claimed a test held it.
  `test_the_null_never_folds_whole_transcripts` now does, by running every scored
  objective against a fold engine whose `mfe` raises.
- **The degradation guard had become a catch-all.** Replacing set equality with
  prefix matching was right; using `"the "` as one of the prefixes was not — four
  characters matching any future sentence starting "the ", inside the test named
  for catching unremarked degradations. Anchored regexes now, each pinned to the
  sentence its source actually emits.
- **The `is_complete` test was a tautology** restating the property's own body. It
  now tests that filling an absence removes its degradation.
- **A dropped assertion** at the origin-spanning site, and a **dead tie-break**
  whose comment claimed work `sorted`'s stability already did.

### Left as a follow-up, deliberately

`RuleSet.gated_out` objectives are named nowhere on the report. On the reference
fixture `b1_five_prime` — the highest-weight SOFT objective — is gated out because
the slot is HEK293 rather than E. coli, which means "does not apply" rather than
"could not be evaluated", so it does not belong in `unavailable` or in
`degradations`. But `solver/catalog.py` records `gated_out` precisely so a report
can say which, and this lane's whole argument is that an absent objective is
indistinguishable from one that was never configured. Surfacing it needs a new
`QcReport` field and a render section; that is a report-contract change S5's CLI
renders, and it belongs in its own PR rather than bolted onto this one.
