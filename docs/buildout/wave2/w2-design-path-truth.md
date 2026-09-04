# W2 — the design path says what it means

*Copy this whole file into a fresh Claude Code session on the BT5 repo.*

**Launch:** permission mode `acceptEdits` · model opus · effort **xhigh** (override) · runs
unattended.
**Do not run this in plan mode.** **You post a design note as the first comment on your
draft PR before you implement** — that is your plan-mode substitute, and it is load-bearing
here because you decide the semantics of `is_complete`.

**Cut your branch only after W1 (`claude/w1-host-reference-sets`) has MERGED.** This is not
a nicety. See "Why you go last".

---

You are finishing the job PR #89 started. `bt5/score/` still exports capability that
`bt5/design/` never calls, and the report still tells every user something that is not true.

## Read this first

Run `/bootstrap`. Then `CLAUDE.md`, then `docs/buildout/wave2/README.md`, then
`docs/decisions/2026-09-02-ranking-increment.md` — that is PR #89's record, it names three
of the five things you are building as deliberate follow-ups, and its "Corrected after
review" section is the standard your own honesty vocabulary will be held to.

Your branch: **`claude/w2-design-path-truth`**, cut from a **freshly fetched** main *after
W1 merges*:

```bash
git remote prune origin
git fetch -q origin main && git checkout -B claude/w2-design-path-truth origin/main
```

## Why you go last

You and W1 **each change the winner**. W1 flips the null from uniform-synonymous to
host-frequency, moving every percentile; your preset weights change the rank key directly.
Landed together, a moved winner cannot be attributed to either — and that attribution is
exactly what the PR template's "scientific impact" section and §7b's owner-merge gate exist
to read.

Branching after W1 lands costs one PR cycle and buys a diff that shows only your delta
against an already-rebaselined `main`.

W3, W4 and W5 are adding five rules in parallel. Their rules land in the catalog you score
against, so rebase rather than assume, and expect `constraint_set_hash` to move.

## What to build

### 1. Presets drive the weights

`design/ranking.py:317` builds the weighted sum from `spec.default_weight`:

```python
weights = {spec.id: float(spec.default_weight) for spec in specs}
```

Meanwhile `score/presets.py` — 19 KB, three presets, `WeightEntry` provenance,
`resolve() -> ResolvedPreset` — is called by **nobody**, and `DesignContext.weights`, a
field that already exists in frozen `core/` with a default factory, is **read by nobody**.
`ResolvedPreset.weights` is even documented as "spec_id → weight, ready for
`DesignContext.weights`".

Wire it: `presets.get(preset_id)` → `resolve()` → `DesignContext.weights` → read by
`weighted_total`, `comparable_totals` and `score_candidate` as a defaulted keyword.

#### The decision this forces, and it needs the owner's sign-off

Preset weights are deliberately sparse, and **the sparseness is the scientific claim**.
`LENTIVIRAL` and `AAV` weight exactly three objectives (`2.F2` 1.0, `2.C1` 0.2, `2.C3` 0.3);
`BACTERIAL` adds `2.B1` at 1.0. Each preset's `rationale` argues its omissions explicitly —
LENTIVIRAL's says the internal-polyA hazard

> carries no weight here precisely because it is not a trade-off: here it is hard, removed
> by repair and proven by the independent validator, which refuses to emit rather than
> pricing a truncated genome against expression.

Today all 13 SOFT rules enter the sum at `default_weight`, which contradicts every one of
those rationales. So this is **not** "presets narrow the sum" — the sum currently in use is
the *wrong* one, and the argued weight set has been sitting unused since M3 landed.
Switching does narrow it, from 13 objectives to 3 or 4.

**The policy — three cases, none of them silent:**

- **`preset_id is None`:** keep `default_weight`, and emit a note naming the preset that
  *would* apply. **Do not auto-apply `preset_for(modality)`.** That would make the ranking
  depend on a value the caller never supplied and cannot see at their own call site, and it
  would silently re-rank every existing `design(modality=LENTIVIRAL)` call. Opt-in is the
  only version where `WeightEntry.note` provenance answers a question the user actually
  asked.
- **The preset names the objective:** its weight wins, with `WeightEntry.note` as the
  provenance the report shows.
- **The preset omits a scored objective:** report it with its percentile and **weight 0.0,
  flagged on the scorecard as outside this preset's claim.** Not dropped, not silently
  defaulted. An objective the user can see that did not move the rank is honest; one that
  vanished is not.

`ResolvedPreset.unimplemented` carries the mirror case (a preset weighting a `brief_ref` no
rule implements) and should reach the rendered report too. An unknown `preset_id` must
**raise**, not become a label.

#### While you are in there: `strain_protocol`

`QcReport.strain_protocol` exists (`report.py:167`), `build_report` accepts it
(`report.py:187`), and `render()` prints it (`report.py:298-301`). `runner.py:632` simply
omits the argument. So this — cited, evidence-backed, already written —

> Propagate in a recA- strain (Stbl3 or NEB Stable) at 30 C. This protects the LTRs… Stbl2
> to Stbl3 alone rescued an HIV vector lost entirely in 0.5 L Stbl2 cultures.

reaches no user. **One keyword argument.**

#### Blocked-by: issue #82

Wiring `resolve()` makes **#82** reachable. `_unscored_enforcement` (`presets.py:168-201`)
has two fallbacks that both resolve to "treat this spec as scored": a spec needing
constructor arguments, and one not exposing `enforcement_for`. Both fall back to the
class-level `enforcement` ClassVar — which is exactly the read that caused #72. It is dead
today by coincidence of the current catalog, not by construction, and **three sessions are
adding rules to that catalog right now.**

Close it with #82's own suggested fix: a `tests/data_integrity` assertion that every
registered spec is no-arg constructible and exposes `gate`/`enforcement_for`. The precedent
is `core/registry.check_engine_calibration` — a rule whose calibration differs *raises*
rather than being skipped, because a skipped rule is a missing constraint nobody sees.

Add a second `data_integrity` assertion: **all three shipped presets `resolve()` without
raising against the live catalog.** That is what stops W3/W4/W5 from landing a rule that
breaks a preset at import — and it is the check the wave README promises them.

### 2. The conflict panel

`runner.py:630` passes `conflicts=()`. Everything else already exists:
`detect_conflicts(evaluations, specs, *, length, circular)` is a pure function
(`score/conflicts.py:144`), the runner already holds evaluations per candidate,
`QcReport.conflicts` is a field, and `render()` prints them (`report.py:264-271`).

This is nearly free and produces real output immediately, because rules genuinely declare
`conflicts_with`: b8↔d1 (NcoI `CCATGG` inside Kozak `GCCACCATGG`), b8↔b9, d8↔e2, e2↔f5,
e3↔{e2,f5}, d2↔d1, plus c1's list. Both halves fire — positional conflicts discovered by
wrap-aware overlap, structural ones declared.

**Ship `hard_versus_soft` beside each conflict. Leave `relaxations = ()` and say so
explicitly in the render.**

`core/result.Relaxation` — "a specific, costed way out of a conflict" — is constructed
nowhere in the engine, and it is **deliberately out of scope**:

- `predicted_cost` can only be filled honestly by re-solving at the relaxed threshold and
  diffing percentiles. A clean solve is ~0.09 s and one needing repair ~2.5 s, against a
  10 s G7 budget already at 7.37 s and about to grow. One relaxation per conflict is a
  second sweep. It does not fit.
- A `predicted_cost` that was not measured **is a predicted number**, in the one field named
  `predicted_*`. That is the vocabulary `CLAUDE.md` bans outright.
- Doing it properly *is* issue #45 — `docs/PLAN.md:505` says a G4 failure makes ε-constraint
  enumeration primary, and #45 is titled "options instead of a conflict at the end".

So: **file the issue with your G7 numbers in it**, so the blocker is recorded as a
measurement rather than a preference. And do not render an empty `relaxations` tuple without
comment — that implies there is no way out. `hard_versus_soft` answers the question the user
actually acts on ("can I slide out of this?") at zero solver cost, and its own docstring
says presenting a hard-over-soft conflict beside genuine trade-offs "invites them to try".

### 3. Retire the biosecurity degradation

`runner.py:783` emits a degradation whenever `screen.status != "clear"`, and
`UNSCREENED = BiosecurityVerdict("not_run", …)` is the default. The owner killed the screen
(PR #87 — see `docs/decisions/2026-09-03-biosecurity-screen-dropped.md`, which W0 landed),
so **nothing will ever return `"clear"`** and `QcReport.is_complete` is pinned `False` by
construction. That is the exact defect PR #89's decision record says it existed to remove.

**Keep the type. Keep the parameter.** `BiosecurityVerdict` is frozen (`core/context.py:97`)
and `design()`'s signature is frozen by the wave's inter-session contract. Removing either
is **MAJOR** under §2a — RFC, deprecation shim, two-window rule, and `pytest tests/contract`
passing without regenerating. It would also destroy the only type in which a `block` is
expressible, and a caller who *has* screened elsewhere can still hand one in.

What changes is only the `not_run` case:

- `_degradations` degrades on **`flag` and `block` only**.
- Add `QcReport.screening: str = ""` rendered in its own **always-present** section: BT5
  does not screen the input protein; DNA-synthesis vendors screen orders before synthesis
  and that is the gate; the verdict carried is `not_run`.
- `block` continues to refuse, via `BiosecurityVerdict.may_proceed`.

**Rejected — delete the degradation outright.** `is_complete` would go True while the report
said *nothing at all* about screening. That is the silent-omission failure the whole
degradation vocabulary exists to prevent (`ObjectiveScore.unavailable_reason`'s docstring
makes the argument in full), and it would render a `flag` or `block` a caller can still
supply invisible.

**Cost: zero `core/` change.** `QcReport` lives in `score/report.py`, which is **not** in
`CORE_MODULES` — `tests/contract/surface.py:47-53` lists exactly
`bt5.core.{types,context,spec,result,services,registry}`. So neither this field nor
`gated_out` is a contract change, and **`regenerate.py` never runs.** Verify that claim
yourself before relying on it.

### 4. Surface `gated_out`

`solver/catalog.py:130,404-418` records the spec ids that do not **apply** in this slot. On
the reference fixture `b1_five_prime` — the highest-weight SOFT objective in any BT5 preset
— is gated out because the slot is HEK293 rather than E. coli. That is "does not apply here",
which is neither `unavailable` ("could not evaluate") nor a degradation ("this run is
incomplete"), and the report says **nothing at all** today. PR #89's record names surfacing
it as a deliberate follow-up.

Add `QcReport.gated_out: tuple[str, ...] = ()` and a render section.

**The pairing that makes this load-bearing rather than decoration.** `gated_out` must **not**
block `is_complete` — a rule that does not apply is an answer, not an absence. **But** the
render must print the gated-out section *above* the completeness line, and a test must
assert that a report with non-empty `gated_out` renders every name.

Without that pairing you have shipped the third permanent lie in a row: the skeleton's
hard-wired `False`, the short-panel degradation's accidental one, and now an
`is_complete: True` printed while the highest-weight objective is absent from every field a
reader looks at. **This is the decision your design note exists to surface.**

### 5. Gate G2's second step, and the deferred findings

- **#74 step 2 / PLAN gate G2.** `tests/data/backbones/real_lenti_pFTMGW_EF177827.gb` — real
  NCBI EF177827.1, 8928 bp, with real LTRs, a real WPRE and `misc_feature` soup — was
  committed by S6 and **no test references it.** Add the end-to-end `design()` run on the
  real map, asserting the annotated GenBank round-trips with features preserved. It lives
  here rather than in a rules session because it is a `tests/design/` file and you already
  own that directory. Only then may G2 be recorded as covered — a hand-built fixture is a
  fixture the gate was written against.
- **#99** — panel `Candidate.construct` values are unannotated, so `design_hash` provenance
  no longer rides on candidates, while `core/result.py`'s `Candidate` docstring implies it
  does. **Pick a reading.** If "constructs stay unannotated now that annotation is a
  per-export step" is intended, the docstring fix is a `core/` change → `/contract-change`,
  and **classification precedes regeneration** because `regenerate.py` writes first and
  returns 0 on every path.
- **#100** — `build_nulls` decides objective availability from `picks[0]` only. Unreachable
  today and it fails in the safe direction (an honest `unavailable` where a real score was
  available). Fix it or close it as won't-fix with the reason — but do not leave it open and
  unaddressed by a PR that touches `build_nulls`.

## Scientific impact

**Non-"none".** Preset weights change the rank key, therefore the winner, therefore the
exported GenBank and the order file. You also change what the app says it did not do. Both
go to the owner under §7b. **You do not self-merge.**

## The design note, posted first

First comment on your own draft PR, before implementing. It must decide, in writing:

1. **`is_complete`'s semantics** — specifically the `gated_out` pairing above. This is the
   one the owner would most want to redirect while it is still cheap.
2. **The three-case weight policy**, and what the scorecard shows for an objective the
   preset omits.
3. **What the screening section says**, verbatim — it is a sentence every user will read on
   every run, forever.
4. Your **cumulative G7 number** (see below).

## G7, and you carry the total

You land last, after W1 re-armed a sweep axis (~2.4 s) and W3–W5 added five rules (each
costing a full catalog pass per repair candidate, and repair is the dominant term). The bar
is `G7_SECONDS = 10.0` at 500 aa, and `main` was at 7.37 s before any of it.

**Never raise the bar. Never mark the test `slow`** — `-m "not slow"` currently deselects
nothing, so the marker's only effect would be that `gates.sh` and CI never run it.
`DEFAULT_SWEEP_STEPS` stays 1. If the budget is exceeded, report the measured number as a
finding: `docs/PLAN.md:508` gives G7's fail consequence verbatim as *"re-allocate budget
before rules multiply cost"*, and this wave is that sentence coming true.

## Acceptance

- **One end-to-end run where `QcReport.is_complete is True`** — never true in this repo's
  history. Proved by a test that *fills* each named absence and watches its degradation
  disappear, **not** by asserting the property's own body. PR #89 shipped exactly that
  tautology and its review caught it; do not repeat it.
- The rendered report shows: preset id, its weights and their provenance; a conflict panel
  with real entries and `hard_versus_soft` per entry, or a stated absence; a `gated_out`
  section above the completeness line; the strain protocol; and **no** biosecurity
  degradation.
- `render()` never emits "clear" for a status that is not `clear`. Pin it with a test.
- An unknown `preset_id` raises rather than becoming a label.
- The degradation-set equality test still pins every degradation by anchored regex — PR
  #89's review found `"the "` used as a prefix match inside the very test named for catching
  unremarked degradations. Do not reintroduce a catch-all.
- `pytest tests/contract` passes **without** `regenerate.py` having been run, which is your
  evidence that `QcReport` is genuinely outside `CORE_MODULES`.
- A decision record under `docs/decisions/`.

## Files

- `packages/engine/src/bt5/design/**`
- `packages/engine/src/bt5/score/**`
- `packages/engine/tests/design/**`, `packages/engine/tests/score/**`
- `tests/data_integrity/**` (carries `approved:oracle-change`)
- `docs/decisions/2026-09-03-*.md`

**Do not** touch `codon/tables.py` or `rules/catalog/c1_cai.py` — W1 owns those and lands
before you. **Do not** touch any rule file: W3, W4 and W5 are live.

## Delegation

- **The G2 real-map fixture's shape** → `Explore` on `tests/data/backbones/`.
- **PR #89's follow-ups and its corrected-after-review list** → `docs-miner` on
  `docs/decisions/2026-09-02-ranking-increment.md`.
- **Gates** → `gate-runner`.
- **A failure a first-pass fix missed** → `debugger`, not another attempt.
- **Never bare-Read `design/runner.py`** (33 KB) or `score/presets.py` (19 KB) — use
  `offset`/`limit`, or delegate.

## Done means

- `bash scripts/gates.sh` → `ALL GATES PASSED`, or G7 fails and you have reported the
  cumulative number as a finding rather than tuned anything.
- Design note posted **before** the implementation commits.
- Issues #82, #99, #100 and #74 each either closed by this PR or restated with the answer
  chosen and the reason.
- The `Relaxation` issue filed, with the G7 measurement in it.
- Scientific impact section says **non-"none"** and names both the ranking change and the
  change to what the app says it did not do.
- Draft PR, `/pre-pr` from the operator, **owner merge — you do not self-merge**.
