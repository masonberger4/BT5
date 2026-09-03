---
name: rule-add
description: Add a rule to bt5/rules/catalog with its paired test, provenance and a resolvable brief_ref.
allowed-tools: Read, Edit, Write, Grep, Glob, Bash, Agent
---

# /rule-add

## Non-negotiables

1. **You are not done without `packages/engine/tests/rules/test_<id>.py`.** Exactly one
   catalog rule ships without a paired test — `d1_restriction_sites.py` — and its own
   line 3 says *"REFERENCE RULE. Copy this file's shape when adding a pattern rule."*
   **Do not copy that file's test story.** Copy `e1_homopolymers.py` together with
   `packages/engine/tests/rules/test_e1_homopolymers.py`.
2. **`brief_ref` must actually resolve.** The contract test only asserts it is non-empty.
   See the resolution procedure below — a `brief_ref` that resolves to nothing means the
   rule's evidence is unauditable.
3. **`SINGLE_PASS` is a downgrade, not the default.** `RepairPolicy.FIXED_POINT` is
   `solver/repair.py:417`'s own default and the `repair` ClassVar carries none of its
   own, so a spec reads `SINGLE_PASS` only because an author typed it. The docstring
   owes the reader why this repair cannot create new instances of what it removes. If
   it can — splice donors are the canonical case — `FIXED_POINT` is **mandatory**
   (CLAUDE.md §3.6), not a judgement call.
4. **There is no registration step.** Rules autodiscover via the `@register` decorator.
   `core/registry.py`: *"Adding a rule edits ZERO shared files, not even an `__init__.py`."*
   Do not add your rule to a list.
5. **`default_weight` must be `0.0` for every hard rule.** The weighted sum only ever
   sees SOFT. Use `steering_weight` if the DP needs nudging (§3.5).

## Which file to copy

- **Scored rule** → `e2_gc_band.py` (line 3: *"REFERENCE RULE. Copy this file's shape
  when adding a scored rule."*)
- **Pattern rule** → `d1_restriction_sites.py` for the *rule* shape — but take the test
  shape from `test_e1_homopolymers.py`, since d1 has none.

## Resolving `brief_ref`

`brief_ref` values look like `2.E4`. **They do not appear literally in `brief.md`** —
`grep -F "2.E4"` returns zero hits. Resolve in two steps:

1. `grep -n '^### 2\.E' docs/research/brief.md` for the section.
2. Find the row beneath it. Two shapes only: a table row `| E4 | ... |`, or a bold
   run-in `**D1 Restriction / Type IIS ...:**`.

Delegate the lookup to `docs-miner` rather than reading the 63 KB file inline. **Check
for supersession**: `brief.md:141` (E4) is struck through and corrected 2026-08-28.
Never encode a struck-through threshold.

## Required ClassVars

`id` (must equal the filename stem) · `version` · `title` · `enforcement` · `evidence` ·
`direction` · `unit` · `citations` (≥1, `https`, labelled) · `last_verified` (ISO date) ·
`weight_provenance` (**required non-empty for SOFT** — prose saying where the default
weight came from) · `default_enabled` (**FOLKLORE ships `False`**; VENDOR_ASSERTED ships
`True`) · `default_weight` · `steering_weight` · `band` (required for `Direction.BAND`,
and `lo < hi`) · `localization` · `repair` · `cost_class` · `conflicts_with` ·
`param_schema` (a JSON Schema object) · `brief_ref` · `engine_calibration`.

Also: never scan the reverse strand yourself — list forward motifs in
`LatticeTerms.forbidden` and let the solver close the set under reverse complement.
A directional scored model must read `slot.strand_of_interest` (`core/spec.py:274`),
never a hard-coded strand 1 (§3.4).

## Verify

```bash
.venv/bin/pytest tests/data_integrity -q
.venv/bin/pytest packages/engine/tests/rules/test_<id>.py -q
```

`tests/data_integrity/test_rule_contract.py` parametrises 11 assertions over your new
spec automatically. Then run `/pre-pr`.
