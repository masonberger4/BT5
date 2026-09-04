---
name: verify-provenance
description: Audit whether rules' citations actually support their thresholds — the standing review rotation.
argument-hint: "[rule-id | --sweep]"
allowed-tools: Read, Grep, Glob, Write, Agent
---

# /verify-provenance

## The only question this asks

**Does the cited source actually support the number in the code?**

`tests/data_integrity/test_rule_contract.py` already asserts, per registered spec:
citations exist, are `https` and labelled; `weight_provenance` non-empty for SOFT;
FOLKLORE defaults off; VENDOR_ASSERTED defaults on; `last_verified` parses; BAND declares
a non-inverted band; `default_weight == 0.0` for every hard rule; `steering_weight == 0.0`
for HARD_LATTICE; HARD_LATTICE declares forbidden motifs; `brief_ref` non-empty;
`param_schema` is a JSON Schema object. **Do not re-run any of that.**

Its docstring names the gap exactly:

> The contract test catches MISSING provenance; it cannot catch WRONG provenance, which
> is why `last_verified` plus a standing review rotation exists alongside it.

`docs/PLAN.md` calls that rotation standing risk #5 — *"the owner reading five rule files
a week against their citations. The one control that is a human habit rather than a CI
job, and the one most likely to lapse."* This skill is that habit.

## Method

**Enumerate the catalog. Never write a rule count into this file.**

```
Glob: packages/engine/src/bt5/rules/catalog/*.py     (ignore __init__.py)
```

That glob is the definition of "all rules". A number written here cannot be, and a stale
one fails silently in the worst direction: this skill said *"audit all 15"* for as long as
the catalog had 25 rules, so a `--sweep` would have called itself complete with ten rules
never audited, and nothing anywhere would have said so. No CI job counts the catalog, and
`rule-auditor` is told what to audit rather than discovering it — so the number here was
the only thing standing between the rotation and a silent short sweep. **If you find a
count in this file, that is the bug: replace it with the glob.**

Delegate to `rule-auditor` (opus/xhigh, project memory). With a rule id, audit that rule.
With `--sweep`, audit every id the glob returns, five at a time, skipping ids its memory
records as audited within the last quarter.

State the enumerated total, the number skipped as recently audited, and the number
actually audited. A sweep that covers less than the catalog is then visible as a number,
which is the property this section exists to guarantee.

`rule-auditor` resolves each `brief_ref` via `docs-miner` — remember `brief_ref` strings
do **not** appear literally in `brief.md`; they are section-qualified. The procedure is
in `.claude/rules/rules-catalog.md`.

## The worked example

`brief.md:141` is E4's row. Its original thresholds are **struck through** and marked
*"corrected 2026-08-28, these thresholds are below the chance floor."* A rule encoding
the struck-through number would pass all 11 CI assertions and be wrong. That is the class
of defect this skill exists to find.

## Return

Per rule: the `brief_ref` and where it resolved, whether the cited text supports the
coded threshold, and one of `SUPPORTED` / `UNSUPPORTED` / `SUPERSEDED` / `UNRESOLVED`.

Open a `--sweep` with the coverage line, so a short sweep cannot read as a complete one:

```
COVERAGE  <globbed> rules enumerated · <n> audited · <n> skipped (within the quarter)
```

Append anything settled to `docs/decisions/` as a new `YYYY-MM-DD-slug.md` file — that is
what `Write` is in the tool list for, and it is the only thing this skill writes. It never
edits a rule, a Spec or a test: finding wrong provenance and changing the number are two
different jobs, and only the first one is this skill's.
