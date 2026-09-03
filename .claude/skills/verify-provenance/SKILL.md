---
name: verify-provenance
description: Audit whether rules' citations actually support their thresholds — the standing review rotation.
argument-hint: "[rule-id | --sweep]"
allowed-tools: Read, Grep, Agent
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

Delegate to `rule-auditor` (opus/xhigh, project memory). With a rule id, audit that rule.
With `--sweep`, audit all 25, five at a time, skipping ids its memory records as audited
within the last quarter.

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
Append anything settled to `docs/decisions/` as a new `YYYY-MM-DD-slug.md` file.
