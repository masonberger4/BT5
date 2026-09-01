---
name: architect
description: Design a cross-lane change, a bt5/core amendment, or a MINOR-vs-MAJOR contract call — at opus and ultracode effort, by orchestrating independent designs. Not for implementing anything, and not for single-lane work including a new rule.
model: opus
effort: xhigh
disable-model-invocation: true
---

# /architect — opus + ultracode

**Ultracode is `xhigh` effort plus standing dynamic-workflow orchestration.** Frontmatter
gives the effort; this body gives the orchestration, and it is mandatory.

## Use this for

- A change that spans two lanes (CLAUDE.md §1 — cross-lane needs an issue first).
- An amendment to `packages/engine/src/bt5/core/**`.
- A MINOR-vs-MAJOR classification you are not certain of.

**Not** for implementing anything, and not for single-lane work — adding a rule is
`/rule-add`, not this.

## The classification that must not be got wrong

Dataclass fields and protocol methods classify in **opposite directions**, deliberately:

- BT5 **constructs** `Breach`, so a new **defaulted** field breaks nobody — MINOR.
- BT5 **implements** `FoldEngine`, so a new protocol method lands on every implementer,
  including every lane's fake — MAJOR.

The underlying question is always **who breaks?** MINOR is a new type, a new defaulted
field, a new enum member, a field that gains a default. MAJOR is a removal, a rename, a
changed annotation or default, a field that loses its default, a changed signature, a
new protocol method. A MAJOR call costs an RFC, a deprecation shim and the two-window
rule — so getting it wrong in either direction is expensive.

Hand the actual regeneration to `/contract-change`; classify here.

## Orchestrate, do not decide alone

Run a `Workflow`:

1. **Fan out** — 3 to 4 independent designs for the same problem, each from a different
   starting assumption. Include one whose only job is "what silently does not work".
2. **Adversarially verify** — pipeline each design into a verifier that checks its
   falsifiable claims against this checkout, citing `file:line`. Refuted claims die.
3. **Synthesize** — merge into one recommendation, saying why each rejected alternative
   lost.

## Return

A plan, a MINOR/MAJOR call with the "who breaks?" reasoning, the alternatives rejected
and why, and the files that would change. Append the settled decision to
`docs/decisions/` (one new file per decision) — what survives compaction is what lives on disk.

**Turn-scoped:** the override reverts on the next prompt.
