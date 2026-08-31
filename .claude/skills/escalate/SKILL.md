---
name: escalate
description: Raise this turn to opus at ultracode effort and orchestrate a workflow instead of solving alone. For a decision a normal attempt has already got wrong.
model: opus
effort: xhigh
disable-model-invocation: true
---

# /escalate — opus + ultracode

**Ultracode is `xhigh` effort plus standing dynamic-workflow orchestration.** The
frontmatter above supplies the effort half. This body supplies the other half, and it is
not optional: on this turn you **orchestrate rather than solve solo**.

## Step 1 — name the decision, in one sentence

Write it out before doing anything else. If you cannot state the specific decision being
escalated in one sentence, stop: you need `debugger` (a failure you cannot explain) or
`/architect` (a design or contract question), not more effort in this window.

## Step 2 — orchestrate

Invoking this skill **is** the opt-in to multi-agent orchestration for this turn, so
build a `Workflow` rather than reasoning alone. The shape that works:

1. **Fan out** — 3 to 4 independent attempts at the *same* question, each through a
   different lens (cost, correctness, will-it-survive-use, what-silently-fails). Give
   each the same brief and let them disagree.
2. **Adversarially verify** — pipeline each attempt straight into a verifier that tries
   to *refute* its falsifiable claims against this repo. Default to refuted when there
   is no evidence: a plausible claim with no file:line behind it is worse than none.
3. **Synthesize** — one agent merges, and a claim that was refuted stays dead.

Use `pipeline()` so verification starts as each attempt lands. Barrier only before the
synthesis, which genuinely needs everything.

## Step 3 — report the disagreement, not just the answer

Say where the lenses disagreed and why the losing view lost. An escalation that returns
a confident consensus it never tested is the failure this skill exists to prevent.

## Scope

**Turn-scoped.** The model and effort override applies to the rest of *this* turn and is
not saved; the next prompt reverts to the session defaults. For a multi-turn piece of
work — a `core/` amendment, say — re-invoke it each turn, or change the session setting.

Do not use this to grind harder on a task that is going fine.
