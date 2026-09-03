---
name: skill-frontmatter-governance
description: disable-model-invocation on /pre-pr was decided open on 2026-09-03 — the flag is gone deliberately, so check the decision record before flagging it; the general lesson is that .claude/skills/** frontmatter governs agent autonomy and is not in CLAUDE.md §2
metadata:
  type: feedback
---

**Settled, do not re-flag.** `/pre-pr` no longer carries `disable-model-invocation: true`.
The owner removed it on 2026-09-03 as a deliberate decision, recorded in
`docs/decisions/2026-09-03-agents-may-attest.md`, which supersedes the open item at
`docs/decisions/2026-09-01-pre-pr-as-ci-gates.md:56`. An agent may run the chain and post
the `/pre-pr <sha>` attestation, under step 10's existing condition: **not** if a gate or
a review came back blocking. `/pre-pr-bypass` stays OWNER-only and stays a human act.

What the attestation proves changed with it: no longer "a human vouched for this commit",
now "the agent's chain ran against this commit and did not come back blocking". That is
weaker on purpose — an attestation only a human can post is not automatable, and automated
CI checks were the goal.

**Why this file used to say the opposite, and what was wrong with it.** It carried a
standing rule to flag removal of this flag as blocking "regardless of what else the PR is
nominally about". Two of its supports did not hold:

- It cited **CLAUDE.md §3** for the `/pre-pr` vs `/pre-pr-bypass` split. §3 is the
  non-negotiable *correctness* rules — genetic code tables, stop codons, `Construct`,
  reverse-complement closure, hard-constraint enforcement, splice fixed-point, RNG
  seeding. It says nothing about either command. The split is argued at
  `.github/workflows/pre-pr-attest.yml:29-31`.
- It read the workflow header's "It does NOT prove the review was done honestly — Claude
  posts the comment" as conditioned on a human having typed the command first. The header
  does not say that; the qualifier was inference presented as quotation.

The finding itself was worth raising — it *was* an undocumented governance change when
first reviewed. What was wrong was asserting a repo rule that did not exist, and making
the verdict unconditional so that it would keep firing after the owner decided.

**How to apply.** `.claude/skills/**` and `.claude/agents/**` frontmatter is not in
CLAUDE.md §2's protected-paths list but does govern agent autonomy, so a change there is
worth *reading* rather than skipping. Before assigning severity:

1. Grep `docs/decisions/` for the flag or skill name. A recorded decision means the
   question is answered — note it and move on. Do not re-litigate an owner decision.
2. If there is no decision record, the finding is that it is **undocumented**, not that it
   is forbidden. Say what the change lets an agent do that it could not before, and ask
   for it to be recorded — that is the actionable gap.
3. Quote repo documents you cite, and check the section number. A confident wrong citation
   costs more than no citation, because it reads as authority.

See also [[ci-workflow-review]].
