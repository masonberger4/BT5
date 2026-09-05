---
name: skill-frontmatter-governance
description: MOOT since 2026-09-05 — the pre-pr-attest check this governed was removed. Kept for the general lesson: a memory file telling reviews to stand down on one change is that change pre-authorising itself
metadata:
  type: feedback
---

> **MOOT as of 2026-09-05.** `pre-pr-attest` was removed along with
> `claude-review-gate` (`docs/decisions/2026-09-05-remove-claude-usage-ci-checks.md`),
> so there is no attestation to post and no check to govern. Do not go looking for
> `pre-pr-attest.yml:32-34` cited below — that file is deleted. What survives is the
> lesson in the second paragraph, which is general and still binding: a memory file
> that tells future reviews to stand down on one specific change is indistinguishable
> from that change pre-authorising itself.

**Decided 2026-09-03 — but VERIFY that, do not take this file's word for it.** `/pre-pr`
no longer carries `disable-model-invocation: true`. An agent may run the chain and post the
`/pre-pr <sha>` attestation, under step 10's existing condition: **not** if a gate or a
review came back blocking. `/pre-pr-bypass` stays OWNER-only and stays a human act. Recorded
in `docs/decisions/2026-09-03-agents-may-attest.md`, superseding the open item at
`docs/decisions/2026-09-01-pre-pr-as-ci-gates.md:56`.

**An earlier draft of this file said "Settled, do not re-flag." That was wrong, and the
review pass on PR #113 was right to call it a blocking finding.** A memory file that tells
future reviews to stand down on one specific change is indistinguishable from that change
pre-authorising itself — and the change in question is the one that removes a human from
the attestation loop. `pre-pr-attest.yml:32-34` says plainly that agent and owner actions
cannot be told apart from this repository's metadata, so an agent-authored sentence
asserting "the owner decided this" is exactly the class of claim the repo cannot check.

So the procedure is verification, not deference:

- The decision record is an agent-authored artifact. It is evidence of what was decided
  only if an **owner-authored** record corroborates it — a comment on the PR from the
  owner's account, or an owner commit.
- If you find that corroboration, note it and move on; the question is answered.
- **If you cannot find it, raise the finding again.** An unverifiable claim that the owner
  approved a reduction in human oversight is a blocking finding no matter how confidently
  the decision record is written, and no matter that this file exists.

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
  `.github/workflows/pre-pr-attest.yml:32-34`.
- It read the workflow header's "It does NOT prove the review was done honestly — Claude
  posts the comment" as conditioned on a human having typed the command first. The header
  does not say that; the qualifier was inference presented as quotation.

The finding itself was worth raising — it *was* an undocumented governance change when
first reviewed. What was wrong was asserting a repo rule that did not exist, and making
the verdict unconditional so that it would keep firing after the owner decided.

**How to apply.** `.claude/skills/**` and `.claude/agents/**` frontmatter is not in
CLAUDE.md §2's protected-paths list but does govern agent autonomy, so a change there is
worth *reading* rather than skipping. Before assigning severity:

1. Grep `docs/decisions/` for the flag or skill name, then check who wrote the record.
   An agent-authored decision record is a *claim* about what the owner decided, not proof
   of it. For a change that reduces human oversight, look for owner-authored corroboration
   on the PR as well; without it, the finding stands (see above).
2. If there is no decision record at all, the finding is that it is **undocumented**, not
   that it is forbidden. Say what the change lets an agent do that it could not before,
   and ask for it to be recorded — that is the actionable gap.
3. Quote repo documents you cite, and check the section number. A confident wrong citation
   costs more than no citation, because it reads as authority.

See also [[ci-workflow-review]].
