## 2026-09-04 — `security-reviewer` is deleted; `boundary-reviewer` takes its checks

**Decided:** `.claude/agents/security-reviewer.md` is removed. `.claude/agents/boundary-reviewer.md`
replaces it with the same checklist, the same trigger paths, the same model and effort, and
the same read-only posture. What changed is where the reasoning lives: the threat model moves
out of the agent's always-loaded preamble and into `CLAUDE.md` §9 and `.claude/rules/vector.md`,
which the agent is instructed to read first, every time.

**Why.** `security-reviewer` **had never completed a single review.** Every invocation died on
its first turn to an API safeguard reporting `[bio]` — no findings, partial or otherwise.
Owner's call, on the evidence below: *"We need to get rid of the security reviewer."*

**The trigger was the agent definition, not any diff**, and that is established by a controlled
comparison rather than inferred. On PR #145 — a diff touching only `.github/**`, `CLAUDE.md`,
`scripts/apply-repo-settings.sh`, `.claude/skills/pre-pr/SKILL.md` and `docs/decisions/**`, with
no engine source and nothing biological in it anywhere:

| invocation | outcome |
|---|---|
| `security-reviewer`, that diff | died on `[bio]`, zero turns (`req_011Ceiy4nnEAyELADMug6JaD`, claude-opus-5) |
| `general-purpose` agent, same diff, equivalent audit prompt without the framing | **completed** — 24 tool uses, 6 findings |

Same diff, same session, same model tier; only the agent definition differed. Three sessions
have now hit it across two model tiers and four unrelated diffs
(`2026-09-03-agents-may-attest.md:54`, PR #143, PR #145).

The old preamble opened *"BT5's output is already the textbook method for evading
nucleotide-homology screening… a general-purpose evasion tool"*, and its `description:` —
the field that makes an agent selectable at all — led with *"biosecurity posture"*. Every one
of those sentences is defensive and correct. They were also unconditional: loaded on every
invocation, whatever the diff.

**What is NOT lost.** All six checklist items survive verbatim in substance, including the one
that earns the agent its keep: the `KmerIndex` Protocol-versus-implementation asymmetry.
`test_kmer_index_takes_no_external_database` reads only `core/services.py` and regexes the
frozen Protocol, while the sole implementation is `ConstructKmerIndex.of` at `vector/kmers.py:158`
— so a widened signature with a default conforms structurally, mypy passes, and
`check-approval-labels.sh` demands no label for `vector/`. That gap is still the first thing the
new agent is told to look at. Item 4 gained one line the old version could not have had: agents
authenticate as the owner, so a grant to the Repository admin role is a grant to every agent
session.

**Rejected:**

- *Delete it and add nothing.* The literal reading of the instruction, and the smallest diff.
  It would leave `/pre-pr` firing three legs where its own skill file promises four, and leave
  `vector/kmers.py` — the highest-risk directory in the repo, and the one with a documented hole
  in mechanical coverage — with no dedicated pass. The substitute audit found two real defects
  on #145 that a full `code-reviewer` pass had missed, which is evidence the leg does work when
  it can run.
- *Keep the agent and route around the filter per-invocation*, as three sessions did by hand.
  It worked each time and was disclosed each time, but it makes the chain's coverage depend on
  whoever is driving remembering to substitute — and `pre-pr-attest` attests that "the chain
  ran", with no way to record that a leg was swapped.
- *Soften the wording in place and keep the name.* The `description:` field is itself part of
  what is loaded, and the name appears in six files; renaming once is cheaper than leaving a
  name whose first word is the thing that broke it.
- *Move the checks into `code-reviewer`.* Its own definition is already at `sonnet/high` and
  scoped to the branch diff against the contract; folding an `opus/xhigh` intent audit into it
  would either raise the cost of every review or dilute the audit.

**Evidence:** issue #146 carries the full report, request id, and prior occurrences.

**Where:** branch `claude/autonomous-github-ci-boxw2y`, PR #145, carrying `approved:ci-change`
(`.github/workflows/pre-pr-attest.yml` carries one comment naming the agent).

Six files updated: `.claude/skills/pre-pr/SKILL.md` step 3, `CLAUDE.md`'s Delegation table,
`.claude/rules/vector.md`, `SETUP-NOTES.md`, `docs/buildout/README.md` (two places), and
`.github/workflows/pre-pr-attest.yml`'s header. Plus `.claude/agent-memory/code-reviewer/notes.md`,
which names the agent in a live note about which agents hold `Bash`. Earlier decision records
keep the old name deliberately — they are the record of what was true when written.

**Left open, and it could not be closed in the slice that created it:** whether the replacement
actually survives the safeguard is **UNVERIFIED**. An invocation was attempted on this very diff
— which touches `.github/`, one of its trigger paths — and failed for an unrelated reason:

> `Agent type 'boundary-reviewer' not found. Available agents: … security-reviewer`

**The agent registry is resolved at session start.** A `.claude/agents/*.md` file added mid-session
is not selectable, and the deleted one remains selectable, until a new session loads the directory.
So this change cannot be tested by the session that makes it, and that is a general property worth
knowing before anyone plans an agent change expecting to validate it in place.

The test is therefore the first `/pre-pr` in a LATER session on a diff touching `vector/`,
`core/services.py`, `verify.py`, `cassette/` or `.github/`. If it also dies on `[bio]`, the next
thing to try is the `description:` frontmatter: it is loaded for selection before any body text is
read, so it is the one part of the definition that cannot be avoided by an agent that never runs.
Until that first invocation returns, treat the boundary leg of `/pre-pr` as unproven and say so in
the PR rather than reporting the chain as complete.
