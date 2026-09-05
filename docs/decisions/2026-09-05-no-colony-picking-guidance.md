## 2026-09-05 — BT5 does not tell the user how many colonies to pick

**Decided:** remove the screening-burden feature entirely. `ScreeningBurden`,
`screening_burden`, `ERROR_FREE_BP`, `ERROR_FREE_LAST_VERIFIED`, `DEFAULT_CONFIDENCE`,
`QcReport.burden`, `build_report`'s `vendor` parameter and the report's "Synthesis
screening" block are gone, along with the `error_free_bp` provenance block. Owner's
call, taken on scope grounds, on PR #131:

> The program does not need to suggest to the user how many colonies to pick. That is
> beyond the scope of the program.

### What it did, and why it went

The report sized colony picking from vendor fidelity: P(perfect clone) = exp(−L/E),
then n = ceil(ln(1−confidence) / ln(1−P)) colonies for 95% confidence. The arithmetic
was right and the vendor figures were sourced. That was never the problem.

The problem is that it is the only number BT5 computed for the user's **bench** rather
than for the **construct**. Everything else in the report is a statement about the
sequence BT5 designed — its rank against a random-synonymous null, which objectives
could not be evaluated, what no codon can fix. "Pick 14 colonies" is a protocol
instruction, downstream of a synthesis order BT5 does not place, at a confidence level
BT5 picked on the user's behalf. A lab decides that from its own cloning history.

### The immediate trigger

PR #131 would have added IDT's published gBlocks figure, the last one missing. Because
`DEFAULT_VENDOR = "idt_gblocks"`, that single value would have switched the number on
for every default-configuration run — the first time BT5 would have told a typical user
a colony count. Reviewing that PR is what surfaced the scope question, and the answer
made the PR moot. #131 is closed unmerged; #56, which asked for the figure, is closed
`not_planned`.

### What changes for the user

- The report no longer carries a "Synthesis screening" section.
- **One forced source of incompleteness is gone. A default `design()` run is still
  `is_complete == False`.** `build_report` used to append a degradation on *every*
  default run, because `DEFAULT_VENDOR` is gBlocks, gBlocks had no figure on file, and
  `is_complete` reads `degradations`. That is gone: `build_report` over a `DesignResult`
  carrying no degradations now returns complete. A real `design()` run does not — it
  still carries the unscreened biosecurity verdict and whatever else `design()` records,
  which is why `tests/design/test_increment.py` continues to assert
  `res.report.is_complete is False`. The change is that BT5 no longer reports itself
  incomplete over a number it has stopped trying to compute; it is not that runs are now
  complete.

  *Corrected 2026-09-05, after this decision merged as #150.* The commit message and
  #150's body both state the stronger claim — "a default run is now `is_complete ==
  True`" — and that is wrong. It was verified against a hand-built `build_report()` call
  rather than a `design()` run, which is exactly the distinction the sentence needed to
  make.
- Nothing about the designed sequences changes. No solver, rule, weight or ranking is
  touched; this removes a reporting line, not a design behaviour.

### How this sits with §0

The honesty posture bans reporting a number nobody measured. This is the adjacent case
and it is worth naming: the burden figure *was* measured and sourced, and BT5 still
should not have shipped it, because a well-sourced number outside the tool's competence
is its own kind of overreach. The old code was scrupulous about the gap — it returned
None for gBlocks rather than borrowing eBlocks' 5000, and said so in a degradation. All
of that care went into a question BT5 should not have been answering.

### Deliberately NOT done

**`brief.md` keeps E11 and its grade-A evidence, and item 11 of the differentiator
list.** The brief is the science record, not the product spec: the literature on
screening burden is unchanged and deleting a finding because of a scope decision would
corrupt the corpus. What changed is BT5's scope, which is recorded here and in
`docs/PLAN.md`, `docs/buildout/wave2/README.md` and `docs/buildout/s4-rules-liabilities.md`
— the three places an agent looks for what to build next.

### Lane note

`rules/_provenance.json` is M4's file and `score/` is M3's, so this change crosses a
lane boundary — the same boundary flagged on #131. It is not separable: the sidecar's
`error_free_bp` block exists only to source `ERROR_FREE_BP`, and leaving it would have
the repo documenting the provenance of a figure no code reads, with the seam test that
kept the two honest deleted alongside the code. Both halves land together, or the
sidecar rots. `ABSENT` stays defined in the sidecar's `_about.policy`: it is a general
provenance status and remains valid for the next absent figure, even though this change
removes its only current use.

### Tests

`packages/engine/tests/rules/test_vendor_error_free.py` is deleted whole. Every test in
it guarded the fidelity figures or the M3/M4 seam that carried them, with **one
exception worth naming**: `test_absent_is_named_in_the_about_block` guarded the
`_about.policy` prose, which this decision deliberately keeps. That prose is therefore
unpinned — nothing now fails if someone deletes the ABSENT definition. The policy text
itself now says so, and says a future ABSENT entry must bring its own test. In
`tests/score/test_report.py`, `TestScreeningBurden` and two further tests go the same
way. This is not §4 suppression: the covered behaviour no longer exists. Two tests were
kept and adapted rather than deleted, because they cover something that survives —
`test_a_run_with_everything_evaluated_is_complete` (which named a vendor only to dodge
the gBlocks degradation) and two guards in `tests/design/test_increment.py` that share one table.
`DEGRADATION_SOURCES` loses its screening-burden entry, which
`test_no_degradation_arrives_unremarked` reads as its closed set of recognised
sentences; and `DEGRADATION_SOURCE_FILES` loses `score/report.py`, which
`test_each_fragment_anchors_its_pattern_to_a_real_source_sentence` reads to check every
fragment still appears in an emitting module. That second check is the one that would
have failed outright, against a sentence no module emits any more.
