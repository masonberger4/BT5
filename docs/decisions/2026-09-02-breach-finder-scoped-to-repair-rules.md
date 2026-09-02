# `breach_finder()` walks only the HARD_REPAIR rules

**Lane:** M1 solver. Cross-lane for the session that found it (S1, design/score),
so it was filed as an issue first and landed on its own branch and PR, per
CLAUDE.md §1.

## What was wrong

`RuleSet.breach_finder()`'s docstring has said since it was written that it is
"SCOPED TO THE HARD_REPAIR RULES", and gave the reason: the callable runs once
per repair candidate, up to `max_candidates` per iteration, so every rule it
touches is paid for hundreds of times per iteration. It even named E8's k-mer
index and B1's fold as the waste to avoid.

The body was `return lambda c: self.findings(c).repairable`, and `findings()`
walks every spec in the set. The docstring described an intent the code never
implemented. Nothing was *wrong* in the output — `.repairable` filters correctly
— every discarded rule was simply evaluated in full and thrown away.

It was cheap enough not to notice until #92 and #93 added rules. It surfaced as
S1's G7 gate: a 500 aa design took 23.5 s against PLAN's 10 s bar, of which
20.2 s was the sweep, which is `optimize()` calls, which is `breach_finder()`.

## Measured

500 aa protein assembled into the 3.1 kb synthetic lentiviral backbone, against
the post-#92/#93 catalog (23 specs, 8 HARD_REPAIR), fold engine absent:

| | per repair candidate |
|---|---|
| unscoped (`findings(c).repairable`) | 65 ms |
| scoped (the eight specs it consumes) | 7 ms |
| **waste** | **~9x** |

Two rules were most of the discarded 57 ms: `f3_inverted_repeats` at ~29 ms and
`e8_kmer_uniqueness` at ~20 ms — the second being the exact rule the docstring
named as the thing not to evaluate here.

## Why the narrowing is behaviour-preserving

`repair_specs()` selects the specs where `enforcement_for(slot)` is HARD_REPAIR
for some active slot. `_enforcement_of` returns only values drawn from that same
set. So a spec outside `repair_specs()` cannot produce a breach that reaches
`.repairable`. The two paths must agree breach for breach, and a test asserts
exactly that rather than trusting the argument.

**With one exception, handled explicitly and not inherited.** `_enforcement_of`
falls back to the `enforcement` ClassVar when there are no active slots, while
`repair_specs()` is empty there because its `any()` is over nothing.
`DesignContext` requires a slot but does not require one to be ENABLED, so a
fully disabled context would silently lose breaches. That branch keeps the full
walk.

## Rejected

- **Caching `findings()` per construct.** Repair generates a fresh candidate each
  call, so the cache would never hit; it would add a keying bug surface for
  nothing.
- **Making the expensive rules cheaper** (`f3`, `e8`). Worth doing on its own
  merits, and orthogonal: they should not be evaluated in this loop at all,
  however fast they get.
- **Changing `findings()` too.** `advise()` reads its `hard_check` and the design
  lane reads its `evaluations` for the scorecard. Both need every spec. It is
  deliberately untouched.
- **A timing assertion in the test.** Flaky on shared CI, and it would not name
  what broke. The tests assert *which specs are evaluated*.

## Follow-on

S1's G7 gate (`test_a_500_residue_design_meets_the_g7_budget`) is red on #89
pending this landing on main. It is not a fix for the whole budget on its own —
it removes ~9x from the dominant term, which is what makes the 10 s bar
reachable.
