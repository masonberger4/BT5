---
name: score-null-pr112
description: Verification techniques from reviewing PR #112 (score/null.py, host-frequency weighted null draw) — empirical RNG checks, caller-chain tracing, PR-description-vs-committed-test gaps
metadata:
  type: project
---

## "this code path never ran in production" claims: chase the caller chain to a filesystem/data mismatch, don't just read the guard
PR #112 claimed the weighted branch of `synonymous_variant`/`null_distribution`
never executed because no host resolves a codon-usage table (issue #98) —
checked out true: `design/runner.py::_host_usage` looks up
`data/codon_usage/{str(HostId)}.json` (e.g. `human.json`), but the files on
disk are named by reference-set, not host (`human_highly_expressed_refseq_w.json`
etc — a different lookup used by C1/C3's `*_REFERENCE_SET` maps). Every host
raises `FileNotFoundError` there today, so the weighted path is dead code.
General pattern: when a PR's "scientific impact: none" / mergeability argument
rests on "this branch is unreachable", grep every production caller (not just
tests) of the function, and check any filename/lookup-key assumption against
`ls` of the actual data directory.

## Inverse-CDF / `bisect_right`-on-`rng.random()` draws: verify empirically, it's cheap and conclusive
For a PR replacing `rng.choice(n, p=weights)` with `bisect_right(cumulative_cdf,
rng.random(), hi=len(cumulative)-1)`: don't just eyeball the math. A
`.venv/bin/python` scratch script (no repo files touched) computing chi-square
over ~100-400k draws both ways, plus checking `rng.bit_generator.state`
equality after one `rng.choice(n, p=p)` vs one `rng.random()` call with the
same seed, settles distribution-equivalence and stream-consumption in under a
minute — far more conclusive than reasoning about numpy's C internals from the
docstring. On #112 they came out bit-identical, stronger than the PR's own
hedged "same distribution, exact stream changes" claim.

## "verified rather than asserted" claims need a matching pinned test
PR #112's *PR description* (check with `gh pr view`, not `git show` — the
squash commit carries only the follow-up line) claimed a 400k-draw chi-square
check, "zero-weight codons drawn 0/60k", and "all-zero family falls back to
uniform" as verification — none of the three was a committed test. The only
weighted-path test used a fixture where every codon has weight 1.0 or 10.0,
never 0.0, so the zero-weight-skip and all-zero-fallback branches were
untested in CI. Flagged as non-blocking test-coverage gap (nothing was
skipped/loosened — coverage was just never added — so not a §4 suppression),
but called out explicitly since the description implied these were locked in.
Merged with the gap open; as of `369429a` the fixture is still `10.0 if
c.endswith("C") else 1.0`, so `null.py`'s all-zero-fallback arm has still never
run.
