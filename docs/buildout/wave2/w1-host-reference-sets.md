# W1 — host reference sets: make the host-frequency null actually engage (#98)

*Copy this whole file into a fresh Claude Code session on the BT5 repo.*

**Launch:** permission mode `acceptEdits` · model opus · effort **xhigh** (override) · runs
unattended.
**Do not run this in plan mode** — three rules sessions are running in parallel and the gate
is the draft PR. **You post a design note as the first comment on your draft PR before you
implement.** That is your plan-mode substitute; see below.

**Wait for W0 to merge before you cut your branch.** You own `codon/tables.py`, which is in
W0's encoding sweep.

---

You are fixing issue **#98**. One lookup keys on the wrong thing, and as a result the
host-frequency null has never engaged for any host in this project's history.

## Read this first

Run `/bootstrap`. Then `CLAUDE.md`, then `docs/buildout/wave2/README.md`, then issue #98 and
`docs/decisions/2026-09-02-ranking-increment.md` — the latter is where the previous session
recorded the mismatch you are fixing, and where the ~2.4 s measurement you will need lives.

Your branch: **`claude/w1-host-reference-sets`**. Cut it from a **freshly fetched** main,
after W0 has merged:

```bash
git remote prune origin
git fetch -q origin main && git checkout -B claude/w1-host-reference-sets origin/main
```

## The bug

`packages/engine/src/bt5/design/runner.py:252-262`:

```python
def _host_usage(host: HostId) -> tuple[CodonUsage | None, str | None]:
    ...
    return FileTableProvider().usage(str(host)), None
```

`FileTableProvider.usage` (`codon/tables.py:186-187`) resolves
`data/codon_usage/{host}.json` — so `human.json`, `hek293.json`, `cho.json`. What ships is
named by **reference set**:

```
data/codon_usage/
  cho_highly_expressed_refseq_w.json
  human_highly_expressed_refseq_w.json
  mouse_highly_expressed_refseq_w.json
  sharp_li_1987_ecoli_w.json
```

So every `HostId` misses, `FileNotFoundError` is caught, the null falls back to
`uniform_synonymous`, `SolveSpace.usage` is always empty, and `live_axes` therefore always
drops the `codon_adaptation` sweep axis. The ranking increment called that axis "provably
dead" — it is dead by a **lookup-key mismatch** as much as by absent data.

`rules/catalog/c1_cai.py:209` already has the correct map, `CAI_REFERENCE_SET`, which is why
CAI works for E. coli and mammals while the null never does.

**Nothing false ships today.** The degradation is reported honestly and the report correctly
says the null is uniform. What is wrong is that a substantial part of what PR #89 built
never actually engages.

## The architectural call, already made

**The map moves to `bt5/codon/tables.py` (M5). Do not re-litigate this; do implement it.**

Add beside `FileTableProvider`:

- `REFERENCE_SET_FOR_HOST: Mapping[HostId, str]` — the `HostId` → reference-set-stem map,
  moved from `c1_cai.py`.
- `FileTableProvider.usage_for_host(host: HostId) -> CodonUsage` — resolves through it and
  raises a **named** error, not `FileNotFoundError`, for a host with no reference set.

`FileTableProvider.usage(name)` keeps its current reference-set-keyed meaning. The addition
is purely additive and breaks nobody.

Then `c1_cai.CAI_REFERENCE_SET` becomes an **import** of `REFERENCE_SET_FOR_HOST`, so the
two can never diverge again. Keep the name — `c1_cai`'s pin test and its
`svc.tables.weights(reference_set, "cai")` path stay untouched.

Then `_host_usage` calls `usage_for_host`. That is the only line you may change in
`design/runner.py`; W2 owns the rest of that file.

**Why here.** `docs/PLAN.md:454` gives M5 "Genetic code tables, **host usage**, CAI/tAI/…".
The provider owns the filename convention, and the bug *is* that the provider keys on
reference set while its callers key on host — the translation table belongs next to the
thing it is a key for. `codon/tables.py` importing `HostId` from `core.context` is fine:
core is the bottom of the graph and has no back-edge.

**Rejected — import `CAI_REFERENCE_SET` from `c1_cai` into `design/`.** It inverts the
dependency against a load-bearing invariant: `core/registry.py` autodiscovers rule files
precisely so adding a rule edits zero shared files. A rule file is a *plugin*; making the
design path import from one means renaming or deleting a rule silently breaks the null. It
is also M11 reaching into M4, which the ownership matrix forbids outright.

**Rejected — `core/context.py`, beside `LOCKED_TRANSLATION_TABLE`.** They look alike and are
not. That is a biological invariant that cannot change without the world changing; this is
an **inventory of which JSON files currently ship in `data/`**, which changes every time
S6-style work adds one. In `core/` every new codon-usage file would need `/contract-change`
plus manifest regeneration plus fixture re-recording, welding `data/`
(`approved:data-change`) to `core/` (`approved:contract-change`) so one file addition needs
two owner labels. Wrong cost curve.

**Note what you are NOT doing:** `c1_cai.py:43` says "codon data arrives through `Services`,
never by importing M5." That constraint is about **data**, and it is why the rule imports a
*constant* rather than gaining a provider call. Do not add a `Services` call to `c1_cai`.

## Scientific impact — read this before you write the PR

**Non-"none", and it is the sharpest one in the wave.** Making the lookup resolve changes
the null from `uniform_synonymous` to `host_frequency` for six hosts, which changes every
percentile, which changes the ranking, which changes the emitted sequence. It also re-arms
the fourth sweep axis, so the *panel itself* changes.

Under `CLAUDE.md` §7b that is an owner-signed change, not a routine follow-up, and your PR's
scientific-impact section must say so in those terms. **You do not self-merge.**

## The design note, posted first

Post a short design note as the **first comment on your own draft PR, before implementing**.
This is the redirect opportunity plan mode would have given, at the point it is still cheap.
It must contain:

1. **The measured G7 time with the fourth axis live.** See the risk section below — this is
   the number the owner most needs before the implementation is sunk.
2. **The before/after percentiles** on the reference fixture, for at least one mammalian
   host.
3. Which hosts gain a real null and which still degrade (`S_CEREVISIAE`, `P_PASTORIS`, `SF9`
   have no reference set — S6 deferred them deliberately, and they must keep reporting
   `unavailable` with a reason, never borrowing another host's table).

## The risk that is probably going to bite you: G7

End-to-end today is **7.37 s** at 500 aa against `G7_SECONDS = 10.0`
(`packages/engine/tests/design/test_timing.py:47`).

You are re-arming the `codon_adaptation` axis. `docs/decisions/2026-09-02-ranking-increment.md`
measured that solve at **~2.4 s** when it merely rediscovered the unsteered design. With a
real usage table it will produce a *different* CDS and may trip repair at **~2.5 s**
instead. 7.37 + 2.4 ≈ **9.8 s**. You are very likely the first session to see that
assertion fail.

**The wrong fixes, pre-refused:**

- **Never raise `G7_SECONDS`.** The bar is the product requirement (`docs/PLAN.md:508`).
- **Never mark the test `slow`.** `.claude/rules/tests.md` records that `-m "not slow"`
  currently deselects nothing because the marker is applied zero times. Marking this test
  would make it the first use of a marker whose only effect is that `gates.sh` and CI never
  run it. A timing gate nothing executes is worse than no timing gate.
- **Never tune the knobs to get under the bar.** `DEFAULT_SWEEP_STEPS` stays 1 — the ranking
  increment measured that density buys nothing (3 distinct designs at 3, 4, 6 and 20 weight
  vectors) and that the only lever that matters is the number of solves.

**The right answer** if four live axes exceed the budget: report the measured number as a
**finding about the budget**. `docs/PLAN.md:508` gives G7's fail consequence verbatim —
*"re-allocate budget before rules multiply cost."* That sentence was written for exactly
this moment. Put the number in the design note and let the owner decide.

Also note W3–W5 are adding five rules in parallel, and every new rule costs a full catalog
pass per repair candidate. **Report your G7 delta in the PR** so the cumulative number is
knowable; W2 lands last and carries the total.

## Acceptance

- `design(hosts=[HostId.HEK293], …)` yields `ObjectiveScore.null_kind == "host_frequency"`
  and **no** usage degradation.
- `live_axes` returns four axes, **and a test proves the fourth axis produces a design the
  other three do not.** This is the inverse of the test the ranking increment wrote when it
  argued the axis was dead — that test asserted sweeping all four axes and sweeping the live
  ones produce the same design set. Find it, and replace it with its opposite. If the fourth
  axis turns out to *still* add no design even with a real table, that is a finding and it
  goes in the decision record, not a silently deleted test.
- A host with no reference set (`S_CEREVISIAE`, `P_PASTORIS`, `SF9`) still degrades honestly
  with a named reason, and **never** borrows another host's table. That failure mode is
  argued in `docs/decisions/2026-09-01-c1-cai-soft-band.md` and is the whole reason the
  `unavailable` path exists.
- `c1_cai`'s existing tests pass unchanged — the re-export must be transparent.
- **The first commit is a characterization commit** pinning today's numbers before behaviour
  changes, so the diff reads as a delta rather than a rewrite.
- A decision record under `docs/decisions/` with the before/after percentile table.

## Files

You own, for this PR only:

- `packages/engine/src/bt5/codon/tables.py`
- `packages/engine/src/bt5/rules/catalog/c1_cai.py` (the re-export, nothing else)
- `packages/engine/src/bt5/design/runner.py` — **`_host_usage` only**, and nothing else in
  that file. W2 owns the rest and lands after you.
- `packages/engine/tests/codon/**`, `packages/engine/tests/rules/test_c1_cai.py`,
  and the `live_axes` test wherever it lives (name it in the PR so W2 knows).
- `docs/decisions/2026-09-03-*.md`

**Do not** touch `data/` — it is a mutex nobody holds this wave, and you need no new file.

## Delegation

- **Locating `live_axes`' existing test and every `_host_usage` caller** → `Explore`.
- **The ~2.4 s measurement and the "provably dead axis" argument** → `docs-miner` on
  `docs/decisions/2026-09-02-ranking-increment.md`.
- **Gates** → `gate-runner`.
- **`/pre-pr` will fire `rule-auditor`** because you touch a rule file. It is opus/xhigh —
  budget for it. Your edit is a re-export, not a threshold change, so say that plainly in
  the PR to keep the audit focused.
- **Never bare-Read `design/runner.py`** (33 KB). Use `offset`/`limit` around line 252.

## Done means

- `bash scripts/gates.sh` → `ALL GATES PASSED`, or G7 fails and you have reported the
  measured number as a finding rather than tuned anything.
- Design note posted as the first comment on the draft PR, **before** the implementation
  commits.
- Scientific impact section says **non-"none"** and explains that percentiles, ranking and
  the panel all move.
- Draft PR, `/pre-pr` from the operator, **owner merge — you do not self-merge**.
- Issue #98 closed by the PR.
