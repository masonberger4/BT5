# Decisions

Settled decisions from working sessions. What survives compaction is what lives on disk —
a decision that exists only in a conversation is gone at the next compaction.

**One file per decision**, named `YYYY-MM-DD-slug.md`. This directory replaced a single
`docs/decisions.md` that every session appended to, which collided exactly as the rules
lane predicted: on 2026-09-01 two concurrent sessions appended to the same tail and PR #79
went `clean` → `dirty`, costing a hand-resolved merge. It is the same fix the catalog
already uses — *"Adding a rule edits ZERO shared files, not even an `__init__.py`"* — for
the same reason.

A second benefit: the old file declared "newest first" and had already stopped obeying it,
because a conflict resolution appended to the end. Filenames sort themselves, so there is
no ordering convention left to violate.

**Scope, against `docs/rfcs/`.** RFCs record *contract amendments*: they are load-bearing
in CI, `check_amendment.py` reads the manifest they correspond to, and a MAJOR change is
unmergeable without one. These files record *session decisions that no gate enforces* —
what was tried, what was rejected, and why. If a decision changes `bt5/core/`, it belongs
in an RFC, and the decision file just points at it.

**Format.**

```
## YYYY-MM-DD — one-line summary
**Decided:** what will happen.
**Rejected:** the alternatives, each with the reason it lost.
**Evidence:** file:line, a command's output, or a measurement.
**Where:** PR / branch / commit, if there is one.
```

Write the `Rejected` section even when it feels obvious. The alternatives that lost, and
why, are the part a future session cannot reconstruct from the code.
