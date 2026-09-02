## 2026-09-02 — The biosecurity screen is dropped; DNA synthesis vendors already screen orders (S2)

**Decided:** Owner determined BT5 does not need its own protein-level biosecurity screen
(`cassette/screen.py`, a `commec`-backed `Screen` protocol with `NullScreen` degradation
and `guard_emission`). DNA synthesis vendors already run their own screening on
incoming orders before synthesis, which is where the actual synthesis risk is gated —
BT5's own screen would have been a redundant, weaker second layer (no reference database
in CI, degrading to `not_run` on every run this repo could exercise). Reverted:
`cassette/screen.py`, `packages/engine/tests/cassette/test_screen.py`, and the prior
decision record `2026-09-02-biosecurity-screen-degradation.md`. `cassette/__init__.py`
restored to its pre-PR-#87 state (envelope exports only).

**Unchanged:** `core/context.BiosecurityVerdict` — frozen, not this lane's to touch, and
out of scope for this reversal regardless. Any future screening work is a fresh decision,
not a resumption of this one.

**Where:** branch `claude/s2-biosecurity-screen`; lane M8, `cassette/` only. PR #87
closed unmerged.
