---
name: cheap-pass
description: Fork one identical mechanical edit across many files onto batch-editor, off the main context.
argument-hint: "<explicit file list> :: <exact before text> -> <exact after text>"
context: fork
agent: batch-editor
disable-model-invocation: true
---

# /cheap-pass

## Abort unless the arguments are complete

`context: fork` injects **this file** into `batch-editor`. It is **not** documented to
carry the conversation history, so the agent knows only what is written here plus
`$ARGUMENTS`.

If `$ARGUMENTS` does not contain **both** an explicit file list **and** the exact before
and after text, stop and reply `ABORT: cheap-pass needs an explicit file list and the
exact before/after text`. Do not guess which files were meant, and do not expand a glob
yourself — a plausible edit to the wrong files is the expensive failure here.

## The request

$ARGUMENTS

## Rules

- Apply the edit **byte-identically** to every listed file.
- Refuse any path under `packages/engine/src/bt5/core/`,
  `packages/engine/src/bt5/verify.py`, `tests/contract/`, `tests/invariants/`,
  `tests/data_integrity/`, `data/`, `benchmarks/` or `.github/` — those need an
  `approved:*` label and a human decision (CLAUDE.md §2).
- Refuse if the edit originates from a `bt5/core/` signature change: that is a contract
  amendment needing `/contract-change`, not a batch edit.
- Report **every** requested file as `APPLIED` or `SKIPPED` with a reason. A file left
  out of the report reads as done; a half-applied refactor that looks complete is the
  failure mode.
- Do not verify by running tests. `/pre-pr` owns verification.
