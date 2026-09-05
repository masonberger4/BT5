---
name: contract-change
description: Classify a bt5/core change MINOR or MAJOR, then regenerate the contract record in the correct order.
disable-model-invocation: true
allowed-tools: Bash, Read, Grep, Agent
---

# /contract-change

## Classify BEFORE regenerating. This order is not stylistic.

`tests/contract/regenerate.py` **writes the manifest and the fixtures before it prints
anything, and returns 0 on every path** — including after its MAJOR warning. So its exit
code cannot tell you whether you owe an RFC, and by the time you read its output your
local baseline is already overwritten.

```bash
git show origin/main:tests/contract/manifest.json > /tmp/baseline.json
.venv/bin/python tests/contract/check_amendment.py /tmp/baseline.json
```

`check_amendment.py` exit codes: **0** clean, MINOR, or an amended MAJOR · **1** drift,
or a MAJOR without an amendment · **2** usage error.

Only once you know the classification, regenerate.

## MINOR — the fast path

A new type, a new **defaulted** field, a new enum member, a field that gains a default.
Nothing that exists stops working.

```bash
.venv/bin/python tests/contract/regenerate.py
```

Commit the updated `tests/contract/manifest.json` and fixtures with your change. The PR
needs **`approved:contract-change`**. That is all.

## MAJOR — the amendment

A removal, a rename, a changed annotation or default, a field that **loses** its default,
a changed signature, a new protocol method.

1. Write `docs/rfcs/NNNN-short-name.md` from `docs/rfcs/TEMPLATE.md`.
2. Ship a **deprecation shim** keeping the old form working, and honour the
   **two-window rule** — the old form stays accepted for two release windows.
3. Regenerate, then bump `contract_version` and add to `amendments` in the manifest:
   `{ "version": 2, "rfc": "docs/rfcs/0001-short-name.md", "summary": "one line" }`
4. **`.venv/bin/pytest tests/contract` must pass WITHOUT regenerating the fixtures.**
   That is the whole point of the shim: a fixture recorded under the old contract has to
   still construct.

`regenerate.py` deliberately does not write the amendment entry for you.

## Regenerating a fixture is not a fix

A fixture records that a value built under the old contract still constructs. Re-recording
it does not make an old caller work — it only stops recording that it doesn't. Same rule
as CLAUDE.md §3.9.

## The direction that catches people

Dataclass fields and protocol methods classify **oppositely**. BT5 *constructs* `Breach`,
so a defaulted field is free. BT5 *implements* `FoldEngine`, so a new method lands on
every implementer including every lane's fake. The question is always **who breaks?**

Uncertain? Run `/architect` — it orchestrates the call rather than guessing.

Full protocol: `docs/rfcs/README.md`. Labels: `approved:contract-change` covers both
`core/` and `tests/contract/`.
