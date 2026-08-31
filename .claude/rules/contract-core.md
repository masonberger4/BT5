---
paths:
  - "packages/engine/src/bt5/core/**"
  - "tests/contract/**"
---

# The frozen contract

Statement and labels are in root `CLAUDE.md` §2 and §2a. Full protocol:
`docs/rfcs/README.md`. This file is the rationale and the order of operations.

## Classify before you regenerate

`tests/contract/regenerate.py` **writes the manifest and fixtures before it prints
anything, and returns 0 on every path** — including after its MAJOR warning. Its exit
code therefore cannot tell you whether you owe an RFC, and by the time you read its
output the local baseline is gone.

```bash
git show origin/main:tests/contract/manifest.json > /tmp/baseline.json
.venv/bin/python tests/contract/check_amendment.py /tmp/baseline.json
```

Exit **0** clean / MINOR / amended MAJOR · **1** drift or unamended MAJOR · **2** usage.

Use `/contract-change`, which sequences this correctly.

## Why fields and methods classify oppositely

The gate asks one question: **who breaks?**

- BT5 **constructs** `Breach`. A new **defaulted** field breaks nobody — every existing
  constructor call still works. **MINOR.**
- BT5 **implements** `FoldEngine`. A new protocol method lands on every implementer,
  including every lane's fake. **MAJOR.**

The same field added *without* a default breaks every caller at once.

## Regenerating a fixture is not a fix

A fixture records that a value built under the old contract still constructs.
Re-recording it does not make an old caller work — it only stops recording that it
doesn't. For a MAJOR change, `pytest tests/contract` must pass **without** regenerating:
that is what the deprecation shim is for, and the old form stays accepted for two
release windows.

## Known deviation — report, do not fix

`CLAUDE.md` §3.1 says the genetic code table is explicit and **never defaulted** (NCBI
table 12 reassigns CTG to Ser; table 4 makes TGA Trp; a wrong table is a silently wrong
protein). Two call sites currently default it:

- `packages/engine/src/bt5/verify.py:308` — `table_id: int = 1`
- `packages/engine/src/bt5/vector/survey.py:150` — `table_id: int = 1`

`solver/pipeline.py:52` does it correctly. `verify.py` is the oracle and needs
`approved:oracle-change`, so this is a finding to raise, not a drive-by fix.

## Scope note

`tests/contract/` shares `core/`'s label rather than having its own: the two change
together, and demanding two labels for one coherent change is friction without signal.
Unprotected, a lane could re-record the manifest and make the freeze agree with whatever
it just broke.
