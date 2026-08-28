# Contract amendments

`packages/engine/src/bt5/core/` is frozen. The freeze is not "core never
changes" — it is that a change to core costs an argument, because after Wave 1
every lane and every one of the ~45 rule files is built on it.

The gate is `contract-freeze` in CI. It classifies your branch against `main`:

| | |
|---|---|
| **MINOR** | a new type, a new **defaulted** field, a new enum member, a field that gains a default |
| **MAJOR** | a removed or renamed anything, a changed annotation or default, a field that **loses** its default, a changed signature, a new protocol method |

The rule underneath is one question: **who breaks?** A field added with a
default breaks nobody — every existing constructor call still works. The same
field added without one breaks every caller at once.

Note that dataclass fields and protocol methods classify in opposite
directions, because the roles are opposite. BT5 *constructs* `Breach`, so a new
defaulted field is free. BT5 *implements* `FoldEngine`, so a new method lands on
every implementer — including every lane's fake.

## MINOR: the fast path

```bash
python tests/contract/regenerate.py
```

Commit the updated `tests/contract/manifest.json` and fixtures with your change.
The PR needs `approved:contract-change`. That is all.

## MAJOR: the amendment

1. Write `docs/rfcs/NNNN-short-name.md` from `TEMPLATE.md`.
2. Ship a **deprecation shim** that keeps the old form working, and the
   **two-window rule**: the old form stays accepted for two release windows.
3. `python tests/contract/regenerate.py`, then bump `contract_version` and add
   to `amendments` in the manifest:

   ```json
   { "version": 2, "rfc": "docs/rfcs/0001-short-name.md", "summary": "one line" }
   ```
4. `pytest tests/contract` must pass **without regenerating the fixtures**.
   That is the point of the shim: a fixture recorded under the old contract has
   to still construct. Regenerating it does not make an old caller work — it
   only stops recording that it doesn't.

`regenerate.py` deliberately does not write the amendment entry for you. A MAJOR
change is supposed to cost a paragraph of thought, and generating that paragraph
would defeat requiring it.
