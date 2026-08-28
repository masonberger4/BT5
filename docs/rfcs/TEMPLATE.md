# RFC NNNN: <short name>

- **Status**: proposed | accepted | superseded
- **Contract version**: <the version this amendment introduces>
- **PR**: <link>

## What breaks

Name every MAJOR change from `contract-freeze`'s output, and for each one say
which existing callers stop working. Be specific — "every rule that constructs a
Breach", not "some callers".

## Why the additive version is not enough

A MAJOR change is only justified when the MINOR form of it is worse, not merely
uglier. If a defaulted field would do, use a defaulted field. If the answer is
"both defaults are wrong", say why, in the terms a reader can check — this is the
argument `Breach.fixable_by_codon_choice` had to make.

## The shim

How the old form keeps working, and for how long. Two release windows minimum.
Name the recorded fixtures that exercise it.

## Scientific impact

What changes about the sequences BT5 produces. If nothing, say "none" — but say
it, because a contract change that silently alters output is the failure this
whole protocol exists to prevent.
