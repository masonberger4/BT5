---
paths:
  - "packages/engine/src/bt5/vector/**"
  - "packages/engine/tests/vector/**"
---

# The vector lane

The busiest directory in the repo — 38 file-touches in the last 200 commits, 3,944 LOC
across 14 files — and the one with the least mechanical protection.

## Biosecurity: the guard does not cover this directory

`CLAUDE.md` §9 bans a `KmerIndex` that accepts an external database, because pointing a
homology minimiser at an arbitrary target turns BT5 into a general-purpose
screening-evasion tool. Constraining the index to the assembled construct is what keeps
it from being one.

**The test enforcing that reads one file, and it is not this one.**
`tests/data_integrity/test_no_expression_claims.py::test_kmer_index_takes_no_external_database`
reads `packages/engine/src/bt5/core/services.py` and regexes the frozen `KmerIndex`
Protocol signature. The only implementation is `ConstructKmerIndex.of` at
`vector/kmers.py:158`.

A `database=` parameter added **here** would:

- leave the frozen Protocol in `core/services.py` untouched, so the test still passes;
- still satisfy structural conformance — a widened signature with a default conforms;
- require no approval label, because `.github/scripts/check-approval-labels.sh` has **no
  rule matching `packages/engine/src/bt5/vector/`**;
- and still satisfy `kmers.py:461`'s
  `_protocol_conformance: type[KmerIndex] = ConstructKmerIndex`. mypy became a required
  CI job in #63, so that assertion is now actually enforced — but a widened signature
  with a default conforms, so it catches a *broken* signature, not a *widened* one.

So `ConstructKmerIndex.of()` takes a `Construct` and nothing else, and that is a rule you
hold yourself. Any change to its signature goes to `boundary-reviewer` via `/pre-pr`.

## Circularity is the recurring defect shape

The construct is circular. Coordinates wrap the origin, and the assembled product has
insert/backbone junctions that neither fragment sees alone. Rules evaluate the whole
circular construct, both strands, including junction- and origin-spanning hits and
including inside the backbone — that is invariant I6 of the oracle. GC windows wrap the
origin too (I7).

When a vector test fails in a way that looks impossible, check the wrap first.

## Strand

Directional models read `slot.strand_of_interest` (`core/spec.py:274`). A hard-coded
strand 1 makes a reverse-oriented lentiviral cassette's polyA and splice analysis exactly
backwards, silently.

## Read these with `offset`/`limit`

`backbone.py` 28.9 KB · `kmers.py` 21.2 KB · `assemble.py` 18.7 KB · `gibson.py` 17.8 KB.
Reading one whole into the main window costs more than the answer. Use `Explore` to
locate and `Read` with a range, or delegate.

## The invariant that matters most

I9 — **every backbone base is byte-identical to the input backbone**. It is the
highest-value assertion in `verify.py`. Anything in `assemble.py` or `remap.py` that
could perturb a backbone base is a defect regardless of what the tests say.
