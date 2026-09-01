## 2026-09-01 — `expand_forbidden` supports IUPAC, with a solver-local expansion table

**Decided:** `LatticeTerms.forbidden` is documented IUPAC (`core/spec.py:196`) but
`expand_forbidden` assumed ACGT, so a degenerate base died as a bare `KeyError` on
`lattice._BASE_INDEX`. Option 1 from #73 — make the engine match the contract — over
option 2 (narrow the docstring to ACGT), because option 2 edits `core/spec.py` and needs
`approved:contract-change` plus an owner merge for what is an M1 bug.

Three parts:

1. **Expansion precedes the reverse-complement closure**, and the order is load-bearing.
   `core.types.reverse_complement` is `str.translate`, which leaves an unmapped character
   ALONE, so revcomp-first does not fail on a degenerate pattern — it half-complements.
   `RGATC` comes back `GATCR`, expanding to {GATCA, GATCG} when the true complements are
   {GATCT, GATCC}: the wrong sequences forbidden, the right ones permitted, silently.
   Pinned by `test_expansion_precedes_the_reverse_complement_closure`.
2. **`MAX_PATTERN_EXPANSION = 1024`, per pattern**, computed from the code widths before
   any string is built. Five fully-ambiguous positions or ten two-fold ones; every
   consensus motif a rule would plausibly declare (`MAGGTRAGT` 4, `GCCRCCATGG` 2, `GGWCC`
   2) sits orders of magnitude below. It excludes the N-spacer interrupted palindromes
   (BstXI `CCANNNNNNTGG` 4,096; XcmI `CCANNNNNNNNNTGG` 262,144) — those are wildcards,
   not ambiguity, and enumerating them is the wrong mechanism. The cost is paid twice
   downstream: an Aho-Corasick state per pattern prefix, then a 64-codon transition row
   per state in `lattice._codon_transitions`, so an all-N 8-mer reads as a hung solver.
3. **A `ValueError` naming the pattern** on a non-IUPAC character, an empty pattern, and
   an over-cap expansion. Eliminating the bare `KeyError` was the point of #73 either way.

**Rejected:**
- *Importing `verify.IUPAC_EXPANSION` / `verify.expand_iupac` into the solver.*
  `tests/data_integrity/test_oracle_independence.py` exists to keep the oracle off every
  lane's code path; sharing the expander defeats it pointing the other way. One
  transposed row would then be invisible, because the design and its check would forbid
  the same wrong set — the differential stops differentiating. The solver keeps its own
  `IUPAC_CODES`, and `test_the_solver_and_the_oracle_agree_on_expansion` asserts the two
  agree on a shared vector of eight patterns, so divergence fails loudly instead of
  silently. Two tables plus one test beats one table and no check.
- *A cap on the TOTAL expanded set rather than per pattern.* The blowup is exponential
  per pattern and merely linear across them; a rule listing fifty pure-ACGT six-cutters
  is fine today and must stay fine. A per-pattern cap is also the only one whose error
  message can name the offender.
- *Option 2 (refuse IUPAC, narrow the docstring).* Protected path, owner merge, and it
  removes a capability the contract advertises rather than supplying it.

**Scientific impact: none.** The only two rules populating `forbidden` today —
`d1_restriction_sites` (six-cutters) and `e1_homopolymers` (A/G runs) — emit pure ACGT,
which expands to itself. `expand_forbidden(["GGTCTC"]) == ("GAGACC", "GGTCTC")` still
holds. No shipped sequence changes.

**Follow-up for the design lane:** #71's `build_catalog` guard asserting lattice terms
are ACGT-only and raising `DesignError` is obsolete once this lands, and can be dropped.

**Evidence:** `packages/engine/src/bt5/solver/reference.py:44` (`IUPAC_CODES`), `:80`
(`MAX_PATTERN_EXPANSION`), `:83` (`expand_iupac`), `:126` (`expand_forbidden`); the three
consumers call it unchanged at `reference.py:185`, `lattice.py:145`, `repair.py:449`, and
are unaffected because expansion preserves pattern LENGTH, which is all
`_creates_forbidden`'s window bound and `repair.py`'s `guard_len` read. `bash
scripts/gates.sh` all eight gates exit 0; `HYPOTHESIS_PROFILE=ci` on `tests/invariants`
and on the two new properties.

**Where:** branch `claude/iupac-forbidden-expansion-orzxz4`, closes #73.
