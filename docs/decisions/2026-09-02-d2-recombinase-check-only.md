## 2026-09-02 — D2 is HARD_CHECK, scans both strands itself, and does not detect Bxb1

**Decided:** `d2_recombinase_sites` detects loxP-family, FRT and Gateway attB1/attB2
sites plus lone half-sites, as `Enforcement.HARD_CHECK`, with its own both-strand
wrap-aware regex scan and no `lattice_terms`.

Three choices worth the record:

### 1. HARD_CHECK, not HARD_LATTICE or HARD_REPAIR

`brief.md:98` grades the row **"H, check-only — 25-48 bp, never arise by chance"**.
`brief.md:90`'s occurrence budget makes "never" quantitative: a 25-mer is ~3e-9 per
1.5 kb. A loxP is in a construct because someone put it there, in a backbone the solver
may not edit, so every breach ships `fixable_by_codon_choice=False`. That is HARD_CHECK's
definition (`core/spec.py:38-41`) and the field's own docstring
(`core/spec.py:154-166`) warns that defaulting it True *"sends the solver after ... none
of which any codon can move -- and it exhausts the mutation space and reports infeasible
on a design that was fine."*

**Rejected:** *HARD_LATTICE with the half-sites in `forbidden`.* The 13 bp arms are pure
ACGT and would expand fine, but making them unreachable by construction is a
mutation-space cost imposed to prevent something that cannot happen by chance, and it
cannot touch the backbone occurrence that is the actual finding.

### 2. It scans both strands itself, and `CLAUDE.md` §3.4 still holds

§3.4 says to list forward motifs in `LatticeTerms.forbidden` and let the solver close
the set. That mechanism is unavailable here twice: the loxP and FRT patterns carry an
`[ACGT]{8}` spacer — 4^8 = 65,536 patterns against `MAX_PATTERN_EXPANSION = 1024`, and
`docs/decisions/2026-09-01-expand-forbidden-iupac.md` already settled that N-spacers are
*"wildcards, not ambiguity, and enumerating them is the wrong mechanism"* — and a
HARD_CHECK rule never reaches the automaton regardless. `d6_non_b_dna._hits` is the
existing precedent for exactly this (a G4 is not a finite motif set either), so this rule
copies its idiom, including the `% n` mapping, rather than inventing a second one.

### 3. Bxb1 is NOT detected, on purpose

`brief.md:99` says *"Bxb1 attB and attP both contain `GGTCTC` (BsaI) — flag this collision
explicitly"* but **gives no attB or attP sequence anywhere in the file**. So the rule
declares `conflicts_with = ("d1_restriction_sites",)` and reports the collision on any
detected site that carries `GGTCTC`, and does not claim to find Bxb1 landing pads.

**Rejected:** *writing the Bxb1 attB/attP sequences from general knowledge.* This is the
tempting one and it is precisely the defect `/verify-provenance` exists to catch — its
question is *"does the cited source actually support the number in the code?"*, and a
sequence recalled rather than cited answers no. A wrong landing-pad sequence would fail
silently in the direction that matters: reporting clean. **Follow-up for the owner:** add
the sequences to `brief.md` with a citation, then extend `SITE_PATTERNS`.

### Two reporting bugs found and fixed while testing

- **Double counting.** loxP's arms are exact reverse complements, so the whole 34 bp
  pattern matches its own reverse complement *at the same coordinates*. Keying the dedup
  on `(span, strand)` reported one physical site as two — telling the user they have twice
  the problem they have. The key is the span alone.
- **Wrapping containment.** The two arms of a complete site are that site described
  again, so they are suppressed — but `outer.start <= inner.start and inner.end <=
  outer.end` is wrong across the origin. A loxP stored as `[25, 59)` on a 46 nt plasmid
  covers residues 25-45 and 0-12, and its first arm is stored as `[0, 13)`, which that
  comparison rejects: the site's own arm came back as a lone half-site, on circular
  constructs only. Containment is now by residue set. Pinned by
  `TestCircular::test_a_wrapping_site_still_swallows_its_own_arms`.

**Evidence:** `brief.md:90, 98, 99`; `core/spec.py:38-41, 154-166`;
`docs/decisions/2026-09-01-expand-forbidden-iupac.md`;
`d6_non_b_dna.py:197-230`. `pytest packages/engine/tests/rules/test_d2_recombinase_sites.py`
29 passed; `pytest tests/data_integrity` 205 passed (194 on the merge base — 11 new
contract assertions for this spec).

**Where:** branch `claude/s4-rules-liabilities`, session S4 of the six-way buildout.
