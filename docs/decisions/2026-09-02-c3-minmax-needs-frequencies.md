## 2026-09-02 — C3 (%MinMax) ships computing nothing, because no shipped table carries usage frequencies

**Decided:** `c3_min_max` ships as `Enforcement.SOFT` / `Direction.BAND`, band
`(-100.0, 0.0)`, `default_weight` 0.3 (matching what all three presets already assign
2.C3), `steering_weight` 0.0, `lattice_terms()` → `None`, window 18 codons.
`MINMAX_REFERENCE_SET` is **empty**, so every one of the nine hosts reports the
objective unavailable with a reason. `ResolvedPreset.unimplemented` is now empty in all
three presets.

**The buildout prompt's premise was wrong, and this is the finding.**
`docs/buildout/s3-rules-translation.md:87` says "C3 %MinMax will hit the same wall" as C1
— unavailable for the seven non-E. coli hosts, computable for the two E. coli ones. It
cannot compute for E. coli either. %MinMax is defined on raw codon usage **frequencies**:
`%Max = 100·Σ(X_ij − X_avg,i)/Σ(X_max,i − X_avg,i)` needs each family's X_avg, X_max and
X_min. `data/codon_usage/` ships one file, Sharp & Li's relative adaptiveness w-index,
and `codon/tables.py:143-146` builds w as `(count + 0.5) / max synonymous (count + 0.5)`
— each family renormalised to its own peak, and the peak discarded.

Per family the differences survive that rescaling exactly: `w_ij − w_avg,i =
(X_ij − X_avg,i) / K_i`, where `K_i` is the family's peak. The pseudocount cancels. But
%MinMax **sums those differences across families**, and the per-family `K_i` that would
make the sum comparable is precisely what normalising to 1.0 threw away. So %MinMax
computed on w is a K-weighted quantity that is not the published metric, and equals it
only if every amino acid's peak count is identical.

**Why the band edges are safe.** Neither is a threshold anyone picked. 0 is Clarke &
Clark's own neutral point — "codon usage equal to the mean of all possible codon choices"
— so a ceiling there says only "do not push above what the host already does". −100 is
the metric's definitional minimum and therefore **cannot be breached**; the floor is inert
by construction, and deliberately so.

**The citation gap is closed by going to the primary source.** brief.md:79 cites nothing
for C3 and gives the window as a range, "z = 17–18 codons (CHARMING default 10)";
PLAN.md:661 pins 18 and attributes it to CodonTransformer. Clarke & Clark 2008 is the
metric's origin and says both "The resulting values are typically averaged over an
18-codon window" and "All results shown in Figures 2–4 used a window size of 18". So 18
is sourced to the paper that measured it rather than to a decision row.

**Rejected:**
- *Computing %MinMax from the w table anyway.* The failure mode
  `docs/decisions/2026-09-01-c1-cai-soft-band.md` already refused, one step over: w is the
  one table on disk, so the substitution would always succeed and always be a different
  statistic wearing this one's name and citation. If the owner prefers a computing-but-
  approximate C3, that is a legitimate call — but it is theirs, and it needs its own unit
  string and a citation that does not say Clarke & Clark.
- *`Direction.LOWER_IS_BETTER` on the mean %Max excursion.* Threshold-free and initially
  attractive, but a monotone "minimise %Max" objective has its optimum at an all-rare-codon
  sequence. C1's floor would normally bound that, except C1 is itself unavailable for the
  seven hosts the mammalian presets use — so the bound evaporates exactly where the presets
  weight C3. A band cannot run away.
- *A band with invented edges* (−30/+30 and similar). brief.md:79 gives C3 no threshold at
  all; any edge would pass all 11 contract assertions and be a number nobody measured.
- *Penalising %Min excursions.* brief.md:84 (C8) asks for ≥80% of native rare-codon clusters
  to be **retained**, and Clarke & Clark's actual finding is that rare codons cluster
  non-randomly at gene termini. A symmetric band would fight both.
- *Adding a frequency table.* `data/**` is S6's mutex and carries `approved:data-change`.
- *Adding a frequency channel to `TableProvider`.* A new protocol method is MAJOR under
  CLAUDE.md §2a and needs an RFC; `core/` is nobody's by default.
- *Leaving 2.C3 unimplemented.* The `unimplemented` docstring says the point is that a
  silently absent objective "reads exactly like a rule nobody has written yet". Shipping the
  rule changes the degradation from "no rule in this build" to "no usage-frequency reference
  set in this build", which is a materially more accurate thing to tell a user weighing a
  ranking — and the arithmetic lands complete and tested.

**Three things found on the way that are not mine to fix.**
1. `solver/catalog.py:151-153` consumes `LatticeTerms.forbidden` and nothing else.
   `codon_weights`, `codon_pair_weights` and `positional` are declared on the dataclass and
   read by no code in the tree. B8 and B9 both wanted `positional` and both ship
   `steering_weight = 0.0` rather than claim a nudge the engine does not perform.
2. `TableProvider.usage` is declared `-> Mapping[str, float]` (`core/services.py:137`) and
   `FileTableProvider.usage` returns a `CodonUsage` dataclass (`codon/tables.py:186`). C3
   codes against the protocol and validates the result, so the divergence surfaces as a
   stated unavailability rather than an AttributeError inside a rule.
3. `.claude/rules/rules-catalog.md:19-21` lists four brief sections — 2.B, 2.D, 2.E, 2.F —
   and omits **2.C**, which is where C3 lives. The documented `brief_ref` resolution
   procedure therefore fails for `2.C3` as written. One-line fix, in a file that governs
   S4's lane too, so it is reported rather than edited into a merge conflict.

**Also in this branch:** B8 (Kozak), B9 (out-of-frame ATG) and B2 (5' structure windows),
each with its paired test. Their own judgement calls are recorded in their module
docstrings and `weight_provenance`, per the catalog's convention: B8's weak/adequate/strong
ordinal is a *reading* of brief.md:68, which defines only "strong"; B9 is HARD_REPAIR with
`FIXED_POINT` because removing one ATG can create another; B2 scores only the proximal
window because brief.md:333 forbids averaging the two.

**Where:** branch `claude/s3-rules-translation`.
