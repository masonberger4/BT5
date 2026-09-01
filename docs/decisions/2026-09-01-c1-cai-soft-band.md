## 2026-09-01 — C1 (CAI): a soft band that reports unavailable for seven of nine hosts

**Decided:** `c1_cai` ships as `Enforcement.SOFT` / `Direction.BAND`, band `(0.70, 0.90)`
per brief.md:77, `default_weight` 0.2 (matching what all three presets already assign
2.C1), `steering_weight` 0.0, `lattice_terms()` → `None`. 2.C1 now resolves in every
preset; only 2.C3 is left in `ResolvedPreset.unimplemented`.

The band's **ceiling is the operative half**. Kudla 2009 measured CAI at r = 0.14 (n.s.)
against 5' folding's r = 0.66 on the same 154 variants, and Welch 2009's deliberately
high-CAI control expressed at ~15% of the best variant. `Direction`'s own docstring names
this rule as the reason `BAND` exists.

`Evidence.CONTESTED`, not `EVIDENCE_BACKED`: the brief grades C1 "A (that it's weak)" —
grade A for the claim that CAI is a *weak* predictor — while Boel 2016 measured codon
content as 3–5x *more* influential than structure on a different readout. BT5 cannot
adjudicate that, and the citations carry both signs.

**Unavailability is the feature, not a caveat.** `data/codon_usage/` holds exactly one
reference set (`sharp_li_1987_ecoli_w.json`) against nine `HostId` values, so C1 computes
for the two E. coli hosts and reports the objective unavailable for HUMAN, HEK293, CHO,
S_CEREVISIAE, P_PASTORIS, SF9 and MOUSE — i.e. unavailable under `lentiviral_hek293` and
`aav_hek293`, a number only under `ecoli_expression`. NaN plus a reason-carrying breach,
copying `b1_five_prime._unavailable`.

**Rejected:**
- *A `codon_weights` lattice term, or a non-zero `steering_weight`.* Both are monotone
  pulls toward maximum CAI — the Tier-A DP maximizes what it is handed — which is exactly
  what "never maximized" (brief.md:73) forbids. A band is not expressible in the automaton
  either: it decides from a bounded codon suffix, and CAI is a geometric mean over the ORF.
- *Adding a mammalian codon usage table so the HEK293 path lights up.* `data/codon_usage/**`
  is protected (`approved:data-change`), and a human highly-expressed reference set is an
  evidence-bearing decision with its own provenance burden. Filed as its own issue.
- *Falling back to the E. coli table for a mammalian host.* It is the one table on disk, so
  the fallback would always succeed and always be wrong — a plausible-looking number
  measuring nothing. This is the failure the unavailable path exists to prevent.
- *Gating on `modality` (B1's choice) or on "host is E. coli".* Gating on host is the
  dangerous one: a lentiviral job propagates in E. coli and expresses in HEK293, so a
  host-keyed rule would find the one host with a w-table and report its CAI as the
  objective for a protein made elsewhere. The gate is `slot.role != "propagation"`.
- *Hard-coding ATG/TGG as the excluded single-codon families.* NCBI table 4 makes TGA a
  second Trp codon, so Trp carries information there. Family size is read from the injected
  table: `len(code.synonymous_codons(aa)) > 1`.
- *Importing `CodonUsage.cai` from M5.* `exp(mean ln w)` is recomputed in the rule so the
  lane keeps no import edge into the codon lane; `Services` is what decouples M4 from M5.
- *Editing `packages/engine/src/bt5/score/`.* The predicted ripple did not occur — the
  `unimplemented` assertions in `test_presets.py` use synthetic `fake_spec` lists, not the
  live registry. Nothing under `score/` was touched.

**Note:** `TableProvider.weights()` keys on the REFERENCE SET, not on `HostId` — every
caller in the tree passes `"sharp_li_1987_ecoli_w"`. C1 therefore carries an explicit
`CAI_REFERENCE_SET` map; a host absent from it has no reference set and reports
unavailable rather than raising `FileNotFoundError`. BL21 shares K-12's entry as a stated
same-species approximation, and the reference set travels in `detail["reference_set"]`.

**Where:** branch `claude/rule-c1-cai-enxfre`.
