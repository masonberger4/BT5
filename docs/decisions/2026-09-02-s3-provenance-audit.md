## 2026-09-02 — provenance audit of the four S3 rules: two citation defects, both fixed

**Decided:** `rule-auditor` ran over `c3_min_max`, `b8_kozak`, `b9_out_of_frame_atg` and
`b2_structure_windows` before the PR left draft, per `/verify-provenance`'s purpose —
catching WRONG provenance, which `tests/data_integrity/test_rule_contract.py` explicitly
cannot. Every threshold returned **SUPPORTED**. Two citations did not, and both are now
corrected on `claude/s3-rules-translation`.

**Defect 1 — C3 cited the wrong paper for terminal enrichment.** `c3_min_max`'s second
citation claimed rare codons are "enriched at the 5' and 3' termini of E. coli genes" and
pointed at Clarke & Clark 2008 (*Rare Codons Cluster*, PLoS ONE 3(10):e3412). That paper
reports non-random clustering; it does not report *where* the clusters sit. The
terminal-enrichment finding is the same group's follow-up — Clarke & Clark 2010, *BMC
Genomics* 11:118 — verified against the source: "rare codon clusters are more likely to
appear at the 5' and 3' ends of E. coli genes, rather than non-terminal positions", with
15 of 26 prokaryotic ORFeomes showing significant 5' enrichment and 12 showing 3'.
Repointed to PMC2833160, year 2010, with a note that it is deliberately a different paper
from the 2008 one above. This also resolves a second finding — citations[0] and [1] had
been the same paper under two URLs, which made the evidence base look broader than it was.

**Defect 2 — B2 attached Cambray's URL to a claim Cambray cannot have made.**
`b2_structure_windows`'s third citation carried the LinearDesign / Expi293F claim
(LinearDesign worst of five schemes for a DNA transgene despite 2-fold lower normalized
MFE, no correlation between MFE and yield) under `https://www.nature.com/articles/nbt.4238`
— Cambray 2018. The claim text is faithful to brief.md:333, but that row carries no URL,
and Cambray 2018 could not discuss LinearDesign (2023) or a 2026 benchmark.
`grep -n "nbt.4238" docs/research/brief.md` returns exactly one line, 19, attached solely
to the 244,000-sequence claim. Repointed to the benchmark itself, named at brief.md:13:
`https://proteininnovation.org/2026/03/...`, year 2026. This is the class of defect the
skill exists to find — it would have passed all 11 contract assertions.

**A claim about B9 that was wrong, and is not a citation.** The PR's scientific-impact
section asserted B9 "changes what the app refuses to build". It does not, and the
enforcement path was checked rather than assumed: an in-frame ATG in the first 50 nt is by
definition Met and therefore forced, so the breach ships `fixable_by_codon_choice=False`;
`solver/repair.py:157` partitions unfixable breaches onto `RepairOutcome.advisory`, and
`verify.py` carries **zero** references to `Enforcement` — it refuses on generic invariants
plus `forbidden` motifs, and B9 publishes none. So the finding is reported, never chased,
and never blocks emission. What remains true and still belongs in front of the owner is
that B9 returns `passes=False` for any protein with a Met at residues 2–16, which is
common.

**Confirmed rather than changed:**
- C3's load-bearing claim that %MinMax cannot be computed from the shipped `w` table
  **holds**, and survived the strongest counter-argument available: if Clarke & Clark's
  `X_ij` were relative-within-family frequency it would be recoverable as
  `w_ij / Σ_j w_ij`, pseudocount and all, and the whole unavailable design would collapse.
  The paper closes it — its Figure 1 tabulates "**absolute** codon frequencies … from the
  entire E. coli genome". Absolute and genome-wide, hence not recoverable. The `K_i` the
  rule's docstring names is `X_max,i`, and `max_j w_ij = 1` by construction destroys it.
- C3's 18-codon window sourced to Clarke & Clark rather than to PLAN.md:661's
  CodonTransformer attribution is a **strict improvement**, both quotes verbatim.
- B2's `PROXIMAL_UPSTREAM/DOWNSTREAM = 30/30` and `DISTAL_START/LENGTH = 30/60` encode
  STR(−30:+30) and STR(+31:+90) correctly on **both** strands — the +1-with-no-zero
  off-by-one was specifically worked and is right.
- B9's literal first-50 reading is what brief.md:69 says: the row qualifies "out-of-frame"
  in its second clause and not its first, one row and one author.
- B9's `FIXED_POINT` is genuinely mandatory. Worked counter-example: `CTC|TGT` (Leu-Cys)
  contains no ATG, but the synonymous `CTA|TGT` creates an out-of-frame ATG across the
  codon boundary — and where the downstream codon is `TGG` (Trp) it cannot be recoded at
  all, so any preceding codon ending in A manufactures one.
- B8's Noderer figures, the purine-class reading of −3 (brief.md:68 writes `R`), and the
  claim that the brief defines only "strong" — all verbatim.
- B9 borrowing B10's effect size (PMC7100133) is legitimate **because the citation
  discloses the borrowing in its own text**. Silently citing it as B9's own evidence would
  not have been.

**Rejected:** *dropping the terminal-enrichment clause from C3 rather than repointing it.*
The finding is real and load-bearing — it is why the band's floor is inert — so citing the
paper that actually reports it is better than deleting the reasoning.

**Not fixed here:** B2 declares `engine_calibration` while coding no kcal/mol threshold of
its own, so `check_engine_calibration` constrains more than the rule strictly needs. Left
deliberately, with the inheritance now stated in the module docstring: the failure it
prevents is silent, and a rule reporting energies weighted alongside B1's is exactly where
a mixed-engine run would go unnoticed.

**Where:** branch `claude/s3-rules-translation`.
