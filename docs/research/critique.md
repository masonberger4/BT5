# Completeness Review — Codon-Optimization App Dossier

The 12 agents produced strong coverage of *sequence science* and *tooling*. The gaps cluster in four places they collectively did not look: **the protein-side input layer**, **the objective-aggregation math that the whole UI promise rests on**, **biosecurity/regulatory**, and **validation without wet lab**. I verified the highest-stakes items by direct fetch (search budget was exhausted by the prior agents).

---

## TIER 1 — Would sink the product

### 1. Biosecurity: this app is, by construction, a screening-evasion tool. Zero agents covered it.

**Why it matters.** Across ~200 searches, the only biosecurity artifact in the entire dossier is "Twist rejects ccdB." Meanwhile the app's core function — produce a functionally identical sequence with maximally different nucleotides — is the textbook method for defeating nucleotide-homology-based synthesis screening. A max-divergence synonymous recode preserves 100% amino-acid identity while dropping nucleotide identity far below typical BLASTN flag thresholds. Agent 4 even specifies `AvoidBlastMatches` and agent 3 specifies whole-construct k-mer uniqueness: pointed at a user-supplied database, those are general-purpose homology minimizers.

Two facts make this urgent rather than theoretical. First, IBBIS states its `commec` screening tool was built partly to "boost resilience to AI-generated sequences, chimeric sequences, and subversion attempts," citing a **2025 *Science* study finding that most AI-generated toxin variants evaded existing screening tools**. Second, agent 11's architecture recommendation is explicitly *fully offline, no telemetry, no phone-home* — the exact design that maximizes misuse attractiveness and gives you zero visibility.

**Recommendation.**
- Bundle or optionally install **`commec` (Common Mechanism)** — MIT licensed, runs entirely locally after a one-time database download, sub-second per sequence, HTML/JSON/CSV output. It screens via HMM biorisk profiles **plus BLASTX protein-level homology**, which is precisely the layer that is robust to your own tool's output.
- **Screen the input protein, before optimization** — not just the output DNA. Protein-level screening is the only screening your product cannot itself defeat.
- **Never ship a "minimize identity to a reference sequence" objective**, and never let repeat/homology-avoidance machinery accept an arbitrary user-supplied target database. Constrain it to the assembled construct only.
- Write the design decision down. For institutional/pharma sales this becomes an asset (IBC and EHS sign-off), not a cost.

### 2. Genetic code table selection is absent entirely — a silent protein-changing bug

**Why it matters.** All 12 agents implicitly assume NCBI translation table 1, yet agent 2's shipped host panel includes organelles and non-standard-code organisms. Verified reassignments that break naive back-translation:

| Table | Context | Reassignment that breaks you |
|---|---|---|
| 12 | Alternative yeast nuclear (*Candida*) | **CTG = Ser, not Leu** |
| 4 | Mycoplasma/Spiroplasma | **TGA = Trp, not stop** |
| 6 | Ciliate nuclear | TAA/TAG = Gln |
| 2 | Vertebrate mitochondrial | ATA = Met; TGA = Trp; AGA/AGG = stop |
| 11 | Bacterial / plant plastid | Chlamydomonas & tobacco chloroplast — agent 2 ships these tables but never flags the code |

Failure mode: user selects a *Candida*-family host, the app emits CTG for Leu, the host translates Ser. The protein is wrong and no assay catches it until months later. Or the app rejects a valid mitochondrial construct as containing an "internal stop."

**Recommendation.** Make the translation table a first-class explicit parameter (`Bio.Data.CodonTable`, 30+ tables), defaulted *and locked* by host selection, printed in the report, and written to the GenBank `/transl_table` qualifier. Property-test `translate(back_translate(p, table), table) == p` for every shipped table. Cheap; eliminates a whole class of silent wrongness.

### 3. Objective normalization is unspecified — and it is the entire UI promise

**Why it matters.** The brief is "user-controlled trade-off weighting." Every agent describes a weighted sum; **none specifies how incommensurable objectives are made comparable.** You are summing kcal/mol (unbounded, negative, length-dependent), CAI in [0,1], motif counts (integers 0–50), GC deviation in percentage points, repetitive-9-mers-per-100bp, and AAV size in kb. Agent 7 quotes DNAChisel's `boost` and GeneOptimizer's `G_q`; neither source says where the numbers came from. Agent 8 offers codonGPT's weights (CAI 2.0, GC 0.2, ΔG 0.05), tuned on one gene.

Unnormalized, the sliders will be dead over most of their range while one term silently dominates. This is the single most likely reason the app feels broken to a user who cannot articulate why.

**Recommendation.** Normalize every objective to a **percentile against an empirical null of 200–500 random synonymous variants of *this* protein in *this* host**. This buys three things at once: unit-free [0,1] scores that sliders act on linearly; an honest report line ("94th percentile for 5′ accessibility among random synonymous variants of this gene"); and automatic length/composition correction, which agent 10 correctly demands for MFE and no one else applied. Agent 1's `codon-bias` package already ships the synonymous-permutation generators. Cost is a few hundred cheap evaluations, once per job.

**Corollary nobody raised: 90% of users will never move a slider.** The trade-off UI is the marketing feature; the **default weight vector is the product**. No agent proposed defaults with provenance. That is the highest-leverage unaddressed design decision in the dossier.

### 4. Repetitive proteins and antibodies: max-CAI is self-defeating, and no agent connected the dots

**Why it matters.** Agent 3 establishes that repeats — not GC — are the top predictor of synthesis failure (the two highest-Gini features in the best published model are longest-repeat length and repetitive-9-mers/100bp). Agent 1 establishes that max-CAI collapses to one codon per amino acid. **These compose catastrophically:** one-codon-per-AA back-translation of any protein with internal repeats produces *perfect* nucleotide repeats.

The proteins people actually express are the repetitive ones: antibodies and scFv/CAR constructs (VH and VL share framework homology; bispecifics and tandem scFvs repeat entire domains), Fc fusions, (GGGGS)ₙ linkers, 2A peptides used twice in one ORF, TALEs, zinc fingers, ankyrin/LRR repeats, collagen, elastin-like polypeptides, His and FLAG tags. **The word "antibody" does not appear anywhere in 12 agents' output**, despite being the dominant recombinant-protein use case.

**Recommendation.** Treat repeat-breaking as a **constraint that overrides codon optimality**, not an objective competing with it. Detect repeated amino-acid segments in the *input protein* before any codon is chosen, then deliberately assign divergent synonymous codons across copies — target ≤85% nucleotide identity between repeat copies and no exact match ≥15–20 bp (agent 6's RecA MEPS floor is 23–27 bp; agent 3's vendor floor is 20 bp). Ship a named "repetitive protein" mode and use an scFv/CAR as a tested example. This is a concrete, demonstrable capability no vendor tool advertises.

---

## TIER 2 — Agents contradicted each other; someone must adjudicate

I resolved the most consequential one by fetching the primary source.

| # | Contradiction | Adjudication |
|---|---|---|
| 1 | **Twist local 50 bp GC window.** Agent 3: High-Complexity trigger is <10% or >90%, and explicitly debunks "35–65%". Agent 4: "local 50 bp GC windows must be 35–65%", stated `confidence: high`, sourced to the Twist FAQ. | **Agent 4 is wrong.** I fetched the live Twist FAQ: the stated trigger is "Local GC content within a 50 bp window falls below 10% or exceeds 90%." Agent 4's numbers are folklore with a fabricated citation. Do not encode them. |
| 2 | Twist homopolymer: ≥14 bp hard reject vs >30 bp High-Complexity | Both appear on the live page simultaneously. Tier-dependent, genuinely unresolved by the vendor. Encode as two tiers, not one rule. |
| 3 | GC "extent"/spread: ≤50 points (agent 3) vs ≤52% (agent 4) | Neither appears on the live FAQ. Both trace to a PDF tech note. Treat as **unverified**; do not make it a hard fail. |
| 4 | **CAI implementation.** Agent 1: pseudocount 1, exclude Met/Trp/stops. Agent 2: Biopython uses pseudocount **0.5** and treats the three stops as a synonymous family. | These produce **different CAI values for the same sequence.** Your reported CAI will not match the vendor's, the literature's, or a reviewer's. Pick one, document it, and ship both as a toggle with the Sharp & Li 1987 *E. coli* table as a regression fixture (agent 2 gives exact values). |
| 5 | %MinMax window: 17 (agents 1, 2) vs 18 (agents 9, 12 — CodonTransformer hard-codes 18) vs 10 (CHARMING default) | Reproducibility issue. Pin one, expose it, print it in the report. |
| 6 | **Structure direction.** Agent 1: maximize 5′ ΔG (less structure). Agent 5/LinearDesign: minimize MFE globally. | Reconciled only by *context* (bacterial 5′ end vs IVT-mRNA body), which agent 10 gets right and agents 1 and 5 each state as a flat rule in opposite directions. **Never expose a single "structure" slider.** Two separate objectives on two separate windows, switched by host/modality. |
| 7 | Codon content vs 5′ folding beyond codon ~16: Boël vs Kudla/Goodman/Cambray | Genuinely unresolved in the literature (agent 1 flags it). Expose as a weighted hypothesis, not a rule. |
| 8 | **Chi site (GCTGGTGG).** Agent 4: folklore, and the primary literature says Chi is *protective* on plasmids. Agent 6: instability determinant (low confidence). Agent 12: in the default scan list. | Agent 4's reading is better sourced. Make it low-priority/optional and label it as such. Cost of avoidance is ~0.05 hits per 1.5 kb, so it is nearly free — but don't present it as evidence-backed. |
| 9 | **Three-way GC conflict.** Agent 6 wants ≥45% GC per 100 nt to suppress cryptic promoters. Agent 5's controlled AAV data shows the 34%-GC natural stuffer *beat* the 44%-GC designed one. Agent 4's CpG depletion drives GC down. Agent 3's vendors reject both extremes. | **Nobody reconciled this.** It may be infeasible for AT-rich proteins. Surface which constraint is binding per window rather than silently steering. |
| 10 | LinearDesign license: agent 7 says no LICENSE file (404); agents 8 and 10 quote `license.txt` verbatim | Agents 8/10 are more credible (specific quotation). Unusable either way — plus a Baidu patent filing. |
| 11 | G-quadruplex avoidance in mammalian mRNA | Agent 10 resolves correctly: on for synthesis feasibility and bacterial CDS, off for mammalian translation (Guo & Bartel: globally unfolded in cells). Agents 3 and 12 apply it flatly. Use agent 10's gating. |
| 12 | Tissue-specific codon tables | Agent 1 recommends shipping them; agent 2 shows the effect is contested (2006 MBE found none; between-tissue variance ≪ between-gene variance). **Agent 2's caveat should win** — opt-in, labeled contested. |

---

## TIER 3 — Missing scientific factors no agent covered

### 5. Introns — the most reliable mammalian expression lever, absent from the dossier
Adding a synthetic/chimeric intron (CMV IE intron A, β-globin/IgG chimeric intron, EF1α intron 1) is among the most dependable ways to raise mammalian transgene expression, and it interacts directly with the app: **agent 4's "remove all cryptic splice sites" rule will happily destroy a deliberately placed intron**, and the app has no concept distinguishing a wanted splice site from an unwanted one. Only agent 5 touches introns, and only via the β-globin RRE case.
**Recommendation.** Parse intron features from the backbone; exclude annotated introns from splice-site removal; offer intron insertion as a cassette option with the size cost surfaced in the AAV/lenti budget.

### 6. Protein-level liabilities created by the app's own tag/linker additions
Agent 12 handles tags and linkers well but misses that **appending a linker can create an N-glycosylation sequon (N-X-S/T, X≠P) across the junction**, silently changing a secreted product's glycoform. Also unaddressed: free cysteines introduced by tags, deamidation/isomerization motifs (NG, DG, DP) at junctions, and N-end-rule/Met-excision consequences of the residue at position 2 (which agent 1 is simultaneously optimizing for *expression* via the codons 3–5 motif — a direct conflict).
**Recommendation.** Run a protein-level liability scan on the *assembled* fusion, not the input protein. Flag junction-created sequons explicitly.

### 7. "It expressed but went into inclusion bodies" — the actual bench complaint
Agent 1 covers co-translational folding via %MinMax, but codon choice is a weak lever for solubility compared to fusion-partner choice (MBP, SUMO, Trx, NusA) — which the app is *already assembling the cassette for*. No agent suggested the app should recommend a solubility tag.

### 8. IVT-mRNA specifics beyond folding
If the app claims IVT support: T7 +1 requirements (G start, GG/GGG strongly preferred), 5′-end constraints for co-transcriptional capping, poly(A) encoding strategy, and — as agent 10 notes — no energy model exists for m1Ψ. Agent 10 flags the last one; the first three are absent.

### 9. "Don't optimize" as a first-class output
Agent 2's 2026 mammalian benchmark (native most consistent), agent 7's vendor coin-flip finding, and agent 9's mushroom-luciferase case all point the same way: for homologous mammalian expression, **the honest default recommendation is often "use the native sequence."** No vendor tool will ever say this. Making it a first-class, reasoned output is a genuine trust differentiator.

### 10. Harmonization mode is unimplementable as specified
Agents 1, 2, 5, and 12 all promote "Harmonize" as a headline mode (agent 2 makes it the mammalian default). But **harmonization requires the source organism's native CDS** — agent 1's own rule is "minimize distance between the %MinMax profile of the designed sequence and that of the native sequence in its native host." The app's stated primary input is a **protein**. With only a protein, there is no source CDS, no profile to match, and the mode is undefined.
**Recommendation.** Require a native CDS to enable harmonization; grey out the mode otherwise; never silently substitute something else. And report agent 2's finding honestly: harmonized was best for 8/18 targets but *the most variable*.

---

## TIER 4 — Implementation questions that will block a developer

| Question | Status |
|---|---|
| How are objectives normalized for the sliders? | **Unanswered** — see Tier 1 #3. Blocking. |
| What are the default weights and where do they come from? | **Unanswered.** Blocking (it is what 90% of users get). |
| How do you present a Pareto front over 5–6 objectives? | **Unanswered.** Agent 8 proposes a 2-D scatter with color for a third — that fails at expression + synthesizability + titer + stability + GC + CpG. Consider parallel-coordinates plus a small ranked candidate panel (agents 1 and 7 both independently recommend emitting 3–8 diverse candidates rather than one answer — take that seriously; it is also the honest response to a coin-flip evidence base). |
| Protein input validation | **Unspecified.** Ambiguity codes (X, B, Z, J), selenocysteine U, pyrrolysine O, lowercase, whitespace, FASTA headers, internal `*`, trailing stop, leading Met present/absent, unicode homoglyphs from PDF copy-paste. Needs a strict validator with actionable errors and an explicit "does your sequence include the initiator Met / the signal peptide?" prompt. |
| Signal-peptide cleavage-site prediction | **Blocked by licensing.** Agent 12 says an N-terminal tag must go *after* the cleavage site, but I confirmed **SignalP 6.0 is academic-download-only; commercial users are told to contact DTU**. DeepTMHMM is similar. Either use a permissive heuristic (n/h/c-region), make SignalP an optional user-installed dependency, or never auto-place an N-terminal tag without explicit user confirmation of the site. |
| Length/scale envelope | **Only specified to ~1 kb** (agent 8). Real targets blow past it: Cas9 1368 aa, Factor VIII 2332 aa (agent 5's own flagship example), full dystrophin 3685 aa, ApoB. Agent 9 notes LinearDesign needed >60 GB at 1450 aa; agent 10 extrapolates RNAfold to 2–7 s at 3 kb. State a ceiling, degrade gracefully (windowed-only structure above N), and test on FVIII and dystrophin explicitly. |
| Circular constructs / origin-spanning features | Agent 8 gives the tripling trick and agent 11 flags 1-based vs 0-based, but the interaction — an insert whose optimal fix lies inside immutable backbone across the origin — has **no documented resolution** (agent 8 lists this as open). |
| Feedback loop | **Absent.** Every agent concludes the evidence is weak and gene-specific; agent 1's best recommendation is per-user calibration. No agent specified the mechanism. Build a local, private, opt-in results log keyed by design hash, plus a calibration view. It is the only path to ever beating a coin flip — and it fits the offline/no-telemetry stance. |
| Design identity | Agent 8 notes DNAChisel has no seed; agent 11 recommends a JSON design record. Neither makes the bench-critical point: **two runs of the same protein producing two different sequences with the same name is how a lab ends up with two tubes and an irreproducible result.** Content-hash the final sequence; print a short hash on the report, the GenBank note, and the ordering FASTA header; make the tube-label string an explicit deliverable. |

---

## TIER 5 — Non-obvious failure modes in real bench use

1. **Silent protein change** from wrong genetic code (Tier 1 #2) or an unhandled ambiguity code.
2. **The user pastes a protein that already has a His-tag** (or already lacks the initiator Met), and the app adds a second one. Agent 12's Minotaor suggestion addresses detection — nobody made it a required pre-flight step.
3. **Backbone GenBank has no annotations, or wrong ones.** Extremely common for files from collaborators and SnapGene auto-annotation. Agent 12's detection ladder ends in a heuristic ("first ATG within 200 bp of the promoter") that will confidently pick the wrong site.
4. **Post-hoc manual editing.** User optimizes in your app, then adds a tag in SnapGene — recreating a BsaI site or breaking frame. Your guarantees evaporate silently. Mitigation: ship a re-validation mode that takes an edited GenBank back and re-checks every constraint.
5. **Copy-paste corruption** into a vendor order form (line breaks, spaces, a truncated selection). Mitigation: the checksum/hash on the ordering file, and a paste-back verifier.
6. **Vendor rejects the order** because your manufacturability model disagreed with theirs, and the user has no recourse path. Mitigation: emit the specific offending window and a one-click "relax this constraint and re-run."
7. **Recoding silently invalidates the user's existing assays** — qPCR primers/probes, siRNA/shRNA targeting the transgene, HDR donor homology arms, and Southern/sequencing primers all stop working. No agent mentioned this. It is a guaranteed, avoidable support ticket. Mitigation: accept a list of existing oligos and warn on loss of match.
8. **Internal ATG in a strong Kozak** produces an N-terminally truncated product that runs at the "wrong" size on a gel and gets misdiagnosed as degradation (agent 4 covers detection; nobody frames the bench symptom).
9. **Over-constrained failure** is the dominant real-world experience with constraint solvers (agent 7 says so explicitly). If the app just says "no solution," it is useless. Agent 8's minimal-conflicting-set + suggested relaxations must be a shipped feature, not a nice-to-have.
10. **The tool wins on metrics and loses at the bench** — because, per agents 7 and 9, optimizers are roughly a coin flip against native and 5–31% of protein-level variance is all any computable feature explains (Cambray). If the UI implies more certainty than that, the first failed construct destroys trust permanently.

---

## TIER 6 — Regulatory, legal, and licensing gaps

**Already found by agents** (consolidate into one attributions/compliance screen): ViennaRNA (no redistribution for a fee), REBASE (CC BY-NC — and inherited by `Bio.Restriction`), LinearDesign (no redistribution + Baidu patent), NUPACK (non-commercial, explicitly forbids GUIs), Codon Statistics DB / D-Tailor / mRNAid / DeepCodon / CodonBERT weights / RiboNN weights / Nucleotide Transformer (all NC), and the GPL cluster (pLannotate, CryptKeeper, OSTIR, Promoter Calculator, RNAstructure, Optimus 5-Prime).

**Newly identified:**
- **SignalP 6.0 / DeepTMHMM: academic-only download; commercial users must contact DTU.** Verified. Blocks the signal-peptide-aware tagging feature.
- **Export control.** Australia Group controls apply to genetic elements of listed agents. A design tool that outputs such a sequence may create obligations for the *user* — a warning is cheap and defensible.
- **Institutional biosafety (IBC) review.** No agent mentioned that constructs may require IBC approval. A "this design contains X" flag helps users, not hinders them.
- **Clinical/GMP context.** If a user designs a gene-therapy plasmid, FDA CMC expectations include full sequence documentation and justification of every element. The app is uniquely positioned to *emit* that documentation package automatically. Nobody suggested it; it is a straightforward premium feature. (Note: I could not retrieve the specific FDA guidance text — the landing page has only metadata. Verify before making claims.)
- **Patent FTO on codon optimization itself.** Agent 7 flagged CureVac/BioNTech GC-enrichment claims as unretrievable (patents.google.com returned 503 throughout). This remains genuinely open and should be reviewed by counsel before any mRNA mode ships.
- **Research Use Only disclaimer / no medical-device claim.** Absent from all 12 agents.

---

## TIER 7 — Validation strategy without wet lab (explicitly asked; no agent addressed it)

Six layers, in order of certainty. The discipline that matters most: **never validate with the same function you optimized with.**

**Layer 1 — Invariants (100% verifiable, no biology required).** Property-based testing over a corpus of ~10,000 UniProt proteins including pathological cases:
- `translate(output, declared_table) == input_protein` for every host × mode × random seed.
- Frame preserved across the full assembled cassette (agent 12's mod-3 invariant).
- **Zero forbidden motifs in the final assembled plasmid, checked by an independent validator written against a different code path** than the optimizer's scorer.
- GenBank/SnapGene round-trip: feature count, every (start, end, strand), every qualifier survives byte-identically (agent 11's golden-file suite).

**Layer 2 — Metamorphic and differential testing.** You don't know the right answer, but you know these relations must hold: adding a redundant constraint must not change the output; permuting constraint order must not change the output (**agent 8 shows DNAChisel fails this** — a ready-made bug class); re-optimizing an already-optimal sequence must be a no-op; a score must improve monotonically in the weight you raised. Differentially: run DNAChisel, CodonTransformer, and your engine on identical input — every disagreement is a bug candidate to triage.

**Layer 3 — Retrospective benchmarking on public *measured* data.** The dossier already names the datasets: Kudla (154 GFP variants, 250-fold range), Goodman (~14,000), Cambray (244,000), Verma (259,134), Boël (6,348), and CodonBERT's bundled mRFP/fungal/*E. coli*/stability sets. Protocol: **hold out by gene, not by variant** (variant-level splits leak). Report Spearman ρ *and* top-k enrichment — "of the 10 sequences my objective ranks highest, what fraction are in the measured top 10%?" That is the decision the user actually makes; ρ is not.

**Layer 4 — Negative controls the field mostly fails.** Your composite objective must (a) rank Kudla's known-bad variants below known-good ones, (b) beat CAI — a low bar the literature says most tools fail, and (c) **beat a random synonymous shuffle**. If it cannot beat random shuffle on a public dataset, ship it labeled honestly rather than quietly.

**Layer 5 — Adversarial inputs.** Poly-Q, poly-A protein, 100% Leu, a 2-residue protein, a 4,000-residue protein, a protein whose every synonymous option creates a BsaI site (proves infeasibility reporting works), a backbone with zero annotations, a feature spanning the origin.

**Layer 6 — The free external oracle nobody proposed.** Submit designs to the **Twist / IDT / GenScript order-entry complexity checkers** and record accept / complex / reject. This is a real, external, quantitative label for the entire manufacturability half of the app, obtainable **at zero cost without ordering anything**. Agent 3 identified reverse-engineering IDT's complexity score as "the single highest-value empirical exercise"; this is the same idea, generalized, and it is the best validation available to you. Build it into CI.

**State plainly what you cannot validate:** absolute expression prediction. Given Cambray's 5–31% variance ceiling, any UI that implies a predicted titer is dishonest. Report ranks, percentiles, and confidence bands.

---

## Suggested sequencing

1. **Now:** genetic code tables; protein input validator; objective normalization; the biosecurity decision (it constrains architecture).
2. **Before any UI work:** default weight vector with written provenance; the multi-objective presentation question; adjudicate the contradiction table above and encode the *sourced* numbers only.
3. **Before v1:** repetitive-protein mode with an scFv test case; vendor-oracle CI; infeasibility reporting; design hashing.
4. **Legal, in parallel (long latency):** ViennaRNA bundling permission, SignalP terms, mRNA-mode patent FTO, Apple Developer enrollment (agent 11).

---

**Sources fetched for this review:**
- [Twist Bioscience Gene Synthesis FAQ](https://www.twistbioscience.com/faq/gene-synthesis) — resolved the 50 bp GC-window contradiction (10–90%, not 35–65%)
- [IBBIS Common Mechanism](https://ibbis.bio/common-mechanism/) and [ibbis-screening/common-mechanism on GitHub](https://github.com/ibbis-screening/common-mechanism) — MIT-licensed, offline, protein-level biosecurity screening; cites the 2025 *Science* AI-variant evasion study
- [SignalP 6.0, DTU Health Tech](https://services.healthtech.dtu.dk/services/SignalP-6.0/) — academic-only download; commercial users must contact DTU
- [List of genetic codes](https://en.wikipedia.org/wiki/List_of_genetic_codes) — NCBI tables 2, 4, 6, 11, 12 reassignments
- [IGSC Harmonized Screening Protocol](https://genesynthesisconsortium.org/harmonized-screening-protocol/) — protocol PDF retrieved but technical thresholds not extractable; still needs manual review
- [FDA CMC guidance for human gene therapy INDs](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/chemistry-manufacturing-and-control-cmc-information-human-gene-therapy-investigational-new-drug) — landing page only; full guidance not parsed, claims above are flagged as unverified