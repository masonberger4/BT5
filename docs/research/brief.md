# Technical Brief: Local Codon-Optimization & Back-Translation Desktop App

**Scope:** decision-ready synthesis of 12 parallel research dossiers into (a) an honest assessment of what sequence-level design actually predicts, (b) a unified constraint/objective model, (c) a context matrix, (d) solver design, (e) data assets, (f) competitive gaps, (g) risks, (h) unresolved contradictions.

---

## 1. State of the art, honestly assessed

### 1.1 The single most important finding: CAI is not an expression objective

Every large, well-controlled synthetic-library experiment since 2009 converges on this. In Kudla 2009 (154 synonymous GFP variants in *E. coli*, 250-fold expression range), predicted folding free energy of the ~41-nt window spanning the start codon explained 44% of variance (r = 0.66; 59% / r = 0.77 in a second promoter system), while whole-mRNA MFE gave r = 0.16 (n.s.) and CAI gave r = 0.14 (n.s.) — https://pmc.ncbi.nlm.nih.gov/articles/PMC3902468/. Welch 2009 (ATUM/DNA2.0, 2×40 synonymous variants, >40-fold range, PLS R² = 0.77) states flatly that "CAI has no value in predicting gene expression," and their deliberately high-CAI control expressed at ~15% of the best variant — https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0007002. Ranaghan 2021 benchmarked nine commercial/academic optimizers and found "a roughly equivalent chance that an algorithm-optimized CDS will increase or diminish recombinant yields," with three tools non-deterministic (one returning 35–39% pairwise codon identity across ten identical submissions) — https://pmc.ncbi.nlm.nih.gov/articles/PMC7893858/.

Two 2026 wet-lab studies extend this to eukaryotes and invert the sign: a Pichia study found CAI **negatively** correlated with titer (coefficient −0.81 for trastuzumab) — https://europepmc.org/article/MED/41701818 — and an 18-glycoprotein / 90-screen Expi293F benchmark concluded "codon optimization to make human proteins, in a human cell line, did not generate increased yields," with native and harmonized constructs most consistent and structure-maximizing LinearDesign worst — https://proteininnovation.org/2026/03/codon-optimization-native-codon-mammalian-protein-expression/.

**Implication:** CAI belongs in the app as a *descriptive statistic and soft band*, never as a maximization target, and never labeled "predicted expression."

### 1.2 What does work: the 5′ region

- **5′ mRNA structure (bacteria).** The −4..+37 window MFE is the strongest single lever, explaining 44–59% of variance (Kudla). Cambray 2018 (244,000 designed sequences, full factorial) found phenotype "dominated by secondary structures and their interactions," localized to STR(−30:+30) and STR(+31:+90) — https://www.nature.com/articles/nbt.4238.
- **N-terminal composition.** Goodman/Church/Kosuri 2013 (>14,000 reporters) showed the N-terminal rare-codon benefit (up to 14-fold, median 4-fold) is *entirely* a reduced-structure effect, not codon rarity — https://europepmc.org/abstract/MED/24072823. A 2023 NAR study reached Pearson r = 0.762 predicting mRFP output from **codons 2–8 alone** — https://academic.oup.com/nar/article/51/5/2363/7016452.
- **Amino-acid identity at codons 3–5.** Verma 2019 (259,134 variants, 9,261 tripeptides) found three orders of magnitude of yield spread; K|N at position 3 + Y|I motif best (GFP score 4.31±0.87), T-V-G strongly detrimental; single-molecule processivity 84% (K-I-H) vs 27% (T-V-G) — https://pmc.ncbi.nlm.nih.gov/articles/PMC6920384/.
- **Kozak, mammalian.** Noderer 2014 measured all 65,536 −6..+5 variants: 12-fold range; −3 purine +58% over −3U; +4G/+5C 24.8% better than +4G/+5A, with cooperativity between +4 and +5 — https://pmc.ncbi.nlm.nih.gov/articles/PMC4299517/. Positions +4/+5 are codon 2, i.e. **designable by codon choice**.

### 1.3 What is folklore or overturned

| Claim | Status |
|---|---|
| Internal Shine–Dalgarno causes translational pausing | **Overturned.** Li/Weissman 2012 was an artifact of selecting 28–42 nt footprints (Mohammad, Green & Buskirk, https://elifesciences.org/articles/42591). Internal SD-driven *initiation* producing truncations is real; score with a TIR model, not a motif ban. |
| Codon-pair-bias deoptimization attenuates via codon pairs | **Largely a CpG/UpA artifact.** Min-E (deoptimized pairs, unchanged dinucleotides) had wild-type fitness; TE vs fitness r = −0.075 n.s. — https://elifesciences.org/articles/04531. Gutman & Hatfield reported over-represented pairs translate *slower* — opposite sign to Coleman 2008. |
| Slow 5′ "translational ramp" improves yield | **Reframed as a spandrel.** Slow-5′ construct gave 67–71% of the fast construct's GFP — https://elifesciences.org/articles/89656. |
| High GC lowers plasmid yield in *E. coli* | **Unsupported.** Documented yield/topology losses trace to non-B DNA (Z-DNA, triplex), toxicity/burden, and nuclease-sensitive tracts. |
| High GC kills AAV titer | **Contradicted by the only controlled dataset.** A computationally designed "inert" stuffer at GC 43.5–44.8% cut yield up to 68% and bioactivity 34–82%; a *lower*-GC (33.8–34.7%) natural UBE3A 3′UTR stuffer of identical length cost neither — https://pmc.ncbi.nlm.nih.gov/articles/PMC12207685/. The liability is repetitiveness/low complexity. |
| Chi-site removal helps mammalian/yeast expression | **Cargo cult.** Primary literature shows Chi is *protective* for plasmids (RecA-mediated survival). Keep as a cheap, low-priority *E. coli*/linear-DNA rule only. |
| rG4 avoidance helps mammalian translation | **Weak.** rG4s are globally unfolded in mammalian cells (median in-cell folding score 0.06) but fold and impair growth in *E. coli* — https://pmc.ncbi.nlm.nih.gov/articles/PMC5367264/. Keep G4 as a hard rule for bacterial CDS and DNA synthesizability, soft for mammalian mRNA. |

### 1.4 The ML picture

Codon language models (CodonTransformer 89.6M/Apache-2.0, CodonBERT 110M, SynCodonLM 102M, EnCodon/DeCodon up to 1B) now beat vendor tools on in-silico metrics and, in the best-controlled wet-lab test, beat both random synonymous variants and vendors: 83.8% of LM variants beat reference vs 5.9% of random; LM median specific productivity 1.20 vs 0.57 (vendor) vs 0.47 (random) — https://doi.org/10.64898/2026.08.11.744178. But the sober number is DeepCodon vs GenScript across 20 proteins: **9 better, 10 indistinguishable, 1 worse** — https://europepmc.org/article/MED/42038710. Prediction ceilings are well-characterized: RiboNN reaches r² = 0.62 for translation efficiency but its own LightGBM feature baseline reaches r = 0.78 vs the CNN's 0.79 — https://www.biorxiv.org/content/10.1101/2024.08.11.607362v2.full; binary high/low expression tops out at AUROC ~0.835; Optimus 5-Prime collapses from r² ≈ 0.93 in-distribution to 0.12 on native HEK293T transcripts. Cambray's ceiling is the honest headline: **all computable design features together explain 5–31% (mean ~14%) of protein-level variance.**

**Verdict:** ML belongs as an optional, pluggable *candidate generator* feeding a deterministic scoring/Pareto layer. Never as the mandatory core, never as a source of hard-constraint guarantees (no LM can guarantee absence of a BsaI site).

---

## 2. The objective/constraint model

Notation: **H** = hard constraint (must reach score 0 or design is refused), **S** = soft objective (weighted). Evidence: **A** = large controlled dataset or physical necessity; **B** = replicated but contested or single-lab; **C** = mechanistic plausibility / vendor assertion only.

### 2.A Correctness invariants (H, evidence A, all contexts)

| ID | Rule |
|---|---|
| A1 | Translated CDS ≡ input protein exactly. Enforce structurally via the mutation space (synonymous codon sets), not by penalty. |
| A2 | Length ≡ 0 mod 3; no internal in-frame stop; exactly one terminal stop unless a C-terminal fusion is requested (then stop removed). |
| A3 | Alphabet ∈ {A,C,G,T}. Reject IUPAC ambiguity in output. |
| A4 | **Cassette frame invariant:** `len(5′ overhang/homology contribution + Kozak spacer + signal peptide + N-tag + N-linker) mod 3 == 0` and `len(C-linker + C-tag) mod 3 == 0`. Signal peptide must be the absolute N-terminus; N-tags go *after* the predicted signal-peptidase site. |

### 2.B Translation initiation / 5′ region

| ID | Type | Rule | Context | Ev |
|---|---|---|---|---|
| B1 | S (highest weight) | Maximize (make less negative) ΔG of window **−4 … +37** relative to the A of ATG. Requires splicing in the real vector 5′UTR. | Bacteria | A |
| B2 | S | Penalize structure in **STR(−30:+30)**; neutral-to-mildly-reward moderate structure in **STR(+31:+90)**. | Bacteria | A |
| B3 | H | Reject if ΔG(5′UTR + codons 1–16) < **−39 kcal/mol** *AND* GC(codons 2–6) > **62%**. Both together (Boël 2016, 6,348 genes) — https://pmc.ncbi.nlm.nih.gov/articles/PMC5054687/ | Bacteria | B |
| B4 | S | Codons 2–6: favor A, avoid G. In the first ~11 codons prefer A/T-ending codons (22/30 A/T-ending have positive initiation log-odds; 22/29 G/C-ending negative). | Bacteria | A |
| B5 | S (opt-in) | Codons 3–5 amino-acid motif: prefer K\|N at position 3 with Y\|I in 3–5; hard-avoid Thr-Val-Gly at 3–5. Only when the N-terminus is designable (tags/linkers). | Bacteria, PURE | A |
| B6 | S | **RBS/TIR model** (Salis): `r ∝ exp(−β·ΔG_tot)`, β = 0.45 ± 0.05 mol/kcal; ΔG_tot = ΔG_mRNA:rRNA + ΔG_start + ΔG_spacing − ΔG_standby − ΔG_mRNA over a 70-nt window (±35 of start). 16S anti-SD = 3′-AUUCCUCCA-5′. ΔG_start = −1.19 (AUG), −0.075 (GUG). Optimal spacer 5 nt. Median ~2.3-fold error — https://pmc.ncbi.nlm.nih.gov/articles/PMC2782888/ | Bacteria | A |
| B7 | H/S | **Internal initiation:** score every internal ATG/GTG/TTG with the TIR model; flag any internal TIR > 10% of the intended start's. Use AGGAGG/GGAGG only as a pre-filter. Do **not** penalize internal SD as a "pausing" signal. | Bacteria | B |
| B8 | S | **Kozak:** target `GCCRCCATGG` (equivalently gccRccATGG). Set +4 = G, +5 = C by codon-2 choice where the residue permits. Score strong (−3 purine AND +4 G) / adequate / weak. | Eukaryote | A |
| B9 | H | No additional ATG in the first 50 nt of CDS; penalize any out-of-frame ATG anywhere that has BOTH −3 purine and +4 G. Ensure optimization does not create an upstream out-of-frame AUG at the UTR/CDS junction. | Eukaryote | A |
| B10 | H (pre-flight, non-editable) | Scan user's 5′UTR for out-of-frame uAUGs. uORFs occur in ~half of human transcripts and typically cut protein 30–80% — https://pmc.ncbi.nlm.nih.gov/articles/PMC7100133/. Warn; cannot be fixed by codon choice. | Eukaryote | A |
| B11 | H | Cap-proximal structure: any hairpin ΔG ≤ −30 kcal/mol with 5′ base within 15 nt of cap = disqualifying; ≤ −50 kcal/mol within 15–60 nt; ≤ −60 kcal/mol beyond ~70 nt (Kozak 1986/1989 — magnitudes approximate, position-dependence solid). | Eukaryote / IVT | B |

### 2.C Codon composition (all S, soft bands, never maximized)

| ID | Rule | Ev |
|---|---|---|
| C1 | **CAI** (Sharp & Li): w_i = count_i / max_synonymous_count computed from a *highly-expressed* reference set (not the genomic table); pseudocount 0.5 before the family max; CAI = exp(mean ln w), excluding ATG/TGG and stops. Target a **band** (e.g. 0.70–0.90, or ±0.1 of host median). Never 1.0. | A (that it's weak) |
| C2 | **tAI / stAI:** W_i = Σ_j (1−s_ij)·tGCN_j; w = W/max(W); zeros → geometric mean of non-zeros. dos Reis s = (0,0,0,0, G:U 0.41, I:C 0.28, I:A 0.9999, U:G 0.68, lysidine C:A 0.89 — **prokaryotes only**). Tuller variant changes G:U to 0.561. Use species-specific s-values (stAIcalc) outside yeast: yeast vector has only median r ≈ 0.5 with species-optimized vectors. | B |
| C3 | **%MinMax:** sliding window z = 17–18 codons (CHARMING default 10). %Max = 100·Σ(X_ij − X_avg,i)/Σ(X_max,i − X_avg,i) for above-average windows; %Min analogous below. Shrink the window at termini rather than dropping positions. | B |
| C4 | **Harmonization objective:** minimize L2 or DTW distance between the design's %MinMax profile in the host and the native sequence's profile in *its* host. Do **not** minimize deviation from a flat max. Per-codon equivalent = DnaChisel `harmonize_rca`. | B |
| C5 | **CFD:** fraction of codons with host relative synonymous frequency < 0.30. Minimize. | B |
| C6 | **ENc/Nc** (Wright): range [20, 61]. Descriptive only — never an objective. | A |
| C7 | **CSC (eukaryote mRNA stability):** per-codon Codon Stability Coefficients (Pearson r between half-life and codon occurrence) as a *separate* objective from elongation terms — https://elifesciences.org/articles/45396 | A |
| C8 | **Rare-codon cluster preservation:** retain ≥80% of native rare-codon clusters when a native CDS is supplied (DeepCodon 80–90%; ICOR/generic <50%). Detect via %MinMax %Min excursions or window-less Maximal Scoring Subsequence. | B |
| C9 | **CPB/CPS:** CPS(AB) = ln[ O(AB) / ((F(A)F(B)/(F(aaA)F(aaB)))·O(aaA-aaB)) ]. **Near-zero weight by default**; expose only for deliberate viral deoptimization, and always paired with explicit CpG/UpA reporting. | C (contested sign) |
| C10 | **Out-of-frame stop density (novel):** maximize stop-codon density in the −1 and +1 frames; constrain max out-of-frame ORF ≤ ~60 codons. Natural human CDSs keep out-of-frame stops at median ~20-codon spacing; 120 therapeutic sequences (incl. approved COVID mRNA vaccines) average 164 aa — https://www.pnas.org/doi/10.1073/pnas.2606609123. Achieved by retaining some T at wobble positions; naturally opposes over-optimization. | B |

### 2.D Motif avoidance — gated by host and application

**Expected-occurrence budget** (encode this and refuse ≤5-mers as hard constraints without override): a non-palindromic k-mer on both strands of length L occurs ≈ 2(L−k+1)/4^k times. For L = 1500: 5-mer ≈ 2.9, 6-mer ≈ 0.73, 7-mer ≈ 0.18, 8-mer ≈ 0.046, 25-mer ≈ 3e−9. **Rule of thumb:** ≥7 nt → hard (cheap); 6-mer → hard only if strongly evidenced; ≤5-mer → soft only.

**D1 Restriction / Type IIS (H, scan BOTH strands — these are non-palindromic):**
BsaI `GGTCTC`(1/5)/`GAGACC`; BsmBI/Esp3I `CGTCTC`(1/5)/`GAGACG`; BbsI/BpiI `GAAGAC`(2/6)/`GTCTTC`; SapI/LguI `GCTCTTC`(1/4, **3-nt** overhang)/`GAAGAGC`; AarI/PaqCI `CACCTGC`(4/8)/`GCAGGTG`; BfuAI/BspMI `ACCTGC`(4/8); BtgZI `GCGATG`(10/14)/`CATCGC`; BsmAI `GTCTC`/`GAGAC`.
Six-cutters: EcoRI GAATTC, BamHI GGATCC, HindIII AAGCTT, XhoI CTCGAG, SalI GTCGAC, XbaI TCTAGA, SpeI ACTAGT, PstI CTGCAG, KpnI GGTACC, SacI GAGCTC, NheI GCTAGC, BglII AGATCT, EcoRV GATATC, NdeI CATATG, NcoI CCATGG, AgeI ACCGGT, MluI ACGCGT, SmaI/XmaI CCCGGG, ApaI GGGCCC, ClaI ATCGAT, AflII CTTAAG, AvrII CCTAGG, BspEI TCCGGA, BsrGI TGTACA, PvuII CAGCTG, StuI AGGCCT, SnaBI TACGTA.
Eight-cutters: NotI GCGGCCGC, AscI GGCGCGCC, PacI TTAATTAA, FseI GGCCGGCC, SbfI CCTGCAGG, PmeI GTTTAAAC, SwaI ATTTAAAT, AsiSI GCGATCGC, SgrAI CRCCGGYG.
**Surface the NcoI CCATGG ⊂ Kozak GCCACCATGG conflict rather than silently breaking one.**

**D2 Recombinase/recombination sites (H, check-only — 25–48 bp, never arise by chance):**
loxP family: regex `ATAACTTCGTATA[ACGT]{8}TATACGAAGTTAT` + revcomp (covers loxP, lox2272, lox5171, loxN, lox511); half-sites `ATAACTTCGTATA` / `TATACGAAGTTAT`. FRT: `GAAGTTCCTATTC[ACGT]{8}GTATAGGAACTTC` + revcomp. Gateway attB1 `ACAAGTTTGTACAAAAAAGCAGGCT`, attB2 `ACCCAGCTTTCTTGTACAAAGTGGT`, shared core `TTTGTACAAA[AG]`. **Bxb1 attB and attP both contain `GGTCTC` (BsaI)** — flag this collision explicitly when a user wants both a landing pad and BsaI-free Golden Gate.

**D3 Splicing (H/S, eukaryotic Pol II contexts only — irrelevant for *E. coli*, yeast heterologous CDS, IVT mRNA):**
5′ donor: score every GT with MaxEntScan 9-mer (3 exonic + 6 intronic); flag >3 bits, hard-constrain >6–8. Literal blacklist: `GGTAAG`, `GGTGAG`, pattern `AN|GT(A/G)AG`, coarse `GTNNG`.
3′ acceptor: score every AG with MaxEntScan 23-mer (20 intronic + 3 exonic); flag >3 bits. Also flag any AG preceded within 5–40 nt by a ≥10-nt window ≥80% pyrimidine and a branch-point-like `YTNAY` 18–40 nt upstream.
**Iterate to a fixed point:** point-mutating one donor activates cryptic donors nearby (A2UCOE, https://jvi.asm.org/content/86/9088). A single-pass "remove motifs" step is unsafe.
**V5 tag special case (H, evidence A):** the standard V5 nucleotide encoding contains `G|GTAAG` and spliced in **17/17** genes tested; 13/17 randomly chosen genes showed aberrant splicing from vector/tag context — https://pmc.ncbi.nlm.nih.gov/articles/PMC9379414/. Recode V5 (and all tags/linkers) to destroy donors. This is a concrete win only a back-translating tool can deliver.

**D4 Polyadenylation (H for lentiviral sense strand, S elsewhere):**
Hard: `AATAAA`, `ATTAAA` (~61.6% and ~15% of dominant hexamers). Soft: AGTAAA, TATAAA, CATAAA, GATAAA, AATATA, AATACA, AATAGA, ACTAAA, AAGAAA, AATGAA. **Escalate to CRITICAL** when a downstream element exists within +10..+40 nt: a GU-rich window containing GTGT/TGTG, or a U-run ≥4 T in a 6-nt window. Optionally score with APARENT2 (MIT).

**D5 Bacterial cryptic transcription (H/S — *E. coli* propagation, applies to ALL constructs regardless of final host):**
(a) hexamer within 1 mismatch of `TTGACA` + 15–19 bp spacer + hexamer within 1 mismatch of `TATAAT`; (b) extended −10 `TGnTATAAT` (needs no −35); (c) **AT-tract rule**: `TATTTAT` or `AATTT` whose 3′ end sits exactly 5 bp upstream of a −10-like hexamer — 103/103 randomized appropriately positioned AT-tracts were active — https://academic.oup.com/nar/article/48/9/4891/5820884; (d) σ38 variant `TATACT`. Score quantitatively with the Salis Promoter Calculator (346 parameters, validated on 17,396 promoters). **Calibration anchor:** the dengue-2 cryptic promoter (−35 `TCAACG` at nt 53, 17-bp spacer, −10 `TTTTTAAT` at nt 72) produced ~10^6.6 mRNA copies/µg uninduced and made the clone uncloneable — https://pmc.ncbi.nlm.nih.gov/articles/PMC3069047/.
**Antisense promoter (novel, no competitor implements):** forbid `ATTATA` (revcomp of TATAAT) on the sense strand. Only 4.76% of 484,741 natural *E. coli* CDSs contain it (2–3× under-represented), yet 77.28% of clean CDSs can silently acquire it by synonymous substitution — https://pmc.ncbi.nlm.nih.gov/articles/PMC13029128/.

**D6 Non-B DNA (H):**
G-quadruplex regex `G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}` both strands; plus G4Hunter (window 25, |score| ≥ 1.2 flag, ≥1.5 severe). *E. coli* mutation rates span 5.5e−5 to 2.7e−10 per cell per generation across G4 variants, with up to 8-fold orientation dependence relative to the fork — https://pmc.ncbi.nlm.nih.gov/articles/PMC10530614/. Note the regex misses ~37% of experimentally detected rG4s.
Z-DNA: alternating purine-pyrimidine, (CG)n ≥ 6 units or (CA/TG)n ≥ 7 units. Ranked *most* destabilizing non-B structure in a cloning vector (single 2004 study — evidence C).
Triplex/H-DNA: homopurine-homopyrimidine mirror repeats ≥ 20 bp.
Telomere repeats: `(TTAGGG){2,}` vertebrate, `(TTTAGGG){2,}` plant, `(TTGGGG){2,}`, `(TTTTGGGG){2,}`.

**D7 Other RNA elements (S):**
ARE — implement **AREScore verbatim**: +1 per `ATTTA`; +1.5 overlapping; +0.75 if 0–3 nt apart; +0.4 if 4–6; +0.2 if 7–9; +0.3 inside an AU-block (opens when a 20-nt window is ≥80% A+T, closes below 55%) — https://journals.plos.org/plosgenetics/article?id=10.1371%2Fjournal.pgen.1002433. Hard-avoid the class-II nonamer `TTATTTATT` / `TTATTTA[TA][TA]`; isolated `ATTTA` is soft only (expected ~2.9 per 1.5 kb).
INS/CRS: **no regex exists** ("no uniformly recognizable sequence commonality"). Implement as a composite proxy (AREScore + 50-nt A+T > 65% + low GC3) and label honestly.
−1 frameshift: `([ACGT])\1\1(AAA|TTT)[ACT]` anchored to codon boundaries; escalate only if a hairpin ΔG ≤ −10 kcal/mol begins 5–9 nt 3′.
Pol III terminator: hard-avoid `TTTT` in any Pol III-transcribed region (T4 minimal; ≥T6 full efficiency). Check antisense `AAAA` for bidirectional contexts.
Chi: `GCTGGTGG`/`CCACCAGC` — *E. coli*/linear-DNA/recombineering only, low priority.
Dam `GATC` / Dcm `CCWGG` — informational, flag when overlapping a chosen restriction site.

**D8 CpG — three separate, separately-toggleable metrics (do not collapse into one slider):**
(a) **TLR9 DNA sensing:** total CpG count + stimulatory hexamers `RRCGYY`, specifically `GTCGTT`, `TTCGTT` (human), `GACGTT` (mouse). CpG-depleted AAVrh32.33 evaded immune detection — https://www.jci.org/articles/view/68205.
(b) **ZAP/KHNYN RNA decay (novel, no competitor exposes this):** sliding 200-nt window; flag any window with ≥14 CpGs at mean inter-CpG spacing ≤14 nt (peak sensitivity 6–14 nt). Spacing ≥32 nt is **not** restricted, and "the magnitude of ZAP-mediated inhibition was not correlated with the number of CpGs introduced" — https://pmc.ncbi.nlm.nih.gov/articles/PMC9519448/. Report the worst window, never the global count.
(c) **Methylation silencing:** CpG islands — Gardiner-Garden (≥200 bp, GC ≥50%, obs/exp ≥0.6) or Takai-Jones strict (≥500 bp, GC ≥55%, ≥0.65). obs/exp = (N_CpG × L)/(N_C × N_G).
**Mechanics:** CpG arises both within codons (CGN Arg, GCG, CCG, TCG, ACG) *and* at codon junctions (codon ends C, next begins G). Evaluate the dinucleotide across boundaries. Warn that full depletion forces AGA/AGG for Arg and can drop GC below vendor floors.

### 2.E Manufacturability (H/S, all contexts, evaluated on the synthesized fragment + adapters)

| ID | Rule | Ev |
|---|---|---|
| E1 | Homopolymers: target ≤9 nt (A/T) and ≤5 nt (G/C). Hard-fail ≥14 (Twist Standard), >15 (GenScript), ≥10 A/T or ≥6 G/C (IDT gBlocks). Twist Complex Genes accepts to 30 bp at a surcharge. The 6-vs-10 asymmetry is the clearest vendor evidence that G/C runs are chemically worse (G-quadruplex aggregation on solid support). | A |
| E2 | Global GC: optimize into 40–60%; warn outside 25–65%; hard-fail outside 25–75%. | A |
| E3 | Windowed GC 50 bp: hard-fail any window <10% or >90% (Twist High-Complexity trigger); warn <25% or >75%. Windowed GC 100 bp: warn outside 25–65% (GenScript GenTitan — the only vendor publishing a windowed rule) — https://www.genscript.com/gentitan-gene-fragments.html | A |
| E4 | **GC variation.** ~~max(GC_50bp) − min(GC_50bp) ≤ 50, target ≤25~~ — **corrected 2026-08-28, these thresholds are below the chance floor.** Random 50% GC DNA has a 50 bp extent of 26.0 at 300 bp, 36.0 at 1.2 kb, 44.0 at 5 kb and 46.0 at 10 kb (200 draws/length), so the ≤25 target fires on every sequence at every length ≥300 bp and ≤50 sits ~4 points above chance at 10 kb. The cause is structural: extent is a RANGE statistic and widens with every window added, whereas `dGC` (SD of GC over 100 bp windows) is a DISPERSION statistic and converges on the binomial floor 100·√(p(1−p)/w) = 5.0 at p=0.5, w=100. So score **`dGC` relative to its own binomial floor** (1.0 = chance; the ratio cancels the composition confound, which matters because synonymous choice alone spans 27–60% achievable GC for one 300 aa protein), and keep extent only to LOCALISE the finding to the two windows worth recoding. Calibration: an ordinary protein reads 1.01, a repetitive one 1.34. Note the 50 bp window is Twist's published window SIZE — Twist publishes no bound on extent or SD at all, which is why this is SOFT. `dGC` outranked global GC in SCP4ssd, but the Synthesis Success Calculator ranks repeats above GC, so the evidence is **contested**, not A. Implemented as `e4_gc_extent`. | B/contested |
| E5 | Direct repeats: hard-fail ≥20 bp (Twist and GenScript both publish 20); warn >12–16 bp; hard-fail >200 bp. Flag repeat clusters at 8–9 bp (Twist names this size). Flag any repeated unit with duplex Tm ≥ 60 °C (catches short-but-GC-rich repeats a length rule misses). | A |
| E6 | **Repeat density:** count repetitive 9-mers per 100 bp, and longest repetitive sequence — the two highest-Gini features of the Synthesis Success Calculator (random forest, 1076 real vendor outcomes, F1 0.928) — https://pubs.acs.org/doi/10.1021/acssynbio.9b00460. **Weight repeats above GC.** | A |
| E7 | STRs (unit 1–6 bp): hard-fail total tract >100 bp; warn >20 bp. Inverted repeats/hairpins: hard-fail tract >100 bp; warn on any 20-bp stem whose revcomp occurs within 200 nt (the IDT rule DnaChisel encodes as `AvoidHairpins(stem_size=20, hairpin_window=200)`). | A/B |
| E8 | **k-mer uniqueness:** every 12-mer (ideally 15-mer) unique across the *whole construct including backbone and vendor adapters*. Single most effective repeat rule — simultaneously protects assembly PCR, Gibson overlaps, and plasmid stability. Twist Gene Fragment adapters: 5′ `CAATCCGCCCTCACTACAACCG`, 3′ `CTACTCTGGCGTCGATGAGGGA`. | A |
| E9 | Length tiers: ≥300 bp minimum (Twist); ≤1,500 IDT eBlocks; ≤3,000 GenScript GenTitan; ≤5,000 Twist Gene Fragments; ≤7,000 Twist Clonal. | A |
| E10 | His-tag: encode 6×His by **alternating CAC/CAT** (`CACCATCACCATCACCAT`). Twist calls this out by name. Generalize: when a homopolymer run is forced by the protein (poly-Lys, poly-Gly), alternate synonymous codons. | A |
| E11 | Screening burden: P(perfect clone) ≈ exp(−L/E) with E = 7,500 bp (Twist) or 5,000 bp (IDT eBlocks, GenScript). Report colonies to pick for 95% confidence. | A |
| E12 | Cost/tier: Gene Fragments ~7¢/base; Clonal ~9¢/base; Complex Clonal 12¢/bp (0.3–1.8 kb) → ~23¢/bp (3.2–7.0 kb). Surface the price consequence of each violation. | B |

### 2.F Plasmid propagation in *E. coli* (H/S — applies to **every** construct that passes through a cloning host)

| ID | Rule | Ev |
|---|---|---|
| F1 | **Direct repeats — hard fail** any perfect direct repeat ≥25 bp anywhere in the assembled plasmid, both strands (RecBCD MEPS = 23–27 bp; RecF pathway 44–90 bp). Warn 20–24 bp; info 15–19. A measured 28-bp repeat pair recombined at 7.8e−7 to 3.1e−5 in **four different recA⁻ strains** — https://www.genoscope.cns.fr/MGE/pubs/Oliveira_Mol_Biotechnol_2008.pdf. recA⁻ is a mitigation, not a fix: below ~200 bp, deletion is RecA-*independent* (slipped-strand/SSA), unaffected by recA/recF/recJ/recO. | A |
| F2 | Near-perfect: flag pairs ≥90% identity over ≥40 bp. Score risk monotone in repeat length, exponentially decaying in spacer length; spacer <3 kb = high risk. | B |
| F3 | **Inverted repeats — hard fail** perfect stem ≥30 bp or total palindrome ≥60 bp; "do not build in a standard host" at perfect palindrome ≥150 bp (SbcCD destroys the replicon). Warn stem 15–29 bp with loop ≤100 nt. Rank by nearest-neighbour ΔG of the **basal 20 bp of the stem** (Sinden's determinant), flag ≤ −20 kcal/mol. IR <100 bp usually stable; 100–1200 bp transiently unstable. | A |
| F4 | If an IR is unavoidable (shRNA, AAV ITR): place ≥438 bp from the origin (ori-proximal ITRs degrade preferentially — *provisional number*), and emit a ΔsbcC-host-at-42 °C protocol recommendation. | B |
| F5 | **AT-window rule:** every 100-nt window ≤55% AT (≥45% GC); none above 60% AT. Toxic horizontally-acquired *E. coli* genes are 63–68% AT vs a non-toxic 55% control — https://www.nature.com/articles/nmicrobiol2016249. **This directly conflicts with vendor GC ceilings** — resolve as a two-sided band 45–60% GC per 100 nt, hard-fail outside 35–65%, and show which side is binding per window. | A |
| F6 | Homopolymers ≤8 nt; dinucleotide repeats ≤5 units; trinucleotide repeats ≤5 units. | B |
| F7 | **Toxicity heuristic:** classify the encoded protein. Predicted TM segments → toxicity probability: 0 TM ≈ 25%, 1 TM ≈ 73%, ≥2 TM ≈ 85–89% (ASKA library, 1589/3956 clones toxic overall) — https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0064893. On trigger, recommend a low-copy backbone. | A |
| F8 | **Copy-number multiplier:** pUC/pMB1* 300–700; pBR322 (rop+) 15–25 / (rop−) 45–60; ColE1 15–20; p15A ~10; pSC101 ~5. Detect the ori from the backbone sequence and use it as a risk multiplier on F1 and F7. | A |
| F9 | **Protocol output** (not sequence design): repeat-tolerant strain (Stbl3 recA13, NEB Stable recA1), 30 °C, ampicillin 50 not 100 µg/mL, harvest late-log, **pick small colonies** (fast growers enrich for deletions), screen ≥4 colonies by diagnostic digest. For IR/ITR: ΔsbcC at 42 °C instead. Stbl2→Stbl3 alone rescued an HIV vector that was lost entirely in 0.5 L Stbl2 cultures — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3563744/ | A |

### 2.G Vector / delivery-specific

**AAV (H):** ssAAV total ITR-to-ITR target 4.0–4.4 kb; warn >4.7; hard-fail >5.0 (nothing above 5.2 kb ever packages — https://pubmed.ncbi.nlm.nih.gov/19904234/). Optimal full-capsid yield window **2.5–3.5 kb**. Yield: 3.5e13 vg at 2.0 kb → 1.4e13 at 4.5 kb (−60%) → 7.0e12 at 5.0 kb (−80%); partial capsids 46–59% at 5.0 kb; **over-filled 25–47% below 2.25 kb** (recommend a stuffer). scAAV: ≤2.4 kb; at ≥3.0 kb 100% of genomes were single-stranded. **Palindrome/IR rule is AAV-critical:** reject IRs with stem ≥20 bp and loop ≤200 nt; flag hairpins ΔG ≤ −20 kcal/mol (engineering guess — no published threshold). AAV-GPseq mapped truncation hotspots exactly to inverted repeats in the CMV enhancer, CB promoter and EGFP ORF — https://www.cell.com/molecular-therapy-family/advances/fulltext/S2329-0501(20)30156-X. **Do not raise GC to "help packaging"** — the only controlled dataset points the other way. Space savings: SV40 late polyA 135 bp vs bGH 225 bp; NRP1 polyA 32 bp; WPRE3/247-bp WPRE vs 600 bp.

**Lentiviral transfer (H):** determine which strand becomes the packaged genome (forward-oriented cassette → sense strand; reverse-oriented β-globin-style → revcomp) and scan **that** strand. Hard-fail sense-strand `AATAAA`/`ATTAAA` — internal polyA raised expression 3–6.5× but cut functional titer **8–9×** with CMV or EF1α (and not at all with β-actin or PSA/Pb) — https://pubmed.ncbi.nlm.nih.gov/18627247/. Hard-fail strong sense-strand splice donors. Size: ≤~8.5 kb LTR-to-LTR; ~3–4× functional titer loss per additional kb. Warn on IRs >20 bp. Offer the forward-oriented+RRE-in-intron design (6× titer, 4–10× LTRC transduction) when the cassette contains a retained intron — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6775231/.

**Adenovirus (H):** genome ≥~27 kb (below 75% of 35.9 kb it rearranges) and ≤37.7 kb (105%); fail >38.9 kb.

**Gamma-retrovirus (H):** ≤~8.2 kb total.

**Producer-cell toxicity (S, high-value, invisible to any pure sequence scan):** if the protein is pro-apoptotic, a cell-cycle regulator, strong TF, membrane channel, protease, fusogen, or secreted clotting factor, warn that producer-cell expression likely dominates titer and recommend (a) tissue-specific promoter, (b) TRiP/TRAP translational repression (tbs = [KAGNN]₁₁ near the initiation codon), or (c) Tet-off. Magnitudes: COX-2 vectors titered ~1000× below GFP; TRAP rescued 600× (EIAV), 30× (CAR), **>150,000×** (adenoviral Bax) — https://pmc.ncbi.nlm.nih.gov/articles/PMC5378976/. Codon-optimized FVIII gave ~10× lentiviral titer, but the mechanism was protein-level (FVIII blocks VSV-G incorporation), not RNA-level.

**IVT mRNA (S, sign-inverted):** after ~codon 15, **minimize** MFE. LinearDesign spike: ~−1500 vs ~−900 kcal/mol, half-life 20.0 h vs 3.9 h (5.1×), 2.3–2.9× protein, 57–128× IgG — https://pmc.ncbi.nlm.nih.gov/articles/PMC10499610/. But keep the 5′UTR and first 30 CDS nt unstructured (LinearDesign excluded the first 5 codons and hand-picked the first 15 nt), and forbid perfect duplexes ≥33 bp (PKR/RIG-I/MDA5). Degradation objective: penalize unpaired uridines — DegScore-style per-linkage scoring (3′ nucleotide identity G≈C < A < U, weighted by unpaired probability) reaches Spearman −0.66 vs −0.50 for MFE — https://pmc.ncbi.nlm.nih.gov/articles/PMC8940940/.

### 2.H Assembly / backbone integration (H)

**Golden Gate:** zero internal recognition sites for the level enzyme (default domestication set: BsaI, BsmBI, BbsI, SapI). Overhang selection: no palindromes (AATT, GATC, TTAA, CATG, GTAC, AGCT, TCGA, ACGT, CCGG, GGCC); no 0% or 100% GC; ≥2 nt difference between any two overhangs *and* between each overhang and the revcomp of every other; cap at 20 four-base junctions (near-perfect fidelity) or 10 three-base SapI junctions (>99%); reject tatapov relative self-annealing <0.4 (hard <0.05) or cross-annealing >0.08. Benchmarks: 35-fragment BsmBI = 71% correct colonies; 13-fragment SapI = 91% — https://pmc.ncbi.nlm.nih.gov/articles/PMC7467295/.
**MoClo/Phytobrick fusion sites (encode verbatim, and note their frame semantics):** A=GGAG, B=TACT, C=AATG, D=AGGT, E=GCTT, F=CGCT, G=TGCC, H=ACTA. `AATG` = A + ATG (supplies the start codon); `AGGT` after a signal peptide reads ...A|GGT = Gly; `TTCG` at a CDS→C-tag junction reads ...T|TCG = Ser. CDS-no-stop must end in frame immediately before TTCG.
**Overlap methods:** Gibson 20–40 bp (default 25), NEBuilder HiFi 15–30, In-Fusion exactly 15 bp. Arms at Tm 48–60 °C, GC 40–60%, not ending in ≥3 identical bases. **Hard check: no exact repeat ≥20 bp shared between insert and any other part of the backbone.**
**Gateway:** ccdB+ CmR destination vector needs a ccdB-survival strain; attB scars add ~8–9 residues in frame. **Do not hard-code attB sequences from secondary sources** — verify against Hartley 2000 and a real pDONR221 GenBank.
**TOPO directional:** insert must begin exactly `CACC-ATG…` (vector GTGG overhang); no extra bases between CACC and ATG.

### 2.I Cassette element library (all back-translated by the same optimizer — the app's differentiator)

2A peptides (aa, incl. GSG): P2A `GSGATNFSLLKQAGDVEENPGP`; T2A `GSGEGRGSLLTCGDVEENPGP`; E2A `GSGQCTNYALLKLAGDVESNPGP`; F2A `GSGVKQTLNFDLLKLAGDVESNPGP`. Skip occurs between the final G and P. Default P2A/T2A; **warn on F2A** (up to ~50% uncleaved). Downstream gene expresses at 5–30% of upstream — put the high-need gene first. **Multiple 2As in one transcript must be different peptides**, and the optimizer must enforce ≤85% nucleotide identity between any two 2A cassettes (identical repeats hamper expression and destabilize lentiviral constructs) — https://pmc.ncbi.nlm.nih.gov/articles/PMC5438344/.
Tags: His6 `HHHHHH`; FLAG `DYKDDDDK`; 3×FLAG `DYKDHDGDYKDHDIDYKDDDDK`; HA `YPYDVPDYA`; c-Myc `EQKLISEEDL`; V5 `GKPIPNPLLGLDST`; Strep-II `WSHPQFEK`; Twin-Strep `SAWSHPQFEKGGGSGGGSGGSAWSHPQFEK`; AviTag `GLNDIFEAQKIEWHE`; ALFA `SRLEEELRRRLTE`; SBP, CBP, TC-tag `CCPGCC`.
Proteases: TEV `ENLYFQ↓G/S`; HRV-3C/PreScission `LEVLFQ↓GP`; Thrombin `LVPR↓GS`; Factor Xa `IEGR↓`; enterokinase `DDDDK↓`; SUMO/Ulp1 cuts after ...GG↓ (native N-terminus, any residue but Pro).
Linkers: `(GGGGS)n`, `(Gly)6/8`, `GSAGSAAGSGEF`, `KESGSVSSEQLAQFRSLD`, XTEN-16 `SGSETPGTSESATPES`; rigid `A(EAAAK)nA` (n=2–5), `A(EAAAK)4ALEA(EAAAK)4A`, `PAPAP`. **Force ≥2 nt differences between consecutive (GGGGS) repeat units** — naive back-translation creates a perfect nucleotide repeat that fails synthesis and recombines.
Stop codons: default `TAA`, emit tandem `TAATAA`. Hard-avoid `TGAC…` and `TGACTAG`. Readthrough TGA > TAG > TAA; basal human UGA contexts 0.9–7.8%, up to 31.3% — https://pmc.ncbi.nlm.nih.gov/articles/PMC12233894/.

---

## 3. Context matrix

### 3.1 By expression host

| Rule family | *E. coli* | Yeast/Pichia | CHO/HEK | Insect (Sf9/Tni) | Plant |
|---|---|---|---|---|---|
| 5′ structure (−4..+37) | **Primary objective** | Moderate | Cap-proximal rules instead | Moderate | Moderate |
| RBS/TIR (B6) | Yes | No | No | No | No |
| Kozak (B8/B9) | No (SD instead) | Yes (`aAaAaAATGTCt`) | **Yes, +4/+5 designable** | Yes | Yes (`acAACAAATGGC`) |
| CAI/tAI weight | Low-moderate (strong bias host) | Moderate | **Very low** (isochore GC, not selection) | Very low | Low |
| Default mode | Structure-first | Structure-first / harmonize | **Native or harmonize** | Native/harmonize | Harmonize |
| Splice/polyA scan | Off | Off (heterologous CDS) | **On** | On | On |
| Cryptic bacterial promoter | On (host) | On (propagation only) | On (propagation only) | On (propagation) | On (propagation) |
| G4 | **Hard** | Hard | Soft | Soft | Soft |
| CSC/mRNA stability | No | No | **Yes, separate objective** | No | No |
| Rare-codon strain gate | pRARE (AGA, AGG, ATA, CTA, GGA, CCC) / pRARE2 (+CGG) — disable penalty if Rosetta declared | — | — | — | — |
| Compartment split | — | — | — | — | **Nuclear vs chloroplast tables are different objects** (C. reinhardtii nuclear GC ~66% vs chloroplast GC3 20.6±7.6%) |

Attach a per-host **evidence-strength badge** derived from genome-median ENc and CAI-vs-ribosomal-protein correlation: high for *E. coli*, *B. subtilis*, *S. cerevisiae*, *K. phaffii*, *C. reinhardtii* nuclear; low for human, mouse, CHO, Sf9, Tni.

### 3.2 By delivery context

| Objective | Plasmid transient | Plasmid stable | Lentiviral | AAV | IVT mRNA | Genome-integrated |
|---|---|---|---|---|---|---|
| Size hard cap | ~10 kb warn | ~10 kb | **8.5 kb** | **4.7 kb ss / 2.4 kb sc** | n/a | HDR arms 500–1000 bp |
| Internal polyA | soft | soft | **HARD FAIL (sense)** | soft (expression, not titer) | n/a | soft |
| Cryptic splice donor | warn | warn | **HARD (titer + safety)** | warn | n/a | **HARD (fusion transcripts)** |
| CDS structure sign | neutral | neutral | neutral | neutral (but no IRs) | **minimize MFE past codon 15** | neutral |
| CpG-ZAP window | off | off | **on** | on | on | off |
| CpG-TLR9 | off | off | on | **on** | on | off |
| CpG methylation/silencing | off | **on** (CMV silences) | on | on | n/a | **on** |
| Inverted repeats | F3 | F3 | F3 + RT recombination | **stem ≥20 bp = hard (surrogate ITR)** | duplex ≥33 bp | F3 |
| Out-of-frame stop density | soft | soft | soft | soft | **on** | soft |
| Producer-cell toxicity | n/a | n/a | **on** | **on** | n/a | n/a |
| Promoter/host silencing warn | — | CMV in iPSC/CHO → EF1α/CAG | same | same | — | same |

### 3.3 The user's compound case: plasmid → virus → transduce → express

This is the case that justifies the whole product, and it means **three constraint sets apply simultaneously to one sequence**:

1. **E. coli propagation (transfer plasmid + packaging plasmids):** all of §2.F — direct repeats ≥25 bp anywhere in the assembled plasmid, inverted repeats (LTRs are themselves long direct repeats; ITRs are 145-bp palindromes), cryptic σ70 promoters and internal SD/ATG pairs driving toxic bacterial expression of a eukaryotic protein (73–89% toxicity for ≥1 TM segment), AT-window ≤55%, copy-number multiplier from the detected ori. **Plus a protocol recommendation:** Stbl3 at 30 °C for lentiviral LTRs, ΔsbcC at 42 °C for AAV ITRs.
2. **Packaging cell (HEK293T):** the transfer plasmid's cassette is transcribed as vector genomic RNA (lentivirus) or replicated by Rep (AAV). Therefore: strand-correct polyA/splice scanning on the *packaged* strand; internal polyA = 8–9× titer loss with a strong internal promoter; cryptic splice acceptors pair with the vector's major splice donor; producer-cell protein toxicity may dominate everything (up to 1000×); AAV genome size and palindrome-driven template switching.
3. **Target cell:** Kozak, uAUG, CSC/mRNA stability, CpG-ZAP density-and-spacing, CpG methylation-driven silencing, tissue codon table (opt-in, contested), and the actual expression objective.

**Conflicts this creates, which the app must surface rather than resolve silently:**
- Cryptic-promoter suppression (*E. coli*) wants GC **up**; synthesis vendors want GC **down**. Two-sided band, per-window "which side is binding" display.
- CpG depletion (target-cell immunogenicity/ZAP) wants GC **down**; the promoter CpG island wants CpG **up** for silencing resistance. Present as a tradeoff, not a default.
- The lentiviral genome and the target-cell mRNA are the **same strand** in a forward-oriented cassette but **opposite strands** in a reverse-oriented one. Cassette orientation is a required input; without it the polyA/splice analysis is exactly backwards.
- LTRs (≥600 bp direct repeats) and ITRs (145-bp palindromes) *violate* the plasmid-stability rules irreducibly. The app must whitelist them as immutable and route the user to strain/temperature mitigation instead of pretending to fix them.

---

## 4. Algorithm recommendation

**Reject** whole-transcript exact structure-aware DP in the interactive loop. LinearDesign takes 11 min exact / 2.7 min at beam 500 for 1,273 aa and >60 GB RAM at 1,450 aa; DERNA takes 6 h at ~1,450 aa — https://pmc.ncbi.nlm.nih.gov/articles/PMC12319323/. LinearDesign is also legally unusable (no OSS license, Baidu patent filing).

**Three-tier decomposition:**

**Tier A — exact Viterbi DP over the codon lattice (all Markov objectives + hard motif avoidance).**
State = (Aho–Corasick automaton state over the forbidden-motif set ∪ revcomps, [previous codon if codon-pair terms are on]). Relaxations = L × |S| × d̄. For L = 1000, |S| ≈ 600 (≈50 six-bp sites plus revcomps), d̄ = 2.94 → **1.76 M relaxations ≈ 10–30 ms fully NumPy-vectorized** over the state axis. Adding codon-pair state multiplies by ≤6 → ~10.6 M ≈ 100 ms. Precompute a codon-level transition table `T[state][64]` once (|S|×64×3 ≈ 115 k steps). This gives **exact optimality** for CAI/tAI/CPB/CFD/CpG/out-of-frame-stop terms and a **by-construction guarantee** that no forbidden motif appears — including motifs created at codon boundaries and at the insert/backbone junction (seed the automaton with the state after consuming the immutable left flank; require a non-accepting final state after the right flank).

**Tier B — DnaChisel-style localized repair for non-Markov constraints** (windowed GC, GC extent, k-mer repeats, hairpins, backbone junction). Copy the calibrated defaults verbatim: exhaustive local search when local mutation-space size < **10,000** variants (m=8 codons at d=3; m=5 at d=6), else guided random with **2 mutations/iteration**, **max 1000 iters**, **stagnation tolerance 100**, **local_extensions (0, 5)** — https://github.com/Edinburgh-Genome-Foundry/DnaChisel. Localization extension rule: motif constraints extend by (len(motif) − 1); windowed-GC by (window − 1); suppress the right-hand extension when two breaches overlap. Circular plasmids: triple to 3L, resolve only constraints overlapping [L, 2L), emit `sequence[L:2L]`, subtract L from coordinates.

**Tier C — structure, evaluated on windows only.** Slide a 100-nt window at 10-nt step; a single-codon mutation invalidates only ~10–20 windows (**~2–5 ms per proposal**, >200 proposals/s single-threaded). Full-ORF fold once at report time. Timing anchors (cubic extrapolation from the published LinearFold benchmark): RNAfold ≈0.24 s at 1 kb, ≈6.5 s at 3 kb; LinearFold ≈0.8 s at 1 kb — **do not use LinearFold below ~2 kb**, the crossover is ~1.5–2 kb. RNAplfold accessibility over 1 kb at L=W=80 is O(n·L²) ≈ tens of ms.

**Multi-objective UI.** Precompute a λ-sweep of ~24 designs by running Tier A with a (|S| × 24) score array in one vectorized pass (**~50–100 ms**, not 24× a single DP). Sliders snap instantly to the nearest precomputed point while an exact re-solve for the precise weights runs in the background. **Warn in code and docs that weighted-sum sweeps recover only *supported* (convex-hull) Pareto points**; offer COSMO-style vector-valued DP with a `vmax` dominance filter (O(L·K) bicriteria, extra O(K·log^(n−2)K) for n≥3) and ε-constraint mode ("maximize CAI subject to GC ≤ 62% and ΔG_5′ ≥ −25"). ε-constraints actually *accelerated* COSMO's search 1.25–2.21× via pruning — https://pmc.ncbi.nlm.nih.gov/articles/PMC7358382/.

**Sampling, not argmax.** Implement forward-filter/backward-sample over the same first-order chain: p(S) ∝ exp(−β·H), H = Σh_codon + ΣJ_codon,codon+1. β = 1 reproduces natural host CAI/CPB/GC statistics; β→∞ is the deterministic argmax — https://pmc.ncbi.nlm.nih.gov/articles/PMC11345917/. Same cost as the Viterbi pass. Expose β as a diversity slider, and **emit a panel of 3–8 diverse candidates**, never one "optimal" sequence — the honest response to a field where optimizer choice is a coin flip and the expression landscape shows negative epistasis.

**10-second budget for a 1,000-codon ORF:** ≤0.2 s preprocessing (mutation space + automaton), ≤0.3 s Tier A, ≤4 s Tier B repair, ≤3 s polish, ≤1 s final validation + one full fold, ~1.5 s slack.

**Determinism contract (non-negotiable):** deterministic tie-break in the DP (lowest codon index); explicit `np.random.default_rng(seed)` threaded through every stochastic stage (DnaChisel has **no seed parameter** and uses the global NumPy RNG — wrap or reimplement); persist {seed, app version, codon-table name+version, weight vector, constraint-set hash} in the output.

**Independent final validator.** Only two mechanisms guarantee constraints: domain restriction (Tier A) and reject-and-repair + full independent revalidation. Penalty weights never guarantee. Re-scan every forbidden motif on both strands of the full assembled plasmid, re-check every GC window, re-check repeats against the backbone, and re-translate. **Refuse to emit on failure**; report the minimal conflicting constraint set (a zero-variant merged segment in the mutation space is a precise infeasibility certificate).

**Architecture:** Tauri v2 + React/TypeScript + a PyInstaller `--onedir` Python sidecar on 127.0.0.1 with an ephemeral port and a per-launch bearer token. Rationale: real filesystem access for GenBank round-tripping (the browser option fails — `showSaveFilePicker` is Chrome/Edge only), the React ecosystem for a custom sequence viewer, an unmodified CPython hosting Biopython/DnaChisel, and a hard OpenAPI process boundary that lets parallel coding sessions work independently. `--onefile` is disqualified: it re-extracts to `_MEI_xxxxxx` on every launch and hides the real child PID from Tauri so `kill()` orphans the server. Run each job in a `ProcessPoolExecutor` worker so cancel == `terminate()` (Python threads cannot be killed). Stream progress via Tauri **Channels**, not events.

---

## 5. Data assets

| Asset | Source | License | Size / notes |
|---|---|---|---|
| Codon usage tables (genomic) | Build at **build time** from NCBI RefSeq CDS FASTA, https://api.ncbi.nlm.nih.gov/datasets/v2alpha/ | NCBI "places no restrictions" | 64 rows ≈ 2 KB each; whole host panel <1 MB. Verified accessions: *E. coli* K-12 MG1655 GCF_000005845.2 (4,290 CDS); BL21(DE3) GCF_000009565.1; *B. subtilis* 168 GCF_000009045.1; *S. cerevisiae* GCF_000146045.2; *K. phaffii* GS115 GCF_000027005.1; CHO GCF_003668045.3 (21,776); human GCF_000001405.40; mouse GCF_000001635.27; *S. frugiperda* GCF_023101765.2; *T. ni* GCF_003590095.1; *C. reinhardtii* GCF_000002595.2. **N. benthamiana has no annotated RefSeq** — use NbLab360/Sol Genomics or substitute *N. tabacum* and label it. |
| Codon + codon-pair + dinucleotide tables | CoCoPUTs/HIVE-CUT, https://dnahive.fda.gov/dna.cgi?cmd=cuts_main | US federal work | ~288 M CDS. **Mirror manually at build time** — the download UI is JS-only with no stable static URLs. |
| Highly-expressed reference sets | PaxDb v5/v6, https://pax-db.org/ | **CC BY 4.0** | Top 5% by protein abundance (min 100 genes). Fallback: cytosolic ribosomal proteins + EF-Tu/EF-G/EF-1α + GAPDH + chaperonins. Last resort: bottom 10% by ENc (genomes >1,000 genes only). |
| Sharp & Li 1987 *E. coli* w-table | Biopython 1.81 `CodonUsageIndices.SharpEcoliIndex`, https://raw.githubusercontent.com/biopython/biopython/biopython-181/Bio/SeqUtils/CodonUsageIndices.py | Biopython/BSD-3 | 60 lines. **Removed in Biopython ≥1.82 — copy it out now** as a CAI regression fixture. |
| tRNA gene copy numbers | GtRNAdb R22, https://gtrnadb.ucsc.edu/ | No explicit license; citation requested (Chan & Lowe NAR 2016) | ~129 KB/genome; ~1.5 MB for 13 hosts. Use `-confidence-set.out` for eukaryotes. **T. ni absent** — run tRNAscan-SE 2.0 (GPL, offline) or fall back to *S. frugiperda* with a caveat. |
| tAI s-values (dos Reis + Tuller) | codon-bias v0.5.0 `tAI_svalues_*.csv`, https://github.com/alondmnt/codon-bias | **MIT** | 9 rows each. Gate lysidine C:A on bacteria/archaea. |
| Species-specific stAI s-vectors | stAIcalc, https://academic.oup.com/bioinformatics/article/33/4/589/2593585 | Free academic; verify | 100 species precomputed. |
| Codon Stability Coefficients | https://elifesciences.org/articles/45396/figures | CC-BY | Static CSVs, multiple human cell lines. |
| Splice-site tables | MaxEntScan (Yeo & Burge 2004); or train your own max-entropy model on GENCODE, https://www.gencodegenes.org/ | GENCODE free; MaxEntScan redistribution ambiguous | Extract only flanking 9-mers/23-mers — a few MB. **Reimplementing from GENCODE sidesteps the license question.** |
| PolyA hexamer priors | PolyA_DB 3 / PolyASite 2.0, https://academic.oup.com/nar/article/46/D1/D315/4561640 | Free academic | 18 hexamers + downstream-element position distributions. |
| Overhang misannealing matrices | tatapov (Potapov 2018, Pryor 2020), https://github.com/Edinburgh-Genome-Foundry/tatapov | **MIT** | Bundle the CSVs so it works offline. |
| Vendor rule profiles | Twist FAQ + DOC-001081; GenScript GenTitan; IDT gBlocks FAQ | Vendor docs, reference only | JSON/TOML with citation URL and "last verified" date. **Assume 12-month staleness** — Twist changed homopolymer 14→30 between 2023 and 2026. |
| Reference element library | Addgene GenBank corpus | Public | AAV2 ITR flip/flop (145 nt), mutant trs-deleted ITR (~130), bGH polyA 225 bp, SV40 late 135, NRP1 32, WPRE 600, WPRE3 247, cPPT/CTS 118/178, RRE, HIV Ψ, 2A peptides, splice-safe recoded V5. |
| Folding engine | **seqfold** (MIT, Rust/PyO3, abi3 wheels, zero deps), https://github.com/Lattice-Automation/seqfold as bundled default; ViennaRNA 2.7.2 as optional user-installed backend | seqfold MIT; ViennaRNA **custom non-OSI: no redistribution for a fee** | seqfold is MFE-only (no partition function, no accessibility, no constraints). ViennaRNA is what every literature threshold was calibrated against — **email TBI Vienna for bundling permission** (the license explicitly invites this). |

**Licensing landmines to avoid:** Codon Statistics Database (CC BY-NC), EMBOSS `data/CODONS/*.cut` (GPL-2+, and vintages are 1994–2005), REBASE (CC BY-NC — inherited by Biopython's `Bio.Restriction` enzyme table; hand-curate ~60 enzymes from vendor catalogues instead), NUPACK (non-commercial, **no GUIs**, paid subscription), LinearDesign (no license + Baidu patent), pLannotate/GenoLIB (GPL-3.0 + SnapGene-derived DB), D-Tailor (CC BY-NC), ATGme (non-commercial), mRNAid (CC BY-NC), SpliceAI models (CC BY-NC), CodonBERT/RiboNN weights/mRNA-LM/DeepCodon/Nucleotide Transformer (all non-commercial). **Permissive and usable:** DnaChisel, DnaCauldron, GoldenHinges, genedom, Geneblocks, tatapov, DnaFeaturesViewer, python_codon_tables (CC0), codon-bias, seqfold, pyahocorasick (BSD-3), OR-Tools CP-SAT (Apache-2.0), pymoo (Apache-2.0), CodonTransformer/ColiFormer (Apache-2.0 / CC BY 4.0), APARENT2 (MIT), @teselagen/ove and seqviz (MIT), Biopython, pydna (BSD).

---

## 6. Prior art gaps

DnaChisel is the reference architecture and its constraint-localization + mutation-space model should be reimplemented faithfully. But it has **no mRNA folding model at all** (`AvoidHairpins` is an exact 20-mer revcomp substring match, not thermodynamics), no synthesizability model, no splice/promoter/polyA/uAUG specs, no vector context, greedy order-dependent objective optimization, linear scalarization only, and no random seed. GeneOptimizer is a strictly 5′→3′ beam-width-1 sliding window (m ≈ 4 codons) that its own paper admits cannot resolve a repeat between the 5′ and 3′ halves except by inserting worse codons downstream, and implements no splice/polyA/TATA/RBS/codon-pair terms — https://pmc.ncbi.nlm.nih.gov/articles/PMC2955205/. GenScript and Twist are closed web tools tied to ordering. Benchling is cloud-only with frequency-table back-translation.

**What we can do that none of them do:**

1. **Optimize in the context of the assembled plasmid.** Splice the candidate into the user's backbone in memory and evaluate all constraints on the *circular* product — junction-spanning restriction sites, GC windows crossing the origin, repeats shared with ITRs/LTRs/WPRE/polyA, homology driving recombination. Every competitor optimizes a free-floating linear CDS.
2. **Auto-detect the Golden Gate destination cassette.** A pair of inward-facing same-enzyme Type IIS sites flanking a drop-out region (ccdB/lacZα/sfGFP/RFP) yields the enzyme, both 4-nt overhangs, and the insertion coordinates with zero user input.
3. **Back-translate the whole cassette, not just the ORF.** Tags, linkers, 2A peptides and signal peptides go through the same optimizer, so P2A comes out domesticated for the user's enzyme, `(GGGGS)3` comes out without a perfect nucleotide repeat, and **V5 comes out without its `G|GTAAG` donor** — a documented 17/17 failure mode no tool addresses.
4. **The ZAP CpG density-and-spacing metric** (≥14 CpGs at ≤14-nt mean spacing in a sliding window), separated from TLR9 hexamer counting and CpG-island detection. Publishable-quality and it changes the answer: a construct can be CpG-rich globally yet ZAP-invisible.
5. **Out-of-frame stop-codon density** as a first-class objective (PNAS 2026), trivially computable, opposes over-optimization.
6. **Three-way context reasoning** (propagation host + packaging cell + target cell) with strand-correct scanning driven by cassette orientation.
7. **An explicit conflict panel** rather than silent resolution: NcoI vs Kozak, Bxb1 att vs BsaI, CpG depletion vs GC floor, Pol III TTTT vs AT-rich proteins, *E. coli* AT-window vs vendor GC ceiling.
8. **Evidence-strength badges per rule** (EVIDENCE-BACKED / VENDOR-ASSERTED / FOLKLORE), with folklore rules default-off, and a per-rule report of how many codons changed and what it cost on every other objective. Given the 20–25% win rate of vendor optimizers, a tool that *quantifies the cost of each constraint* is more defensible than one that silently applies twelve.
9. **User-supplied codon tables** built in-app from their own high-expressing clone FASTA — dissolving the BL21-vs-K-12-vs-W3110, CHO-K1-vs-S-vs-DG44, GS115-vs-CBS7435 strain-mismatch problem entirely.
10. **A calibration view:** paste sequences with measured expression and see which metric actually correlates *in your system*. Published R² for codon metrics ranges 0.03–0.23 while 5′ models reach r ≈ 0.76 — per-user calibration beats any fixed weighting we could ship.
11. **Protocol output derived from the scan:** strain, temperature, antibiotic concentration, colony-picking guidance, expected Golden Gate fidelity, colonies to pick for 95% confidence.

---

## 7. Risks and known-hard problems

**Scientific.** (i) The ceiling is real: ~14% of protein-level variance explained by all computable features (Cambray), and negative epistasis makes the landscape rugged enough that free methods match premium commercial ones — https://doi.org/10.1101/2025.06.03.657573. Never report a predicted expression number; report rank/percentile with a confidence band. (ii) Over-optimization has documented costs: CAT lost ~20% specific activity when 16 rare codons were sped up; SufI folding was impaired; codon-optimized luciferase made more protein but dramatically less light; codon-optimized FIX adopts a distinct conformation. (iii) No published dataset cleanly quantifies the titer cost of a cryptic splice donor inside a therapeutic ORF — every quantified case is confounded. Say so. (iv) No published ΔG threshold defines when a DNA hairpin truncates AAV genomes; any number we encode (−20 kcal/mol) is an engineering guess. (v) tAI in mammalian hosts rests on tGCN as a proxy, tissue-specific tRNA pools, and U34 wobble modification state that no static table captures; direct tRNA-seq methods disagree with each other more than with tGCN (nano- vs hydro-tRNAseq ρ = 0.182).

**Technical.** (i) Windowed GC in the automaton DP requires ~17 codons of history — naive state augmentation is intractable; either Lagrangian relaxation on a multiplier, a bounded running-deviation state, or accept Tier-B repair. No literature precedent. (ii) Over-constrained infeasibility is the dominant real-world UX; DnaChisel raises `NoSolutionError` constantly. Report the minimal conflicting set and offer specific relaxations. (iii) Absolute ΔG values differ between seqfold and RNAfold; silently applying a ViennaRNA-calibrated threshold (−39, −30/−50/−61) to seqfold output is the most likely correctness bug in the whole feature. Pin engine, parameter set, dangles model and temperature in a provenance record and gate thresholds on the engine. (iv) GenBank 1-based inclusive vs Biopython 0-based half-open, plus circular origin-spanning `join()` features, plus `SeqRecord` slicing silently dropping partially-overlapping features and **all** `.annotations`/`.dbxrefs` — write an explicit interval remapper and a golden round-trip test suite over real Addgene plasmids. (v) SnapGene `.dna` is **read-only in every open implementation** (Biopython handles only packets 0x00/0x05/0x06/0x0A) — import `.dna`, export `.gb`, never promise round-trip. (vi) macOS notarization of a hardened-runtime CPython needs `com.apple.security.cs.allow-unsigned-executable-memory` and `disable-library-validation`; prove this path in week one with a hello-world sidecar. (vii) Tauri has no delta updates — every release re-downloads the full ~100 MB Python payload; make updates user-initiated.

**Legal.** ViennaRNA bundling requires written permission (routinely granted, long latency — start immediately). LinearDesign is patent-encumbered: any codon-lattice × folding-DP implementation, including a from-scratch CDSfold-style one, needs an FTO review. REBASE's NC clause is inherited through `Bio.Restriction`.

---

## 8. Contradictions and disagreements

1. **Boël 2016 vs Kudla/Goodman/Cambray.** Boël (6,348 genes, T7-driven, soluble-protein readout) finds codon content ~3–5× more influential than folding beyond codon ~16; the others (host RNAP, single-gene libraries, fluorescence readout) find structure dominant throughout. Unresolved whether the discrepancy is T7 vs σ70, many-genes vs many-variants-of-one-gene, or solubility vs fluorescence. **Consequence:** we cannot principledly set the relative weight of 5′-structure vs codon terms past codon 17. Expose both as weighted hypotheses.

2. **Structure sign flips between contexts and nobody has reconciled them.** Kudla says minimize 5′ structure (bacteria, expression); LinearDesign says maximize global structure (IVT mRNA, half-life, 5.1× validated); the Expi293F glycoprotein benchmark says LinearDesign was the *worst* of five schemes for a DNA transgene despite 2-fold lower normalized MFE, with **no correlation between MFE and yield**. These are three different objectives on three different molecules. Model them as separate, sign-switched, window-specific terms and never average them into one "structure" slider.

3. **Codon-pair effect sign.** Coleman 2008 says under-represented pairs slow translation (attenuation); Gutman & Hatfield say *over*-represented pairs translate slower; Simmonds shows the poliovirus phenotype is CpG/UpA, not pairs. CPB gets near-zero default weight.

4. **G-quadruplex in mammalian cells.** Reporter assays show 55–85% repression from 5′UTR rG4s; transcriptome-wide in-cell probing shows they are globally unfolded (median folding score 0.06). Same motif, two independent flags: default-on for bacterial CDS and DNA synthesizability, default-off for mammalian translation.

5. **Tissue-specific codon tables.** A 2023 Genome Biology study validated tissue-matched eGFP/mCherry ratios in HEK293T vs A549 (p<0.05, AUC>0.70 for kidney/breast/lung/rectum/tonsil); a 2006 MBE analysis found no evidence for tissue-specific translational adaptation and attributed the signal to isochore GC. Between-tissue variance is far smaller than between-gene variance. **Opt-in, clearly labeled contested.**

6. **ZAP in engineered vectors.** Strong effects in replicating viruses (Takata, Ficarelli); one 2021 study reports "minimal impact of ZAP on lentiviral vector production and transduction efficiency." Present the metric, not a mandate.

7. **Twist's own documentation contradicts itself:** the FAQ still lists homopolymer ≥14 bp as a hard reject while the 2026 Complex Genes product explicitly accepts to 30 bp. The 14 figure is now the Standard-tier line, not a reject line. Similarly, Twist's codon optimizer targets GC 25–65% and max homopolymer <10 — tighter than its own acceptance envelope. **Optimize to the tool targets, accept up to the acceptance envelope.**

8. **The widely repeated "Twist rejects >65% GC over 50 bp" is not supported by any Twist document.** Twist's published 50-bp window bound is 10–90%; 65% appears only as a *global* recommendation. Several numeric thresholds circulating in this field (GENEWIZ/Azenta and Telesis Bio publish none at all) are unsourced.

9. **Best-codon disagreement in *E. coli*.** Welch's charged-tRNA-derived favorites (Ser AGC, Thr ACG, Leu TTG — explicitly *not* CTG) vs Boël's 6AA/31C-FO sets have never been reconciled experimentally, and it is unclear whether either transfers across strains.

10. **Harmonization.** Best for 8/18 mammalian targets but also the *most variable* — consistent with either a real target-dependent benefit or noise. No large pre-registered harmonization-vs-native trial exists.

11. **Method disagreement is itself the finding.** Across 10 tools on the same proteins, GC spanned 51–64% (*E. coli* insulin) and 48–66% (CHO adalimumab HC), with ΔG from −156 to −104 kcal/mol; tools modeling codon context showed "no significant improvement." This is the strongest argument for shipping a Pareto panel with visible provenance rather than a single opaque answer.