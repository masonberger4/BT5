# BT5 — Codon Optimization & Back-Translation App

## Context

You want a locally-run app that takes a protein sequence plus a vector backbone you already have,
and emits a DNA coding sequence optimized across competing objectives — protein expression, DNA
synthesizability, viral titer, plasmid stability, and avoidance of the problem motifs molecular
biology actually cares about (Chi sites, telomere repeats, cryptic promoters, splice sites, polyA
signals). The output must work *in your backbone*, and the UI must let you steer the trade-offs
rather than hand you one opaque answer.

`masonberger4/BT5` is empty apart from a README, so this is greenfield. The plan is sized to be
executed by simultaneous Claude Code sessions, each owning one lane and opening its own PR.
**This document is the sole input to those sessions.**

### Citation policy

Every factual claim in this plan carries a link to its primary source. Rules inherit that
requirement: `Spec.citations` is a **required tuple** and CI rejects any rule with an empty one.
Where the literature disagrees, both sides are cited and the rule is badged `contested`. Where a
number comes from a vendor rather than a paper it is badged `vendor_asserted` and carries a
`last_verified` date, because these drift — Twist moved its homopolymer limit from 14 to 30 bp
between 2023 and 2026 ([Twist FAQ](https://www.twistbioscience.com/faq/gene-synthesis)).

### Research basis

Two research fan-outs and a judged design pass (24 agents, ~3.6M tokens) produced four documents.
Copy them into `docs/research/` in the foundation PR. They currently sit in
`/tmp/claude-0/-home-user-BT5/f8a714a6-08d8-56fc-a44d-5336c2583b6a/scratchpad/`.

| File | Contents |
|---|---|
| `brief.md` (60KB) | Constraint/objective model with exact motifs and thresholds, context matrix, solver design, data assets with licenses, competitive gaps, contradictions |
| `critique.md` (27KB) | Completeness critique: mandatory Tier-1 gaps, adjudicated contradictions, validation strategy |
| `githubPlan.md` (121KB) | Full CI/ruleset/merge-gate design with committable YAML |
| `proposal_C.md` + `verdicts.json` | The winning architecture and the three judges' corrections |

---

## What "backbone" means — the central data model decision

This governs which bases the optimizer may touch, and everything else follows from it.

**The CDS (designable — every base back-translated):** everything between the start and stop codon,
as one continuous reading frame. The gene of interest **plus** all in-frame cassette elements — N-
and C-terminal tags, linkers, 2A peptides, signal peptides, protease sites. These are part of the
ORF and go through the same optimizer.

This is a genuine differentiator, because naively back-translated cassette elements are a documented
failure source. The standard V5 tag encoding contains a `G|GTAAG` splice donor and spliced in
**17/17** genes tested, with 13/17 randomly chosen genes showing aberrant splicing from vector/tag
context ([PMC9379414](https://pmc.ncbi.nlm.nih.gov/articles/PMC9379414/)). A naive `(GGGGS)₃` linker
is a perfect nucleotide repeat that fails synthesis and recombines. Two copies of the same 2A peptide
in one ORF is a perfect direct repeat — the literature explicitly requires different 2A peptides and
≤85% nucleotide identity between them
([PMC5438344](https://pmc.ncbi.nlm.nih.gov/articles/PMC5438344/)).

**The backbone (fixed — never edited, but always scored against):** everything else. Promoter,
5'UTR, 3'UTR (**including WPRE**), polyA signal, origin of replication, selection marker, LTRs, ITRs,
insulators.

> **WPRE is backbone, not CDS.** It sits after the stop codon and is never translated, so by the
> definition above it is fixed context. It still matters three ways: it is worth ~2–10× expression,
> transgene-, promoter- and vector-independently
> ([Donello/Salk](https://www.salk.edu/wp-content/uploads/2016/04/WPRE-Donello-RD9436-NCD-FY2016.pdf),
> [Nature Gene Therapy](https://www.nature.com/articles/3302979)); at ~600 bp (WPRE3: 247 bp) it eats
> the AAV/lenti size budget; and it is a repeat liability the insert must not duplicate.

**The backbone is not passive.** Three ways fixed sequence changes what the optimizer must do:

1. **The 5'UTR is required input for the highest-weight objective.** The −4…+37 window that explains
   44–59% of expression variance in bacteria *spans the UTR/CDS junction*
   ([Kudla 2009, PMC3902468](https://pmc.ncbi.nlm.nih.gov/articles/PMC3902468/)) — it cannot be
   computed from the CDS alone. Independently, synthetic 5'UTRs produce translation differences over
   **two orders of magnitude** in mammalian lines and up to 3.5× titer changes, and the ranking of
   UTRs is **cell-type dependent**
   ([NAR 2020](https://academic.oup.com/nar/article/48/20/e119/5922799),
   [CHO study](https://www.authorea.com/doi/full/10.22541/au.173150120.01568890/v1)).
   **Consequence:** with no annotated 5'UTR the app degrades honestly and says the 5'-structure
   objective is unavailable — it must never silently fold the CDS alone and report a number as if the
   UTR were there. Gate G6 enforces this.
2. **uAUGs and uORFs in the user's own 5'UTR are real and unfixable by codon choice.** uORFs occur in
   ~half of human transcripts and typically cut protein 30–80%
   ([PMC7100133](https://pmc.ncbi.nlm.nih.gov/articles/PMC7100133/)). Report-only (`HARD_CHECK`).
3. **Backbone repeats set the insert's repeat budget.** Plasmid stability and Gibson assembly both
   depend on uniqueness across the *whole* construct, so a 20 bp sequence fine in isolation is
   disqualifying if the backbone already contains it.

---

## The finding that reshapes the product

- **CAI does not predict expression and must not be the objective.**
  [Kudla 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC3902468/) (154 synonymous GFP variants,
  250-fold expression range): 5′ folding free energy explained 44% of variance (r = 0.66; 59% in a
  second promoter system); CAI gave r = 0.14, not significant; whole-mRNA MFE r = 0.16, not
  significant. [Welch 2009](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0007002)
  states flatly that "CAI has no value in predicting gene expression," with their deliberately
  high-CAI control expressing at ~15% of the best variant.
  [Ranaghan 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC7893858/) benchmarked nine optimizers:
  "a roughly equivalent chance that an algorithm-optimized CDS will increase or diminish recombinant
  yields," with three tools non-deterministic (one returning 35–39% pairwise codon identity across
  ten identical submissions).
- **In mammalian cells it is often worse than doing nothing.** An 18-glycoprotein / 90-screen Expi293F
  benchmark concluded codon optimization of human proteins in a human cell line "did not generate
  increased yields," with native and harmonized most consistent
  ([Protein Innovation 2026](https://proteininnovation.org/2026/03/codon-optimization-native-codon-mammalian-protein-expression/)).
  A Pichia study found CAI *negatively* correlated with titer, −0.81 for trastuzumab
  ([EuropePMC 41701818](https://europepmc.org/article/MED/41701818)).
- **The ceiling is low.** [Cambray 2018](https://www.nature.com/articles/nbt.4238) (244,000 designed
  sequences, full factorial): all computable design features together explain 5–31% (mean ~14%) of
  protein-level variance.
- **Several standard rules are folklore.** Internal Shine-Dalgarno "pausing" was an artifact of
  selecting 28–42 nt ribosome footprints
  ([Mohammad, Green & Buskirk, eLife](https://elifesciences.org/articles/42591)). Codon-pair-bias
  attenuation is largely a CpG/UpA effect — deoptimized pairs with unchanged dinucleotides had
  wild-type fitness ([eLife 04531](https://elifesciences.org/articles/04531)). The slow-5′ "ramp"
  gave 67–71% of the fast construct's GFP ([eLife 89656](https://elifesciences.org/articles/89656)).
  "High GC kills AAV titer" is contradicted by the only controlled dataset: a designed stuffer at
  GC 43.5–44.8% cut yield up to 68%, while a *lower*-GC (33.8–34.7%) natural stuffer of identical
  length cost neither yield nor bioactivity — the liability is repetitiveness
  ([PMC12207685](https://pmc.ncbi.nlm.nih.gov/articles/PMC12207685/)).

**Consequences, structural rather than cosmetic:**

1. The app **never reports a predicted expression number.** Ranks, percentiles against a
   random-synonymous null, confidence bands. A CI grep over the *generated* OpenAPI bans field names
   matching `/predict|titer|yield|expression_level/`.
2. Every rule carries an **evidence badge** — `evidence_backed`, `contested`, `vendor_asserted`,
   `folklore`. Only `folklore` defaults off. (Collapsing vendor-asserted into folklore would ship
   every manufacturability rule disabled, which is backwards — vendors enforce those at order time.)
3. **"Use the native sequence" is a first-class output** — `DesignResult.native_baseline` is a field,
   not a UI afterthought.
4. The app emits a **panel of 3–8 genuinely different candidates**, not one "optimal" sequence.
5. **The defensible value is not better expression prediction.** It is guaranteeing hard constraints,
   doing so *in the context of the assembled plasmid*, and quantifying what each constraint costs.

### Where the differentiation actually is

Every competing tool optimizes a free-floating linear CDS. DNA Chisel has no thermodynamic folding
model at all — its `AvoidHairpins` is an exact 20-mer reverse-complement substring match
([DnaChisel](https://github.com/Edinburgh-Genome-Foundry/DnaChisel)) — and GeneArt's GeneOptimizer is
a strictly 5′→3′ beam-width-1 sliding window whose own paper admits it cannot resolve a repeat
between the 5′ and 3′ halves ([PMC2955205](https://pmc.ncbi.nlm.nih.gov/articles/PMC2955205/)).

1. Constraint evaluation on the assembled circular product including backbone junctions.
2. **Whole-CDS back-translation** — fixing the V5 donor, the `(GGGGS)ₙ` repeat, the duplicate-2A repeat.
3. **Gibson-aware uniqueness, on a 2-D repeat-risk surface.** No exact repeat ≥20 bp shared between
   insert and any other part of the construct — repeats and high GC are the documented Gibson
   misassembly cause, and no two junctions may share a homology sequence
   ([Addgene](https://blog.addgene.org/plasmids-101-gibson-assembly); NEB arms 15–40 bp, Tm > 48 °C).
   Risk is scored over **(repeat length × spacer distance)**, not a flat length cutoff, because
   RecA-independent recombination is strongly proximity-sensitive. This matters most in exactly the
   LVV/AAV workflow: `recA⁻` strains suppress only the >200–300 bp RecA-dependent pathway, so the
   15–100 bp repeats codon choice controls are the app's responsibility, not the strain's (see Q4a).
4. The **ZAP CpG density-and-spacing metric** — flag any 200-nt window with ≥14 CpGs at mean
   inter-CpG spacing ≤14 nt; spacing ≥32 nt is not restricted, and inhibition magnitude is *not*
   correlated with total CpG count ([PMC9519448](https://pmc.ncbi.nlm.nih.gov/articles/PMC9519448/)).
   Kept separate from TLR9 hexamer counting ([JCI 68205](https://www.jci.org/articles/view/68205))
   and CpG-island detection. Report the worst window, never a global count.
5. **Three-way context reasoning** — E. coli propagation + packaging cell + target cell — evaluated
   simultaneously, never merged, with strand-correct scanning driven by cassette orientation.
6. An explicit **conflict panel** (NcoI CCATGG ⊂ Kozak GCCACCATGG; CpG depletion vs vendor GC floor;
   E. coli AT-window vs vendor GC ceiling).
7. **Per-user calibration** — paste sequences with measured expression, see which metric correlates
   *in your system*. Published R² for codon metrics spans 0.03–0.23 while 5′ models reach r ≈ 0.76
   ([NAR 2023](https://academic.oup.com/nar/article/51/5/2363/7016452)), so per-user calibration beats
   any fixed weighting we could ship.

---

## Locked decisions

**1. Stack — Python core + local web UI.** Python 3.11 engine (numpy/biopython/ViennaRNA) behind a
FastAPI server on 127.0.0.1; React/TypeScript UI in the browser at localhost. Chosen for scientific
library access and the cleanest split for parallel sessions. **Explicitly designed for easy later
transfer to a desktop app** — the process boundary is a hard OpenAPI contract, so a later Tauri v2
wrapper with a PyInstaller `--onedir` sidecar is a packaging change, not a rewrite. Install via `uv`.

**2. v1 covers all four delivery contexts:** mammalian + viral (lenti/AAV), E. coli, yeast/insect, IVT
mRNA. **Deep curation** for **human/HEK293, CHO, and E. coli** only; other organisms are
codon-table-level behind the same interface.

**3. Trade-off UI — presets + sliders + candidate gallery.**

**4. Licensing — personal use, best available tools. You have ViennaRNA bundling permission**, so
ViennaRNA 2.7.2 is the bundled default. This is material: every literature threshold — the Boël
−39 kcal/mol dual gate ([PMC5054687](https://pmc.ncbi.nlm.nih.gov/articles/PMC5054687/)), the
cap-proximal −30/−50/−60 ladder — was calibrated against ViennaRNA, so we apply them directly rather
than transferring across engines, which the research named the most likely correctness bug in the
feature. `FoldEngine` stays an interface; there is no cross-engine calibration risk to manage.

**5. Inputs — GenBank, SnapGene, FASTA.** GenBank with annotation-driven detection (default path);
SnapGene `.dna` read-only import (every open implementation is read-only — Biopython handles only
packets 0x00/0x05/0x06/0x0A — so import `.dna`, export `.gb`, never promise round-trip);
FASTA/paste with the insertion point marked manually.

**6. Output — annotated construct + vendor order CSV.** Optimized CDS, complete assembled GenBank,
QC report, and **a simple CSV for bulk upload to IDT or another vendor**, with complexity rules
pre-checked. **No assembly design tooling.**

**7. Assembly context — Gibson, not Golden Gate.** Type IIS avoidance (BsaI/BsmBI/BbsI/SapI) remains
available as *optional, default-off* rules; Golden Gate destination auto-detection is dropped. What
Gibson needs is junction and repeat uniqueness, now the better-justified core rule.

**8. Buildout — independent modules, one PR each.**

Verified: Python 3.11.15, Node 22.22, Rust 1.94, uv 0.8.17. `ViennaRNA` 2.7.2 ships prebuilt wheels;
`biopython`, `snapgene-reader`, `pydna`, `dnachisel`, `seqfold` all install cleanly.

### The order file — exact IDT format (confirmed)

Verified against the IDT **eBlocks 96-well plate upload template** you supplied
(`eblocksplateuploadtemplate96.xlsx`). It is simpler than expected:

- **Format is `.xlsx`, not `.csv`**, for plate orders. The exporter must write a real workbook.
- **The sheet name is the plate name** — the template ships a single sheet literally named
  `Plate Name`, which the user renames. So plate identity travels as sheet metadata, not a column.
- **Exactly three columns:** `A = Well Position`, `B = Name`, `C = Sequence`.
- **96 rows pre-filled with well positions in row-major order** — `A1…A12, B1…B12, … H1…H12`.

Implementation: emit the workbook with wells assigned row-major, `Name` carrying the construct name
plus the short design hash (so the tube label traces back to the run that produced it), and `Sequence`
as bare ACGT. Keep a plain `Name,Sequence` CSV path for non-plate/tube orders and other vendors.
Vendor profiles live in a data file with a `last_verified` date, since these templates drift.

Screening burden is reported alongside: P(perfect clone) ≈ exp(−L/E), E ≈ 7,500 bp (Twist) or
5,000 bp (IDT eBlocks), so the report states how many colonies to pick for 95% confidence.

---

## Mandatory correctness items

1. **Genetic code tables.** Table 12 (*Candida*) reassigns **CTG = Ser, not Leu**; table 4 makes
   TGA = Trp; table 2 makes AGA/AGG stops. Wrong table ⇒ silently wrong protein no assay catches for
   months. The table lives on `TranslationUnit` (not per-segment — one ribosome reads one ORF under
   one code, and a per-segment field is incoherent *and* passes verification while wrong), is a
   **required field with no default**, is cross-validated against
   `HostProfile.locked_translation_table_id`, printed in the report, and written to `/transl_table`.
   Property-test `translate(back_translate(p, t), t) == p` over all 30+ NCBI tables.

2. **Objective normalization is the entire UI promise.** Normalize every objective to a percentile
   against 200–500 random synonymous variants of this protein. The null must be computed **on the
   assembled construct** (otherwise the percentile is against a distribution that never contained a
   backbone) and **uses windowed folding only** — see Q2 below for the arithmetic.

3. **The default weight vector is the product.** 90% of users never move a slider. Every soft rule
   carries a CI-enforced non-empty `weight_provenance`; presets carry a `rationale`.

4. **Repetitive proteins.** Max-CAI collapses to one codon per amino acid, and one-codon-per-AA
   back-translation of a repeat-containing protein produces *perfect* nucleotide repeats — and
   repetitive 9-mers per 100 bp plus longest-repeat length are the two highest-importance features in
   the best published synthesis-success model (random forest over 1,076 real vendor outcomes,
   F1 0.928 — [ACS Synth Biol](https://pubs.acs.org/doi/10.1021/acssynbio.9b00460)), as well as the
   main Gibson misassembly cause. The proteins people actually express are the repetitive ones:
   antibodies, scFv/CAR, Fc fusions, `(GGGGS)ₙ`, duplicate 2A peptides, TALEs, zinc fingers, His tags.
   Detect repeats in the **input protein before any codon is chosen** and force divergent synonymous
   assignment across copies — ≤85% nucleotide identity, no exact match ≥15–20 bp. A constraint that
   overrides codon optimality. trastuzumab HC and an anti-CD19 scFv-CAR are in `benchmarks/panel.json`
   from PR #0. Encode His6 as alternating `CACCATCACCATCACCAT`, which Twist calls out by name.

5. **Biosecurity.** BT5's core function is the textbook method for defeating nucleotide-homology
   synthesis screening, so **screen the input protein before optimization** — protein-level screening
   is the one layer this tool cannot itself defeat. `DesignRequest.screen_token` is a **required,
   non-defaulted field** verified against `sha256(protein)`; the pipeline refuses on
   `verdict == "block"`. `NullScreen` reports `status="not_run"`, **never** `"clear"`. **No "minimize
   identity to a reference" objective anywhere**, and `KmerIndex.of()` takes a `Construct` and nothing
   else, with a CI grep asserting no constructor accepts an external database.
   Ships **on by default** — see Q10. Document in `docs/biosecurity.md`.

6. **Harmonization requires a native CDS.** `ProteinInput.native_cds` is required to enable it; grey
   out otherwise; never silently substitute. Report honestly that harmonized was best for 8/18
   mammalian targets but also the *most variable*.

7. **Introns.** There are no introns in the CDS. Splice scanning applies to the CDS to remove
   *cryptic* donors/acceptors arising from codon choice. The exception is an intron deliberately
   placed to aid expression (normally in the 5'UTR, therefore backbone): if annotated,
   `SegmentKind.ANNOTATED_INTRON` marks it immutable **and exempt from scanning**. The same mechanism
   (`WHITELISTED_REPEAT`) covers LTRs and ITRs, which irreducibly violate the repeat rules.

---

## Architecture

Winner of the judged design pass: **risk-first vertical slice** (unanimous), grafted with the
domain-layered proposal's rule-catalog conventions and the contract-first proposal's port
record/replay and amendment protocol.

**The thesis:** six of the eight hardest problems push on one interface — the shape of "a thing that
evaluates a rule against a sequence." Circularity, immutable regions, three simultaneous contexts,
strand-of-interest, per-window conflict attribution, minimal-conflicting-set reporting, and objective
normalization are *all* properties of that type. If it is wrong, every lane's rule code is wrong with
it, and it surfaces in week four. So: build a working, ugly, one-host, three-rule, one-plasmid product
in week one that exercises all of them, run numeric gates, then freeze.

### Frozen interfaces — `packages/engine/src/bt5/core/`

```python
@dataclass(frozen=True, slots=True)
class Interval:
    """Half-open [start, end), 0-based, CONSTRUCT coordinates.
    end > length means the interval wraps the origin. Exactly ONE representation;
    GenBank join() normalises into it and re-splits on export."""

    start: int
    end: int
    strand: Strand = 1


@dataclass(frozen=True)
class Construct:
    """The ONLY thing a rule is ever evaluated against. No API anywhere in BT5
    evaluates a rule against a bare str — which is what makes 'evaluate on the
    assembled circular plasmid' the default rather than something to remember."""

    sequence: str
    topology: Topology
    features: tuple[Feature, ...]
    editable: tuple[Interval, ...]  # THE CDS. Complement is backbone, immutable by construction.
    codon_map: tuple[Interval, ...]  # one per codon, may wrap, may be discontiguous
    annotations: Mapping[str, str]  # SeqRecord silently drops these
    provenance: Provenance
```

`editable` **is the CDS/backbone boundary made structural**; invariant I9 proves it byte-for-byte.

**The rule protocol** carries `id`, `version`, `evidence`, `citations` (**a tuple** with
label/url/year/sign — a single URL makes the badge dishonest on contested rules resting on two sources
with opposite signs), `last_verified`, `weight_provenance`, `parameters` + `param_schema` + `unit` +
`band`, `default_enabled`, `conflicts_with`, `cost_class`, `enforcement_by_modality`, and:

```python
def gate(self, slot: ContextSlot) -> bool: ...
def evaluate(self, c: Construct, ctx: DesignContext, svc: Services) -> Evaluation: ...
def lattice_terms(self, ctx: DesignContext) -> LatticeTerms | None: ...  # opt-in to Tier A
```

`lattice_terms()` decouples the solver lane from every rule lane: a rule ships Tier-B-only on day 3 and
upgrades to Tier-A on day 20 **without the solver agent touching a file.**

**Corrections the judges required — each load-bearing:**

| # | Correction |
|---|---|
| **Enforcement enum** | `HARD_LATTICE / HARD_REPAIR / HARD_CHECK / SOFT / REPORT_ONLY` replaces `Severity{HARD,SOFT,INFO}`, with a `pipeline.py` assertion that a HARD rule can **never** be routed into the weighted sum. The only mechanism that structurally enforces "guaranteed by construction or by repair-plus-independent-validation, never by penalty weight." Three values cannot express uAUG-in-the-user's-5'UTR (real, unfixable) or ENc (must never be an objective). |
| **`Direction` incl. `BAND`** | Raw score in native units plus explicit `unit` and `Direction`. "Lower is always better" cannot represent CAI's 0.70–0.90 band, GC 40–60%, or the two-sided 45–60% AT-window — and a monotone sum drives CAI to 1.0, the exact failure the research refutes. Carries `binding_side` per window. |
| **`RepairPolicy.FIXED_POINT`** | Cryptic splice removal **must iterate to convergence** with an iteration cap. Point-mutating one donor activates cryptic donors nearby (the A2UCOE case, [JVI](https://jvi.asm.org/content/86/9088)); single-pass repair ships a construct whose donors were "removed" into new donors, and the validator passes it because the specific 9-mer is gone. All three proposals got this wrong. |
| **Strand** | Revcomp closure is scoped to `LatticeTerms.forbidden` **only**. Directional scored models — MaxEntScan, the Salis promoter calculator, the TIR/RBS model ([PMC2782888](https://pmc.ncbi.nlm.nih.gov/articles/PMC2782888/)), polyA hexamer+downstream-element — are not revcomp-symmetric and must read `packaged_strand`. Static `strands=(1,)` is AST-banned. For a reverse-oriented lentiviral cassette the packaged genome is the reverse complement, so hard-coding forward makes polyA/splice analysis exactly backwards. |
| **`FoldEnergy` value type** | Carries `dg_kcal_mol` + engine + version + param_set + temperature + dangles so a ΔG never travels without provenance, including into goldens. |
| **`fixable_by_codon_choice`** | On every `Breach`, routing unfixable findings (uAUG in the UTR, TM-segment toxicity, AAV size overflow) to the advisor rather than letting the solver chase them into spurious infeasibility. |
| **`LocalizationPolicy` as data** | `MOTIF_LEN_MINUS_1 / WINDOW_MINUS_1 / WHOLE_SCOPE / PAIRED_SEGMENTS / FIXED_POINT`, consumed generically by Tier B — otherwise ~45 rule agents each hand-roll DnaChisel's localization and each get the overlapping-breach suppression rule wrong independently. |
| **`enforcement_by_modality`** | Internal polyA is HARD on the lentiviral sense strand (it raised expression 3–6.5× but cut functional titer **8–9×** — [PubMed 18627247](https://pubmed.ncbi.nlm.nih.gov/18627247/)) and soft elsewhere; cryptic donors HARD for lentiviral and genome-integrated, warn-only for plasmid; G-quadruplex hard for bacterial CDS, soft for mammalian mRNA (rG4s are globally *unfolded* in mammalian cells, median folding score 0.06 — [PMC5367264](https://pmc.ncbi.nlm.nih.gov/articles/PMC5367264/) — but fold and impair growth in E. coli, where G4 variants span mutation rates 5.5e−5 to 2.7e−10 per cell per generation, [PMC10530614](https://pmc.ncbi.nlm.nih.gov/articles/PMC10530614/)). One frozen level per rule is wrong in every job. |

### The oracle — `verify_construct()`

Hand-written, in `src/` not `tests/`, **import-banned from lane modules** so it cannot share a code
path with the scorer. Called on every `optimize()`.

```
I1 alphabet   I2 frame   I3 round trip   I4 initiator   I5 stops
I6 forbidden motifs on the CIRCULAR construct, pattern set closed under revcomp,
   including junction- and origin-spanning hits, including inside the backbone
I7 GC band global + windowed, windows wrapping the origin
I8 homopolymer / repeat ceiling across the whole construct
I9 every backbone interval is byte-identical to the input backbone      <-- highest value
I10 cassette frame invariant across the assembled CDS
```

**I9 makes "the optimizer silently edited your backbone" a raised exception rather than a shipped
plasmid.** It is the mechanical enforcement of the CDS/backbone boundary.

### Solver

- **Tier A — exact Viterbi DP over the codon lattice.** State = Aho-Corasick automaton over the
  forbidden set ∪ revcomps (+ previous codon for codon-pair terms). 1000 codons × ~600 states ≈ 1.76M
  relaxations, **10–30 ms** NumPy-vectorized. Junction guarantee: seed the automaton with the state
  after consuming the backbone left flank, require a non-accepting state after the right flank.
- **Tier B — DnaChisel-style localized repair.** Calibrated defaults verbatim: exhaustive local search
  below 10,000 local variants, else guided random, 2 mutations/iter, max 1000 iters, stagnation 100
  ([DnaChisel](https://github.com/Edinburgh-Genome-Foundry/DnaChisel)). Circular: triple to 3L, resolve
  constraints overlapping [L, 2L), emit `sequence[L:2L]`.
- **Tier C — structure on windows.** 100-nt window at 10-nt step; a single-codon mutation invalidates
  ~10–20 windows (~2–5 ms/proposal). ViennaRNA ≈0.24 s at 1 kb, ≈6.5 s at 3 kb; do not use LinearFold
  below ~2 kb. **The 5' window must be folded with the real backbone 5'UTR spliced on.** Whole-transcript
  exact structure-aware DP is rejected for the interactive loop — LinearDesign takes 11 min exact for
  1,273 aa and >60 GB RAM at 1,450 aa; DERNA takes 6 h
  ([PMC12319323](https://pmc.ncbi.nlm.nih.gov/articles/PMC12319323/)).

**What is and is not guaranteed by construction.** The Tier-A automaton carries a small summary of
recent bases, so anything decidable from the last few bases is guaranteed absolutely: every forbidden
motif, every homopolymer/G-run limit, every restriction or recombinase site, including ones created
across codon boundaries, on the reverse strand, at the CDS/backbone junction, or spanning the origin.
**Windowed GC *content* is the sole exception** — deciding whether a codon pushes a 50 bp window past
its GC bound requires the G+C count over the previous ~17 codons, and enumerating that history is
combinatorially intractable with no literature precedent. See Q3.

**Three separate reasons GC is constrained, which must not be collapsed into one slider:**

| Concern | What it constrains | Why | Direction |
|---|---|---|---|
| Manufacturability | Windowed and global GC *content* | Extreme local GC breaks phosphoramidite coupling and assembly PCR. Twist's published 50 bp bound is **10–90%**; the widely repeated "35–65%" has no vendor source ([Twist FAQ](https://www.twistbioscience.com/faq/gene-synthesis)). Global target 40–60%. | two-sided band |
| Synthesis + replication | G-runs specifically | IDT allows ~10 A/T but only ~6 G/C consecutively. G-quadruplexes foul solid-phase synthesis and stall replication forks — G4 variants span mutation rates 5.5e−5 to 2.7e−10 per cell per generation ([PMC10530614](https://pmc.ncbi.nlm.nih.gov/articles/PMC10530614/)). | lower |
| Immunology / silencing | **CpG dinucleotides — a different quantity from GC content** | A construct can be GC-rich yet CpG-poor. Three distinct mechanisms: TLR9 sensing ([JCI](https://www.jci.org/articles/view/68205)), ZAP RNA decay (**spacing**-dependent, not count-dependent — [PMC9519448](https://pmc.ncbi.nlm.nih.gov/articles/PMC9519448/)), methylation silencing. | lower |
| E. coli propagation | Windowed AT | Cryptic AT-rich promoters drive toxic leaky expression; toxic horizontally-acquired genes run 63–68% AT vs a 55% non-toxic control ([Nat Microbiol](https://www.nature.com/articles/nmicrobiol2016249)). | **raises GC** |

The last row opposes the first three. The app resolves this as a two-sided band and **displays which
side is binding per window** rather than silently steering.

**Sampling, not argmax.** Forward-filter/backward-sample over the same chain, p(S) ∝ exp(−β·H); β = 1
reproduces natural host CAI/CPB/GC statistics, β→∞ is the argmax
([PMC11345917](https://pmc.ncbi.nlm.nih.gov/articles/PMC11345917/)). Same cost as the Viterbi pass.

**Determinism.** Deterministic DP tie-break; `np.random.default_rng(seed)` threaded everywhere —
**DnaChisel has no seed parameter and uses the global NumPy RNG**, so wrap or reimplement.
`PYTHONHASHSEED=0` in CI plus a grep banning global RNG in `src/`. Content-hash every design; the short
hash goes on the report, the GenBank note, and the CSV `Name` field.

**Infeasibility.** A zero-variant merged segment in the mutation space is a **proof**: emit
`InfeasibilityCertificate` with the minimal conflicting spec set, protein span, and ranked
`Relaxation`s. DnaChisel raises `NoSolutionError` constantly, and "no solution" is a useless product.

---

## Lanes

| Lane | Directory | Owns |
|---|---|---|
| **M1 solver** | `bt5/solver/` | Mutation space, Aho-Corasick, Tier-A DP, β-sampler, Tier-B repair, infeasibility certificates |
| **M2 vector** | `bt5/vector/` | GenBank/SnapGene/FASTA I/O, `IntervalRemapper`, **CDS/backbone partitioning**, insertion-site + UTR detection, construct assembly, `KmerIndex`, Gibson junction uniqueness |
| **M3 score** | `bt5/score/` | Null model, percentile normalization, weights, presets + provenance, λ-sweep + ε-constraint gallery, conflict detection, QC report, design hash, **GenBank + order-CSV export** |
| **M4 rules** | `bt5/rules/catalog/` | ~45 rules, **one file per rule**, `<brief_id>_<slug>.py` so "did we implement D5(c)?" is an `ls`. Includes vendor complexity pre-checks. |
| **M5 codon** | `bt5/codon/` | Genetic code tables, host usage, CAI/tAI/stAI/%MinMax/CFD/ENc/CSC, harmonization, user-supplied tables |
| **M6 structure** | `bt5/structure/` | ViennaRNA `FoldEngine`, windowed ΔG with incremental invalidation, **UTR-aware 5' terms**, IVT terms, DegScore |
| **M7 packaging** | `packaging/` | **Lane zero.** `uv` install, ViennaRNA, `commec` DB flow; later the Tauri sidecar + macOS notarization |
| **M8 cassette** | `bt5/cassette/` | Protein validator, genetic-code selection, tag/linker/2A/protease library, repetitive-protein analysis, protein liability scan on the assembled fusion, biosecurity screening |
| **M9 server** | `packages/server/` | FastAPI on 127.0.0.1, `ProcessPoolExecutor` jobs, cancellation, SSE progress |
| **M10 web** | `apps/web/` | React/TS, sequence viewer, sliders, gallery, conflict panel, CSV download |

### Why nobody waits

1. **Registry autodiscovery** via `pkgutil.walk_packages`. Adding a rule edits **zero** shared files.
   No committed catalog file — both losing proposals reintroduced one and it collides on every PR.
2. **`lattice_terms()`** decouples solver ⊥ rules.
3. **Services injection** decouples rules from tables, folding and k-mer indexing.
4. **Entry-point provider resolution** — a lane goes stub → real by editing one line in its own
   directory. All three proposals left `pipeline.py` owner-serialized; that is a guaranteed bottleneck.
5. **Manifest-driven UI.** `GET /capabilities` serves the catalogue from the live registry.
6. **Everything shared is pre-registered in PR #0**: every dependency (no lane adds one — `uv.lock` is
   the least-mergeable file in the repo); every benchmark metric name at `blocking: false`; the `data/`
   line **split** (protein-changing data gated in root `data/`, liability-flagging data in-lane with a
   `_provenance.json` sidecar); **per-lane golden directories**; fixtures namespaced by consuming lane;
   **per-lane CI path filters**.
7. **`ReferenceSolver`** — deliberately dumb, exhaustive over ≤12-codon windows, slow, never wrong.
   Unblocks every rules lane from minute one *and* survives as the permanent differential oracle.
8. **Amendment protocol.** MINOR on a fast path; MAJOR requires an RFC plus a deprecation shim with a
   `model_validator(mode="before")`, a two-window rule, and `test_backward_compat` re-parsing every
   recorded fixture inside the amendment PR. `contract-freeze` is a required CI job.

---

## Build sequence

**PR #0 (Contract):** `bt5/core/**`; `verify_construct.py`; all `.github/**`; `CLAUDE.md`; labels;
`data/**` + `MANIFEST.sha256`; `benchmarks/{panel,tolerances}`; `tests/invariants/**`; `pyproject.toml`
with every dependency; `uv.lock`; `brief_index.toml`; 2–3 hand-written **reference rule files**;
`ReferenceSolver`; `docs/research/**`. Pushed to `main` **before the rulesets exist**.

**PR #1 (Skeleton):** one host (E. coli K-12, table 11), three rules, one real Addgene lentiviral
GenBank **with an annotated 5'UTR and an origin-spanning feature**, ViennaRNA behind `FoldEngine`, one
preset, real λ-sweep, real 200-variant null, one report, one UI page — protein → validated → screened →
CDS planned → spliced into the circular backbone → mutation space over the CDS only → Tier-A DP →
Tier-B repair → independent `verify_construct` → normalized scorecard → 5-candidate gallery → annotated
GenBank + order CSV, **under 10 s.** Doubles as the ruleset rehearsal. ~5–7k hand-written lines.

### Week-one gates — each is a number

| Gate | Measures | Pass | Fail ⇒ |
|---|---|---|---|
| **G0** | `Spec`/`Evaluation` carries localized, attributable, context-gated circular results. **Validated against one rule from each of the four hardest families** — fixed-point splice removal, Salis TIR over a window crossing into the backbone 5'UTR, harmonization, and a two-sided band. | 0 additions to `Evaluation`/`Breach` | protocol is wrong; iterate before wave 1 |
| **G1** | Tier-A DP, 1000 codons, ~600 states | ≤50 ms | >500 ms ⇒ drop Tier A, re-plan around repair-only, say so in the README |
| **G2** | Addgene lentiviral GenBank round-trips byte-identically incl. origin-spanning `join()`; forbidden sites planted across the CDS/backbone junction **and** across the origin both caught; I9 catches a deliberate backbone edit | all three | the coordinate model is wrong |
| **G3** | 200-variant windowed-fold null on the assembled construct, 500 aa | <2 s | normalization becomes report-time-only, recorded in `Provenance.degradations` |
| **G4** | Dense λ-sweep + ε-constraint → greedy max-min pick of 5 | pairwise codon distance ≥15% | ε-constraint enumeration with explicit diversity constraints becomes primary, *before* M10 builds a UI on it |
| **G5** | `junction_trap` adversarial fixture | raises with ≤3 spec ids, interval ≤10 codons | redesign the mutation-space merge |
| **G6** | 5' ΔG computed with the real backbone 5'UTR differs measurably from CDS-alone, and the no-UTR path degrades rather than reporting a number | both | the UTR coupling is not wired; the highest-weight objective is wrong |
| **G7** | End-to-end skeleton, 500 aa | ≤10 s | re-allocate budget before rules multiply cost |
| **G8** | `uv` install on a clean machine incl. ViennaRNA and the <1 GB commec biorisk DB | clean install | packaging is the real critical path |

**G2 and G4 should scare you most.** G2 because a wrong coordinate model invalidates everything built
on it; G4 because its failure invalidates a *product* decision, not a technical one.

### Waves

**Wave 1 — 3 agents:** M1 solver, M2 vector, M3 score — the only lanes that can still force a `core/`
change. Three, not eight, because a `core/` change costs three rebases rather than ten. `core/` freezes
for real when wave 1 lands. **Watch for:** any amendment changing the *shape* of `Evaluation`/`Breach`
rather than adding a field. Two in a week means G0 was insufficient and the freeze should slip.

**Wave 2:** M4, M5, M6, M8, M9, M10 against a frozen core. M7 runs throughout. **If G2 fails**, M4/M6/M8
are largely coordinate-independent — start those against recorded fixtures and hold M1/M2/M3.

**Concurrency — RATIFIED: cap at 5 open PRs, 10 lanes rotating.** A free personal account gets 20
concurrent job slots; a Python PR consumes ~12, so 8 open PRs would be ~4–5× oversubscribed and agents
would queue behind each other with no way to tell "slow" from "hung". Ten lanes remain as ten *ownership
boundaries* — that is what keeps merge conflicts near zero — but only **5 non-draft PRs are open at
once**. Agents work in drafts (which skip `e2e` and `benchmark-gate`) and mark ready when they believe
they are done. No infra change and no cost. Cut order if it still bites: CodeQL off PRs →
benchmark-gate to ready-only → transfer to a free org for merge queue. **Do not turn `strict` on** —
that is a livelock worse than the conflicts it prevents.

---

## GitHub setup

Public repo, personal free account, `main` default, admin rights: unlimited Actions minutes and full
rulesets, but **no merge queue** (organization-only).

1. **Required approvals = 0.** GitHub forbids self-approval and your agents authenticate as you, so any
   non-zero count is a permanent deadlock. Sign-off is re-implemented as a status check.
2. **One aggregating required check per workflow.** A path-filtered workflow that never triggers never
   reports, and a required check that never reports blocks the PR forever with no error.
3. **Loose status checks (`strict: false`).** With N parallel PRs, strict is O(N²) CI.

**Settings:** squash-merge only; auto-delete head branches; allow auto-merge. Actions workflow
permissions **read-only**, and **"Allow GitHub Actions to create and approve pull requests" UNCHECKED**
— with it on, an agent that can edit `.github/workflows` can approve its own PR. Secret scanning +
**push protection** (separate toggle, off by default). CodeQL **advanced setup only** — default setup
silently disables a committed workflow. Zero classic branch protection rules.

**Ruleset `main-protection`:** squash-only, no force-push, no deletion, linear history, required
conversation resolution, `required_approving_review_count: 0`, `strict: false`, required check pinned to
the GitHub Actions app id, and **`bypass_actors` empty, including you** — rulesets do *not* implicitly
exempt admins the way legacy branch protection did. Emergency escape is *not* a bypass actor (bypass is
per-*actor*, so it disables every rule): set enforcement to `disabled` for five minutes and flip back.

**Solo-owner sign-off:** `required_approving_review_count: 1`, `require_last_push_approval` and
`require_code_owner_review` are all traps that make the repo mathematically unmergeable. Instead: CI as
the machine gate, required conversation resolution as the checklist, and
`anthropics/claude-code-action/agent-approval-check` as a status check accepting `/approve <head-sha>`
from the PR author. Land it **not required** for a few days first.

### 🔴 The benchmark gate ships scientifically backwards — fix in PR #0

The base CI plan ships `cai: direction: higher_is_better, blocking: true` plus an invariant
`test_cai_beats_random_backtranslation`. These are **blocking required checks**, and the evidence above
establishes CAI has no predictive value and must be a *band*. As written, every agent optimizing toward
green CI is trained by the merge gate itself to maximize the one metric the science refutes — and any PR
correctly moving CAI into a band reads as a regression and fails. **Change `cai` to `direction: band`
and delete the CAI-beats-random invariant before any lane starts.**

### Other CI traps

- `PYTHONHASHSEED: "0"`; grep banning global RNG in `src/`.
- **CAI is meaningless without a pinned reference set.** Version-pin it in the baseline, hard-fail when
  it moves, and track a **held-out** set — if optimized CAI rises but held-out CAI does not, the
  optimizer is fitting the table rather than improving the sequence.
- **Never put ΔG in a byte-exact snapshot.** One runner architecture, compare with tolerance.
- **ViennaRNA version policy — stay current without breaking comparability.** Pinning is for
  correctness, not conservatism: energy parameters determine every ΔG, so a silent version change
  shifts every baseline and makes old results incomparable to new ones. So: pin exactly
  (`viennarna==2.7.2`, current as of 2025-12-30; cadence is roughly annual — 2.7.0 Oct 2024, 2.7.1 and
  2.7.2 Dec 2025); let Renovate open a PR the day a new version ships so you are never more than one
  PR behind; treat that PR as a **scientific change, not a dependency bump** (label
  `approved:algorithm-change`, regenerate the benchmark baseline); and have the comparability guard
  **hard-fail any benchmark comparison across differing `viennarna_version`** so nobody compares
  numbers taken with two different rulers. `FoldEnergy` carries `engine_version` and `param_set`, so
  results stay interpretable after an upgrade. Note the software version and the energy parameter set
  (Turner 2004) are separate — most releases do not change ΔG at all, and the guard reports which did.
- Biopython's `Seq.translate()` **silently truncates** a non-multiple-of-3 sequence (only a
  `BiopythonWarning`), so a frame-length bug passes a naive round-trip test. Check length *before*
  translating; run pytest with `-W error::BiopythonWarning`.
- Forbidden motifs are created **across codon boundaries and on the reverse strand** — the two bug
  classes AI agents produce most reliably, both invisible to any per-codon test.
- Banned-schema-field grep against the **generated** OpenAPI.

### The eight things that silently disable this entire system

1. `paths:` in `on:` on any workflow owning a required check → permanent deadlock, no error.
2. `if: success()` or no `if:` on the gate job → gate **skipped** when a dependency fails, and skipped
   **satisfies** a required check. Failing tests merge.
3. `if: always()` with no `needs.*.result` inspection → gate reports success while everything burns.
4. A new job not added to `required-checks.needs` → runs, can fail, unenforced.
5. Yourself in `bypass_actors` → every rule off for the only person who merges.
6. `required_approving_review_count: 1` or `require_last_push_approval` → unmergeable, forcing (5).
7. Required check left as "any source" instead of pinned to the Actions app → forgeable green.
8. CodeQL default setup enabled alongside the advanced workflow → silently disables it.

### Standing rules for CLAUDE.md

Stay in your lane. Never hand-resolve a lockfile — rebase and regenerate. **Hypothesis counterexample on
an unrelated lane's PR:** file it as a `tests/fixtures/regressions/` fixture plus an issue; if it
reproduces on the merge base, the owner merges that PR manually — **never weaken the property**.
Suppression is not a fix. `--snapshot-update` is not a fix.

---

## Verification

**Never validate with the same function you optimized with.**

1. **Invariants**: `translate(output, declared_table) == input_protein` for every host × mode × seed;
   frame preserved across the assembled CDS; **zero forbidden motifs in the final assembled plasmid,
   checked by the import-banned independent validator**; GenBank round-trip; I9 backbone byte-identity.
2. **Metamorphic**: a redundant constraint must not change output; permuting constraint order must not
   change output (**DnaChisel fails this** — a ready-made bug class); re-optimizing an optimal sequence
   is a no-op. **Differential** against `ReferenceSolver` (permanent) and DnaChisel (nightly).
3. **Retrospective benchmarking** on public *measured* data — Kudla (154 GFP variants),
   [Goodman 2013](https://europepmc.org/abstract/MED/24072823) (~14,000 reporters),
   Cambray (244,000), [Verma 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6920384/) (259,134),
   Boël (6,348). **Hold out by gene, not by variant** (variant splits leak). Report Spearman ρ **and
   top-k enrichment** — "of the 10 sequences my objective ranks highest, what fraction are in the
   measured top 10%?" is the decision the user makes; ρ is not.
4. **Negative controls the field mostly fails**: beat CAI, and beat a random synonymous shuffle. The
   same 200–500 null variants double as this control at zero extra cost.
5. **Adversarial**: poly-Q, 100% Leu, a 2-residue protein, FVIII (2332 aa) and dystrophin (3685 aa) so
   graceful degradation is exercised in CI, a protein whose every synonymous option creates a forbidden
   site, a backbone with zero annotations, a backbone with **no annotated 5'UTR**, an origin-spanning
   feature.
6. **Vendor-complexity oracle** — paste designs into the Twist/IDT complexity checkers and record
   accept/complex/reject. **Manual periodic calibration, not CI**: automating authenticated commercial
   order forms is scraping against ToS and breaks silently on any form change.

**State plainly what cannot be validated: absolute expression prediction.** Given the ~14% ceiling
(Cambray), any UI implying a predicted titer is dishonest.

---

## Open questions — resolutions

| # | Question | Resolution |
|---|---|---|
| **1** | Does the candidate gallery diversify? | **Design resolved; G4 confirms.** [Das & Dennis 1997](https://link.springer.com/article/10.1007/BF01197559) proves weighted sums have two failure modes: non-convex Pareto regions are unreachable, *and* an even spread of weights gives an uneven, clustered spread of points because weight tracks the curve's slope. **So never select the gallery by evenly-spaced weights** — sweep λ densely, then pick 3–8 by greedy max-min distance in *sequence* space. ε-constraint is promoted from fallback to co-primary; it also *accelerated* search 1.25–2.21× via pruning ([PMC7358382](https://pmc.ncbi.nlm.nih.gov/articles/PMC7358382/)). |
| **2** | Is the empirical null affordable on the assembled construct? | **Resolved by design.** Naive is fatal: ViennaRNA ~0.24 s/kb × 200 variants ≈ 48 s. But the highest-weight term is a *windowed* ~90 nt fold, not a full-transcript fold. **The null uses windowed folding only**; whole-transcript MFE becomes report-time-only with a smaller (n≈50) null. G3 confirms <2 s. |
| **3** | Windowed GC guarantee | **RATIFIED: steer → repair → refuse.** Everything decidable from the last few bases — forbidden motifs, restriction and recombinase sites, homopolymer and G-run limits — is guaranteed absolutely by the Tier-A automaton, including across codon boundaries, on the reverse strand, at the CDS/backbone junction and spanning the origin. Windowed GC *content* is the sole exception, because it needs the G+C count over the previous ~17 codons and enumerating that history is intractable (no literature precedent). Block-state DP (600 × 52 = 31.2k states, ~92M relaxations, 0.5–1.5 s vs ~10–30 ms) was rejected as the default because *block-wise* GC is a strictly weaker surrogate — a window straddling two blocks can still violate — so it costs ~50× for a guarantee weaker than the actual rule. **Implementation:** Lagrangian steering term in Tier A, `HARD_REPAIR` in Tier B, and the independent validator **refuses to emit** if any window is still out of band, reporting the minimal conflicting set instead. A bad sequence never reaches the user; worst case the app declines and explains. Block-state stays available behind an opt-in flag for a final design. |
| **4** | Backbone repeats vs plasmid yield; 3'UTR beyond WPRE | **Materially advanced — and my earlier rule was too crude.** Two thresholds were conflated: RecBCD MEPS (23–27 bp) is the floor for RecA-*dependent* recombination, but below ~200 bp deletion is RecA-*independent* (slipped-strand / single-strand annealing) and unaffected by `recA⁻`, with recA-dependence rising only above ~300 bp ([Springer](https://link.springer.com/article/10.1007/BF00290109), [PMC5426353](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5426353/), [PNAS](https://www.pnas.org/doi/10.1073/pnas.111008398)). **Proximity is a first-class variable** — recA-independent recombination is strongly distance-sensitive, and inserting sequence between repeats suppresses it. **So repeat risk is a 2-D surface over (length, spacer), not a length threshold.** 3'UTR beyond WPRE: no clean answer found; ships `contested`, weight 0, non-blocking. |
| **4a** | *You noted `recA⁻` is standard for LVV and AAV plasmid cloning — which sharpens this considerably.* | **`recA⁻` protects the long repeats and does essentially nothing for the short ones.** Stbl3 (`recA13`) / NEB Stable are the right strains and they matter: Stbl2→Stbl3 alone rescued an HIV vector that was lost entirely in 0.5 L Stbl2 cultures ([PMC3563744](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3563744/)). But that protection is against the RecA-*dependent* pathway, which needs >~200–300 bp homology — i.e. the LTRs and ITRs themselves. **The repeats a codon optimizer creates or removes are 15–100 bp, squarely inside the RecA-independent regime the strain does not suppress.** Measured directly: a 28 bp repeat pair still recombined at 7.8e−7 to 3.1e−5 in **four different `recA⁻` strains** ([Oliveira 2008](https://www.genoscope.cns.fr/MGE/pubs/Oliveira_Mol_Biotechnol_2008.pdf)). **Consequence for BT5:** short-repeat avoidance is the app's job and cannot be delegated to the strain, so the repeat rules carry high weight by default in the LVV/AAV presets — and the report must not tell a user that `recA⁻` covers them. The strain recommendation stays in the protocol output (Stbl3 at 30 °C for LTRs; ΔsbcC at 42 °C for ITR palindromes) but is presented as covering the *long* repeats only. |
| **5** | Boël vs Kudla past codon ~17 | **Not resolvable.** [Boël 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5054687/) (6,348 genes, T7, soluble-protein readout) finds codon content 3–5× more influential; Kudla/Goodman/Cambray (host RNAP, single-gene libraries, fluorescence) find structure dominant. Unresolved whether the discrepancy is T7 vs σ70, many-genes vs many-variants, or solubility vs fluorescence. Ship both as separately-weighted hypotheses with `weight_provenance` naming the disagreement; the calibration view lets a user resolve it in *their* system. A limit to display, not a bug to fix. |
| **6** | CAI convention | **Resolved.** [Sharp & Li 1987](https://pubmed.ncbi.nlm.nih.gov/3547335/) defines CAI as the geometric mean of relative adaptiveness **excluding start and stop codons**. Pin that: exclude ATG/TGG/stops, pseudocount 0.5 for zero-count codons, reference set from highly-expressed genes ([PaxDb](https://pax-db.org/), CC BY 4.0). Ship the alternative convention as a toggle, print which was used, and keep the Sharp & Li E. coli table as a regression fixture — **copy it out of Biopython 1.81 now, it was removed in 1.82**. |
| **7** | %MinMax window | **Resolved by decision.** Pin 18 (CodonTransformer's value; CHARMING defaults to 10, other sources say 17). Expose in `param_schema`, print in the report. |
| **8** | IDT/Twist order-file headers | **Resolved — you supplied the template.** IDT eBlocks 96-well plate upload is `.xlsx` with a single sheet whose **name is the plate name**, three columns (`Well Position`, `Name`, `Sequence`), and 96 pre-filled wells in row-major order `A1…A12, B1…B12, … H1…H12`. See "The order file" above. Twist still needs a template confirmed at implementation; the vendor-profile mechanism covers it. |
| **9** | Patent FTO on codon-lattice × folding DP | **Off the critical path.** Baidu's LinearDesign patent filing is real, as is [CN117660445A](https://patents.google.com/patent/CN117660445A/en) on adjacent-codon optimization. **But you specified personal use only**, so this is not a build-time blocker — it becomes a gate on any future decision to distribute or commercialize. I previously put FTO review in week one; that was wrong for your context. Tier C stays post-hoc windowed polish anyway, for performance reasons. |
| **10** | `commec` database size | **Resolved, and better than I said.** The biorisk-only tier is **under 1 GB and runs on a laptop**; only the full protein/nucleotide similarity search needs ~600 GB ([commec-databases](https://github.com/ibbis-bio/commec-databases/), [IBBIS FAQ](https://ibbis.bio/our-work/common-mechanism/faq/)). So biosecurity screening ships **on by default** with the <1 GB HMM biorisk tier, and the 600 GB full-homology tier is an explicit opt-in. My earlier "multi-GB, possibly prohibitive" was wrong in the direction that matters. |

### What could still kill this

1. **The interface is wrong and G0 passes anyway.** Mitigated by scoping G0 to the four hardest rule
   families, but a mode nobody has implemented can still surprise.
2. **Wave 0 becomes Wave 0-through-3.** ~5–7k owner-serial lines plus eight gates, with no branch for
   "G2 fails." Mitigation: start the coordinate-independent lanes against recorded fixtures.
3. **"Use the native sequence" is right often enough that the product feels pointless.** This is the
   trust differentiator, but the value must then rest on the **mechanical** half — manufacturability,
   propagation stability, junction-spanning guarantees, whole-CDS back-translation, the V5 donor, repeat
   divergence. **Watch for:** any week where more rule-lane effort goes into expression families than
   into motif/manufacturability/stability families.
4. **Concurrency starvation makes ten lanes theatre.** Five open non-draft PRs maximum.
5. **Scientific correctness rots quietly as the rule count grows.** The contract test catches *missing*
   provenance; it cannot catch *wrong* provenance. Mitigation is `last_verified`, quarterly vendor
   recalibration, and the owner reading five rule files a week against their citations. The one control
   that is a human habit rather than a CI job, and the one most likely to lapse.
6. **Packaging is the real critical path and nobody notices until week eight.** G8 exists to surface it
   and is the gate most likely to be skipped because it feels like plumbing. It is not plumbing.
