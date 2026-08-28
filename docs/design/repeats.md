# Repeats — what we are actually protecting against, and what the user gets to choose

Status: design proposal. Research completed 2026-08-28.
Supersedes nothing; refines brief rows 2.E5–E8 and 2.F1–F3, all of which are
already implemented. Read `docs/research/brief.md` §2.E and §2.F first.

BT5 already has eight repeat-adjacent rules. They were each built against a
correct local reading of the brief, but they were built one at a time, and the
brief states thresholds without stating **which failure each one is protecting
against**. This document supplies that, because the research turned up two
things that the per-rule view cannot see: the distance term has **opposite
signs** in two of the four failure modes, and the length ranges of those two
modes **do not overlap**.

---

## 1. Four failures, not one

"Repeat" names four mechanically unrelated failures. They differ in the length
that matters, in whether distance helps or hurts, and — decisively — in whether
codon choice controls them at all.

| # | Stage | Mechanism | Length that matters | Distance | Codon-controllable? |
|---|---|---|---|---|---|
| 1 | **Synthesis** | oligo mispriming, assembly-PCR misannealing | **8–100 bp**, clusters of 8–9 bp count | clustering matters, absolute distance does not | **Yes** |
| 2 | **Plasmid propagation** | RecA-**independent** SSA / slipped-strand | **~14–100 bp** (sole contributor below ~100 bp) | **closer = worse** | **Yes** |
| 2b | Plasmid propagation | RecA-**dependent** homologous recombination | **>200–300 bp** | farther = more recA-dependent | No — LTR/ITR scale |
| 3 | **Packaging / reverse transcription** | RT template switching between the two co-packaged genomes | **≥~114 bp, steep above ~350 bp** | **farther = worse** | Only at architectural scale |
| 4 | **AAV plasmid prep** | ITR palindrome collapse during bacterial amplification | ITR scale | — | No — ITRs are `WHITELISTED_REPEAT` |

Rows 2 and 3 are the ones that matter, and they disagree with each other.

### 1a. The sign reversal

In the plasmid, RecA-independent deletion is **proximity-sensitive in the
ordinary direction**: putting sequence between two repeats suppresses it. Bi &
Liu measured this directly — "increasing the distance separating the homologous
regions preferentially inhibits the recA-independent recombination", while
"shortening of the homology preferentially inhibits recA-dependent
recombination"
([J Mol Biol 235:414–423](https://pubmed.ncbi.nlm.nih.gov/8289271/)).

In the packaged lentiviral genome the relationship **inverts**. Reverse
transcriptase deletes the sequence between direct repeats by template switching,
and the further apart the repeats, the more reliably it does so: for a 701 bp
direct repeat, deletion ran 37% at 799 bp of separation, 82% at 1,407 bp, and
**>90% at every distance above ~1,545 bp**, out to 4,175 bp
([J Virol 73:7923–7932](https://journals.asm.org/doi/10.1128/jvi.73.10.7923-7932.1999)).
Intermolecular template switching between the two co-packaged genomes happens
**~3–30 times per infection** — more often than base substitution.

`f2_near_perfect_repeats.py` currently implements one decay, `SPACER_DECAY_BP =
3000`, in the plasmid direction, and applies it in every modality. **On a
lentiviral construct that term has the wrong sign.**

### 1b. …but the lengths do not overlap, which defuses it

The same J Virol study tested shorter repeats: 114 bp, 225 bp and 349 bp gave
only **6–9% deletion regardless of distance**. Template switching needs
substantial 3' homology before it becomes efficient.

Codon choice controls repeats of roughly **15–100 bp**. That is entirely below
the RT template-switching regime and entirely inside the RecA-independent
plasmid regime.

**So the practical conclusion is not "add a modality-dependent sign to F2".** It
is:

> Codon-level repeat control protects **synthesis** and **plasmid propagation**.
> It does **not** protect the packaged viral genome, and BT5 must not imply that
> it does.

This matters for how the viral presets are justified. Today `_REPEAT_NOTE` in
`bt5/score/presets.py` justifies the elevated viral repeat weight by the
recA⁻/short-repeat argument — which is correct, and is a **propagation**
argument. It should not acquire a packaging argument it cannot support.

### 1c. Where the two regimes DO intersect — and it is a real case

There is one overlap, and it is exactly the case `docs/PLAN.md` already names as
a differentiator: **architectural repeats inside the ORF**. Two copies of the
same 2A peptide, a duplicated binding domain, a tandem scFv, a repeated tag —
these are 100–700 bp, which is squarely in the template-switching regime, *and*
they are inside the CDS, so whole-CDS back-translation can diverge them.

The documented case is a multi-shRNA lentiviral vector whose cassettes were
~300 bp of which **250 bp (83%) was repeated**: deletion ran 2–36% across
combinations and **100% for the 4-shRNA construct**, and the authors placed the
event **after viral production, during or after transduction** — RT
recombination in the transduced cell, not plasmid instability
([PMC2775741](https://pmc.ncbi.nlm.nih.gov/articles/PMC2775741/)).

That is the finding that earns BT5 a rule the competition does not have. It is
also the one place where the inverted distance sign is real, because these
repeats are long enough for template switching to work.

---

## 2. Vendor limits: the published numbers are not the acceptance boundary

Verified against live vendor pages on **2026-08-28**.

**Twist**, high-complexity scoring triggers
([FAQ](https://www.twistbioscience.com/faq/gene-synthesis/what-do-scoring-results-my-gene-mean)):

| Trigger | Threshold |
|---|---|
| Direct repeats | **>200 bp** |
| Short tandem repeats | **>100 bp** |
| Hairpin-forming repeats | **>100 bp** |
| Homopolymer (any single-nucleotide run) | **>30 bp** |
| Local GC in a 50 bp window | **<10% or >90%** |

Hard rules: homopolymers ≥14 bp, no CcdB. Recommended GC 25–65%. Separately,
Twist's design guidance recommends **avoiding direct repeats longer than 12–16
bp** and **avoiding clusters of short repeats at 8–9 bp**.

**IDT**: GC below 25% or above 75% is problematic; homopolymers of ≥10 A/T or
≥6 G/C; and IDT publishes a numeric **complexity score with a hard cutoff — any
sequence scoring above 10 cannot be synthesised**
([IDT gBlocks FAQ](https://www.idtdna.com/pages/support/faqs/what-types-of-sequence-motifs-should-be-avoided-when-ordering-gblocks-gene-fragments-)).

### The gap that matters

Read literally, Twist's published trigger tolerates a **200 bp** direct repeat.
Actual behaviour is nothing like that. Eight constructs whose repeated elements
ran **20 bp to 81 bp** were *all* rejected by IDT, with complexity scores of
**30 to 139** against a cutoff of 10 — a 133 bp construct with a 2× repeat
scored 30.2, and a 280 bp construct with an 8× repeat scored 139.2. Several also
failed Twist's screen. The authors had to abandon commercial synthesis and
assemble from ≤80 bp synthons by Golden Gate
([PMC10949351](https://pmc.ncbi.nlm.nih.gov/articles/PMC10949351/)).

**The published thresholds are labels on a complexity tier, not the accept/reject
boundary.** A rule that encodes Twist's 200 bp figure as "the limit" would pass
designs that get rejected at order time by an order of magnitude.

This vindicates the thresholds already in `e5_synthesis_repeats.py` (flag at 12
bp, hard at 20 bp, severe at 200 bp) and means the 200 bp figure should keep its
current role — the *severe* tier, not the constraint. It also means:

- **`vendor_asserted` is too weak a badge for the 12–20 bp numbers.** They are
  now backed by measured vendor *outcomes* (PMC10949351, and the 1,076-outcome
  random forest already cited at 2.E6), not just vendor copy. The badge should
  say which of the two a number rests on.
- **The IDT complexity score is a published, thresholded, vendor-authoritative
  scalar.** BT5 cannot compute it (proprietary), but it is the right shape for
  the manual calibration loop the plan already schedules under
  "vendor-complexity oracle".

---

## 3. The quantitative model we should adopt

`f1_direct_repeats.py` and `f2_near_perfect_repeats.py` currently implement the
(length × spacer) risk surface with hand-chosen tiers — `INFO_BP = 15`,
`WARN_BP = 20`, `HARD_BP = 25`, and an exponential spacer decay with a
3 kb constant. Those tiers are defensible individually; the surface joining them
is invented.

It does not need to be. Oliveira et al. published **a fitted non-linear function
predicting recombination frequency from repeat length and intervening-sequence
length jointly**, built on collected deletion data covering direct repeats of
**14–856 bp** and intervening sequences of **0–3,872 bp**, and parameterised by
the strain's recA genotype. More than **92% of predictions fall within ±5-fold**
of experiment
([Plasmid 60:159–165](https://pubmed.ncbi.nlm.nih.gov/18647618/)).

That is the exact object `docs/PLAN.md` Q4 asks for — "repeat risk is a 2-D
surface over (length, spacer), not a length threshold" — already fitted,
already validated, and already covering the recA⁻ case that the LVV/AAV workflow
runs in.

Two properties make it a good fit for BT5 specifically:

1. It is **parameterised by recA genotype**, so the "recA⁻ protects the long
   repeats and does nothing for the short ones" claim stops being a narrative in
   `_REPEAT_NOTE` and becomes a computed difference between two curves.
2. Its accuracy is **±5-fold**, which is honest and which BT5 is already built to
   express. A ±5-fold band renders as a confidence band and a rank, never as a
   number — the same discipline the expression objectives are held to.

**Caveat, and it is a real one:** the functional form is behind a paywall and I
could not extract it. Adopting this is a task with a research step, not a
copy-paste. See work item R3.

---

## 4. What the user gets to choose

The ask: the user should be able to say what repeat length they are comfortable
with in the sequence BT5 hands them.

### The parameter

```
max_repeat_bp: int   # no exact repeat >= this length anywhere in the
                     # assembled construct, either strand, including
                     # backbone junctions and the origin
```

- **Default 20.** It is where two vendors' published hard-fail sits, it is below
  the RecBCD MEPS floor of 23–27 bp, and it is above the point where the
  measured rejection data starts.
- **Scope: the assembled construct.** That is the molecule the user ends up
  holding. BT5 still reports the *fragment* scope separately, because the vendor
  screens the fragment plus its adapters and will reject on that regardless of
  what the assembled plasmid looks like. The two scopes already exist and
  already disagree deliberately (`bt5/rules/fragment.py`); the knob must not
  collapse them.
- **Enforcement: `HARD_REPAIR`.** Not `HARD_LATTICE` — k-mer uniqueness is not
  decidable from a bounded suffix, since it depends on the entire prefix, so the
  Aho-Corasick automaton cannot guarantee it. Repair plus the independent
  validator, which refuses to emit, is the mechanism (§3.5 of `CLAUDE.md`).

### The part that is not a slider

A repeat-length floor is not freely choosable, and a UI that presents it as one
is lying. A `(GGGGS)₃` linker is 45 bp drawn from two amino acids; poly-Gly has
four codons and poly-Lys has two. Below some length **this protein cannot be
back-translated at all** under this genetic code, and the honest response is to
say so before the user picks, not to hand them an `InfeasibilityCertificate`
after.

So the knob ships with a **computed floor**:

1. Compute the protein's own repeat structure *before* codon choice — which
   `bt5/cassette/` already does for the repetitive-protein rule.
2. From that, derive the minimum achievable `max_repeat_bp` for this protein
   under this table.
3. Render the slider with that floor as a hard stop, annotated with **which
   protein feature sets it** ("the (GGGGS)₃ linker at residues 118–132").

That last string is the whole value. "You cannot go below 26 bp" is a wall;
"you cannot go below 26 bp *because of the GGGGS linker*" tells the user their
options are to accept it, shorten the linker, or use a different one — a protein
design decision BT5 is uniquely positioned to surface, since it is the only tool
that looks at the assembled construct and the protein together.

### And it must show what it costs

Tightening `max_repeat_bp` spends the sequence's freedom. That freedom is what
every other objective is competing for, so the knob needs the same treatment as
every other trade-off in BT5: as the user drags it down, show the percentile
cost to the other objectives, from the same λ-sweep machinery
(`bt5/score/gallery.py`) the candidate gallery already uses.

---

## 5. Work items

Ordered. Lane owner in brackets; cross-lane items need an issue first per
`CLAUDE.md` §1.

| # | Item | Lane |
|---|---|---|
| **R1** | Give `f2_near_perfect_repeats.py` a modality-aware distance term. Plasmid sign (closer = worse) stays the default; the packaged-genome sign inverts. Gate the inverted term on repeat length ≥100 bp so it only fires where template switching is real. | M4 |
| **R2** | New rule: **in-ORF architectural repeats**, ≥100 bp at ≥90% identity within the CDS, read on the packaged strand. This is the duplicate-2A / tandem-domain case, the one place codon choice reaches the RT regime. Report-only until R3 lands a calibrated risk number. | M4 |
| **R3** | Obtain Oliveira et al. 2008 (*Plasmid* 60:159–165), extract the fitted function, and replace F1/F2's invented surface. If the paper cannot be obtained, refit against the same published deletion dataset and say so in `weight_provenance`. Carries `approved:algorithm-change`. | M4 + M3 |
| **R4** | `max_repeat_bp` parameter, `HARD_REPAIR`, assembled-construct scope, default 20. | M4 + M1 |
| **R5** | Computed feasibility floor for `max_repeat_bp`, with the protein feature that sets it named in the response. | M8 |
| **R6** | Slider + floor + cost curve in the UI. Depends on R4, R5. | M10 |
| **R7** | Correct `_REPEAT_NOTE` in `bt5/score/presets.py` so the viral repeat weight is justified by propagation and synthesis only, and state plainly that codon-level repeat control does not protect the packaged genome. | M3 |
| **R8** | Split the `vendor_asserted` badge so a number backed by measured vendor *outcomes* is distinguishable from one backed by vendor *copy*. The 12–20 bp synthesis thresholds are the former. | M4 + core |

### Verification debt

Several numbers above reached me through search summaries and page extraction
rather than my reading the primary text — specifically the Bi & Liu 1994
distance/homology asymmetry, the "RecA-independent is the sole contributor below
100 bp" claim, and the Oliveira function's form and coefficients. They are good
enough to plan against and **not good enough to become rule constants**. Pull the
PDFs before any of R1, R2 or R3 hard-codes a number, and record `last_verified`
per `Spec.citations` when you do.

The Twist and IDT thresholds in §2, and every figure from PMC10949351, PMC2775741
and J Virol 73:7923, were read from the sources directly and are quotable as-is.
