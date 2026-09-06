## 2026-09-06 — D5 carries two enforcement mechanisms, one per sub-rule shape

**Decided:** `d5_cryptic_transcription` ships with `enforcement = HARD_LATTICE` as its
ClassVar and `enforcement_for(slot) -> HARD_REPAIR` for every slot, because D5's five
sub-rules split cleanly into two shapes and BT5 gives a `Spec` one of each switch:

| Part (`brief.md:110-112`) | Shape | Mechanism |
|---|---|---|
| (b) extended −10 `TGnTATAAT` | fixed IUPAC string | `LatticeTerms.forbidden` — unreachable by construction |
| (d) σ38 variant `TATACT` | fixed string | `LatticeTerms.forbidden` |
| (e) antisense `ATTATA` | fixed string | `LatticeTerms.forbidden` |
| (a) −35 + 15–19 bp spacer + −10 | geometry | `evaluate()` → Tier-B repair |
| (c) AT-tract 3′ end exactly 5 bp upstream of a −10-like hexamer | geometry | `evaluate()` → Tier-B repair |

The two switches are genuinely independent rather than in conflict:
`solver/catalog.py:149` reads the **ClassVar** when collecting forbidden patterns for the
Tier-A automaton (and those same patterns reach the independent validator), while
`solver/catalog.py:302` (`repair_specs`) and `_enforcement_of` read **`enforcement_for`**
when routing breaches to Tier-B. `d4_internal_polya:93,156-161` already ships the same
split in the other direction — a SOFT ClassVar with a HARD_REPAIR `enforcement_for`.

`repair = FIXED_POINT`, not inherited and not copied: mutating one element of a promoter
creates −10-like hexamers nearby, and a new hexamer 5 bp downstream of an AT-tract that
was previously harmless is a **new** promoter a single pass would never re-scan for. That
is exactly the condition CLAUDE.md §3.6 names.

`gate()` returns True for **every** modality, including `IVT_MRNA`.

**Rejected:**

- *ClassVar `HARD_REPAIR` with a non-zero `steering_weight`.* Would let the DP be nudged
  away from the (a)/(c) geometry, but `solver/catalog.py:149` collects forbidden patterns
  from HARD_LATTICE rules **only**, so (b)/(d)/(e) would silently never reach the automaton
  *or* the oracle. Losing an absolute guarantee plus an independent proof, on the half that
  carries the disaster anchor and the novel finding, to buy a search hint. The cost of the
  choice actually made is real and is stated in the rule's docstring:
  `test_only_lattice_rules_skip_steering` forces `steering_weight = 0.0`, so repair does
  all of the geometry work with no DP steering.
- *ClassVar `HARD_LATTICE` and `enforcement_for` returning it too* (the `d1` shape). The
  geometry breaches would then match neither arm of `findings()`'s router
  (`solver/catalog.py:184-189`), so every one would be discovered, reported into the
  scorecard, and silently dropped by the only thing that could fix it.
- *`HARD_CHECK` for the geometry.* Wrong by the enum's own definition — HARD_CHECK is
  "real, but not fixable by codon choice", and a cryptic promoter built by synonymous
  substitution is precisely fixable by codon choice.
- *Splitting D5 into two registered rules.* `_index_by_ref` (`score/presets.py:120`) raises
  when two rules claim one `brief_ref`, and both halves are `2.D5`.
- *Reimplementing the Salis Promoter Calculator* (346 parameters, validated on 17,396
  promoters, `brief.md:111`). It is a model, not a threshold; an in-PR reimplementation
  would be a second, unvalidated model wearing a validated one's citation.
  `_scored_path_unavailable` reports its absence at magnitude 0.0 instead — `d3_splicing`'s
  pattern for its missing MaxEntScan model. Issue filed.
- *Gating off `IVT_MRNA`,* which `d3`, `d4` and `d6` all do. Copying them would have been
  wrong: those three read the **transcript**, and an in-vitro transcript is neither
  replicated nor spliced. D5 reads the **template plasmid during *E. coli* propagation**,
  which an IVT mRNA construct still passes through. `brief.md:110` — "applies to ALL
  constructs regardless of final host"; `brief.md:209` marks the cryptic bacterial
  promoter "On (propagation)" in every column.
- *Working around the reverse-complement closure.* `ATTATA` **is** the revcomp of
  `TATAAT`, and both the solver and `verify.find_motifs:187` close the pattern set, so
  listing `ATTATA` necessarily forbids a bare sense-strand `TATAAT` as well. Kept rather
  than engineered around: a sense `TATAAT` is the −10 box that (a), (b) and (c) are each
  built around, and the cost is small — a 6-mer occurs ~0.73× per 1.5 kb, so the pair costs
  about one site per 1.5 kb. The 1-mismatch variants that (a) and (c) turn on stay legal,
  so both geometry halves still do real work. Pinned by
  `test_a_sense_strand_tataat_is_also_caught_by_the_closure`.
- *Declaring `e2_gc_band` in `conflicts_with`.* Suggested by the wave plan, checked, and
  dropped: e2's band is (0.28, 0.80) over the whole construct — far too wide for a 6-mer
  substitution to bind against. `d8_cpg_depletion` and `f5_at_window` **are** declared:
  D5 removes AT-rich hexamers and so raises local GC against f5's 60% ceiling, while d8's
  CpG fixes lower GC and walk back toward what D5 forbids.

**Evidence:**

- Sub-rule text and every threshold verbatim: `docs/research/brief.md:110-112`. No `~~`,
  "corrected" or "superseded" marker anywhere in the block, unlike E4 (`brief.md:141`).
- The two independent switches: `solver/catalog.py:139-154` (ClassVar → automaton, and the
  docstring's own reason for the restriction) against `solver/catalog.py:184-189, 296-305`
  (`enforcement_for` → repair routing).
- The `RulePolicy` window is `getattr(spec, "window", 50)` and must be an `int`
  (`solver/catalog.py:286-290`), so `self.window` is 6 + 19 + 6 = 31, the widest
  architecture, letting the repair window reach either element.
- Variant enumeration is 19 strings at 1 mismatch and 154 at 2, verified directly; the hot
  path is 38 C-level `str.find` scans per candidate rather than a per-position character
  loop, because `evaluate` runs once per repair candidate.
- Boundary behaviour verified before the test was written: spacers 15/17/19 fire and 14/20
  do not; an AT-tract at exactly 5 bp fires and at 4 or 6 does not.
- The calibration anchor is **not** a test vector, and saying so matters: the dengue-2 −35
  `TCAACG` is 3 mismatches from `TTGACA`, so part (a) does not match it at the brief's
  1-mismatch budget (nor at 2). Its −10 half does — `TTTAAT` is 1 mismatch from `TATAAT`.
  The anchor is evidence that cryptic promoters cause disasters, not a claim about this
  rule's recall, and no test asserts otherwise.

**Where:** branch `claude/next-build-plan-cjmqj6`; wave-2 slice W3
(`docs/buildout/wave2/w3-d5-cryptic-transcription.md`).

---

## 2026-09-06 — the brief's "5 bp upstream" is a position offset, and D5 read it as a gap

**Decided:** `AT_TRACT_GAP = 4`, meaning **four bases between** the AT-tract and the −10-like
hexamer — not five. `brief.md:111` says the tract's "3′ end sits exactly 5 bp upstream of a
−10-like hexamer", which reads naturally as five intervening bases and is how D5 first
shipped. It is a **position** offset: Warman & Grainger place the tract at promoter
positions −23..−17 and the −10 at −12..−7, and −17 to −12 is five positions but four bases
in between. Their Table 1 primers settle it without arithmetic — `tatttat` + `TGAC` +
`tataat`.

The consequence of the original reading was not a near-miss. Run against both constructs
the 103/103 result was measured on, the rule returned **zero** AT-tract breaches, and fired
instead on a one-base-shifted class nobody has tested. The rule cited a result it could not
reproduce, and the paired test asserted the wrong offset was correct in both directions.

**Rejected:**

- *Keeping 5 because the brief says 5.* The brief is a summary; the rule-auditor's whole
  question is whether the cited source supports the number, and here it does not. Following
  the summary over the source would have shipped a rule that misses its own calibration set.
- *Widening to a window of 4-5, or 3-6.* "Appropriately positioned" is the entire finding —
  the paper's own control is a randomized tract at a **fixed** offset. A window would report
  the AT-richness of ordinary sequence, which is the noise this sub-rule exists to avoid.
- *Reporting both offsets at different magnitudes.* Same objection, plus it would invent a
  strength ordering no source measured.

**Evidence:** reproduced before the constant was touched — `TATTTAT`+`TGAC`+`TATAAT` and
`AATTT`+`TGAC`+`TATAAT` both scored 0 `at_tract` breaches at gap 5 and 1 at gap 4.
`test_the_papers_own_constructs_are_caught` now pins the primers themselves, so the citation
and the code cannot drift apart again.

**Also settled in the same audit:**

- The **dengue-2 anchor is not a test vector**, and nothing may imply it is. Its −35
  `TCAACG` is 3 mismatches from `TTGACA` (above even the schema's maximum budget of 2), and
  its real spacer is 13 bp — `brief.md:111`'s "17-bp spacer" does not survive the paper's own
  coordinates (nt 53 → nt 72). The paper also does not contain the word "uncloneable", and it
  showed the −35 **not** essential, so `sign="qualifies"`. Two tests pin the non-detection
  with the reason, because an unpinned negative reads as recall the rule does not have.
- **`Evidence.CONTESTED`, not `EVIDENCE_BACKED`.** `brief.md:110` marks D5 "H/S" and attaches
  no A/B/C letter (legend at `brief.md:46`), so the badge is the rule's own call. The part
  enforced most aggressively — `ATTATA`, unreachable by construction — rests on a
  computational study with no transcription assay.
- **`brief.md:209` cannot support the IVT gate** and was dropped from `gate()`'s docstring:
  it is the by-expression-host table and has no IVT mRNA column. `brief.md:110`'s "applies to
  ALL constructs regardless of final host" carries the decision alone.

**Where:** branch `claude/next-build-plan-cjmqj6`, commit `18a3163`; PR #153.
