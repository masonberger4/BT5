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
