## 2026-09-02 — D8 CpG: two of the three named TLR9 hexamers are not `RRCGYY`

**Decided:** `d8_cpg_depletion` ships `Enforcement.SOFT`, `Evidence.CONTESTED`, with the
three metrics `brief.md:128` demands as three independent toggles, and it scans the named
TLR9 hexamers **unioned with** the `RRCGYY` consensus rather than as a filter applied to
it.

### The finding: "specifically" does not mean "for example"

`brief.md:129` reads: *"total CpG count + stimulatory hexamers `RRCGYY`, specifically
`GTCGTT`, `TTCGTT` (human), `GACGTT` (mouse)."* The natural reading — that the three are
instances of the consensus, so scanning `RRCGYY` covers them — **is wrong**, and wrong in
the silent direction:

```
R = A or G

GTCGTT ->  G  T  C  G  T  T     position 2 is T.   NOT RRCGYY
TTCGTT ->  T  T  C  G  T  T     position 1 is T.   NOT RRCGYY
GACGTT ->  G  A  C  G  T  T                        matches
```

Both **human** motifs fail the consensus; only the mouse one matches. A rule that scans
`RRCGYY` and then labels hits by membership in the named set — which is what this rule
did in its first draft, and which passed a smoke test on CpG-dense sequence because the
ZAP and island metrics fired — reports **zero** human TLR9 hexamers on a construct full
of them. It was caught only because the first three TLR9 tests asserted a specific
magnitude and got an empty list.

This is not an error in the brief: the classic human CpG-B motifs really are `GTCGTT` and
`TTCGTT`, and they really are not the `RRCGYY` consensus. The brief is listing two things,
not one thing and its examples.

**Decided:** `TLR9_MOTIF` is the alternation of the three named hexamers *and* the
consensus. `TLR9_CONSENSUS_ONLY` is kept as a separate compiled pattern purely so the
test can assert the non-membership that makes the union necessary
(`TestTlr9::test_the_named_human_hexamers_are_not_rrcgyy`) — otherwise a later
simplification back to "just scan the consensus" would look like a tidy-up.

**Rejected:** *scanning only the three named hexamers.* The consensus is a real, separate
part of the row and catches motifs the three do not.

### Species attribution is not guessed

`brief.md:129` attributes hexamers to human and mouse only. `HOST_SPECIES` therefore maps
`HUMAN`/`HEK293` → human and `MOUSE` → mouse, and **nothing else**. A hexamer escalates
only when its documented species matches the slot's host; in CHO, or for an unattributed
`RRCGYY` hit, it is reported at the general magnitude.

**Rejected:** *mapping CHO to human as "the nearest mammal".* CHO is hamster, the brief
attributes nothing to it, and an invented attribution is the defect `/verify-provenance`
looks for.

### SOFT, and CONTESTED, on the evidence's own terms

`brief.md:128`'s header carries **no H/S marker** — unlike D3's `(H/S, ...)`, D4's
`(H for lentiviral sense strand, S elsewhere)` or D2's `(H, check-only)`. Absent a hard
grade the rule stays a weighted objective and never refuses a construct.

The badge is `CONTESTED` because the ZAP arm's own primary source cuts both ways:
`brief.md:130` quotes *"the magnitude of ZAP-mediated inhibition was not correlated with
the number of CpGs introduced."* That is carried as a third `Citation` with
`sign="refutes"` — the field exists for exactly this, and `Citation`'s docstring names
"ZAP CpG" as one of the three cases it was added for. It is also why the ZAP metric
reports **the worst 200-nt window and never a global count**, as `brief.md:130` instructs:
the count is the quantity the source says does not predict the effect.

**Rejected:** *`EVIDENCE_BACKED`.* It would badge a contested quantity as settled in the
one place the badge is supposed to be load-bearing.

### The weight prices a certain cost against a contested benefit

`default_weight = 0.35`, below `d4_internal_polya`'s 0.7. `brief.md:133` records that full
CpG depletion forces AGA/AGG for Arg and can drop GC below vendor floors — a *certain*
manufacturability loss bought with a *contested* immunological gain. `e2_gc_band` already
names `d8_cpg_depletion` in its own `conflicts_with`; this rule now names it back.

### Metric (a) is a count as well as a motif scan

`brief.md:129` says "total CpG count **+** stimulatory hexamers". The count ships as one
breach at magnitude 0.0 with `fixable_by_codon_choice=False` — a measurement rather than a
finding, so it stays out of `passes`, out of the weighted sum's meaning and off the
solver's target list. Same channel `d3_splicing` uses to report that a scan did not run.

### Performance, and a discrepancy worth an issue

`RuleSet.findings` (`solver/catalog.py:176-190`) calls `spec.evaluate` for **every** spec
in `self.specs` and filters by enforcement afterwards. Its `breach_finder` docstring
(`catalog.py:200-206`) says the callable is *"SCOPED TO THE HARD_REPAIR RULES, which is a
performance decision"* — that describes the intent, not the loop: `build_rule_set`
(`catalog.py:363`) puts every gated-in spec in `specs`, so a **SOFT rule is still
evaluated once per candidate, up to 256 per repair iteration**.

That is why `_zap` sweeps its 200-nt window with two pointers over the CpG position list
instead of re-slicing per window. It is also worth an **M1 issue**: either the docstring
overstates the scoping or `breach_finder` should build from `repair_specs()`. Not this
lane's edit.

**Evidence:** `brief.md:128-133`; `core/spec.py` `Citation` docstring;
`solver/catalog.py:176-190, 200-206, 363`; `e2_gc_band.py:151`. D8 evaluates a 5 kb
circular construct in 2.7 ms. `pytest packages/engine/tests/rules/test_d8_cpg_depletion.py`
33 passed.

**Where:** branch `claude/s4-rules-liabilities`, session S4 of the six-way buildout.
