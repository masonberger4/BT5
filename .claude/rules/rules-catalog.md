---
paths:
  - "packages/engine/src/bt5/rules/**"
  - "packages/engine/tests/rules/**"
---

# Writing a rule

Statements of the non-negotiable rules live in root `CLAUDE.md` §3. This file is the
rationale and the mechanics.

## Resolving `brief_ref` — it does not resolve by literal search

Every `Spec` carries `brief_ref`, "section id in `docs/research/brief.md`". The values
look like `2.E4`, `2.B1`, `2.D1`. **None of those strings appears in `brief.md`** —
`grep -F "2.E4"` returns zero hits for all 25. The reference is *section-qualified* and
resolves in two steps:

1. **Section:** `grep -n '^### 2\.E' docs/research/brief.md`.
   The four sections are `### 2.B` (translation initiation / 5′), `### 2.D` (motif
   avoidance), `### 2.E` (manufacturability), `### 2.F` (plasmid propagation).
2. **Row**, beneath that heading. Ids appear in exactly two shapes:
   - a **table row** — `| E4 | **GC variation.** ... |`
   - a **bold run-in** — `**D1 Restriction / Type IIS (H, scan BOTH strands ...):**`

So `2.E4` = section `### 2.E`, row `E4`. Delegate the lookup to `docs-miner`; `brief.md`
is 63 KB and reading it inline costs more than the answer.

### Superseded rows

`brief.md:141` — E4 — has its original thresholds struck through with `~~...~~` and
marked **"corrected 2026-08-28, these thresholds are below the chance floor."** A rule
encoding the struck-through number would pass all 11 contract assertions and still be
wrong. Always scan a resolved row for `~~`, `corrected`, `superseded` or `provisional`.

## What CI already checks — do not re-derive it

`tests/data_integrity/test_rule_contract.py` parametrises these over every registered
spec: citations exist, are `https`, are labelled · `weight_provenance` non-empty for
SOFT · FOLKLORE defaults off · VENDOR_ASSERTED defaults on · `last_verified` parses ·
BAND declares a non-inverted band · `default_weight == 0.0` for every hard rule ·
`steering_weight == 0.0` for HARD_LATTICE · HARD_LATTICE declares forbidden motifs ·
`brief_ref` non-empty · `param_schema` is a JSON Schema object. Plus unique, sorted ids
and no self-conflict.

Its docstring names its own limit: *"catches MISSING provenance; it cannot catch WRONG
provenance."* That gap is what `rule-auditor` and `/verify-provenance` are for.

## Reverse strand: list motifs, do not scan (CLAUDE.md §3.4)

Put forward motifs in `LatticeTerms.forbidden`; the solver closes the set under reverse
complement. Hand-rolling the scan is how a junction-spanning or origin-spanning hit gets
missed.

But **directional scored models are not revcomp-symmetric** — MaxEntScan, Salis TIR, the
promoter calculator, the polyA downstream element. Those must read
`slot.strand_of_interest` (`core/spec.py:274`). Hard-coding strand 1 makes a
reverse-oriented lentiviral cassette's polyA and splice analysis exactly backwards, and
nothing fails loudly when it does.

## Repair policy: `SINGLE_PASS` is a downgrade

`RepairPolicy.FIXED_POINT` is the default of `solver/repair.py:417`. 22 of the 25
current catalog rules explicitly declare `SINGLE_PASS`, so each of those is an opt-out;
`d3_splicing`, `b9_out_of_frame_atg` and `f5_at_window` keep `FIXED_POINT`.

Point-mutating one cryptic splice donor activates cryptic donors nearby. A single pass
ships a construct whose donors were removed *into* new donors — and the validator passes
it, because the specific 9-mer it was told to avoid is gone. If your repair can create
new instances of what it removes, `FIXED_POINT` is mandatory (§3.6). Otherwise say in the
docstring why single-pass is sufficient.

## Rules take a `Construct`, never a string (§3.3)

That is what makes junction-spanning, origin-spanning and reverse-strand hits impossible
to miss rather than something you remembered to handle.

## Mechanics

- **No registration step.** `@register` autodiscovers. `core/registry.py`: *"Adding a
  rule edits ZERO shared files, not even an `__init__.py`."* There is deliberately no
  committed catalog list — it collided on every PR at the highest-volume lane.
- **One file per rule**, named `<brief_id>_<slug>.py`, and `id` equals the filename stem.
- **Reference rules:** `e2_gc_band.py` for a scored rule, `d1_restriction_sites.py` for a
  pattern rule. **`d1` has no paired test** — it is the only catalog rule without one, so
  copy `e1_homopolymers.py` + `test_e1_homopolymers.py` for the test shape.
- **`engine_calibration`** must name the engine a kcal/mol threshold was measured on.
  `registry.register` raises `CalibrationMismatchError` otherwise: applying a
  ViennaRNA-calibrated number to another engine's output fails silently.
