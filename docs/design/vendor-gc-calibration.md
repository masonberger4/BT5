# Vendor GC calibration — measured, not read off a webpage

Status: **complete — both vendors, all 18 probes**. Measured 2026-08-28.

`docs/PLAN.md` schedules a **vendor-complexity oracle** — "paste designs into the
Twist/IDT complexity checkers and record accept/complex/reject", as *manual
periodic calibration, not CI*. This is that loop, run for the first time.

It exists because the published numbers have already been caught being wrong
twice. `docs/design/repeats.md` §2 records the first: Twist's published trigger
tolerates a 200 bp direct repeat while constructs with 20–81 bp repeats were
rejected outright. This document records the second, for GC — and the gap is
just as wide.

---

## Method

18 sequences, each **500 bp**, homopolymers capped at 3, and every one certified
against BT5's own homopolymer / repeat / STR / k-mer / inverted-repeat rules
before submission, so a verdict can only be about GC. Two ladders:

- **LOC_** — one interior **50 bp** window swept 5→95% GC on a **byte-identical**
  50% background. Isolates any *windowed* rule.
- **GLB_** — uniform whole-sequence GC 20→80%, exact to the percent (25.0, 65.0,
  75.0 all hit exactly). Isolates the *global* rule.

Generator and panel: scratch, regenerable. IDT confirmed the construction — it
independently reported "overall GC content of the sequence is 20% / 25% / 70% /
75% / 80%", matching the designed values exactly.

---

## IDT — gBlocks plate-entry complexity checker

### The tiers, from the UI itself

| score | verdict |
|---|---|
| < 12 | clean (green) |
| **12 – 24** | **"Accepted – Moderate Complexity"** — *"we will attempt this order"* |
| **≥ 24** | **"Denied – High Complexity"** — *"prevent manufacturing … as a fragment"* |

### Results

| probe | GC | verdict | score |
|---|---|---|---|
| LOC_w50_gc05 … gc95 (all 8) | one 50 bp window at 4–96% | **all green** | <12 |
| GLB_gc20 | 20% | **Denied** | 64.6 |
| GLB_gc25 | 25% | **Denied** | 27.7 |
| GLB_gc30 / 40 / 50 / 60 / 65 | 30–65% | green | <12 |
| GLB_gc70 | 70% | Accepted, Moderate | 14.2 |
| GLB_gc75 | 75% | Accepted, Moderate | 21.2 |
| GLB_gc80 | 80% | **Denied** | 28.2 |

### The global-GC penalty is exactly linear, and now solved

Both high-GC probes produced **one finding only** — "overall GC content" — making
this a clean single-variable measurement:

```
score = 1.40 × GC%  −  83.8          (for GC above ~60%)
```

| GC% | predicted | observed |
|---|---|---|
| 60 | 0.2 | green ✓ |
| 65 | 7.2 | green ✓ |
| 70 | 14.2 | **14.2** |
| 75 | 21.2 | **21.2** |
| 80 | 28.2 | **28.2** |

Three exact hits and two consistent greens. The 75% point was **predicted before
it was measured** and came back to the decimal.

Which converts the tiers into GC thresholds:

| | global GC |
|---|---|
| penalty onset | **59.9%** |
| Moderate (score 12) | **68.4%** |
| **Denied (score 24)** | **77.0%** |

### The rules IDT states in its own remediation text

The checker does not merely score — it names its targets:

| rule | IDT's wording | so the spec is |
|---|---|---|
| global GC, high | *"Redesign to reduce the GC content below 62%"* | target ≤ 62% |
| global GC, low | *"Redesign to increase the GC content above 32%"* | target ≥ 32% |
| **windowed GC** | *"a window of **100 bases** starting at base 117 with a GC content of 15%. Redesign this region to have a GC content **greater than 30%**"* | **100 bp window, floor 30%** |
| **terminal GC** | *"GC content of 18.3% in the **terminal 60 bases of the 3′ end**. Redesign the end to have a GC content **greater than 24%**"* | **60 bp 3′ terminus, floor 24%** |
| end proximity | *"a region of low GC content on or near the 5′ / 3′ end"* | scored 10 each, threshold unstated |

### Three findings that change what BT5 should encode

**1. IDT's GC window is 100 bp, not 50.** E2 uses 50 bp. This is why all eight
LOC probes passed: a 50 bp window at 4% GC, diluted by 50 bp of 50% background,
is ~27% across 100 bp — barely under the 30% floor, scoring a few points and
staying green. A 50 bp rule is the wrong geometry for IDT and will flag
sequences IDT does not care about.

**2. The rule is asymmetric — low GC is policed locally, high GC only globally.**

- `GLB_gc20` (64.6): a 100 bp window finding, the terminal-60 finding, *both*
  end findings, and global. Global contributed only **9 of 64.6**.
- `GLB_gc80` (28.2): **one** finding, global GC. No window finding. No terminal
  finding, despite being 80% GC throughout.

This independently confirms SCP4ssd's conclusion that *"local fragments with low
GC content might have a more important impact than fragments with high GC
content"* — an ML model trained on synthesis outcomes and a production order
checker agreeing from opposite directions. It also means a **two-sided windowed
band is the wrong shape**: IDT enforces a windowed *floor* and no windowed
ceiling at all.

**3. Twist and IDT do not share a geometry.** Twist publishes 10–90% over
**50 bp**; IDT enforces a 30% floor over **100 bp** plus a 60 bp terminal rule.
A single windowed GC rule cannot represent both, which is the strongest evidence
yet for issue #43's per-vendor evaluation (V3) — the same sequence has to be
scored twice, under two window sizes.

---

## Twist — "Analyze your gene sequences for manufacturability"

| probe | GC | verdict |
|---|---|---|
| LOC_w50_gc05 … gc95 (all 8) | one 50 bp window at 4–96% | **all Standard** |
| GLB_gc20 | 20% | **Not Accepted** |
| GLB_gc25 | 25% | **Not Accepted** |
| GLB_gc30 … gc80 (all 8) | 30–80% | **all Standard** |

Twist accepts **80% GC as Standard** — not "complex", Standard — where IDT
denies the same sequence. Its upper bound was never reached by this panel.

### Stated rules

- low GC: *"Increasing GC to **> 25%** will be optimal for success"*
- repeat density: *"More than **15%** of your sequence is composed of small
  repeats (**9bp or longer**)… break up repeats, perhaps by varying your codon
  usage"* — independent corroboration of E6's `KMER_BP = 9`
- positional: "Problematic area" spans, *"repeats or extreme high/low GC"*

### Two things the data refused to confirm

**The repeat metric is not what I reconstructed.** Reading "15% of your sequence
composed of repeats" as *positions covered by any 9-mer occurring twice* gives
25.2% for `GLB_gc20` (flagged) but **19.0% for `GLB_gc80`, which passed as
Standard**. So that is not their metric. The data bounds it only loosely — the
trigger sits between 19.0% and 25.2% *in those units*, or the definition differs
entirely. Recorded as unresolved rather than guessed.

**The "Problematic area" detector is not a windowed GC threshold.** For
`GLB_gc25` the flagged spans were lower-GC than the rest (16.7% / 24.3% against
26.5%), which fits. For `GLB_gc20` it **inverts**: the *unflagged* remainder was
lower-GC (17.0%) and far more repetitive (45%) than the flagged spans. Sliding
50 bp windows do not separate them either (flagged 18–30%, unflagged 16–30%).
Two probes cannot resolve it, and it is not reducible to a windowed bound.

---

## The combined result

| global GC | Twist | IDT |
|---|---|---|
| 20% | Not Accepted | Denied (64.6) |
| 25% | Not Accepted | Denied (27.7) |
| 30 – 65% | Standard | clean |
| 70% | Standard | Moderate (14.2) |
| 75% | Standard | Moderate (21.2) |
| 80% | **Standard** | **Denied** (28.2) |
| **one 50 bp window, 4–96% GC** | **Standard ×8** | **green ×8** |

Three conclusions, each resting on all 18 probes at both vendors:

**1. A single extreme 50 bp window is irrelevant to both vendors.** Sixteen of
sixteen passed, spanning 4% to 96% local GC. Neither vendor's verdict moved.

**2. Global GC is what gates, and both reject the same low end.** ≤25% is
refused by both; ≥30% is accepted by both.

**3. The vendors diverge at the high end, sharply.** Twist ships 80% as
Standard; IDT denies above ~77%. A single GC ceiling cannot serve both, which is
[#43](https://github.com/masonberger4/BT5/issues/43) V3 (per-vendor evaluation)
arriving as a measurement rather than an argument.

---

## What this says about E2 as shipped — X7, settled

`e2_gc_band` ships **40–60% GC over a 50 bp window at `HARD_REPAIR`**, so the
validator refuses to emit outside it. Against 18 probes at two vendors that is
wrong in every dimension at once:

| | E2 as shipped | measured |
|---|---|---|
| **quantity** | windowed GC | **global GC** — the windowed probes all passed |
| **window** | 50 bp | **no vendor gates on a 50 bp window at all** |
| **shape** | two-sided band | asymmetric, and vendor-specific at the top |
| **bound** | 40–60% | both accept **30–65%** clean; Twist to **80%** |
| **enforcement** | **HARD**, refuses to emit | Twist calls 80% *Standard* |

E2 refuses to emit sequences **both vendors would manufacture without comment**.
Its 60% ceiling is 20 points below Twist's demonstrated tolerance, applied to a
geometry neither vendor uses, at the one enforcement level that blocks output.

### The proposal

Split the number into the two things it is doing:

- **Hard bound → global GC, per vendor.** Reject below **28%** (both vendors
  refuse 25%, accept 30%). Ceiling is vendor-specific: **77%** for IDT (from the
  fitted denial threshold), none demonstrated for Twist. This belongs in the
  vendor profile of
  [#43](https://github.com/masonberger4/BT5/issues/43) V1, not in a module
  constant.
- **40–60% → steering target only.** It stays a design preference and keeps
  pulling on codon choice through `steering_weight = 0.5`, which is what
  actually selects codons. It stops being a gate.
- **Keep a windowed floor, at IDT's geometry.** IDT states one: **100 bp window,
  GC > 30%**. That is the only windowed rule any vendor was observed to apply,
  and it is a floor with no ceiling — so the two-sided windowed band goes.

That change makes a `(GGGGS)₆` construct emittable, which is the failure X1
surfaced: its 62% window floor is inside every measured vendor tolerance and
outside only E2's invented one.

## Caveats

- **Length is fixed at 500 bp.** IDT's score may scale with length; every number
  here is calibrated at 500 bp and should not be extrapolated to a 3 kb order
  without re-running the ladder.
- **The low-GC global points are partly confounded.** `GLB_gc20`'s denial was
  driven mainly by a 15% *100 bp window*, not by its 20% global GC — the panel
  was flattened per 50 bp block, so sliding 100 bp windows still vary. The high
  side is clean (a single global finding both times), which is why the fitted
  function is quoted for GC above 60% only. A sharp low-side global number needs
  a panel flattened at 100 bp.
- **The penalty onset at 59.9% is an extrapolation** from three points all above
  62%. It is consistent with the 60% and 65% greens but is not directly measured.
- **This is one product line** (gBlocks plate entry). eBlocks and Twist Gene
  Fragments may differ; the panel is cheap to re-run against each.
