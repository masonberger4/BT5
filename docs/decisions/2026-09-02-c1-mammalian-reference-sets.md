## 2026-09-02 — C1 gains four hosts; C3 gains none, and that is the same fact twice

**Decided:** `CAI_REFERENCE_SET` in `c1_cai.py` now maps HUMAN, HEK293, MOUSE and CHO
onto the reference sets S6 shipped in #90, alongside the two E. coli entries. C1 computes
for **five of nine** hosts instead of two. `MINMAX_REFERENCE_SET` in `c3_min_max.py`
stays **empty**. The same three tables cause both outcomes.

This is the follow-up S6 explicitly declined and handed over:
`docs/decisions/2026-09-02-s6-host-data-and-real-backbone.md` lists under *Rejected* —
"**Wiring `CAI_REFERENCE_SET` myself.** That map lives in `c1_cai.py` (S3's engine file);
the data lane ships data only."

**HEK293 → the human set, and it is S6's mapping rather than an invention here.** That
record's own table reads "`human_highly_expressed_refseq_w` … serves hosts HUMAN,
HEK293", with the reason: "HEK293 → human is the load-bearing one: both shipped mammalian
presets (`lentiviral_hek293`, `aav_hek293`) key on HEK293, so a human set is what lets C1
score instead of report `unavailable` for the presets the walking skeleton actually runs."
It is also the approximation `c1_cai.py` already makes one taxon up — BL21 shares K-12's
entry — and HEK293 is a *Homo sapiens* cell line, so this is the same move with less
distance in it: codon usage is a property of the organism's translational machinery, not
of the cell line. Stated rather than hidden, as before: the set used travels in every
breach's `detail["reference_set"]`.

**The band did not transfer, and the first draft of this change shipped it anyway.**
`rule-auditor` caught it before the PR left draft. brief.md:77 offers its numbers as an
example — "Target a band (e.g. 0.70–0.90, or ±0.1 of host median)" — and they were
calibrated on E. coli. Measured on the shipped tables, with a composition-neutral random
synonymous encoding as the chance baseline:

| table | chance CAI | the 0.70 floor sits |
|---|---|---|
| `sharp_li_1987_ecoli_w` | 0.238 | **0.46 above chance** |
| `human_highly_expressed_refseq_w` | 0.656 | 0.04 above chance |
| `mouse_highly_expressed_refseq_w` | 0.633 | 0.07 above chance |
| `cho_highly_expressed_refseq_w` | 0.660 | 0.04 above chance |

Mammalian bias is weak (brief.md:206, "isochore GC, not selection"), so the tables are
near-flat and the same floor loses nearly all its discriminating power. `CAI_BAND` now
holds a band per host, and the two halves are treated differently because they behave
differently:

- **The floor is 0.0 for the weak-bias hosts and can never bind.** The obvious repair —
  rescale the floor to preserve E. coli's headroom — was tried and rejected: it puts
  human's floor at **0.864**, *above* where a native human CDS sits, so C1 would flag
  native sequence as "rare codons across the ORF" and hand the optimizer pressure to
  raise its CAI. That is the opposite of what the evidence says. brief.md:206 grades the
  CAI weight "very low" for CHO/HEK with default mode "Native or harmonize"; brief.md:215
  marks per-host evidence "low for human, mouse, CHO, Sf9, Tni"; brief.md:13's Expi293F
  benchmark found optimization did not increase yields; and a 2026 Pichia study found CAI
  *negatively* correlated with titer (−0.81 for trastuzumab). Nothing here claims a
  low-CAI mammalian CDS is worse than a native one.
- **The ceiling does transfer and stays operative**, scaled to each host's own
  chance-to-1.0 headroom so it means the same thing: E. coli's 0.90 is 0.8687 of its
  headroom, and every other ceiling is that same fraction of its own (human/HEK293
  0.9548, mouse 0.9519, CHO 0.9553). Max-CAI collapse is a *mechanical* failure — it
  drives each amino acid onto one codon and manufactures perfect direct repeats — and
  that is true of any organism. **E. coli's pair is unchanged at exactly (0.70, 0.90).**

`TestBandCalibration` re-derives every constant from the shipped tables, so they cannot
drift from the data they came from, and pins that the mammalian floor is inert while the
mammalian ceiling still bites.

**What the earlier "measured, not assumed" claim was worth.** It said the four hosts
return 0.807 / 0.788 / 0.801 against E. coli's 0.763, "all inside the band". The
arithmetic was right and the conclusion was hollow: `IN_BAND` is a five-codon repeat
chosen to land in-band *for E. coli*, and landing inside a band whose floor is four
hundredths above chance is close to unavoidable. It measured that the tables are
distinct, which is a different proposition from the band being valid.

**C3 stays dark, and the reason is structural.** All four shipped tables — including the
three new ones — are relative-adaptiveness w-indices: `w = (count + 0.5) / family_max`,
so every family is renormalised to its own peak and the peak is discarded. %MinMax sums
per-family differences ACROSS families, and the `K_i` that would make that sum comparable
is exactly what the normalisation throws away. More w-tables cannot fix it; only a
frequency table can. `docs/decisions/2026-09-02-c3-minmax-needs-frequencies.md` predicted
this exact case and it held.

**Also fixed: a mis-citation that had been on `main` since before this branch.**
`c1_cai.py` carried "Nine benchmarked commercial optimizers were a coin flip against
native sequence, and all computable design features together explain 5-31% …" under a
single URL, `nature.com/articles/nbt.4238` — Cambray 2018. The first half is **Ranaghan
2021** (`PMC7893858`), which a 2018 paper cannot have reported, and C1 cited Ranaghan
nowhere else. Split into two citations, each under the paper that made its claim. Found
by the `rule-auditor` sweep on the S3 rules PR and deliberately left out of that PR as
out of scope; this is its own change, in the same file, under one review.

**Rejected:**
- *Mapping S_CEREVISIAE, P_PASTORIS or SF9 to anything.* S6 deferred them deliberately —
  no shipped preset consumes them and their RefSeq symbol/transcript coverage is
  materially messier. They stay absent from the map and keep reporting `unavailable`,
  which is the honest state, not an oversight. A test now pins that this is on purpose.
- *Computing %MinMax from the new w-tables because there are more of them now.* Three
  more of the wrong shape is still the wrong shape.
- *Rescaling the floor to preserve E. coli's headroom.* Measured at 0.864 for human —
  above native sequence — so it would manufacture the exact pressure the evidence
  says not to apply. Rejected on the number, not on principle.
- *brief.md:77's own second formulation, "±0.1 of host median".* The brief offers two
  band constructions and this one is already per-host, so it deserved measuring rather
  than being passed over. Two things sink it. First, the input is not on disk: a "host
  median" means the median CAI of that host's own native genes, and `data/codon_usage/`
  ships w-indices only, so the quantity the brief names cannot be computed here without
  inventing it. Second, the shape of the answer is wrong wherever that median lands in
  the plausible range. The nearest computable stand-in — the expected CAI of a random
  synonymous encoding — is 0.656 on the human table, and ±0.1 around it puts the
  **ceiling at 0.756**: 0.10 above the score a random encoding gets for free, and 0.20
  below the 0.9548 at which max-CAI collapse actually begins. A ceiling that close to
  chance cannot separate the mechanical failure it exists to catch from ordinary
  sequence — this file's own `IN_BAND` fixture scores 0.807 against the human table
  and would breach it. That is a claim about one computable point, not about every
  median the construction might use: a native median well above ~0.855 would put the
  ±0.1 ceiling back near 0.955 and the two constructions would agree. The objection
  that does not depend on where the median lands is the first one — the median of a
  host's native genes is not in `data/codon_usage/`, so this band cannot be built
  here at all without inventing its input. Rejected on that, with the computed point
  as evidence that guessing the input would not be harmless.
- *Leaving the E. coli band on the mammalian hosts and documenting the caveat.*
  A disclosed wrong threshold is still a live 0.2-weighted objective pushing on
  native sequence; a comment does not stop a solver.
- *Giving HEK293 its own table.* There is no HEK293-specific codon usage to have; the
  cell line's translational machinery is human.
- *Folding this into the S3 rules PR (#92).* It was already reviewed and attested at a
  fixed SHA; adding an unrelated change would have invalidated both.

**Note on the guard that replaced a weaker one.** The old
`test_only_the_e_coli_hosts_have_a_reference_set_in_this_build` hard-coded the map's
contents, so it could only ever fail by going stale. It is replaced by a check that every
mapped host's table actually exists on disk — which is the failure that matters (a map
entry without a file turns an honest `unavailable` into a `FileNotFoundError` from inside
a rule) and which keeps working as the map grows.

**Still open, and not fixed here.** brief.md:215 asks for a per-host evidence-strength
badge — high for E. coli, low for human/mouse/CHO/Sf9/Tni. `Spec.evidence` and
`default_weight` are single ClassVars with no per-host channel, so C1 ships one
`CONTESTED` badge and one 0.2 weight across hosts whose evidence the brief grades
differently. Expressing that needs a `core/` change (a new field or a per-host accessor),
which is MAJOR under CLAUDE.md §2a and belongs to an RFC, not to this branch. The band is
the half that was expressible here; the weight is not.

**Where:** branch `claude/s3-c1-mammalian-hosts`.

## Follow-up in the same PR: the override sentinel, and three surfaces the per-host band made stale

A delta review of the band change found one blocking defect and a set of surfaces that
still described a single band. All are fixed in this PR rather than deferred, because
each is a statement about what the rule refuses.

**Blocking — `cai_min`/`cai_max` used value-equality as their "unset" sentinel.**
`_lo_override = cai_min != BAND_LO` cannot tell "omitted" from "explicitly set to
0.70". On the floor that was a silent no-op; on the ceiling it **loosened** the limit
the caller asked for. `CodonAdaptationIndex(cai_max=0.90)` on a HEK293 job is a caller
deliberately pinning the tighter published anti-max-CAI ceiling, and it resolved to
0.9548 — so a HEK293 CDS at CAI 0.95 passed a rule explicitly configured to refuse it,
and the only way to obtain 0.90 was to perturb the value to 0.8999999. That is the
"silently permitted higher CAI" failure this rule exists to prevent, reached through
the rule's own parameter. Fixed with a real `None` sentinel, which is what
`e2_gc_band.py:189-190` already does per vendor and what C1's own comment claimed to
be copying. Two consequences handled with it: the range guards now run only over
supplied values, and `_band_for` validates `lo < hi` on the **composed** pair, because
once one side can be supplied alone the old accidental protection against an inverted
band is gone (`cai_min=0.93` against E. coli's 0.90 ceiling would otherwise have
reported a max-CAI sequence as "rare codons").

**No fallback band for an unmapped host.** `_band_for` returned `(BAND_LO, BAND_HI)`
for a host absent from `CAI_BAND`, which would reinstate this whole bug one host at a
time: the next host to gain a reference set would be scored against E. coli's band
with no error and no test failure. It now returns `None` and the caller reports the
objective unavailable, and a test pins `CAI_BAND.keys() == CAI_REFERENCE_SET.keys()`.
The one exception is a caller who supplied *both* bounds — then the band is theirs and
needs no host calibration.

**The `band` ClassVar is the loosest envelope, computed.** It still read
`(0.70, 0.90)` while four of the six scoring hosts use something else. It is now
`(min lo, max hi)` over `CAI_BAND` = `(0.0, 0.9553)`, computed rather than transcribed
so adding a host cannot make it a lie, and `weight_provenance` says it is an envelope
and not the gate — `e2_gc_band`'s convention, which `solver/catalog.py:272` states
outright ("Read off the INSTANCE, never `Spec.band`").

**`param_schema` advertised defaults it cannot honour.** `"default": 0.70` / `0.90`
are E. coli's. A form materializing them showed a HEK293 job the wrong band and, if it
posted its own displayed values back, silently became an override. Both `default` keys
are dropped and the descriptions now name `CAI_BAND` — again e2's answer to the
identical problem one axis over.

**The multi-slot tiebreak was degenerate under an inert floor.** Among in-band slots
the rule ranked by distance from the band's *midpoint*, which for `(0.0, 0.9548)` is
0.477 — below where any real mammalian CDS sits. The mammalian slot was therefore
always the "least interesting" one and a two-slot report handed itself to the E. coli
slot, discarding the CAI the job is actually about. Now decided by **`slot.role`**, as a declared policy rather than a derived number,
after three attempts to derive it from the CAIs each failed in the same way. Each
produced a *fixed preference for one host across a large region of realizable CAI*
while presenting as a neutral comparator:

| comparator | who wins, and where |
|---|---|
| distance to the band's **middle** | E. coli, always — the middle of (0.0, 0.9548) is 0.477, below where any real mammalian CDS sits |
| **raw** distance to nearest live edge | E. coli, for every human CAI below 0.8548 |
| distance / **band width** | the mammal, *unconditionally* |
| distance in **chance-to-1.0 headroom** | E. coli, for every human CAI below 0.9096 |

The width version inverted the preference rather than removing it: an in-band E. coli
slot scores at most 0.5 while a mammalian slot scores `(hi-cai)/hi` < 0.5 for any CAI
above 0.477, which chance alone (0.656) already clears. It also rested on a quantity that
is not a measurement — HEK293's band is 0.9548 wide only because the inert floor is
*written* as `0.0`, and writing it as chance (0.656) flips the result.

Headroom fixed that encoding dependence and was still a fixed preference, in a way that
was worse than what it replaced: E. coli carries **two** live edges in a narrow band, so
its distance is bounded above by 0.1313, while a weak-bias host carries **one** (its floor
is inert) and its distance is unbounded. Headroom multiplies the mammalian distance by
2.904 and E. coli's by 1.313 — a constant 2.212× tilt — so the dominance region *grew*,
from human CAI 0.8548 up to 0.9096. The commit that shipped it rejected raw distance in
its own comment as "the same pathology, narrower" and then shipped the same pathology,
wider.

**The common failure is structural, and no rescaling addresses it.** "Nearest a binding
edge" is not a symmetric predicate across slots whose hosts have different *numbers* of
live edges: a min over two edges is systematically smaller than a min over one drawn from
a wider interval. Every rescaling moved the dominance region instead of removing it,
which is why the fourth attempt stops rescaling.

So the policy is stated. CAI is a statement about **translation** — which is already why
`gate` keys on `role` and not on `host`, and why propagation slots are excluded outright —
and the producer slot is where the protein is actually made. "The CAI this job is about"
is the producer's, by the same argument that excludes propagation. Within one role the
higher CAI wins, which is the side the operative ceiling is on and is deterministic.
Breaches are unaffected and still outrank every in-band slot, ordered by deviation: that
ordering is about safety, not relevance, and a breach must never be hidden behind the
producer's comfortable number.

The tests changed shape with it. The previous set asserted five specific winners, of which
only three could distinguish the new comparator from the old one and *none* could
distinguish headroom from raw distance — reverting to raw would have left the suite green.
The set now includes the two cases where the rejected comparators disagreed (human 0.88
and 0.70 against E. coli 0.80), a role-swap on an identical construct proving the answer
follows the role and not the host, and a case where a breaching target must still outrank
an in-band producer.

**Recorded, not fixed:** `Breach.magnitude` is a raw CAI deviation and is not
comparable across hosts — human's chance-to-1.0 headroom is 0.344 against E. coli's
0.762, so a full max-CAI collapse yields |1 − ceiling| = 0.045 on HEK293 against 0.100
on E. coli, a factor of **2.2**, which is exactly that headroom ratio.
Normalizing by headroom would be more honest within C1 but would misrank C1 against
every other rule in `score/conflicts.py`, which compares unnormalized magnitudes
repo-wide. The non-comparability is now stated where `rank` is defined rather than
left to be inferred.

**Left to its owner:** `d3_splicing.py:584` cross-references `c1_cai.py:490-525` for
the `_unavailable` pattern; that method has moved (those lines are `lattice_terms`
now) and the range was already stale before this change. `rules/catalog/` is this
lane, but `d*` is another session's file in it, so the fix is flagged rather than made
— it should cite `CodonAdaptationIndex._unavailable` by name, not by line range.
