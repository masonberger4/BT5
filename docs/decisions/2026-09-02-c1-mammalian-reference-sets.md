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

**Measured, not assumed.** Against the in-band fixture the four new hosts return CAI
0.807 (human, HEK293), 0.788 (mouse) and 0.801 (CHO) versus E. coli's 0.763 — distinct
numbers from distinct tables, all inside the (0.70, 0.90) band. A test asserts a
mammalian CAI is *not* the E. coli number, which is the failure the unavailable path
existed to prevent and which only becomes testable now that a real alternative exists.

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

**Where:** branch `claude/s3-c1-mammalian-hosts`.
