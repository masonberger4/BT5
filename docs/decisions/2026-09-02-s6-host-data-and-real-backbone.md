## 2026-09-02 — S6: mammalian codon-usage reference sets, a real lentiviral map, and the empty genetic_codes path

Buildout session S6 (`data/` lane). Closes #78 (data half) and #74 (acquisition half);
opens a follow-up for the `CAI_REFERENCE_SET` wiring S6 may not make.

### Decided — mammalian codon-usage reference sets (#78)

Shipped three highly-expressed reference sets under `data/codon_usage/`, keyed by
reference set (not `HostId`, per `FileTableProvider.usage` → `{key}.json`):

| file / key | organism | serves hosts |
|---|---|---|
| `human_highly_expressed_refseq_w` | *Homo sapiens* (9606) | HUMAN, HEK293 |
| `mouse_highly_expressed_refseq_w` | *Mus musculus* (10090) | MOUSE |
| `cho_highly_expressed_refseq_w` | *Cricetulus griseus* (10029) | CHO |

HEK293 → human is the load-bearing one: both shipped mammalian presets
(`lentiviral_hek293`, `aav_hek293`) key on HEK293, so a human set is what lets C1 score
instead of report `unavailable` for the presets the walking skeleton actually runs.

**Method.** `w` is built from a fixed, auditable panel of ~100 canonical highly-expressed
genes — cytosolic ribosomal proteins (large + small subunit), translation elongation
factors, core glycolytic enzymes, and chaperonins — extending the fallback reference set
named in `brief.md:282`. Resolution is **gene-first and verified**: for each symbol the
build finds the NCBI Gene record whose *official* symbol equals it (never an alias), links
to that gene's RefSeq RNA (curated `NM_` preferred over predicted `XM_`), and — after
fetching — re-checks that the transcript's `/gene` equals the requested symbol before
counting. A symbol that resolves to no gene, or to a transcript for a different gene, is
**dropped and disclosed** (`unresolved_symbols` / `mismatched_symbols`), never mapped to a
paralog. Each accession is counted once. The CDS's in-frame sense codons are counted
(terminal stop and non-ACGT CDS excluded) and `w = (count + 0.5) / family_max` per family
— **identical to `bt5.codon.tables.CodonUsage.from_counts`** (verified: matches to 6 dp on
random counts). Genetic code table 1 (nuclear) for all three. The committed build script
`data/codon_usage/build_reference_set.py` is the reproducible recipe; each JSON's
`_provenance` lists every contributing accession *and its verified gene*, the `NM_`/`XM_`
split, unresolved/mismatched symbols, rejected records, and the source `sha256`.

**A first pass got this wrong, and `rule-auditor` caught it — which is why it was run
before the PR, not after.** The original resolver searched `{symbol}[gene]` in nuccore and
took the first hit without verifying identity, so human `RPL10` silently resolved to
RPL15's transcript and three mouse symbols collapsed onto paralogs — inflating
`genes_contributing` and double-counting a handful of transcripts. The gene-first, verified
resolver above is the fix; the tables here are the regenerated, verified output.

Result: **human 99/99** genes (all `NM_`, 25 846 codons); **mouse 96** (all `NM_`, 24 747
codons; 3 unresolved — `RPS3A`, `TUBB`, `FTL`, whose mouse orthologs carry different
*official* symbols `Rps3a1`/`Tubb5`/`Ftl1`, so the strict matcher refuses to guess);
**CHO 91** (20 `NM_` + 71 predicted `XM_`, 23 903 codons; 4 unresolved, and 1 mismatch —
`RPS19` → an uncharacterized `LOC` model — correctly dropped). All load through
`FileTableProvider` with 61 sense codons, no stops, `ATG`/`TGG` = 1.0, an https citation,
and every row's symbol equal to its transcript's gene. Known mammalian bias reproduces
(Leu `CTG` = 1.0, Ala `GCC` = 1.0, `GCG` rare ≈ 0.15–0.17).

### Rejected — codon-usage

- **PaxDb v5/v6 top-5%-by-abundance as the primary source (`brief.md:282`).** PaxDb is the
  brief's *primary*, the curated panel its *fallback*. Reaching PaxDb's number requires
  mapping abundance IDs (STRING/Ensembl) to CDS nucleotide sequences — a multi-step
  pipeline whose mis-mappings produce exactly the "plausible-looking number measuring
  nothing" this project refuses, and one a reviewer cannot audit from the shipped file.
  The curated highly-expressed panel is the brief's own sanctioned fallback, is the
  classical CAI reference method (Sharp & Li built theirs the same way), and is auditable
  down to the accession. Chosen deliberately over the primary for that auditability.
- **Falling back to the E. coli table for a mammalian host.** Already rejected in
  `2026-09-01-c1-cai-soft-band.md`; restated because it is the failure the whole
  `unavailable` path exists to prevent. Not done.
- **Shipping S_CEREVISIAE, P_PASTORIS, SF9 now.** No shipped preset consumes them, and
  their RefSeq symbol/transcript coverage is materially messier than the vertebrate hosts
  (non-standard symbols, sparse curated transcripts). Rushing them would trade the one
  thing this lane cannot trade — provenance integrity — for breadth no preset needs.
  Deferred to the follow-up; the build script generalizes to them by taxid.
- **Wiring `CAI_REFERENCE_SET` myself.** That map lives in `c1_cai.py` (S3's engine file);
  the data lane ships data only. Follow-up issue instead (see below).
- **A separate `data/codon_usage/_provenance.json`.** The shipped `sharp_li` file embeds
  its provenance inline; matched that shape rather than introducing a sidecar.

**Caveat carried on CHO:** 71 of 91 contributing transcripts are RefSeq *predicted*
(`XM_`) models, because *C. griseus* RefSeq is largely gene-prediction based. The counts
are made explicit in the file's `_provenance` (`predicted_xm_transcripts: 71`). Predicted
housekeeping CDSs are reliable at the codon level and 23 903 codons is ample, but a
reviewer who wants CHO on curated transcripts only should say so — it is a deliberate,
disclosed choice, not an oversight.

### Decided — a real lentiviral map (#74)

Committed `tests/data/backbones/real_lenti_pFTMGW_EF177827.gb`: NCBI GenBank **EF177827.1**
(pFTMGW, a FUGW-derived third-generation lentiviral transfer vector), 8928 bp, circular,
**unmodified**. It carries the real-world annotation the synthetic fixtures smooth over —
real LTRs, a real WPRE, misc_feature soup, an eGFP CDS with a depositor's quirky
`/transl_table=11`. Verified to round-trip losslessly through `bt5.vector`: sequence,
topology and every feature location/key/qualifier preserved, writer idempotent (PLAN G2's
semantic round-trip).

### Rejected — real map

- **An Addgene-authored map kept out of the repo behind an env var (#74's own proposal).**
  That workaround exists to dodge Addgene redistribution ambiguity. An NCBI GenBank deposit
  of the same vector class is public domain, so the fixture can be committed and checked
  in-repo — strictly better provenance than a file that isn't there.
- **Hand-adding a 5'UTR feature or re-linearizing to force an origin-spanning `join()`.**
  Both would corrupt the map's entire value as an *unmodified, independent* parse target
  and put fabricated annotation on a real accession. That geometry is already covered,
  cleanly, by `synthetic_lenti_ef1a.gb`; PLAN #1's origin-spanning bar stays satisfied
  there. This fixture's job is annotation realism, and it is left exactly as deposited.
- **Committing a second, richer map (PZ267690.1, a 31-feature pLenti clone).** Also
  round-trips losslessly and stresses the parser harder, but the brief asks for one real
  map and EF177827 matches #74's named requirements (real LTRs + WPRE) most directly. Noted
  here as the ready alternative if a heavier annotation-quality fixture is wanted later.

### Reported — `data/genetic_codes/` (do NOT fabricate)

`CLAUDE.md` §2 protects `data/genetic_codes/**`, but the directory **does not exist** and
nothing reads it: the genetic code comes from Biopython at runtime
(`Bio.Data.CodonTable.unambiguous_dna_by_id`, `codon/tables.py:64,181`); a repo-wide grep
for `genetic_codes` finds only the `CLAUDE.md` §2 line itself. No files were invented to
satisfy the path — the protection is defensive/aspirational and correctly left empty.

### Follow-up

Opened an issue for S3 to add HUMAN/HEK293/MOUSE/CHO to `CAI_REFERENCE_SET` in
`c1_cai.py` and update its pin test, naming the reference-set keys and the hosts each
serves. Until then C1 keeps reporting `unavailable` for these hosts, which is correct
behaviour, not a regression.

**Where:** branch `claude/s6-host-data`. PR carries `approved:data-change`; goes to the
owner (a protected path is sign-off, not a self-merge licence — §7b).
