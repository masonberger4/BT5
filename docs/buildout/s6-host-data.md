# S6 — Host data and a real backbone

*Copy this whole file into a fresh Claude Code session on the BT5 repo.*

**Launch:** permission mode `default` · model opus · effort high · **needs a person
within reach** — every edit you make lands under `data/`, which prompts, and every PR
you open needs the `approved:data-change` label.
**Do not run this in plan mode** — five other sessions are running in parallel and the
gate is the draft PR, not a plan approval.

---

Your entire output is evidence. Getting the provenance wrong is the failure mode, not
getting the file format wrong.

## Read this first

Run `/bootstrap` — a fresh checkout has no `.venv`, and `gates.sh` exit **10** means
BROKEN, not a code failure. Then `CLAUDE.md`, then `docs/buildout/README.md`.

Your branch: **`claude/s6-host-data`**, cut from `main` (green at `628e130`).

## The situation

Exactly one codon-usage table ships:

```
data/codon_usage/sharp_li_1987_ecoli_w.json
data/vendors/templates/idt_eblocks_plate_96.xlsx
data/vendors/templates/_provenance.json
```

`data/genetic_codes/` **does not exist**, despite `CLAUDE.md` §2 protecting it — the
genetic code comes from Biopython at runtime. Worth confirming and reporting rather
than "fixing" by inventing files.

## What to build, in priority order

1. **Mammalian codon-usage reference sets (issue #78).** C1 (CAI) reports
   `unavailable` for all seven non-E. coli hosts because no reference set ships. S3 is
   building C3 (%MinMax) in parallel and will hit the same wall.

   `docs/decisions.md` already records why the obvious shortcuts were rejected, and
   you are bound by that reasoning:
   - *Falling back to the E. coli table for a mammalian host* — it is the one table on
     disk, so the fallback would always succeed and always be wrong: a
     plausible-looking number measuring nothing. This is the failure the `unavailable`
     path exists to prevent.
   - A human highly-expressed reference set is **an evidence-bearing decision with its
     own provenance burden**. That burden is now yours.

   `TableProvider.weights()` keys on the **reference set**, not on `HostId` — every
   caller passes `"sharp_li_1987_ecoli_w"`. Confirm how `FileTableProvider` resolves a
   key before choosing filenames.

2. **A real Addgene lentiviral map (issue #74).** `docs/PLAN.md:490-495`'s v1 bar wants
   "one real Addgene lentiviral GenBank with an annotated 5'UTR and an origin-spanning
   feature". Only synthetic fixtures exist —
   `tests/data/backbones/synthetic_lenti_ef1a.gb` and `synthetic_mcs_ef1a.gb`. The
   synthetic ones are sufficient for correctness but not for annotation quality, which
   is what PLAN's **G2** gate measures (lossless GenBank round-trip including an
   origin-spanning `join()`).

   PLAN says G2 should scare you most, "because a wrong coordinate model invalidates
   everything built on it". A real map is how you find out.

3. **Whatever `data/genetic_codes/` should be**, if anything. Report first; do not
   invent files to satisfy a path in `CLAUDE.md` §2.

## Provenance is the deliverable

Every file you add carries where it came from, when, and what it licenses.
`data/vendors/templates/_provenance.json` is the shape already in use — read it before
you write yours.

For a codon-usage table specifically: the organism, the reference gene set, the paper
or database release it derives from, the retrieval date, and enough detail that
`rule-auditor` can later answer *does the cited source actually support this?* A table
without that is not usable evidence no matter how correct its numbers are.

## Files

**You own:** `data/**` and `tests/data/backbones/**`.

**Never touch any engine source.** Not `rules/`, not `codon/`, not `design/`, not
`score/` — no `packages/engine/src/**` at all. Also never `core/`, `verify.py`,
`.github/`, `pyproject.toml`, `tests/contract/`, `tests/invariants/`,
`tests/data_integrity/`.

`data/` is a **global mutex** and you hold it. No other session writes there.

## Your contract with the other five

- **You ship data only.** Wiring a new reference set into `CAI_REFERENCE_SET` is an
  edit to `c1_cai.py`, which is **S3's file**. Do not make it.
- Instead: when your data lands, **open a follow-up issue** naming the reference sets
  and the hosts they serve, so S3 can wire them. If S6 merges before S3 gets to C1,
  S3 will pick it up directly.
- Rules report `unavailable` for a missing reference set. That is correct behaviour and
  your PR does not need to change it.

## Delegation

- `docs-miner` — what provenance a shipped reference set must carry, and PLAN's G2
  round-trip bar. Never read `PLAN.md` (58 KB) or `brief.md` (63 KB) inline.
- **`rule-auditor`** (opus/xhigh) on the acquired tables **before** the PR, not after.
  Its one question — *does the cited source actually support the number in the code?* —
  is exactly this session's question, so run it as a design check rather than waiting
  for review.
- `Explore` — how `FileTableProvider` keys a reference set, and every caller of
  `TableProvider.weights()`.
- `gate-runner` — gates before each push. `tests/data_integrity/` will exercise your
  files; you may **not** edit that directory to accommodate them.

## Done means

- `bash scripts/gates.sh` reaches `ALL GATES PASSED`. Exit 10 = `/bootstrap`; pytest
  exits 2, 3, 4, 5 are BROKEN, and 5 is never success.
- Every added file has provenance in the established shape.
- A real Addgene lentiviral map round-trips losslessly, origin-spanning `join()`
  included (PLAN's G2).
- No `packages/engine/src/**` file is in your diff. Check with
  `git diff --name-only origin/main`.
- A follow-up issue exists for the `CAI_REFERENCE_SET` wiring S3 will do.
- The PR is **open as a draft** and carries **`approved:data-change`**.
- You appended a `docs/decisions.md` entry: decided, **rejected** and why, with
  evidence. Newest-first, so expect a conflict and keep both entries.

**Do not self-merge.** `data/` is a protected path under `CLAUDE.md` §2, and §7b is
explicit that the `approved:*` label is sign-off on the change, **not** a licence to
merge it unreviewed. Goes to the owner.
