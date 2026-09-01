# S3 — Rules: translation and expression (2.B, 2.C)

*Copy this whole file into a fresh Claude Code session on the BT5 repo.*

**Launch:** permission mode `acceptEdits` · model opus · effort high (both the repo
default, so no override needed) · runs unattended.
**Do not run this in plan mode** — five other sessions are running in parallel and the
gate is the draft PR, not a plan approval.

---

You are filling the expression half of the rule catalog. Sixteen of roughly forty-five
planned rules exist; ten of eleven 2.B rules and nine of ten 2.C rules do not.

## Read this first

Run `/bootstrap` — a fresh checkout has no `.venv`, and `gates.sh` exit **10** means
BROKEN, not a code failure. Then `CLAUDE.md`, then `.claude/rules/rules-catalog.md`
(which governs your lane specifically), then `docs/buildout/README.md`.

Your branch: **`claude/s3-rules-translation`**, cut from `main` (green at `628e130`).

## What to build, in priority order

1. **C3 — %MinMax. Start here.** All three presets already carry
   `WeightEntry("2.C3", 0.3, _NATIVE_NOTE)` (`score/presets.py:328,354,388`) and C3 is
   the only id left in `ResolvedPreset.unimplemented`. The preset machinery is
   *waiting* for this rule; shipping it closes a declared gap rather than opening a
   new one.
2. **B8 — Kozak.**
3. **B9 — out-of-frame ATG.**
4. **B2 — structure penalty**, STR(−30:+30)/(+31:+90).

Take them one PR at a time if they are large, or batch the small ones — but every rule
ships with its paired test in the same PR.

## How to add a rule here

**Use `/rule-add`.** It is the skill built for this: it scaffolds the rule with its
paired test, its provenance and a resolvable `brief_ref`, and it knows the eleven
`tests/data_integrity/test_rule_contract.py` assertions your rule must satisfy. Hand-
writing a rule file means rediscovering those assertions one CI cycle at a time.

**Pull the threshold verbatim before you type it.** Use `docs-miner` to get the exact
`brief.md` row — never read `brief.md` (63 KB) inline, and never paraphrase a number.
At least one row in the brief is struck through and superseded (E4's extent
thresholds, corrected 2026-08-28 as below the chance floor), so "the brief says" is
not good enough; you need the current row with its line number.

## The rules that are not negotiable

From `CLAUDE.md` §3, with the rationale in `.claude/rules/rules-catalog.md`:

- **The genetic code table is explicit and never defaulted.** NCBI table 12 reassigns
  CTG to Ser; table 4 makes TGA Trp. A wrong table is a silently wrong protein no
  assay catches for months.
- **Never emit a codon that is also a stop in the target table.** Tables 27 and 28 make
  TGA both Trp and a stop.
- **Never evaluate a rule against a bare string.** Rules take a `Construct` — that is
  what makes junction-, origin-spanning and reverse-strand hits impossible to miss.
- **Never scan the reverse strand yourself for motif rules.** List forward motifs in
  `LatticeTerms.forbidden` and let the solver close the set. Directional scored models
  are *not* revcomp-symmetric and must read `slot.strand_of_interest`. B9
  (out-of-frame ATG) is directional — this matters to you.
- **Hard constraints are never enforced by a penalty weight.** `HARD_LATTICE`,
  `HARD_REPAIR` or `HARD_CHECK`, with `default_weight` 0.0 for all three;
  `steering_weight` nudges the DP.
- **Seed every RNG explicitly** with `np.random.default_rng(seed)`. Global
  `np.random.*` and any stdlib `random` import are banned and CI greps for both.

Every SOFT rule explains its default weight in `weight_provenance`. Every rule carries
citations, an evidence badge and `last_verified`.

## A host-data caveat

C1 (CAI) reports `unavailable` for all seven non-E. coli hosts because only
`data/codon_usage/sharp_li_1987_ecoli_w.json` ships (issue #78). **C3 %MinMax will hit
the same wall.** That is correct behaviour, not a bug to route around: report
`unavailable` with a reason rather than falling back to the E. coli table, which would
always succeed and always be wrong — a plausible-looking number measuring nothing.
`docs/decisions.md` records this reasoning for C1; follow it.

S6 is acquiring mammalian reference sets in parallel. **If S6 has merged by the time
you get there**, wire the new sets into `CAI_REFERENCE_SET` in `c1_cai.py` (your file).
If not, leave it and open a follow-up issue.

## Files

**You own:** `packages/engine/src/bt5/rules/catalog/b*.py` and `c*.py`, plus
`packages/engine/tests/rules/test_b*.py` and `test_c*.py`.

**Never touch:** `d*`, `e*`, `f*` rule files (S4 owns those), `rules/vendors.py`,
`rules/_provenance.json`, `rules/fragment.py`, `rules/exempt.py` (S4's), `score/`,
`design/` (S1's), `solver/`, `vector/`, `cassette/`, `core/`, `verify.py`, `.github/`,
`pyproject.toml`, `data/`, `tests/contract/`, `tests/invariants/`,
`tests/data_integrity/`.

**Rule registration stays autodiscovery.** `core/registry.py` walks
`bt5.rules.catalog` with `pkgutil`; `rules/__init__.py` and `rules/catalog/__init__.py`
are both empty. Adding a rule edits **zero** shared files — which is exactly why you
and S4 can work at the same time. Do not introduce a hand-maintained rule list; that
would create the one collision this whole design avoids.

## Delegation

- **`/rule-add`** per rule. This is the main event.
- **`docs-miner`** for each `brief.md` row, verbatim, before any threshold is typed.
- **`rule-auditor`** (opus/xhigh) fires on your PR automatically — `/pre-pr` runs it iff
  the diff changes a Spec's `citations`, `weight_provenance`, `enforcement`,
  `last_verified` or a threshold, which describes every PR you will open. Run
  **`/verify-provenance <rule-id>`** early, while the number is still cheap to change,
  rather than discovering the problem at PR time.
- `Explore` — how an existing SOFT scored rule is shaped. `c1_cai.py` is the closest
  model for C3 and is the most recently reviewed rule in the tree.
- `gate-runner` — gates before each push.

No file in your lane exceeds 20 KB.

## Done means

- `bash scripts/gates.sh` reaches `ALL GATES PASSED`. Exit 10 = `/bootstrap`; pytest
  exits 2, 3, 4, 5 are BROKEN, and 5 is never success.
- Every rule you added has a paired test at
  `packages/engine/tests/rules/test_<id>.py`. `/pre-pr` checks this mechanically;
  `d1_restriction_sites` is the one known pre-existing gap and anything else is
  blocking.
- `rule-auditor` returns SUPPORTED for every threshold you shipped.
- `ResolvedPreset.unimplemented` is empty once C3 lands.
- The PR is **open as a draft**.
- You appended a `docs/decisions.md` entry: decided, **rejected** and why, with
  evidence. Newest-first, so expect a conflict and keep both entries.

**Do not self-merge.** Adding a rule changes the sequences the app produces, so your
scientific impact is non-"none" and `CLAUDE.md` §7b sends it to the owner. Fill in the
PR template's scientific-impact section properly: what is the evidence, and what does
enforcing this rule cost on the other objectives?
