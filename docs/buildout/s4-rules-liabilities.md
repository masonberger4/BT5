# S4 — Rules: sequence liabilities (2.D, 2.E, 2.F)

*Copy this whole file into a fresh Claude Code session on the BT5 repo.*

**Launch:** permission mode `acceptEdits` · model opus · **effort xhigh** (override;
the repo default is high) · runs unattended.
**Do not run this in plan mode** — five other sessions are running in parallel and the
gate is the draft PR, not a plan approval.

---

You are filling the liability half of the rule catalog, and the first item on your list
is the hardest correctness problem in this buildout.

## Read this first

Run `/bootstrap` — a fresh checkout has no `.venv`, and `gates.sh` exit **10** means
BROKEN, not a code failure. Then `CLAUDE.md`, then `.claude/rules/rules-catalog.md`
(which governs your lane) and `docs/design/repeats.md`. Then
`docs/buildout/README.md`.

Your branch: **`claude/s4-rules-liabilities`**. Cut it yourself, from a
**freshly fetched** main — do not reuse a branch someone made earlier, and do not
assume main is where it was when this prompt was written:

```bash
git remote prune origin
git fetch -q origin main && git checkout -B claude/s4-rules-liabilities origin/main
```

`git fetch --prune origin main` does **not** clear a stale ref — pruning is bounded by
the refspec, so that form prunes only `origin/main` (CLAUDE.md §7a).

## What to build, in priority order

1. **D3 — splicing. Start here, and budget for it.** `CLAUDE.md` §3.6 mandates
   `RepairPolicy.FIXED_POINT` for splice-site removal:

   > A single pass ships a construct whose donors were removed *into* new donors, and
   > the validator passes it.

   **There is no splice rule in the catalog at all.** That rule was written to govern
   code nobody has written yet — you are writing it. MaxEntScan scoring, fixed-point
   iteration, and the V5-tag special case the brief documents.
2. **D2 — recombinase sites** (loxP / FRT / Gateway / Bxb1).
3. **D8 — CpG**, three separately-toggleable metrics (TLR9, ZAP/KHNYN, methylation
   islands).
4. **E3 — windowed GC** bands at 50/100 bp.
5. **F5 — the AT-window rule.**
6. **Issue #56** — no error-free length on file for gBlocks, so the default report
   cannot state a screening burden. `rules/_provenance.json` is your file.

Take them one PR at a time if large; every rule ships with its paired test in the same
PR.

## How to add a rule here

**Use `/rule-add`.** It scaffolds the rule with its paired test, its provenance and a
resolvable `brief_ref`, and it knows the eleven
`tests/data_integrity/test_rule_contract.py` assertions your rule must satisfy.

**Pull the threshold verbatim before you type it**, via `docs-miner`. Never paraphrase
a number, and never read `brief.md` (63 KB) inline. **Your lane contains the known
superseded row:** `brief.md:141` struck through E4's extent thresholds on 2026-08-28 as
*below the chance floor* — random 50% GC DNA has a 50 bp extent of 26.0 at 300 bp
rising to 46.0 at 10 kb, so a fixed cutoff is wrong at every length but one. The live
metric is `dGC` against its own binomial floor, evidence grade **B/contested**. If you
audit or extend `e4_gc_extent.py`, check it against the corrected formulation, not the
struck-through numbers.

## The rules that are not negotiable

From `CLAUDE.md` §3:

- **Never evaluate a rule against a bare string.** Rules take a `Construct` — that is
  what makes junction-, origin-spanning and reverse-strand hits impossible to miss.
  Your lane is where this bites hardest: repeats and inverted repeats wrap the origin.
- **Never scan the reverse strand yourself for motif rules.** List forward motifs in
  `LatticeTerms.forbidden` and let the solver close the set. Directional scored models
  are *not* revcomp-symmetric and must read `slot.strand_of_interest` — D3 splicing is
  directional, and so is D2.
- **Hard constraints are never enforced by a penalty weight.** `HARD_LATTICE`
  (guaranteed by the automaton), `HARD_REPAIR` (repair plus the independent validator,
  which refuses to emit) or `HARD_CHECK` (real but unfixable by codon choice).
  `default_weight` 0.0 for all three; `steering_weight` nudges the DP.
- **Splice-site removal must use `RepairPolicy.FIXED_POINT`** — see above. This is
  yours.
- **The genetic code table is explicit and never defaulted**; never emit a codon that
  is also a stop in the target table.
- **Seed every RNG explicitly** with `np.random.default_rng(seed)`. Global
  `np.random.*` and stdlib `random` are banned; CI greps for both.

## Files

**You own:** `packages/engine/src/bt5/rules/catalog/d*.py`, `e*.py`, `f*.py`;
`packages/engine/src/bt5/rules/vendors.py`, `rules/_provenance.json`,
`rules/fragment.py`, `rules/exempt.py`; and
`packages/engine/tests/rules/test_{d,e,f}*.py` plus the vendor tests.

**Never touch:** `b*` and `c*` rule files (S3 owns those), `score/`, `design/` (S1's),
`solver/`, `vector/`, `cassette/`, `codon/`, `structure/`, `core/`, `verify.py`,
`.github/`, `pyproject.toml`, `data/`, `tests/contract/`, `tests/invariants/`,
`tests/data_integrity/`.

**`packages/engine/tests/rules/conftest.py` is READ-ONLY for you**, and for the other
rules session too. It holds the shared helpers every rule test uses — `construct()`,
`wrapping_construct()`, `slot()`, `context()` and the `services` fixture. Both rules
sessions run at once, so an edit there is the one collision this split cannot absorb.
Need a helper it does not have? **Define it in your own test file.** A local helper
costs a few lines; a shared-conftest edit costs a merge conflict in the file every rule
test imports. If a helper genuinely belongs to both lanes, open an issue and let it land
after both sessions merge.

**Rule registration stays autodiscovery.** `core/registry.py` walks
`bt5.rules.catalog` with `pkgutil`; both `rules/__init__.py` and
`rules/catalog/__init__.py` are empty, so adding a rule edits **zero** shared files —
which is why you and S3 can work simultaneously. Do not introduce a hand-maintained
rule list.

## The cross-lane problem in D3

`RepairPolicy.FIXED_POINT` is enforced by the **solver's** repair seam, and you may not
edit `solver/`. So D3 is cross-lane by construction: your rule declares the policy,
the solver honours it.

**Use `/architect` for this** before writing D3 — it runs opus at ultracode and
orchestrates independent designs rather than answering alone. If the seam turns out to
need a solver change, that is an issue for the M1 lane, not an edit you make.

## Delegation

- **`/architect`** for D3's repair-policy seam, before implementing.
- **`/rule-add`** per rule.
- **`docs-miner`** for each `brief.md` row, verbatim, before any threshold is typed.
- **`rule-auditor`** (opus/xhigh) fires on your PR automatically via `/pre-pr` — it
  runs iff a Spec's `citations`, `weight_provenance`, `enforcement`, `last_verified`
  or a threshold changed, which is every PR you will open. Run
  **`/verify-provenance <rule-id>`** early instead of finding out at PR time.
- `gate-runner` — gates before each push.
- **`debugger`** if a fix you already tried did not work. It is explicitly the second
  attempt, not the first.

**Over the 20 KB bare-Read limit, in your lane:** `rules/vendors.py` (25.7 KB) and
`catalog/e4_gc_extent.py` (19.1 KB). Use `offset`/`limit`, or delegate to `Explore`.

## Done means

- `bash scripts/gates.sh` reaches `ALL GATES PASSED`. Exit 10 = `/bootstrap`; pytest
  exits 2, 3, 4, 5 are BROKEN, and 5 is never success.
- Every rule you added has a paired test at
  `packages/engine/tests/rules/test_<id>.py`.
- D3 removes splice donors to a **fixed point**, with a test that fails under a single
  pass — the exact failure `CLAUDE.md` §3.6 describes.
- `rule-auditor` returns SUPPORTED for every threshold you shipped.
- The PR is **open as a draft**.
- You added a decision file at `docs/decisions/2026-XX-XX-<slug>.md`: decided,
  **rejected** and why, with evidence. One file per decision — never append to a
  shared one.

- **`/pre-pr` is run by the operator, not by you.** It is
  `disable-model-invocation: true`, so a session cannot self-invoke it and must not
  replicate its steps by other means. Ask for it when the branch is ready.
- **The attestation is posted last.** After `/pre-pr` and after the final push, comment
  the full 40-character head SHA on the PR:

  ```
  /pre-pr <head-sha>
  ```

  The advisory `pre-pr-attest` check reads that comment. An attestation names **one**
  commit, and pushing again makes it stale on purpose — a review of the previous tree
  says nothing about this one. Never attest a SHA that was not just reviewed; the whole
  value is that the claim is on the record. If a gate or review came back blocking and
  you are pushing anyway, do **not** attest — say so in the PR and let the check stay
  red. Only the owner may waive it, with `/pre-pr-bypass <head-sha>`.

**Do not self-merge.** Adding a rule changes the sequences the app produces — non-"none"
scientific impact, so `CLAUDE.md` §7b sends it to the owner. Fill in the PR template's
scientific-impact section properly: the evidence, and what enforcing the rule costs on
the other objectives.
