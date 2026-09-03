# Code-reviewer notes for BT5

## Enforcement is per-slot, not per-class (score/presets.py)
`Spec.enforcement` (ClassVar) is only the FLOOR. The real routing is
`enforcement_for(slot) -> Enforcement`, which can escalate SOFT to HARD_REPAIR/
HARD_CHECK depending on `slot.modality` (e.g. `d4_internal_polya`: SOFT floor,
HARD_REPAIR on LENTIVIRAL/AAV/GENOME_INTEGRATED). Any guard in `score/` that
enforces CLAUDE.md §3.5 ("hard constraints never carry a weight") must ask
`enforcement_for` per slot admitted by the preset's `modality`, not read the
ClassVar alone — reading the ClassVar alone is exactly the bug fixed in #72
(LENTIVIRAL/AAV both shipped `WeightEntry("2.D4", ...)` past a ClassVar-only
guard). Watch for this pattern recurring anywhere a preset/objective assembly
step decides "is this rule scorable" — check it asks per-slot.

## `_slots_admitted_by`-style enumeration is deliberately over-broad
`packages/engine/src/bt5/score/presets.py::_slots_admitted_by` builds every
`(host, table_id) x SlotRole` pair for a fixed modality, not just the "real"
host for that modality. Intentional (documented in-file): a guard should err
toward refusing a weight rather than guessing a representative host.
Re-check this if a future rule starts gating on `slot.host`/`slot.role`
instead of only `slot.modality`, since over-broad enumeration could then
cause spurious refusals.

## `spec()` instantiation inside a guard (try/except TypeError)
`_unscored_enforcement` tries `spec()` to call `enforcement_for` per slot; if
the spec needs required constructor args it falls back to the ClassVar floor
only. Dead code today (every catalog rule's `__init__` has all-default args)
but a silent hole if a future rule adds a required constructor arg while
relying on `enforcement_for` to escalate. Worth a non-blocking mention if it
recurs, not blocking on its own.

## Test inversions aren't automatically suppression
Flipping an assertion from "preset DOES weight X" to "does NOT weight X" is
legitimate when the old assertion was pinning the bug itself and prod code
changed first to fix it, especially with a new invariant test class pinning
the property going forward. Suppression loosens a check to dodge a
still-present failure; distinguish the two by whether prod code actually
changed.

## docs/decisions/ scientific-impact flagging
A lane PR that changes ranking/weight behavior in a shipped preset should say
so explicitly under "Scientific impact" (not "none"), triggering §7b's
owner-merge requirement. Treat its absence on a weight/ranking-affecting diff
as a finding.

## A thin adapter's "reportable errors" tuple is an allowlist, not a catch-all
`bt5/cli.py`'s `_REPORTABLE_ERRORS` misses raise sites reachable from `bt5
design` (e.g. `FileTableProvider.genetic_code()`'s bare `ValueError`,
Biopython `SeqIO.read()` on a malformed `--backbone` file) — both escape as
raw tracebacks, though exit code stays nonzero. Non-blocking UX gap, not a
CLAUDE.md violation. Check the tuple against what callees one level of
transitive dependency down actually raise, not just the caller's docstring.

## Reviewing prompt/doc diffs (docs/buildout/*.md session prompts)
This repo's docs culture is precise — file:line citations, quoted strings,
exported-symbol lists and SHAs have consistently checked out exactly against
source. Recurring gap classes worth the verification budget:
- diff each prompt's "never touch" list against the others and against
  `ls packages/engine/src/bt5/` (9 subdirs), not just the README table;
- check a shared fixture one level up from a glob-partitioned test split is
  declared read-only by every session whose glob doesn't cover it;
- diff a decision doc's plural-attribution sentences ("S2 and S4 do X")
  against each named prompt individually — a plural claim applying to only
  one file is the recurring failure shape;
- cross-check issue/PR numbers only against `git log --oneline` squash-merge
  `(#N)` trails (no GitHub API here); say so when a citation has no trail.

## "read-only agent" is policy here, not mechanism — do not claim otherwise
`rule-auditor` holds `Agent` (to resolve `brief_ref` via `docs-miner`), which
also reaches `batch-editor` (holds `Edit`). Do NOT write this up as "the only
agent whose read-only-ness is not mechanical" — `code-reviewer`, `debugger`,
`docs-miner`, `gate-runner`, `security-reviewer` all hold `Bash`, which writes
via `sed -i` too. The `Agent` grant widens an existing surface, it doesn't
create a new class. `tools: Agent(docs-miner)` (parameterised) is dropped
silently by this CLI, leaving no `Agent` tool at all (fails closed, doesn't
scope). An agent holding `Bash` can reach the `claude` CLI on PATH regardless
of its tool list. Before calling a control "mechanical", name the mechanism
and check it holds. When a diff corrects a fact, grep the whole file for the
old value — a stale citation can survive in an example block below the fix.

## Rule-count/repair-count prose in `.claude/*.md` is a recurring drift point — verify, don't assume loss of a number is a defect
This repo has repeatedly patched `.claude/agents/`, `.claude/rules/`, `.claude/skills/`
prose that hardcoded "25 catalog rules" or "22 of 25 declare SINGLE_PASS" style counts
(574ea0e stopgap 15->25; #103 replaced a hardcoded count with an enumeration; a later
diff removed the "25" count and the "d3_splicing/b9_out_of_frame_atg/f5_at_window"
FIXED_POINT enumeration entirely, replacing with structural explanations). When review-
ing such a diff: (1) grep `git log --all -p -- <path>` for any historical number the
prose cites as "how long this bug existed" before calling it fabricated — the true prior
state can be one commit further back than a three-dot diff shows; (2) verify every
line:number and mechanism citation in the *replacement* text against live source (e.g.
`core/spec.py:231`, `solver/repair.py:417`, `solver/catalog.py`'s `policies()`) rather
than trusting the count it removed was accurate — sometimes the old enumeration was
itself already stale (e.g. a "four brief.md sections" list that omitted `2.C` even
though `c1_cai`/`c3_min_max` already cited it); (3) losing a concrete example (which
rules currently declare `FIXED_POINT`) in exchange for a rule that never goes stale is a
legitimate trade, not a diligence failure, as long as a canonical example survives
elsewhere (e.g. "splice donors are the canonical case" stays in `rule-add/SKILL.md`).
`.claude/agents/`, `.claude/rules/`, `.claude/skills/` are not in CLAUDE.md's lane table
or its §2 protected-path list — editing their prose needs no lane issue and no `approved:*` label.

## "this code path never ran in production" claims: verify by chasing the caller chain to a filesystem/data mismatch, not just by reading the guard
PR #112 (`score/null.py`, host-frequency weighted null draw, merged as `10456d8`)
claimed the weighted
branch of `synonymous_variant`/`null_distribution` never executed because no host
resolves a codon-usage table (issue #98) — checked out true: `design/runner.py::
_host_usage` looks up `data/codon_usage/{str(HostId)}.json` (e.g. `human.json`,
`e_coli_k12.json`), but the files actually on disk are named by reference-set, not
host (`human_highly_expressed_refseq_w.json`, `sharp_li_1987_ecoli_w.json`, etc. —
a DIFFERENT lookup used by C1/C3's `CAI_REFERENCE_SET`/`MINMAX_REFERENCE_SET`
maps). `_host_usage` is the only production caller feeding `null_distribution`, so
every host raises `FileNotFoundError` there today and the weighted path is dead
code. General pattern: when a PR's "scientific impact: none" / agent-mergeability
argument rests on "this branch is unreachable", grep every production caller (not
just tests) of the function that would reach it, and check any filename/lookup-key
assumption against `ls` of the actual data directory — a two-path split (one
correct lookup key elsewhere in the codebase, one stale one in the code under
review) is exactly the shape that makes an "unreachable" claim true today but
fragile to a future data-file rename.

## Inverse-CDF / `bisect_right`-on-`rng.random()` draws: verify empirically, it's cheap and conclusive
For a PR replacing `rng.choice(n, p=weights)` with `bisect_right(cumulative_cdf,
rng.random(), hi=len(cumulative)-1)`: this is provably the same categorical
distribution (bisect_right's tie-breaking toward higher indices means a
zero-width interval from a zero-weight option can never be selected, and the
`hi=len-1` cap only guards the rounding edge case where cumulative floating-point
division leaves the last entry a hair under 1.0). Don't just eyeball it — a
`.venv/bin/python` scratch script (no repo files touched) computing chi-square
over ~100-400k draws both ways, plus checking `rng.bit_generator.state` equality
after a single `rng.choice(n, p=p)` vs a single `rng.random()` call with the same
seed, settles distribution-equivalence and stream-consumption questions in under a
minute and is far more conclusive than reasoning about numpy's C internals from
the docstring alone. In this PR they came out bit-identical for the tested weight
vectors, stronger than the PR's own hedged "same distribution, exact stream
changes" claim.

## "verified rather than asserted" claims need a matching pinned test
PR #112's *PR description* claimed a 400k-draw chi-square check and "zero-weight
codons drawn 0/60k" and "all-zero family falls back to uniform" as verification
(check with `gh pr view`, not `git show` — the squash commit `10456d8` carries only
the follow-up line, so grepping the commit body for these will find nothing) — none of
the three is a committed test in `tests/score/test_null.py`; the only weighted-path
test (`test_weights_actually_bias_the_sampling`) uses a fixture (`usage` in
`tests/score/conftest.py`) where every codon has weight 1.0 or 10.0, never 0.0, so
the zero-weight-skip and all-zero-fallback branches of `weight_table()` are
completely untested in CI. Flag this as non-blocking test-coverage gap (not a
CLAUDE.md §4 suppression — nothing was skipped/loosened, coverage was just never
added), but call it out explicitly since the PR description implies these were
checked and locked in. #112 merged with the gap open; as of `369429a` the fixture
is still `10.0 if c.endswith("C") else 1.0` (`tests/score/conftest.py:27`), so
`running > 0` is always true and `null.py:133`'s `else None` arm has still never
run — a ready-made first test for whoever turns host steering on.
