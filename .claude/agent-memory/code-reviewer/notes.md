# Code-reviewer notes for BT5

## Enforcement is per-slot, not per-class (score/presets.py)
`Spec.enforcement` (ClassVar) is only the FLOOR. Real routing is
`enforcement_for(slot) -> Enforcement`, which can escalate SOFT to HARD_REPAIR/
HARD_CHECK depending on `slot.modality` (e.g. `d4_internal_polya`: SOFT floor,
HARD_REPAIR on LENTIVIRAL/AAV/GENOME_INTEGRATED). Any guard enforcing CLAUDE.md
§3.5 ("hard constraints never carry a weight") must ask `enforcement_for` per
slot admitted by the preset's `modality`, not read the ClassVar alone — that
was the #72 bug (LENTIVIRAL/AAV both shipped `WeightEntry("2.D4", ...)` past a
ClassVar-only guard).

**`gate is None` WIDENS the probe, it does not grant a free pass.**
`_unscored_enforcement`: `if gate is not None and not gate(slot): continue`. A
missing `gate` short-circuits the `and` to True, so `continue` never fires and
*every* slot gets checked against `enforcement_for` — strictly more likely to
find a barring Enforcement, not less. Only a missing/None `enforcement_for` is
the real silent-permission-to-weight hole (guard returns `None` = "scored").
`gate` is still worth asserting separately — Spec Protocol conformance, and a
non-callable one raises `TypeError` inside `resolve()` — just not for the
"silent weight" reason. Confirmed live against `presets.py`/`core/spec.py` on
0598071 (#129); an earlier commit on that branch had conflated the two
rationales, correctly flagged, then fixed.

`_slots_admitted_by` deliberately enumerates every `(host, table_id) x
SlotRole` for a modality rather than guessing one representative host — errs
toward refusing a weight. Re-check if a rule starts gating on `slot.host`/
`slot.role` instead of only `slot.modality`.

`_unscored_enforcement` falls back to the ClassVar floor if `spec()` raises
`TypeError` (required constructor arg) — dead code today, every catalog rule
is default-constructible. Pinned per-spec since #82 by
`tests/data_integrity/test_rule_contract.py::TestRuleContract::
test_is_probeable_per_slot`; check that test rather than re-deriving.

## Near-duplicate `resolve()`-over-`PRESETS` guards: check the buildout plan before flagging as accidental
`packages/engine/tests/score/test_presets.py::TestShippedWeights::
test_no_shipped_preset_weights_a_hard_rule_in_this_build` (since #72) already
does `for preset in PRESETS: resolve(preset)`. PR #129 added a functionally
identical `test_every_shipped_preset_resolves_against_the_live_catalog` to
`tests/data_integrity/test_rule_contract.py` — looked like accidental
duplication into a §2 protected path (`approved:oracle-change`) for free
coverage the score lane's own suite already had. It wasn't: `docs/buildout/
wave2/w2-design-path-truth.md` ("Blocked-by: issue #82") explicitly specs this
exact second `data_integrity` assertion — "all three shipped presets
`resolve()` without raising against the live catalog... the check the wave
README promises" W3/W4/W5 (rule-adding sessions in OTHER lanes, who would have
no reason to read `score/`'s own test file). Real, non-vacuous, and
deliberately placed for defense-in-depth across lanes that can't be expected
to check each other's test suites. Grep `docs/buildout/**` for the issue
number before calling this kind of thing redundant.

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
raw tracebacks, though exit code stays nonzero. Non-blocking UX gap. Check the
tuple against what callees one level down actually raise, not just the
caller's docstring.

## Reviewing prompt/doc diffs (docs/buildout/*.md session prompts)
This repo's docs culture is precise — file:line citations, quoted strings,
exported-symbol lists and SHAs have consistently checked out exactly against
source. Recurring gap classes worth the verification budget: diff each
prompt's "never touch" list against the others and against `ls
packages/engine/src/bt5/`; check a shared fixture one level up from a
glob-partitioned test split is declared read-only by every session whose glob
doesn't cover it; diff a decision doc's plural-attribution sentences ("S2 and
S4 do X") against each named prompt individually; cross-check issue/PR numbers
only against `git log --oneline` squash-merge `(#N)` trails (no GitHub API
here) and say so when a citation has no trail.

## "read-only agent" is policy here, not mechanism — do not claim otherwise
`rule-auditor` holds `Agent` (to resolve `brief_ref` via `docs-miner`), which
also reaches `batch-editor` (holds `Edit`). Don't write this up as "the only
agent whose read-only-ness is not mechanical" — `code-reviewer`, `debugger`,
`docs-miner`, `gate-runner`, `security-reviewer` all hold `Bash`, which writes
too. `tools: Agent(docs-miner)` (parameterised) is dropped silently by this
CLI (fails closed). Before calling a control "mechanical", name the mechanism
and check it holds. When a diff corrects a fact, grep the whole file for the
old value — a stale citation can survive in an example block below the fix.

## Rule-count/repair-count prose in `.claude/*.md` drifts — verify, don't assume a removed number is a defect
This repo has repeatedly patched hardcoded counts ("25 catalog rules", "22 of
25 declare SINGLE_PASS") in `.claude/agents|rules|skills`. When reviewing such
a diff: grep `git log --all -p -- <path>` for the number's history before
calling it fabricated; verify every line:number/mechanism citation in the
*replacement* text against live source rather than trusting the removed count
was accurate (it sometimes wasn't); losing a concrete example in exchange for
a rule that never goes stale is a legitimate trade if a canonical example
survives elsewhere. `.claude/agents|rules|skills` are outside CLAUDE.md's lane
table and §2 protected-path list — editing their prose needs no lane issue and
no `approved:*` label.

See also [[score-null-pr112]] for verification techniques (empirical RNG
checks, caller-chain tracing, PR-description-vs-committed-test gaps) from the
`score/null.py` weighted-draw review.
