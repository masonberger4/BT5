## 2026-09-02 — The biosecurity screen degrades to `not_run`, never `clear` (S2)

**Context.** `core/context.BiosecurityVerdict` was already frozen with
`status: Literal["not_run","clear","flag","block"]` and the docstring rule "status is
never 'clear' when screening did not run". `commec>=0.2` was already declared in the
`screen` extra. This lane makes the screen real behind a protocol and gets the
degradation right — the part CI actually runs, because `/bootstrap` installs
`dev,fold,export` and **not** `screen`, so commec is absent in CI.

**Decided:**

- **A `Screen` protocol in `cassette/screen.py`, mirroring `FoldEngine`.** Concrete
  `CommecScreen` (commec-backed) and `NullScreen` (reports `not_run`). `load_screen()`
  returns a `CommecScreen` only when commec is installed AND a database is configured,
  otherwise a `NullScreen` carrying the specific reason. Unlike `load_fold_engine`,
  `load_screen` never returns `None`: a screen must always produce a verdict, and the
  honest verdict when no real screen can run is `not_run`, which `NullScreen` embodies.

- **The degradation is fail-safe by construction.** `_status_for` maps commec's outcome
  vocabulary through an explicit table and **defaults every unrecognised outcome to
  `not_run`**. So commec absent, no database, a timeout, a non-zero exit, an unparseable
  result, or a future outcome word BT5 does not know — all land on `not_run`, never
  `clear`. The load-bearing test asserts the negative: `test_only_explicit_clear_reads_clear`
  (Hypothesis over arbitrary text) proves nothing maps to `clear` except commec's own
  "Clear"; `test_no_failure_mode_reads_clear` proves no exception path does either.

- **commec outcome → BT5 status mapping** (the one decision here that changes what the app
  REFUSES to build; grounded in commec's documented Clear/Warning/Flag vocabulary):
  `Clear → clear`, `Warning → flag` (surfaced, still emits), `Flag → block` (regulated
  pathogen match, refuses to emit). Conservative: commec's strongest outcome is the only
  one that refuses emission.

- **`database_version` recorded on every real verdict.** Read from commec's JSON, falling
  back to a version stamp beside the databases. A verdict from `NullScreen` carries
  `None`, because nothing produced it.

- **`guard_emission(verdict)` honours the block path** by raising `BiosecurityBlockedError`
  when `not verdict.may_proceed`. An exception, not a boolean, so a caller cannot forget
  to check it. It enforces only the one refusal the frozen type defines and has no
  argument that could turn a block into a pass.

- **The commec invocation is an injectable seam (`CommecScreen.runner`).** The
  status-mapping and degradation logic — the part with the correctness risk — is pure and
  fully tested without commec installed. The commec-specific argv/JSON handling
  (`_run_commec_cli`, `_parse_commec_output`) is quarantined into clearly-marked functions
  CI never reaches, written defensively so anything unexpected raises and becomes
  `not_run`.

**Rejected:**

- *`load_screen` returning `None` like `load_fold_engine`.* For folding, `None` is fine —
  a rule reports its objective `unavailable`. For screening the semantics differ: we
  always screen, and a screen that could not run must still produce a `not_run` verdict a
  report renders. `NullScreen` makes that explicit and un-forgettable; `None` would push
  the "what do we say when it didn't run" decision onto every caller, which is exactly how
  a silent `clear` creeps in.

- *Hard-coding commec's full JSON schema.* commec's exact output shape is a
  version-dependent detail, and reasoning about it is reasoning about commec internals —
  out of this lane. `_parse_commec_output` reads only the two fields BT5 needs
  (recommendation, database version) via a small key search and treats any shape it does
  not recognise as `not_run`, not `clear`.

- *Encoding `strict_biosecurity` policy in `guard_emission`.* `DesignContext` carries
  `strict_biosecurity: bool = True`, and item 4 of the brief ties it to the block path.
  But whether a `not_run` should ALSO fail-closed under strict is an **emit-policy**
  decision, and emission happens in `design/runner.py` — **S1's lane, which this session
  may not edit**. The walking skeleton currently emits with `not_run` under
  `strict_biosecurity=True`, so encoding "strict ⇒ not_run refuses" here would break it.
  `guard_emission` therefore enforces only the unambiguous block refusal and never weakens
  it; the strict-vs-advisory and fail-closed-on-not_run semantics are left to the runner.
  **Open question for S1/owner** (raised in the PR design note): does `strict_biosecurity`
  mean not_run fails closed, or does it downgrade a block to advisory? A knob that lets a
  block be bypassed is a posture weakening and would need owner sign-off.

- *Protein input validation (brief item 5).* Deferred. It is orthogonal to the hazard
  screen, `cassette/envelope.py` already rejects an amino acid with no codon under the
  table, and the managing session narrowed this session to the screen plumbing (items
  1–4). Adding it now would widen the diff without serving the failure this lane exists to
  prevent. A follow-up issue if v1 needs it.

**Scientific impact — not "none".** This changes what the app refuses to build: a `block`
verdict now refuses emission via `guard_emission`. But the wiring into the design path is
S1's, and in CI the screen still degrades to `not_run` (commec absent), so no sequence
this session can produce today is newly refused — the mechanism exists, the enforcement
point is S1's runner. Because it touches what the app REFUSES to build, the owner merges
under `CLAUDE.md` §7b; this session does not self-merge.

**Where:** branch `claude/s2-biosecurity-screen`; lane M8, `cassette/screen.py`,
`cassette/__init__.py`, `tests/cassette/test_screen.py` only.
