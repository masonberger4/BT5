## 2026-09-04 — every text I/O site names UTF-8, and a data_integrity test stands in for PLW1514

**Decided:** add `encoding="utf-8"` to all 35 file-I/O call sites #102 enumerated, plus 7
`subprocess` text-mode sites the issue's enumeration did not count, and pin the property
with `tests/data_integrity/test_text_io_declares_encoding.py` — an AST scan of every
git-tracked `.py` file, in two tests (file I/O, subprocess decode).

**Rejected:**

- *Enable ruff's `PLW1514` instead.* It is the right long-term mechanism and is one entry
  in `[tool.ruff.lint] select`. That entry lives in `pyproject.toml`, which CLAUDE.md §2
  protects and for which **no `approved:*` label exists** — `protect_paths.py` says
  "Open an issue instead." The test is the stand-in until the owner enables the rule; its
  docstring says to delete the file rather than maintain both. `PLW1514` would still not
  cover the subprocess half, which has no ruff equivalent.
- *A new `.ruff.toml` to add the rule without touching `pyproject.toml`.* Ruff reads one
  config file, so a root `.ruff.toml` would **replace** the `[tool.ruff]` block entirely
  and silently drop `line-length`, `target-version`, the 13-item `select` and the
  per-file-ignores. A worse outcome than the omission it fixes.
- *A new top-level `tests/policy/` directory.* CI runs exactly four pytest invocations
  (`ci.yml:157,169,194,214`): `tests/invariants`, `tests/data_integrity`, `tests/contract`,
  `packages/engine/tests`. A fifth directory would never run — a gate that looks installed
  and gates nothing, which is what `.claude/rules/tests.md` catalogues. `data_integrity` is
  the right home anyway: it already holds the repo-wide source scans
  (`test_no_expression_claims`, `test_oracle_independence`).
- *`-X warn_default_encoding` as the gate itself.* It is the authority and was used as the
  verification (below), but it only reports lines that **execute**. Most of these sites are
  error paths, so as a gate it would pass a tree full of omissions. Turning it on also
  needs `pyproject.toml`.
- *Stopping at the 35 sites #102 counted.* `subprocess.run(..., text=True)` resolves
  `encoding=None` through `locale.getencoding()` in `Popen.__init__` for exactly the reason
  `open()` does. Leaving it would have left `.claude/hooks/push_gate.py` and
  `statusline.py` decoding git output as cp1252 on Windows, and would have made the new
  guard file an offender against its own property.
- *`errors="replace"` on the user-facing reads (`cli.py`, `vector/io.py`).* Strict is the
  point: a FASTA that is not UTF-8 should fail identically on every platform, not decode
  differently per machine. Behaviour on Linux is unchanged either way.

**Evidence:**

- 35 file sites found by AST walk, matching #102's count exactly; 0 remain after the fix.
- `python -X warn_default_encoding -W error::EncodingWarning -m pytest tests packages/engine/tests`
  → 8 failures before, all 7 subprocess sites plus the new guard's own `git ls-files` call.
  After: **2116 passed**, with the single third-party frame excluded
  (`-W ignore::EncodingWarning:openpyxl.worksheet._writer`; `openpyxl/worksheet/_writer.py:36`
  opens a text `NamedTemporaryFile` with no encoding — upstream's defect, not ours, and our
  xlsx cell values are ASCII).
- Guard probed both ways: it fires on `open`, `Path.open`, `read_text`, `write_text`,
  `NamedTemporaryFile` and subprocess text mode; it stays quiet on `mode="rb"`, `os.open`,
  `**kwargs` passthrough and calls that already name an encoding.
- `bash scripts/gates.sh` → ALL GATES PASSED.
- The three CI scripts touched still run clean: `check-workflow-gate.py`,
  `check-hook-commands.py` (probes the live guard), `test-attestation-matcher.py`.

**Where:** branch `claude/github-issues-pykbvy`, closes #102.
