---
name: text-io-encoding
description: PR #102/bbffdb8 pattern — encoding="utf-8" additions and the AST-scan guard that pins them
metadata:
  type: project
---

**bbffdb8 (issue #102) added `encoding="utf-8"` to 35 file-I/O sites + 7 subprocess
text-mode sites, plus `tests/data_integrity/test_text_io_declares_encoding.py`.** Reviewed
clean: every site checked was a plain `read_text`/`write_text`/`subprocess.run(text=True)`
call with no mode change, `score/order.py`'s `newline=""` CSV writer kept its kwarg, and
`text=True` combined with `encoding="utf-8"` is not a new crash path on Linux CI (locale
there already resolves to UTF-8, and default `errors` stays `strict` either way — the
`-X warn_default_encoding` flag the PR cites flags the *omission itself*, not an actual
decode-behavior difference).

**If `tests/data_integrity/test_text_io_declares_encoding.py` is touched again**, check
these are still true (weakening any one is a silent hole, not a visible test loosening):
`NON_TEXT_OPEN_OWNERS` still excludes `os`/`zipfile`/`gzip`/etc so a real binary `.open()`
isn't force-fed `encoding=`; `_opens_in_binary_mode` still checks literal `mode=`/`"...b..."`
before flagging; `_forwards_kwargs` still exempts `**kwargs` passthroughs (a real gap the
file documents rather than hides); `MIN_FILES_SCANNED` floor still exists so a broken
`git ls-files -z` can't report a clean scan of zero files.

**This diff touched files under FOUR different §2 protected-path prefixes in one PR**
(`.github/scripts/**`, `data/codon_usage/**`, `tests/contract/**`, `tests/data_integrity/**)
— each needs its own `approved:*` label (`ci-change`, `data-change`, `contract-change`,
`oracle-change` respectively) even though issue #102 pre-authorized the cross-lane scope.
The issue authorizes the *lane* crossing (§1); it does not substitute for the per-path
labels (§2), which are separate sign-off. Check both independently on any PR this broad.
