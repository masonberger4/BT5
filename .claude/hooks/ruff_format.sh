#!/usr/bin/env bash
# PostToolUse(Edit|Write): format the edited Python file so the `ruff format --check .`
# CI gate stops failing on whitespace.
#
# Uses .venv/bin/ruff ONLY. The ruff on PATH is a different version from the pinned
# `ruff>=0.4` in the venv, and formatting with the wrong version produces a diff CI
# then rejects -- which is the exact failure this hook exists to prevent.
#
# Silent no-op when there is no venv (a fresh checkout has none), when the file is not
# Python, or when it is outside the source and test trees. Never touches .claude/ or docs/.
set -u

INPUT="$(cat)"
RUFF=".venv/bin/ruff"
[ -x "$RUFF" ] || exit 0

FILE=$(printf '%s' "$INPUT" | python3 -c '
import json,sys
try: d=json.load(sys.stdin)
except Exception: sys.exit(0)
ti=d.get("tool_input") or {}
p=ti.get("file_path") if isinstance(ti,dict) else None
print(p or "")
' 2>/dev/null) || exit 0

[ -n "$FILE" ] || exit 0
case "$FILE" in *.py) ;; *) exit 0 ;; esac
[ -f "$FILE" ] || exit 0

REL="${FILE#"$PWD"/}"
case "$REL" in
  packages/engine/src/*|packages/engine/tests/*|tests/*) ;;
  *) exit 0 ;;
esac

"$RUFF" format "$FILE" >/dev/null 2>&1 || true
exit 0
