#!/usr/bin/env bash
# SessionStart: report the environment. NEVER installs anything.
#
# A fresh checkout has no .venv and bt5's third-party deps are absent, so every session
# otherwise rediscovers this by running a command that fails in a way that looks like a
# code defect. Installing here instead would pull viennarna and biopython -- minutes,
# silently, on every session start -- so this only reports.
set -u
cat <<'HDR'
BT5 environment
HDR

if [ -x .venv/bin/python ]; then
  PYV=$(.venv/bin/python -c 'import sys;print("%d.%d.%d"%sys.version_info[:3])' 2>/dev/null || echo "?")
  if .venv/bin/python -c "import numpy, Bio, hypothesis, bt5" 2>/dev/null; then
    echo "  venv     ready (.venv, python $PYV)"
  else
    MISSING=$(.venv/bin/python - <<'PY' 2>/dev/null || echo "unknown"
import importlib.util as u
print(", ".join(m for m in ("numpy","Bio","hypothesis","bt5") if u.find_spec(m) is None))
PY
)
    echo "  venv     INCOMPLETE (missing: $MISSING) -- run /bootstrap"
  fi
else
  echo "  venv     MISSING -- run /bootstrap before any gate"
  echo "           Bare pytest/ruff/mypy resolve to an interpreter without numpy;"
  echo "           bare pytest exits 4 on a conftest import error, which looks like"
  echo "           a code failure and is not one."
fi

BRANCH=$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "?")
DIRTY=$(git status --porcelain 2>/dev/null | wc -l | tr -d ' ')
echo "  branch   $BRANCH ($DIRTY uncommitted)"
echo "  gates    bash scripts/gates.sh   (or /pre-pr for the full chain)"
echo "  note     Hypothesis runs the dev profile (50 examples) locally; CI runs ci (200)."
echo "           A property passing here has had a quarter of CI's search."
exit 0
