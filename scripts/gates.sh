#!/usr/bin/env bash
# The local gate chain, in the order CI runs it.
#
# Every gate runs INDEPENDENTLY. They are deliberately not chained with `&&`:
# chaining stops at the first failure and hides the rest, so a caller fixes one
# thing, re-runs, and discovers the next -- one CI cycle per gate. Here you get
# the whole picture in one pass.
#
# Each gate prints `GATE <name> EXIT <n>`. Exit codes are NOT all equivalent:
#   0  pass
#   1  fail (tests failed, lint findings)
#   2  BROKEN -- interrupted or collection error, i.e. the environment
#   3  BROKEN -- internal error
#   4  BROKEN -- usage error (bad args or ini)
#   5  BROKEN -- no tests collected; never read this as success
# This script's own exit is 10 when there is no usable venv, 1 if any gate was
# non-zero, 0 only when every gate passed.
#
# mypy is included and is LOCAL-ONLY: no CI job runs it, and it is the only
# thing that checks kmers.py's KmerIndex conformance assertion.
set -u

VENV_PY=".venv/bin/python"
FAILED=0
export PYTHONHASHSEED=0

if [ ! -x "$VENV_PY" ]; then
  echo "GATE environment EXIT 10"
  echo "No .venv found. Run /bootstrap:"
  echo "  uv venv --python 3.11 .venv"
  echo '  uv pip install --python .venv/bin/python -e ".[dev,fold,export]"'
  echo "Bare pytest/ruff/mypy resolve to an interpreter without numpy and will mislead you."
  exit 10
fi

if ! "$VENV_PY" -c "import numpy, Bio, hypothesis, bt5" 2>/dev/null; then
  echo "GATE environment EXIT 10"
  echo "The venv exists but is incomplete. Missing one of: numpy, Bio, hypothesis, bt5."
  echo 'Run /bootstrap: uv pip install --python .venv/bin/python -e ".[dev,fold,export]"'
  exit 10
fi
echo "GATE environment EXIT 0"

run_gate() {
  # run_gate <name> <command...>
  local name="$1"; shift
  echo "--- $name"
  "$@"
  local rc=$?
  echo "GATE $name EXIT $rc"
  [ "$rc" -ne 0 ] && FAILED=1
  return 0
}

run_gate ruff-check      .venv/bin/ruff check . --output-format=concise
run_gate ruff-format     .venv/bin/ruff format --check .
run_gate mypy            .venv/bin/mypy
run_gate invariants      .venv/bin/pytest tests/invariants -q -p no:randomly --tb=short
run_gate data-integrity  .venv/bin/pytest tests/data_integrity -q --tb=short
run_gate contract        .venv/bin/pytest tests/contract -q --tb=short
run_gate engine-tests    .venv/bin/pytest packages/engine/tests -q -m "not slow" --tb=short --maxfail=10

echo
if [ "$FAILED" -eq 0 ]; then
  echo "ALL GATES PASSED"
  echo "NOTE mypy is local-only; no CI job runs it."
  exit 0
fi
echo "ONE OR MORE GATES FAILED -- see the GATE ... EXIT lines above."
echo "NOTE exit 2/3/4/5 mean BROKEN (environment or usage), not a test failure."
exit 1
