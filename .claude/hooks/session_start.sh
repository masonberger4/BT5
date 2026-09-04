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

# Guardrail liveness. Two probes, because the two failure classes are different.
# The SHAPE grep catches the config regressing on a box where the interpreter is fine
# (i.e. it works on Linux, where the python3-stub bug cannot reproduce). The END-TO-END
# probe runs the real hook through the real launcher, so it survives a rename of either.
#
# `pwd -W` is not cosmetic. Under Git Bash $PWD is an MSYS path (/c/Users/...), which a
# native Windows Python resolves to the nonexistent C:\c\Users\... -- the walk up to
# `.git` then finds nothing and probe B reports a false PARTIAL. `pwd -W` yields the
# native path; it fails on Linux/macOS, where plain `pwd` is already correct.
if grep -qE '"command": *"(python|python3|py|node|npx|sh|ruby|perl) ' .claude/settings.json; then
  echo "  hooks    MISROUTED -- a hook in settings.json invokes a bare interpreter name."
  echo "           See .claude/hooks/run_py.sh: that is how this went silently dead before."
fi
_P=$(pwd -W 2>/dev/null || pwd)
_gp() {  # $1 = the cwd to claim; the guard must answer regardless of it
  printf '{"tool_name":"Write","tool_input":{"file_path":"%s/pyproject.toml"},"cwd":"%s"}' \
    "$_P" "$1" | bash .claude/hooks/run_py.sh --required .claude/hooks/protect_paths.py 2>/dev/null
}
_A=$(_gp "$_P"); _B=$(_gp "$(dirname "$_P")")
case "$_A$_B" in
  *'"ask"'*'"ask"'*)
    echo "  hooks    ok (protect_paths answers; python = $(bash .claude/hooks/run_py.sh --which))" ;;
  *'"ask"'*)
    echo "  hooks    PARTIAL -- the guard answers for this cwd but NOT for an off-root cwd."
    echo "           A session started outside this worktree could write to protected"
    echo "           paths unprompted. See repo_root_for() in protect_paths.py." ;;
  *)
    echo "  hooks    DEAD -- protect_paths.py returned no decision for pyproject.toml."
    echo "           It is the ONLY thing stopping a Write to verify.py, core/**,"
    echo "           tests/contract/**, data/**, .github/**, pyproject.toml and uv.lock."
    echo "           Treat CLAUDE.md section 2 as UNENFORCED; bash .claude/verify-setup.sh" ;;
esac
[ -n "${BT5_PYTHON-}" ] && echo "  hooks    BT5_PYTHON override active: $BT5_PYTHON"

BRANCH=$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "?")
DIRTY=$(git status --porcelain 2>/dev/null | wc -l | tr -d ' ')
echo "  branch   $BRANCH ($DIRTY uncommitted)"
echo "  gates    bash scripts/gates.sh   (or /pre-pr for the full chain)"
echo "  note     Hypothesis runs the dev profile (50 examples) locally; CI runs ci (200)."
echo "           A property passing here has had a quarter of CI's search."
exit 0
