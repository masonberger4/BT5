#!/usr/bin/env bash
# Launch a .claude/hooks/*.py hook under a Python that is PROVED to work, or fail in a
# way somebody sees.
#
# THE FAILURE THIS EXISTS TO AVOID. Every hook here used to be registered as
# `python3 <script>`. On Windows `python3` resolves to the Microsoft Store stub: it
# exits 49 and writes NOTHING to stdout. Claude Code blocks a PreToolUse hook only on
# exit 2, so 49 is a "non-blocking error" -- rendered to the human as one red line,
# never entering the model's context, tool call proceeds. Hundreds of such records
# accumulated in this repo's session logs while protect_paths.py -- the only mechanical
# enforcement of CLAUDE.md section 2 for the `Write` tool and for pyproject.toml /
# uv.lock, which .github/scripts/check-approval-labels.sh says verbatim it does
# "deliberately NOT" enforce -- was completely off.
#
# EXISTENCE CHECKS CANNOT DETECT IT. `command -v python3` SUCCEEDS on the stub; it is a
# real executable file. So does `which`, `type -p` and `[ -x ]`. The only test that
# separates a live interpreter from a corpse is RUNNING it and reading its OUTPUT.
# `-V` is the probe: on Python 3 it prints "Python 3.x.y" on STDOUT; the stub prints
# nothing; Python 2 prints to STDERR. One cheap spawn rejects both wrong answers.
#
# NO NAME IS CORRECT EVERYWHERE. python3 is dead on Windows and is the only name
# guaranteed on the Linux web container. `python` is absent from many Linux images, and
# on Windows it works only by PATH ordering (a second, dead WindowsApps `python` stub
# sits one entry behind the good one). `py` exists only on Windows. `.venv/bin/python`
# is absent from a fresh checkout and from the claude-review CI job, which runs no
# bootstrap, and `uv venv` never creates `bin/python3` on Windows at all. Hence: probe
# at hook-run time, ordered for speed, correct regardless of order.
#
# TIERS. The tier is a property of the REGISTRATION SITE, not of the script.
#   --required    no interpreter -> exit 2. Claude Code BLOCKS the tool call and puts
#                 the message in front of the MODEL. Only for guardrails whose loss is
#                 covered by nothing else: protect_paths.py on Edit|Write, and
#                 push_gate.py on the MCP pull-request tools.
#   --optional    no interpreter -> exit 0 + stderr. For hooks whose loss costs
#                 convenience, not a contract. This is what keeps `Bash` ALIVE on a
#                 Python-less box, so the machine stays repairable from inside a session.
#   --statusline  no interpreter -> a banner on STDOUT, which IS the status line.
#                 statusLine failures are otherwise the most silent path in the whole
#                 system: debug log plus once-per-session telemetry, no UI at all.
#   --which       print the resolved interpreter, or exit 1. The single resolver shared
#                 by ruff_format.sh, verify-setup.sh and apply-repo-settings.sh, so the
#                 verifier cannot drift into testing a different chain than the hooks use.
#
# NO OFF SWITCH. There is no BT5_HOOKS_DISABLE and no BT5_HOOKS_FAIL_OPEN. If this picks
# the wrong interpreter, point BT5_PYTHON at the right one -- that is a fix. A committed
# lever that restores the silent no-op is CLAUDE.md section 4 by the back door.
#
# BT5_PYTHON IS AUTHORITATIVE. If it is set, it is the ONLY candidate. It does not fall
# back, because silently running a different interpreter than the one an operator named
# is the same class of bug this file exists to remove -- and because that makes the
# negative control in verify-setup.sh a one-liner. Bash stays alive under --optional, and
# session_start.sh reports the override on every session, so a stale value is loud and
# repairable, never a brick.
set -u

usage() {
  printf 'usage: run_py.sh --required|--optional|--statusline <script.py> [args...]\n' >&2
  printf '       run_py.sh --which\n' >&2
  exit 2
}

MODE="${1-}"
[ -n "$MODE" ] || usage
shift
SCRIPT=""
case "$MODE" in
  --which) ;;
  --required|--optional|--statusline)
      SCRIPT="${1-}"
      [ -n "$SCRIPT" ] || usage
      shift ;;
  *) usage ;;
esac

# A candidate is usable only if it RUNS and PRINTS "Python 3." on stdout.
# `command -v` / `[ -x ]` are used ONLY as a free NEGATIVE filter -- they cannot be
# trusted positively (see the header), and nothing here does.
usable() {
  case "$1" in
    */*) [ -x "$1" ] || return 1 ;;
    *)   command -v "$1" >/dev/null 2>&1 || return 1 ;;
  esac
  case "$("$1" -V 2>/dev/null </dev/null)" in
    "Python 3."*) return 0 ;;
    *)            return 1 ;;
  esac
}

# Ordering is a PERFORMANCE choice only; the probe makes it safe either way. Windows
# first-hits `python` and never spawns the ~280ms dead stub; elsewhere `python3` first so
# a stale `python` can never shadow a live `python3`. $OS is set by Windows itself in
# every process environment; $OSTYPE is a bash build detail (it reads "cygwin" on this
# machine's Git Bash, not "msys") and is only a backstop.
if [ "${OS-}" = "Windows_NT" ]; then
  CHAIN="python python3 py .venv/bin/python"
else
  case "${OSTYPE-}" in
    msys*|cygwin*|win*) CHAIN="python python3 py .venv/bin/python" ;;
    *)                  CHAIN="python3 python .venv/bin/python" ;;
  esac
fi

PY=""
BAD_OVERRIDE=""
if [ -n "${BT5_PYTHON-}" ]; then
  if usable "$BT5_PYTHON"; then PY="$BT5_PYTHON"; else BAD_OVERRIDE="$BT5_PYTHON"; fi
else
  # Unquoted ON PURPOSE: word splitting is wanted, and every token is space-free.
  # BT5_PYTHON is handled above precisely so a path with spaces never reaches here.
  for c in $CHAIN; do
    if usable "$c"; then PY="$c"; break; fi
  done
fi

if [ "$MODE" = "--which" ]; then
  [ -n "$PY" ] || exit 1
  printf '%s\n' "$PY"
  exit 0
fi

# exec: byte-exact passthrough of stdin, stdout, stderr and exit status. No retry loop,
# so a hook is never executed twice and its real stdout is never discarded. A hook that
# runs and crashes exits non-zero WITH its traceback, correctly attributed to the hook
# instead of being misreported as a dead interpreter.
if [ -n "$PY" ] && [ -f "$SCRIPT" ]; then
  exec "$PY" "$SCRIPT" "$@"
fi

# Failure paths only. Drain the payload first: exiting without reading stdin sends
# SIGPIPE to whoever is writing the tool_input, and this fires exactly on a large Write.
[ -t 0 ] || cat >/dev/null 2>&1

if [ -n "$BAD_OVERRIDE" ]; then
  MSG="BT5 hooks: BT5_PYTHON is set to '$BAD_OVERRIDE', which did not print 'Python 3.x'
  for -V. That override is authoritative: run_py.sh does NOT quietly fall back to the
  normal chain. Fix or unset BT5_PYTHON in .claude/settings.local.json -> \"env\"."
elif [ -z "$PY" ]; then
  MSG="BT5 hooks: no working Python 3 (tried: $CHAIN -- each RUN, not merely looked up).
  '$SCRIPT' did NOT run. On Windows 'python3' is normally a Microsoft Store stub that
  exits 49 with empty stdout, which Claude Code treats as a NON-blocking hook failure --
  i.e. silently. Install Python 3.11, or set BT5_PYTHON in
  .claude/settings.local.json -> \"env\" to a working interpreter."
else
  MSG="BT5 hooks: '$SCRIPT' does not exist (cwd: $PWD), so the hook did NOT run.
  .claude/settings.json references a script that is not on disk. Fix the path or the
  settings entry; .github/scripts/check-hook-commands.py fails CI on exactly this."
fi

case "$MODE" in
  --statusline)
    printf 'BT5  HOOKS DEAD - guardrails OFF (bash .claude/verify-setup.sh)\n'
    printf '%s\n' "$MSG" >&2
    exit 0 ;;
  --required)
    printf '%s\n' "$MSG" >&2
    printf 'This tool call is BLOCKED, because a guardrail that cannot run must not be\n' >&2
    printf 'indistinguishable from a guardrail that said yes. Bash and Read still work,\n' >&2
    printf 'so this is repairable from inside the session.\n' >&2
    exit 2 ;;
  *)
    printf '%s\n' "$MSG" >&2
    exit 0 ;;
esac
