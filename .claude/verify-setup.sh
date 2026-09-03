#!/usr/bin/env bash
# Verify the routing configuration is actually wired up.
#
# The failures this catches are all SILENT ones: an agent file missing `description`
# is skipped with no error, a `paths:` glob matching nothing means the rules file
# never loads, and a protected-path hook that matches nothing exits 0 and looks
# installed. None of these announce themselves.
set -u -o pipefail
cd "$(dirname "$0")/.." || exit 1
FAIL=0
ok()   { printf '  ok    %s\n' "$1"; }
bad()  { printf '  FAIL  %s\n' "$1"; FAIL=1; }

# PREFLIGHT. Everything below probes a hook by capturing its stdout, and for seven of
# those checks EMPTY stdout was the pass condition. A dead interpreter produces empty
# stdout, so this script used to report those seven as `ok` while the hooks were not
# running at all -- it was blind to the exact defect it exists to catch. Resolve the
# interpreter through the same launcher the hooks use, and refuse to run at all if
# there is none, because every result after that point would be vacuous.
HOOK_PY="bash .claude/hooks/run_py.sh"
PY=$($HOOK_PY --which) || {
  echo "  FAIL  no working Python 3. Every check below would be VACUOUS, so this stops here."
  echo "        Every python-launched hook in settings.json is a SILENT NO-OP, and"
  echo "        CLAUDE.md section 2 protected paths are UNENFORCED on this machine."
  echo; echo "SOME CHECKS FAILED"; exit 1; }
ok "interpreter $PY ($("$PY" -V 2>/dev/null))"

# The two frontmatter sections need PyYAML, and they SKIP without it -- which is a
# second, independent fail-open: with no PyYAML every agent, skill and rule file goes
# unverified and the script still says ALL CHECKS PASSED. The hook probes below must
# keep using $HOOK_PY (never drift from what the hooks really run), but these two may
# use any interpreter that can import yaml, and the venv always can.
YAML_PY=""
for c in .venv/bin/python "$PY"; do
  if [ -n "$c" ] && "$c" -c "import yaml" >/dev/null 2>&1; then YAML_PY="$c"; break; fi
done
if [ -z "$YAML_PY" ]; then
  NAG=$(ls .claude/agents/*.md 2>/dev/null | wc -l | tr -d ' ')
  NSK=$(ls -d .claude/skills/*/ 2>/dev/null | wc -l | tr -d ' ')
  NRU=$(ls .claude/rules/*.md 2>/dev/null | wc -l | tr -d ' ')
  bad "no interpreter can import yaml -- $NAG agents, $NSK skills and $NRU rules go UNVERIFIED. Run /bootstrap."
  YAML_PY="$PY"
else
  ok "yaml interpreter $YAML_PY"
fi

# Native path: under Git Bash $PWD is /c/Users/..., which a native Windows Python
# resolves to a nonexistent C:\c\Users\... -- that would make the off-root probe below
# report a false failure. `pwd -W` is Git Bash only and fails everywhere else.
P=$(pwd -W 2>/dev/null || pwd)

echo "== frontmatter =="
"$YAML_PY" - <<'PY' || FAIL=1
import glob, sys
try:
    import yaml
except ImportError:
    print("  SKIP  PyYAML not installed"); sys.exit(0)

bad = False
total_desc = 0
for path in sorted(glob.glob(".claude/agents/*.md") + glob.glob(".claude/skills/*/SKILL.md")):
    # glob returns OS separators; on Windows that is a backslash, which made the
    # startswith() below never match and pinned the token total at 0 -- another check
    # that passed by measuring nothing.
    path = path.replace("\\", "/")
    text = open(path).read()
    if not text.startswith("---\n"):
        print(f"  FAIL  {path}: no frontmatter (must start with ---)"); bad = True; continue
    fm = yaml.safe_load(text.split("---\n", 2)[1])
    missing = [k for k in ("name", "description") if not (fm or {}).get(k)]
    if missing:
        print(f"  FAIL  {path}: missing {missing} -- file is silently skipped"); bad = True; continue
    eff = fm.get("effort")
    if eff is not None and eff not in ("low", "medium", "high", "xhigh", "max") and not isinstance(eff, int):
        print(f"  FAIL  {path}: effort={eff!r} is not a valid value"); bad = True; continue
    if path.startswith(".claude/agents/"):
        total_desc += len(fm["description"])
    print(f"  ok    {path}  ({fm.get('model','inherit')}/{fm.get('effort','inherit')})")
tokens = total_desc // 4
if tokens > 15000:
    print(f"  FAIL  agent descriptions ~{tokens} tokens (over the 15000 threshold)"); bad = True
else:
    print(f"  ok    agent descriptions ~{tokens} tokens (warning threshold 15000)")
sys.exit(1 if bad else 0)
PY

echo "== rules globs match real files =="
"$YAML_PY" - <<'PY' || FAIL=1
import glob, subprocess, sys
try:
    import yaml
except ImportError:
    print("  SKIP  PyYAML not installed"); sys.exit(0)
bad = False
for path in sorted(glob.glob(".claude/rules/*.md")):
    text = open(path).read()
    fm = yaml.safe_load(text.split("---\n", 2)[1]) if text.startswith("---\n") else {}
    for pattern in (fm or {}).get("paths", []):
        n = len(subprocess.run(["git", "ls-files", pattern], capture_output=True, text=True).stdout.split())
        if n == 0:
            print(f"  FAIL  {path}: glob {pattern!r} matches 0 files -- rule never loads"); bad = True
        else:
            print(f"  ok    {path}: {pattern} -> {n} files")
sys.exit(1 if bad else 0)
PY

echo "== hooks are executable =="
for f in .claude/hooks/*; do
  [ -x "$f" ] && ok "$f" || bad "$f not executable"
done
[ -x scripts/gates.sh ] && ok "scripts/gates.sh" || bad "scripts/gates.sh not executable"

echo "== negative control (the verifier must be able to FAIL) =="
# Without this, every probe below can pass vacuously. BT5_PYTHON is authoritative and
# does not fall back, which is what makes this a one-liner rather than a PATH stunt.
PAY=$(printf '{"tool_name":"Write","tool_input":{"file_path":"%s/pyproject.toml"},"cwd":"%s"}' "$P" "$P")
NRC=0
printf '%s' "$PAY" | BT5_PYTHON=/nonexistent/python \
  $HOOK_PY --required .claude/hooks/protect_paths.py >/dev/null 2>&1 || NRC=$?
[ "$NRC" -eq 2 ] && ok "a dead interpreter BLOCKS (rc=2)" \
  || bad "rc=$NRC with a dead interpreter -- this verifier is blind to the original bug"

echo "== hook behaviour =="
probe_bash() { printf '{"tool_name":"Bash","tool_input":{"command":%s}}' "$1" | $HOOK_PY --optional .claude/hooks/compact_output.py; }
[ -n "$(probe_bash '"pytest tests/invariants"')" ] && ok "rewrites a bare pytest" || bad "did NOT rewrite a bare pytest"
for c in '"pytest -q && mypy"' '"mypy"' '"pytest -v"' '"pytest --maxfail=1"' '"pytest tests/contract"' '"python tests/contract/regenerate.py"'; do
  # rc AND emptiness, separately. `[ -z ]` alone made a dead interpreter the pass case.
  out=$(probe_bash "$c"); rc=$?
  if [ "$rc" -eq 0 ] && [ -z "$out" ]; then ok "stands down on $c"
  else bad "stands down on $c -- rc=$rc out=${out:0:60}"; fi
done
D=$(printf '{"tool_name":"Edit","tool_input":{"file_path":"%s/packages/engine/src/bt5/core/types.py"},"cwd":"%s"}' "$P" "$P" | $HOOK_PY --required .claude/hooks/protect_paths.py)
case "$D" in *'"ask"'*) ok "asks on an ABSOLUTE protected path" ;; *) bad "protected-path hook matched nothing (the silent no-op)" ;; esac
D2=$(printf '%s' "$PAY" | $HOOK_PY --required .claude/hooks/protect_paths.py)
case "$D2" in *'"ask"'*) ok "asks on pyproject.toml (the pair CI does not cover)" ;; *) bad "no decision for pyproject.toml" ;; esac
# The OFF-ROOT probe: the only one here that fails on a pre-repo_root_for() hook. Every
# probe that uses the same value for cwd and file_path tests the one configuration that
# always worked, and would attest "the guard answers" on a machine where it is blind.
W=$(printf '{"tool_name":"Write","tool_input":{"file_path":"%s/pyproject.toml"},"cwd":"%s"}' "$P" "$(dirname "$P")" | $HOOK_PY --required .claude/hooks/protect_paths.py)
case "$W" in *'"ask"'*) ok "asks when the session cwd is not the repo root" ;;
  *) bad "SILENT PASS with an off-root cwd -- a main-repo session editing a worktree file writes to protected paths unprompted. See repo_root_for() in protect_paths.py." ;; esac
B=$(printf '{"tool_name":"Edit","tool_input":{"file_path":"%s/README.md"},"cwd":"%s"}' "$P" "$P" | $HOOK_PY --required .claude/hooks/protect_paths.py); brc=$?
if [ "$brc" -eq 0 ] && [ -z "$B" ]; then ok "allows a benign path"; else bad "benign path: rc=$brc out=${B:0:60}"; fi
# Assert the ARTIFACT, not just the exit status: statusline.py returns 0 without
# printing on unparseable JSON, so `>/dev/null && ok` could pass while nothing rendered.
S=$(echo '{"model":{"display_name":"Opus"},"effort":{"level":"high"}}' | $HOOK_PY --statusline .claude/hooks/statusline.py)
case "$S" in *Opus*) ok "statusline renders" ;; *) bad "statusline printed: ${S:-<empty>}" ;; esac
echo garbage | $HOOK_PY --optional .claude/hooks/compact_output.py >/dev/null 2>&1 && ok "hooks fail open on garbage" || bad "hook errored on garbage input"

echo "== settings.json =="
[ -s .claude/settings.json ] || bad "settings.json is missing or empty"
"$PY" -c "import json;json.load(open('.claude/settings.json'))" && ok "valid JSON" || bad "invalid JSON"
grep -q '"Edit(.claude/\*\*)"' .claude/settings.json && ok "Edit(.claude/**) present" || bad "Edit(.claude/**) missing -- hooks/skills become uneditable"
grep -qE '"command": *"(python|python3|py|node|npx|sh|ruby|perl) ' .claude/settings.json && bad "a hook invokes a bare interpreter name (see .claude/hooks/run_py.sh)" || ok "every hook goes through run_py.sh"
grep -q '"Write(.claude/\*\*)"' .claude/settings.json && bad "Write(path) rule present (never consulted; warns at startup)" || ok "no Write(path) rule"
grep -q '"ultracode"' .claude/settings.json && bad "session-wide ultracode set (every turn would fan out)" || ok "ultracode not set session-wide"

echo "== no invalid effort value in frontmatter =="
# Only real frontmatter lines: an `effort:` at the start of a line in an agent or
# skill file. Prose mentioning the invalid value (this script included) is fine.
if [ ! -d .claude/agents ] || [ ! -d .claude/skills ]; then
  bad "'.claude/agents' or '.claude/skills' is missing -- the scan below would pass vacuously"
elif grep -rn --include='*.md' -E '^effort: *ultracode' .claude/agents .claude/skills >/dev/null 2>&1; then
  bad "an agent or skill declares the ultracode effort, which the schema rejects"
else ok "no agent or skill declares an invalid effort"; fi

echo "== CLAUDE.md =="
L=$(wc -l < CLAUDE.md)
[ "$L" -lt 200 ] && ok "CLAUDE.md $L lines (< 200)" || bad "CLAUDE.md $L lines (>= 200)"
N=$(sed -n '/## 3\./,/## 4\./p' CLAUDE.md | grep -cE '^[0-9]\. \*\*')
[ "$N" -eq 7 ] && ok "all 7 correctness statements present" || bad "only $N of 7 correctness statements"
for s in "## Stack" "## Delegation" "## Compact instructions"; do
  grep -q "^$s" CLAUDE.md && ok "$s present" || bad "$s missing"
done

echo
[ "$FAIL" -eq 0 ] && echo "ALL CHECKS PASSED" || echo "SOME CHECKS FAILED"
exit "$FAIL"
