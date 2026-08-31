#!/usr/bin/env bash
# Verify the routing configuration is actually wired up.
#
# The failures this catches are all SILENT ones: an agent file missing `description`
# is skipped with no error, a `paths:` glob matching nothing means the rules file
# never loads, and a protected-path hook that matches nothing exits 0 and looks
# installed. None of these announce themselves.
set -u
cd "$(dirname "$0")/.." || exit 1
FAIL=0
ok()   { printf '  ok    %s\n' "$1"; }
bad()  { printf '  FAIL  %s\n' "$1"; FAIL=1; }

echo "== frontmatter =="
python3 - <<'PY' || FAIL=1
import glob, sys
try:
    import yaml
except ImportError:
    print("  SKIP  PyYAML not installed"); sys.exit(0)

bad = False
total_desc = 0
for path in sorted(glob.glob(".claude/agents/*.md") + glob.glob(".claude/skills/*/SKILL.md")):
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
print(f"  ok    agent descriptions ~{total_desc // 4} tokens (warning threshold 15000)")
sys.exit(1 if bad else 0)
PY

echo "== rules globs match real files =="
python3 - <<'PY' || FAIL=1
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

echo "== hook behaviour =="
probe_bash() { printf '{"tool_name":"Bash","tool_input":{"command":%s}}' "$1" | python3 .claude/hooks/compact_output.py; }
[ -n "$(probe_bash '"pytest tests/invariants"')" ] && ok "rewrites a bare pytest" || bad "did NOT rewrite a bare pytest"
for c in '"pytest -q && mypy"' '"mypy"' '"pytest -v"' '"pytest --maxfail=1"' '"pytest tests/contract"' '"python tests/contract/regenerate.py"'; do
  [ -z "$(probe_bash "$c")" ] && ok "stands down on $c" || bad "rewrote $c and must not"
done
P="$PWD"
D=$(printf '{"tool_name":"Edit","tool_input":{"file_path":"%s/packages/engine/src/bt5/core/types.py"},"cwd":"%s"}' "$P" "$P" | python3 .claude/hooks/protect_paths.py)
case "$D" in *'"ask"'*) ok "asks on an ABSOLUTE protected path" ;; *) bad "protected-path hook matched nothing (the silent no-op)" ;; esac
B=$(printf '{"tool_name":"Edit","tool_input":{"file_path":"%s/README.md"},"cwd":"%s"}' "$P" "$P" | python3 .claude/hooks/protect_paths.py)
[ -z "$B" ] && ok "allows a benign path" || bad "protected-path hook fired on README.md"
echo '{"model":{"display_name":"Opus"},"effort":{"level":"high"}}' | python3 .claude/hooks/statusline.py >/dev/null && ok "statusline renders" || bad "statusline failed"
echo garbage | python3 .claude/hooks/compact_output.py >/dev/null 2>&1 && ok "hooks fail open on garbage" || bad "hook errored on garbage input"

echo "== settings.json =="
python3 -c "import json;json.load(open('.claude/settings.json'))" && ok "valid JSON" || bad "invalid JSON"
grep -q '"Write(.claude/\*\*)"' .claude/settings.json && bad "Write(path) rule present (never consulted; warns at startup)" || ok "no Write(path) rule"
grep -q '"ultracode"' .claude/settings.json && bad "session-wide ultracode set (every turn would fan out)" || ok "ultracode not set session-wide"

echo "== no invalid effort value in frontmatter =="
# Only real frontmatter lines: an `effort:` at the start of a line in an agent or
# skill file. Prose mentioning the invalid value (this script included) is fine.
if grep -rn --include='*.md' -E '^effort: *ultracode' .claude/agents .claude/skills >/dev/null 2>&1; then
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
