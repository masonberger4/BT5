#!/usr/bin/env bash
# Protected paths require their approval label.
#
# CLAUDE.md section 2 lists files that must never change without a matching
# `approved:*` label. Until now that was prose: no CI job read the labels, and
# the labels did not exist in the repository at all, so the gate could neither
# pass nor fail. This makes it mechanical.
#
# In a script rather than inline YAML so it can be run against fixtures. A gate
# that has never been executed against a case it should REJECT is a gate nobody
# has tested, and this one exists precisely to reject things.
#
# Usage:  check-approval-labels.sh <changed-files-file> <labels-json>
# Exits 0 when every protected path touched carries its label, 1 otherwise.
set -euo pipefail

CHANGED_FILE="${1:?usage: $0 <changed-files-file> <labels-json>}"
LABELS_JSON="${2:-[]}"
CHANGED="$(cat "$CHANGED_FILE")"
missing=0

has_label() {
  # Match "name" as a complete JSON string so `approved:ci-change` cannot be
  # satisfied by a label that merely contains it as a prefix.
  printf '%s' "$LABELS_JSON" | grep -qF "\"$1\""
}

require() {
  local label="$1"
  shift
  local hit=""
  for pattern in "$@"; do
    hit="$(printf '%s\n' "$CHANGED" | grep -E "$pattern" || true)"
    [ -n "$hit" ] && break
  done
  [ -z "$hit" ] && return 0
  local first
  first="$(printf '%s\n' "$hit" | head -1)"
  if has_label "$label"; then
    printf '  ok       %-28s (%s)\n' "$label" "$first"
  else
    printf '  MISSING  %-28s (%s)\n' "$label" "$first"
    echo "::error::$label is required: this pull request changes $first"
    missing=1
  fi
}

echo "Protected-path approval check"
require "approved:oracle-change" \
  '^packages/engine/src/bt5/verify\.py$' \
  '^tests/(invariants|data_integrity)/'
require "approved:contract-change" '^packages/engine/src/bt5/core/'
require "approved:algorithm-change" '^benchmarks/(baseline\.json|tolerances\.yaml)$'
require "approved:data-change" '^data/(genetic_codes|codon_usage)/'
require "approved:ci-change" '^\.github/'

if [ "$missing" -ne 0 ]; then
  cat <<'EOF'

One or more protected paths changed without their approval label. Add the label
named above to this pull request; the check re-runs on `labeled`.

If the label does not exist in the repository yet it must be created once:
  approved:oracle-change  approved:contract-change  approved:algorithm-change
  approved:data-change    approved:ci-change

pyproject.toml and uv.lock are protected by CLAUDE.md section 2 but no label is
named for them, so they are deliberately NOT enforced here.
EOF
  exit 1
fi
echo "  all protected paths carry their label (or none were touched)"
