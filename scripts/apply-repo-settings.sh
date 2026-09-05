#!/usr/bin/env bash
# Apply the BT5 repository settings and the main-branch ruleset.
#
# Run this from a machine with `gh` authenticated as the repository owner.
# It cannot be run from a Claude Code web session: that session's proxy blocks
# every repo-scoped GitHub REST path, and the GitHub MCP server exposes no
# ruleset tool.
#
# ORDER MATTERS. Repo settings are applied before the ruleset, because the
# ruleset requires the squash merge method -- if the repository has squash
# disabled when the ruleset lands, EVERY merge blocks.
set -euo pipefail

REPO="${1:-masonberger4/BT5}"
RULESET_FILE="$(dirname "$0")/../.github/rulesets/main-protection.json"

command -v gh >/dev/null || { echo "error: gh CLI not found" >&2; exit 1; }
gh auth status >/dev/null 2>&1 || { echo "error: gh is not authenticated" >&2; exit 1; }

echo "==> 1/4 repository settings"
gh api -X PATCH "repos/$REPO" \
  -F allow_squash_merge=true \
  -F allow_merge_commit=false \
  -F allow_rebase_merge=false \
  -f squash_merge_commit_title=PR_TITLE \
  -f squash_merge_commit_message=PR_BODY \
  -F allow_auto_merge=true \
  -F delete_branch_on_merge=true \
  -F allow_update_branch=true \
  -F has_wiki=false \
  -F has_projects=false \
  -F has_issues=true >/dev/null
echo "    squash-only, auto-merge on, head branches auto-deleted"

echo "==> 2/4 Actions permissions"
# Read-only default token. "can_approve_pull_request_reviews=false" is
# security-critical: with it enabled, anything able to edit .github/workflows
# could write a workflow that approves its own pull request.
gh api -X PUT "repos/$REPO/actions/permissions/workflow" \
  -f default_workflow_permissions=read \
  -F can_approve_pull_request_reviews=false >/dev/null
echo "    default token read-only; Actions cannot approve PRs"

# Fork pull requests must not run workflows unreviewed. This guards
# `required-checks`, which as of 2026-09-05 is the ONLY required context.
#
# WHY. A required status check is matched by CONTEXT NAME plus integration id,
# and 15368 is the generic "GitHub Actions" app -- not this workflow, not this
# job. A `pull_request`-triggered workflow runs the file FROM THE PULL REQUEST'S
# OWN BRANCH. So a pull request that adds `.github/workflows/anything.yml` with a
# job named `required-checks` produces a second check run of that exact name,
# from that exact app, on that exact SHA -- and it can be made to succeed
# trivially. Nothing in check-workflow-gate.py sees it: that script reads the
# workflows the repo already ships, not ones a pull request introduces.
#
# GitHub's default here is `first_time_contributors`, which lets a RETURNING
# outside contributor's workflows run unreviewed. `all_external_contributors`
# requires a human to approve every fork run, which is what actually stands
# between that name collision and a forged green. Dropping `pre-pr-attest` from
# the ruleset narrowed the target but did not remove it.
gh api -X PUT "repos/$REPO/actions/permissions/fork-pr-contributor-approval" \
  -f approval_policy=all_external_contributors >/dev/null
echo "    every fork-PR workflow run needs human approval"

echo "==> 3/4 discovering the GitHub Actions app id"
# Pinning the required check to its producing app is what stops any token with
# statuses:write from posting a forged green. Discover it rather than assume:
# a WRONG id means the check never matches and the pull request blocks forever.
APP_ID="$(gh api "repos/$REPO/commits/main/check-runs" \
  --jq '[.check_runs[] | select(.name=="required-checks") | .app.id] | first' 2>/dev/null || true)"
if [ -z "$APP_ID" ] || [ "$APP_ID" = "null" ]; then
  echo "    could not find a 'required-checks' run on main; falling back to the" >&2
  echo "    documented GitHub Actions app id 15368. Verify with:" >&2
  echo "      gh api repos/$REPO/commits/main/check-runs --jq '.check_runs[] | {name, app: .app.id}'" >&2
  APP_ID=15368
fi
echo "    integration_id = $APP_ID"

echo "==> 4/4 creating the main-protection ruleset"
TMP="$(mktemp)"
trap 'rm -f "$TMP"' EXIT
PY_BIN=$(bash .claude/hooks/run_py.sh --which) || { echo "no working Python 3 -- see .claude/hooks/run_py.sh" >&2; exit 1; }
"$PY_BIN" - "$RULESET_FILE" "$APP_ID" > "$TMP" <<'PY'
import json, sys
spec = json.load(open(sys.argv[1]))
for rule in spec["rules"]:
    if rule["type"] == "required_status_checks":
        for check in rule["parameters"]["required_status_checks"]:
            check["integration_id"] = int(sys.argv[2])
json.dump(spec, sys.stdout)
PY

EXISTING="$(gh api "repos/$REPO/rulesets" --jq '[.[] | select(.name=="main-protection") | .id] | first' 2>/dev/null || true)"
if [ -n "$EXISTING" ] && [ "$EXISTING" != "null" ]; then
  echo "    updating existing ruleset $EXISTING (PUT, not PATCH)"
  gh api -X PUT "repos/$REPO/rulesets/$EXISTING" --input "$TMP" >/dev/null
else
  gh api -X POST "repos/$REPO/rulesets" --input "$TMP" >/dev/null
fi

echo
echo "==> verification"
gh api "repos/$REPO/rulesets" --jq '.[] | "  \(.name): \(.enforcement)"'
echo "  legacy branch protection (expect 404):"
gh api "repos/$REPO/branches/main/protection" >/dev/null 2>&1 \
  && echo "    WARNING: legacy protection exists and stacks most-restrictive-wins; remove it" \
  || echo "    none (good)"
echo
echo "Done. Open a throwaway PR and confirm the merge box shows: one required"
echo "check, squash-only, no approval requirement, and NO bypass available to you."
