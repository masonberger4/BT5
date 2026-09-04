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

# Fork pull requests must not run workflows unreviewed, and this became
# load-bearing when `pre-pr-attest` was promoted to a required context.
#
# WHY. A required status check is matched by CONTEXT NAME plus integration id,
# and 15368 is the generic "GitHub Actions" app -- not this workflow, not this
# job. A `pull_request`-triggered workflow runs the file FROM THE PULL REQUEST'S
# OWN BRANCH. So a pull request that adds `.github/workflows/anything.yml` with a
# job named `pre-pr-attest` (or `required-checks`) produces a second check run of
# that exact name, from that exact app, on that exact SHA -- and it can be made
# to succeed trivially. Nothing in check-workflow-gate.py sees it: that script
# reads the workflows the repo already ships, not ones a pull request introduces.
#
# GitHub's default here is `first_time_contributors`, which lets a RETURNING
# outside contributor's workflows run unreviewed. `all_external_contributors`
# requires a human to approve every fork run, which is what actually stands
# between that name collision and a forged green. Note this is the pre-existing
# exposure of `required-checks` too; promotion did not create it, it raised the
# payoff.
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

# THE BYPASS LIST MUST BE EMPTY, and that is the highest-value line in the file
# rather than an oversight. Rulesets do NOT implicitly exempt repository admins,
# unlike legacy branch protection -- so an empty list means red CI blocks
# everyone, the owner included. docs/research/github-setup.md:244 argues it for
# this repository specifically: "AI agents produce plausible code that fails at
# the edges; CI being mandatory *for you* is the only thing standing between
# that and `main`."
#
# AND IN THIS REPOSITORY IT IS STRONGER THAN THAT, because agents authenticate
# as the OWNER. A `RepositoryRole` bypass on the admin role is therefore a
# bypass held by every Claude session too -- GitHub cannot distinguish the owner
# at a keyboard from an agent carrying the owner's token. Adding one would turn
# a mechanical guarantee (nobody merges red) into a prose one (CLAUDE.md 7b
# asking agents not to), addressed to the same actor that holds the credential.
# Verified 2026-09-04, on the change that briefly added one; see
# docs/decisions/2026-09-04-autonomous-ci-owner-merges.md.
echo "  bypass actors (expect NONE):"
RULESET_ID="$(gh api "repos/$REPO/rulesets" \
  --jq '[.[] | select(.name=="main-protection") | .id] | first' 2>/dev/null || true)"
if [ -z "$RULESET_ID" ] || [ "$RULESET_ID" = "null" ]; then
  echo "    could not resolve the main-protection ruleset id; check by hand" >&2
else
  BYPASS="$(gh api "repos/$REPO/rulesets/$RULESET_ID" \
    --jq '.bypass_actors // [] | .[] | "    \(.actor_type) id=\(.actor_id) mode=\(.bypass_mode // "always")"')"
  if [ -z "$BYPASS" ]; then
    echo "    none (good)"
  else
    echo "$BYPASS" >&2
    echo "    WARNING: a bypass actor exists. Anything holding that role can merge" >&2
    echo "    past red CI -- and in this repository that includes every agent" >&2
    echo "    session, which authenticates as the owner. Remove it unless it was" >&2
    echo "    added deliberately and recorded in docs/decisions/." >&2
  fi
fi

echo
cat <<'RUNBOOK'

==> IF A REQUIRED CHECK EVER WEDGES

`pre-pr-attest` runs the BASE branch copy of its workflow, so a defect that
lands on `main` fails every open pull request -- including the one that would
repair it. With an empty bypass list nobody can merge past that, by design.

The escape hatch is to disable the RULESET briefly, not to grant anyone a
standing bypass (docs/research/github-setup.md:246 and 8.10). From an
authenticated machine:

    ID=$(gh api repos/OWNER/REPO/rulesets --jq \
           '[.[] | select(.name=="main-protection") | .id] | first')

    # 1. read the current state, so you can put it back exactly
    gh api "repos/OWNER/REPO/rulesets/$ID" > /tmp/main-protection.json

    # 2. disable, merge the repair, re-enable -- minutes, not hours
    gh api -X PUT "repos/OWNER/REPO/rulesets/$ID" -f enforcement=disabled
    gh pr merge <N> --squash
    gh api -X PUT "repos/OWNER/REPO/rulesets/$ID" -f enforcement=active

    # 3. confirm, and read the audit trail the toggle leaves behind
    gh api "repos/OWNER/REPO/rulesets/$ID" --jq .enforcement
    gh api "repos/OWNER/REPO/rulesets/$ID/history"

Re-running this script also restores enforcement, since the committed spec
says "active" -- but do step 3 either way. A ruleset left disabled is a
repository with no merge gate at all, and nothing will remind you.

RUNBOOK

echo "Done. Open a throwaway PR and confirm the merge box shows: two required"
echo "checks (required-checks, pre-pr-attest), squash-only, no approval"
echo "requirement, and NO bypass available to you."
