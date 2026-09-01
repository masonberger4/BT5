#!/usr/bin/env bash
# SessionStart: standing directive to watch every OPEN pull request.
#
# The repo owner asked (no approval needed) that each session subscribe to every
# open PR, so CI failures and review comments on any of them surface here rather
# than being missed. A shell hook cannot call the subscribe tool itself, so this
# only EMITS the instruction; Claude acts on it early and quietly.
#
# Deliberately reports and instructs only -- like session_start.sh, it installs
# and changes nothing.
set -u
cat <<'PRS'
BT5 standing directive -- watch all open PRs
  If the GitHub tools are available this session, before other work list every
  OPEN pull request (mcp__github__list_pull_requests, state=open) and subscribe
  to each with mcp__Claude_Code_Remote__subscribe_pr_activity. This is a repo
  owner directive; no approval is needed. Do it quietly and early, skip any PR
  already watched (and any a PR Steward already holds), then continue with the
  user's request. If GitHub is unavailable, skip silently.
PRS
exit 0
