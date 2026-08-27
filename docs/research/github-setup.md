# GitHub Setup Plan — `masonberger4/BT5`

**Context assumed:** public repo, personal free account, `masonberger4` sole human, many Claude Code sessions opening PRs in parallel. Monorepo: `packages/engine` (Python 3.11 — numpy / biopython / ViennaRNA), `packages/server` (FastAPI), `apps/web` (React + TS). Toolchain: `uv` + `pnpm`.

**The three decisions everything else follows from:**

1. **Required approvals = 0.** GitHub forbids self-approval. Any non-zero count on a solo repo is a permanent deadlock. Human sign-off is re-implemented as a *status check* (`agent-approval-check`, `/approve <sha>`), not as a review.
2. **One aggregating required check per workflow.** A path-filtered workflow that never triggers never reports, and a required check that never reports blocks the PR forever with no error. Every workflow that owns a required check runs unconditionally and short-circuits *inside* jobs.
3. **Loose status checks (`strict: false`).** Merge queue is organization-only. Strict mode with N parallel agent PRs is O(N²) CI and a livelock. Loose + squash-only + a post-merge `main-is-green` safety net is the correct trade.

---

## 1. Repo settings checklist

Grouped by the Settings sidebar. `gh` one-shot script at the end of the section.

### Settings → General

| Setting | Value | Note |
|---|---|---|
| Default branch | `main` | Confirm. |
| Features → Wikis | **off** | |
| Features → Projects | **off** | |
| Features → Issues | **on** | Needed by `main-is-green` and nightly failure reporting. |
| Features → Discussions | off | |
| Pull Requests → Allow merge commits | **off** | |
| Pull Requests → Allow squash merging | **on** | |
| Pull Requests → Allow rebase merging | **off** | |
| Squash commit title | **Pull request title** (`PR_TITLE`) | |
| Squash commit message | **Pull request description** (`PR_BODY`) | |
| Always suggest updating PR branches | **on** | Gives you a one-click "Update branch" *without* requiring up-to-date. This is the loose-mode companion. |
| Allow auto-merge | **on** | Your merge-queue substitute. |
| Automatically delete head branches | **on** | 8 agents × branches otherwise. |
| Allow forking | **on** (leave default) | Public repo. |

⚠️ Repo settings and the ruleset must agree on merge methods — if the repo disables a method the ruleset requires, **all merges block**.

### Settings → Collaborators and teams

Leave **empty**. Do not adopt the second-account workaround; `agent-approval-check` solves the same problem without a second login. (Second account is the only path to a *genuine* `required_approving_review_count: 1` if you ever want it.)

### Settings → Branches

Zero classic branch protection rules. Verify: `gh api repos/masonberger4/BT5/branches/main/protection` should 404. Rulesets and classic rules **stack, most-restrictive-wins**, and debugging the union is miserable.

### Settings → Rules → Rulesets

Two rulesets — `main-protection` (§2) and `release-tags` (§2). Bypass list **empty** on both.

### Settings → Actions → General

| Setting | Value | Why |
|---|---|---|
| Actions permissions | Allow all actions and reusable workflows | Fine for public; tighten to an allowlist if agents start adding third-party actions. |
| Fork PR workflows from outside collaborators | **Require approval for all external contributors** | Public repo, drive-by PRs will happen. |
| Workflow permissions | **Read repository contents and packages permissions** | Grant writes per-workflow in YAML. |
| Allow GitHub Actions to create and approve pull requests | **UNCHECKED** | 🔴 Security-critical. If on, an agent that can edit `.github/workflows` can write a workflow that approves its own PR. Off by default on new personal repos — *verify it stayed off*. |
| Artifact and log retention | 30 days | Public repos cap at 90. Archive anything you need for reproducing an optimization run. |

### Settings → Code security (Advanced Security)

| Setting | Value |
|---|---|
| Secret scanning | **on** (free on public) |
| **Push protection** | **on** — 🔴 disabled by default at repo level even when secret scanning is on. Verify the toggle. |
| Dependabot alerts | **on** |
| Dependabot security updates | **on** |
| Code scanning → CodeQL | **Advanced setup only** (committed `codeql.yml`) |

⚠️ Never enable CodeQL **default setup** alongside the advanced workflow — default setup *overrides* existing configuration, disables your workflow file, and blocks CodeQL API uploads. Advanced setup is required here because you want the `actions` language (scans your own workflows) and `security-extended`.

### Settings → Webhooks / Deploy keys / Environments / Secrets

- **Webhooks:** none. Do not let an agent add one.
- **Deploy keys:** none. Agents use a PAT or GitHub App so access is revocable per identity.
- **Environments:** empty (desktop app, no deploy target). Later: a `release` environment for a PyPI/npm token — never repo-level.
- **Secrets → Actions:** minimum. `CLAUDE_CODE_OAUTH_TOKEN` (or `ANTHROPIC_API_KEY`) only. On a public repo secrets are not exposed to fork workflows, but any workflow you merge can exfiltrate them.

### Settings → Copilot / Codespaces / Pages / Moderation / Planning

- **Copilot code review:** skip. It **can never block a merge** ("Copilot always leaves a *Comment* review… will not block merging changes") and costs 13 premium requests per review — ~23 reviews/month on Copilot Pro. You may qualify for free Pro as an OSS maintainer; even then, don't build safety on it.
- **Codespaces / Pages:** defaults.
- **Moderation:** interaction limits off; know where they are if the public repo attracts spam.

### 🚫 Unavailable on a free personal account

| Feature | Status |
|---|---|
| **Merge queue** | Org-owned repos only. This is the big one — it is the canonical fix for N parallel PRs. |
| **Ruleset "Evaluate" (dry-run) mode + Rule Insights** | Team/Enterprise only. You cannot dry-run a ruleset — create it Active and fix forward. |
| **"Restrict who can push to matching branches"** | Org-owned repos only. |
| **"Required reviewers" (per-path team reviewers)** | Explicitly "unavailable for user-owned repositories." |
| **Organization-level rulesets** | Org only. *Repository* rulesets work fine on free personal public repos — don't misread the docs. |
| **Delegated bypass for push protection** | Org-oriented. As sole owner you can always bypass push protection with a reason. |

> **Strongly consider transferring BT5 into a free GitHub Organization you own.** It costs nothing, redirects are automatic, and it unlocks merge queue — which is the architecturally correct answer to 8 parallel agent PRs. Do it *before* you accumulate forks/stars. Every workflow below already carries `merge_group:` so the move is a ruleset edit and nothing else.

### One-shot settings script

```bash
#!/usr/bin/env bash
set -euo pipefail
REPO=masonberger4/BT5

# General → merges & features
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
  -F has_issues=true

# Actions → read-only default token; workflows may NOT create/approve PRs
gh api -X PUT "repos/$REPO/actions/permissions/workflow" \
  -f default_workflow_permissions=read \
  -F can_approve_pull_request_reviews=false

# Verify the Actions app id used for pinning required checks (expect 15368)
gh api "repos/$REPO/commits/main/check-runs" --jq '.check_runs[] | {name, app_id: .app.id}'

# Confirm no legacy branch protection exists
gh api "repos/$REPO/branches/main/protection" >/dev/null 2>&1 \
  && echo "WARNING: legacy protection present — delete it" \
  || echo "OK: no legacy protection"
```

Secret scanning, push protection, and Dependabot are toggled in the UI (Settings → Code security).

---

## 2. The ruleset

### `.github/rulesets/main-protection.json`

```json
{
  "name": "main-protection",
  "target": "branch",
  "enforcement": "active",
  "bypass_actors": [],
  "conditions": {
    "ref_name": {
      "include": ["~DEFAULT_BRANCH"],
      "exclude": []
    }
  },
  "rules": [
    { "type": "deletion" },
    { "type": "non_fast_forward" },
    {
      "type": "pull_request",
      "parameters": {
        "required_approving_review_count": 0,
        "dismiss_stale_reviews_on_push": true,
        "require_code_owner_review": false,
        "require_last_push_approval": false,
        "required_review_thread_resolution": true,
        "allowed_merge_methods": ["squash"]
      }
    },
    {
      "type": "required_status_checks",
      "parameters": {
        "strict_required_status_checks_policy": false,
        "do_not_enforce_on_create": true,
        "required_status_checks": [
          { "context": "required-checks",   "integration_id": 15368 },
          { "context": "codeql-passed",     "integration_id": 15368 },
          { "context": "dependency-review", "integration_id": 15368 }
        ]
      }
    }
  ]
}
```

Phase 3 adds two more contexts (see §7):

```json
{ "context": "agent-approval-check", "integration_id": 15368 },
{ "context": "claude-review-gate",   "integration_id": 15368 }
```

### `.github/rulesets/release-tags.json`

```json
{
  "name": "release-tags",
  "target": "tag",
  "enforcement": "active",
  "bypass_actors": [],
  "conditions": { "ref_name": { "include": ["refs/tags/v*"], "exclude": [] } },
  "rules": [
    { "type": "deletion" },
    { "type": "non_fast_forward" },
    { "type": "update" }
  ]
}
```

### Creating them

```bash
REPO=masonberger4/BT5

gh api -X POST "repos/$REPO/rulesets" --input .github/rulesets/main-protection.json
gh api -X POST "repos/$REPO/rulesets" --input .github/rulesets/release-tags.json

# Verify
gh api "repos/$REPO/rulesets" --jq '.[] | {id, name, target, enforcement}'
gh api "repos/$REPO/rules/branches/main" --jq '.[].type'

# Update later — note PUT, not PATCH
# gh api -X PUT "repos/$REPO/rulesets/<ID>" --input .github/rulesets/main-protection.json
# Audit trail
# gh api "repos/$REPO/rulesets/<ID>/history"
```

Raw REST equivalent:

```
POST https://api.github.com/repos/masonberger4/BT5/rulesets
Accept: application/vnd.github+json
Authorization: Bearer <token with repo admin>
X-GitHub-Api-Version: 2022-11-28
Body: <main-protection.json>
```

### What it enforces, in plain English

- **Nothing reaches `main` except through a pull request.** No direct pushes, from you or from any agent.
- **`main` cannot be force-pushed or deleted.** A confused Claude Code session cannot rewrite history and destroy other agents' merged work.
- **The only merge button is Squash.** Each module change lands as one commit — trivial to revert, tractable to bisect. Agent branches' "fix lint / try again" noise never enters `main`.
- **Every review conversation must be resolved before the merge button unlocks.** This is the one review-flavored control that actually works for a solo owner — you comment, the agent fixes, the thread closes. A speed bump, not a wall (you can resolve your own threads), but it stops "merged and forgot the comment."
- **Three status checks must be green, and they must come from GitHub Actions specifically.** Pinning `integration_id` matters: without it, "any person or integration with write access can set status check states" — a token with `statuses:write` could post a green `required-checks` without CI ever running.
- **Branches do NOT have to be up to date.** Deliberate. See §8.
- **`v*` tags are immutable.** A published optimization result stays reproducible.
- **The bypass list is empty — including for you.** 🔴 This is the single highest-value setting here. Rulesets do **not** implicitly exempt repo admins (the opposite of legacy branch protection). Red CI genuinely blocks the one human who can click Merge. AI agents produce plausible code that fails at the edges; CI being mandatory *for you* is the only thing standing between that and `main`.

**Emergency escape hatch — do not add a bypass actor.** Bypass is per-*actor*, not per-*rule*: granting yourself admin bypass to get past one annoying rule also lets you merge red CI, force-push `main`, and delete the branch. Instead create a third ruleset named `emergency-bypass` with `"enforcement": "disabled"`, identical targeting, and no rules; when you truly need to bypass, set `main-protection` to `disabled` for five minutes and flip it back. That leaves an audit trail in `/rulesets/{id}/history`.

### The solo-owner self-approval problem — solved explicitly

**The constraint:** GitHub states flatly that "pull request authors cannot approve their own pull requests." No setting, permission, or API call changes this. Your agents authenticate as you, so every PR is authored by you.

**What that rules out:**

- `required_approving_review_count: 1` → mathematically unmergeable. You'd add yourself to the bypass list, and a bypass used on 100% of PRs is not a control — it silently disables *every* rule in the ruleset for the only person who merges.
- `require_last_push_approval: true` → demands approval from someone other than the last pusher. You are always the last pusher. Guaranteed lockout. Second most common way solo owners brick their own repo.
- `require_code_owner_review: true` → code owners need write access, so you're the only possible owner, and you're the author. No-op trap.

**What we do instead — a three-layer substitute:**

| Layer | Mechanism | What it buys |
|---|---|---|
| 1. Machine gate | `required_status_checks` with an empty bypass list | CI is genuinely mandatory. Nothing merges red. |
| 2. Conversation gate | `required_review_thread_resolution: true` | An enforceable checklist on agent PRs. |
| 3. **Human sign-off gate** | `anthropics/claude-code-action/agent-approval-check` as a **status check**, `required_approvals: 1` | A deliberate, timeline-visible human act that a solo owner *can* perform. |

Layer 3 is the real answer. The sub-action ("the same gate Anthropic runs internally on every agent-authored PR") counts approvals itself and accepts a **`/approve <head-sha>`** comment from the PR author — it exists precisely because "the PR author can't approve their own PR in GitHub's UI." It only engages on PRs containing agent-authored commits, **re-arms on every push** (a new head SHA stales prior approvals), verifies approvers actually have write access, treats >100-commit PRs as agent-authored (commit list unverifiable), and **fails closed** on any error.

Your merge ritual becomes: read the diff → comment `/approve <sha>` → click Merge (or Enable auto-merge). One deliberate act per PR, which is exactly the checkpoint that a solo-owner repo otherwise lacks.

**Land this in Phase 2, not Phase 0** — it depends on the Claude GitHub App being installed, and you want to watch it behave before it can block you.

---

## 3. Workflow files

### `.github/workflows/ci.yml`

The main PR workflow. **Owns exactly one required check: `required-checks`.** No workflow-level `paths:` filter anywhere — that is non-negotiable.

```yaml
name: CI

on:
  pull_request:
    types: [opened, synchronize, reopened, ready_for_review, labeled, unlabeled]
  push:
    branches: [main]
  merge_group:          # pre-wired for a future org move + merge queue
  workflow_dispatch:

# DELIBERATELY no `paths:` filter. A workflow that never triggers emits no
# check run, and a required check that never reports blocks the PR forever
# with no error and no timeout. All path filtering happens INSIDE jobs.

permissions:
  contents: read
  pull-requests: write   # sticky benchmark comment

concurrency:
  group: ci-${{ github.event.pull_request.number || github.ref }}
  cancel-in-progress: ${{ github.event_name == 'pull_request' }}

env:
  PYTHON_VERSION: "3.11"      # load-bearing: ViennaRNA ships no cp312/cp313 wheels
  NODE_VERSION: "22"
  PNPM_VERSION: "10"
  UV_VERSION: "0.12.6"
  PYTHONHASHSEED: "0"         # set/dict order must not vary between runs
  FORCE_COLOR: "1"

jobs:
  # ------------------------------------------------------------------
  # Path filtering + the label path-guard. Cheap, always runs, fails OPEN.
  # ------------------------------------------------------------------
  changes:
    name: changes
    runs-on: ubuntu-latest
    timeout-minutes: 5
    outputs:
      python:  ${{ steps.decide.outputs.python }}
      web:     ${{ steps.decide.outputs.web }}
      e2e:     ${{ steps.decide.outputs.e2e }}
      engine:  ${{ steps.decide.outputs.engine }}
      python_os: ${{ steps.decide.outputs.python_os }}
    steps:
      - uses: actions/checkout@v7
        with: { fetch-depth: 0 }

      - name: Filter changed paths
        id: filter
        if: github.event_name == 'pull_request'
        uses: dorny/paths-filter@v4
        with:
          filters: |
            shared: &shared
              - '.github/workflows/ci.yml'
              - 'pyproject.toml'
              - 'uv.lock'
            engine:
              - *shared
              - 'packages/engine/**'
              - 'data/**'
              - 'benchmarks/**'
            python:
              - *shared
              - 'packages/engine/**'
              - 'packages/server/**'
              - 'data/**'
              - 'benchmarks/**'
            web:
              - *shared
              - 'apps/web/**'
              - 'package.json'
              - 'pnpm-lock.yaml'
              - 'pnpm-workspace.yaml'
            e2e:
              - *shared
              - 'packages/server/**'
              - 'apps/web/**'
              - 'tests/e2e/**'

      # FAIL OPEN. An empty output must mean "run it", never "skip it" --
      # a skip is swallowed by the gate and reads as green.
      - name: Decide what to run
        id: decide
        shell: bash
        env:
          EVENT:      ${{ github.event_name }}
          F_PY:       ${{ steps.filter.outputs.python }}
          F_WEB:      ${{ steps.filter.outputs.web }}
          F_E2E:      ${{ steps.filter.outputs.e2e }}
          F_ENGINE:   ${{ steps.filter.outputs.engine }}
          FULL_LABEL: ${{ contains(github.event.pull_request.labels.*.name, 'ci-full') }}
        run: |
          set -euo pipefail
          if [ "$EVENT" = "pull_request" ] && [ "$FULL_LABEL" != "true" ]; then
            py="${F_PY:-true}"; web="${F_WEB:-true}"
            e2e="${F_E2E:-true}"; engine="${F_ENGINE:-true}"
            os='["ubuntu-latest"]'
          else
            py=true; web=true; e2e=true; engine=true
            os='["ubuntu-latest","macos-latest","windows-latest"]'
          fi
          {
            echo "python=$py"; echo "web=$web"; echo "e2e=$e2e"
            echo "engine=$engine"; echo "python_os=$os"
          } >> "$GITHUB_OUTPUT"

      - name: Load-bearing paths require an owner label
        if: github.event_name == 'pull_request'
        env:
          LABELS: ${{ toJSON(github.event.pull_request.labels.*.name) }}
          BASE:   ${{ github.event.pull_request.base.sha }}
          HEAD:   ${{ github.event.pull_request.head.sha }}
        run: python .github/scripts/path_guard.py

  # ------------------------------------------------------------------
  python-quality:
    name: python-quality
    needs: changes
    if: needs.changes.outputs.python == 'true'
    runs-on: ubuntu-latest
    timeout-minutes: 10
    steps:
      - uses: actions/checkout@v7
      - uses: astral-sh/setup-uv@v10
        with:
          version: ${{ env.UV_VERSION }}
          python-version: ${{ env.PYTHON_VERSION }}
          enable-cache: true
          cache-dependency-glob: "uv.lock"

      # --locked, NOT --frozen. --frozen uses the lockfile without checking it
      # matches pyproject.toml, so an agent that adds a dep and forgets
      # `uv lock` gets a silently wrong env and a green CI.
      - name: Sync (fails if uv.lock is stale)
        run: uv sync --locked --all-packages --all-extras --dev

      # Per-SHA key + prefix restore-key so the cache always re-saves.
      # A fixed key never re-saves after the first hit and goes stale silently.
      - uses: actions/cache@v6
        with:
          path: .mypy_cache
          key: mypy-${{ runner.os }}-py${{ env.PYTHON_VERSION }}-${{ hashFiles('uv.lock') }}-${{ github.sha }}
          restore-keys: |
            mypy-${{ runner.os }}-py${{ env.PYTHON_VERSION }}-${{ hashFiles('uv.lock') }}-

      - name: Ruff lint
        run: uv run ruff check --output-format=github .
      - name: Ruff format check
        if: ${{ !cancelled() }}
        run: uv run ruff format --check --diff .
      - name: Mypy (strict)
        if: ${{ !cancelled() }}
        run: uv run mypy .

      # Agents suppress rather than fix. Specific codes are fine; blanket is not.
      - name: Ban blanket suppressions
        if: ${{ !cancelled() }}
        run: |
          if git grep -nE '(#\s*type:\s*ignore(\s|$))|(#\s*noqa(\s|$))|(#\s*pragma:\s*no cover)' -- 'packages/**/*.py'; then
            echo "::error::Use a specific code (# type: ignore[arg-type], # noqa: E501) and justify it."
            exit 1
          fi

      - name: Ban global RNG state in src/
        if: ${{ !cancelled() }}
        run: |
          if grep -rnE '(^|[^.\w])(np|numpy)\.random\.(seed|rand|randn|choice|randint)|(^|[^.\w])random\.seed\(' packages/*/src/; then
            echo "::error::Use np.random.default_rng(seed) / random.Random(seed). Global RNG breaks reproducibility."
            exit 1
          fi

      - name: Verification must default to ON
        if: ${{ !cancelled() }}
        run: |
          grep -q '_verify: bool = True' packages/engine/src/bt5/optimizer.py || {
            echo "::error::optimize() must verify by default"; exit 1; }

      - name: Prune uv cache
        if: always()
        run: uv cache prune --ci

  # ------------------------------------------------------------------
  python-tests:
    name: python-tests (${{ matrix.os }})
    needs: changes
    if: needs.changes.outputs.python == 'true'
    runs-on: ${{ matrix.os }}
    timeout-minutes: 20
    strategy:
      fail-fast: false     # a Windows-only break must not be hidden by macOS
      matrix:
        os: ${{ fromJSON(needs.changes.outputs.python_os) }}
    steps:
      - uses: actions/checkout@v7
      - uses: astral-sh/setup-uv@v10
        with:
          version: ${{ env.UV_VERSION }}
          python-version: ${{ env.PYTHON_VERSION }}
          enable-cache: true
          cache-dependency-glob: "uv.lock"
      - run: uv sync --locked --all-packages --all-extras --dev

      # Fails in 2s with a clear message if a wheel is missing, instead of
      # 8 minutes into a source build of a large C library.
      - name: Smoke-import native deps
        run: uv run python -c "import numpy, RNA, Bio; print(numpy.__version__, RNA.__version__)"

      - name: Unit tests + data integrity + goldens
        run: >-
          uv run pytest packages/*/tests tests/unit tests/data_integrity tests/goldens
          -m "not slow"
          -n logical --dist loadfile
          --snapshot-warn-unused
          --cov --cov-report=xml --cov-report=term-missing:skip-covered
          --cov-fail-under=85
          --durations=15 -q -W error::Bio.BiopythonWarning

      - name: Upload coverage
        if: ${{ !cancelled() && matrix.os == 'ubuntu-latest' }}
        uses: actions/upload-artifact@v7
        with:
          name: coverage-xml
          path: coverage.xml
          if-no-files-found: error
          retention-days: 7

      - name: Prune uv cache
        if: always()
        run: uv cache prune --ci

  # ------------------------------------------------------------------
  invariants:
    name: invariants
    needs: changes
    if: needs.changes.outputs.engine == 'true'
    runs-on: ubuntu-latest
    timeout-minutes: 20
    steps:
      - uses: actions/checkout@v7
      - uses: astral-sh/setup-uv@v10
        with:
          version: ${{ env.UV_VERSION }}
          python-version: ${{ env.PYTHON_VERSION }}
          enable-cache: true
          cache-dependency-glob: "uv.lock"
      - run: uv sync --locked --all-packages --all-extras --dev

      # Replays previously-found counterexamples first. Without this a bug
      # found yesterday is rediscovered only by luck.
      - uses: actions/cache@v6
        with:
          path: .hypothesis
          key: hypothesis-${{ github.run_id }}
          restore-keys: hypothesis-

      - name: Property-based invariants
        env: { HYPOTHESIS_PROFILE: ci }
        run: uv run pytest tests/invariants -q -W error::Bio.BiopythonWarning

  # ------------------------------------------------------------------
  goldens-not-hand-edited:
    name: goldens-not-hand-edited
    needs: changes
    if: needs.changes.outputs.engine == 'true'
    runs-on: ubuntu-latest
    timeout-minutes: 15
    steps:
      - uses: actions/checkout@v7
      - uses: astral-sh/setup-uv@v10
        with:
          version: ${{ env.UV_VERSION }}
          python-version: ${{ env.PYTHON_VERSION }}
          enable-cache: true
          cache-dependency-glob: "uv.lock"
      - run: uv sync --locked --all-packages --all-extras --dev
      - name: Regenerate goldens from scratch and diff
        run: |
          cp -r tests/goldens/__snapshots__ /tmp/committed
          rm -rf tests/goldens/__snapshots__
          uv run pytest tests/goldens --snapshot-update -q
          if ! diff -ru /tmp/committed tests/goldens/__snapshots__ > /tmp/g.diff; then
            {
              echo "### Committed goldens do not match freshly generated output"
              echo "A snapshot file appears to have been hand-edited."
              echo '```diff'; head -c 60000 /tmp/g.diff; echo '```'
            } >> "$GITHUB_STEP_SUMMARY"
            exit 1
          fi

  # ------------------------------------------------------------------
  benchmark-gate:
    name: benchmark-gate
    needs: changes
    if: needs.changes.outputs.engine == 'true'
    runs-on: ubuntu-latest    # folding dG is compared on ONE architecture only
    timeout-minutes: 30
    steps:
      - uses: actions/checkout@v7
      - uses: astral-sh/setup-uv@v10
        with:
          version: ${{ env.UV_VERSION }}
          python-version: ${{ env.PYTHON_VERSION }}
          enable-cache: true
          cache-dependency-glob: "uv.lock"
      - run: uv sync --locked --all-packages --all-extras --dev
      - name: Run fast panel and compare to committed baseline
        run: |
          uv run python benchmarks/run_panel.py \
            --tier fast \
            --out benchmarks/results.json \
            --compare benchmarks/baseline.json \
            --summary "$GITHUB_STEP_SUMMARY"
      - uses: actions/upload-artifact@v7
        if: always()
        with: { name: benchmark-results, path: benchmarks/results.json, retention-days: 14 }
      - name: Post diff table as a sticky PR comment
        if: always() && github.event_name == 'pull_request'
        uses: marocchino/sticky-pull-request-comment@v2
        with:
          header: benchmark
          path: ${{ env.GITHUB_STEP_SUMMARY }}

  # ------------------------------------------------------------------
  web:
    name: web
    needs: changes
    if: needs.changes.outputs.web == 'true'
    runs-on: ubuntu-latest
    timeout-minutes: 12
    steps:
      - uses: actions/checkout@v7
      # pnpm MUST precede setup-node: setup-node shells out to `pnpm store path`.
      - uses: pnpm/action-setup@v6
        with: { version: "${{ env.PNPM_VERSION }}" }
      - uses: actions/setup-node@v7
        with:
          node-version: ${{ env.NODE_VERSION }}
          cache: pnpm
          cache-dependency-path: pnpm-lock.yaml
      - run: pnpm install --frozen-lockfile
      - name: Typecheck
        run: pnpm -r exec tsc --noEmit
      - name: ESLint
        if: ${{ !cancelled() }}
        run: pnpm exec eslint . --max-warnings=0
      - name: Prettier
        if: ${{ !cancelled() }}
        run: pnpm exec prettier --check .
      - name: Vitest
        if: ${{ !cancelled() }}
        run: pnpm exec vitest run --coverage --reporter=dot --reporter=github-actions

  # ------------------------------------------------------------------
  e2e:
    name: e2e
    needs: changes
    if: needs.changes.outputs.e2e == 'true'
    runs-on: ubuntu-latest
    timeout-minutes: 20
    steps:
      - uses: actions/checkout@v7
      - uses: pnpm/action-setup@v6
        with: { version: "${{ env.PNPM_VERSION }}" }
      - uses: actions/setup-node@v7
        with:
          node-version: ${{ env.NODE_VERSION }}
          cache: pnpm
          cache-dependency-path: pnpm-lock.yaml
      - uses: astral-sh/setup-uv@v10
        with:
          version: ${{ env.UV_VERSION }}
          python-version: ${{ env.PYTHON_VERSION }}
          enable-cache: true
          cache-dependency-glob: "uv.lock"
      - run: pnpm install --frozen-lockfile
      - run: uv sync --locked --all-packages

      # Key on the Playwright version: Playwright pins exact browser builds.
      - name: Resolve Playwright version
        id: pw
        run: |
          v=$(pnpm --filter @bt5/web exec playwright --version | tr -dc '0-9.')
          echo "version=$v" >> "$GITHUB_OUTPUT"
      - uses: actions/cache@v6
        id: pw-cache
        with:
          path: ~/.cache/ms-playwright
          key: playwright-${{ runner.os }}-${{ steps.pw.outputs.version }}-chromium

      - name: Install browser + OS deps
        if: steps.pw-cache.outputs.cache-hit != 'true'
        run: pnpm --filter @bt5/web exec playwright install --with-deps chromium

      # The cache holds ~/.cache/ms-playwright but NOT the apt libs in /usr/lib.
      # Skipping this on a cache hit yields a browser that cannot launch.
      - name: Install OS deps only
        if: steps.pw-cache.outputs.cache-hit == 'true'
        run: pnpm --filter @bt5/web exec playwright install-deps chromium

      - run: pnpm --filter @bt5/web exec playwright test
      - uses: actions/upload-artifact@v7
        if: ${{ !cancelled() }}
        with:
          name: playwright-report
          path: apps/web/playwright-report/
          retention-days: 7

  # ==================================================================
  # THE GATE. The only name from this workflow in the ruleset.
  #
  #  * `if: always()` is REQUIRED. With success() (or no if:) the gate is
  #    SKIPPED when a dependency fails -- and GitHub accepts a skipped
  #    required check as satisfied. Failing tests would merge.
  #  * always() ALONE IS NOT ENOUGH. A job with no failing steps reports
  #    success. The explicit needs.*.result step below is the enforcement.
  #  * `skipped` is intentionally allowed -- that is what makes in-job path
  #    filtering safe -- which is why `changes` is asserted separately and
  #    written to fail open.
  #  * EVERY NEW JOB MUST BE ADDED TO needs: OR IT IS UNENFORCED.
  # ==================================================================
  required-checks:
    name: required-checks
    if: ${{ always() }}
    needs:
      - changes
      - python-quality
      - python-tests
      - invariants
      - goldens-not-hand-edited
      - benchmark-gate
      - web
      - e2e
    runs-on: ubuntu-latest
    timeout-minutes: 5
    steps:
      - name: Show dependency results
        env: { RESULTS: ${{ toJSON(needs) }} }
        run: echo "$RESULTS"

      - name: Fail if anything failed, was cancelled, or the gate lost its inputs
        if: >-
          ${{ needs.changes.result != 'success'
              || contains(needs.*.result, 'failure')
              || contains(needs.*.result, 'cancelled') }}
        run: |
          echo "::error::One or more required CI jobs did not succeed."
          exit 1

      - name: Success
        run: echo "All required jobs succeeded or were legitimately skipped."
```

### `.github/workflows/codeql.yml`

Same aggregator pattern, so the required name never changes when the language matrix does.

```yaml
name: CodeQL

on:
  pull_request:
    branches: [main]
  push:
    branches: [main]
  merge_group:
  schedule:
    - cron: "27 4 * * 1"
  workflow_dispatch:

permissions:
  contents: read

concurrency:
  group: codeql-${{ github.event.pull_request.number || github.ref }}
  cancel-in-progress: ${{ github.event_name == 'pull_request' }}

jobs:
  pick:
    name: pick-languages
    runs-on: ubuntu-latest
    timeout-minutes: 5
    outputs:
      langs: ${{ steps.decide.outputs.langs }}
    steps:
      - uses: actions/checkout@v7
      - uses: dorny/paths-filter@v4
        id: filter
        if: github.event_name == 'pull_request'
        with:
          filters: |
            python:  ['packages/engine/**','packages/server/**','**/*.py']
            js:      ['apps/web/**','**/*.ts','**/*.tsx','**/*.js']
            actions: ['.github/workflows/**','.github/actions/**']
      - id: decide
        shell: bash
        env:
          EVENT: ${{ github.event_name }}
          P: ${{ steps.filter.outputs.python }}
          J: ${{ steps.filter.outputs.js }}
          A: ${{ steps.filter.outputs.actions }}
        run: |
          set -euo pipefail
          if [ "$EVENT" != "pull_request" ]; then
            echo 'langs=["python","javascript-typescript","actions"]' >> "$GITHUB_OUTPUT"; exit 0
          fi
          out=()
          [ "${P:-true}" = "true" ] && out+=('"python"')
          [ "${J:-true}" = "true" ] && out+=('"javascript-typescript"')
          [ "${A:-true}" = "true" ] && out+=('"actions"')
          [ ${#out[@]} -eq 0 ] && out+=('"actions"')
          echo "langs=[$(IFS=,; echo "${out[*]}")]" >> "$GITHUB_OUTPUT"

  analyze:
    name: analyze (${{ matrix.language }})
    needs: pick
    runs-on: ubuntu-latest
    timeout-minutes: 30
    permissions:
      security-events: write
      packages: read
      actions: read
      contents: read
    strategy:
      fail-fast: false
      matrix:
        language: ${{ fromJSON(needs.pick.outputs.langs) }}
    steps:
      - uses: actions/checkout@v7
      - uses: github/codeql-action/init@v4
        with:
          languages: ${{ matrix.language }}
          build-mode: none          # valid for python and javascript-typescript
          queries: security-extended
      - uses: github/codeql-action/analyze@v4
        with:
          category: "/language:${{ matrix.language }}"

  codeql-passed:
    name: codeql-passed        # <- the required check name
    if: ${{ always() }}
    needs: [pick, analyze]
    runs-on: ubuntu-latest
    timeout-minutes: 5
    steps:
      - if: >-
          ${{ needs.pick.result != 'success'
              || contains(needs.*.result, 'failure')
              || contains(needs.*.result, 'cancelled') }}
        run: |
          echo "::error::CodeQL analysis did not succeed."
          echo '${{ toJSON(needs) }}'
          exit 1
      - run: echo "CodeQL clean."
```

### `.github/workflows/dependency-review.yml`

The one third-party check that **blocks by default** rather than commenting.

```yaml
name: Dependency Review

on:
  pull_request:
    branches: [main]
  merge_group:

permissions:
  contents: read
  pull-requests: write

jobs:
  dependency-review:
    name: dependency-review       # <- the required check name
    runs-on: ubuntu-latest
    timeout-minutes: 10
    steps:
      - uses: actions/checkout@v7
      # merge_group has no PR to diff; report success so the required check reports.
      - if: github.event_name == 'pull_request'
        uses: actions/dependency-review-action@v4
        with:
          fail-on-severity: moderate
          comment-summary-in-pr: on-failure
      - run: echo "dependency review complete"
```

⚠️ Dependency review is free on public repos only. If BT5 ever goes private this required check silently becomes unavailable and blocks every PR.

### `.github/workflows/main-is-green.yml`

The safety net that pays for running loose status checks.

```yaml
name: main-is-green

on:
  push:
    branches: [main]
  workflow_dispatch:

permissions:
  contents: read
  issues: write

concurrency:
  group: main-is-green
  cancel-in-progress: false      # never cancel a post-merge validation

env:
  PYTHONHASHSEED: "0"

jobs:
  full-suite:
    runs-on: ubuntu-latest
    timeout-minutes: 45
    steps:
      - uses: actions/checkout@v7
      - uses: astral-sh/setup-uv@v10
        with:
          version: "0.12.6"
          python-version: "3.11"
          enable-cache: true
          cache-dependency-glob: "uv.lock"
      - uses: pnpm/action-setup@v6
        with: { version: "10" }
      - uses: actions/setup-node@v7
        with: { node-version: "22", cache: pnpm }
      - run: uv sync --locked --all-packages --all-extras --dev
      - run: pnpm install --frozen-lockfile

      # NO path filters. This is the semantic-conflict detector: the whole
      # point is to catch two PRs that were each green against an older main.
      - name: Full Python suite
        run: uv run pytest -n logical --dist loadfile -q -W error::Bio.BiopythonWarning
      - name: Full invariants (deep profile)
        env: { HYPOTHESIS_PROFILE: ci }
        run: uv run pytest tests/invariants -q
      - name: Full benchmark panel vs baseline
        run: |
          uv run python benchmarks/run_panel.py --tier full \
            --out benchmarks/results.json --compare benchmarks/baseline.json \
            --summary "$GITHUB_STEP_SUMMARY"
      - name: Web build + typecheck
        run: pnpm -r exec tsc --noEmit && pnpm -r build

      - name: Open an issue if main is broken
        if: failure()
        uses: actions/github-script@v7
        with:
          script: |
            const sha = context.sha.slice(0,7);
            const { data: existing } = await github.rest.issues.listForRepo({
              owner: context.repo.owner, repo: context.repo.repo,
              state: 'open', labels: 'broken-main'
            });
            const url = `${context.serverUrl}/${context.repo.owner}/${context.repo.repo}/actions/runs/${context.runId}`;
            if (existing.length) {
              await github.rest.issues.createComment({
                owner: context.repo.owner, repo: context.repo.repo,
                issue_number: existing[0].number,
                body: `Still broken at \`${sha}\`: ${url}`
              });
            } else {
              await github.rest.issues.create({
                owner: context.repo.owner, repo: context.repo.repo,
                title: `main is broken at ${sha}`,
                body: `Post-merge full suite failed.\n\n${url}\n\nLikely a semantic conflict between two PRs that were each green against an older main.`,
                labels: ['broken-main','priority:high']
              });
            }
```

### `.github/workflows/nightly.yml`

```yaml
name: Nightly

on:
  schedule:
    - cron: "0 7 * * *"
  workflow_dispatch:      # scheduled workflows auto-disable after 60 days of
                          # repo inactivity -- always keep a manual trigger

permissions:
  contents: read
  issues: write

concurrency:
  group: nightly
  cancel-in-progress: false

env:
  PYTHONHASHSEED: "0"

jobs:
  full-matrix:
    name: full suite (${{ matrix.os }})
    runs-on: ${{ matrix.os }}
    timeout-minutes: 60
    strategy:
      fail-fast: false
      matrix:
        os: [ubuntu-latest, macos-latest, windows-latest]
    steps:
      - uses: actions/checkout@v7
      - uses: astral-sh/setup-uv@v10
        with: { version: "0.12.6", python-version: "3.11", enable-cache: true, cache-dependency-glob: "uv.lock" }
      - run: uv sync --locked --all-packages --all-extras --dev
      - run: uv run python -c "import numpy, RNA, Bio"
      - run: uv run pytest -n logical --dist loadfile --cov --cov-fail-under=85 -q

  deep-invariants:
    name: deep invariants
    runs-on: ubuntu-latest
    timeout-minutes: 90
    steps:
      - uses: actions/checkout@v7
      - uses: astral-sh/setup-uv@v10
        with: { version: "0.12.6", python-version: "3.11", enable-cache: true, cache-dependency-glob: "uv.lock" }
      - run: uv sync --locked --all-packages --all-extras --dev
      - uses: actions/cache@v6
        with:
          path: .hypothesis
          key: hypothesis-nightly-${{ github.run_id }}
          restore-keys: hypothesis-
      - env: { HYPOTHESIS_PROFILE: deep }
        run: uv run pytest tests/invariants -q

  differential:
    name: differential vs DNA Chisel
    runs-on: ubuntu-latest
    timeout-minutes: 45
    steps:
      - uses: actions/checkout@v7
      - uses: astral-sh/setup-uv@v10
        with: { version: "0.12.6", python-version: "3.11", enable-cache: true, cache-dependency-glob: "uv.lock" }
      - run: uv sync --locked --all-packages --all-extras --dev --extra diff
      - run: uv run pytest tests/differential -q -m slow

  mutation:
    name: mutation score (correctness spine)
    runs-on: ubuntu-latest
    timeout-minutes: 90
    steps:
      - uses: actions/checkout@v7
      - uses: astral-sh/setup-uv@v10
        with: { version: "0.12.6", python-version: "3.11", enable-cache: true, cache-dependency-glob: "uv.lock" }
      - run: uv sync --locked --all-packages --all-extras --dev
      - run: uv run mutmut run --max-children 4 || true
      - name: Enforce score floor
        run: |
          uv run mutmut export-cicd-stats > mutmut.json
          uv run python - <<'PY'
          import json, os, sys
          s = json.load(open("mutmut.json"))
          killed, total = s["killed"], s["killed"] + s["survived"]
          score = 100.0 * killed / max(total, 1)
          open(os.environ["GITHUB_STEP_SUMMARY"], "a").write(
              f"## Mutation score: {score:.1f}% ({killed}/{total})\n")
          sys.exit(0 if score >= 90 else 1)
          PY

  lock-drift:
    name: lockfile drift (upstream breakage early warning)
    runs-on: ubuntu-latest
    timeout-minutes: 20
    steps:
      - uses: actions/checkout@v7
      - uses: astral-sh/setup-uv@v10
        with: { version: "0.12.6", python-version: "3.11" }
      - run: uv lock --upgrade
      - run: uv sync --locked --all-packages --all-extras --dev
      - run: uv run pytest -q -m "not slow"

  audit:
    name: security audit
    runs-on: ubuntu-latest
    timeout-minutes: 20
    steps:
      - uses: actions/checkout@v7
      - uses: astral-sh/setup-uv@v10
        with: { version: "0.12.6", python-version: "3.11" }
      - run: uv run --with pip-audit pip-audit --strict
      - uses: pnpm/action-setup@v6
        with: { version: "10" }
      - uses: actions/setup-node@v7
        with: { node-version: "22", cache: pnpm }
      - run: pnpm install --frozen-lockfile
      - run: pnpm audit --audit-level=high

  gate-integrity:
    name: gate integrity
    runs-on: ubuntu-latest
    timeout-minutes: 5
    steps:
      - uses: actions/checkout@v7
      - name: Every ci.yml job must be reachable from required-checks
        run: python .github/scripts/assert_gate_complete.py

  e2e-all-browsers:
    name: playwright (${{ matrix.browser }} ${{ matrix.shardIndex }}/2)
    runs-on: ubuntu-latest
    timeout-minutes: 45
    strategy:
      fail-fast: false
      matrix:
        browser: [chromium, firefox, webkit]
        shardIndex: [1, 2]
    steps:
      - uses: actions/checkout@v7
      - uses: pnpm/action-setup@v6
        with: { version: "10" }
      - uses: actions/setup-node@v7
        with: { node-version: "22", cache: pnpm }
      - run: pnpm install --frozen-lockfile
      - run: pnpm --filter @bt5/web exec playwright install --with-deps ${{ matrix.browser }}
      - run: >-
          pnpm --filter @bt5/web exec playwright test
          --project=${{ matrix.browser }}
          --shard=${{ matrix.shardIndex }}/2

  report:
    name: report nightly failure
    if: failure()
    needs: [full-matrix, deep-invariants, differential, mutation, lock-drift, audit, gate-integrity, e2e-all-browsers]
    runs-on: ubuntu-latest
    steps:
      - uses: actions/github-script@v7
        with:
          script: |
            await github.rest.issues.create({
              owner: context.repo.owner, repo: context.repo.repo,
              title: `Nightly failed ${new Date().toISOString().slice(0,10)}`,
              body: `${context.serverUrl}/${context.repo.owner}/${context.repo.repo}/actions/runs/${context.runId}`,
              labels: ['nightly-failure']
            })
```

### `.github/workflows/agent-approval-check.yml` (Phase 2)

```yaml
# Requires an explicit human sign-off on any PR containing agent-authored
# commits. BOTH triggers run the workflow file from the BASE branch, so a PR
# cannot edit this check to approve itself.
#
# Solo-owner path (GitHub will not let you Approve your own PR in the UI):
#   comment  /approve <current-head-sha>
name: agent-approval-check

on:
  pull_request_target:
    types: [opened, synchronize, reopened, ready_for_review]
  issue_comment:
    types: [created]

# NEVER use pull_request_review here -- it runs from the merge ref, so the PR
# under review could edit the check to approve itself.

permissions:
  contents: read
  pull-requests: write
  statuses: write

jobs:
  check:
    if: github.event_name != 'issue_comment' || github.event.issue.pull_request
    runs-on: ubuntu-latest
    timeout-minutes: 10
    steps:
      # No checkout of the PR head. pull_request_target runs privileged.
      - uses: anthropics/claude-code-action/agent-approval-check@main
        with:
          required_approvals: 1
          agent_emails: noreply@anthropic.com
          agent_logins: claude[bot],claude-code[bot]
          excluded_approvers: dependabot[bot]
          protected_bases: main
```

### `.github/workflows/claude-review.yml` (Phase 3)

```yaml
name: Claude Review

on:
  pull_request:
    types: [opened, synchronize, reopened, ready_for_review]

concurrency:
  group: claude-review-${{ github.event.pull_request.number }}
  cancel-in-progress: true      # a burst of agent pushes pays for one review

permissions:
  contents: read
  pull-requests: write
  id-token: write

jobs:
  claude-review-gate:
    name: claude-review-gate      # <- the required check name
    if: github.event.pull_request.draft == false
    runs-on: ubuntu-latest
    timeout-minutes: 25
    steps:
      - uses: actions/checkout@v7
        with: { fetch-depth: 0 }

      - id: review
        uses: anthropics/claude-code-action@v1
        with:
          claude_code_oauth_token: ${{ secrets.CLAUDE_CODE_OAUTH_TOKEN }}
          use_sticky_comment: true
          prompt: |
            REPO: ${{ github.repository }}
            PR:   ${{ github.event.pull_request.number }}

            Reviewing an AI-agent-authored PR for BT5: a locally-run codon
            optimization / protein back-translation desktop app. Python 3.11
            engine (numpy, biopython, ViennaRNA), FastAPI local server,
            React + TypeScript front end, one monorepo.

            The PR branch is checked out. Read CLAUDE.md first.

            Review priorities, highest first:

            1. SCIENTIFIC CORRECTNESS
               - Codon tables must match the declared NCBI translation table id.
               - Back-translation must round-trip: translate(back_translate(p)) == p.
               - Reading frame, start/stop handling, ambiguity codes (N,R,Y),
                 selenocysteine / pyrrolysine edge cases.
               - CAI / GC / CpG / rare-codon math must match the docstring
                 formula. Flag any uncited hard-coded constant.
            2. NUMERICAL SAFETY
               - numpy dtype/overflow, divide-by-zero, log(0), windowed
                 off-by-one. Any non-deterministic iteration order or unseeded RNG.
            3. INTERFACE DRIFT
               - This PR owns ONE module. Verify it did not silently change a
                 shared Pydantic model, a FastAPI route signature, or a TS type
                 another module consumes.
            4. SECURITY
               - CORS restricted to localhost, no shell=True with user input,
                 no path traversal in upload/export, no eval/pickle of user data.
            5. TESTS
               - New logic needs tests. Flag tests that assert on the
                 implementation rather than on the biology.

            Post findings as inline comments; post a short summary with
            `gh pr comment`. Then return the JSON verdict.

            Set blocking=true ONLY for defects you are confident are real:
            wrong science, a crash, a security hole, a broken shared interface.
            Style, naming and "consider refactoring" are NEVER blocking.
            If you could NOT complete the review (ran out of turns, could not
            read the diff), set blocking=true and say so.

          claude_args: |
            --model claude-sonnet-5
            --max-turns 30
            --allowedTools "mcp__github_inline_comment__create_inline_comment,Bash(gh pr comment:*),Bash(gh pr diff:*),Bash(gh pr view:*)"
            --json-schema '{"type":"object","properties":{"blocking":{"type":"boolean"},"summary":{"type":"string"},"blocking_findings":{"type":"array","items":{"type":"string"}}},"required":["blocking","summary","blocking_findings"],"additionalProperties":false}'

      # THIS is what makes the review blocking. The action's own `conclusion`
      # output means "did Claude run", not "did the review pass".
      - name: Gate on Claude verdict
        if: always()
        env: { OUT: "${{ steps.review.outputs.structured_output }}" }
        run: |
          set -euo pipefail
          if [ -z "${OUT:-}" ]; then
            echo "::error::No verdict (crash, timeout, or --max-turns). Failing closed."
            exit 1
          fi
          echo "$OUT" | jq .
          if [ "$(echo "$OUT" | jq -r '.blocking')" = "true" ]; then
            echo "::error::Claude found blocking issues:"
            echo "$OUT" | jq -r '.blocking_findings[]'
            exit 1
          fi
          echo "Claude review passed: $(echo "$OUT" | jq -r '.summary')"
```

### `.github/workflows/claude.yml` (`@claude` mentions)

```yaml
name: Claude

on:
  issue_comment:
    types: [created]
  pull_request_review_comment:
    types: [created]
  issues:
    types: [opened, assigned]

jobs:
  claude:
    if: |
      (github.event_name == 'issue_comment' && contains(github.event.comment.body, '@claude')) ||
      (github.event_name == 'pull_request_review_comment' && contains(github.event.comment.body, '@claude')) ||
      (github.event_name == 'issues' && contains(github.event.issue.body, '@claude'))
    runs-on: ubuntu-latest
    timeout-minutes: 30
    permissions:
      contents: write
      pull-requests: write
      issues: write
      id-token: write
    steps:
      - uses: actions/checkout@v7
        with: { fetch-depth: 0 }
      - uses: anthropics/claude-code-action@v1
        with:
          claude_code_oauth_token: ${{ secrets.CLAUDE_CODE_OAUTH_TOKEN }}
          claude_args: |
            --model claude-sonnet-5
            --max-turns 40
```

⚠️ Leave `allowed_bots` at its empty default. On a public repo `*` lets external Apps invoke this action with prompts they control.

### `.github/workflows/regen-goldens.yml`

Makes updating goldens a deliberate, reviewable ritual instead of `--snapshot-update` on an agent's laptop.

```yaml
name: Regenerate goldens

on:
  workflow_dispatch:
    inputs:
      reason:
        description: "Why is the algorithm output changing?"
        required: true

permissions:
  contents: write
  pull-requests: write

env:
  PYTHONHASHSEED: "0"

jobs:
  regen:
    runs-on: ubuntu-latest
    timeout-minutes: 45
    steps:
      - uses: actions/checkout@v7
      - uses: astral-sh/setup-uv@v10
        with: { version: "0.12.6", python-version: "3.11", enable-cache: true, cache-dependency-glob: "uv.lock" }
      - run: uv sync --locked --all-packages --all-extras --dev
      - run: |
          rm -rf tests/goldens/__snapshots__
          uv run pytest tests/goldens --snapshot-update -q
          uv run python benchmarks/run_panel.py --tier full --update-baseline
      - uses: peter-evans/create-pull-request@v7
        with:
          branch: chore/regen-goldens
          title: "chore: regenerate goldens and benchmark baseline"
          labels: approved:algorithm-change
          body: |
            Regenerated by workflow_dispatch.

            **Reason:** ${{ inputs.reason }}

            Review the snapshot and baseline diffs below as a *scientific* change.
            Every sequence here was validated by `verify_solution()` before being
            written, so nothing biologically invalid can be laundered in.
          commit-message: "chore: regenerate goldens and benchmark baseline"
```

---

## 4. Other `.github` files

### `.github/CODEOWNERS`

```
# Routing + documented module ownership boundaries. The RULE
# "Require review from Code Owners" is deliberately NOT enabled --
# code owners need write access, so you are the only possible owner,
# you are also the author, and authors cannot self-approve. Enabling it
# would guarantee an unmergeable PR for every file.

*                               @masonberger4

# Load-bearing paths. Changes here also trip the path-guard label check.
/.github/                       @masonberger4
/packages/engine/src/bt5/verify.py  @masonberger4
/tests/invariants/              @masonberger4
/benchmarks/baseline.json       @masonberger4
/benchmarks/tolerances.yaml     @masonberger4
/data/                          @masonberger4

# Module boundaries -- point each agent session at its own directory.
/packages/engine/               @masonberger4
/packages/server/               @masonberger4
/apps/web/                      @masonberger4
```

### `.github/scripts/path_guard.py`

```python
#!/usr/bin/env python3
"""Load-bearing paths require an explicit owner label.

There is no second reviewer to require, so the enforceable substitute is:
touching the correctness spine, the goldens, the reference data, or CI
itself requires a label the owner applies after reading the diff. The label
appears in the PR timeline as an explicit human act -- which is exactly the
reviewable checkpoint a solo-owner repo otherwise lacks.
"""

from __future__ import annotations
import fnmatch, json, os, subprocess, sys

RULES: dict[str, str] = {
    "packages/engine/src/bt5/verify.py": "approved:oracle-change",
    "tests/invariants/*": "approved:oracle-change",
    "tests/goldens/__snapshots__/*": "approved:algorithm-change",
    "benchmarks/baseline.json": "approved:algorithm-change",
    "benchmarks/tolerances.yaml": "approved:algorithm-change",
    "benchmarks/panel.json": "approved:algorithm-change",
    "data/*": "approved:data-change",
    ".github/*": "approved:ci-change",
}


def main() -> int:
    labels = set(json.loads(os.environ["LABELS"]))
    changed = subprocess.check_output(
        ["git", "diff", "--name-only", os.environ["BASE"], os.environ["HEAD"]],
        text=True,
    ).split()

    needed: set[str] = set()
    for pattern, label in RULES.items():
        prefix = pattern.rstrip("*")
        for f in changed:
            if fnmatch.fnmatch(f, pattern) or (pattern.endswith("*") and f.startswith(prefix)):
                needed.add(label)

    missing = sorted(needed - labels)
    if missing:
        print(
            f"::error::This PR touches load-bearing paths and needs label(s): "
            f"{missing}. The repository owner applies these after reviewing the diff."
        )
        for pattern, label in RULES.items():
            if label in missing:
                print(f"  {label}  <- {pattern}")
        return 1
    print("path-guard: ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

### `.github/scripts/assert_gate_complete.py`

Closes the "agent adds a job and forgets to wire it into the gate" hole.

```python
#!/usr/bin/env python3
"""Every job in ci.yml must be reachable from `required-checks.needs`.

An agent that adds a job without adding it to needs: creates a job that can
fail while the gate never looks at it -- the PR merges green. This runs
nightly and on any .github change.
"""

from __future__ import annotations
import sys, yaml, pathlib

wf = yaml.safe_load(pathlib.Path(".github/workflows/ci.yml").read_text())
jobs = set(wf["jobs"]) - {"required-checks"}
needs = set(wf["jobs"]["required-checks"]["needs"])

orphans = sorted(jobs - needs)
if orphans:
    print(f"::error::ci.yml jobs not enforced by required-checks: {orphans}")
    sys.exit(1)

stale = sorted(needs - jobs)
if stale:
    print(f"::error::required-checks needs nonexistent jobs: {stale}")
    sys.exit(1)

gate = wf["jobs"]["required-checks"]
if "always()" not in str(gate.get("if", "")):
    print("::error::required-checks must use `if: always()`")
    sys.exit(1)
if not any("needs.*.result" in str(s.get("if", "")) for s in gate["steps"]):
    print("::error::required-checks has no needs.*.result inspection step")
    sys.exit(1)

print(f"gate integrity: ok ({len(jobs)} jobs enforced)")
```

### `.github/pull_request_template.md`

```markdown
## What changed

<!-- One paragraph. Which module does this PR own? -->

**Module owned by this PR:** `packages/engine` | `packages/server` | `apps/web` | `infra`

## Scientific impact

- [ ] No change to optimizer output for any golden case
- [ ] Output changes — golden/baseline diff justified below

<!-- If output changed, explain the SCIENCE, not the code. What is now more
     correct, and how do you know? Attach the benchmark diff table. -->

## Interface contract

- [ ] No shared Pydantic model, FastAPI route signature, or exported TS type changed
- [ ] Shared interface changed — listed below, with the consumers I checked

## Checklist

- [ ] `uv lock` / `pnpm install --lockfile-only` regenerated if deps changed
- [ ] New logic has tests that assert on the **biology**, not the implementation
- [ ] Any RNG is explicitly seeded (`np.random.default_rng(seed)`)
- [ ] No blanket `# type: ignore`, `# noqa`, `# pragma: no cover`, `// eslint-disable`
- [ ] If this is a bugfix: a fixture exists under `tests/data/regressions/issue_NNNN_*`
      that **fails on the parent commit**
- [ ] I did not modify `.github/**`, `verify.py`, `tests/invariants/**`,
      `data/**`, `benchmarks/baseline.json`, or `**/__snapshots__/**`
      (these require an owner label)

## Owner sign-off

<!-- Owner: after reading the diff, comment `/approve <head-sha>` -->
```

### `.github/ISSUE_TEMPLATE/bug_report.yml`

```yaml
name: Bug report
description: Something is wrong with the software
labels: ["bug"]
body:
  - type: input
    id: version
    attributes: { label: BT5 version / commit SHA }
    validations: { required: true }
  - type: dropdown
    id: area
    attributes:
      label: Area
      options: [engine (optimizer), engine (metrics), server (FastAPI), web (UI), packaging]
    validations: { required: true }
  - type: textarea
    id: repro
    attributes:
      label: Minimal reproduction
      description: Protein sequence (or accession), constraint set, and seed.
      render: json
    validations: { required: true }
  - type: textarea
    id: expected
    attributes: { label: Expected vs actual }
    validations: { required: true }
  - type: input
    id: platform
    attributes: { label: OS + Python version }
```

### `.github/ISSUE_TEMPLATE/scientific_correctness.yml`

```yaml
name: Scientific correctness
description: The output is syntactically fine but biologically wrong
labels: ["bug", "scientific-correctness", "priority:high"]
body:
  - type: markdown
    attributes:
      value: |
        This is the highest-severity class of bug in BT5. A wrong sequence
        that looks right ships silently. Please give us enough to build a
        permanent regression fixture.
  - type: textarea
    id: protein
    attributes: { label: Input protein (FASTA or UniProt accession) }
    validations: { required: true }
  - type: textarea
    id: constraints
    attributes: { label: Constraints JSON (including seed and table_id), render: json }
    validations: { required: true }
  - type: textarea
    id: output
    attributes: { label: Produced DNA }
    validations: { required: true }
  - type: dropdown
    id: violation
    attributes:
      label: Which invariant is violated?
      options:
        - Back-translation does not round-trip
        - Premature or missing stop codon
        - Forbidden motif present (plus strand)
        - Forbidden motif present (minus strand)
        - Forbidden motif created across a codon junction
        - GC constraint violated (global)
        - GC constraint violated (window)
        - Non-deterministic output for a fixed seed
        - Wrong codon table for the declared organism
        - Silent best-effort instead of InfeasibleConstraints
        - Other
    validations: { required: true }
  - type: textarea
    id: evidence
    attributes:
      label: Independent verification
      description: Output of an independent check (Biopython translate, DNA Chisel, benchling…).
```

### `.github/ISSUE_TEMPLATE/config.yml`

```yaml
blank_issues_enabled: true
contact_links:
  - name: Security vulnerability
    url: https://github.com/masonberger4/BT5/security/advisories/new
    about: Report privately via a security advisory, not a public issue.
```

### `.github/dependabot.yml`

Scope Dependabot to **github-actions only**. GitHub's supported-ecosystem list does not include `uv`, so Dependabot cannot regenerate `uv.lock` — and with `uv sync --locked` in CI, a PR that bumps `pyproject.toml` without the lock fails immediately. Renovate owns the dependency updates.

```yaml
version: 2
updates:
  - package-ecosystem: "github-actions"
    directory: "/"
    schedule:
      interval: "weekly"
      day: "monday"
      time: "06:00"
    open-pull-requests-limit: 2
    labels: ["dependencies", "ci", "approved:ci-change"]
    groups:
      github-actions:
        patterns: ["*"]
```

### `renovate.json`

```json
{
  "$schema": "https://docs.renovatebot.com/renovate-schema.json",
  "extends": ["config:recommended", ":semanticCommits"],
  "schedule": ["* 0-6 * * 1"],
  "prConcurrentLimit": 3,
  "prHourlyLimit": 2,
  "rebaseWhen": "conflicted",
  "lockFileMaintenance": { "enabled": true, "schedule": ["* 0-6 * * 1"] },
  "packageRules": [
    {
      "matchManagers": ["pep621"],
      "matchDepTypes": ["dev", "dependency-groups"],
      "groupName": "python dev deps",
      "automerge": true
    },
    {
      "matchManagers": ["npm"],
      "matchDepTypes": ["devDependencies"],
      "groupName": "web dev deps",
      "automerge": true
    },
    {
      "description": "ViennaRNA: no cp312/cp313 wheels, and the energy parameter set moves dG. Bump by hand with a baseline regeneration.",
      "matchPackageNames": ["ViennaRNA", "viennarna"],
      "enabled": false
    },
    {
      "description": "numpy/biopython changes can move numeric output.",
      "matchPackageNames": ["numpy", "biopython"],
      "automerge": false,
      "addLabels": ["ci-full", "needs-baseline-check"]
    },
    {
      "matchUpdateTypes": ["major"],
      "automerge": false,
      "addLabels": ["major-update", "ci-full"]
    },
    {
      "matchPackageNames": ["@playwright/test", "playwright"],
      "groupName": "playwright",
      "addLabels": ["ci-full"]
    }
  ],
  "vulnerabilityAlerts": { "labels": ["security"], "automerge": true }
}
```

### `CLAUDE.md` — the contract every AI session must follow

```markdown
# BT5 — contract for AI coding sessions

You are one of several agents working on this repository in parallel. The
repository owner is the only human and cannot review everything closely.
The rules below are what make parallel work safe. Follow them exactly.

## 0. What this project is

BT5 is a locally-run codon-optimization / protein back-translation desktop
app for working scientists. A wrong answer that looks right is the worst
possible outcome — worse than a crash, worse than an ugly UI, worse than
missing a deadline. Optimize for being right.

## 1. Module ownership — stay in your lane

    packages/engine/   Python 3.11 optimizer, codon tables, metrics
    packages/server/   FastAPI local server
    apps/web/          React + TypeScript UI
    tests/             invariants, goldens, data_integrity, differential, e2e
    benchmarks/        metric regression panel
    data/              codon tables, motif DB

Your session owns ONE of these. Do not edit files in another module. If your
change requires a change in another module, STOP and open an issue describing
the interface change instead of making it. Cross-module PRs are the main
source of merge conflicts and semantic breakage between parallel agents.

## 2. Files you must NEVER modify

    packages/engine/src/bt5/verify.py     the correctness oracle
    tests/invariants/**                   the properties the oracle is asserted through
    tests/goldens/**/__snapshots__/**     pinned optimizer output
    benchmarks/baseline.json              pinned metric baseline
    benchmarks/tolerances.yaml            regression thresholds
    data/**                               codon tables, motif databases
    .github/**                            CI, rulesets, this contract's enforcement

Every one of these is label-gated: CI fails your PR unless the owner has
applied the matching `approved:*` label. Do not ask for the label as a
shortcut. If you believe one of these files is wrong, open an issue with your
evidence and let the owner decide.

**Never run `pytest --snapshot-update`.** Never regenerate
`benchmarks/baseline.json`. If goldens need regenerating, that is the owner
running the `Regenerate goldens` workflow.

## 3. Correctness rules

- `optimize()` calls `verify_solution()` as a post-condition, by default,
  always. Do not add a code path that bypasses it. `_verify=False` exists only
  for benchmarking internals and must never be reachable from the public API
  or the FastAPI surface.
- The optimizer is a **trichotomy**: it either returns a fully compliant
  sequence or raises `InfeasibleConstraints`. There is no third option.
  Returning a "best effort" sequence with a logged warning is a P0 bug, not a
  feature. If you find yourself writing `logger.warning("could not satisfy")`,
  raise instead.
- Never guess at biology. Reject non-canonical residues with `ValueError`;
  do not silently map `X` to anything.
- Forbidden motifs must be checked on **both strands** and across **codon
  junctions**. These are the two bug classes agents produce most reliably.
- Any randomness uses `np.random.default_rng(seed)` or `random.Random(seed)`.
  `np.random.seed(...)` and `random.seed(...)` are banned by CI — global RNG
  state destroys reproducibility across parallel modules.
- Any hard-coded numeric constant needs a source citation in a comment.

## 4. Suppression is not a fix

Banned by CI: blanket `# type: ignore`, blanket `# noqa`,
`# pragma: no cover`, `// eslint-disable` without a rule name.
Specific, justified forms (`# type: ignore[arg-type]  # ViennaRNA has no stubs`)
are fine. If a check is failing, fix the code. If you genuinely believe the
check is wrong, say so in the PR body — do not silence it.

## 5. Dependencies and lockfiles

- Python: edit `pyproject.toml` then run `uv lock`. CI runs
  `uv sync --locked` and will fail on a stale lock.
- Web: `pnpm install` then commit `pnpm-lock.yaml`. CI runs
  `--frozen-lockfile`.
- **Never hand-resolve a lockfile merge conflict.** Rebase on `main` and
  regenerate (`uv lock`, `pnpm install --lockfile-only`). A three-way merge of
  a lockfile is garbage.
- Dependency changes go in their own small PR, merged first.

## 6. Bugfix PRs

Every bugfix adds a fixture under `tests/data/regressions/issue_NNNN_<slug>/`
containing the protein, the constraints, and a README explaining the bug.
The fixture must **fail on the parent commit**. State in the PR body that you
verified this. A bug without a regression fixture will ship again.

## 7. Branching and PRs

- Branch: `agent/<module>/<short-slug>` (e.g. `agent/engine/junction-motifs`).
- Branch from the latest `main`. Do not branch from another agent's branch.
- One module, one concern, one PR. Squash merge is the only merge method.
- Fill in the PR template honestly, especially the "Scientific impact" and
  "Interface contract" sections.
- The PR title becomes the commit message on `main`. Write it for a human
  reading `git log` in six months.

## 8. CI

- Local before pushing: `uv run ruff check . && uv run ruff format --check .
  && uv run mypy . && uv run pytest -m "not slow" -q`
- `required-checks` is an aggregator. If you add a job to `.github/workflows/ci.yml`
  you must add it to `required-checks.needs` — but you are not allowed to edit
  `.github/**` anyway, so this is a note about why, not permission.
- If a required check is red, the PR does not merge. Not for the owner either.
  There is no bypass. Fix it.
- If CI fails for a reason unrelated to your diff (a Hypothesis counterexample
  in a module you did not touch, a flaky runner), do NOT weaken the test. Say so
  in the PR, file an issue with the counterexample, and let the owner triage.

## 9. Never

- Add a `paths:` filter to `on:` in any workflow (it permanently deadlocks a
  required check).
- Add a bypass actor to any ruleset.
- Enable "Allow GitHub Actions to create and approve pull requests".
- Commit a secret, an API key, or a `.env`. Push protection will block you and
  the incident is yours.
- Use `pull_request_target` in any workflow that checks out PR head code.
- Weaken an assertion, a threshold, or a tolerance to make CI pass.
```

---

## 5. The merge gate table

| Required check | Workflow | What it proves | Typical wall clock | Phase |
|---|---|---|---|---|
| **`required-checks`** | `ci.yml` | Aggregate: nothing in the CI workflow failed or was cancelled, and the path-filter job itself succeeded. The **only** name from `ci.yml` in the ruleset, so adding/renaming/matrix-ing jobs never touches repo settings. | +10s over its slowest dependency | 0 |
| ↳ `changes` | `ci.yml` | Which modules the PR touches (fails open), **and** that load-bearing paths carry their owner label. | ~20s | 0 |
| ↳ `python-quality` | `ci.yml` | ruff lint + format, mypy `--strict`, no blanket suppressions, no global RNG, `optimize()` verifies by default. | ~60s warm | 0 |
| ↳ `python-tests` | `ci.yml` | Unit + data-integrity + golden tests, ≥85% coverage, native wheels importable. On PRs: ubuntu only. | ~90s | 0 |
| ↳ `invariants` | `ci.yml` | Hypothesis properties: round-trip-or-raise trichotomy, determinism under a fixed seed, junction motifs, constraint monotonicity. **The single most valuable gate.** | ~90s (200 examples) | 0 |
| ↳ `goldens-not-hand-edited` | `ci.yml` | Committed snapshots match freshly regenerated output — i.e. nobody hand-edited a snapshot file to fake green. | ~60s | 0 |
| ↳ `benchmark-gate` | `ci.yml` | CAI, held-out CAI, GC band, GC windows, forbidden hits (hard 0), homopolymer, 5′ ΔG did not regress beyond direction-aware, noise-aware tolerance on the fast panel. | ~4–8 min | 0 |
| ↳ `web` | `ci.yml` | `tsc --noEmit`, eslint `--max-warnings=0`, prettier, vitest + coverage thresholds. | ~90s | 0 |
| ↳ `e2e` | `ci.yml` | Playwright chromium against the real FastAPI server. | ~2.5 min | 1 |
| **`codeql-passed`** | `codeql.yml` | Aggregate SAST over python, javascript-typescript, and **`actions`** (scans your own workflows for injection and over-permissioned tokens), `security-extended`. | ~4–7 min | 0 |
| **`dependency-review`** | `dependency-review.yml` | No new dependency carries a ≥moderate advisory. The one third-party gate that **fails by default** instead of commenting. | ~20s | 0 |
| **`agent-approval-check`** | `agent-approval-check.yml` | A human deliberately vouched for the current head SHA (`/approve <sha>`). Re-arms on every push. Fails closed. **This is the review requirement, re-implemented as a status check because GitHub forbids self-approval.** | ~20s | 2 |
| **`claude-review-gate`** | `claude-review.yml` | An independent Claude pass found no confirmed scientific-correctness, crash, security, or shared-interface defect — enforced via `--json-schema` → `structured_output` → `exit 1`, failing closed on an incomplete review. | ~3–6 min | 3 |

**Full-suite wall clock on a Python PR:** roughly 8–10 minutes, dominated by `benchmark-gate` and `codeql-passed`. On a web-only PR: ~5 minutes (Python jobs skip and report `skipped`, which satisfies the gate).

**Not required, deliberately:**

- **CodeRabbit** — free forever on public repos, genuinely different blind spots from Claude. Install it, leave it **comment-only**. If the double-commenting becomes noise, drop CodeRabbit, never Claude — Claude holds the veto.
- **Copilot code review** — cannot block, period. Skipped.
- **Nightly jobs** (deep Hypothesis, mutation score, DNA Chisel differential, OS matrix, lock drift, audit, gate integrity) — they open issues, they don't block PRs.

---

## 6. Scientific-correctness gates

The organizing principle: **the correctness oracle lives in production code, not in tests.** If verification lives in `tests/`, the agent that writes a broken optimizer also writes the test and the bug ships. If it lives in `src/` and runs on every call, a wrong sequence raises in the user's app, in every property example, in every golden test, and in every benchmark row — from *one* definition that is label-gated against modification.

### `packages/engine/src/bt5/verify.py` — the oracle

```python
"""The correctness oracle. PRODUCTION code, called by optimize() on every run.

This module independently RE-DERIVES the translation using Biopython's NCBI
tables rather than reusing the optimizer's own codon table -- so a corrupted
bundled table breaks the round trip instead of being consistently wrong on
both sides.

Changing this file requires the `approved:oracle-change` label.
"""

from __future__ import annotations

import re
from Bio.Data import CodonTable
from Bio.Seq import Seq

from bt5.constraints import Constraints
from bt5.seqio import gc_fraction, iupac_regex, revcomp

DNA = frozenset("ACGT")


class VerificationError(AssertionError):
    """The optimizer produced a sequence that does not encode the input."""


class InfeasibleConstraints(Exception):
    """No sequence satisfies the given constraints. The ONLY legal failure."""


def verify_solution(protein: str, cons: Constraints, dna: str, *, table_id: int = 1) -> None:
    problems: list[str] = []

    # I1 -- alphabet
    bad = sorted(set(dna) - DNA)
    if bad:
        problems.append(f"non-ACGT characters: {bad}")

    # I2 -- reading frame. Checked BEFORE translating: Biopython's translate()
    # on a non-multiple-of-3 emits a BiopythonWarning and silently truncates,
    # so a frame-length bug can pass a naive round-trip test.
    want = 3 * len(protein)
    if len(dna) != want:
        problems.append(f"len {len(dna)} != 3*len(protein) = {want}")

    # I3 -- THE round trip, via an independent re-translation.
    if not problems:
        back = str(Seq(dna).translate(table=table_id, cds=False))
        if back != protein:
            i = next(
                (k for k, (a, b) in enumerate(zip(back, protein)) if a != b),
                min(len(back), len(protein)),
            )
            problems.append(
                f"back-translation mismatch at residue {i}: got {back[i : i + 1]!r} "
                f"want {protein[i : i + 1]!r} (codon {dna[i * 3 : i * 3 + 3]!r})"
            )

    # I4 -- initiator
    if protein.startswith("M") and cons.require_atg_start and dna[:3] != "ATG":
        problems.append(f"initiator codon {dna[:3]!r}, expected ATG")

    # I5 -- stops: terminal required if declared, interior always forbidden
    stops = set(CodonTable.unambiguous_dna_by_id[table_id].stop_codons)
    if protein.endswith("*") and dna[-3:] not in stops:
        problems.append(f"terminal codon {dna[-3:]!r} is not a stop")
    interior = [dna[i : i + 3] for i in range(0, max(0, len(dna) - 3), 3)]
    premature = [i for i, c in enumerate(interior) if c in stops]
    if premature:
        problems.append(f"premature stop at codon index {premature[:5]}")

    # I6 -- forbidden motifs on BOTH STRANDS. Junction-spanning and
    # minus-strand motifs are the two bug classes agents produce most often,
    # and both are invisible to any per-codon test.
    rc = revcomp(dna)
    for motif in cons.forbidden_motifs:
        rx = iupac_regex(motif)
        for strand, seq in (("+", dna), ("-", rc)):
            hits = [m.start() for m in re.finditer(f"(?={rx})", seq)]
            if hits:
                problems.append(f"forbidden motif {motif!r} on {strand} strand at {hits[:5]}")

    # I7 -- GC band, global and windowed, only when declared hard
    if cons.gc_is_hard:
        gc = gc_fraction(dna)
        if not (cons.gc_min <= gc <= cons.gc_max):
            problems.append(f"global GC {gc:.3f} outside [{cons.gc_min},{cons.gc_max}]")
        if cons.gc_window:
            w = cons.gc_window
            gcs = [1 if b in "GC" else 0 for b in dna]
            run = sum(gcs[:w])
            for i in range(len(dna) - w + 1):  # rolling sum, O(n) not O(n*w)
                if i:
                    run += gcs[i + w - 1] - gcs[i - 1]
                frac = run / w
                if not (cons.gc_min <= frac <= cons.gc_max):
                    problems.append(f"window GC {frac:.3f} at offset {i} outside band")
                    break

    # I8 -- homopolymer / repeat ceiling
    if cons.avoid_repeats_len:
        m = re.search(r"(A{%d,}|C{%d,}|G{%d,}|T{%d,})" % ((cons.avoid_repeats_len,) * 4), dna)
        if m:
            problems.append(f"homopolymer run of {len(m.group(1))} at offset {m.start()}")

    if problems:
        raise VerificationError(
            "invalid encoding of input protein:\n  - " + "\n  - ".join(problems)
        )
```

Wired in as a default-on post-condition:

```python
# packages/engine/src/bt5/optimizer.py
def optimize(protein: str, cons: Constraints, *, _verify: bool = True) -> str:
    """Return a DNA sequence encoding `protein` under `cons`.

    Raises InfeasibleConstraints if no compliant sequence exists.
    NEVER returns a partially-compliant "best effort" sequence.
    """
    dna = _search(protein, cons)
    if _verify:
        verify_solution(protein, cons, dna, table_id=cons.table_id)
    return dna
```

The `python-quality` job greps for `_verify: bool = True` so the default cannot be flipped, and `path_guard.py` requires `approved:oracle-change` to touch `verify.py` at all.

### `tests/strategies.py`

```python
from hypothesis import strategies as st
from bt5.constraints import Constraints

CANONICAL_AA = "ACDEFGHIKLMNPQRSTVWY"  # no B/Z/J/X, no U/O, no '*'
SINGLE_CODON = "MW"  # zero degrees of freedom
MOTIF_POOL = [
    "GAATTC",
    "GGATCC",
    "AAGCTT",
    "CTGCAG",
    "GTCGAC",
    "TCTAGA",
    "GCGGCCGC",
    "AAAAAAAA",
    "GGGGGGGG",
    "TTTTTTTT",
]


@st.composite
def proteins(draw, min_size=1, max_size=250, start_met=None, with_stop=None):
    core = draw(st.text(alphabet=CANONICAL_AA, min_size=min_size, max_size=max_size))
    if start_met is None:
        start_met = draw(st.booleans())
    if with_stop is None:
        with_stop = draw(st.booleans())
    return ("M" if start_met else "") + core + ("*" if with_stop else "")


@st.composite
def degenerate_proteins(draw):
    """Poly-M/W: the optimizer has no choices at all. Excellent shrink target."""
    return "M" + draw(st.text(alphabet=SINGLE_CODON, min_size=1, max_size=60))


@st.composite
def constraint_sets(draw, feasible_bias=False):
    lo = draw(st.integers(35, 45) if feasible_bias else st.integers(25, 50))
    hi = draw(st.integers(lo + 20, 70) if feasible_bias else st.integers(lo + 10, 75))
    return Constraints(
        forbidden_motifs=tuple(
            draw(
                st.lists(
                    st.sampled_from(MOTIF_POOL), max_size=2 if feasible_bias else 4, unique=True
                )
            )
        ),
        gc_min=lo / 100,
        gc_max=hi / 100,
        gc_window=draw(st.sampled_from([None, 50, 100])),
        gc_is_hard=draw(st.booleans()),
        avoid_repeats_len=draw(st.sampled_from([None, 8, 12, 15])),
        organism=draw(st.sampled_from(["e_coli_k12", "h_sapiens", "s_cerevisiae", "cho"])),
        seed=draw(st.integers(0, 2**32 - 1)),
    )
```

### `tests/conftest.py`

```python
import os
from hypothesis import HealthCheck, Phase, settings

settings.register_profile("dev", max_examples=25, deadline=None)
settings.register_profile(
    "ci", max_examples=200, deadline=None, suppress_health_check=[HealthCheck.too_slow]
)
settings.register_profile(
    "deep",
    max_examples=5000,
    deadline=None,
    phases=list(Phase),
    suppress_health_check=[HealthCheck.too_slow],
)
settings.load_profile(os.getenv("HYPOTHESIS_PROFILE", "dev"))

# Do NOT set derandomize=True on PRs. For correctness invariants any failure
# is a real bug, never a flake, and derandomizing throws away exploration.
```

### `tests/invariants/test_round_trip.py` — the trichotomy

```python
import statistics
import pytest
from hypothesis import HealthCheck, assume, given, settings, target

from bt5 import optimize
from bt5.verify import InfeasibleConstraints, verify_solution
from tests.strategies import constraint_sets, degenerate_proteins, proteins
from tests.helpers import default_constraints, near_miss_count


@given(protein=proteins(), cons=constraint_sets())
@settings(
    max_examples=200,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)
def test_round_trip_or_declared_infeasible(protein, cons):
    """THE property. Either raise, or be completely correct.

    Silent 'best effort' partial compliance is the #1 AI failure mode in
    constraint-solving code and must be structurally impossible.
    """
    try:
        dna = optimize(protein, cons)
    except InfeasibleConstraints:
        return  # declaring failure is always allowed
    verify_solution(protein, cons, dna, table_id=cons.table_id)


@given(protein=degenerate_proteins(), cons=constraint_sets())
@settings(max_examples=100, deadline=None)
def test_zero_degrees_of_freedom(protein, cons):
    """Poly-M/W: M and W have exactly one codon each. The optimizer has no
    choices, so it must either produce that unique sequence or declare
    infeasibility -- never invent a codon."""
    try:
        dna = optimize(protein, cons)
    except InfeasibleConstraints:
        return
    verify_solution(protein, cons, dna, table_id=cons.table_id)


@given(protein=proteins(), cons=constraint_sets())
@settings(max_examples=100, deadline=None)
def test_determinism_under_fixed_seed(protein, cons):
    try:
        a = optimize(protein, cons)
    except InfeasibleConstraints:
        return
    assert optimize(protein, cons) == a, "optimizer is not deterministic for a fixed seed"


def test_seed_actually_varies_output_for_stochastic_strategies():
    """Anti-cheat: an optimizer that ignores the seed passes the determinism
    test trivially. Only meaningful for stochastic strategies."""
    from bt5.constraints import Constraints

    p = read_fasta("tests/data/proteins/P42212.faa")
    outs = {
        optimize(p, Constraints(strategy="sampled", organism="e_coli_k12", seed=s))
        for s in range(8)
    }
    assert len(outs) > 1, "seed has no effect; is it wired through?"


@given(protein=proteins(), cons=constraint_sets())
@settings(max_examples=100, deadline=None)
def test_tightening_constraints_never_breaks_looser_ones(protein, cons):
    """A solution to a superset of constraints must satisfy the subset."""
    loose = cons.without_motifs()
    try:
        tight = optimize(protein, cons)
    except InfeasibleConstraints:
        return
    verify_solution(protein, loose, tight, table_id=cons.table_id)


@given(protein=proteins(min_size=30), cons=constraint_sets())
@settings(max_examples=200, deadline=None)
def test_junction_motifs(protein, cons):
    """Motifs created ACROSS codon boundaries are the #1 real bug class.
    target() steers Hypothesis toward Hamming-distance-1 near misses."""
    assume(cons.forbidden_motifs)
    try:
        dna = optimize(protein, cons)
    except InfeasibleConstraints:
        return
    target(float(near_miss_count(dna, cons.forbidden_motifs)))
    verify_solution(protein, cons, dna, table_id=cons.table_id)


@pytest.mark.parametrize("ch", list("BZJXUO") + ["m", " ", "-", "1", "*"])
def test_non_canonical_residues_rejected_not_guessed(ch):
    with pytest.raises(ValueError):
        optimize("MA" + ch + "KL", default_constraints())


@pytest.mark.parametrize("seeds", [[11, 22, 33, 44, 55, 66, 77]])
def test_cai_beats_random_backtranslation_on_every_seed(seeds):
    """Invariant-strength on every seed; distributional on the mean.
    Asserting a quality claim on ONE seed produces a test that fails one run
    in twenty and gets deleted by the next agent."""
    from bt5.metrics import cai
    from bt5.constraints import Constraints

    p = read_fasta("tests/data/proteins/P42212.faa")
    naive = [cai(random_backtranslate(p, seed=s), "e_coli_k12") for s in seeds]
    opt = [
        cai(optimize(p, Constraints(organism="e_coli_k12", seed=s)), "e_coli_k12") for s in seeds
    ]
    assert all(o > n for o, n in zip(opt, naive)), list(zip(opt, naive))
    assert statistics.fmean(opt) >= 0.80
    assert statistics.pstdev(opt) < 0.05, f"seed spread {statistics.pstdev(opt):.3f}"
```

### `tests/goldens/test_golden_sequences.py`

```python
import json
from pathlib import Path

import pytest

from bt5 import optimize
from bt5.constraints import Constraints
from bt5.metrics import deterministic_metrics
from bt5.verify import verify_solution
from tests.helpers import read_fasta

CASES = json.loads((Path(__file__).parent / "cases.json").read_text())


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_golden_sequence(case, snapshot):
    protein = read_fasta(Path(__file__).parent / "proteins" / case["fasta"])
    cons = Constraints(**case["constraints"], seed=case["seed"])
    dna = optimize(protein, cons)

    # 1. ORACLE FIRST. A snapshot can be wrong; the oracle is label-gated.
    #    This is what stops `--snapshot-update` from laundering a
    #    biologically invalid sequence into a "passing" baseline.
    verify_solution(protein, cons, dna, table_id=cons.table_id)

    # 2. Byte-exact pin of the optimizer's DECISION.
    assert dna == snapshot

    # 3. Deterministic metrics only. NO folding energy here -- ViennaRNA
    #    floats are parameter-set and platform sensitive; dG lives in the
    #    benchmark gate with a tolerance.
    assert {k: round(v, 4) for k, v in deterministic_metrics(dna).items()} == snapshot(
        name="metrics"
    )
```

Trustworthiness comes from the paired CI job `goldens-not-hand-edited` (§3), which deletes `__snapshots__`, regenerates from scratch, and diffs — catching the *other* way agents fake green.

### `tests/data_integrity/test_codon_tables.py`

Semantic validation beats checksums, because an agent can regenerate a checksum in the same commit but cannot regenerate biology.

```python
import hashlib
import json
from pathlib import Path

import pytest
from Bio.Data import CodonTable

from bt5.tables import ALL_ORGANISMS, load_codon_table


def test_manifest_matches_files():
    for line in Path("data/MANIFEST.sha256").read_text().splitlines():
        digest, rel = line.split("  ", 1)
        assert hashlib.sha256(Path("data", rel).read_bytes()).hexdigest() == digest, rel


@pytest.mark.parametrize("org", ALL_ORGANISMS)
def test_codon_table_is_structurally_valid(org):
    t = load_codon_table(org)
    ncbi = CodonTable.unambiguous_dna_by_id[t.table_id]

    codons = [c.codon for aa in t.by_aa.values() for c in aa]
    assert len(codons) == 64 and len(set(codons)) == 64, "not a partition of 64 codons"

    for aa, entries in t.by_aa.items():
        assert entries, f"{aa} has no codons"
        total = sum(e.frequency for e in entries)
        assert abs(total - 1.0) < 1e-6, f"{aa} frequencies sum to {total}"
        assert all(e.frequency >= 0 for e in entries)

    # THE decisive check. Our table must AGREE with NCBI on what each codon
    # means. This cannot be regenerated to match a corrupted table.
    for aa, entries in t.by_aa.items():
        for e in entries:
            expected = "*" if e.codon in ncbi.stop_codons else ncbi.forward_table[e.codon]
            assert expected == aa, (
                f"{e.codon} maps to {aa}, NCBI table {t.table_id} says {expected}"
            )

    assert set(t.by_aa) == set("ACDEFGHIKLMNPQRSTVWY*")


def test_motif_db_entries_are_valid_iupac_and_have_provenance():
    db = json.loads(Path("data/motifs.json").read_text())
    for name, rec in db.items():
        assert set(rec["pattern"]) <= set("ACGTRYSWKMBDHVN"), name
        assert rec["source"].startswith("https://"), f"{name} lacks provenance"
        assert rec.get("both_strands") in (True, False), name
```

⚠️ Semantic validation only covers what you thought to encode. It will *not* catch an agent swapping E. coli K-12 frequencies for B-strain frequencies — both are structurally valid. That is covered by the CAI-reference-set version pin in the benchmark baseline, which makes any usage-table swap fail the comparability check.

### `benchmarks/run_panel.py` — the metric regression gate

Invariants prove the sequence is *legal*; they say nothing about whether it is *good*. An agent refactoring the search can keep every invariant green while silently dropping CAI from 0.86 to 0.61.

```python
#!/usr/bin/env python
"""Metric regression gate. Exits non-zero on regression."""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
import time
from pathlib import Path

import yaml

from bt5 import optimize
from bt5.constraints import Constraints
from bt5.metrics import (
    cai,
    cai_reference_version,
    count_motifs,
    five_prime_dg,
    gc_fraction,
    gc_window_extremes,
    longest_homopolymer,
    mfe_total,
    vienna_version,
)
from bt5.verify import verify_solution
from tests.helpers import read_fasta

ROOT = Path(__file__).resolve().parent
PANEL = json.loads((ROOT / "panel.json").read_text())
TOL = yaml.safe_load((ROOT / "tolerances.yaml").read_text())["metrics"]
SEEDS = [11, 22, 33, 44, 55]  # FROZEN. Changing them IS an algorithm change.


def measure_once(protein: str, cons: Constraints) -> dict[str, float]:
    t0 = time.perf_counter()
    dna = optimize(protein, cons)
    rt = time.perf_counter() - t0
    verify_solution(protein, cons, dna, table_id=cons.table_id)  # gate the gate
    lo, hi = gc_window_extremes(dna, window=cons.gc_window or 50)
    return {
        "cai": cai(dna, cons.organism),
        "cai_heldout": cai(dna, cons.organism, reference="heldout"),
        "gc_percent": gc_fraction(dna),
        "gc_window_min": lo,
        "gc_window_max": hi,
        "forbidden_hits": count_motifs(dna, cons.forbidden_motifs, both_strands=True),
        "homopolymer_max": longest_homopolymer(dna),
        "five_prime_dg": five_prime_dg(dna, window=45),
        "mfe_total": mfe_total(dna),
        "runtime_s": rt,
    }


def measure(entry: dict) -> dict:
    protein = read_fasta(ROOT / entry["fasta"])
    runs = [
        measure_once(
            protein,
            Constraints(
                **entry["constraints"],
                organism=entry["organism"],
                table_id=entry["table_id"],
                seed=s,
            ),
        )
        for s in SEEDS
    ]
    out = {}
    for k in runs[0]:
        vals = [r[k] for r in runs]
        out[k] = {
            "mean": statistics.fmean(vals),
            "stdev": statistics.pstdev(vals) if len(vals) > 1 else 0.0,
            "min": min(vals),
            "max": max(vals),
        }
    # runtime uses the median -- one slow runner must not move it
    out["runtime_s"]["mean"] = statistics.median(r["runtime_s"] for r in runs)
    # hard per-protein wall-clock ceiling; runner variance dwarfs real deltas,
    # so this is the runtime gate, not a percentage from baseline.
    if out["runtime_s"]["mean"] > entry["budget_s"]:
        out["_budget_exceeded"] = True
    return out


def compare(base: dict, cur: dict):
    rows, failures = [], []
    for cid, c in cur["cases"].items():
        if c.get("_budget_exceeded"):
            failures.append(f"{cid}: exceeded wall-clock budget")
        b = base["cases"].get(cid)
        if b is None:
            rows.append((cid, "-", "", "", "", "NEW", "no baseline"))
            continue
        for metric, spec in TOL.items():
            bm, cm = b[metric]["mean"], c[metric]["mean"]
            noise = max(b[metric]["stdev"], c[metric]["stdev"])
            delta = cm - bm
            verdict, why = "ok", ""
            if spec["direction"] == "exact":
                if cm != spec["value"]:
                    verdict, why = "FAIL", f"must equal {spec['value']}"
            else:
                # A move in the GOOD direction is never a failure.
                sign = -1 if spec["direction"] in ("higher_is_better",) else 1
                regression = sign * delta
                allowed = spec.get("abs")
                if allowed is None:
                    allowed = abs(bm) * spec["rel"] + spec.get("floor", 0.0)
                allowed = max(allowed, 2.0 * noise)  # never fail inside noise
                if regression > allowed:
                    verdict = "FAIL" if spec.get("blocking", True) else "warn"
                    why = f"regressed {regression:.4g} > allowed {allowed:.4g}"
            rows.append((cid, metric, f"{bm:.4g}", f"{cm:.4g}", f"{delta:+.4g}", verdict, why))
            if verdict == "FAIL":
                failures.append(f"{cid}.{metric}: {why}")
    return rows, failures


def markdown(rows, failures, verbose: bool) -> str:
    head = "## Benchmark panel\n\n" + (
        "All metrics within tolerance.\n\n"
        if not failures
        else f"**{len(failures)} regression(s) beyond tolerance.**\n\n"
    )
    tbl = [
        "| case | metric | baseline | current | delta | verdict | note |",
        "|---|---|--:|--:|--:|---|---|",
    ]
    for r in rows:
        if r[5] == "ok" and not verbose:
            continue  # surface movement only
        tbl.append("| " + " | ".join(str(x) for x in r) + " |")
    return head + "\n".join(tbl) + "\n"


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", default="fast", choices=["fast", "full"])
    ap.add_argument("--out", default="benchmarks/results.json")
    ap.add_argument("--compare", default="benchmarks/baseline.json")
    ap.add_argument("--summary")
    ap.add_argument("--update-baseline", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()

    entries = [e for e in PANEL if a.tier == "full" or e["tier"] == "fast"]
    sha = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    cur = {
        "generated_from_sha": sha,
        "seeds": SEEDS,
        "cai_reference_set_version": cai_reference_version(),
        "viennarna_version": vienna_version(),
        "cases": {e["id"]: measure(e) for e in entries},
    }
    Path(a.out).write_text(json.dumps(cur, indent=2, sort_keys=True) + "\n")

    if a.update_baseline:
        Path(a.compare).write_text(json.dumps(cur, indent=2, sort_keys=True) + "\n")
        sys.exit(0)

    base = json.loads(Path(a.compare).read_text())

    # Comparability guard: CAI is meaningless without a pinned reference set,
    # and ViennaRNA dG depends on the energy parameter version. A swap of
    # either makes every number shift and reads as an algorithmic "win".
    for k in ("cai_reference_set_version", "viennarna_version"):
        if base[k] != cur[k]:
            print(
                f"::error::{k} changed ({base[k]} -> {cur[k]}); metrics are not "
                "comparable. Regenerate the baseline in its own "
                "`approved:algorithm-change` PR."
            )
            sys.exit(1)

    rows, failures = compare(base, cur)
    md = markdown(rows, failures, a.verbose)
    print(md)
    if a.summary:
        with open(a.summary, "a") as fh:
            fh.write(md)
    sys.exit(1 if failures else 0)
```

### `benchmarks/tolerances.yaml`

```yaml
# direction: higher_is_better | lower_is_better | exact | band
# A move in the GOOD direction is never a failure.
# Every band is widened to at least 2*stdev across seeds, so the gate can
# never fail inside measurement noise.
metrics:
  cai:
    direction: higher_is_better
    abs: 0.02
    blocking: true

  cai_heldout:            # guards against fitting the usage table itself
    direction: higher_is_better
    abs: 0.03
    blocking: true

  gc_percent:
    direction: band
    band: [0.40, 0.60]
    abs: 0.01
    blocking: true

  gc_window_max:
    direction: lower_is_better
    abs: 0.03
    blocking: true

  gc_window_min:
    direction: higher_is_better
    abs: 0.03
    blocking: true

  forbidden_hits:         # hard gate, zero tolerance
    direction: exact
    value: 0
    blocking: true

  homopolymer_max:
    direction: lower_is_better
    abs: 0
    blocking: true

  five_prime_dg:          # kcal/mol; LESS negative = less structure = better
    direction: higher_is_better
    abs: 1.5
    blocking: true

  mfe_total:
    direction: higher_is_better
    abs: 30.0
    blocking: false       # whole-transcript MFE is noisy; informational

  runtime_s:              # runner variance dominates -> warn only.
    direction: lower_is_better   # Hard ceilings live in panel.json budget_s.
    rel: 0.30
    floor: 0.5
    blocking: false
```

### `benchmarks/panel.json` (excerpt)

```json
[
  {
    "id": "gfp_ecoli",
    "tier": "fast",
    "fasta": "../tests/data/proteins/P42212.faa",
    "source": "UniProt P42212 v2, retrieved 2026-08-01",
    "organism": "e_coli_k12",
    "table_id": 11,
    "budget_s": 20,
    "cai_floor": 0.80,
    "constraints": {
      "forbidden_motifs": ["GAATTC", "GGATCC", "AAGCTT", "GCGGCCGC"],
      "gc_min": 0.40, "gc_max": 0.60,
      "gc_window": 50, "gc_is_hard": true,
      "avoid_repeats_len": 12,
      "require_atg_start": true
    }
  },
  {
    "id": "cas9_human",
    "tier": "full",
    "fasta": "../tests/data/proteins/Q99ZW2.faa",
    "source": "UniProt Q99ZW2 v1, retrieved 2026-08-01",
    "organism": "h_sapiens",
    "table_id": 1,
    "budget_s": 120,
    "cai_floor": 0.78,
    "constraints": {
      "forbidden_motifs": ["GAATTC", "GGATCC", "AAGCTT", "TCTAGA", "AAAAAAAA"],
      "gc_min": 0.45, "gc_max": 0.65,
      "gc_window": 100, "gc_is_hard": true,
      "avoid_repeats_len": 15,
      "require_atg_start": true
    }
  }
]
```

### `benchmarks/baseline.json` shape

```json
{
  "generated_from_sha": "…",
  "seeds": [11, 22, 33, 44, 55],
  "cai_reference_set_version": "ecoli_heg_2019-r3",
  "cai_reference_set_sha256": "…",
  "viennarna_version": "2.7.2",
  "energy_params": "rna_turner2004",
  "cases": { "gfp_ecoli": { "cai": { "mean": 0.8631, "stdev": 0.0042 } } }
}
```

### `tests/differential/test_dnachisel_agrees.py` (nightly)

Do **not** compare two optimizers' sequences — they legitimately differ on almost every codon, so the test either never passes or gets weakened to meaninglessness. The valuable differential is on the *checking* side.

```python
import pytest
from dnachisel import (
    AvoidPattern,
    DnaOptimizationProblem,
    EnforceGCContent,
    EnforceTranslation,
    Location,
)

from bt5 import optimize
from bt5.constraints import Constraints
from tests.helpers import read_fasta
from benchmarks.run_panel import PANEL


@pytest.mark.slow
@pytest.mark.parametrize("case", PANEL, ids=lambda c: c["id"])
def test_dnachisel_confirms_our_constraints(case):
    """Cross-validate OUR CHECKER against an independent, peer-reviewed
    re-implementation of the same predicates."""
    protein = read_fasta(case["fasta"])
    cons = Constraints(
        **case["constraints"], organism=case["organism"], table_id=case["table_id"], seed=7
    )
    dna = optimize(protein, cons)

    constraints = [EnforceTranslation(location=Location(0, len(dna), 1))]
    constraints += [AvoidPattern(m) for m in cons.forbidden_motifs]
    if cons.gc_is_hard:
        constraints.append(
            EnforceGCContent(mini=cons.gc_min, maxi=cons.gc_max, window=cons.gc_window)
        )
    problem = DnaOptimizationProblem(sequence=dna, constraints=constraints, logger=None)
    evals = problem.constraints_evaluations()
    assert evals.all_evaluations_pass(), evals.to_text()
```

### Test corpus layout

```
tests/data/
  proteins/                        # each with a sidecar .json: accession, version, retrieval date
    P42212.faa           GFP, A. victoria                    238 aa
    X5DSL3.faa           mCherry                             236 aa
    Q99ZW2.faa           SpCas9, S. pyogenes                1368 aa  (size/runtime canary)
    trastuzumab_hc.faa   IgG1 heavy chain                   ~450 aa  (repeat-prone CDRs)
    P04637.faa           human p53                           393 aa
    gc_rich_human.faa    >65% GC native gene                          (band-conflict case)
  adversarial/
    poly_MW.faa          zero degrees of freedom
    poly_A_prone.faa     forces homopolymer runs
    junction_trap.faa    residue pairs whose every codon pair spells EcoRI
    infeasible.json      provably unsolvable -> must raise InfeasibleConstraints
    single_residue.faa   "M"
    titin_fragment.faa   10k aa  (nightly only -- memory/runtime ceiling)
  regressions/                     # one dir per fixed bug, named for the issue
    issue_0041_minus_strand_bamhi/{protein.faa,constraints.json,README.md}
    issue_0057_stop_codon_dropped_when_no_trailing_star/…
```

Every bugfix PR adds a `regressions/issue_NNNN_*` fixture that fails on the parent commit — enforced by the PR template and by CLAUDE.md §6. This is what stops the same bug shipping twice when next month's agent has no memory of this month's incident.

### `pyproject.toml` — the sections the gates depend on

```toml
[project]
requires-python = "==3.11.*"   # load-bearing: ViennaRNA has no cp312/cp313 wheels

[project.optional-dependencies]
dev = [
  "pytest>=8", "pytest-cov", "pytest-xdist",
  "hypothesis>=6.100", "syrupy>=6",
  "biopython>=1.83", "numpy>=1.26", "pyyaml",
  "ruff==0.16.4", "mypy>=1.11", "mutmut>=3",
]
bench = ["viennarna==2.7.2"]   # PINNED: energy parameter set moves every dG
diff  = ["dnachisel>=3.2"]     # nightly differential only

[tool.ruff]
target-version = "py311"
line-length = 100

[tool.ruff.lint]
select = ["E","F","W","I","N","UP","B","A","C4","SIM","ARG","PTH","RUF","PGH","TID","ANN","S","NPY","PL"]
# PGH003 = blanket type-ignore, PGH004 = blanket noqa

[tool.ruff.lint.per-file-ignores]
"tests/**" = ["S101", "ANN", "PLR2004"]

[tool.mypy]
python_version = "3.11"
strict = true                   # includes warn_unused_ignores
warn_unreachable = true
enable_error_code = ["ignore-without-code", "redundant-expr", "truthy-bool"]

[[tool.mypy.overrides]]
module = ["RNA.*", "dnachisel.*"]   # no type stubs shipped
ignore_missing_imports = true

[tool.pytest.ini_options]
addopts = "-ra --strict-markers --strict-config"
testpaths = ["packages/engine/tests", "packages/server/tests", "tests"]
markers = [
  "slow: long-running folding / property tests (nightly only)",
  "differential: requires the `diff` extra",
]
filterwarnings = ["error::Bio.BiopythonWarning"]

[tool.coverage.run]
branch = true
source = ["packages/engine/src", "packages/server/src"]

[tool.mutmut]
# Deliberately scoped to the deterministic correctness spine. Mutating the
# search heuristics produces a flood of equivalent mutants (a changed weight
# yields a different-but-valid sequence, killing nothing), which trains
# everyone to ignore the report.
paths_to_mutate = [
  "packages/engine/src/bt5/verify.py",
  "packages/engine/src/bt5/codon_table.py",
  "packages/engine/src/bt5/translate.py",
  "packages/engine/src/bt5/constraints/motifs.py",
  "packages/engine/src/bt5/metrics/cai.py",
]
tests_dir = ["tests/invariants", "tests/data_integrity"]
```

### `.pre-commit-config.yaml`

Local convenience only. **Do not add a required `pre-commit run --all-files` CI job** — a fresh clone has no hooks installed, agents work in ephemeral checkouts and essentially never run them, and re-running the same tools in CI adds ~40s for zero signal plus a second source of truth that drifts.

```yaml
default_install_hook_types: [pre-commit, commit-msg]
repos:
  - repo: https://github.com/pre-commit/pre-commit-hooks
    rev: v6.0.0
    hooks:
      - id: check-yaml
        args: ["--unsafe"]        # needed for GitHub Actions ${{ }} tags
      - id: check-toml
      - id: check-merge-conflict
      - id: end-of-file-fixer
      - id: trailing-whitespace
      - id: check-added-large-files
        args: ["--maxkb=512"]
      - id: mixed-line-ending
        args: ["--fix=lf"]        # keeps the Windows matrix row honest

  # NOTE: the hook id was renamed from `ruff` to `ruff-check`.
  # Linter before formatter: ruff's fixes can produce code needing reformatting.
  - repo: https://github.com/astral-sh/ruff-pre-commit
    rev: v0.16.4
    hooks:
      - id: ruff-check
        args: [--fix]
      - id: ruff-format

  - repo: https://github.com/astral-sh/uv-pre-commit
    rev: 0.12.6
    hooks:
      - id: uv-lock
```

---

## 7. Ordering

The constraint: **parallel agent sessions start the moment the foundation lands**, so anything not in place by then is retrofitted onto a moving target.

### PR #0 — "Foundation" (you write this by hand, merged directly or via a single PR before the ruleset goes live)

Everything an agent could later weaken must exist **before** any agent has write access to the repo.

**Ship:**

| Item | Why it must be first |
|---|---|
| Monorepo skeleton: `packages/engine`, `packages/server`, `apps/web`, `pyproject.toml`, `uv.lock`, `pnpm-workspace.yaml`, `pnpm-lock.yaml` | Nothing else can be built or tested. |
| **`packages/engine/src/bt5/verify.py`** with real invariants I1–I8 | The oracle must predate the optimizer. Every other gate is derived from it. |
| `optimize()` stub calling `verify_solution()`, `InfeasibleConstraints` defined | Establishes the trichotomy contract before any agent writes search code. |
| `tests/invariants/` + `tests/strategies.py` + `tests/conftest.py` (Hypothesis profiles) | The properties are the spec. |
| `tests/data/proteins/` (≥4 real proteins with provenance sidecars) + `adversarial/` | Agents cannot invent a good corpus. |
| `data/` codon tables + `MANIFEST.sha256` + `tests/data_integrity/` | Reference data must be correct before anything consumes it. |
| `benchmarks/{panel.json,tolerances.yaml,run_panel.py}` + an initial `baseline.json` | The panel needs a real baseline; generate it from the stub optimizer and regenerate once the real one lands. |
| `.github/workflows/ci.yml` (all jobs + `required-checks`) | |
| `.github/workflows/{codeql,dependency-review,main-is-green,nightly,regen-goldens}.yml` | |
| `.github/scripts/{path_guard.py,assert_gate_complete.py}` | |
| `.github/CODEOWNERS`, PR template, issue templates, `dependabot.yml`, `renovate.json` | |
| **`CLAUDE.md`** | The contract must exist before session one. |
| `.pre-commit-config.yaml`, `.gitattributes` (`* text=auto eol=lf`) | |
| Labels: `approved:oracle-change`, `approved:algorithm-change`, `approved:data-change`, `approved:ci-change`, `ci-full`, `broken-main`, `nightly-failure`, `scientific-correctness` | `path_guard.py` references them by name. |

**Then, in order:**

1. Run the settings script (§1) + toggle secret scanning, push protection, Dependabot in the UI.
2. Push the foundation to `main` **before** creating the ruleset.
3. Open **one throwaway PR** and confirm CI is green and `required-checks` reports. You cannot dry-run a ruleset on this plan, so this is your only rehearsal.
4. Create both rulesets (§2). Open a **second** throwaway PR and confirm the merge box says exactly what you expect — three green required checks, squash-only, no approval requirement, no bypass available to you.
5. Only now start agent sessions.

### Phase 1 — first week, in parallel with early agent work

- `e2e` job + Playwright config (needs a real server and UI first). Add it to `required-checks.needs` in the same PR.
- Golden test cases + first `__snapshots__` (needs a real optimizer). Add `goldens-not-hand-edited` once snapshots exist.
- Install CodeRabbit (comment-only, zero config on public repos).
- Install the Claude GitHub App via `/install-github-app`; store `CLAUDE_CODE_OAUTH_TOKEN` (subscription-billed — the single largest cost lever for a solo dev). Land `claude.yml` (`@claude` mentions) — no required check involved.

### Phase 2 — week two, once you've felt the merge loop

- `agent-approval-check.yml`. Land it **not required** for a few days, confirm `/approve <sha>` works and that the status re-arms on push, then add `agent-approval-check` to the ruleset.
- `claude-review.yml`. Land it **not required**. Watch the false-positive rate on `blocking=true` for a week.
- Install Renovate (GitHub App) once `uv.lock` and `pnpm-lock.yaml` are stable enough that lockfile-maintenance PRs aren't constant conflict.

### Phase 3 — week three or four

- Promote `claude-review-gate` to a required check once you trust the veto scope.
- Consider CodeQL results as a ruleset rule (severity threshold) rather than a status check.
- Pin all third-party actions by commit SHA with `# vN` comments (supply-chain hygiene on a public repo). Renovate's `pinDigests` maintains them.
- Enable mutation-score floor enforcement (nightly, non-blocking until the spine is stably ≥90%).

### Phase 4 — when the pain justifies it

- **Transfer to a free organization** and enable merge queue. Add the `merge_queue` rule, flip `strict` back on, and the queue handles up-to-dateness for you. 🔴 Every workflow above already carries `merge_group:` — omitting it inside a queue reproduces the never-reports deadlock in a place that is much harder to diagnose.

```json
{ "type": "merge_queue",
  "parameters": { "merge_method": "SQUASH", "grouping_strategy": "ALLGREEN",
                  "max_entries_to_build": 5, "min_entries_to_merge": 1,
                  "min_entries_to_merge_wait_minutes": 5, "max_entries_to_merge": 5,
                  "check_response_timeout_minutes": 60 } }
```

### Never do

Add a bypass actor. Enable `require_last_push_approval`. Enable `require_code_owner_review`. Enable `required_signatures` (Claude Code sessions sign nothing by default; the push is rejected mid-task with a message agents handle badly). Enable CodeQL default setup. Enable "Allow GitHub Actions to create and approve pull requests". Add a `paths:` filter to `on:` in any workflow owning a required check.

---

## 8. Failure modes at 8 concurrent AI PRs

### 8.1 🔴 CI concurrency starvation — the binding constraint

**What happens.** Free personal accounts get **20 concurrent jobs** (5 max on macOS). A Python PR runs ~6 CI jobs + ~4 CodeQL jobs + dependency-review + claude-review ≈ 12. Eight concurrent PRs demand ~96 job slots against 20. You are 4–5× oversubscribed. Queue depth means a 10-minute pipeline becomes a 40-minute wait, agents idle, and you can no longer tell "slow" from "hung."

**Mitigations, in order of impact:**

1. **`cancel-in-progress: true` on PR events** (already in every workflow). Agents push corrective commits in bursts; without this each push stacks a full pipeline.
2. **Aggressive in-job path filtering** (already in place). A web-only PR runs 3 jobs, not 12.
3. **Move CodeQL off PRs** if the queue is choking: change `codeql.yml` to `push: [main] + schedule` only and drop `codeql-passed` from required checks. Costs you pre-merge SAST; buys back ~4 jobs × 8 PRs = 32 slots. Do this first if you have to cut something.
4. **`benchmark-gate` on `engine` paths only** (already gated) — it is the longest job.
5. **OS matrix off PRs entirely** (already: PRs are ubuntu-only; the full matrix is nightly or `ci-full`-labeled). macOS is separately capped at 5 concurrent jobs even on Pro.
6. **Agents open PRs as drafts** and mark ready when they believe they're done. `claude-review` skips drafts; you can add `if: github.event.pull_request.draft == false` to `e2e` and `benchmark-gate` too.
7. **Cap agent parallelism at 4–5 sessions.** This is the honest answer. Eight simultaneous agents on a free personal account is above the hardware.
8. **Move to an org** — doesn't raise the job cap, but merge queue eliminates the re-run cascades that make the cap bite.

### 8.2 🔴 Lockfile conflict storm

**What happens.** `uv.lock` and `pnpm-lock.yaml` are the highest-churn, least-mergeable files in the repo. With 8 agents, several will touch dependencies; the first merge conflicts all others; a three-way merge of a lockfile is garbage; an agent that "resolves" it by hand produces a lockfile that satisfies neither `--locked` nor reality.

**Mitigations.** CLAUDE.md §5 forbids hand-resolving — rebase and regenerate. Dependency changes go in dedicated small PRs, merged first, before feature work. `uv sync --locked` fails loudly on a stale lock rather than producing a silently-wrong environment. Renovate is capped at `prConcurrentLimit: 3` / `prHourlyLimit: 2` so bot PRs don't compete with agents for the same 20 slots.

### 8.3 🔴 Semantic conflict on `main` (the price of loose checks)

**What happens.** The engine agent renames a Pydantic field; the web agent writes a client against the old name. Both PRs are green against their own (older) `main`. Git merges cleanly. Runtime breaks. With 8 PRs merging in a day this is not hypothetical — it is scheduled.

**Mitigations.** `main-is-green.yml` runs the *entire* suite with no path filters on every push to `main` and opens (or comments on) a `broken-main` issue within minutes. "Always suggest updating pull request branches" gives you a one-click update for any PR you judge semantically stale. CLAUDE.md §1 forbids cross-module edits and requires an issue for interface changes; the Claude review prompt makes "interface drift" priority 3 explicitly. The real fix is merge queue — an org transfer.

**Do not** respond by turning `strict` on. With 8 open PRs, every merge invalidates the other 7, each must be updated (new SHA) and re-run full CI; at ~8 minutes a suite you burn ~an hour of CI per merge and the slowest agent's PR may never reach the front. That is the livelock, and it is worse than the semantic conflicts.

### 8.4 Aggregator gate rot

**What happens.** An agent adds a job to `ci.yml` and doesn't add it to `required-checks.needs`. The job runs, it can fail, and the gate never looks at it — the PR merges green. Agents editing `ci.yml` do this constantly. There is no visible signal.

**Mitigations.** `.github/**` is label-gated by `path_guard.py` (`approved:ci-change`), so an agent cannot silently edit CI at all. CODEOWNERS routes it to you. `assert_gate_complete.py` runs nightly and parses `ci.yml`, failing if any job is unreachable from `needs`, if the gate doesn't use `always()`, or if it has no `needs.*.result` inspection step. Treat any `.github/**` diff as owner-review-required, full stop.

### 8.5 Hypothesis finds a pre-existing bug on an unrelated agent's PR

**What happens.** With 8 PRs × 200 examples each, Hypothesis explores a lot of space. A counterexample surfaces on the web agent's PR for a latent engine bug the web agent didn't cause. Blocking that PR indefinitely stalls the parallel model; the agent's natural next move is to weaken the assertion.

**Mitigation — a standing rule, written into CLAUDE.md §8.** File the counterexample as a fixture under `tests/data/regressions/` plus an issue. If the failure reproduces on the merge base, let that PR proceed (you merge it manually after confirming). Never weaken the property. The `.hypothesis` cache keys the counterexample so it replays first from then on.

### 8.6 Benchmark baseline thrash

**What happens.** Two agents both legitimately improve CAI. The first merges and the baseline is regenerated; the second's PR now compares against a moved baseline and either reads as a regression or gets its own regeneration — and now nobody knows which numbers are real.

**Mitigations.** The baseline is only ever regenerated by the `regen-goldens` workflow, which opens its own PR labeled `approved:algorithm-change`. Agents cannot touch `benchmarks/baseline.json` (label-gated). The comparability guard (`cai_reference_set_version`, `viennarna_version`) hard-fails rather than silently comparing incomparable numbers. Tolerances are direction-aware, so an *improvement* is never a failure — two independent improvements both pass without a baseline change.

### 8.7 Claude review cost and turn truncation

**What happens.** 8 PRs × several pushes each × a full review = real money, and `--max-turns` truncates rather than fails, so a large PR can produce a *partial* review that still reports `blocking: false`.

**Mitigations.** `--model claude-sonnet-5` ($2/$10 per MTok vs Opus 5 at $5/$25). `concurrency` + `cancel-in-progress` so a burst pays for one review. Drafts skipped. `CLAUDE_CODE_OAUTH_TOKEN` bills against a Pro/Max subscription instead of metered API. The prompt explicitly instructs: *"If you could not complete the review, set blocking=true and say so"*, and the gate step fails closed on empty `structured_output`. Watch actual spend in the Anthropic Console for the first week — a Sonnet review of a moderate diff lands in the low tens of cents to ~$1, scaling with how much of the repo Claude reads.

### 8.8 Auto-merge silently disarms

**What happens.** You enable auto-merge on 5 agent PRs and walk away. Auto-merge is disabled automatically if the base branch changes or if someone without write permission pushes to the head branch. You come back to PRs sitting green and unmerged.

**Mitigations.** Auto-merge is a convenience, not a guarantee — check the PR list, don't assume. Note also the "Enable auto-merge" button only appears while the PR is *actually blocked*; if CI already finished green you just merge normally. And auto-merge merges without re-validating against the new `main` (that's merge queue's job) — which is exactly the residual risk `main-is-green` covers.

### 8.9 The path-guard label is self-applicable

**What happens.** Agents authenticate as you and therefore have write access, so in principle an agent can apply its own `approved:oracle-change` label.

**Mitigations.** Two options. (a) Accept it: the point of the label is to create a *visible* review checkpoint in the PR timeline, not an unbypassable wall — and you are the only one who can merge anyway. A self-applied label is anomalous and obvious in the timeline. (b) Have agents work from **forks**: fork PRs get a read-only `GITHUB_TOKEN`, cannot modify workflow files in a way that runs, and cannot apply labels. Option (b) is materially stronger and worth adopting if you scale past ~5 agents.

### 8.10 Everything blocks because CI itself broke

**What happens.** A flaky runner, an expired token, an upstream action yanked — and with an empty bypass list you cannot merge anything.

**Mitigation.** That is the intended cost, and it's the right one. The escape hatch is the `emergency-bypass` pattern from §2: set `main-protection` to `"enforcement": "disabled"`, merge, flip it straight back. It leaves an entry in `/rulesets/{id}/history`. **Never** solve it by adding yourself as a bypass actor — bypass is per-actor, not per-rule, and one permanent bypass silently disables every rule for the only person who merges.

---

### Quick reference: the eight things that silently disable this entire system

1. `paths:` in `on:` on any workflow that owns a required check → permanent deadlock, no error, no timeout.
2. `if: success()` (or no `if:`) on the gate job → gate gets *skipped* when a dependency fails, and skipped satisfies a required check.
3. `if: always()` on the gate with no `needs.*.result` inspection step → gate reports success while everything under it burns.
4. A new job not added to `required-checks.needs` → runs, can fail, unenforced.
5. Yourself in `bypass_actors` → every rule off for the only person who merges.
6. `required_approving_review_count: 1` or `require_last_push_approval: true` → mathematically unmergeable, forcing (5).
7. Required check source left as "any source" instead of pinned to the Actions app → any token with `statuses:write` can forge green.
8. "Allow GitHub Actions to create and approve pull requests" turned on → an agent that can edit `.github/workflows` can approve its own PR.