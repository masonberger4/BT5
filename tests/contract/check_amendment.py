"""Classify this branch's contract changes against a baseline. CI's half of the gate.

    python tests/contract/check_amendment.py <baseline-manifest.json>

The local test cannot do this. Once you regenerate the manifest, the recorded
contract and the live code agree by construction, and the severity of what you
just did is only visible against the version you started from -- which in CI is
main's manifest, and nowhere else.

The decision itself lives in `surface.review()` so it can be tested against the
cases it exists to reject. This file is the I/O around it: read the baseline,
print what was found, choose an exit code.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import surface  # noqa: E402

PROTOCOL = (
    "MAJOR means something that exists today stops working. docs/PLAN.md's\n"
    "amendment protocol: an RFC, a deprecation shim with a\n"
    '`model_validator(mode="before")`, the two-window rule, and\n'
    "test_backward_compat re-parsing every recorded fixture IN THIS PR.\n"
    "Add docs/rfcs/<name>.md, bump contract_version, and record the amendment."
)


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print("usage: check_amendment.py <baseline-manifest.json>")
        return 2

    live = surface.extract()
    current = json.loads(surface.MANIFEST_PATH.read_text())

    # The manifest on this branch must already match the live code. The local
    # test says so too, but repeating it here means CI never classifies against
    # a manifest that was never regenerated.
    drift = surface.diff(current["surface"], live)
    if drift:
        print("::error::the recorded contract does not match bt5.core on this branch")
        for change in drift:
            print(f"  {change}")
        print("Run `python tests/contract/regenerate.py`.")
        return 1

    baseline_path = Path(argv[1])
    if not baseline_path.exists() or not baseline_path.read_text().strip():
        print("No contract manifest on the base branch: this commit IS the freeze.")
        print(f"Recording {sum(len(v) for v in live.values())} entries as contract version 1.")
        return 0

    verdict = surface.review(json.loads(baseline_path.read_text()), current, live)

    if not verdict.changes:
        print("Contract unchanged.")
        return 0

    print(f"{len(verdict.changes)} contract change(s) against the base branch:")
    for change in verdict.changes:
        print(f"  {change}")

    if not verdict.major:
        print("\nAll MINOR: nothing that exists stops working. Fast path, no RFC needed.")
        return 0

    if verdict.ok:
        amendment = verdict.amendment or {}
        print(
            f"\n{len(verdict.major)} MAJOR change(s), amended by contract version "
            f"{current['contract_version']} ({amendment.get('rfc')})"
        )
        print(f"  {amendment.get('summary')}")
        return 0

    print(f"\n::error::{len(verdict.major)} MAJOR contract change(s) without a complete amendment")
    for problem in verdict.problems:
        print(f"  - {problem}")
    print(f"\n{PROTOCOL}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
