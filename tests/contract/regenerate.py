"""Re-record the frozen contract: `python tests/contract/regenerate.py`.

Run this in the PR that changes `core/`. It rewrites the manifest and every
fixture from the live code, and prints what changed and at what severity, so
you find out you owe an RFC before CI tells you.

It deliberately does NOT write the amendment entry for you. A MAJOR change is
supposed to cost a paragraph of thought; generating that paragraph would defeat
the point of requiring it.
"""

from __future__ import annotations

import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import fixtures  # noqa: E402
import surface  # noqa: E402


def main() -> int:
    live = surface.extract()

    if surface.MANIFEST_PATH.exists():
        manifest = json.loads(surface.MANIFEST_PATH.read_text(encoding="utf-8"))
        changes = surface.diff(manifest.get("surface", {}), live)
    else:
        manifest = {
            "frozen_at": date.today().isoformat(),
            "contract_version": 1,
            "amendments": [],
        }
        changes = ()

    manifest["surface"] = live
    surface.MANIFEST_PATH.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    fixtures.record()

    entries = sum(len(v) for v in live.values())
    print(f"recorded {entries} contract entries across {len(live)} modules")
    print(f"recorded {len(fixtures.specimens())} fixtures")

    if not changes:
        print("no change to the contract surface")
        return 0

    major = surface.majors(changes)
    print(f"\n{len(changes)} change(s), {len(major)} MAJOR:")
    for change in changes:
        print(f"  {change}")

    if major:
        print(
            "\nMAJOR changes need an amendment: bump contract_version, add an "
            "entry to `amendments` naming an RFC under docs/rfcs/, and ship the "
            "deprecation shim. CI checks this against main."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
