#!/usr/bin/env python3
"""Status line: which tier this turn actually ran on, plus context pressure.

The point of the routing config is that different work runs at different cost. Without
this you cannot see which tier a turn used, so a misrouted session looks exactly like a
correctly routed one. Also surfaces whether the venv exists, because every gate in this
repo depends on it and its absence is otherwise invisible until a command fails oddly.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys


def branch() -> str:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            capture_output=True,
            text=True,
            encoding="utf-8",
            timeout=3,
        )
    except (OSError, subprocess.SubprocessError):
        return "-"
    return out.stdout.strip() or "-"


def main() -> int:
    try:
        data = json.load(sys.stdin)
    except Exception:
        return 0
    if not isinstance(data, dict):
        return 0

    model = (data.get("model") or {}).get("display_name") or "?"
    effort = (data.get("effort") or {}).get("level") or ""
    used = (data.get("context_window") or {}).get("used_percentage")

    parts = [f"{model} {effort}".strip(), branch()]

    if isinstance(used, (int, float)):
        filled = max(0, min(10, int(used // 10)))
        bar = "#" * filled + "." * (10 - filled)
        parts.append(f"ctx {bar} {used:.0f}%")

    parts.append("venv" if os.access(".venv/bin/python", os.X_OK) else "NO VENV")
    print("  ".join(parts))
    return 0


if __name__ == "__main__":
    sys.exit(main())
