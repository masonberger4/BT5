"""The oracle must not share a code path with anything it validates.

Three callers of one pure function is zero independence. bt5/verify.py may import
only from bt5.core (types and result) and the standard library / Biopython.
"""

from __future__ import annotations

import ast
from pathlib import Path

VERIFY = Path(__file__).resolve().parents[2] / "packages" / "engine" / "src" / "bt5" / "verify.py"
ALLOWED_PREFIXES = (
    "bt5.core.types",
    "bt5.core.result",
    "Bio",
    "__future__",
    "collections",
    "typing",
)


def test_oracle_imports_no_lane_module() -> None:
    tree = ast.parse(VERIFY.read_text(encoding="utf-8"))
    bad: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            mod = node.module
            if mod.startswith("bt5.") and not mod.startswith(ALLOWED_PREFIXES):
                bad.append(mod)
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.startswith("bt5.") and not alias.name.startswith(ALLOWED_PREFIXES):
                    bad.append(alias.name)
    assert not bad, (
        f"bt5/verify.py must stay independent of the lanes it validates, but imports: {bad}"
    )


def test_oracle_does_not_import_the_rules_registry() -> None:
    src = VERIFY.read_text(encoding="utf-8")
    assert "from bt5.core.registry" not in src
    assert "bt5.rules" not in src
