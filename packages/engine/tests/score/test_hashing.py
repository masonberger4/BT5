"""The content hash on the tube label.

Two runs producing two different sequences under one name is how a lab ends up
with two tubes and an irreproducible result.
"""

from __future__ import annotations

from bt5.score import HASH_LENGTH, design_hash


def test_it_is_deterministic_across_calls() -> None:
    assert design_hash("ATGAAACCC") == design_hash("ATGAAACCC")


def test_different_designs_hash_differently() -> None:
    assert design_hash("ATGAAACCC") != design_hash("ATGAAACCG")


def test_it_is_stable_across_processes() -> None:
    """Not Python's hash(), which is salted per process -- a value that changes
    between runs cannot identify anything. CI sets PYTHONHASHSEED=0 for exactly
    this class of bug, and relying on that would be relying on a CI setting."""
    import subprocess
    import sys

    code = (
        "import sys; sys.path.insert(0, 'packages/engine/src');"
        "from bt5.score import design_hash; print(design_hash('ATGAAACCC'))"
    )
    out = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        encoding="utf-8",
        check=True,
        env={"PATH": "/usr/bin:/bin"},
    )
    assert out.stdout.strip() == design_hash("ATGAAACCC")


def test_the_context_changes_the_hash() -> None:
    a = design_hash("ATGAAACCC", context={"table_id": 11})
    b = design_hash("ATGAAACCC", context={"table_id": 1})
    assert a != b


def test_context_key_order_does_not_change_the_hash() -> None:
    """Otherwise the hash identifies the caller's iteration order as much as
    the design."""
    a = design_hash("ATG", context={"seed": 1, "table_id": 11})
    b = design_hash("ATG", context={"table_id": 11, "seed": 1})
    assert a == b


def test_case_does_not_change_the_hash() -> None:
    assert design_hash("atgaaaccc") == design_hash("ATGAAACCC")


def test_it_is_short_enough_to_write_on_a_tube() -> None:
    h = design_hash("ATGAAACCC")
    assert len(h) == HASH_LENGTH
    assert all(c in "0123456789abcdef" for c in h)
