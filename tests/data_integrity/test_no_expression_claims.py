"""BT5 must never claim to predict expression.

All computable design features together explain only 5-31% (mean ~14%) of
protein-level variance, and nine benchmarked optimizers were a coin flip against
native. Any field name implying a predicted titer or yield is dishonest, so the
vocabulary is banned at the source level and (once the server lands) against the
GENERATED OpenAPI schema, which catches what grepping source misses.
"""

from __future__ import annotations

import re
from pathlib import Path

ENGINE = Path(__file__).resolve().parents[2] / "packages" / "engine" / "src"
BANNED = re.compile(
    r"\b(predicted_expression|expression_score|titer_prediction|predicted_titer|"
    r"fold_improvement|expression_level|predicted_yield)\b"
)


def test_no_banned_prediction_vocabulary_in_engine_source() -> None:
    offenders: list[str] = []
    for path in ENGINE.rglob("*.py"):
        for n, line in enumerate(path.read_text().splitlines(), 1):
            if BANNED.search(line) and "BANNED" not in line:
                offenders.append(f"{path.relative_to(ENGINE)}:{n}: {line.strip()}")
    assert not offenders, "banned prediction vocabulary:\n" + "\n".join(offenders)


def test_kmer_index_takes_no_external_database() -> None:
    """BIOSECURITY. Pointing a homology minimiser at an arbitrary target database
    turns BT5 into a general-purpose screening-evasion tool. KmerIndex.of() must
    accept a Construct and nothing else."""
    src = (ENGINE / "bt5" / "core" / "services.py").read_text()
    match = re.search(r"def of\(cls, ([^)]*)\)", src)
    assert match, "KmerIndex.of not found"
    params = match.group(1)
    assert "Construct" in params, "KmerIndex.of must take a Construct"
    for forbidden in ("database", "db", "reference", "fasta", "blast"):
        assert forbidden not in params.lower(), (
            f"KmerIndex.of must not accept {forbidden!r}: constraining it to the "
            f"assembled construct is what keeps BT5 from being a screening-evasion tool"
        )
