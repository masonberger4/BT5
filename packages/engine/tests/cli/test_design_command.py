"""`bt5 design` end to end: the CLI is a thin argv-to-design() adapter, and these
tests exist to pin that adapter, not to re-test `design()` itself (that is
`tests/design/test_skeleton.py`'s job).
"""

from __future__ import annotations

import csv
import io
import subprocess
import sys
from pathlib import Path

import pytest
from bt5.cli import main
from bt5.vector import read_genbank

ROOT = Path(__file__).resolve().parents[4]
MCS_PATH = ROOT / "tests" / "data" / "backbones" / "synthetic_mcs_ef1a.gb"

#: Same 140-residue protein `tests/design/test_skeleton.py` uses, so a failure
#: here and a failure there are directly comparable.
PROTEIN = (
    "MKLVTAAFERSKSVQNYVVSTKDSPLYYLRKWVRSGYKFDCEEVGLREHQGPAATYTPTQAIWRLTLPSPLL"
    "NVDVWQNSCKSLQHTASWKKHRFGLFTLVISPLIRLGEVASLCGLCEHTATSEVKVCPIDCLQSPTSF"
)

_BASE_ARGS = [
    "design",
    "--backbone",
    str(MCS_PATH),
    "--protein",
    PROTEIN,
    "--table-id",
    "1",
    "--modality",
    "lentiviral",
    "--host",
    "hek293",
    "--seed",
    "0",
]


def test_design_writes_genbank(tmp_path: Path) -> None:
    out_gb = tmp_path / "out.gb"
    exit_code = main([*_BASE_ARGS, "--out-genbank", str(out_gb), "--quiet"])

    assert exit_code == 0
    assert out_gb.exists()
    # Round-trips through the same reader real users will hand it back to BT5.
    parsed = read_genbank(out_gb)
    assert parsed.length > 0


def test_design_writes_order_csv(tmp_path: Path) -> None:
    out_gb = tmp_path / "out.gb"
    out_csv = tmp_path / "order.csv"
    exit_code = main(
        [*_BASE_ARGS, "--out-genbank", str(out_gb), "--out-order", str(out_csv), "--quiet"]
    )

    assert exit_code == 0
    assert out_csv.exists()
    rows = list(csv.reader(io.StringIO(out_csv.read_text())))
    assert rows[0] == ["Name", "Sequence"]
    assert len(rows) == 2
    name, sequence = rows[1]
    assert name
    assert set(sequence) <= set("ACGT")


def test_design_default_genbank_path_is_relative(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    exit_code = main([*_BASE_ARGS, "--quiet"])

    assert exit_code == 0
    assert (tmp_path / "design.gb").exists()


def test_design_rejects_protein_without_initiator(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    bad_args = list(_BASE_ARGS)
    bad_args[bad_args.index("--protein") + 1] = "AKLVTAAFERSKS"
    exit_code = main([*bad_args, "--out-genbank", str(tmp_path / "out.gb"), "--quiet"])

    assert exit_code == 1
    captured = capsys.readouterr()
    assert "bt5:" in captured.err
    assert "initiator" in captured.err


def test_design_rejects_unknown_host() -> None:
    bad_args = list(_BASE_ARGS)
    bad_args[bad_args.index("--host") + 1] = "not_a_real_host"

    with pytest.raises(SystemExit) as excinfo:
        main(bad_args)
    assert excinfo.value.code == 2


def test_design_requires_table_id() -> None:
    args = [a for a in _BASE_ARGS if a not in {"--table-id", "1"}]
    with pytest.raises(SystemExit) as excinfo:
        main(args)
    assert excinfo.value.code == 2


def test_design_accepts_multiple_hosts(tmp_path: Path) -> None:
    args = [*_BASE_ARGS, "--host", "human", "--out-genbank", str(tmp_path / "out.gb"), "--quiet"]
    exit_code = main(args)
    assert exit_code == 0


def test_design_reads_protein_from_fasta_file(tmp_path: Path) -> None:
    fasta = tmp_path / "protein.fasta"
    fasta.write_text(f">example\n{PROTEIN}\n")
    out_gb = tmp_path / "out.gb"
    args = [
        "design",
        "--backbone",
        str(MCS_PATH),
        "--protein-file",
        str(fasta),
        "--table-id",
        "1",
        "--modality",
        "lentiviral",
        "--host",
        "hek293",
        "--out-genbank",
        str(out_gb),
        "--quiet",
    ]
    exit_code = main(args)
    assert exit_code == 0
    assert out_gb.exists()


def test_python_dash_m_bt5_runs(tmp_path: Path) -> None:
    """`python -m bt5` is wired to the same CLI, not just the console script."""
    out_gb = tmp_path / "out.gb"
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "bt5",
            *_BASE_ARGS,
            "--out-genbank",
            str(out_gb),
            "--quiet",
        ],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert completed.returncode == 0, completed.stderr
    assert out_gb.exists()
