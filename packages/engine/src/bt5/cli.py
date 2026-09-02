"""`bt5` -- the command-line entry point over `bt5.design.design()`.

Before this, running BT5 meant writing Python: importing `design()`, building a
`VectorBackbone` and a `Modality`/`HostId` pair by hand, and writing the returned
GenBank text out yourself. This module is that same call, wired to argv, so a
user who never opens a Python REPL can still get an annotated construct out of a
protein and a backbone.

`table_id` is required with no default here for the same reason `design()` never
defaults it (CLAUDE.md section 3.1): a wrong genetic code table is a silently
wrong protein no assay catches for months, and a CLI flag with a convenient
default would put that mistake one keystroke away.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from bt5.core.context import HostId, Modality
from bt5.core.result import InfeasibleConstraints, VerificationError
from bt5.design import DesignError, design
from bt5.score.order import OrderEntry, OrderError, entry_name, write_csv
from bt5.vector import VectorError, read_genbank

#: Errors that name a remedy in their own message -- reported cleanly on stderr
#: rather than as a Python traceback, which is not something a CLI user can act on.
_REPORTABLE_ERRORS = (
    DesignError,
    InfeasibleConstraints,
    VerificationError,
    VectorError,
    OrderError,
    OSError,
    # An unrecognised --table-id (FileTableProvider.genetic_code) and a
    # malformed --backbone file (Bio.SeqIO's own parse errors) both surface
    # as ValueError; a raw traceback for either is not something a CLI user
    # can act on the way a clean "bt5: ..." line is.
    ValueError,
)


def _read_protein(args: argparse.Namespace) -> str:
    if args.protein is not None:
        return str(args.protein).strip()
    text = Path(args.protein_file).read_text()
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bt5",
        description="Back-translate a protein and codon-optimize it in the context "
        "of an assembled construct.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    design_parser = subparsers.add_parser(
        "design", help="Design one protein into one backbone, verified end to end."
    )
    design_parser.add_argument(
        "--backbone", required=True, type=Path, help="GenBank vector map to insert into."
    )
    protein_source = design_parser.add_mutually_exclusive_group(required=True)
    protein_source.add_argument(
        "--protein", help="Protein sequence (one-letter code), starting with the initiator M."
    )
    protein_source.add_argument(
        "--protein-file",
        type=Path,
        help="File containing the protein sequence (plain text or single-record FASTA).",
    )
    design_parser.add_argument(
        "--table-id",
        required=True,
        type=int,
        help="NCBI genetic code table id. Required, never defaulted: the wrong table "
        "silently mistranslates the protein.",
    )
    design_parser.add_argument(
        "--modality",
        required=True,
        choices=[m.value for m in Modality],
        help="How the construct is delivered.",
    )
    design_parser.add_argument(
        "--host",
        dest="hosts",
        action="append",
        required=True,
        choices=[h.value for h in HostId],
        metavar="HOST",
        help="Expression host; repeat --host for more than one. Choices: "
        + ", ".join(h.value for h in HostId),
    )
    design_parser.add_argument("--seed", type=int, default=0, help="RNG seed (default 0).")
    design_parser.add_argument(
        "--site-label", default=None, help="Label of the insertion site to use, if ambiguous."
    )
    design_parser.add_argument("--preset-id", default=None, help="Vendor preset id for the report.")
    design_parser.add_argument(
        "--max-candidates", type=int, default=256, help="Solver candidate cap (default 256)."
    )
    design_parser.add_argument(
        "--out-genbank",
        type=Path,
        default=Path("design.gb"),
        help="Where to write the annotated GenBank (default ./design.gb).",
    )
    design_parser.add_argument(
        "--out-order",
        type=Path,
        default=None,
        help="Also write a Name,Sequence order CSV for the designed CDS.",
    )
    design_parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress the rendered QC report on stdout; still prints file paths written.",
    )
    design_parser.set_defaults(func=_cmd_design)

    return parser


def _cmd_design(args: argparse.Namespace) -> int:
    protein = _read_protein(args)
    backbone = read_genbank(args.backbone)
    hosts = tuple(HostId(h) for h in args.hosts)

    result = design(
        backbone=backbone,
        protein=protein,
        table_id=args.table_id,
        modality=Modality(args.modality),
        hosts=hosts,
        site_label=args.site_label,
        seed=args.seed,
        preset_id=args.preset_id,
        max_candidates=args.max_candidates,
    )

    args.out_genbank.write_text(result.genbank)
    print(f"wrote {args.out_genbank}")

    if args.out_order is not None:
        candidate = result.result.candidates[0]
        entry = OrderEntry(
            name=entry_name(backbone.name, candidate.design_hash),
            sequence=candidate.cds,
        )
        write_csv([entry], args.out_order)
        print(f"wrote {args.out_order}")

    if not args.quiet:
        print(result.rendered)
        for note in result.notes:
            print(f"note: {note}")

    return 0


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except _REPORTABLE_ERRORS as exc:
        print(f"bt5: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
