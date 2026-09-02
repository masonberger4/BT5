"""Build a highly-expressed codon-usage reference set for a host from NCBI RefSeq.

BT5 computes CAI (rule C1) as relative adaptiveness ``w`` against a *highly-expressed*
reference set, never against a whole-genome table (docs/research/brief.md:77, :282).
Only the E. coli set (Sharp & Li 1987) shipped, so C1 reported ``unavailable`` for every
eukaryotic host (issue #78).

This script reproduces the reference set from a fixed, auditable panel of canonical
highly-expressed genes -- cytosolic ribosomal proteins, translation elongation factors,
core glycolytic enzymes, and chaperonins -- an extension of the fallback chain named in
brief.md:282. Each gene symbol is resolved to its RefSeq transcript *gene-first* (via the
NCBI Gene database, confirming the official symbol, then linking to the gene's RefSeq RNA),
and every fetched transcript is re-verified to actually belong to its symbol before its
codons are counted. A symbol that cannot be resolved to a transcript whose ``/gene`` equals
that symbol contributes nothing and is disclosed -- never silently mapped to a paralog.

``w`` is built exactly as ``bt5.codon.tables.CodonUsage.from_counts`` does (pseudocount 0.5
before the family max). Provenance is the deliverable: the emitted JSON records every
contributing accession and its verified gene, so an auditor can confirm the cited source
supports each number.

Run: ``python build_reference_set.py <host_key> <taxid> <organism_label> <out.json>``
NCBI E-utilities, no API key: capped at 3 requests/second here.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from datetime import date
from io import StringIO

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord

# --- The highly-expressed panel -------------------------------------------------
# Canonical, constitutively high-abundance genes shared across vertebrate hosts.
# These extend the brief's fallback set (brief.md:282): cytosolic ribosomal proteins +
# elongation factors + GAPDH + chaperonins, widened with core glycolytic enzymes and
# cytoskeletal genes that sit in the top abundance tier of every proteome measured.
PANEL: dict[str, list[str]] = {
    "ribosomal_large_subunit": [
        "RPL3",
        "RPL4",
        "RPL5",
        "RPL6",
        "RPL7",
        "RPL7A",
        "RPL8",
        "RPL9",
        "RPL10",
        "RPL10A",
        "RPL11",
        "RPL12",
        "RPL13",
        "RPL13A",
        "RPL14",
        "RPL15",
        "RPL18",
        "RPL18A",
        "RPL19",
        "RPL21",
        "RPL23",
        "RPL23A",
        "RPL24",
        "RPL26",
        "RPL27",
        "RPL27A",
        "RPL28",
        "RPL30",
        "RPL31",
        "RPL32",
        "RPL34",
        "RPL35",
        "RPL35A",
        "RPL36",
        "RPL37",
        "RPL37A",
        "RPL38",
        "RPLP0",
        "RPLP1",
        "RPLP2",
    ],
    "ribosomal_small_subunit": [
        "RPS2",
        "RPS3",
        "RPS3A",
        "RPS4X",
        "RPS5",
        "RPS6",
        "RPS7",
        "RPS8",
        "RPS9",
        "RPS10",
        "RPS11",
        "RPS12",
        "RPS13",
        "RPS14",
        "RPS15",
        "RPS15A",
        "RPS16",
        "RPS17",
        "RPS18",
        "RPS19",
        "RPS20",
        "RPS23",
        "RPS24",
        "RPS25",
        "RPS27",
        "RPS27A",
        "RPS28",
        "RPS29",
        "RPSA",
    ],
    "elongation_factors": ["EEF1A1", "EEF2", "EEF1B2", "EEF1G", "EEF1D"],
    "glycolytic": [
        "GAPDH",
        "ACTB",
        "PGK1",
        "ENO1",
        "PKM",
        "LDHA",
        "ALDOA",
        "TPI1",
        "PGAM1",
    ],
    "chaperonins": [
        "HSPA8",
        "HSP90AA1",
        "HSP90AB1",
        "CCT2",
        "CCT3",
        "CCT4",
        "CCT5",
        "TCP1",
        "HSPD1",
    ],
    "other_abundant": ["PABPC1", "YWHAZ", "TUBB", "TUBA1B", "VIM", "FTL", "FTH1"],
}

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_last_call = [0.0]


def _throttled_get(url: str) -> str:
    # NCBI allows 3 requests/second without an API key.
    wait = 0.34 - (time.monotonic() - _last_call[0])
    if wait > 0:
        time.sleep(wait)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=45) as resp:  # noqa: S310
                _last_call[0] = time.monotonic()
                return resp.read().decode("utf-8")
        except Exception:  # noqa: BLE001
            if attempt == 3:
                raise
            time.sleep(2**attempt)
    raise RuntimeError("unreachable")


def _esearch(db: str, term: str) -> list[str]:
    url = (
        f"{EUTILS}/esearch.fcgi?db={db}&retmode=json&retmax=25&tool=bt5"
        f"&term={urllib.parse.quote(term)}"
    )
    return json.loads(_throttled_get(url)).get("esearchresult", {}).get("idlist", [])


def _esummary(db: str, uids: list[str]) -> dict[str, object]:
    url = f"{EUTILS}/esummary.fcgi?db={db}&retmode=json&tool=bt5&id={','.join(uids)}"
    return json.loads(_throttled_get(url)).get("result", {})


def _refseq_rna_uids(gene_uid: str) -> list[str]:
    url = (
        f"{EUTILS}/elink.fcgi?dbfrom=gene&db=nuccore&retmode=json&tool=bt5"
        f"&linkname=gene_nuccore_refseqrna&id={gene_uid}"
    )
    linksets = json.loads(_throttled_get(url)).get("linksets", [])
    for ls in linksets:
        for db in ls.get("linksetdbs", []):
            if db.get("linkname") == "gene_nuccore_refseqrna":
                return list(db.get("links", []))
    return []


def resolve_accession(symbol: str, taxid: str) -> str | None:
    """Resolve a gene symbol to a RefSeq mRNA accession, gene-first.

    Finds the NCBI Gene record whose OFFICIAL symbol equals ``symbol`` (never an alias),
    links to that gene's RefSeq RNA transcripts, and returns the canonical accession
    (curated NM_ preferred over predicted XM_; lowest accession among isoforms). Returns
    None if the organism has no gene with that official symbol or no RefSeq RNA for it.
    The transcript is re-verified against its ``/gene`` after fetching (see main), so a
    wrong-gene resolution can never enter the counts silently.
    """
    org = f"txid{taxid}[Organism:exp]"
    gene_ids = _esearch("gene", f"{symbol}[sym] AND {org} AND alive[prop]")
    if not gene_ids:
        return None
    summ = _esummary("gene", gene_ids)
    gene_uid = next(
        (
            uid
            for uid in gene_ids
            if isinstance(summ.get(uid), dict)
            and str(summ[uid].get("name", "")).upper() == symbol.upper()  # type: ignore[union-attr]
        ),
        None,
    )
    if gene_uid is None:
        return None
    rna_uids = _refseq_rna_uids(gene_uid)
    if not rna_uids:
        return None
    ssum = _esummary("nuccore", rna_uids)
    uids = ssum.get("uids", []) if isinstance(ssum, dict) else []
    accs = [
        str(ssum[u]["accessionversion"])  # type: ignore[index]
        for u in uids
        if isinstance(ssum.get(u), dict) and "accessionversion" in ssum[u]  # type: ignore[operator]
    ]
    nm = sorted(a for a in accs if a.startswith("NM_"))
    xm = sorted(a for a in accs if a.startswith("XM_"))
    if nm:
        return nm[0]
    if xm:
        return xm[0]
    return None


def fetch_genbank(accessions: list[str]) -> str:
    ids = ",".join(accessions)
    url = f"{EUTILS}/efetch.fcgi?db=nuccore&rettype=gb&retmode=text&tool=bt5&id={ids}"
    return _throttled_get(url)


def cds_codons(record: SeqRecord) -> tuple[str, str, list[str]] | None:
    """Return (accession, gene symbol, in-frame codons minus the terminal stop).

    Rejects any record whose CDS is not a clean multiple of three, contains an
    ambiguous base, or does not end in a stop -- a codon count is only evidence if the
    reading frame is real. The gene is read from the CDS ``/gene`` (fallback: the gene
    feature) so the caller can verify the transcript belongs to the requested symbol.
    """
    cds = next((f for f in record.features if f.type == "CDS"), None)
    if cds is None:
        return None
    gene = ""
    for feat in (cds, *(f for f in record.features if f.type == "gene")):
        vals = feat.qualifiers.get("gene")
        if vals:
            gene = str(vals[0])
            break
    seq = str(cds.location.extract(record).seq).upper()
    if len(seq) % 3 or set(seq) - set("ACGT"):
        return None
    codons = [seq[i : i + 3] for i in range(0, len(seq), 3)]
    if len(codons) < 2:
        return None
    std = CodonTable.unambiguous_dna_by_id[1]
    if codons[-1] not in std.stop_codons:
        return None
    body = codons[:-1]
    if any(c in std.stop_codons for c in body):
        return None
    return record.id, gene, body


def build_w(counts: dict[str, int]) -> dict[str, float]:
    """Relative adaptiveness, identical to CodonUsage.from_counts (pseudocount 0.5)."""
    std = CodonTable.unambiguous_dna_by_id[1]
    families: dict[str, list[str]] = {}
    for codon, aa in sorted(std.forward_table.items()):
        if codon in std.stop_codons:
            continue
        families.setdefault(aa, []).append(codon)
    w: dict[str, float] = {}
    for codons in families.values():
        adjusted = {c: counts.get(c, 0) + 0.5 for c in codons}
        peak = max(adjusted.values())
        for c, v in adjusted.items():
            w[c] = round(v / peak, 6)
    return w


def main() -> None:
    host_key, taxid, organism, out_path = sys.argv[1:5]

    # 1. Resolve each panel symbol to a candidate accession (gene-authoritative).
    symbol_of: dict[str, str] = {}  # symbol -> category
    candidate: dict[str, str] = {}  # symbol -> accession
    unresolved: list[str] = []
    for category, symbols in PANEL.items():
        for symbol in symbols:
            symbol_of[symbol] = category
            acc = resolve_accession(symbol, taxid)
            if acc:
                candidate[symbol] = acc
                print(f"  {symbol:10s} -> {acc}", file=sys.stderr)
            else:
                unresolved.append(symbol)
                print(f"  {symbol:10s} -> (no RefSeq gene/transcript)", file=sys.stderr)

    # 2. Fetch each unique candidate accession once, then verify gene identity.
    unique_accs = sorted(set(candidate.values()))
    gb_text = fetch_genbank(unique_accs)
    records = list(SeqIO.parse(StringIO(gb_text), "genbank"))
    parsed: dict[str, tuple[str, list[str]]] = {}  # accession -> (gene, codons)
    rejected: list[str] = []
    for record in records:
        result = cds_codons(record)
        if result is None:
            rejected.append(record.id)
            continue
        acc, gene, body = result
        parsed[acc] = (gene, body)
        parsed[acc.split(".")[0]] = (gene, body)

    # 3. Count each verified accession ONCE. A symbol whose transcript's /gene does not
    #    match the symbol is a mismatch: dropped and disclosed, never counted.
    counts: Counter[str] = Counter()
    used: list[dict[str, object]] = []
    mismatched: list[dict[str, str]] = []
    counted_accs: set[str] = set()
    for symbol, acc in sorted(candidate.items()):
        entry = parsed.get(acc) or parsed.get(acc.split(".")[0])
        if entry is None:
            rejected.append(acc)
            continue
        gene, body = entry
        if gene.upper() != symbol.upper():
            mismatched.append({"symbol": symbol, "accession": acc, "actual_gene": gene})
            print(f"  MISMATCH {symbol} -> {acc} is {gene}", file=sys.stderr)
            continue
        if acc in counted_accs:
            continue
        counted_accs.add(acc)
        counts.update(body)
        used.append(
            {
                "symbol": symbol,
                "category": symbol_of[symbol],
                "accession": acc,
                "gene": gene,
                "codons": len(body),
            }
        )

    w = build_w(dict(counts))
    total_codons = sum(counts.values())
    src_sha = hashlib.sha256(gb_text.encode("utf-8")).hexdigest()

    doc = {
        "_provenance": {
            "name": f"{organism} highly-expressed codon usage (relative adaptiveness w)",
            "organism": organism,
            "ncbi_taxid": taxid,
            "reference_gene_set": (
                "Canonical highly-expressed panel: cytosolic ribosomal proteins "
                "(large + small subunit), translation elongation factors, core "
                "glycolytic enzymes, and chaperonins -- an extension of the fallback "
                "reference set of docs/research/brief.md:282 (ribosomal proteins + "
                "elongation factors + GAPDH + chaperonins), widened with further core "
                "glycolytic and other high-abundance genes; resolved to RefSeq transcripts."
            ),
            "method": (
                "Each gene symbol was resolved gene-first: the NCBI Gene record whose "
                "OFFICIAL symbol equals the symbol (never an alias) was found, then linked "
                "to its RefSeq RNA (curated NM_ preferred over predicted XM_). Every "
                "fetched transcript's CDS /gene was re-verified to equal the requested "
                "symbol before counting; a mismatch is dropped and listed in "
                "mismatched_symbols, never mapped to a paralog. In-frame sense codons were "
                "counted (terminal stop and any non-ACGT CDS excluded), each accession "
                "once, and w = (count + 0.5) / family_max per family, identical to "
                "bt5.codon.tables.CodonUsage.from_counts."
            ),
            "citation": "O'Leary NA et al. RefSeq. Nucleic Acids Res 2016;44(D1):D733-45",
            "citation_url": "https://pubmed.ncbi.nlm.nih.gov/26553804/",
            "source": "NCBI RefSeq via E-utilities (esearch/esummary/elink + efetch, db=nuccore)",
            "source_url": "https://www.ncbi.nlm.nih.gov/refseq/",
            "license": "Public domain (U.S. Government work; NCBI RefSeq is not copyrighted)",
            "retrieved": date.today().isoformat(),
            "genes_requested": sum(len(v) for v in PANEL.values()),
            "genes_contributing": len(used),
            "curated_nm_transcripts": sum(1 for r in used if str(r["accession"]).startswith("NM_")),
            "predicted_xm_transcripts": sum(
                1 for r in used if str(r["accession"]).startswith("XM_")
            ),
            "codons_counted": total_codons,
            "sense_codons": len(w),
            "source_sha256": src_sha,
            "genetic_code_table": 1,
            "unresolved_symbols": unresolved,
            "mismatched_symbols": mismatched,
            "rejected_records": rejected,
            "contributing_accessions": sorted(
                used, key=lambda r: (str(r["category"]), str(r["symbol"]))
            ),
            "note": (
                "61 sense codons; the three stops are excluded because CAI is a "
                "geometric mean over sense codons only. ATG and TGG are 1.0 by "
                "definition. This is a highly-expressed reference set, NOT a genomic "
                "average: swapping it changes every CAI value."
            ),
        },
        "w": w,
    }
    with open(out_path, "w") as fh:
        json.dump(doc, fh, indent=2)
        fh.write("\n")
    print(
        f"{host_key}: {len(used)} genes, {total_codons} codons, "
        f"{len(unresolved)} unresolved, {len(mismatched)} mismatched, "
        f"{len(rejected)} rejected -> {out_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
