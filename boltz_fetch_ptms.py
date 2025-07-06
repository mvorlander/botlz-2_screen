#!/usr/bin/env python
"""
uniprot_helper.py

CLI utility to fetch curated PTM sites from the EBI Proteins PTM
endpoint and emit wrapper-ready tokens, e.g.:

    $ python uniprot_helper.py Q13838
    Q13838 38:Phosphoserine
    Q13838 123:Phosphothreonine

Options let you filter by PTM type (substring match) and cap the number
of sites to avoid combinatorial explosions.
"""

from __future__ import annotations
from typing import Iterable, List
import sys, argparse, json, requests
import csv, pathlib

# --------------------------------------------------------------------------- #
# Constants                                                                   #
# --------------------------------------------------------------------------- #
_PTMS_ENDPOINT = (
    "https://www.ebi.ac.uk/proteins/api/proteomics/ptm/{acc}?format=json"
)
_PTM_KEYS = {"MOD_RES", "LIPID", "CARBOHYD", "DISULFID"}  # UniProt feature types

# --------------------------------------------------------------------------- #
# CCD code mapping & sequence endpoint                                        #
# --------------------------------------------------------------------------- #
_CCD_MAP = {
    "phosphoserine":    "SEP",
    "phosphothreonine": "TPO",
    "phosphotyrosine":  "PTR",
    "n6-acetyllysine":  "ALY",
    "selenocysteine":   "SEC",
    "selenium-methionine": "MSE",
}

# --------------------------------------------------------------------------- #
# PTM description + residue → CCD code mapping                                #
# --------------------------------------------------------------------------- #
def _resolve_ccd(description: str, residue: str) -> str:
    """
    Map a human‑readable PTM *description* plus the one‑letter *residue*
    into a CCD code that Boltz‑2 understands.  Falls back to a cleaned
    description if the PTM is unknown.

    Parameters
    ----------
    description : str
        PTM description from the API (e.g. "Phosphoserine", "Acetylation").
    residue : str
        One‑letter residue at that position (uppercase).

    Returns
    -------
    str
        CCD code (e.g. "SEP") or the cleaned description with spaces→_
    """
    d = description.lower()
    if "phospho" in d:
        return {"S": "SEP", "T": "TPO", "Y": "PTR"}.get(residue.upper(), "PHOS")
    if "acetyl" in d and residue.upper() == "K":
        return "ALY"
    # fall back to explicit mapping table or sanitised description
    return _CCD_MAP.get(d, description.replace(" ", "_"))

# --------------------------------------------------------------------------- #
# Helper: Format EBI syntax fragment for TSV column 1                         #
# --------------------------------------------------------------------------- #
def _format_ebi(peptide_seq: str, ptm_type: str, rel_pos: int) -> str:
    """Human‑readable EBI syntax fragment for TSV column 1."""
    return f"{peptide_seq}:{rel_pos}:{ptm_type}"

# --------------------------------------------------------------------------- #
# Helper: Robust int conversion for nested JSON values                        #
# --------------------------------------------------------------------------- #
def _to_int(x):
    """Robustly convert nested JSON number representations to int or None."""
    if x is None:
        return None
    if isinstance(x, int):
        return x
    if isinstance(x, str) and x.isdigit():
        return int(x)
    if isinstance(x, dict):
        return _to_int(x.get("value") or x.get("position") or x.get("begin") or x.get("start"))
    return None
_SEQ_ENDPOINT = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"

from functools import lru_cache

@lru_cache(maxsize=256)
def _fetch_sequence(acc: str) -> str:
    """Return canonical UniProt sequence or '' on error."""
    try:
        r = requests.get(_SEQ_ENDPOINT.format(acc=acc), timeout=15)
        r.raise_for_status()
        return "".join(r.text.splitlines()[1:])  # drop FASTA header
    except Exception:
        return ""
# --------------------------------------------------------------------------- #
# Core helper                                                                 #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# Core helper: fetch_ptm_records                                             #
# --------------------------------------------------------------------------- #
def fetch_ptm_records(
    acc: str,
    *,
    keep_types: Iterable[str] | None = None,
    cap: int | None = None,
    debug: bool = False,
) -> List[tuple]:
    """
    Return a list of PTM records for *acc*.

    Each record is a tuple: (ebi_syntax, boltz_notation, source)

    Parameters
    ----------
    acc : str
        UniProt accession (e.g. "Q13838").
    keep_types : Iterable[str] | None
        Case-insensitive substrings; keep only PTMs whose description
        contains **any** of them.  None ➜ no filter.
    cap : int | None
        Stop after *cap* sites (safety valve).

    Returns
    -------
    list[tuple]
        Each tuple: (ebi_syntax, boltz_notation, source)
    """
    dbg: list[str] = [] if debug else None
    url = _PTMS_ENDPOINT.format(acc=acc)
    r = requests.get(url, headers={"Accept": "application/json"}, timeout=15)
    r.raise_for_status()
    data = r.json()
    if dbg is not None:
        dbg.append(f"# Raw JSON top‑level keys: {list(data[0].keys() if isinstance(data, list) else data.keys())}")

    features: list[dict] = []
    if isinstance(data, dict):
        features.extend(data.get("features", []))
    elif isinstance(data, list):
        for entry in data:
            if isinstance(entry, dict):
                features.extend(entry.get("features", []))

    if not features:
        return []

    wanted = {k.lower() for k in keep_types} if keep_types else None
    records: list[tuple] = []
    seen: set = set()
    seq = _fetch_sequence(acc)
    for feat in features:
        if not isinstance(feat, dict):
            continue
        nested = feat.get("ptms")
        if not isinstance(nested, list):
            nested = []
        nested.append(feat)
        for ptm in nested:
            if not isinstance(ptm, dict):
                continue
            descr = (
                ptm.get("name")
                or ptm.get("modification")
                or feat.get("description")
                or ""
            )
            if not descr or (wanted and not any(k in descr.lower() for k in wanted)):
                continue
            rel = _to_int(ptm.get("position") or ptm.get("begin") or ptm.get("start"))
            if rel is None:
                continue
            pep_info = feat.get("peptide")
            if not isinstance(pep_info, dict):
                pep_info = {}
            pep_start = _to_int(
                (pep_info.get("start") if isinstance(pep_info, dict) else None)
                or (pep_info.get("begin") if isinstance(pep_info, dict) else None)
                or feat.get("begin")
                or feat.get("start")
            )
            abs_pos = pep_start + rel - 1 if pep_start else rel
            pep_seq = pep_info.get("sequence") if isinstance(pep_info, dict) else None
            if pep_seq and rel <= len(pep_seq):
                res = pep_seq[rel - 1]
            elif seq and 1 <= abs_pos <= len(seq):
                res = seq[abs_pos - 1]
            else:
                res = "X"
            if dbg is not None:
                dbg.append(
                    f"{descr!r:<20} rel={rel:<4} pep_start={pep_start} "
                    f"abs={abs_pos:<4} res={res}  → { _resolve_ccd(descr, res) }"
                )
            ccd = _resolve_ccd(descr, res)
            peptide_seq = pep_seq or ""
            ebi_str = _format_ebi(peptide_seq, descr, rel)
            bolt_str = f"{acc} {abs_pos}:{ccd}"
            src_str  = ",".join(ptm.get("sources", []))
            record   = (ebi_str, bolt_str, src_str)
            dedup_key = (bolt_str, src_str)
            if dedup_key not in seen:
                records.append(record)
                seen.add(dedup_key)
                if cap and len(records) >= cap:
                    if dbg is not None:
                        print("\n".join(dbg), file=sys.stderr)
                    return records
    if dbg is not None:
        print("\n".join(dbg), file=sys.stderr)
    return records

# --------------------------------------------------------------------------- #
# Command-line interface                                                      #
# --------------------------------------------------------------------------- #
def _cli() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Fetch PTM sites for one or many UniProt accessions and write:\n"
            "  • a TSV  (EBI-syntax | Boltz-notation | data source)\n"
            "  • a plain list of Boltz tokens (default <output>.boltz.txt)\n"
            "Use -t/--ptm-types to filter, --verbose to echo the table."
        )
    )
    ap.add_argument(
        "accession",
        nargs="?",
        help="UniProt accession(s), comma-separated (e.g. Q13838,P05067)",
    )
    ap.add_argument(
        "-f", "--file",
        type=str,
        help="File with UniProt IDs (one per line)."
    )
    ap.add_argument(
        "-n",
        "--max-sites",
        type=int,
        metavar="N",
        help="Stop after N PTM sites (safety valve)",
    )
    ap.add_argument(
        "-o",
        "--output",
        type=str,
        default="ptms.tsv",
        help="TSV output filename (default: ptms.tsv)",
    )
    ap.add_argument(
        "--verbose",
        action="store_true",
        help="Also print TSV output to stdout",
    )
    ap.add_argument("--debug", action="store_true",
                    help="Print raw feature→token mapping to stderr")
    ap.add_argument(
        "--ptm-types", "-t",
        metavar="PTM",
        nargs="+",
        help=(
            "Limit the fetch to specific PTM *keywords* (case-insensitive, space-separated).\n"
            "Recognised keywords include the most common UniProt annotations:\n"
            "  • phospho    (Phosphorylation)\n"
            "  • acetyl     (Lysine acetylation / N-terminal acetylation)\n"
            "  • sumo       (SUMOylation)\n"
            "  • ubiquitin  (Ubiquitinylation)\n"
            "  • glyco      (N-/O-linked Glycosylation)\n"
            "  • methyl     (Lys/Arg methylation)\n"
            "  • lipid      (Lipidation, e.g. myristoyl, palmitoyl)\n"
            "\n"
            "Examples:\n"
            "  # phospho sites only\n"
            "  uniprot_helper.py Q13838 -t phospho\n"
            "  \n"
            "  # phospho + acetyl\n"
            "  uniprot_helper.py Q13838 -t phospho acetyl\n"
            "  \n"
            "  # batch mode from file, but keep only ubiquitin & SUMO marks\n"
            "  uniprot_helper.py -f ids.txt -t ubiquitin sumo\n"
        )
    ),
    ap.add_argument(
        "--boltz-output", type=str,
        help="Write a second file containing ONLY the Boltz tokens; "
            "default is <output>.boltz.txt")
    args = ap.parse_args()


    print("⏳  Starting UniProt PTM fetch …", file=sys.stderr)
    # Gather accessions from file and/or arg
    accessions = set()
    if args.accession:
        #print(f"  • {acc}", file=sys.stderr)
        for acc in args.accession.split(","):
            acc = acc.strip()
            if acc:
                accessions.add(acc)
    if args.file:
        with open(args.file) as fh:
            for line in fh:
                acc = line.strip()
                if acc and not acc.startswith("#"):
                    accessions.add(acc)
    if not accessions:
        ap.error("No accessions provided. Use positional arg or -f/--file.")

    all_records = []
    try:
        for acc in accessions:
            print(f"  • {acc}", file=sys.stderr)
            recs = fetch_ptm_records(
                acc,
                keep_types=args.ptm_types,
                cap=args.max_sites,
                debug=args.debug,
            )
            all_records.extend(recs)
    except Exception as exc:
        sys.exit(f"❌  {exc}")

    # Write TSV
    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["ebi_syntax", "boltz_notation", "source"])
        writer.writerows(all_records)
        # ---------- write plain-text Boltz tokens ---------------------------
        boltz_path = pathlib.Path(args.boltz_output
                                if args.boltz_output
                                else args.output).with_suffix(".boltz.txt")

        with boltz_path.open("w") as fh_b:
            for _, boltz_not, _ in all_records:
                fh_b.write(boltz_not + "\n")

        print(f"✅  Wrote {len(all_records)} PTMs to {args.output}", file=sys.stderr)
        print(f"✅  Wrote Boltz tokens   to {boltz_path}", file=sys.stderr)
        
        
    if args.verbose:
        writer = csv.writer(sys.stdout, delimiter="\t")
        writer.writerows(all_records)

if __name__ == "__main__":
    _cli()