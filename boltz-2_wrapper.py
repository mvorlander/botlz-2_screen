#!/usr/bin/env python3
"""
boltz-2_wrapper.py  –  Run Boltz‑2 locally *or* via Slurm
=========================================================

This wrapper now supports **all Boltz‑2 prediction modes** and the new
*inline PTM* notation that works for both *bait* and *screen* sequences.

Supported prediction flavours
-----------------------------
* **Structure** (default)  
* **Affinity**  (automatic insertion with `--affinity L`)  
* **Diffusion with potentials** (`--boltz "--use_potentials"`)  
* **High‑throughput screening** (`--bait / --screen`)  

Inline PTM syntax
-----------------
Add one or more *position:CCD* tokens **after** the sequence token:

* `Q13838 38:SEP` → Ser‑38 is phosphorylated (SEP) in chain *Q13838*  
* `P05067 15:CSO 42:MSE` → two modifications in the same chain  
* When several chains are listed, a bare `pos:CCD` applies to the **nearest
  preceding** sequence token.  

Ligands / ions
--------------
Any 3‑letter all‑caps token is treated as a CCD ligand or ion.
SMILES strings or `.smi` files (`SMILES  NAME`) are also supported.

Quick Examples
--------------
Structure on Slurm  
    ./boltz-2_wrapper.py  complex.fa

Structure on the current GPU node  
    ./boltz-2_wrapper.py  P69905,P68871  --local

Affinity (ligand chain “L”)  
    ./boltz-2_wrapper.py  complex.fa --affinity L

Potential‑steered, 300 sampling steps  
    ./boltz-2_wrapper.py  complex.yaml  --local \\
          --boltz "--use_potentials --sampling_steps 300"

Screening examples
------------------
**What you need**

1. *Bait specification*  (`--bait …`)
   • A **single file** that lists every
     *constant* entity you want in all complexes:

     ─ UniProt IDs  
     ─ paths to FASTA / YAML / SMILES files  
     ─ ligand / ion CCD tokens (three-letter)  
     ─ optional PTMs written inline, e.g. `Q13838 38:SEP`


IMPORTANT: ordering of  the bait file
--------------------------------
For every job Boltz numbers chains strictly in the order they appear
inside the generated YAML. Therefore within the bait.txt, always order entities like this

    ▸ proteins / RNA / DNA  (in the exact order you want)
    ▸ ligands  (ATP, GTP, …)
    ▸ simple ions  (MG, ZN, CA ...)
 

2. *Screen specification*  (`--screen …`)
   • A plain-text file containing one target per line*
   • Each line may be a UniProt ID, FASTA path, or raw sequence and may
     include PTMs the same way: `Q86W42 101:SEP`.

Chain-map file
--------------
To give human-readable labels in the analysis stage, create a plain-text
*chain-map* file that assigns names to your prediction targets. ATTENTION:  
For the mapping file, keep in mind that the chain ID of the bait entity will
depend on its type; if is a protein or nucelic aid type entity, and your bait 
contains both protein/NAs AND ligands, the screen wntity will be inserted BETWEEN
the bait protein/NA entities and any Ligand entities. 


Concrete examples
-----------------
Example 1 –  Screening a pre-defined bait complex of protein, RNA, ATP, and MG ion against proteins and RNAs.     

  contents of bait.txt (Keep this order!):
    Q13838          #uniprot of bait protein
    ./RNA.fa        #path to file containing bait RNA sequence 
    ATP             #Ligand, only CCD ligand codes are accepted
    MG              #MG ion

    
  contents of screen.txt:
    Q09161          #uniprot screen target. 1
    Q13838 38:SEP   #uniprot screen target, 2 phosphorylated on Ser2
    /RNA2.fa        #path to file containing screen RNA sequence 
    
   contents of chain_mapping.txt 
    0=UAP56         #name of bait protein listed first in bait.txt
    1=baitRNA       #name of bait protein listed first in bait.txt
    2=target        #placeholder for the target that will be screened, since we are screening protein/RNA types, the final prediction will insert the chain BEFORE any ligands in the prediction, so it gets position 3
    3=ATP           #ligands before ions
    4=Mg2+          #ions last
    
  Submission command:
    ./boltz-2_wrapper.py --bait bait.txt --screen screen.txt --chain_mapping chain_mapping.txt --n MyScreenName

Example 2 – A complex of two bait proteins (second protein phosphorylated on Ser 38) versus 2 RNAs
  contents of  bait_proteins.txt:
    Q09161
    Q13838 38:SEP
    
  contents of screen_rnas.txt:
    /data/RNA1.fa
    /data/RNA2.fa
  
  contents of chain_mapping.txt 
    0=UAP56         #name of bait protein listed first in bait.txt
    1=baitRNA       #name of bait protein listed second in bait.txt
    2=target        #placeholder for the target that will be screened, since we are screening protein/RNA types, the final prediction will insert the chain BEFORE any ligands in the prediction, so it gets position 3

  Submission command
    ./boltz-2_wrapper.py --bait bait_proteins.txt --screen screen_rnas.txt -n prot_vs_rna

Example 3 – protein + RNA + bait versus different ligands
  contents of bait_combo.list:
    P69905
    bait_rna.fa
    MG
    
  contents of  screen_proteins.txt
    ATP
    ADP
    AMP
    
  contents of chain_mapping.txt 
    0=UAP56         #name of bait protein listed first in bait.txt
    1=baitRNA       #name of bait protein listed second in bait.txt
    2=targetLigand  #placeholder for the target that will be screened, since we are screening protein/RNA types, the final prediction will insert the chain AFTER any proteins.RNAs in the prediction, so it gets position 3 (index 2)
    3=MG            #Ion form bait list

---------------------------------------------------------------------------
"""

import argparse, io, pathlib, re, shlex, shutil, subprocess, sys, tempfile, textwrap, time
from concurrent.futures import ThreadPoolExecutor, as_completed
import yaml, requests
from requests.exceptions import RequestException
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from datetime import datetime
import json       # ← add this line
import re   # already imported once, keep only one import

def _clean_seq(seq: str) -> str:
    """Return *seq* stripped of ALL whitespace & line-breaks."""
    return re.sub(r"\s+", "", seq)

# ------------------------------------------------------------------
# Absolute directory where this wrapper script resides.
# Lets us reference boltz_analysis.py regardless of CWD.
# ------------------------------------------------------------------
WRAPPER_DIR = pathlib.Path(__file__).resolve().parent# 
#--------------------------------------------------------------------------- #
# Helpers                                                                     #
# --------------------------------------------------------------------------- #

# ------------------------------------------------------------------
# Cluster GPU inventory (Plaschka lab, July 2025)
# ------------------------------------------------------------------
# key  = Slurm --constraint value
# cap  = usable VRAM in **GB** (before safety deduction)
#   g1 = P100‑12 GB     g2 = V100‑32 GB
#   g3 = RTX 6000‑24 GB g4 = A100‑40 GB
_GPU_CAPACITY = {
    "g1": 12,
    "g3": 24,
    "g2": 32,
    "g4": 40,
}
# Small safety cushion so borderline cases are skipped instead of OOMing
_SAFETY_MARGIN_GB = 2      # e.g. 40‑2 = 38 GB max for A100‑40 GB


# --------------------------------------------------------------------------- #
# FASTA validation helper                                                     #
# --------------------------------------------------------------------------- #
FASTA_EXT = {".fa", ".fasta", ".faa", ".rna", ".dna"}

def _validate_fasta(path: pathlib.Path) -> None:
    """
    Ensure *path* exists and is a well-formed FASTA.

    • file must exist  
    • first non-blank line must start with '>'  

    The script aborts with a clear message if the check fails.
    """
    if not path.is_file():
        sys.exit(f"❌ FASTA file not found: {path}")
        
        
    with path.open() as fh:
        for ln in fh:
            if ln.strip():                     # first non-empty line
                if not ln.startswith(">"):
                    sys.exit(f"❌ {path} is not a valid FASTA – "
                             "missing leading '>' header.")
                break
# ---------------------------------------------------------------------------
# chain-map helper (shared by wrapper & analysis)
# ---------------------------------------------------------------------------
def load_chain_map(path: pathlib.Path | None) -> dict[str, str]:
    """
    Read a chain-mapping file.

    Numeric keys may be **0‑ or 1‑based**; both are accepted transparently.
    • Keys may be positional indices (1-based: \"1\", \"2\", …) **or**
      explicit tokens (UniProt accession, ligand CCD, fasta basename).
    • Values are the user-friendly names shown in plots etc.

    Lines starting with ‘#’ are ignored.
    """
    mapping: dict[str, str] = {}
    if not path or not path.is_file():
        return mapping
    with path.open() as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            if "=" in ln:
                k, v = ln.split("=", 1)
            else:
                k, v = ln.split(None, 1)
            # add upper‑case and 5‑char aliases so that truncated IDs
            # (Boltz chain‑ids are ≤ 5 chars) still find a mapping
            ku = k.strip().upper()
            mapping[ku] = v.strip()
            if len(ku) > 5:
                mapping[ku[:5]] = v.strip()
    return mapping

def total_residues(path: pathlib.Path) -> int:
    if path.suffix.lower() in {".yml", ".yaml"}:
        doc = yaml.safe_load(path.read_text())
        return sum(len(ent["sequence"])
                for item in doc["sequences"]
                for ent in item.values()
                if list(item.keys())[0] in {"protein", "rna", "dna"})
    return sum(len(r.seq) for r in SeqIO.parse(str(path), "fasta"))

def count_ligands(yaml_path: pathlib.Path) -> int:
    """Return number of ligand / ion entries in a Boltz YAML."""
    if yaml_path.suffix.lower() not in {".yml", ".yaml"}:
        return 0
    doc = yaml.safe_load(yaml_path.read_text())
    return sum(1 for item in doc.get("sequences", [])
                 if list(item.keys())[0] == "ligand")

# --------------------------------------------------------------------------- #
# YAML quick statistics (residues, chains, ligands, modifications)            #
# --------------------------------------------------------------------------- #
def yaml_stats(yaml_path: pathlib.Path) -> dict[str, int]:
    """Return basic counts used for memory / GPU heuristics."""
    stats = dict(residues=0, chains=0, ligands=0, mods=0)
    if yaml_path.suffix.lower() not in {".yml", ".yaml"}:
        return stats
    doc = yaml.safe_load(yaml_path.read_text())
    for item in doc.get("sequences", []):
        kind, seq_dict = next(iter(item.items()))
        if kind in {"protein", "rna", "dna"}:
            stats["chains"] += 1
            stats["residues"] += len(seq_dict["sequence"])
            stats["mods"] += len(seq_dict.get("modifications", []))
        elif kind == "ligand":
            stats["ligands"] += 1
    return stats

# --------------------------------------------------------------------------- #
# GPU class & extra‑flags heuristic                                           #
# --------------------------------------------------------------------------- #
#
# --------------------------------------------------------------------------- #
# GPU / host‑memory heuristic + auto‑flags                                      #
# --------------------------------------------------------------------------- #

def gpu_and_flags(residues: int, chains: int, ligands: int, mods: int):
    """
    Return
       (constraint_str, host_mem_string, extra_flags_list, est_vram_gb)

    • If the job is too big for the largest GPU we have, *constraint_str*
      is **None** – the caller can skip it.
    """
    # ---------- adaptive Boltz flags -----------------------------------
    flags: list[str] = []
    if residues >= 2500 or chains >= 6:
        flags += ["--sampling_steps", "40"]
    elif residues >= 1800 or chains >= 4:
        flags += ["--sampling_steps", "60"]
    if ligands:
        flags += ["--no_kernels"]
    if mods >= 20:
        flags += ["--recycles", "1"]

    # helpers to read our own override back
    def _par(opt: str, default: int) -> int:
        try:
            return int(flags[flags.index(opt) + 1])
        except ValueError:
            return default
    steps    = _par("--sampling_steps", 100)
    recycles = _par("--recycles", 3)

    # ---------- crude VRAM estimator (good ±15 %) ----------------------
    est_vram_gb = (
        0.0000035 * residues ** 2 +
        0.015 * residues +
        0.5 * ligands +
        2
    ) * (steps / 100) * (recycles / 3)

    # ---------- pick smallest GPU that fits ---------------------------
    if est_vram_gb <= _GPU_CAPACITY["g1"] - _SAFETY_MARGIN_GB:
        gpu = "g1"; host_mem_gb = 70
    elif est_vram_gb <= _GPU_CAPACITY["g3"] - _SAFETY_MARGIN_GB:
        gpu = "g3"; host_mem_gb = 90
    elif est_vram_gb <= _GPU_CAPACITY["g2"] - _SAFETY_MARGIN_GB:
        gpu = "g2"; host_mem_gb = 120
    elif est_vram_gb <= _GPU_CAPACITY["g4"] - _SAFETY_MARGIN_GB:
        gpu = "g4"; host_mem_gb = 140
    else:
        # would still OOM on A100‑40 → tell caller to skip
        return None, None, None, est_vram_gb

    # round host RAM up to next 1000 MB
    host_mem_mb = int((host_mem_gb + 1) * 1000)

    constraint_map = {
        "g1": "g4|g2|g3|g1",
        "g3": "g4|g3",
        "g2": "g4|g2",
        "g4": "g4",
    }
    return constraint_map[gpu], f"{host_mem_mb}M", flags, est_vram_gb
    
# --------------------------------------------------------------------------- #
# Sequence classification helper                                              #
# --------------------------------------------------------------------------- #

def classify_sequence(seq: str) -> str:
    """Return 'dna', 'rna', or 'protein' based on alphabet heuristic."""
    s = set(seq.upper())
    dna_lets = {"A", "C", "G", "T", "N"}
    rna_lets = {"A", "C", "G", "U", "N"}
    if s.issubset(dna_lets):
        return "dna"
    if s.issubset(rna_lets):
        return "rna"
    return "protein"

# --------------------------------------------------------------------------- #
# YAML generation                                                             #
# --------------------------------------------------------------------------- #

# ---------------- FASTA‑header token extraction --------------------------- #
def _header_token(header: str) -> str:
    """
    Return a reasonable accession/token from a FASTA header.

    * UniProt convention “sp|Q9VKM1|PIWI_DROME …”  →  Q9VKM1
    * Anything else         →  first part before whitespace or ‘|’
    """
    parts = header.split("|")
    if len(parts) >= 3 and parts[0] in {"sp", "tr"}:
        return parts[1]
    return re.split(r"[\\s|:_-]+", header, 1)[0]

UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{}.fasta"

def _fetch_single(uid: str, timeout: int):
    """Return (uid, text_or_None) for a single UniProt ID."""
    try:
        r = requests.get(UNIPROT_FASTA_URL.format(uid), timeout=timeout)
        if r.ok:
            return uid, r.text
    except RequestException:
        pass
    return uid, None

def _records_from_uniprot(ids, *, timeout: int = 10, retries: int = 2, workers: int = 4):
    """
    Parallel download of UniProt FASTA records with progress, timeout and retry.
    """
    ids = [uid.strip() for uid in ids if uid.strip()]
    if not ids:
        return

    remaining = ids
    round_n = 1
    while remaining and round_n <= retries:
        print(f"⏳  UniProt round {round_n}/{retries} for {len(remaining)} ID(s)…", flush=True)
        new_remaining = []
        with ThreadPoolExecutor(max_workers=min(workers, len(remaining))) as pool:
            futures = {pool.submit(_fetch_single, uid, timeout): uid for uid in remaining}
            for fut in as_completed(futures):
                uid, fasta = fut.result()
                if fasta:
                    print(f"  ✓ {uid}")
                    yield from SeqIO.parse(io.StringIO(fasta), "fasta")
                else:
                    print(f"  ✗ {uid}")
                    new_remaining.append(uid)
        remaining = new_remaining
        round_n += 1

    if remaining:
        sys.stderr.write("⚠ UniProt fetch failed for: " + ", ".join(remaining) + "\n")

# --------------------------------------------------------------------------- #
# Ligand/ion token detection helper                                           #
# --------------------------------------------------------------------------- #
def _token_to_entry(token: str):
    """
    Return tuple(kind, entry_dict) based on token heuristics.

    * protein/rna/dna  -> handled later via sequence
    * 3‑letter all‑caps → CCD ligand
    * SMILES string    → ligand with smiles
    * .smi file path   → each line 'SMILES  NAME'
    """
    p = pathlib.Path(token)
    lig_entries = []
    if p.is_file() and p.suffix.lower() in {".smi", ".smiles"}:
        for line in p.read_text().splitlines():
            if not line.strip(): continue
            parts = line.strip().split()
            smi, name = parts[0], (parts[1] if len(parts) > 1 else f"LIG{len(lig_entries)+1}")
            lig_entries.append(("ligand", {"id": name, "smiles": smi}))
        return lig_entries   # special multi‑return
    # recognise 1‑to‑3‑letter (or digit) CCD/ion codes, e.g. MG, NA, ZN2
    if re.fullmatch(r"[A-Z0-9]{1,3}", token):
        return [("ligand", {"id": token, "ccd": token})]
    # very simple SMILES heuristic (contains lowercase + atoms)
    if any(c in token for c in "=#") or re.search(r"[a-z]", token):
        return [("ligand", {"id": f"LIG{token[:3]}", "smiles": token})]
    return []  # not a ligand token

def make_yaml(src: str,
              msa_paths: list[pathlib.Path],
              chain_names: dict[str, str]) -> pathlib.Path:
    """
    Build YAML from mixed inputs: UniProt IDs, FASTA files, or a single YAML.
    Multiple items can be comma- or whitespace‑separated.
    """
    # If the string contains a comma / whitespace it is surely a mixed token list
    if re.search(r"[,\s]", src):
        tokens = re.split(r"[,\s]+", src.strip())
    else:
        p = pathlib.Path(src)
        try:
            if p.is_file():
                # Treat sequence / SMILES files as single tokens
                suf = p.suffix.lower()
                if suf in FASTA_EXT.union({".fastq"}):
                    _validate_fasta(p)
                    tokens = [str(p)]
                elif suf in {".smi", ".smiles"}:
                    tokens = [str(p)]
                else:
                    tokens = [src]
            else:
                tokens = [src]
        except OSError:
            # extremely long token cannot be treated as a file path
            tokens = [src]

    print(f"[{datetime.now():%H:%M:%S}] Building YAML from mixed source “{src}”")
    records = []
    seq_expect = 0          # number of sequences we expect to parse
    extra_ligands = []
    modifications_raw = []      # tuples (chain_index, position, CCD)

    for tok in tokens:
        if not tok:
            continue
        # -------------------------------------------------------------- #
        # Residue‑level modification tokens                              #
        #   • chain:pos:CCD   → applies to explicit chain index (legacy) #
        #   • pos:CCD         → applies to the most‑recent sequence      #
        # -------------------------------------------------------------- #
        m_mod_long  = re.fullmatch(r"(?i)(\d+):(\d+):([A-Z0-9]{3})", tok)
        m_mod_short = re.fullmatch(r"(?i)(\d+):([A-Z0-9]{3})", tok)

        if m_mod_long:
            chain_idx = int(m_mod_long.group(1))      # 1‑based index
            mod_pos   = int(m_mod_long.group(2))
            mod_ccd   = m_mod_long.group(3).upper()
            modifications_raw.append((chain_idx, mod_pos, mod_ccd))
            continue

        if m_mod_short:
            if not records:
                print(f"⚠  Ignoring modification '{tok}' – no preceding sequence token.")
                continue
            chain_idx = len(records)                  # last parsed sequence
            mod_pos   = int(m_mod_short.group(1))
            mod_ccd   = m_mod_short.group(2).upper()
            modifications_raw.append((chain_idx, mod_pos, mod_ccd))
            continue
        path = pathlib.Path(tok)
        if path.is_file() and path.suffix.lower() in FASTA_EXT:
            _validate_fasta(path)
            fa_recs = list(SeqIO.parse(str(path), "fasta"))
            seq_expect += len(fa_recs)
            records.extend(fa_recs)
        else:
            # ligand / ion tokens
            lig_try = _token_to_entry(tok)
            if lig_try:
                extra_ligands.extend(lig_try)
                continue

            # skip lone FASTA headers
            if tok.startswith(">"):
                continue

            # raw sequence pasted inline → create a SeqRecord
            if re.fullmatch(r"[ACGUTNRYSWKMBDHV]+", tok.upper()) and len(tok) > 20:
                seq_id = f"SEQ{len(records)+1}"
                seq_expect += 1
                records.append(SeqRecord(Seq(tok), id=seq_id, description=""))
                continue
            fetched = list(_records_from_uniprot([tok]))
            seq_expect += len(fetched)
            records.extend(fetched)

    if not records and not extra_ligands:
        sys.exit("❌ No valid records found from input")

    if len(records) != seq_expect:
        sys.exit(f"❌ Parsed {len(records)} sequences but expected "
                 f"{seq_expect}. Check FASTA headers or file paths.")

    fixed, ids_in_order, used, types = [], [], set(), []

    for rec in records:
        kind = classify_sequence(str(rec.seq))
        raw_token = _header_token(rec.id)
        mapped = chain_names.get(raw_token.upper()) if chain_names else None

        # ---------- friendly & chain IDs ---------------------------------
        friendly_name = mapped if mapped else raw_token  # full descriptive label

        token_clean   = re.sub(r"\W+", "", raw_token)
        base_token    = token_clean.upper() or f"SEQ{len(used)+1}"
        cid_base      = base_token[:5]                 # ✱ Boltz requires ≤ 5 chars
        cid           = cid_base
        n = 1
        while cid in used:
            cid = f"{cid_base[:4]}{n}"  # ensure uniqueness but keep ≤5
            n  += 1
        used.add(cid)

        # rewrite record header (5‑char id | kind) so MSAs match
        rec.id, rec.description = f"{cid}|{kind}", ""
        fixed.append(rec)
        ids_in_order.append(cid)
        types.append(kind)
        # store friendly mapping for later use when writing YAML
        rec.friendly = friendly_name

    protein_count = sum(1 for t in types if t == "protein")
    if msa_paths and len(msa_paths) != protein_count:
        sys.exit("❌  --msa count must match number of protein chains")

    yaml_struct = {"version": 1, "sequences": []}
    for pos, (cid, rec) in enumerate(zip(ids_in_order, fixed), 1):
        raw_token = rec.id.split("|")[0]          # original identifier
        # precedence: explicit token (full or 5‑char)  →  positional index
        explicit_key = raw_token.upper()
        explicit_key5 = explicit_key[:5]
        friendly = (
            chain_names.get(explicit_key)        # full token
            or chain_names.get(explicit_key5)    # 5‑char alias
            or chain_names.get(str(pos))         # 1‑based index
            or chain_names.get(str(pos - 1))     # 0‑based index
            or raw_token
        )
        kind = classify_sequence(str(rec.seq))          # protein / rna / dna

        entry = {"id": cid, "sequence": str(rec.seq)}
        if kind == "protein" and msa_paths:
            entry["msa"] = str(msa_paths.pop(0))

        yaml_struct["sequences"].append({kind: entry})  # keep true kind
        
    # append ligand/ion entries, applying chain-map if supplied, and keep IDs at max 5 chars
    for kind, entry in extra_ligands:
        if len(entry["id"]) > 5:
            entry["id"] = entry["id"][:5]
        yaml_struct["sequences"].append({kind: entry})

    # ------------------------------------------------------------------ #
    # per‑chain residue modifications                                    #
    # ------------------------------------------------------------------ #
    for idx, pos, ccd in modifications_raw:
        if 1 <= idx <= len(yaml_struct["sequences"]):
            seq_block = yaml_struct["sequences"][idx - 1]   # ordered entry
            key       = next(iter(seq_block))
            seq_dict  = seq_block[key]
            seq_dict.setdefault("modifications", []).append(
                {"position": pos, "ccd": ccd}
            )
        else:
            print(f"⚠  modification chain index {idx} out of range – skipped")

    # build a concise but unique filename stem that includes ligand IDs too
    ligand_ids = [entry["id"] for kind, entry in extra_ligands]
    stem_tokens = list(dict.fromkeys(ids_in_order + ligand_ids))  # preserve order, dedup
    stem_ids = _safe_name("_".join(stem_tokens[:4]), 50)
    out = pathlib.Path(f"boltz_{stem_ids}.yaml")
    counter = 1
    while out.exists():
        out = pathlib.Path(f"boltz_{stem_ids}_{counter}.yaml")
        counter += 1
    out.write_text(yaml.safe_dump(yaml_struct, sort_keys=False))
    print(f"✅  YAML written → {out.resolve()}")
    return out.resolve()

# --------------------------------------------------------------------------- #
# AlphaFold-Server JSON helpers                                               #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# AlphaFold-Server JSON converter                                             #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# AlphaFold-Server JSON converter  (stable, proven)                           #
# --------------------------------------------------------------------------- #
def yaml_to_af3(yaml_dict: dict, job_name: str) -> list[dict]:
    """
    Translate a Boltz YAML document into AlphaFold-Server JSON.

    • Returns *[job_dict]* – a list with a single job.
    • Keys match the April-2025 AlphaFold-Server README exactly.
    • Seeds are left empty so the server chooses its own.
    """
    ion_set = {
        "LI", "NA", "K", "MG", "CA", "MN", "FE", "CO", "NI",
        "CU", "ZN", "CL", "BR", "IOD", "IOD3", "IOD5"
    }

    def _ptm(entry):
        return [
            {"ptmType": f"CCD_{m['ccd'].upper()}",
             "ptmPosition": int(m["position"])}
            for m in entry.get("modifications", [])
        ]

    def _na_mods(entry):
        return [
            {"modificationType": f"CCD_{m['ccd'].upper()}",
             "basePosition": int(m["position"])}
            for m in entry.get("modifications", [])
        ]

    seqs: list[dict] = []
    for item in yaml_dict.get("sequences", []):
        kind, data = next(iter(item.items()))

        # ---------------- sequence-bearing kinds ----------------
        if kind in {"protein", "dna", "rna"}:
            seq_clean = _clean_seq(data["sequence"])
            base = {"sequence": seq_clean, "count": 1}

            if kind == "protein":
                mods = [{"ptmType":      f"CCD_{m['ccd'].upper()}",
                        "ptmPosition":  int(m["position"])}
                        for m in data.get("modifications", [])]
                if mods:
                    base["modifications"] = mods
                seqs.append({"proteinChain": base})

            elif kind == "dna":
                mods = [{"modificationType": f"CCD_{m['ccd'].upper()}",
                        "basePosition":     int(m["position"])}
                        for m in data.get("modifications", [])]
                if mods:
                    base["modifications"] = mods
                seqs.append({"dnaSequence": base})

            else:  # rna
                mods = [{"modificationType": f"CCD_{m['ccd'].upper()}",
                        "basePosition":     int(m["position"])}
                        for m in data.get("modifications", [])]
                if mods:
                    base["modifications"] = mods
                seqs.append({"rnaSequence": base})

        # ---------------- ligand / ion ------------------------------------
    

        elif kind == "ligand":
            ccd = (data.get("ccd") or data.get("id") or "").upper()
            # use ion_set defined at the top of the function
            if ccd in ion_set:
                seqs.append({"ion": {"ion": ccd, "count": 1}})
            else:
                seqs.append({"ligand": {"ligand": f"CCD_{ccd}", "count": 1}})

    job = {
        "name": job_name,
        "modelSeeds": [],   # let server generate a random seed
        "sequences": seqs,
        "dialect": "alphafoldserver",
        "version": 1
    }
    return [job]

# ---- central hub -----------------------------------------------------------
def _link_af3_json(json_path: pathlib.Path, root_dir: pathlib.Path) -> None:
    hub = root_dir / "AF3_JSON"
    hub.mkdir(exist_ok=True)
    link = hub / json_path.name
    if link.exists():
        return
    try:
        link.symlink_to(json_path.resolve())
    except (OSError, NotImplementedError):
        import shutil
        shutil.copy2(json_path, link)

# --------------------------------------------------------------------------- #
# Slurm template (Apptainer)                                                  #
# --------------------------------------------------------------------------- #

SLURM_TEMPLATE = """#!/bin/bash
#SBATCH --job-name=__JOB__
#SBATCH --output=__OUTDIR__/__JOB___%j.out
#SBATCH --error=__OUTDIR__/__JOB___%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=__PART__
#SBATCH --constraint=__CONSTRAINT__
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=__MEM__

set -euo pipefail
export PYTORCH_MATMUL_PRECISION=medium
export TORCH_MATMUL_ALLOW_TF32=1
INPUT="$1"
OUTDIR="$2"
EXTRA="$3"

echo "[run]  $(date)  GPU on $(hostname)"
echo "Input : $INPUT"
echo "Out   : $OUTDIR"

mkdir -p "$OUTDIR"

apptainer exec --nv \\
  --bind "$INPUT:$INPUT" \\
  /groups/plaschka/shared/software/boltz-2/boltz2:develop \\
  boltz predict "$INPUT" \\
        --out_dir "$OUTDIR" \\
        --accelerator gpu \\
        --use_msa_server $EXTRA
# ------------------------------------------------------------------
# Flatten output – move boltz_results_* one level up
# ------------------------------------------------------------------
RES_SUB=$(find "$OUTDIR" -maxdepth 1 -type d -name "boltz_results_*" | head -n 1 || true)
if [ -n "$RES_SUB" ]; then
    shopt -s dotglob
    mv "$RES_SUB"/* "$OUTDIR"/
    rmdir "$RES_SUB"
fi

"""

# ---------- job‑array & analysis templates -------------------------------
ARRAY_TEMPLATE = """#!/bin/bash
#SBATCH --job-name={root}
#SBATCH --array=1-{n}
# (initial temporary log files in root/slurm – will be re‑directed below)
#SBATCH --output={root}/slurm/array_%A_%a.tmp.out
#SBATCH --time=01:00:00
#SBATCH --partition=g
#SBATCH --constraint=__CONSTRAINT__
#SBATCH --gres=gpu:1
#SBATCH --mem={mem}

set -euo pipefail
export PYTORCH_MATMUL_PRECISION=medium
export TORCH_MATMUL_ALLOW_TF32=1
# ------------------------------------------------------------------
# Re‑load line “SLURM_ARRAY_TASK_ID” from jobs.list and redirect logs
# into the final per‑job directory so *.out/.err sit beside the YAML.
# ------------------------------------------------------------------
IFS=$'\\t' read -r YAML OUTDIR FLAGS MEM CONSTRAINT < <(
      sed -n "${{SLURM_ARRAY_TASK_ID}}p" {root}/jobs.list)

# Re‑direct *all* subsequent stdout / stderr to that directory
mkdir -p "${{OUTDIR}}"
exec > "${{OUTDIR}}/slurm_${{SLURM_JOB_ID}}.out" 2>&1

# Restore FLAGS placeholder
[ "${{FLAGS}}" = "-" ] && FLAGS=""

echo "[run]  $(date) on $(hostname)  — array ${{SLURM_ARRAY_TASK_ID}} / ${{SLURM_ARRAY_TASK_MAX}}"
echo "yaml  : $YAML"
echo "outdir: $OUTDIR"
echo "flags : $FLAGS"

apptainer exec --nv --bind "$YAML:$YAML" \\
  /groups/plaschka/shared/software/boltz-2/boltz2:develop \\
  boltz predict "$YAML" \\
        --out_dir "$OUTDIR" \\
        --accelerator gpu \\
        --use_msa_server $FLAGS
        
# ------------------------------------------------------------------
# Flatten output – move boltz_results_* one level up
# ------------------------------------------------------------------
RES_SUB=$(find "$OUTDIR" -maxdepth 1 -type d -name "boltz_results_*" | head -n 1 || true)
if [ -n "$RES_SUB" ]; then
    shopt -s dotglob
    mv "$RES_SUB"/* "$OUTDIR"/
    rmdir "$RES_SUB"
fi        
"""

ANALYSIS_TEMPLATE = """#!/bin/bash
#SBATCH --job-name={root}_ana
#SBATCH --output={root}/slurm/analysis_%j.out
#SBATCH --error={root}/slurm/analysis_%j.err
#SBATCH --partition=c
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --dependency=afterok:{array_id}

module purge
python {wrapper_dir}/boltz_analysis.py {root} --no-labels {chain_flag}
"""

def write_slurm(job: str, yaml_path: pathlib.Path, outdir: pathlib.Path, flags: str, mem_req: str, part: str, constraint: str) -> pathlib.Path:
    script_txt = SLURM_TEMPLATE.replace("__JOB__", job)
    script_txt = script_txt.replace("__MEM__", mem_req).replace("__PART__", part)
    script_txt = script_txt.replace("__CONSTRAINT__", constraint)
    script_txt = script_txt.replace("__OUTDIR__", str(outdir))
    script_path = outdir / f"{job}.slurm"
    with open(script_path, "w") as fh:
        fh.write(script_txt)
    print(f"📝  Slurm script written: {script_path}")
    return script_path

# --------------------------------------------------------------------------- #
# _parse_list helper                                                          #
# ----------    ----------------------------------------------------------------- #
def _parse_list(spec: str) -> list[str]:
    """
    Return a list of meaningful tokens extracted from *spec*.

    * If *spec* is a comma/whitespace‑separated string → split and clean.
    * If *spec* is a path to a FASTA / SMILES file → keep it as ONE token
      so the file’s internal records aren’t treated as separate targets.
    * If *spec* is a path to a generic *.txt list → return one cleaned
      token per non‑empty, non‑comment line.

    Comment lines start with ‘#’ after optional leading whitespace and are
    ignored.  Blank lines are ignored, too.
    """
    p = pathlib.Path(spec)

    # ---------- Direct file reference -----------------------------------
    if p.is_file():
        # Treat sequence / SMILES files as single tokens
        if p.suffix.lower() in {
            ".fa", ".fasta", ".faa", ".fastq",
            ".rna", ".dna", ".smi", ".smiles"
        }:
            return [str(p)]

        # Ordinary text list: one meaningful line ⇒ one token
        clean: list[str] = []
        for ln in p.read_text().splitlines():
            ln_clean = ln.strip()
            if not ln_clean or ln_clean.startswith("#"):
                continue          # skip blanks & comments
            clean.append(ln_clean)
        return clean

    # ---------- Comma / whitespace separated string ---------------------
    return [tok for tok in re.split(r"[,\s]+", spec) if tok]

# --------------------------------------------------------------------------- #
# sanity‑check helper                                                         #
# --------------------------------------------------------------------------- #
def _assert_token_ok(token: str) -> None:
    """
    Abort with a clear message if *token* references a file that
    does not exist or is obviously malformed/empty.

    • For paths ending in recognised FASTA/SEQ extensions – check file
      existence **and** verify a leading '>' header via _validate_fasta().
    • For YAML / YML paths – require the file to exist and contain the
      keyword 'sequences' (cheap 1‑kB read).
    • For .smi/.smiles – just check that the file exists and is non‑empty.
    Everything else (uniProt IDs, CCD tokens, raw sequences) is allowed
    through unchanged.
    """
    p = pathlib.Path(token)
    if p.is_file():
        suf = p.suffix.lower()
        if suf in FASTA_EXT.union({".fastq"}):
            _validate_fasta(p)             # will exit on error
        elif suf in {".yml", ".yaml"}:
            if p.stat().st_size == 0:
                sys.exit(f"❌ YAML file is empty: {p}")
            with p.open("r", encoding="utf-8", errors="ignore") as fh:
                head = fh.read(2048)
            if "sequences" not in head:
                sys.exit(f"❌ YAML seems malformed (no 'sequences' key): {p}")
        elif suf in {".smi", ".smiles"}:
            if p.stat().st_size == 0:
                sys.exit(f"❌ SMILES file is empty: {p}")
    elif p.expanduser().is_absolute() or "/" in token or "\\" in token:
        # looks like a path but does not exist
        sys.exit(f"❌ Referenced path not found: {token}")
# --------------------------------------------------------------------------- #
# main                                                                        #
# --------------------------------------------------------------------------- #

def _safe_name(token: str, maxlen: int = 15) -> str:
    """Return filesystem‑safe, truncated version of token."""
    return re.sub(r'[^A-Za-z0-9]+', '_', token)[:maxlen] or "X"

def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent(__doc__))
    ap.add_argument("input", nargs="?", help="FASTA / YAML / UniProt IDs (omit when using --screen)")
    ap.add_argument("--out", "-o")
    ap.add_argument("--msa", "-m", help="comma-separated list of MSA files")
    ap.add_argument("--no-kernels", action="store_true",
                    help="adds --no_kernels to boltz")
    ap.add_argument("--local", action="store_true",
                    help="Run directly on this node (no Slurm); needs a visible GPU.")
    ap.add_argument("--affinity", metavar="CHAIN_ID",
                    help="Inject affinity prediction block with the given "
                         "ligand chain‑id (e.g. L). YAML must already contain "
                         "a ligand entry with that id.")
    ap.add_argument("--boltz", default="",
                    help="Extra flags passed verbatim to `boltz predict` "
                         "(quoted as one string).")
    ap.add_argument("--bait", help="Comma-separated list or file with IDs/FASTA that stay constant in every screening pair.")
    ap.add_argument("--screen", metavar="FILE",
                    help="Text file: each line a UniProt ID or FASTA path to be paired with the bait(s). "
                         "Triggers screening mode: one job per line.")
    ap.add_argument("--screen-name", "-n",
                    help="Short label to prefix screening output/jobs. "
                         "If omitted, the wrapper uses sanitized target tokens.")
    ap.add_argument("--chain_map", "-c", metavar="FILE",
                    help="Optional mapping file: TOKEN  NEW_ID   "
                         "(one entry per line, tokens match UniProt IDs, "
                         "FASTA basenames, or ligand CCDs).")
    args = ap.parse_args()

    # ------------------------------------------------------------------
    # Chain map loading
    # ------------------------------------------------------------------
    chain_map_path = pathlib.Path(args.chain_map).expanduser() if args.chain_map else None
    chain_names: dict[str, str] = load_chain_map(chain_map_path)

    # ------------------------------------------------------------------ #
    # Sanity‑check CLI combination                                       #
    # ------------------------------------------------------------------ #
    if args.input is None and not args.screen:
        ap.error("the following arguments are required: input "
                 "(unless --screen is supplied)")

    # ---- screening mode ----------------------------------------------------
    if args.screen:
        if not args.bait:
            sys.exit("❌  --screen needs --bait")
        bait_list   = _parse_list(args.bait)
        targets     = _parse_list(args.screen)
        # ----------------------------------------------
        # Basic sanity checks before generating scripts
        # ----------------------------------------------
        if not bait_list:
            sys.exit("❌ Bait specification yielded NO usable tokens.")
        if not targets:
            sys.exit("❌ Screen specification yielded NO targets.")

        for tok in bait_list + targets:
            _assert_token_ok(tok)

        root_label  = args.screen_name or datetime.now().strftime("screen_%Y%m%d_%H%M%S")
        root_out    = pathlib.Path(args.out or f"boltz_{root_label}").resolve()
        root_out.mkdir(parents=True, exist_ok=True)

        # --- unified layout ------------------------------------------------
        inputs_dir  = root_out / "inputs"
        results_dir = root_out / "results"
        slurm_root  = root_out / "slurm"
        for _d in (inputs_dir, results_dir, slurm_root):
            _d.mkdir(parents=True, exist_ok=True)
            
        # Keep the raw input specs for provenance
        def _cp_if_file(path_spec):
            if not path_spec:
                return
            p = pathlib.Path(path_spec).expanduser()
            if p.is_file():
                try:
                    shutil.copy2(p, inputs_dir / p.name)
                except shutil.SameFileError:
                    pass

        _cp_if_file(args.bait)
        _cp_if_file(args.screen)
        _cp_if_file(args.chain_map)
        
                    
        print(f"📊  Screening {len(targets)} targets vs bait set ({len(bait_list)} entries)")
        jobs_list = []
        for idx, tgt in enumerate(targets, 1):
            index_str = f"{idx:03d}"
            combined   = ",".join(bait_list + [tgt])

            # Use mapped chain names for job folder naming
            def mapped_token(token: str) -> str:
                """Return user‑friendly name if present in chain_names, else the raw token."""
                return chain_names.get(token.upper(), token)
            bait_label   = _safe_name("_".join([mapped_token(tok) for tok in bait_list]))
            target_label = _safe_name(mapped_token(tgt))
            if args.screen_name:
                job_label = f"{index_str}_{args.screen_name}_{target_label}"
            else:
                job_label = f"{index_str}_{bait_label}_{target_label}"

            sub_out = results_dir / job_label            
            sub_out.mkdir(parents=True, exist_ok=True)

            yaml_src = make_yaml(combined, [], chain_names)
            yaml_dst = sub_out / f"{job_label}.yaml"
            shutil.move(str(yaml_src), yaml_dst)
            print(f"🔄  YAML moved → {yaml_dst}")
            
            # ---- create AF3 JSON ------------------------------------------------
            try:
                af3_path  = yaml_dst.with_suffix("").with_name(yaml_dst.stem + "_for_AF3.json")
                created   = not af3_path.exists()
                if created:
                    yobj     = yaml.safe_load(yaml_dst.read_text())
                    af3_json = yaml_to_af3(yobj, job_label)
                    af3_path.write_text(json.dumps(af3_json, indent=2))
                _link_af3_json(af3_path, root_out)
                status = "written" if created else "ready"
                print(f"   AF3 json {status}: {af3_path.name} – upload to AlphaFold-Server")
            except Exception as exc:
                print(f"⚠  AF3 json FAILED for {yaml_dst.name}: {exc}")
                
            stats = yaml_stats(yaml_dst)
            constraint, mem_req, auto_flags, vram_est = gpu_and_flags(
                stats["residues"], stats["chains"], stats["ligands"], stats["mods"])
            # ---------- skip jobs that won’t fit on any GPU ----------
            if constraint is None:                 # gpu_and_flags() could not place it
                skipped_path = root_out / "skipped_jobs.txt"
                skipped_path.open("a").write(
                    f"{yaml_dst}\tneeds ≈{vram_est:.1f} GB VRAM (skip)\n")
                print(f"⏭  SKIPPED – too large for any GPU: {job_label} "
                    f"(≈{vram_est:.1f} GB VRAM)")
                shutil.rmtree(sub_out, ignore_errors=True)   # clean up empty folder
                continue
                # ---------- NEW: honour --local in screening mode ----------
            if args.local:
                run_local(yaml_dst, sub_out, auto_flags)
                print(f"🏃  Local prediction finished: {job_label}")
                continue            # skip Slurm generation for this target
            flags = auto_flags[:]
            if args.no_kernels:
                flags += ["--no_kernels"]
            if args.boltz:
                flags += shlex.split(args.boltz)

            flags_str = " ".join(flags).strip() or "-"   # “-” marks “no extra flags”

            # ------------------------------------------------------------------
            # save a full per‑job Slurm script (for provenance / reruns)
            # ------------------------------------------------------------------
            part = "g"
            _ = write_slurm(
                job_label,               # job name
                yaml_dst,                # YAML input
                sub_out,                 # output dir (script lives here)
                flags_str if flags_str != "-" else "",  # extra flags
                mem_req,                 # memory request
                part,                    # partition (always 'g')
                constraint               # GPU constraint
            )

            jobs_list.append((yaml_dst, sub_out, flags_str, mem_req, constraint))
        if args.local:
            print("✅  All screening targets were run locally – no Slurm submission.")
            print(f"📂  AF3_JSON hub: {(root_out / 'AF3_JSON').resolve()}")
            return
        # ---------------- write jobs.list --------------------------------
        list_path = root_out / "jobs.list"
        with open(list_path, "w") as fh:
            for y, d, f, m, c in jobs_list:
                fh.write(f"{y}\t{d}\t{f}\t{m}\t{c}\n")

        max_mem = max(int(m[:-1]) for *_, m, _ in jobs_list)  # strip 'M'
        # Use the constraint of the first job (all jobs should use the same heuristic)
        constraint = jobs_list[0][4] if jobs_list else "g4|g3|g2|g1"
        array_slurm = root_out / "array.slurm"
        array_slurm.write_text(ARRAY_TEMPLATE.format(
            root=str(root_out), n=len(jobs_list), mem=f"{max_mem}M").replace("__CONSTRAINT__", constraint))

        # submit array and capture ID
        res = subprocess.run(["sbatch", str(array_slurm)],
                             check=True, capture_output=True, text=True)
        array_id = res.stdout.strip().split()[-1]
        print(f"🗄  array job {array_id} ({len(jobs_list)} tasks)")

        # analysis job
        # ----- analysis job -------------------------------------------------
        analysis_slurm = root_out / "analysis.slurm"
        wrapper_dir = pathlib.Path(__file__).parent
        chain_flag = f"--chain-map {args.chain_map}" if args.chain_map else ""
        ana_script = ANALYSIS_TEMPLATE.format(
            root=str(root_out),
            array_id=array_id,
            wrapper_dir=WRAPPER_DIR,
            chain_flag=chain_flag
        )

        analysis_slurm.write_text(ana_script)
        subprocess.run(["sbatch", str(analysis_slurm)], check=True)

        print("✅  Screening dispatch complete.")
        return

    msa_paths = [pathlib.Path(p).expanduser() for p in args.msa.split(",")] if args.msa else []
    # single‑run: validate input token / path
    _assert_token_ok(args.input)
    yaml_path = make_yaml(args.input, msa_paths, chain_names)
        # ---------------- companion AF-Server JSON (single-run) --------------
    try:
        af3_path  = yaml_path.with_suffix("").with_name(yaml_path.stem + "_for_AF3.json")
        created   = not af3_path.exists()
        if created:
            yobj     = yaml.safe_load(yaml_path.read_text())
            af3_json = yaml_to_af3(yobj, yaml_path.stem)
            af3_path.write_text(json.dumps(af3_json, indent=2))
        _link_af3_json(af3_path, yaml_path.parent.parent)  # root dir
        status = "written" if created else "ready"
        print(f"📝  AF3 json {status}: {af3_path.name} – upload to AlphaFold-Server")
    except Exception as exc:
        print(f"⚠  AF3 json FAILED: {exc}")
    print(f"📂  AF3_JSON hub: {(yaml_path.parent.parent / 'AF3_JSON').resolve()}")
    # ------------------------------------------------------------------ #
    # Optionally inject affinity property                                #
    # ------------------------------------------------------------------ #
    if args.affinity:
        with yaml_path.open() as fh:
            ydoc = yaml.safe_load(fh)
        if "properties" not in ydoc:
            ydoc["properties"] = []
        # check if already present
        if not any("affinity" in p for p in ydoc["properties"]):
            ydoc["properties"].append(
                {"affinity": {"binder": args.affinity}}
            )
            yaml_path.write_text(yaml.safe_dump(ydoc, sort_keys=False))
            print(f"🔧  Affinity block injected with binder={args.affinity}")
        else:
            print("ℹ️  Affinity property already present; skipping injection.")
    stats = yaml_stats(yaml_path)
    constraint, mem_req, auto_flags, vram_est = gpu_and_flags(
        stats["residues"], stats["chains"], stats["ligands"], stats["mods"])
    
    if constraint is None:                      # won’t fit on A100‑40
        # record oversized job and abort gracefully
        skipped_path = pathlib.Path("skipped_jobs.txt")
        skipped_path.open("a").write(
            f"{yaml_path}\tneeds ≈{vram_est:.1f} GB VRAM (skip)\n")
        print("⏭  SKIPPED – too large for any GPU:", yaml_path.name)
        return

    
    flags = auto_flags[:]
    if args.no_kernels:
        flags += ["--no_kernels"]
    if args.boltz:
        flags += shlex.split(args.boltz)

    part = "g"
    outdir = pathlib.Path(args.out or f"boltz_results_{yaml_path.stem}").resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    if args.local:
        run_local(yaml_path, outdir, flags)
        print("✅  Local run finished.")
        return

    slurm_script = write_slurm("boltz", yaml_path, outdir, " ".join(flags), mem_req, part, constraint)
    cmd = ["sbatch", str(slurm_script), str(yaml_path), str(outdir), " ".join(flags)]
    print("# Residues =", stats["residues"], "| flags:", " ".join(flags) or "none")
    print("Submitting:", " ".join(cmd))

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Submitting job via sbatch …")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode:
        sys.stderr.write(result.stderr)
        sys.exit(result.returncode)
    print("✅  Slurm job", result.stdout.strip())
    print("🕒  You can monitor with: squeue -j", result.stdout.split()[-1])


def run_local(yaml_path: pathlib.Path, outdir: pathlib.Path, flags: list[str]) -> None:
    """Run boltz predict inside Apptainer directly on this node."""
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "apptainer", "exec", "--nv",
        "--bind", f"{yaml_path}:{yaml_path}",
        "/groups/plaschka/shared/software/boltz-2/boltz2:develop",
        "boltz", "predict", str(yaml_path),
        "--out_dir", str(outdir),
        "--accelerator", "gpu",
        "--use_msa_server", *flags
    ]
    print("Running:", " ".join(shlex.quote(c) for c in cmd))
    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()
