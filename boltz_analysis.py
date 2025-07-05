#!/usr/bin/env python3
"""
boltz_screen_plotly.py
----------------------

• Creates a `summary_metrics.csv`
• Builds an interactive HTML dashboard for scalar metrics
• Generates one HTML heat-map (PAE / PDE / PLDDT etc.) per job
• By default the script makes iptm swarm plots, PAE heat-maps and a
  scatter dashboard.  Use --all-plots to generate every available plot.

Static images: SVG always, PNG if `cairosvg` installed.

Run inside your screen root:

    python boltz_screen_plotly.py .

Outputs
-------
summary_metrics.csv
plots/
   scalar_dashboard.html
   001_test_combo_Q09161_pae.html
   002_test_combo_Q13838_pde.html
   ...

Note: if you run `./boltz_analysis.py` ensure the file has executable permission
(`chmod +x`) or invoke with `python boltz_analysis.py .`.
"""

# Imports
import argparse, json, pathlib, re, textwrap
import subprocess, yaml, shlex, re as _re
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
from plotly.subplots import make_subplots

import functools, time

import random
import requests

# ------------------------------------------------------------------
# ChimeraX helper – write a ready‑to‑open .cxc per prediction folder
# ------------------------------------------------------------------
CXC_TEMPLATE = textwrap.dedent("""\
open predictions/*/*.cif
open predictions/*/pae*.npz structure last-opened
rainbow last-opened chains palette bupu
interfaces last-opened
hide
show ligand
show nucleic
show c
style (protein|nucleic|solvent) & @@draw_mode=0 stick

cartoon style modeh def arrows t arrowshelix f arrowscale 1.5 wid 2 thick 1 sides 12 div 20
cartoon style ~(nucleic|strand) x round
cartoon style (nucleic|strand) x rect
cartoon style nucleic x round width 4 thick 4
nucleotides stubs
lighting shadows false
lighting simple
camera ortho

lighting soft intensity 0.2 fillIntensity 0.5 ambientIntensity 1
lighting soft

#bg and model colors
set bgColor white

nucleotides  fill
style nucleic  stick

graphics silhouettes true

cartoon style  arrowScale 1.5
""")

def _write_cxc(job_dir: pathlib.Path):
    """Create {job}_analysis.cxc inside *job_dir* if not present."""
    cxc_path = job_dir / f"CHIMERAX_{job_dir.name}_analysis.cxc"
    try:
        if not cxc_path.exists():
            cxc_path.write_text(CXC_TEMPLATE)
    except OSError:
        pass

# ---------------------------------------------------------------------------
# UniProt helper -------------------------------------------------------------
# ---------------------------------------------------------------------------
_UNIPROT_LABEL_CACHE: dict[str, str] = {}


def uniprot_label(uid: str) -> str:
    """
    Return a concise “GENE – description” label for a UniProt accession.
    Falls back to the bare accession on errors.  Results are cached.
    """
    if uid in _UNIPROT_LABEL_CACHE:
        return _UNIPROT_LABEL_CACHE[uid]

    label = uid  # fallback
    try:
        r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uid}.txt",
                         timeout=8)
        if r.status_code == 200:
            txt = r.text
            # description line
            desc = ""
            for line in txt.splitlines():
                if line.startswith("DE   RecName:"):
                    desc = line.partition(":")[2].strip().rstrip(".")
                    break
            gene = uid
            m = re.search(r"GN   Name=([A-Za-z0-9_-]+)", txt)
            if m:
                gene = m.group(1)
            label = f"{gene} – {desc[:60]}" if desc else gene
    except Exception:
        pass

    _UNIPROT_LABEL_CACHE[uid] = label
    return label

# --- Slurm string helpers ----------------------------------------------
import math

# Global cache for slurm_stats
_SLURM_STAT_CACHE = {}

def _mem_to_mb(s: str) -> float | None:
    """Convert strings like '81056M' or '8.67G' to MiB (float)."""
    if not s:
        return None
    m = _re.match(r"([\d.]+)\s*([KMG]?)", s)
    if not m:
        return None
    val, unit = float(m.group(1)), m.group(2).upper()
    factor = {'K': 1/1024, 'M': 1, 'G': 1024}.get(unit, 1)
    return val * factor


def _hms_to_min(t: str) -> float | None:
    """Convert HH:MM:SS wall‑time strings to minutes (float)."""
    if not t or t.count(":") != 2:
        return None
    h, m, s = map(int, t.split(":"))
    return h * 60 + m + s / 60

# Simple random jitter helper
def random_jitter(n: int, width: float = 0.35) -> np.ndarray:
    """Return *n* random offsets in [-width, +width]."""
    rng = random.Random(0)          # deterministic
    return np.array([rng.uniform(-width, width) for _ in range(n)])

from collections import defaultdict
def beeswarm_offsets(values: np.ndarray, radius: float = 0.12,
                     max_points: int = 800) -> np.ndarray:
    """
    Return symmetric x‑offsets (beeswarm).  O(N log N) with grid bins.

    For very large N (> max_points) fall back to random jitter to avoid
    pathological runtimes.
    """
    n = len(values)
    if n > max_points:
        return random_jitter(n, width=radius*3)

    order   = np.argsort(values)
    offsets = np.zeros(n)

    # binning on y to reduce neighbour checks
    bin_h = radius             # vertical bin height
    bins  = defaultdict(list)  # bin_index -> indices placed in that bin

    for idx in order:
        y = values[idx]
        bin_idx = int(y // bin_h)
        neigh   = bins[bin_idx-1] + bins[bin_idx] + bins[bin_idx+1]

        x = 0.0
        while any((abs(y - values[j]) < radius) and (abs(x - offsets[j]) < radius)
                  for j in neigh):
            x = -x if x > 0 else abs(x) + radius

        offsets[idx] = x
        bins[bin_idx].append(idx)

    return offsets

# use clean white layout globally
px.defaults.template = "simple_white"

VERBOSE = False
_warned_kaleido = False    # ensure only one Kaleido warning is printed
def vprint(*a, **kw):
    if VERBOSE:
        print(*a, **kw)

#
# Ensure static image export (requires kaleido)
try:
    import kaleido  # noqa: F401
except ModuleNotFoundError:
    print("⚠  'kaleido' not installed – PNG export will be skipped.")

# -- Optional PNG export ----------------------------------------------------
# Kaleido (for SVG) works without extra libs, but PNG via cairosvg needs
# a system lib 'cairo'.  If that is missing we silently disable PNG once.

HAVE_CAIROSVG = False
try:
    import cairosvg  # noqa: F401
    HAVE_CAIROSVG = True
except Exception as _e:   # ImportError or OSError (missing libcairo)
    # one concise warning, printed exactly once
    print(f"ℹ️  PNG export disabled (cairosvg unavailable: {_e.__class__.__name__})")

def _save(fig, html_path: pathlib.Path, png_scale: int = 2):
    """
    Write interactive HTML.
    Try Kaleido → SVG → PNG (via cairosvg); if Kaleido missing, only HTML.
    """
    fig.write_html(html_path, include_plotlyjs="cdn")

    # Attempt Kaleido for static images
    try:
        svg_data = fig.to_image(format="svg", scale=png_scale)  # Kaleido backend
    except Exception:
        # print the Kaleido‑missing notice only once
        global _warned_kaleido
        if not _warned_kaleido:
            print("ℹ️  Static SVG/PNG export skipped – Kaleido unavailable.")
            _warned_kaleido = True
        return

    svg_path = html_path.with_suffix(".svg")
    svg_path.write_bytes(svg_data)

    if HAVE_CAIROSVG:
        cairosvg.svg2png(bytestring=svg_data,
                        write_to=str(html_path.with_suffix(".png")))

# --------------------------------------------------------------------------- #
# discovery                                                                   #
# --------------------------------------------------------------------------- #
FAILED_MARKER = "failed"   # global constant

def walk_predictions(root: pathlib.Path):
    """
    Yield tuples (job_name, status, conf_json_or_None, matrices_dict, yaml_path_or_None)
      • status == "ok"      → run finished, confidence JSON present
      • status == "failed"  → no confidence JSON (OOM, runtime error, etc.)
    """
    def is_numeric_square(arr: np.ndarray) -> bool:
        return (np.issubdtype(arr.dtype, np.number)
                and arr.ndim == 2
                and arr.shape[0] == arr.shape[1])

    base_dir = root / "results" if (root / "results").is_dir() else root
    for job_dir in sorted(base_dir.glob("*_*")):
        if not job_dir.is_dir():
            # skip stray files like slurm_*.err/out etc.
            continue
        conf_jsons = list(job_dir.rglob("confidence_*json"))

        # collect yaml path (if any)
        yaml_candidates = list(job_dir.glob("*.yaml"))
        yaml_path = yaml_candidates[0] if yaml_candidates else None

        # ----- failed job --------------------------------------------------
        if not conf_jsons:
            vprint(f"→ found *FAILED* folder (no confidence JSON): {job_dir.name}")
            # Preserve any YAML we may have so the caller can still extract
            # chain counts / residue numbers etc. for correlation analysis
            yield job_dir.name, FAILED_MARKER, None, {}, yaml_path
            continue

        # ----- successful job ---------------------------------------------
        conf_json = conf_jsons[0]
        matrices: dict[str, list[np.ndarray]] = {}
        for npz in job_dir.rglob("*.npz"):
            try:
                with np.load(npz, allow_pickle=False) as raw:
                    first = raw[raw.files[0]]
                if first.dtype.fields is not None or first.ndim != 2 or first.shape[0] != first.shape[1]:
                    continue
            except Exception:
                continue
            key = npz.stem.split("_")[0]
            matrices.setdefault(key, []).append(first)
            vprint(f"   collected matrix {key} {first.shape} for {job_dir.name}")

        vprint(f"→ found prediction folder: {job_dir.name}")
        yield job_dir.name, "ok", conf_json, matrices, yaml_path

def load_conf(path):
    with open(path) as fh:
        return json.load(fh)

# --------------------------------------------------------------------------- #
# helper to shorten job names                                                  #
# --------------------------------------------------------------------------- #
def short_label(job: str) -> str:
    """
    Return 'IDX_Target' from a job string like
    '014_uap56_atp_mg_vs_top_50_Q53F19'  ->  '014_Q53F19'
    """
    parts = job.split("_")
    if len(parts) < 2:
        return job
    idx  = parts[0]
    tgt  = parts[-1]
    return f"{idx}_{tgt}"

# --------------------------------------------------------------------------- #
# Slurm stats helper                                                          #
# --------------------------------------------------------------------------- #
def slurm_stats(job_dir: pathlib.Path) -> dict[str, str]:
    """
    Fetch Slurm run‑time statistics for *job_dir*:

     1. Look for ``slurm_*.out`` *inside the job directory*.
     2. If none found (common when Slurm stdout is written at the
        screen‑root), scan the **parent directory** for a ``slurm_*.out``
        file mentioning the job‑folder name in its contents.
     3. Extract the JOBID (or JOBID_TASK) from the file name and query
        either::

            jobinfo <JOBID>
        or, if that fails,

            sacct  -j <JOBID>  --parsable2 \
                   --format JobID,Elapsed,MaxRSS,Reserved,State

    The function returns a *dict* which may include:

        • ``"Mem reserved"``         (e.g. ``81056M``)
        • ``"Full Max Mem usage"``   (e.g. ``8.67G``)
        • ``"Used walltime"``        (``HH:MM:SS``)
        • ``"State"``                (``COMPLETED``, ``OOM``, …)

    If no statistics can be obtained an **empty dict** is returned.
    """
    # ------------------------------------------------------------- #
    # fast path: return cached result if we've already seen job id
    # ------------------------------------------------------------- #
    def _from_cache(jid: str):
        cached = _SLURM_STAT_CACHE.get(jid)
        if cached is not None:
            return cached.copy()
        return None

    # ------------------------------------------------------------------
    # Extra heuristic (failed‑jobs case):
    # If no Slurm file lives *inside* the job directory, try to match
    #   slurm_<JOBID>_<TASK>.out
    # in the parent folder by comparing the *array index* (first token
    # of the job‑folder name, e.g. "041"  in  "041_jobname_…").
    # ------------------------------------------------------------------
    folder_idx = None
    m_idx = _re.match(r"^(\d+)[_A-Za-z]", job_dir.name)
    if m_idx:
        folder_idx = m_idx.group(1).lstrip("0")  # strip leading zeros

    # 1) locate slurm output file belonging to this job directory
    out_files = sorted(job_dir.glob("slurm_*[0-9].out"))

    # If nothing inside the folder, look at parent directory – this is
    # the layout produced by boltz‑2_wrapper at the moment.
    if not out_files:
        for p in job_dir.parent.glob("slurm_*[0-9].out"):
            try:
                with open(p, "r") as fh:
                    head = fh.read(2000)        # scan a small chunk
                if job_dir.name in head:
                    out_files = [p]
                    break
            except OSError:
                continue
        # Extra heuristic: match by array index if job name not found
        if not out_files and folder_idx:
            # Try to match by the numeric array index if job name not present
            for p in job_dir.parent.glob(f"slurm_*_{folder_idx}.out"):
                out_files = [p]
                break
    # --- additional search one level higher and in root/slurm -------
    if not out_files:
        parent2 = job_dir.parent.parent           # usually the screen root
        # direct slurm_* files in the root folder
        for p in parent2.glob("slurm_*[0-9].out"):
            try:
                with open(p, "r") as fh:
                    head = fh.read(2000)
                if job_dir.name in head:
                    out_files = [p]
                    break
            except OSError:
                continue

        # dedicated root/slurm directory
        if not out_files:
            slurm_root = parent2 / "slurm"
            if slurm_root.is_dir():
                for p in slurm_root.glob("slurm_*[0-9].out"):
                    try:
                        with open(p, "r") as fh:
                            head = fh.read(2000)
                        if job_dir.name in head:
                            out_files = [p]
                            break
                    except OSError:
                        continue
                if not out_files and folder_idx:
                    for p in slurm_root.glob(f"slurm_*_{folder_idx}.out"):
                        out_files = [p]
                        break
    # If still nothing was found bail out early – avoids IndexError
    if not out_files:
        return {}
    # 2) extract JOBID (array task IDs look like 12345_7)
    m = _re.search(r"slurm_([0-9]+(?:_[0-9]+)?)", out_files[0].name)
    if not m:
        return {}
    jobid_full = m.group(1)
    jobid_base = jobid_full.split("_")[0]

    cached = _from_cache(jobid_base)
    if cached is not None:
        return cached

    def _run(cmd: list[str]) -> list[str]:
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, check=True)
            return res.stdout.splitlines()
        except (FileNotFoundError, subprocess.CalledProcessError):
            return []

    lines = _run(["jobinfo", jobid_full]) or _run(["jobinfo", jobid_base])

    # fallback ‑‑ sacct
    if not lines:
        lines = _run(["sacct", "-j", jobid_base, "--parsable2",
                      "--format", "JobID,Elapsed,MaxRSS,Reserved,State"])

    if not lines:
        return {}

    info: dict[str, str] = {}

    for ln in lines:
        # jobinfo style "Key  :  value"
        if " : " in ln:
            k, v = [x.strip() for x in ln.split(":", 1)]
            k = {"Max Mem (Node/step)": "Full Max Mem usage",
                 "Reserved walltime":   "Reserved walltime"}.get(k, k)
            info[k] = v
            continue

        # sacct (pipe‑separated)
        if "|" in ln and ln.startswith(jobid_base):
            _, elapsed, maxrss, reserved, state = (ln.split("|") + [""]*5)[:5]
            if elapsed:
                info["Used walltime"] = elapsed.replace("-", ":")
            if maxrss:
                info["Full Max Mem usage"] = maxrss
            if reserved:
                info["Mem reserved"] = reserved
            info["State"] = state

    # add to cache (use base id so task ids hit same entry)
    _SLURM_STAT_CACHE[jobid_base] = info

    return info

def yaml_summary(yaml_path: pathlib.Path) -> dict[str, int]:
    # keep counters for every possible sequence type plus overall stats
    out = {
        "seqs": 0, "residues": 0, "mods": 0,          # overall
        "protein": 0, "rna": 0, "dna": 0,             # chain‑type counts
        "ligands": 0                                  # number of ligand/ion entries
    }
    if not yaml_path or not yaml_path.is_file():
        return out
    y = yaml.safe_load(yaml_path.read_text())
    for item in y.get("sequences", []):
        k, d = next(iter(item.items()))
        if k == "ligand":
            out["ligands"] += 1
            continue
        if k in ("protein", "rna", "dna"):
            out["seqs"] += 1
            # ignore unexpected keys gracefully
            if k not in out:
                continue
            out[k] += 1
            out["residues"] += len(d["sequence"])
            out["mods"] += len(d.get("modifications", []))
    return out

# --------------------------------------------------------------------------- #
# plotting helpers                                                            #
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------------------- #
# per‑metric vertical dot plot                                                #
# --------------------------------------------------------------------------- #

def plot_metric_dot(df: pd.DataFrame, metric: str, out_html: pathlib.Path,
                    style: str = "strip", annotate: bool = True):
    if metric not in df.columns:
        return

    d = df.sort_values(metric, ascending=False).reset_index(drop=True).copy()
    # ensure numeric values
    d[metric] = pd.to_numeric(d[metric], errors="coerce")
    d = d.dropna(subset=[metric])
    if d.empty:
        return
    d["dummy"] = ""   # single x category

    import time
    t0 = time.time()
    if style == "swarm" and len(d) <= 800:
        # Only run costly beeswarm for reasonably sized sets
        d["offset"] = beeswarm_offsets(d[metric].to_numpy(), radius=0.12)
    else:  # strip – or auto-downgraded swarm when >800 points
        d["offset"] = random_jitter(len(d), width=0.35)
    if VERBOSE and style == "swarm":
        vprint(f"      swarm offsets computed in {time.time()-t0:.3f}s for {len(d)} dots [{metric}]")

    fig = go.Figure(go.Scatter(
        x=d["offset"],
        y=d[metric],
        mode="markers",
        marker=dict(size=8, line=dict(width=0), color="mediumpurple"),
        text=d["hover_label"],
        hovertemplate="%{text}<br>%{y:.3f}<extra></extra>",
    ))
    fig.update_xaxes(visible=False, showticklabels=False)
    fig.update_layout(
        title=f"{metric} – {style} plot",
        showlegend=False,
        template="simple_white",
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=40, r=40, t=80, b=40),
    )

    # optional annotations: top-5
    if annotate:
        top_n = min(5, len(d))
        for idx in range(top_n):
            fig.add_annotation(
                x=d.loc[idx, "offset"],
                y=d.loc[idx, metric],
                text=d.loc[idx, "hover_label"],
                showarrow=True, arrowhead=1, ax=40, ay=0, font_size=10
            )

    _save(fig, out_html)
    vprint(f"   dot plot {metric} done → {out_html.name}")
# --------------------------------------------------------------------------- #
# combined multi‑metric dot plot                                              #
# --------------------------------------------------------------------------- #
def plot_combined_dots(df: pd.DataFrame, out_html: pathlib.Path,
                       style: str = "strip",
                       top_n: int | None = None,
                       annotate: bool = True):
    metrics = [c for c in df.columns if c != "job"]
    if not metrics:
        return

    # Determine a safe vertical spacing that respects Plotly limits:
    nrows = len(metrics)
    if nrows > 1:
        # Plotly requires vertical_spacing ≤ 1/(rows‑1)
        max_vs = (1 / (nrows - 1)) - 1e-4
        vspace = min(0.06, max_vs)
    else:
        vspace = 0.0

    fig = make_subplots(
        rows=nrows,
        cols=1,
        shared_xaxes=False,
        vertical_spacing=vspace,
        subplot_titles=metrics,
    )

    palette = px.colors.qualitative.Dark24

    for r, metric in enumerate(metrics, 1):
        d = df.sort_values(metric, ascending=False).reset_index(drop=True).copy()
        d[metric] = pd.to_numeric(d[metric], errors="coerce")
        d = d.dropna(subset=[metric])
        if d.empty:
            continue

        if style == "swarm":
            d["jitter"] = beeswarm_offsets(d[metric].to_numpy(), radius=0.12)
        else:
            d["jitter"] = random_jitter(len(d), width=0.6)
        colours = ["#d9d9d9"] * len(d)
        sel = range(len(d)) if top_n is None else range(min(top_n, len(d)))
        for rank, idx in enumerate(sel):
            colours[idx] = palette[rank % len(palette)]

        fig.add_trace(
            go.Scatter(
                x=d["jitter"],
                y=d[metric],
                mode="markers",
                marker=dict(color=colours, size=7),
                text=d["hover_label"],
                hovertemplate="%{text}<br>%{y:.3f}<extra></extra>",
                showlegend=False,
            ),
            row=r, col=1,
        )
        # annotate top‑N
        for idx in sel:
            if annotate:
                fig.add_annotation(
                    x=d.loc[idx, "jitter"],
                    y=d.loc[idx, metric],
                    text=d.loc[idx, "hover_label"],
                    showarrow=True, arrowhead=1, ax=40, ay=0, font_size=9,
                    row=r, col=1,
                )

        fig.update_xaxes(visible=False, row=r, col=1)
        fig.update_yaxes(title_text=metric, row=r, col=1)

    fig.update_layout(
        height=300*len(metrics)+100,
        title="Boltz‑2 screening metrics – combined dot plots",
        margin=dict(l=60, r=40, t=80, b=40),
        template="simple_white",
        plot_bgcolor="white",
        paper_bgcolor="white",
    )
    _save(fig, out_html)
    vprint("   combined dot plot written:", out_html.name)
    print("✓ combined dot plot:", out_html.name)

def plot_matrix(mat, title, out_html, typ=""):
    """
    Heat-map plot helper.

    • Uses **RdBu_r** for PAE (blue = high confidence, red = low)
    • Viridis for everything else.
    """
    mat = np.asarray(mat, dtype=float)
    cmap = "RdBu_r" if typ.lower() == "pae" else "Viridis_r"   # invert Viridis
    fig = px.imshow(mat,
                    color_continuous_scale=cmap,
                    aspect="auto",
                    origin="lower",
                    title=title)
    fig.update_layout(height=500, width=500,
                      margin=dict(l=40,r=40,t=80,b=40))
    _save(fig, out_html)

# --------------------------------------------------------------------------- #
# interactive scatter dashboard                                               #
# --------------------------------------------------------------------------- #
def plot_interactive_scatter(df: pd.DataFrame,
                             out_html: pathlib.Path,
                            filter_ok: bool = True,
                            init_x: str | None = None,
                            init_y: str | None = None,
                            init_color: str | None = None,
                            init_size: str | None = None):
    """
    Build a Plotly dashboard with 3 dropdowns:
      • X‑axis metric
      • Y‑axis metric
      • colour‑by metric
    Only numeric cols are offered.
    """
    # ------------------------------------------------------------------
    # Ensure we only plot successful predictions (status == "ok")
    # Even though the caller already passes ok_df, an inline check makes
    # the function robust when reused elsewhere.
    # ------------------------------------------------------------------
    if filter_ok and "status" in df.columns:
        df = df[df["status"] == "ok"].copy()
    numeric_cols = [c for c in df.columns
                    if c not in ("job", "status")
                    and pd.api.types.is_numeric_dtype(df[c])
                    and df[c].notna().any()]   # keep only columns that have real data

    if len(numeric_cols) < 2:
        return   # nothing to plot

    def _pick(name, fallback_idx):
        return name if (name in numeric_cols) else (numeric_cols[fallback_idx]
                                                    if len(numeric_cols) > fallback_idx else None)
    x0 = _pick(init_x, 0)
    y0 = _pick(init_y, 1)
    c0 = _pick(init_color, 2)
    s0 = _pick(init_size, 3)

    fig = go.Figure()

    if s0:
        base_size = df[s0]
        size_ref = 2.0 * base_size.max() / 40.0**2   # Plotly’s area scaling
    else:
        base_size = 10
        size_ref = 1

    marker_dict = dict(
        size=base_size,
        sizemode="area",
        sizeref=size_ref,
        sizemin=4,
        colorscale="Viridis_r",
        color=df[c0] if c0 else "mediumpurple",
        showscale=bool(c0),
        colorbar=dict(title=c0) if c0 else None,
    )

    fig.add_trace(go.Scatter(
        x=df[x0], y=df[y0],
        mode="markers",
        marker=marker_dict,
        text=df["hover_label"],
        hovertemplate="%{text}<br>X=%{x:.3f}<br>Y=%{y:.3f}<extra></extra>"
    ))

    # --- build dropdowns ---------------------------------------------------
    def make_buttons(kind: str):
        """
        kind:
            'x' – X-axis
            'y' – Y-axis
            'c' – colour
            's' – size
        Each button updates the corresponding trace attributes.
        """
        buttons = []
        for col in numeric_cols:
            if kind == "x":      # X-axis
                args = [{"x": [df[col]]},
                        {"xaxis.title.text": col}]
            elif kind == "y":    # Y-axis
                args = [{"y": [df[col]]},
                        {"yaxis.title.text": col}]
            elif kind == "c":    # colour
                args = [{
                    "marker.color": [df[col]],
                    "marker.colorscale": ["Viridis_r"],
                    "marker.showscale": [True],
                    "marker.colorbar.title": [col]
                }]
            else:                # size
                args = [{"marker.size": [df[col]],
                        "marker.sizeref": [2.*df[col].max()/40.0**2]}]
            buttons.append(dict(label=col, method="update", args=args))
        return buttons

    fig.update_layout(
        updatemenus=[
            dict(buttons=make_buttons("x"), direction="down", x=0.00, y=1.2,
                pad={"r": 8, "t": 10}, showactive=True, bgcolor="white",
                font=dict(size=11), xanchor="left"),
            dict(buttons=make_buttons("y"), direction="down", x=0.12, y=1.2,
                pad={"r": 8, "t": 10}, showactive=True, bgcolor="white",
                font=dict(size=11), xanchor="left"),
            dict(buttons=make_buttons("c"), direction="down", x=0.24, y=1.2,
                pad={"r": 8, "t": 10}, showactive=True, bgcolor="white",
                font=dict(size=11), xanchor="left"),
            dict(buttons=make_buttons("s"), direction="down", x=0.36, y=1.2,
                pad={"r": 8, "t": 10}, showactive=True, bgcolor="white",
                font=dict(size=11), xanchor="left"),
        ],
        annotations=[
            dict(text="X:", x=-0.04, y=1.2, xref="paper", yref="paper",
                 showarrow=False, font=dict(size=12)),
            dict(text="Y:", x=0.08, y=1.2, xref="paper", yref="paper",
                 showarrow=False, font=dict(size=12)),
            dict(text="Size:", x=0.32, y=1.2, xref="paper", yref="paper",
                showarrow=False, font=dict(size=12)),
        ],
        title="Interactive scatter plot – choose metrics",
        height=600,
        margin=dict(l=60, r=40, t=100, b=40),
        template="simple_white",
        plot_bgcolor="white",
        paper_bgcolor="white",
    )
    _save(fig, out_html)
    vprint("   scatter dashboard written:", out_html.name)
    print("✓ scatter dashboard:", out_html)

# --------------------------------------------------------------------------- #
# main                                                                        #
# --------------------------------------------------------------------------- #
def main():
    pa = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__doc__))
    pa.add_argument("screen_root", help="Directory produced by boltz-2_wrapper.py")
    pa.add_argument("--labels", dest="labels", action="store_true",
                    help="Annotate top‑N dots with static text "
                         "(index_target). Use --no-labels to disable.")
    pa.set_defaults(labels=True)
    pa.add_argument("--no-labels", dest="labels", action="store_false")
    pa.add_argument("-v", "--verbose", action="store_true",
                    help="Print progress messages")
    pa.add_argument(
        "--all-plots", action="store_true",
        help=("Generate *all* dot-plots and matrix heat-maps. "
              "If omitted, the script makes only iptm-related swarm plots, "
              "PAE matrices, and the interactive scatter dashboard."))
    pa.add_argument("--style", choices=["strip", "swarm"], default="swarm",
                    help="Dot arrangement for per‑metric plots: 'strip' (random jitter) or 'swarm' (symmetric beeswarm, default).")
    pa.add_argument("--chain-map", metavar="FILE",
                    help="Optional text file that maps numerical chain indices to human‑readable "
                         "names.  Format: one mapping per line, either '0=RNA' or '0 RNA'.  "
                         "These names replace occurrences of 'chain0', 'chain1', … in metric "
                         "column titles and in the interactive scatter dashboard. (--color-map is an alias.)")
    pa.add_argument("--color-map", "-c", dest="chain_map", metavar="FILE",
                    help="Alias for --chain-map. Optional text file mapping numerical chain indices to human‑readable names.")
    args = pa.parse_args()
    all_plots = args.all_plots
    
    global VERBOSE
    VERBOSE = args.verbose
    root = pathlib.Path(args.screen_root).expanduser().resolve()
    vprint("Gathering data and making plots, please be patient...")

    # ------------------------------------------------------------------
    # optional custom chain names
    # ------------------------------------------------------------------
    chain_names: dict[str, str] = {}
    if args.chain_map:
        cmap_path = pathlib.Path(args.chain_map).expanduser()
        if cmap_path.is_file():
            with open(cmap_path) as fh:
                for ln in fh:
                    ln = ln.strip()
                    if not ln or ln.startswith("#"):
                        continue
                    if "=" in ln:
                        idx, name = ln.split("=", 1)
                    else:
                        idx, name = ln.split(None, 1)
                    idx = idx.strip()
                    chain_names[f"chain{idx}"] = name.strip()
        else:
            print(f"⚠  chain‑map file not found: {cmap_path}")

    # Helper to map chain tokens to friendly names if available
    def mapped_token(token: str) -> str:
        """Return friendly chain name if provided, else the original token."""
        key = token.upper()
        return chain_names.get(key, token)


    plots_dir = root / "plots"
    plots_dir.mkdir(exist_ok=True)
    dot_dir = plots_dir / "dot_plots"
    mat_dir = plots_dir / "matrix_plots"
    dot_dir.mkdir(exist_ok=True)
    mat_dir.mkdir(exist_ok=True)

    df_rows: list[dict] = []

    # ------------------------------------------------------------------- #
    # iterate through all prediction folders                              #
    # ------------------------------------------------------------------- #
    for job, status, conf_json, mats, yaml_file in walk_predictions(root):
        # absolute path of this prediction folder
        job_dir_path = (root / "results" / job) if (root / "results").is_dir() else (root / job)
        _write_cxc(job_dir_path)
        if status == FAILED_MARKER:
            parts = job.split("_")
            idx  = parts[0] if parts else ""
            acc  = parts[-1] if parts else ""
            row = {
                "job": job,
                "status": FAILED_MARKER,
                "hover_label": f"{idx} – {uniprot_label(acc)}"
            }
            # --- Slurm statistics (convert to numeric like for successful jobs) ---
            sstats = slurm_stats((root / job))
            row["mem_reserved"]  = _mem_to_mb(sstats.get("Mem reserved", ""))
            row["max_mem_used"]  = _mem_to_mb(sstats.get("Full Max Mem usage", ""))
            row["walltime_min"]  = _hms_to_min(sstats.get("Used walltime", ""))

            # YAML may still exist even for failed runs (e.g. early OOM)
            if yaml_file and yaml_file.is_file():
                row.update(yaml_summary(yaml_file))

            df_rows.append(row)
            continue

        conf = load_conf(conf_json)
        acc = job.split("_")[-1]
        idx = job.split("_")[0]
        row = {"job": job,
            "status": "ok",
            "hover_label": f"{idx} – {uniprot_label(acc)}"}
        row.update({k: v for k, v in conf.items() if isinstance(v, (int, float))})
        # per-chain metrics
        for cid, val in conf.get("chains_ptm", {}).items():
            row[f"chain{cid}_ptm"] = val
        for c1, sub in conf.get("pair_chains_iptm", {}).items():
            for c2, val in sub.items():
                row[f"iptm_chain_{c1}_vs_chain_{c2}"] = val
        # --- Slurm stats -------------------------------------------------
        sstats = slurm_stats((root / job))  
        row["mem_reserved"]   = _mem_to_mb(sstats.get("Mem reserved", ""))   # MiB
        row["max_mem_used"]  = _mem_to_mb(sstats.get("Full Max Mem usage", ""))
        row["walltime_min"]  = _hms_to_min(sstats.get("Used walltime", ""))
        # yaml details
        if yaml_file:
            ysum = yaml_summary(yaml_file)
            row.update(ysum)
        df_rows.append(row)

        # --- matrix plots (PAE, PDE, etc.) ------------------------------- #
        groups: dict[tuple[str, tuple[int, int]], list[np.ndarray]] = {}
        for typ, arrs in mats.items():
            for arr in arrs:
                if arr.dtype.fields is None and arr.ndim == 2 and arr.shape[0] == arr.shape[1]:
                    groups.setdefault((typ, arr.shape), []).append(arr)
        for (typ, shape), arrs in groups.items():
            if not all_plots and typ.lower() != "pae":
                continue
            matrix = np.mean(arrs, axis=0) if len(arrs) > 1 else arrs[0]
            out_name = f"{job}_{typ}.html"
            plot_matrix(matrix, f"{job} – {typ.upper()}",
                mat_dir / out_name, typ)

    if not df_rows:
        print("No predictions found in", root)
        return

    df = pd.DataFrame(df_rows).sort_values("job")
    # ------------------------------------------------------------------
    # Ensure Slurm‑derived numeric columns are real numbers.
    # Empty strings ("") would make the whole column dtype=object and
    # then these columns get ignored for the interactive scatter menu.
    # Convert to float, coercing blanks to NaN.
    # ------------------------------------------------------------------
    for col in ("mem_reserved", "max_mem_used", "walltime_min"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df.to_csv(root / "summary_metrics.csv", index=False)
    print("✏  summary_metrics.csv written.")
    analytics_dir = root / "analytics"
    analytics_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(analytics_dir / "slurm_metrics.csv", index=False)
    print("✏  slurm_metrics.csv written.")

    # --- Job‑run statistics ------------------------------------------------
    total_jobs  = len(df)
    failed_jobs = (df["status"] == FAILED_MARKER).sum()
    ok_jobs     = total_jobs - failed_jobs
    print(f"🧮  Job summary: {total_jobs} total — {ok_jobs} ok, {failed_jobs} failed.")

    # For failed predictions we capture Slurm metrics separately
    if failed_jobs:
        failed_df = df[df["status"] == FAILED_MARKER].copy()

        # Ensure relevant numeric columns exist (add blanks if missing)
        needed_cols = ["job", "mem_reserved", "max_mem_used",
                       "walltime_min", "seqs", "residues", "mods"]
        for col in needed_cols:
            if col not in failed_df.columns:
                failed_df[col] = np.nan

        # Save a dedicated CSV with the extra Slurm information
        failed_df.to_csv(analytics_dir / "failed_slurm_metrics.csv", index=False)
        print(f"✏  failed_slurm_metrics.csv written ({len(failed_df)} rows).")

        # Quick aggregate stats – handy when diagnosing OOM / runtime problems
        mem_vals = failed_df["mem_reserved"].dropna()
        if not mem_vals.empty:
            print(f"   Failed mem_reserved: mean {mem_vals.mean():.0f} MiB — "
                  f"max {mem_vals.max():.0f} MiB")
        wt_vals = failed_df["walltime_min"].dropna()
        if not wt_vals.empty:
            print(f"   Failed walltime   : mean {wt_vals.mean():.1f} min — "
                  f"max {wt_vals.max():.1f} min")

        # ------------------------------------------------------------------
        # Scatter dashboard exclusively for failed jobs
        # ------------------------------------------------------------------
        if not failed_df.empty:
            # Remove columns that are entirely NaN to avoid empty dropdown entries
            failed_df = failed_df.dropna(axis=1, how='all')
            failed_plot = plots_dir / "failed_scatter_dashboard.html"
            plot_interactive_scatter(failed_df, failed_plot,
                                     filter_ok=False)
            print("✓ failed‑jobs scatter dashboard:", failed_plot)

    ok_df = df[df["status"] == "ok"].copy()
    # ------------------------------------------------------------------
    # apply custom chain labels to dataframe columns
    # ------------------------------------------------------------------
    if chain_names:
        new_cols = {}
        for col in ok_df.columns:
            # rename chainX_ptm
            m = _re.match(r"^chain(\d+)_ptm$", col)
            if m:
                idx = m.group(1)
                if f"chain{idx}" in chain_names:
                    new_cols[col] = f"{chain_names[f'chain{idx}']}_ptm"
                    continue
            # rename iptm_chain_i_vs_chain_j
            m = _re.match(r"^iptm_chain_(\d+)_vs_chain_(\d+)$", col)
            if m:
                i, j = m.groups()
                repl_i = chain_names.get(f"chain{i}", f"chain{i}")
                repl_j = chain_names.get(f"chain{j}", f"chain{j}")
                new_cols[col] = f"iptm_{repl_i}_vs_{repl_j}"
        if new_cols:
            ok_df.rename(columns=new_cols, inplace=True)
    if ok_df.empty:
        print("No successful predictions – skipping plots.")
        return

    # plot_scalar_dashboard(ok_df, plots_dir / "scalar_dashboard.html")
    if all_plots:
        metric_list = [c for c in ok_df.columns if c not in ("job", "status")]
    else:
        # quick, informative default – iptm-based metrics only
        metric_list = [c for c in ok_df.columns
                       if c.lower().startswith("iptm")]

    for m in metric_list:
        plot_metric_dot(ok_df, m, dot_dir / f"{m}_{args.style}.html",
                        style=args.style, annotate=args.labels)
    if all_plots:
        plot_combined_dots(ok_df, plots_dir / "all_metrics_dots.html",
                           style=args.style, top_n=None, annotate=args.labels)
    vprint("   combined dot plot finished")

    # interactive scatter dashboard
    
    plot_interactive_scatter(
        ok_df,
        plots_dir / "scatter_dashboard.html",
        init_x="confidence_score",
        init_y="complex_plddt",
        init_size="iptm",
        init_color="complex_iplddt"
    )

    vprint("📊  All plots in", plots_dir)
    print("📊  All plots in", plots_dir)
    import sys
    sys.exit(0)

if __name__ == "__main__":
    main()