# project_layout.py
from pathlib import Path

ROOT = Path.cwd()            # or Path("/absolute/path/to/root")
DIRS = {
    "analytics": ROOT / "analytics",
    "inputs":    ROOT / "inputs",
    "plots":     ROOT / "plots",
    "results":   ROOT / "results",
    "slurm":     ROOT / "slurm",
}

def result_dir(job_name: str) -> Path:
    """Return path for an individual job under results/."""
    return DIRS["results"] / job_name