from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd
import yaml


def test_analysis_main_writes_summary_and_handles_failed_jobs(tmp_path, monkeypatch, analysis_mod):
    root = tmp_path / "screen_root"
    ok_dir = root / "results" / "001_demo_Q1"
    failed_dir = root / "results" / "002_demo_Q2"
    ok_dir.mkdir(parents=True)
    failed_dir.mkdir(parents=True)

    yaml_doc = {
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": "MSTNPKPQR", "modifications": [{"position": 2, "ccd": "SEP"}]}},
            {"rna": {"id": "B", "sequence": "ACGUACGU"}},
            {"ligand": {"id": "ATP", "ccd": "ATP"}},
        ],
    }
    (ok_dir / "001_demo_Q1.yaml").write_text(yaml.safe_dump(yaml_doc, sort_keys=False))
    (failed_dir / "002_demo_Q2.yaml").write_text(yaml.safe_dump(yaml_doc, sort_keys=False))

    conf = {
        "iptm": 0.82,
        "ptm": 0.74,
        "chains_ptm": {"0": 0.91},
        "pair_chains_iptm": {"0": {"1": 0.61}},
    }
    (ok_dir / "confidence_model_0.json").write_text(json.dumps(conf))
    np.savez(ok_dir / "pae_model_0.npz", arr=np.array([[1.0, 2.0], [2.0, 1.0]]))
    (ok_dir / "slurm_12345.out").write_text("ok\n")
    (failed_dir / "slurm_12346.out").write_text("failed\n")

    def fake_run(cmd, *args, **kwargs):
        class Result:
            def __init__(self, stdout):
                self.stdout = stdout

        if cmd[0] == "jobinfo":
            return Result(
                "Mem reserved : 8000M\n"
                "Max Mem (Node/step) : 7.5G\n"
                "Used walltime : 00:10:00\n"
            )
        raise AssertionError(f"unexpected subprocess call: {cmd}")

    monkeypatch.setattr(analysis_mod, "uniprot_label", lambda uid: uid)
    monkeypatch.setattr(analysis_mod.subprocess, "run", fake_run)
    monkeypatch.setattr(analysis_mod, "_save", lambda _fig, path, png_scale=2: path.write_text("stub"))
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["boltz_analysis.py", str(root), "--no-labels"])

    try:
        analysis_mod.main()
    except SystemExit as exc:
        assert exc.code == 0

    summary = pd.read_csv(root / "summary_metrics.csv")
    assert set(summary["status"]) == {"ok", "failed"}
    ok_row = summary.loc[summary["status"] == "ok"].iloc[0]
    failed_row = summary.loc[summary["status"] == "failed"].iloc[0]

    assert ok_row["seqs"] == 2
    assert ok_row["ligands"] == 1
    assert ok_row["mods"] == 1
    assert ok_row["iptm"] == 0.82
    assert failed_row["residues"] == 17

    assert (root / "analytics" / "slurm_metrics.csv").is_file()
    assert (root / "plots").is_dir()
    assert (ok_dir / "CHIMERAX_001_demo_Q1_analysis.cxc").is_file()
