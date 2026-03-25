from __future__ import annotations

import json
import sys

import yaml


def test_screening_submission_builds_expected_files(tmp_path, monkeypatch, wrapper_mod):
    bait_fa = tmp_path / "bait.fa"
    bait_fa.write_text(">bait\nMSTNPKPQRITFVKSHFSRQDILD\n")

    rna_fa = tmp_path / "target.rna"
    rna_fa.write_text(">rna\nACGUACGUACGUACGUACGUACGU\n")

    bait_txt = tmp_path / "bait.txt"
    bait_txt.write_text(f"{bait_fa}\nMG\n")

    screen_txt = tmp_path / "screen.txt"
    screen_txt.write_text(f"{rna_fa}\nATP\n")

    chain_map = tmp_path / "chain_map.txt"
    chain_map.write_text("0=Bait\n1=Target\n2=Mg\n")

    calls = []

    class FakeResult:
        def __init__(self, stdout):
            self.stdout = stdout

    def fake_run(cmd, *args, **kwargs):
        calls.append(cmd)
        if len(calls) == 1:
            return FakeResult("Submitted batch job 1001\n")
        return FakeResult("Submitted batch job 1002\n")

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(wrapper_mod.subprocess, "run", fake_run)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "boltz-2_wrapper.py",
            "--bait",
            str(bait_txt),
            "--screen",
            str(screen_txt),
            "--chain_map",
            str(chain_map),
            "-n",
            "screen_demo",
        ],
    )

    wrapper_mod.main()

    root = tmp_path / "boltz_screen_demo"
    assert (root / "jobs.list").is_file()
    assert (root / "array.slurm").is_file()
    assert (root / "analysis.slurm").is_file()
    assert len(list((root / "results").glob("*"))) == 2

    jobs = (root / "jobs.list").read_text().strip().splitlines()
    assert len(jobs) == 2
    assert "#SBATCH --array=1-2" in (root / "array.slurm").read_text()
    assert "screen_demo_ana" in (root / "analysis.slurm").read_text()
    assert calls[0][0] == "sbatch"
    assert calls[1][0] == "sbatch"

    af3_jsons = sorted((root / "AF3_JSON").glob("*.json"))
    assert len(af3_jsons) == 2

    first_yaml = next((root / "results").glob("*/*.yaml"))
    doc = yaml.safe_load(first_yaml.read_text())
    assert doc["version"] == 1


def test_screening_local_mode_skips_sbatch(tmp_path, monkeypatch, wrapper_mod):
    bait_fa = tmp_path / "bait.fa"
    bait_fa.write_text(">bait\nMSTNPKPQRITFVKSHFSRQDILD\n")
    screen_fa = tmp_path / "screen.fa"
    screen_fa.write_text(">screen\nMQQQQQQQQQQQQQQQQQQQQQQQ\n")

    bait_txt = tmp_path / "bait.txt"
    bait_txt.write_text(f"{bait_fa}\n")
    screen_txt = tmp_path / "screen.txt"
    screen_txt.write_text(f"{screen_fa}\n")

    observed = []

    def fake_local(yaml_path, outdir, flags):
        observed.append((yaml_path, outdir, tuple(flags)))

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(wrapper_mod, "run_local", fake_local)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "boltz-2_wrapper.py",
            "--bait",
            str(bait_txt),
            "--screen",
            str(screen_txt),
            "--local",
            "-n",
            "local_demo",
        ],
    )

    wrapper_mod.main()

    assert len(observed) == 1
    assert observed[0][0].is_file()
    assert not (tmp_path / "boltz_local_demo" / "jobs.list").exists()


def test_screening_all_skipped_exits_cleanly(tmp_path, monkeypatch, wrapper_mod):
    bait_fa = tmp_path / "bait.fa"
    bait_fa.write_text(">bait\nMSTNPKPQRITFVKSHFSRQDILD\n")
    screen_fa = tmp_path / "screen.fa"
    screen_fa.write_text(">screen\nMQQQQQQQQQQQQQQQQQQQQQQQ\n")

    bait_txt = tmp_path / "bait.txt"
    bait_txt.write_text(f"{bait_fa}\n")
    screen_txt = tmp_path / "screen.txt"
    screen_txt.write_text(f"{screen_fa}\n")

    def fake_gpu(*_args, **_kwargs):
        return None, None, None, 999.0

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(wrapper_mod, "gpu_and_flags", fake_gpu)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "boltz-2_wrapper.py",
            "--bait",
            str(bait_txt),
            "--screen",
            str(screen_txt),
            "-n",
            "skip_demo",
        ],
    )

    wrapper_mod.main()

    root = tmp_path / "boltz_skip_demo"
    skipped = root / "skipped_jobs.txt"
    assert skipped.is_file()
    assert "needs" in skipped.read_text()
    assert not (root / "jobs.list").exists()
