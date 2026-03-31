from __future__ import annotations

import json
import os
import pathlib
import sys

import pytest
import yaml


def test_classify_sequence_modalities(wrapper_mod):
    assert wrapper_mod.classify_sequence("ACGTNN") == "dna"
    assert wrapper_mod.classify_sequence("ACGUNN") == "rna"
    assert wrapper_mod.classify_sequence("MSTNPKPQR") == "protein"


def test_make_yaml_supports_mixed_modalities(tmp_path, monkeypatch, wrapper_mod):
    rna_fa = tmp_path / "bait.rna"
    rna_fa.write_text(">rna1\nACGUACGUACGU\n")

    smi = tmp_path / "ligand.smi"
    smi.write_text("CCO ETH\n")

    seqio, seq, seqrecord = wrapper_mod._bio()

    def fake_records(ids, **_kwargs):
        for uid in ids:
            yield seqrecord(seq("MSTNPKPQRITF"), id=f"sp|{uid}|TEST", description="")

    monkeypatch.setattr(wrapper_mod, "_records_from_uniprot", fake_records)
    monkeypatch.chdir(tmp_path)

    out = wrapper_mod.make_yaml(
        f"P11111 2:SEP,{rna_fa},ATP,{smi},ACGTACGTACGTACGTACGTACGTA",
        [],
        {},
    )

    doc = yaml.safe_load(out.read_text())
    kinds = [next(iter(item.keys())) for item in doc["sequences"]]
    assert kinds == ["protein", "rna", "dna", "ligand", "ligand"]

    protein = doc["sequences"][0]["protein"]
    assert protein["modifications"] == [{"position": 2, "ccd": "SEP"}]
    assert len(protein["id"]) <= 5
    assert doc["sequences"][2]["dna"]["sequence"] == "ACGTACGTACGTACGTACGTACGTA"
    assert doc["sequences"][3]["ligand"]["ccd"] == "ATP"
    assert doc["sequences"][4]["ligand"]["smiles"] == "CCO"


def test_yaml_to_af3_maps_all_sequence_kinds(wrapper_mod):
    payload = {
        "version": 1,
        "sequences": [
            {"protein": {"id": "P1", "sequence": "MSTN", "modifications": [{"position": 2, "ccd": "SEP"}]}},
            {"rna": {"id": "R1", "sequence": "ACGU", "modifications": [{"position": 1, "ccd": "PSU"}]}},
            {"dna": {"id": "D1", "sequence": "ACGT", "modifications": [{"position": 3, "ccd": "5MC"}]}},
            {"ligand": {"id": "ATP", "ccd": "ATP"}},
            {"ligand": {"id": "MG", "ccd": "MG"}},
        ],
    }

    af3 = wrapper_mod.yaml_to_af3(payload, "job1")
    seqs = af3[0]["sequences"]

    assert seqs[0]["proteinChain"]["modifications"][0]["ptmType"] == "CCD_SEP"
    assert seqs[1]["rnaSequence"]["modifications"][0]["modificationType"] == "CCD_PSU"
    assert seqs[2]["dnaSequence"]["modifications"][0]["modificationType"] == "CCD_5MC"
    assert seqs[3] == {"ligand": {"ligand": "CCD_ATP", "count": 1}}
    assert seqs[4] == {"ion": {"ion": "MG", "count": 1}}


def test_gpu_and_flags_cover_size_and_modifiers(wrapper_mod):
    constraint, mem, flags, est = wrapper_mod.gpu_and_flags(1200, 6, 2, 25)
    assert constraint is not None
    assert "--sampling_steps" in flags
    assert "--recycling_steps" in flags
    assert est > 0

    constraint, mem, flags, est = wrapper_mod.gpu_and_flags(300, 1, 0, 0)
    assert constraint is not None
    assert mem.endswith("M")
    assert est > 0

    constraint, mem, flags, est = wrapper_mod.gpu_and_flags(5000, 8, 4, 0)
    assert constraint is None


def test_array_constraint_intersects_job_classes(wrapper_mod):
    assert wrapper_mod.array_constraint(["g4|g3", "g4|g2|g3|g1", "g4|g2"]) == "g4"
    assert wrapper_mod.array_constraint(["g4|g3", "g4|g3"]) == "g4|g3"


def test_analysis_dependency_and_retry_logic_present(wrapper_mod):
    assert "dependency=afterany:{array_id}" in wrapper_mod.ANALYSIS_TEMPLATE
    assert "[retry] transient-looking failure" in wrapper_mod.ARRAY_TEMPLATE
    assert "OOM-like failure detected; not retrying." in wrapper_mod.ARRAY_TEMPLATE
    assert "PIPESTATUS[0]" in wrapper_mod.ARRAY_TEMPLATE
    assert "BOLTZ_ANALYSIS_APPTAINER_IMAGE" in wrapper_mod.ANALYSIS_TEMPLATE
    assert "--no-mount hostfs" in wrapper_mod.ANALYSIS_TEMPLATE
    assert '/bin/bash --noprofile --norc -lc' in wrapper_mod.ANALYSIS_TEMPLATE
    assert 'export PATH="__BIN_DIR__:$PATH"; exec python -I "$@"' in wrapper_mod.ANALYSIS_TEMPLATE
    assert "__BIN_DIR__" in wrapper_mod.ANALYSIS_TEMPLATE
    assert "[requeue] transient runtime failure before prediction" in wrapper_mod.ARRAY_TEMPLATE
    assert "from boltz.main import cli" in wrapper_mod.ARRAY_TEMPLATE
    assert "ExcNodeList" in wrapper_mod.ARRAY_TEMPLATE
    assert "#SBATCH --requeue" in wrapper_mod.ARRAY_TEMPLATE
    assert 'unset PYTHONPATH PYTHONHOME PYTHONUSERBASE' in wrapper_mod.ARRAY_TEMPLATE
    assert 'unset APPTAINERENV_PYTHONPATH APPTAINERENV_PYTHONHOME APPTAINERENV_PYTHONUSERBASE' in wrapper_mod.ARRAY_TEMPLATE
    assert 'export APPTAINERENV_PYTHONPATH="$BOLTZ_CONTAINER_SITEPKGS"' in wrapper_mod.ARRAY_TEMPLATE
    assert 'export APPTAINERENV_TMPDIR="$RUNTIME_TMP"' in wrapper_mod.ARRAY_TEMPLATE
    assert 'export APPTAINERENV_PYTHONDONTWRITEBYTECODE=1' in wrapper_mod.ARRAY_TEMPLATE
    assert '--home "$RUNTIME_HOME"' in wrapper_mod.ARRAY_TEMPLATE
    assert 'JOB_DIR="$(dirname "$YAML")"' in wrapper_mod.ARRAY_TEMPLATE
    assert '--bind "$JOB_DIR:$JOB_DIR"' in wrapper_mod.ARRAY_TEMPLATE
    assert '--no-mount hostfs --nv' in wrapper_mod.ARRAY_TEMPLATE
    assert '${PYTHONPATH:+:$PYTHONPATH}' not in wrapper_mod.ARRAY_TEMPLATE
    assert 'python -I - <<' in wrapper_mod.ARRAY_TEMPLATE
    assert 'unset PYTHONPATH PYTHONHOME PYTHONUSERBASE' in wrapper_mod.ANALYSIS_TEMPLATE
    assert 'unset APPTAINERENV_PYTHONPATH APPTAINERENV_PYTHONHOME APPTAINERENV_PYTHONUSERBASE' in wrapper_mod.ANALYSIS_TEMPLATE
    assert 'export APPTAINERENV_TMPDIR="$RUNTIME_TMP"' in wrapper_mod.ANALYSIS_TEMPLATE
    assert 'export APPTAINERENV_PYTHONDONTWRITEBYTECODE=1' in wrapper_mod.ANALYSIS_TEMPLATE
    assert '--home "$RUNTIME_HOME"' in wrapper_mod.ANALYSIS_TEMPLATE
    assert '[preflight] analysis runtime import check failed:' in wrapper_mod.ANALYSIS_TEMPLATE
    assert '[requeue] transient runtime failure before analysis' in wrapper_mod.ANALYSIS_TEMPLATE
    assert 'import pandas as pd' in wrapper_mod.ANALYSIS_TEMPLATE
    assert '${PYTHONPATH:+:$PYTHONPATH}' not in wrapper_mod.ANALYSIS_TEMPLATE


def test_container_defaults_point_at_current_image(wrapper_mod):
    assert str(wrapper_mod.CONTAINER_IMAGE).endswith("/containers/current")
    assert str(wrapper_mod.ANALYSIS_CONTAINER_IMAGE).endswith("/containers/current")


def test_shell_wrappers_default_to_current_image():
    repo_root = pathlib.Path(__file__).resolve().parents[1]
    screen_text = (repo_root / "boltz-screen.sh").read_text()
    fetch_text = (repo_root / "boltz-fetch-ptms.sh").read_text()
    analysis_text = (repo_root / "boltz-analysis.sh").read_text()

    assert 'DEFAULT_IMAGE="$ROOT_DIR/containers/current"' in screen_text
    assert 'PREP_IMAGE="${BOLTZ_PREPARE_IMAGE:-$DEFAULT_IMAGE}"' in screen_text
    assert 'unset PYTHONPATH PYTHONHOME PYTHONUSERBASE' in screen_text
    assert 'unset APPTAINERENV_PYTHONPATH APPTAINERENV_PYTHONHOME APPTAINERENV_PYTHONUSERBASE' in screen_text
    assert 'runtime_base_default()' in screen_text
    assert 'export APPTAINERENV_TMPDIR="$RUNTIME_TMP"' in screen_text
    assert 'export APPTAINERENV_PYTHONDONTWRITEBYTECODE=1' in screen_text
    assert '--home "$RUNTIME_HOME"' in screen_text
    assert 'exec python -I "$@"' in screen_text
    assert 'DEFAULT_IMAGE="$ROOT_DIR/containers/current"' in fetch_text
    assert 'PREP_IMAGE="${BOLTZ_PREPARE_IMAGE:-$DEFAULT_IMAGE}"' in fetch_text
    assert 'unset PYTHONPATH PYTHONHOME PYTHONUSERBASE' in fetch_text
    assert 'unset APPTAINERENV_PYTHONPATH APPTAINERENV_PYTHONHOME APPTAINERENV_PYTHONUSERBASE' in fetch_text
    assert 'runtime_base_default()' in fetch_text
    assert 'export APPTAINERENV_TMPDIR="$RUNTIME_TMP"' in fetch_text
    assert 'export APPTAINERENV_PYTHONDONTWRITEBYTECODE=1' in fetch_text
    assert '--home "$RUNTIME_HOME"' in fetch_text
    assert 'exec python -I "$@"' in fetch_text
    assert 'ANALYSIS_IMAGE="${BOLTZ_ANALYSIS_APPTAINER_IMAGE:-${BOLTZ_APPTAINER_IMAGE:-$DEFAULT_IMAGE}}"' in analysis_text
    assert 'unset PYTHONPATH PYTHONHOME PYTHONUSERBASE' in analysis_text
    assert 'unset APPTAINERENV_PYTHONPATH APPTAINERENV_PYTHONHOME APPTAINERENV_PYTHONUSERBASE' in analysis_text
    assert 'runtime_base_default()' in analysis_text
    assert 'TARGET_BIND=()' in analysis_text
    assert '[ -n "${1:-}" ] && [ "${1#-}" = "$1" ]' in analysis_text
    assert 'cmd=(' in analysis_text
    assert '[ "${#TARGET_BIND[@]}" -gt 0 ]' in analysis_text
    assert 'export APPTAINERENV_TMPDIR="$RUNTIME_TMP"' in analysis_text
    assert 'export APPTAINERENV_PYTHONDONTWRITEBYTECODE=1' in analysis_text
    assert '--home "$RUNTIME_HOME"' in analysis_text
    assert 'cmd+=("${TARGET_BIND[@]}")' in analysis_text
    assert 'exec python -I "$@"' in analysis_text
    assert '${PYTHONPATH:+:$PYTHONPATH}' not in analysis_text


def test_apptainer_env_strips_host_python_settings(tmp_path, monkeypatch, wrapper_mod):
    monkeypatch.setattr(wrapper_mod, "CONTAINER_SITEPKGS", tmp_path / "sitepkgs")
    wrapper_mod.CONTAINER_SITEPKGS.mkdir()
    monkeypatch.setenv("PYTHONPATH", "/tmp/user/python")
    monkeypatch.setenv("PYTHONHOME", "/tmp/user/home")
    monkeypatch.setenv("PYTHONUSERBASE", "/tmp/user/base")
    monkeypatch.setenv("APPTAINERENV_PYTHONPATH", "/tmp/apptainer/python")
    monkeypatch.setenv("APPTAINERENV_PYTHONHOME", "/tmp/apptainer/home")
    monkeypatch.setenv("APPTAINERENV_PYTHONUSERBASE", "/tmp/apptainer/base")

    env = wrapper_mod._apptainer_env()

    assert env["APPTAINERENV_PYTHONPATH"] == str(wrapper_mod.CONTAINER_SITEPKGS)
    assert env["APPTAINERENV_PYTHONNOUSERSITE"] == "1"
    assert "PYTHONPATH" not in env
    assert "PYTHONHOME" not in env
    assert "PYTHONUSERBASE" not in env
    assert "APPTAINERENV_PYTHONHOME" not in env
    assert "APPTAINERENV_PYTHONUSERBASE" not in env


def test_parse_list_and_assert_token_validation(tmp_path, wrapper_mod):
    fasta = tmp_path / "x.fa"
    fasta.write_text(">x\nMSTNPKPQR\n")

    txt = tmp_path / "list.txt"
    txt.write_text("# comment\nA\n\nB\n")

    assert wrapper_mod._parse_list(str(fasta)) == [str(fasta)]
    assert wrapper_mod._parse_list(str(txt)) == ["A", "B"]

    wrapper_mod._assert_token_ok(str(fasta))
    with pytest.raises(SystemExit, match="Referenced path not found"):
        wrapper_mod._assert_token_ok(str(tmp_path / "missing.fa"))


def test_main_single_run_affinity_and_slurm(tmp_path, monkeypatch, wrapper_mod):
    fasta = tmp_path / "complex.fa"
    fasta.write_text(">protA\nMSTNPKPQRITFVKSHFSRQ\n")

    calls = []

    class FakeResult:
        def __init__(self, stdout="Submitted batch job 12345\n", returncode=0, stderr=""):
            self.stdout = stdout
            self.returncode = returncode
            self.stderr = stderr

    def fake_run(cmd, *args, **kwargs):
        calls.append(cmd)
        return FakeResult()

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(wrapper_mod.subprocess, "run", fake_run)
    monkeypatch.setattr(sys, "argv", ["boltz-2_wrapper.py", str(fasta), "--affinity", "L"])

    wrapper_mod.main()

    yaml_path = next(tmp_path.glob("boltz_*.yaml"))
    doc = yaml.safe_load(yaml_path.read_text())
    assert doc["properties"] == [{"affinity": {"binder": "L"}}]

    slurm_path = pathlib.Path(calls[0][1])
    assert slurm_path.is_file()
    assert "Submitted batch job 12345" in FakeResult().stdout


def test_main_local_passes_potential_flags(tmp_path, monkeypatch, wrapper_mod):
    fasta = tmp_path / "complex.fa"
    fasta.write_text(">protA\nMSTNPKPQRITFVKSHFSRQ\n")
    observed = {}

    def fake_local(yaml_path, outdir, flags):
        observed["yaml"] = yaml_path
        observed["outdir"] = outdir
        observed["flags"] = flags

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(wrapper_mod, "run_local", fake_local)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "boltz-2_wrapper.py",
            str(fasta),
            "--local",
            "--boltz",
            "--use_potentials --sampling_steps 300",
        ],
    )

    wrapper_mod.main()

    assert observed["yaml"].is_file()
    assert "--use_potentials" in observed["flags"]
    assert "300" in observed["flags"]
