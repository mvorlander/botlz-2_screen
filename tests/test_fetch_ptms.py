from __future__ import annotations

import sys


class FakeResponse:
    def __init__(self, *, json_data=None, text="", status_code=200):
        self._json = json_data
        self.text = text
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"http {self.status_code}")

    def json(self):
        return self._json


def test_fetch_ptm_records_maps_and_filters(monkeypatch, ptm_mod):
    ptm_mod._fetch_sequence.cache_clear()

    def fake_get(url, *args, **kwargs):
        if "proteomics/ptm" in url:
            return FakeResponse(
                json_data={
                    "features": [
                        {
                            "description": "Phosphoserine",
                            "begin": 10,
                            "peptide": {"start": 10, "sequence": "ASDF"},
                            "ptms": [
                                {
                                    "name": "Phosphoserine",
                                    "position": 2,
                                    "sources": ["PRIDE"],
                                }
                            ],
                        },
                        {
                            "description": "Acetylation",
                            "begin": 20,
                            "peptide": {"start": 20, "sequence": "KAAA"},
                            "ptms": [
                                {
                                    "name": "Acetylation",
                                    "position": 1,
                                    "sources": ["PSP"],
                                }
                            ],
                        },
                    ]
                }
            )
        if "uniprotkb" in url and url.endswith(".fasta"):
            return FakeResponse(text=">sp|Q1|TEST\nMSTNPKPQRITFVKSHFSRQ")
        raise AssertionError(f"unexpected URL {url}")

    monkeypatch.setattr(ptm_mod.requests, "get", fake_get)

    recs = ptm_mod.fetch_ptm_records("Q1", keep_types=["phospho", "acetyl"])
    assert ("ASDF:2:Phosphoserine", "Q1 11:SEP", "PRIDE") in recs
    assert ("KAAA:1:Acetylation", "Q1 20:ALY", "PSP") in recs


def test_fetch_ptm_cli_writes_outputs(tmp_path, monkeypatch, ptm_mod):
    monkeypatch.setattr(
        ptm_mod,
        "fetch_ptm_records",
        lambda acc, **kwargs: [("pep:2:Phosphoserine", f"{acc} 11:SEP", "PRIDE")],
    )
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        sys,
        "argv",
        ["boltz_fetch_ptms.py", "Q1,Q2", "-o", "ptms.tsv"],
    )

    ptm_mod._cli()

    assert (tmp_path / "ptms.tsv").is_file()
    assert (tmp_path / "ptms.boltz.txt").is_file()
    assert "Q1 11:SEP" in (tmp_path / "ptms.boltz.txt").read_text()
