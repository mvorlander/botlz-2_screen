from __future__ import annotations

import importlib.util
import pathlib
import uuid

import pytest


REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]


def _load_module(path: pathlib.Path):
    name = f"testmod_{path.stem.replace('-', '_')}_{uuid.uuid4().hex}"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def wrapper_mod():
    return _load_module(REPO_ROOT / "boltz-2_wrapper.py")


@pytest.fixture()
def analysis_mod():
    return _load_module(REPO_ROOT / "boltz_analysis.py")


@pytest.fixture()
def ptm_mod():
    return _load_module(REPO_ROOT / "boltz_fetch_ptms.py")
