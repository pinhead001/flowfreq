"""Provisional (flowfreq-developed, non-USGS) equation sets.

Oregon regression update plan, task 1
(``docs/superpowers/plans/2026-10-12-oregon-regression-update.md``): a provisional
set lives under ``provisional/`` and is returned only on opt-in, so it can never
change the default ``load_state``.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pytest

from flowfreq.regression import available_states, load_state, provisional_states
from flowfreq.regression.library import DATA_DIR, PROVISIONAL_SUBDIR

EQ = {
    "region_code": "R1",
    "aep": 0.01,
    "intercept": 2.0,
    "variables": [{"code": "DRNAREA"}],
    "coefficients": [0.7],
    "citation": {"publication": "Synthetic", "table": "T1"},
}


def _write(path: Path, status: str, intercept: float, state: str = "OR") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "schema_version": 1,
        "state": state,
        "status": status,
        "equations": [{**EQ, "intercept": intercept}],
    }
    path.write_text(json.dumps(data), encoding="utf-8")


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    _write(tmp_path / "OR.json", "verified", 2.0)
    _write(tmp_path / PROVISIONAL_SUBDIR / "OR.json", "provisional", 3.0)
    return tmp_path


def test_default_lookup_ignores_the_provisional_set(data_dir):
    lib = load_state("OR", data_dir=data_dir)
    assert lib.status == "verified" and not lib.provisional and lib.advisory == ""
    assert lib.equation("R1", 0.01).intercept == 2.0


def test_opt_in_returns_the_provisional_set(data_dir, caplog):
    with caplog.at_level(logging.WARNING, logger="flowfreq.regression.library"):
        lib = load_state("OR", data_dir=data_dir, allow_provisional=True)
    assert lib.provisional and lib.equation("R1", 0.01).intercept == 3.0
    assert lib.advisory.startswith("PROVISIONAL (not USGS)")
    assert "PROVISIONAL" in caplog.text


def test_opt_in_without_a_provisional_set_returns_the_published_one(tmp_path):
    _write(tmp_path / "OR.json", "verified", 2.0)
    lib = load_state("OR", data_dir=tmp_path, allow_provisional=True)
    assert lib.status == "verified" and lib.advisory == ""


def test_provisional_status_cannot_replace_the_published_file(tmp_path):
    _write(tmp_path / "OR.json", "provisional", 3.0)
    with pytest.raises(ValueError, match="belongs under provisional/"):
        load_state("OR", data_dir=tmp_path)


@pytest.mark.parametrize("status, state", [("verified", "OR"), ("provisional", "WA")])
def test_provisional_file_must_say_so(tmp_path, status, state):
    _write(tmp_path / PROVISIONAL_SUBDIR / "OR.json", status, 3.0, state=state)
    with pytest.raises(ValueError, match="must have status 'provisional'"):
        load_state("OR", data_dir=tmp_path, allow_provisional=True)


def test_listings_keep_the_two_apart(data_dir):
    assert available_states(data_dir) == {"OR": "verified"}
    assert provisional_states(data_dir) == ["OR"]
    assert provisional_states(data_dir / "nowhere") == []


def test_packaged_provisional_sets_are_valid():
    for code in provisional_states():
        lib = load_state(code, allow_provisional=True)
        assert lib.provisional and lib.equations
        assert available_states()[code] != "provisional"
    assert all(s != "provisional" for s in available_states(DATA_DIR).values())
