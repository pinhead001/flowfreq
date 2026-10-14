"""USGS-published equation sets that NSS/StreamStats does not serve yet.

Status ``awaiting_nss`` lives under ``awaiting_nss/`` and, like a provisional set, is
returned only on opt-in (``allow_awaiting_nss=True``), so it never changes the default
``load_state``. Its advisory says it is USGS-published, not flowfreq-developed.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pytest

from flowfreq.regression import (
    available_states,
    awaiting_nss_states,
    load_state,
    provisional_states,
)
from flowfreq.regression.library import AWAITING_NSS_SUBDIR, PROVISIONAL_SUBDIR, STATUSES

EQ = {
    "region_code": "R1",
    "aep": 0.01,
    "intercept": 2.0,
    "variables": [{"code": "DRNAREA"}],
    "coefficients": [0.7],
    "citation": {"publication": "Synthetic", "table": "T1"},
}
REPORT = "Doe, J., 2026, Synthetic report: U.S. Geological Survey SIR 2026-0000"


def _write(path: Path, status: str, intercept: float, state: str = "CO") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "schema_version": 1,
        "state": state,
        "status": status,
        "source_reports": [{"citation": REPORT}],
        "equations": [{**EQ, "intercept": intercept}],
    }
    path.write_text(json.dumps(data), encoding="utf-8")


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    _write(tmp_path / "CO.json", "verified", 2.0)
    _write(tmp_path / AWAITING_NSS_SUBDIR / "CO.json", "awaiting_nss", 4.0)
    return tmp_path


def test_status_is_defined():
    assert "awaiting_nss" in STATUSES


def test_default_lookup_ignores_the_awaiting_nss_set(data_dir):
    lib = load_state("CO", data_dir=data_dir)
    assert lib.status == "verified" and not lib.awaiting_nss and lib.advisory == ""
    assert lib.equation("R1", 0.01).intercept == 2.0


def test_provisional_opt_in_does_not_return_it(data_dir):
    lib = load_state("CO", data_dir=data_dir, allow_provisional=True)
    assert lib.status == "verified"


def test_opt_in_returns_it_with_a_usgs_published_warning(data_dir, caplog):
    with caplog.at_level(logging.WARNING, logger="flowfreq.regression.library"):
        lib = load_state("CO", data_dir=data_dir, allow_awaiting_nss=True)
    assert lib.awaiting_nss and not lib.provisional
    assert lib.equation("R1", 0.01).intercept == 4.0
    assert lib.advisory.startswith("NOT YET IN STREAMSTATS (USGS-published)")
    assert REPORT in lib.advisory
    assert "not yet served by NSS/StreamStats" in lib.advisory
    assert "load_state('CO')" in lib.advisory
    # Never the flowfreq-developed wording.
    assert "flowfreq" not in lib.advisory and "not USGS" not in lib.advisory
    assert "NOT YET IN STREAMSTATS" in caplog.text


def test_usgs_set_wins_over_a_provisional_one_when_both_opted_in(data_dir):
    _write(data_dir / PROVISIONAL_SUBDIR / "CO.json", "provisional", 3.0)
    both = load_state("CO", data_dir=data_dir, allow_provisional=True, allow_awaiting_nss=True)
    assert both.awaiting_nss and both.equation("R1", 0.01).intercept == 4.0
    assert load_state("CO", data_dir=data_dir, allow_provisional=True).provisional


def test_opt_in_without_a_file_returns_the_published_set(tmp_path):
    _write(tmp_path / "CO.json", "verified", 2.0)
    lib = load_state("CO", data_dir=tmp_path, allow_awaiting_nss=True)
    assert lib.status == "verified" and lib.advisory == ""


def test_status_cannot_replace_the_published_file(tmp_path):
    _write(tmp_path / "CO.json", "awaiting_nss", 4.0)
    with pytest.raises(ValueError, match="belongs under awaiting_nss/"):
        load_state("CO", data_dir=tmp_path)


@pytest.mark.parametrize(
    "status, state", [("verified", "CO"), ("provisional", "CO"), ("awaiting_nss", "UT")]
)
def test_awaiting_nss_file_must_say_so(tmp_path, status, state):
    _write(tmp_path / AWAITING_NSS_SUBDIR / "CO.json", status, 4.0, state=state)
    with pytest.raises(ValueError, match="must have status 'awaiting_nss'"):
        load_state("CO", data_dir=tmp_path, allow_awaiting_nss=True)


def test_awaiting_nss_status_is_refused_under_provisional(tmp_path):
    _write(tmp_path / PROVISIONAL_SUBDIR / "CO.json", "awaiting_nss", 4.0)
    with pytest.raises(ValueError, match="must have status 'provisional'"):
        load_state("CO", data_dir=tmp_path, allow_provisional=True)


def test_listings_keep_the_sets_apart(data_dir):
    assert available_states(data_dir) == {"CO": "verified"}
    assert awaiting_nss_states(data_dir) == ["CO"]
    assert provisional_states(data_dir) == []
    assert awaiting_nss_states(data_dir / "nowhere") == []


def test_packaged_awaiting_nss_sets_are_valid():
    assert "CO" in awaiting_nss_states()
    for code in awaiting_nss_states():
        lib = load_state(code, allow_awaiting_nss=True)
        assert lib.awaiting_nss and lib.equations and lib.source_reports
        assert available_states()[code] not in ("provisional", "awaiting_nss")
