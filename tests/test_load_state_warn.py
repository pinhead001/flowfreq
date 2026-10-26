"""``load_state(..., warn=False)``: opt-in sets load without the logged advisory (#214).

A coverage sweep over the 41 ``nss_snapshot`` jurisdictions used to log 41 warnings. With
``warn=False`` the advisory is still set on the library, and logged at DEBUG instead.
"""

from __future__ import annotations

import importlib.util
import json
import logging
from pathlib import Path

import pytest

from flowfreq.regression.library import NSS_SNAPSHOT_SUBDIR, EquationsUnavailable, load_state

LOGGER = "flowfreq.regression.library"
REPO_ROOT = Path(__file__).resolve().parent.parent

EQ = {
    "region_code": "GC2015",
    "aep": 0.01,
    "intercept": 5.0,
    "variables": [{"code": "DRNAREA"}],
    "coefficients": [0.7],
    "citation": {"publication": "Synthetic", "table": "T1"},
}


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    path = tmp_path / NSS_SNAPSHOT_SUBDIR / "TN.json"
    path.parent.mkdir(parents=True)
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "state": "TN",
                "status": "nss_snapshot",
                "source_reports": [{"citation": "Synthetic report"}],
                "equations": [EQ],
                "retrieved_utc": "2026-10-25T17:20:40+00:00",
            }
        ),
        encoding="utf-8",
    )
    return tmp_path


def _records(caplog: pytest.LogCaptureFixture, level: int) -> list[logging.LogRecord]:
    return [r for r in caplog.records if r.name == LOGGER and r.levelno == level]


def test_default_still_warns(data_dir: Path, caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.DEBUG, logger=LOGGER):
        lib = load_state("TN", data_dir=data_dir, allow_nss_snapshot=True)
    warnings = _records(caplog, logging.WARNING)
    assert len(warnings) == 1 and warnings[0].getMessage() == lib.advisory


def test_warn_false_keeps_advisory_and_logs_at_debug(
    data_dir: Path, caplog: pytest.LogCaptureFixture
) -> None:
    with caplog.at_level(logging.DEBUG, logger=LOGGER):
        lib = load_state("TN", data_dir=data_dir, allow_nss_snapshot=True, warn=False)
    assert lib.nss_snapshot and lib.advisory.startswith("NSS SNAPSHOT (not report-verified)")
    assert _records(caplog, logging.WARNING) == []
    debug = _records(caplog, logging.DEBUG)
    assert [r.getMessage() for r in debug] == [lib.advisory]


def test_warn_is_keyword_only(data_dir: Path) -> None:
    with pytest.raises(TypeError):
        load_state("TN", data_dir, False, False, True, False)  # type: ignore[misc]


def test_warn_false_without_any_file_still_raises(tmp_path: Path) -> None:
    with pytest.raises(EquationsUnavailable, match="no regression library file"):
        load_state("TN", data_dir=tmp_path, allow_nss_snapshot=True, warn=False)


def test_coverage_report_does_not_warn(caplog: pytest.LogCaptureFixture) -> None:
    spec = importlib.util.spec_from_file_location(
        "gen_cov_warn", REPO_ROOT / "tools" / "gen_regression_coverage.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        mod.render()
    assert _records(caplog, logging.WARNING) == []
