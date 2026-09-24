"""Tests for flowfreq.catalog and the generated coverage report."""

from __future__ import annotations

import importlib.util

import pandas as pd
import pytest

from flowfreq.catalog import COLUMNS, SEED_PATH, catalog_row, load_catalog, select
from tests.fixtures.paths import REPO_ROOT


def test_seed_loads_with_all_columns():
    df = load_catalog(SEED_PATH)
    assert tuple(df.columns) == COLUMNS and len(df) == 3
    assert df["site_no"].str.len().min() >= 8


def _write(tmp_path, body, header="site_no,site_name,state,drainage_area_sqmi\n"):
    p = tmp_path / "c.csv"
    p.write_text(header + body)
    return p


@pytest.mark.parametrize(
    "body, match",
    [
        ("12345678,A,WA,10\n12345678,B,WA,10\n", "duplicate"),
        ("1234,A,WA,10\n", "digits"),
        ("12345678,A,ZZ,10\n", "unknown state"),
        ("12345678,A,WA,0\n", "positive"),
    ],
)
def test_validation(tmp_path, body, match):
    with pytest.raises(ValueError, match=match):
        load_catalog(_write(tmp_path, body))


def test_missing_required_column(tmp_path):
    with pytest.raises(ValueError, match="required"):
        load_catalog(_write(tmp_path, "12345678,A\n", header="site_no,site_name\n"))


def test_catalog_row_and_select():
    peaks = pd.DataFrame(
        {
            "water_year": [1990, 1991, 1995],
            "peak_flow_cfs": [1.0, 2.0, 3.0],
            "qualification_code": ["", "6", ""],
        }
    )
    row = catalog_row("12345678", "wa", peaks, site_name="Test", drainage_area_sqmi=5.0)
    assert (row["n_peaks"], row["first_water_year"], row["last_water_year"]) == (3, 1990, 1995)
    assert row["regulation_class"] == "regulated" and row["state"] == "WA"
    df = pd.DataFrame(
        [row, {**row, "site_no": "87654321", "regulation_class": "no_code_evidence", "n_peaks": 30}]
    )
    assert len(select(df, state="WA", exclude_regulated=True)) == 1
    assert len(select(df, min_years=10)) == 1


def test_regression_coverage_doc_is_current():
    """docs/REGRESSION_COVERAGE.md must match what the generator produces now."""
    spec = importlib.util.spec_from_file_location(
        "gen_cov", REPO_ROOT / "tools" / "gen_regression_coverage.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    committed = (REPO_ROOT / "docs" / "REGRESSION_COVERAGE.md").read_text(encoding="utf-8")
    assert committed == mod.render(), "stale: run python tools/gen_regression_coverage.py"
