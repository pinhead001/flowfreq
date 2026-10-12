"""The Washington SIR 2016-5118 regression dataset (tools/build_wa_gls_dataset.py).

Literals are the report's own: Table 5 counts and Table 8 values, the latter as
already pinned in tests/test_regression_washington.py (PR #65).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

CSV = Path(__file__).resolve().parents[1] / "data" / "regression_study" / "wa_2016_5118_sites.csv"


@pytest.fixture(scope="module")
def sites() -> pd.DataFrame:
    return pd.read_csv(CSV, dtype={"site_no": str})


def test_one_row_per_table5_gage(sites):
    assert len(sites) == 641
    assert sites["site_no"].is_unique
    assert set(sites["region"].dropna().astype(int)) == {1, 2, 3, 4}


def test_used_counts_by_region(sites):
    used = sites[sites["used"]].groupby(sites["region"].astype("Int64")).size().to_dict()
    # Region 4: the 139 calibration gages of PR #65.
    assert used == {1: 93, 2: 89, 3: 142, 4: 139}


@pytest.mark.parametrize(
    "site, aep, line",
    [
        ("12024400", 0.01, (7290, 3840, 6840, 3720, 12600)),
        ("12024400", 0.002, (10100, 5150, 9210, 4280, 19800)),
        ("12020800", 0.002, (18500, 4950, 16100, 7690, 33700)),
    ],
)
def test_table8_lines(sites, site, aep, line):
    r = sites.set_index("site_no").loc[site]
    got = tuple(r[f"{k}_{aep}"] for k in ("q", "qreg", "qw", "ci_lo", "ci_hi"))
    assert got == line


def test_table8_covers_washington_only(sites):
    # Idaho and Oregon gages used in the regression have no Table 8 rows.
    missing = sites[sites["used"] & sites["q_0.01"].isna()]
    assert set(missing["state"]) == {"ID", "OR"}
    assert len(missing) == 51


def test_recovered_lp3_reproduces_station_quantiles(sites):
    used = sites[sites["used"] & sites["q_0.01"].notna()]
    # Three-significant-figure rounding allows ~0.5 %; the exceptions are listed.
    off = set(used.loc[used["fit_max_rel_err"] > 0.006, "site_no"])
    # 14017200's printed 1-percent value (1674) breaks its own LP3 pattern
    # (neighbours 1330 and 3170): a Table 8 typo. 12144500 is at 1.4 %.
    assert off == {"14017200", "12144500"}
    assert np.isfinite(used[["mu", "sd", "skew_w"]].to_numpy()).all()
