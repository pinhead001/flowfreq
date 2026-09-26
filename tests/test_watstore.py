"""Tests for flowfreq.watstore, the WATSTORE (``I ASCI``) peak-file reader."""

from __future__ import annotations

import math

import pandas as pd
import pytest

from flowfreq.watstore import WATSTORE_COLUMNS, parse_watstore, read_watstore, separate_codes
from tests.fixtures.paths import SKIP_REASON, TESTDATA_AVAILABLE
from tests.fixtures.paths import testdata_path as _testdata

# Copied from vendor/peakfqr/inst/testdata/wymt_ffa_2022A_WATSTORE.TXT
# (06177720.00), plus one October peak and one gage-height-only record.
SAMPLE = """\
Z06177720.00                   USGS
H06177720.00    4731561051406003030021ST1006000214.7          2520
N06177720.00    West Fork Sullivan Creek near Richey MT
306177720.00    1972        103         7Bm    5.45  Bm
306177720.00    1974       0.00          Bm
306177720.00    19780322   10.0           2    4.75   1
306177720.00    19791015    250          48    3.00
306177720.00    19800601                       2.00
"""


class TestSeparateCodes:
    @pytest.mark.parametrize(
        "raw, expected",
        [("7Bm", "7,Bm"), ("  4Bm ", "4,Bm"), ("125", "1,2,5"), ("", ""), ("Bd", "Bd")],
    )
    def test_forms(self, raw, expected):
        assert separate_codes(raw) == expected


class TestParse:
    def test_sample(self):
        df = parse_watstore(SAMPLE)
        assert tuple(df.columns) == WATSTORE_COLUMNS
        assert df.attrs["station_names"] == {
            "06177720.00": "West Fork Sullivan Creek near Richey MT"
        }
        assert list(df["site_no"].unique()) == ["06177720.00"]
        # Year only: the year is the water year.
        first = df.iloc[0]
        assert (first["water_year"], first["peak_flow_cfs"]) == (1972, 103.0)
        assert first["qualification_code"] == "7,Bm"
        assert first["peak_date"] == pd.Timestamp(1972, 1, 1)
        # October onward belongs to the next water year.
        oct_row = df[df["peak_flow_cfs"] == 250.0].iloc[0]
        assert oct_row["water_year"] == 1980
        assert oct_row["qualification_code"] == "4,8"
        # A zero is a real value; a blank discharge is NaN (gage height only).
        assert df.loc[df["water_year"] == 1974, "peak_flow_cfs"].iloc[0] == 0.0
        assert math.isnan(df.iloc[-1]["peak_flow_cfs"])

    def test_duplicates_removed_and_logged(self, caplog):
        line = "306177720.00    19780322   10.0           2    4.75   1\n"
        df = parse_watstore(line + line)
        assert len(df) == 1
        assert "duplicate" in caplog.text

    def test_peak_without_year_raises(self):
        with pytest.raises(ValueError, match="no year"):
            parse_watstore("306177720.00            10.0\n")

    def test_unknown_record_logged(self, caplog):
        df = parse_watstore("X something\n306177720.00    1978       10.0\n")
        assert len(df) == 1
        assert "unrecognised" in caplog.text


@pytest.mark.skipif(not TESTDATA_AVAILABLE, reason=SKIP_REASON)
def test_reads_vendored_wymt_file():
    df = read_watstore(_testdata("wymt_ffa_2022A_WATSTORE.TXT"))
    powder = df[df["site_no"] == "06326500.00"]
    assert len(powder) == 85
    assert powder["water_year"].min() == 1938 and powder["water_year"].max() == 2022
    assert len(df.attrs["station_names"]) == df["site_no"].nunique()
