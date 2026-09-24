"""Tests for flowfreq.psf against the vendored WY/MT PeakFQ specification file."""

from __future__ import annotations

import csv

import pytest

from flowfreq.psf import parse_psf, read_psf
from tests.fixtures.paths import SKIP_REASON, TESTDATA_AVAILABLE
from tests.fixtures.paths import testdata_path as _testdata

needs_testdata = pytest.mark.skipif(not TESTDATA_AVAILABLE, reason=SKIP_REASON)

SAMPLE = """\
I ASCI peaks.txt
O ConfInterval 0.9
' a comment line
Station 06177820.00
     PCPT_Thresh 1974 2019 0 1E+20  DEFAULT
     PCPT_Thresh 2015 2015 12 1E+20  PEAK < STATED VALUE
     Interval 2015 0 12  PEAK < STATED VALUE
     SkewOpt Weighted
     GenSkew -0.397702
     SkewSE 0.64
     LOType MGBT
     Urb/Reg Yes
Station 06185500.10
     PCPT_Thresh 1942 2022 0 1E+20
     SkewOpt Station
"""


class TestParse:
    def test_sample(self):
        psf = parse_psf(SAMPLE)
        assert psf.confidence_interval == 0.9
        assert list(psf.stations) == ["06177820.00", "06185500.10"]
        st = psf.stations["06177820.00"]
        assert st.site_no == "06177820"
        assert st.skew_option == "Weighted"
        assert st.generalized_skew == pytest.approx(-0.397702)
        assert st.generalized_skew_mse == pytest.approx(0.4096)
        assert st.include_urban_regulated
        assert st.analysis_years == (1974, 2019)
        assert st.thresholds[1].comment == "PEAK < STATED VALUE"
        assert (st.intervals[0].year, st.intervals[0].lower, st.intervals[0].upper) == (
            2015,
            0.0,
            12.0,
        )
        other = psf.stations["06185500.10"]
        assert other.generalized_skew is None and not other.include_urban_regulated
        assert other.thresholds[0].comment == ""

    def test_field_before_station_raises(self):
        with pytest.raises(ValueError, match="before any Station"):
            parse_psf("SkewOpt Weighted\n")

    def test_duplicate_station_raises(self):
        with pytest.raises(ValueError, match="duplicate"):
            parse_psf("Station 1\nStation 1\n")

    def test_malformed_threshold_raises(self):
        with pytest.raises(ValueError, match="malformed"):
            parse_psf("Station 1\n PCPT_Thresh 1990 x 0 1E+20\n")

    def test_unknown_and_deprecated_fields_warn(self, caplog):
        parse_psf("Station 1\n BegYear 1990\n Bogus 3\n Peak 1990 100\n")
        assert "Deprecated PSF field" in caplog.text
        assert "Unknown PSF field" in caplog.text
        assert "Deprecated PSF keyword 'Peak'" in caplog.text


@needs_testdata
class TestWymtFile:
    @pytest.fixture(scope="class")
    def psf(self):
        return read_psf(_testdata("wymt_ffa_2022A.psf"))

    def test_counts(self, psf):
        assert len(psf.stations) == 24
        assert sum(len(s.thresholds) for s in psf.stations.values()) == 53
        assert sum(len(s.intervals) for s in psf.stations.values()) == 10

    def test_regional_skew_matches_peakfq_output(self, psf):
        """GenSkew/SkewSE in the spec match RegSkew/RegMSEG in peakfq's own EXPinfo output."""
        with open(_testdata("wymt_ffa_2022A_EXPinfo_7_4.csv"), newline="") as fh:
            info = {r["site_no"]: r for r in csv.DictReader(fh)}
        checked = 0
        for sid, st in psf.stations.items():
            if st.skew_option != "Weighted" or sid not in info:
                continue
            assert st.generalized_skew == pytest.approx(float(info[sid]["RegSkew"]), abs=5e-4)
            assert st.generalized_skew_mse == pytest.approx(float(info[sid]["RegMSEG"]), abs=5e-3)
            checked += 1
        assert checked > 0
