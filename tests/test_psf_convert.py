"""Tests for flowfreq.psf_convert: .psf station -> siteQT rows -> Bulletin17C.

The rules under test are peakfq 8.1.0's (``vendor/peakfqr/R/readInputs.R::siteQT``
and ``vendor/peakfqr/R/main.R``). Numeric parity against peakfq's own output is
in ``tests/fortran_parity/test_psf_vs_golden.py`` and
``tests/fortran_parity/test_live_psf_convert.py``.
"""

from __future__ import annotations

import csv
import logging

import numpy as np
import pandas as pd
import pytest

from flowfreq.bulletin17c import Bulletin17C
from flowfreq.peak_codes import Q_MAX, Q_MIN
from flowfreq.psf import parse_psf
from flowfreq.psf_convert import (
    PEAKFQ_AEPS,
    PEAKFQ_AEPS_EXTENDED,
    STATION_SKEW_MSE_SENTINEL,
    StationInputs,
    UnsupportedSpecError,
    convert_peak_frame,
    convert_psf,
    convert_station,
    read_psf_peaks,
    station_rows,
)
from tests.fixtures.paths import SKIP_REASON, TESTDATA_AVAILABLE
from tests.fixtures.paths import testdata_path as _testdata

needs_testdata = pytest.mark.skipif(not TESTDATA_AVAILABLE, reason=SKIP_REASON)

YEARS = list(range(2001, 2016))
FLOWS = [float(100 + 37 * i % 211) for i in range(len(YEARS))]


def _spec(body: str, header: str = ""):
    psf = parse_psf(f"{header}Station S1\n{body}")
    return psf.stations["S1"], psf


def _peaks(codes=None, years=YEARS, flows=FLOWS) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "water_year": years,
            "peak_flow_cfs": flows,
            "qualification_code": codes or [""] * len(years),
        }
    )


BASE = "PCPT_Thresh 2001 2015 0 1E+20\nSkewOpt Station\nLOType MGBT\n"


def _row(rows, year):
    return next(r for r in rows if r.year == year)


class TestStationRows:
    def test_plain_systematic(self):
        spec, _ = _spec(BASE)
        rows = station_rows(spec, _peaks())
        assert [r.year for r in rows] == YEARS
        assert all(r.is_exact and r.tl == Q_MIN and r.tu == Q_MAX for r in rows)
        assert all(r.dtype == 0 for r in rows)

    def test_codes_follow_sitqt(self):
        codes = [""] * len(YEARS)
        codes[1], codes[2], codes[3], codes[4] = "4", "8", "3", "6"
        codes[5] = "7"
        spec, _ = _spec(BASE)
        rows = station_rows(spec, _peaks(codes))
        assert (_row(rows, 2002).ql, _row(rows, 2002).qu) == (Q_MIN, FLOWS[1])
        assert (_row(rows, 2003).ql, _row(rows, 2003).qu) == (FLOWS[2], Q_MAX)
        # Removed peaks carry no information, so the year is dropped.
        assert 2004 not in {r.year for r in rows}
        assert 2005 not in {r.year for r in rows}
        assert _row(rows, 2006).dtype == 1

    def test_urban_regulated_kept_by_spec(self):
        codes = [""] * len(YEARS)
        codes[4] = "C"
        spec, _ = _spec(BASE + "Urb/Reg Yes\n")
        rows = station_rows(spec, _peaks(codes))
        assert _row(rows, 2005).is_exact

    def test_gap_year_censored_at_threshold_or_dropped(self, caplog):
        years = [y for y in YEARS if y not in (2004, 2010)]
        flows = [f for y, f in zip(YEARS, FLOWS) if y not in (2004, 2010)]
        spec, _ = _spec(BASE + "PCPT_Thresh 2004 2004 500 1E+20\n")
        rows = station_rows(spec, _peaks(years=years, flows=flows))
        gap = _row(rows, 2004)
        assert (gap.ql, gap.qu, gap.tl, gap.tu) == (Q_MIN, 500.0, 500.0, Q_MAX)
        # 2010 sits under the zero default threshold: missing data, dropped.
        assert 2010 not in {r.year for r in rows}
        assert "zero perception threshold" in caplog.text

    def test_later_threshold_wins(self):
        spec, _ = _spec(BASE + "PCPT_Thresh 2003 2004 50 1E+20\nPCPT_Thresh 2004 2004 70 1E+20\n")
        rows = station_rows(spec, _peaks())
        assert (_row(rows, 2003).tl, _row(rows, 2004).tl) == (50.0, 70.0)

    def test_interval_and_peak_lines_override_data(self, caplog):
        spec, _ = _spec(BASE + "Peak 2002 999\nInterval 2002 0 12\nPeak 2003 5\n")
        rows = station_rows(spec, _peaks())
        assert (_row(rows, 2002).ql, _row(rows, 2002).qu) == (Q_MIN, 12.0)
        assert _row(rows, 2003).ql == _row(rows, 2003).qu == 5.0
        assert "overwritten by Interval" in caplog.text

    def test_missing_threshold_period_drops_its_years(self):
        spec, _ = _spec(BASE + "PCPT_Thresh 2003 2004 1E+20 1E+20 MISSING\n")
        rows = station_rows(spec, _peaks())
        assert {2003, 2004}.isdisjoint(r.year for r in rows)

    def test_zero_flow_is_exact_at_qmin(self):
        flows = list(FLOWS)
        flows[0] = 0.0
        spec, _ = _spec(BASE)
        row = _row(station_rows(spec, _peaks(flows=flows)), 2001)
        assert row.ql == row.qu == Q_MIN and row.flow == 0.0

    def test_peaks_outside_period_dropped(self):
        spec, _ = _spec("PCPT_Thresh 2005 2010 0 1E+20\nSkewOpt Station\nLOType MGBT\n")
        assert [r.year for r in station_rows(spec, _peaks())] == list(range(2005, 2011))

    def test_filters_multi_site_frame(self):
        spec, _ = _spec(BASE)
        both = pd.concat([_peaks().assign(site_no="S1"), _peaks().assign(site_no="S2")])
        assert len(station_rows(spec, both)) == len(YEARS)

    def test_keep_no_info(self):
        codes = [""] * len(YEARS)
        codes[0] = "3"
        spec, _ = _spec(BASE)
        rows = station_rows(spec, _peaks(codes), keep_no_info=True)
        assert _row(rows, 2001).tl == Q_MAX

    @pytest.mark.parametrize(
        "body, match",
        [
            ("SkewOpt Station\n", "no PCPT_Thresh"),
            (BASE + "Interval 1990 0 5\n", "not in analysis range"),
            (BASE + "Peak 2003 -8888\n", "Interval YYYY 0 1E"),
            (BASE + "Peak 2003 -5\n", "Negative discharge"),
            (
                "PCPT_Thresh 2001 2005 0 1E+20\nPCPT_Thresh 2010 2015 0 1E+20\n",
                "missing perception",
            ),
        ],
    )
    def test_errors(self, body, match):
        spec, _ = _spec(body)
        with pytest.raises(ValueError, match=match):
            station_rows(spec, _peaks())

    def test_duplicate_year_raises(self):
        spec, _ = _spec(BASE)
        with pytest.raises(ValueError, match="Duplicated water year"):
            station_rows(spec, _peaks(years=[2001, 2001, 2002], flows=[1.0, 2.0, 3.0]))

    def test_no_peaks_raises(self):
        spec, _ = _spec(BASE)
        with pytest.raises(ValueError, match="No peak flow data"):
            station_rows(spec, _peaks().iloc[0:0])


class TestSettings:
    @pytest.mark.parametrize(
        "skew, expected",
        [
            ("SkewOpt Station\n", ("station", None, None, 0.0, STATION_SKEW_MSE_SENTINEL)),
            (
                "SkewOpt Weighted\nGenSkew -0.3\nSkewSE 0.5\n",
                ("weighted", -0.3, 0.25, -0.3, 0.25),
            ),
            (
                "SkewOpt Generalized\nGenSkew -0.3\nSkewSE 0.5\n",
                ("generalized", -0.3, -0.25, -0.3, -0.25),
            ),
            (
                "SkewOpt regional\nGenSkew -0.3\nSkewSE 0.5\n",
                ("generalized", -0.3, -0.25, -0.3, -0.25),
            ),
        ],
    )
    def test_skew_encoding(self, skew, expected):
        spec, _ = _spec("PCPT_Thresh 2001 2015 0 1E+20\nLOType MGBT\n" + skew)
        s = convert_station(spec, _peaks())
        got = (s.skew_option, s.regional_skew, s.regional_skew_mse, s._r_g, s._r_g_mse)
        assert got == pytest.approx(expected) if None not in expected else got == expected

    @pytest.mark.parametrize(
        "skew, match",
        [
            ("", "no SkewOpt"),
            ("SkewOpt Sideways\n", "Invalid skew option"),
            ("SkewOpt Weighted\nSkewSE 0.5\n", "regional skew coefficient"),
            ("SkewOpt Weighted\nGenSkew -0.3\n", "standard error"),
            ("SkewOpt Weighted\nGenSkew -0.3\nSkewSE 0\n", "not positive"),
        ],
    )
    def test_skew_errors(self, skew, match):
        spec, _ = _spec("PCPT_Thresh 2001 2015 0 1E+20\nLOType MGBT\n" + skew)
        with pytest.raises(ValueError, match=match):
            convert_station(spec, _peaks())

    @pytest.mark.parametrize(
        "lo, method, threshold",
        [
            ("LOType MGBT\n", "MGBT", None),
            ("LOType FIXED\nLoThresh 150\n", "FIXED", 150.0),
            ("LOType NONE\n", "NONE", Q_MIN),
        ],
    )
    def test_low_outlier(self, lo, method, threshold):
        spec, _ = _spec("PCPT_Thresh 2001 2015 0 1E+20\nSkewOpt Station\n" + lo)
        s = convert_station(spec, _peaks())
        assert (s.low_outlier_method, s.user_low_outlier_threshold) == (method, threshold)
        # emafitpr's gbthrsh0 encoding: <= -6 runs MGBT, otherwise a log10 threshold.
        gb = s.emafit_arrays().gbthrsh0
        if method == "MGBT":
            assert gb == -99.0
        else:
            assert gb == pytest.approx(np.log10(threshold))

    @pytest.mark.parametrize(
        "lo, match",
        [("", "no LOType"), ("LOType FIXED\n", "low outlier threshold"), ("LOType X\n", "option")],
    )
    def test_low_outlier_errors(self, lo, match):
        spec, _ = _spec("PCPT_Thresh 2001 2015 0 1E+20\nSkewOpt Station\n" + lo)
        with pytest.raises(ValueError, match=match):
            convert_station(spec, _peaks())

    def test_weight_option(self):
        spec, _ = _spec(BASE + "WeightOpt erl\n")
        assert convert_station(spec, _peaks()).weight_option == "ERL"
        bad, _ = _spec(BASE + "WeightOpt XYZ\n")
        with pytest.raises(ValueError, match="weighting option"):
            convert_station(bad, _peaks())

    def test_output_options(self, caplog):
        header = (
            "O ConfInterval 0.95\nO Extended YES\nO EMA YES\nO File out.PRT\n"
            "O Plot Position 0.3\nO Mystery 1\n"
        )
        spec, psf = _spec(BASE, header)
        with caplog.at_level(logging.INFO):
            s = convert_station(spec, _peaks(), psf)
        assert s.confidence == 0.95
        assert s.aeps == PEAKFQ_AEPS_EXTENDED
        assert s.method == "ema"
        assert "O File out.PRT" in caplog.text and "no flowfreq equivalent" in caplog.text
        assert "Plot Position" in caplog.text
        assert "unknown output option" in caplog.text
        assert "only reads the spelling" in caplog.text  # "Extended", not "EXTENDED"

    def test_defaults_without_psf(self):
        spec, _ = _spec(BASE)
        s = convert_station(spec, _peaks())
        assert (s.confidence, s.aeps, s.method) == (0.90, PEAKFQ_AEPS, "ema")

    def test_invalid_confidence_defaults(self, caplog):
        spec, psf = _spec(BASE, "O ConfInterval 1.5\n")
        assert convert_station(spec, _peaks(), psf).confidence == 0.90
        assert "Invalid confidence" in caplog.text

    def test_ema_no_is_mom(self, caplog):
        spec, psf = _spec(BASE, "O EMA NO\n")
        s = convert_station(spec, _peaks(), psf)
        assert s.method == "mom"
        assert "not PeakFQ 7's" in caplog.text
        assert s.unsupported_reasons("native") == []
        assert s.unsupported_reasons("fortran")

    def test_ema_invalid_raises(self):
        spec, psf = _spec(BASE, "O EMA MAYBE\n")
        with pytest.raises(ValueError, match="O EMA"):
            convert_station(spec, _peaks(), psf)

    def test_short_record_warns(self, caplog):
        spec, _ = _spec("PCPT_Thresh 2001 2005 0 1E+20\nSkewOpt Station\nLOType MGBT\n")
        convert_station(spec, _peaks())
        assert "fewer than 10 rows" in caplog.text


class TestStationInputs:
    def test_kwargs_build_bulletin17c(self):
        spec, _ = _spec(BASE)
        s = convert_station(spec, _peaks())
        kw = s.bulletin17c_kwargs("native")
        assert list(kw["water_years"]) == YEARS
        assert kw["regional_skew"] is None and kw["historical_peaks"] is None
        assert isinstance(s.to_bulletin17c(), Bulletin17C)

    def test_run_reports_at_psf_aeps(self):
        spec, psf = _spec(BASE, "O ConfInterval 0.95\n")
        results = convert_station(spec, _peaks(), psf).run("native")
        assert np.allclose(results.quantiles["aep"], PEAKFQ_AEPS)
        assert len(results.confidence_limits) == len(PEAKFQ_AEPS)

    def test_historic_peak_and_period(self):
        # A pre-systematic historical period with one code-7 flood: the case
        # both engines express, so both round trips must reproduce siteQT.
        years = [1990] + YEARS
        flows = [5000.0] + FLOWS
        codes = ["7"] + [""] * len(YEARS)
        spec, _ = _spec(
            "PCPT_Thresh 1985 2015 0 1E+20\nPCPT_Thresh 1985 2000 3000 1E+20\n"
            "SkewOpt Station\nLOType MGBT\n"
        )
        s = convert_station(spec, _peaks(codes, years, flows))
        assert s.historical_peaks == ((1990, 5000.0),)
        assert s.n_historic == 1
        assert s.unsupported_reasons("native") == []
        assert s.unsupported_reasons("fortran") == []

    def test_interval_is_unsupported(self):
        codes = [""] * len(YEARS)
        codes[3] = "8"
        spec, _ = _spec(BASE)
        s = convert_station(spec, _peaks(codes))
        for engine in ("native", "fortran"):
            with pytest.raises(UnsupportedSpecError, match="greater-than"):
                s.bulletin17c_kwargs(engine)

    def test_zero_flow_unsupported_natively(self):
        flows = list(FLOWS)
        flows[2] = 0.0
        spec, _ = _spec(BASE)
        s = convert_station(spec, _peaks(flows=flows))
        assert any("zero flows" in r for r in s.unsupported_reasons("native"))
        assert s.unsupported_reasons("fortran") == []

    def test_generalized_skew_unsupported_natively(self):
        spec, _ = _spec(
            "PCPT_Thresh 2001 2015 0 1E+20\nLOType MGBT\nSkewOpt Generalized\n"
            "GenSkew -0.2\nSkewSE 0.5\n"
        )
        s = convert_station(spec, _peaks())
        assert any("Generalized" in r for r in s.unsupported_reasons("native"))
        assert s.unsupported_reasons("fortran") == []

    def test_non_hwn_weighting_unsupported(self):
        spec, _ = _spec(
            "PCPT_Thresh 2001 2015 0 1E+20\nLOType MGBT\nSkewOpt Weighted\n"
            "GenSkew -0.2\nSkewSE 0.5\nWeightOpt INV\n"
        )
        s = convert_station(spec, _peaks())
        assert all(
            any("WeightOpt" in r for r in s.unsupported_reasons(e)) for e in ("native", "fortran")
        )

    def test_unknown_engine_raises(self):
        spec, _ = _spec(BASE)
        with pytest.raises(ValueError, match="engine"):
            convert_station(spec, _peaks()).unsupported_reasons("gpu")

    def test_emafit_arrays_are_log10_rows(self):
        spec, _ = _spec(BASE)
        s = convert_station(spec, _peaks())
        arrays = s.emafit_arrays()
        assert np.allclose(10**arrays.ql, FLOWS)
        assert arrays.n == len(YEARS) and arrays.n_censored == 0
        assert isinstance(s, StationInputs)


class TestConvertPeakFrame:
    def test_no_codes_matches_plain_bulletin17c(self):
        """Uncoded contiguous data: the same arguments, so the same fit."""
        s = convert_peak_frame(_peaks(), regional_skew=-0.1, regional_skew_mse=0.3)
        assert (s.skew_option, s.regional_skew, s.regional_skew_mse) == ("weighted", -0.1, 0.3)
        assert s.low_outlier_method == "MGBT"
        kw = s.bulletin17c_kwargs("native")
        assert list(kw["peak_flows"]) == FLOWS and list(kw["water_years"]) == YEARS
        via = Bulletin17C(**kw).run_analysis()
        plain = Bulletin17C(
            np.array(FLOWS), water_years=np.array(YEARS), regional_skew=-0.1, regional_skew_mse=0.3
        ).run_analysis()
        assert (via.mean_log, via.std_log, via.skew_used) == (
            plain.mean_log,
            plain.std_log,
            plain.skew_used,
        )

    def test_codes_applied(self):
        codes = [""] * len(YEARS)
        codes[0], codes[4] = "7", "6,Bd"
        s = convert_peak_frame(_peaks(codes), low_outlier_threshold=10.0)
        assert s.historical_peaks == ((2001, FLOWS[0]),)
        assert 2005 not in s.water_years
        assert (s.low_outlier_method, s.user_low_outlier_threshold) == ("FIXED", 10.0)
        kept = convert_peak_frame(_peaks(codes), include_urban_regulated=True)
        assert 2005 in kept.water_years

    def test_generalized_and_station(self):
        gen = convert_peak_frame(_peaks(), regional_skew=-0.1, regional_skew_mse=-0.3)
        assert (gen.skew_option, gen.regional_skew_mse) == ("generalized", -0.3)
        assert convert_peak_frame(_peaks()).skew_option == "station"

    def test_zero_mse_raises(self):
        with pytest.raises(ValueError, match="positive skew standard error"):
            convert_peak_frame(_peaks(), regional_skew=-0.1, regional_skew_mse=0.0)

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="no water years"):
            convert_peak_frame(_peaks().iloc[0:0])


class TestAnalyzeGageWiring:
    """``flowfreq.analyze_gage``'s peak-frame -> Bulletin17C step."""

    def test_default_ignores_codes_but_says_so(self, caplog):
        from flowfreq import _gage_analysis

        codes = [""] * len(YEARS)
        codes[2] = "6"
        b17c = _gage_analysis(_peaks(codes), "S1", method="ema")
        assert b17c.n == len(YEARS)
        assert "code 6: 1" in caplog.text

    def test_default_silent_without_acted_on_codes(self, caplog):
        from flowfreq import _gage_analysis

        _gage_analysis(_peaks(["2"] * len(YEARS)), "S1", method="ema")
        assert "qualification codes" not in caplog.text

    def test_apply_peak_codes(self):
        from flowfreq import _gage_analysis

        codes = [""] * len(YEARS)
        codes[-1] = "6"
        b17c = _gage_analysis(_peaks(codes), "S1", method="ema", apply_peak_codes=True)
        assert b17c.n == len(YEARS) - 1

    def test_mid_record_removal_is_missing_data(self):
        """siteQT drops a removed year; the native engine, given the gap, must
        build the same rows -- which the round trip checks rather than assumes."""
        from flowfreq import _gage_analysis

        codes = [""] * len(YEARS)
        codes[2] = "6"
        b17c = _gage_analysis(_peaks(codes), "S1", method="ema", apply_peak_codes=True)
        b17c.run_analysis()
        assert b17c.results.n_peaks == len(YEARS) - 1

    def test_apply_peak_codes_refuses_censored_peaks(self):
        from flowfreq import _gage_analysis

        codes = [""] * len(YEARS)
        codes[2] = "4"
        with pytest.raises(UnsupportedSpecError, match="2003"):
            _gage_analysis(_peaks(codes), "S1", method="ema", apply_peak_codes=True)

    @pytest.mark.parametrize(
        "kwargs, match",
        [
            ({"method": "mom"}, "method='ema'"),
            ({"method": "ema", "historical_peaks": [(1990, 5.0)]}, "historical_peaks"),
        ],
    )
    def test_apply_peak_codes_errors(self, kwargs, match):
        from flowfreq import _gage_analysis

        with pytest.raises(ValueError, match=match):
            _gage_analysis(_peaks(), "S1", apply_peak_codes=True, **kwargs)


class TestReadPsfPeaks:
    def test_no_peak_file_raises(self):
        _, psf = _spec(BASE)
        with pytest.raises(ValueError, match="names no peak file"):
            read_psf_peaks(psf, "x.psf")

    def test_rdb_not_implemented(self):
        _, psf = _spec(BASE, "I RDB peaks.txt\n")
        with pytest.raises(NotImplementedError, match="convert_station"):
            read_psf_peaks(psf, "x.psf")

    def test_missing_file_raises(self, tmp_path):
        _, psf = _spec(BASE, "I ASCI nothing_here.txt\n")
        with pytest.raises(FileNotFoundError):
            read_psf_peaks(psf, tmp_path / "x.psf")


@needs_testdata
class TestWymt:
    """Every station in peakfq's own WY/MT test file converts, and agrees with peakfq."""

    @pytest.fixture(scope="class")
    def stations(self):
        return convert_psf(_testdata("wymt_ffa_2022A.psf"))

    @pytest.fixture(scope="class")
    def info(self):
        with open(_testdata("wymt_ffa_2022A_EXPinfo_7_4.csv"), newline="") as fh:
            return {r["site_no"]: r for r in csv.DictReader(fh)}

    def test_all_stations_convert(self, stations, info):
        assert set(stations) == set(info)

    def test_site_info_matches_peakfq(self, stations, info):
        """peakfq's own ``Site Info`` test (``tests/testthat/test-moments.R``)."""
        for sid, s in stations.items():
            row = info[sid]
            got = (s.begin_year, s.end_year, s.n_historic, s.skew_option)
            want = (
                int(row["BegYear"]),
                int(row["EndYear"]),
                int(row["HistPeaks"]),
                row["SkewOption"].lower(),
            )
            assert got == want, sid

    def test_gaged_peak_count_matches_peakfq(self, stations, info):
        """``GagedPeaks`` counts every non-historic row that holds a recorded peak.

        Censored ones included (code 4 or an ``Interval``); gap years filled
        against a perception threshold excluded.
        """
        from flowfreq.watstore import read_watstore

        peaks = read_watstore(_testdata("wymt_ffa_2022A_WATSTORE.TXT"))
        for sid, s in stations.items():
            recorded = set(peaks.loc[peaks["site_no"] == sid, "water_year"])
            n = sum(1 for r in s.rows if r.dtype == 0 and r.year in recorded)
            assert n == int(info[sid]["GagedPeaks"]), sid

    def test_every_station_runs_through_some_route(self, stations):
        """Only 06328100 (upper threshold 407 cfs in 2021) needs fortran_reference()."""
        needs_direct = sorted(
            sid for sid, s in stations.items() if s.unsupported_reasons("fortran")
        )
        assert needs_direct == ["06328100.00"]

    def test_skew_settings(self, stations):
        powder = stations["06326500.00"]
        assert powder.regional_skew == -0.2002
        assert powder.regional_skew_mse == pytest.approx(0.64**2)
        assert powder.low_outlier_method == "MGBT"
        assert powder.confidence == 0.9
        station_skew = stations["06185500.11"]
        assert station_skew.regional_skew is None and station_skew.regional_skew_mse is None
