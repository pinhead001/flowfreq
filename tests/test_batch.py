"""Tests for flowfreq.batch (multi-site batch processing).

`analyze_sites` is the one function here that talks to the network, via
`flowfreq.usgs.fetch_nwis_batch`. `flowfreq.batch` imports that name directly
(``from .usgs import fetch_nwis_batch``), so tests patch
``flowfreq.batch.fetch_nwis_batch`` rather than hitting NWIS.

Since 0.4.0, ``B17CEngine.fit`` -- and therefore every function here, all of
which sit on top of it -- uses the Bulletin 17C Eq. 7-2 unbiased station skew
rather than the pre-0.4.0 biased population coefficient (CHANGELOG.md). The
skew test below pins agreement with the ``Bulletin17C`` path rather than a
hardcoded number, so a regression to the old estimator is caught here too.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from flowfreq.batch import analyze_sites, batch_summary_table, run_multi_site
from flowfreq.bulletin17c import Bulletin17C
from flowfreq.core import PeakRecord
from flowfreq.engine import STANDARD_RETURN_PERIODS
from tests.fixtures.big_sandy import SYSTEMATIC_PEAKS


@pytest.fixture(scope="module")
def records():
    years = sorted(SYSTEMATIC_PEAKS)
    return [PeakRecord(year=y, flow=SYSTEMATIC_PEAKS[y]) for y in years]


class TestRunMultiSite:
    def test_fits_every_site(self, records):
        results = run_multi_site({"good": records})
        assert set(results["good"]) == {"params", "n", "quantiles", "ci"}
        assert results["good"]["n"] == len(records)

    def test_default_return_periods(self, records):
        results = run_multi_site({"site": records})
        assert set(results["site"]["quantiles"]) == set(STANDARD_RETURN_PERIODS)
        assert set(results["site"]["ci"]) == set(STANDARD_RETURN_PERIODS)

    def test_custom_return_periods(self, records):
        results = run_multi_site({"site": records}, return_periods=[10, 100])
        assert set(results["site"]["quantiles"]) == {10, 100}
        assert set(results["site"]["ci"]) == {10, 100}

    def test_multiple_sites_are_independent(self, records):
        results = run_multi_site({"full": records, "short": records[:15]})
        assert results["full"]["n"] == len(records)
        assert results["short"]["n"] == 15
        assert results["full"]["params"] != results["short"]["params"]

    def test_insufficient_data_reports_error_without_raising(self, records):
        """A site with fewer than 3 flows must not blow up the whole batch."""
        results = run_multi_site({"bad": records[:1]})
        assert "error" in results["bad"]
        assert "3" in results["bad"]["error"]

    def test_error_site_does_not_block_other_sites(self, records):
        results = run_multi_site({"good": records, "bad": records[:1]})
        assert "error" not in results["good"]
        assert "error" in results["bad"]

    def test_empty_input_returns_empty_dict(self):
        assert run_multi_site({}) == {}


class TestAnalyzeSites:
    def test_delegates_site_nos_and_workers_to_fetch(self, records, monkeypatch):
        captured = {}

        def fake_fetch(sites, workers=6):
            captured["sites"] = sites
            captured["workers"] = workers
            return ({"03606500": records}, {})

        monkeypatch.setattr("flowfreq.batch.fetch_nwis_batch", fake_fetch)
        analyze_sites(["03606500"], workers=3)
        assert captured["sites"] == ["03606500"]
        assert captured["workers"] == 3

    def test_returns_analysis_and_fetch_errors_as_a_tuple(self, records, monkeypatch):
        monkeypatch.setattr(
            "flowfreq.batch.fetch_nwis_batch",
            lambda sites, workers=6: ({"good": records}, {"bad_site": "404"}),
        )
        analysis, fetch_errors = analyze_sites(["good", "bad_site"])
        assert "good" in analysis
        assert set(analysis["good"]) == {"params", "n", "quantiles", "ci"}
        assert fetch_errors == {"bad_site": "404"}

    def test_fetch_failures_are_absent_from_analysis_results(self, monkeypatch):
        """A site that failed to fetch has no records to analyze at all."""
        monkeypatch.setattr(
            "flowfreq.batch.fetch_nwis_batch",
            lambda sites, workers=6: ({}, {"bad_site": "404 not found"}),
        )
        analysis, fetch_errors = analyze_sites(["bad_site"])
        assert analysis == {}
        assert fetch_errors == {"bad_site": "404 not found"}

    @pytest.mark.xfail(strict=True, reason="real bug: see docstring")
    def test_real_fetch_output_shape_is_analyzable(self):
        """Pins a genuine interface break between usgs.fetch_nwis_batch and
        batch.run_multi_site, discovered while writing these tests.

        ``usgs.fetch_nwis_peaks`` (and therefore ``fetch_nwis_batch``, which
        just fans it out) returns plain dicts --
        ``{"year": ..., "flow": ..., "source": "USGS"}`` -- per its own
        annotated return type ``List[Dict]``. But ``run_multi_site`` feeds
        those records straight to ``B17CEngine.fit``, which reads ``r.flow``
        as an *attribute*, not a key, because its real contract is
        ``List[PeakRecord]``.

        So every real call to ``analyze_sites`` fails for every site with
        ``'dict' object has no attribute 'flow'`` -- caught by
        ``run_multi_site``'s broad ``except Exception`` and reported as
        ``{"error": ...}`` rather than raising, so nothing surfaces the
        failure short of reading the output. This is not a precision or
        tolerance issue, and not fixed here: the fix means deciding whether
        ``usgs.fetch_nwis_batch`` should return ``PeakRecord`` objects or
        ``run_multi_site``/``analyze_sites`` should adapt dicts, and
        ``usgs.py`` is outside this lane.
        """
        dict_records = [
            {"year": int(y), "flow": float(f), "source": "USGS"}
            for y, f in sorted(SYSTEMATIC_PEAKS.items())
        ]
        import flowfreq.batch as batch_mod

        batch_mod.fetch_nwis_batch = lambda sites, workers=6: (
            {"03606500": dict_records},
            {},
        )
        try:
            analysis, _ = analyze_sites(["03606500"])
        finally:
            import importlib

            importlib.reload(batch_mod)
        assert "error" not in analysis["03606500"]


class TestBatchSummaryTable:
    def test_columns_include_requested_return_periods(self, records):
        results = run_multi_site({"site": records})
        df = batch_summary_table(results, return_periods=(10, 100))
        assert list(df.columns) == [
            "Site",
            "n",
            "Mean (log)",
            "Std (log)",
            "Skew",
            "Q10",
            "Q100",
        ]

    def test_default_return_periods(self, records):
        results = run_multi_site({"site": records})
        df = batch_summary_table(results)
        assert {"Q10", "Q50", "Q100"} <= set(df.columns)

    def test_one_row_per_site(self, records):
        results = run_multi_site({"a": records, "b": records[:20]})
        df = batch_summary_table(results)
        assert len(df) == 2
        assert set(df["Site"]) == {"a", "b"}

    def test_values_match_run_multi_site(self, records):
        results = run_multi_site({"site": records})
        df = batch_summary_table(results, return_periods=(100,))
        row = df.iloc[0]
        mu, sigma, skew = results["site"]["params"]
        assert row["Mean (log)"] == pytest.approx(mu)
        assert row["Std (log)"] == pytest.approx(sigma)
        assert row["Skew"] == pytest.approx(skew)
        assert row["Q100"] == pytest.approx(results["site"]["quantiles"][100])

    def test_skew_matches_bulletin17c_unbiased_estimator(self, records):
        """The Skew column comes from B17CEngine.fit, which moved from the
        biased population coefficient to Bulletin 17C's unbiased Eq. 7-2
        estimator in 0.4.0. Pin agreement with the Bulletin17C path so a
        regression to the pre-0.4.0 numbers is caught here too.
        """
        years = sorted(SYSTEMATIC_PEAKS)
        flows = [SYSTEMATIC_PEAKS[y] for y in years]
        b = Bulletin17C(peak_flows=np.array(flows), water_years=np.array(years))
        b.run_analysis(method="mom")

        results = run_multi_site({"site": records})
        df = batch_summary_table(results)
        assert df.loc[0, "Skew"] == pytest.approx(b.results.skew_station, rel=1e-6)

    def test_error_sites_get_an_error_column_instead_of_stats(self, records):
        results = run_multi_site({"bad": records[:1]})
        df = batch_summary_table(results)
        assert "Error" in df.columns
        assert "n" not in df.columns or pd.isna(df.loc[df["Site"] == "bad", "n"]).all()
        assert "3" in df.loc[df["Site"] == "bad", "Error"].iloc[0]

    def test_mixed_success_and_error_rows(self, records):
        results = run_multi_site({"good": records, "bad": records[:1]})
        df = batch_summary_table(results)
        good_row = df[df["Site"] == "good"].iloc[0]
        bad_row = df[df["Site"] == "bad"].iloc[0]
        assert pd.notna(good_row["n"])
        assert pd.isna(good_row["Error"])
        assert pd.notna(bad_row["Error"])

    def test_empty_results_returns_empty_dataframe(self):
        df = batch_summary_table({})
        assert df.empty
