"""Live verification of ``flowfreq.workflow.compare_engines``.

Design doc section 7's last correctness-plan row: "on all four parity sites,
``max_quantile_deviation_pct`` must be under the tolerances already measured
in the goldens." Runs the whole stack end to end -- native fit, live
``emafitpr`` call, ``Bulletin17C.validate()``/``FrequencyComparator`` -- on
each of the four registered parity sites.

Three of the four reproduce TODO.md P3's own measured quantile tolerances
almost exactly (Big Sandy 0.06%, Powder River 0.10%, and 12363000 0.106% --
which is also the exact number design doc section 5's own usage example
quotes: "0.106 for 12363000 (at the 500-yr)"). Cains Coulee's 9.7% and
overall FAIL are the known, already-documented ``skew_weighted`` residual
(TODO.md P3, ``tests/fortran_parity/test_wymt_vs_golden.py``'s standing
``xfail(strict=True)``) surfacing through this tool -- not a defect in the
comparison itself, and exactly the kind of real disagreement section 1 of
the design doc says this feature exists to surface.
"""

from __future__ import annotations

import pytest

pytest.importorskip(
    "flowfreq.peakfqr",
    reason="Fortran extension not built; run python build_fortran/build.py "
    "(see docs/FORTRAN_UPLOAD.md)",
    # See tests/fortran_parity/test_fortran_oracles.py's exc_type comment: pytest
    # 9.1 narrowed importorskip's default exc_type to ModuleNotFoundError, but this
    # module deliberately re-raises as ImportError.
    exc_type=ImportError,
)

pytestmark = pytest.mark.requires_fortran


def _big_sandy_report():
    from flowfreq.workflow import compare_engines
    from tests.fixtures.big_sandy import (
        HISTORICAL_PEAKS,
        REGIONAL_SKEW,
        REGIONAL_SKEW_SD,
        SYSTEMATIC_PEAKS,
        THRESHOLDS,
    )

    return compare_engines(
        peak_flows=list(SYSTEMATIC_PEAKS.values()),
        water_years=list(SYSTEMATIC_PEAKS.keys()),
        regional_skew=REGIONAL_SKEW,
        regional_skew_se=REGIONAL_SKEW_SD,
        historical_peaks=list(HISTORICAL_PEAKS.items()),
        perception_thresholds={(t["start"], t["end"]): t["lower"] for t in THRESHOLDS},
        site_name="Big Sandy River at Bruceton, TN (USGS 03606500)",
    )


def _wymt_report(site_no, label):
    from flowfreq.workflow import compare_engines
    from tests.fixtures.wymt_peaks import load_site

    site = load_site(site_no)
    years = sorted(site.peaks)
    return compare_engines(
        peak_flows=[site.peaks[y] for y in years],
        water_years=years,
        regional_skew=site.regional_skew,
        regional_skew_se=site.regional_skew_mse**0.5,
        site_name=label,
    )


def _site_12363000_report():
    from flowfreq.workflow import compare_engines
    from tests.fixtures.site_12363000 import REGIONAL_SKEW, REGIONAL_SKEW_SD, SYSTEMATIC_PEAKS

    years = sorted(SYSTEMATIC_PEAKS)
    return compare_engines(
        peak_flows=[SYSTEMATIC_PEAKS[y] for y in years],
        water_years=years,
        regional_skew=REGIONAL_SKEW,
        regional_skew_se=REGIONAL_SKEW_SD,
        site_name="USGS 12363000",
    )


class TestBigSandy:
    def test_passes_and_stays_under_the_measured_tolerance(self):
        report = _big_sandy_report()
        assert report.comparison.passed
        # TODO.md P3: "every quantile to <= 0.06%". A little headroom above
        # the measured value keeps this from flapping on a different
        # gfortran build's last-ulp noise, while still catching a real
        # regression an order of magnitude off.
        assert report.max_quantile_deviation_pct < 0.15

    def test_markdown_renders(self):
        report = _big_sandy_report()
        md = report.to_markdown()
        assert md.startswith("# Engine comparison: Big Sandy")
        assert "## Quantiles (cfs)" in md
        assert "## Confidence intervals (cfs)" in md


class TestPowderRiver:
    def test_passes_and_stays_under_the_measured_tolerance(self):
        try:
            report = _wymt_report("06326500.00", "Powder River")
        except FileNotFoundError as exc:
            pytest.skip(str(exc))
        assert report.comparison.passed
        # TODO.md P3: "quantiles <= 0.10%".
        assert report.max_quantile_deviation_pct < 0.2


class TestSite12363000:
    def test_passes_and_stays_under_the_measured_tolerance(self):
        report = _site_12363000_report()
        assert report.comparison.passed
        # design doc section 5's own worked example: "0.106 for 12363000
        # (at the 500-yr)".
        assert report.max_quantile_deviation_pct < 0.2


class TestCainsCoulee:
    """MGBT-censored site; the engines agree now.

    Its ``skew_weighted`` was 0.058 off peakfq (a standing ``xfail(strict=True)``)
    until the native at-site skew MSE followed ``emafit.f:707``'s switch to the
    Bulletin 17B formula when MGBT finds low outliers. Both tests here were the
    comparison-level view of that gap and flipped with it.
    """

    def test_quantile_deviation_stays_under_the_measured_bound(self):
        """Measured about 0.001% (was up to 9.7% before the B17B skew-MSE switch)."""
        try:
            report = _wymt_report("06327450.00", "Cains Coulee")
        except FileNotFoundError as exc:
            pytest.skip(str(exc))
        assert report.max_quantile_deviation_pct < 0.1

    def test_overall_comparison_passes(self):
        try:
            report = _wymt_report("06327450.00", "Cains Coulee")
        except FileNotFoundError as exc:
            pytest.skip(str(exc))
        assert report.comparison.passed


class TestCodedRecord:
    """A real coded record (#30): both engines get the same siteQT rows.

    USGS 01426500 from the vendored HU02 WATSTORE file: 50 code 6 (regulated)
    peaks that peakfq removes and one code 7 historic peak. ``peak_code_kwargs``
    checks each engine rebuilds ``siteQT``'s rows before either is run, so the
    comparison is like for like.
    """

    def test_coded_record_passes(self):
        from flowfreq.watstore import read_watstore
        from flowfreq.workflow import compare_engines
        from tests.fixtures.paths import TESTDATA_AVAILABLE, testdata_path

        if not TESTDATA_AVAILABLE:
            pytest.skip("peakfqr reference test data not present")
        frame = read_watstore(testdata_path("extra_tests/HU02_WATSTORE.txt"))
        site = frame[frame["site_no"].astype(str).str.strip() == "01426500"]
        report = compare_engines(
            peak_flows=site["peak_flow_cfs"].to_numpy(),
            water_years=site["water_year"].to_numpy(),
            station_skew_only=True,
            peak_codes=site["qualification_code"].tolist(),
            site_name="USGS 01426500",
        )
        assert report.native.n_peaks == 52
        assert report.comparison.passed
        assert report.max_quantile_deviation_pct < 0.2
