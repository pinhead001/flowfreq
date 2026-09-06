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
    """The one site with a known, already-documented residual.

    Its ``skew_weighted`` gap (0.058 skew units) is the standing
    ``xfail(strict=True)`` in
    ``tests/fortran_parity/test_wymt_vs_golden.py::TestCainsCouleeCensored::
    test_weighted_skew_matches``. ``compare_engines`` surfacing the same
    known gap as an overall FAIL is correct behaviour, not a bug in the
    comparison -- so that FAIL is asserted here too, `xfail(strict=True)`
    the same way, rather than skipped or worked around. If that residual is
    ever resolved, this xfail should flip at the same time as the one in
    ``test_wymt_vs_golden.py``.
    """

    def test_quantile_deviation_stays_under_the_measured_bound(self):
        """Not xfailed: TODO.md P3 already measured this at up to 9.7%, and
        that bound is not itself in question -- only whether the overall
        comparison *passes* is."""
        try:
            report = _wymt_report("06327450.00", "Cains Coulee")
        except FileNotFoundError as exc:
            pytest.skip(str(exc))
        assert report.max_quantile_deviation_pct < 15.0

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "Known, documented residual (TODO.md P3, "
            "test_wymt_vs_golden.py::test_weighted_skew_matches): skew_weighted is "
            "0.058 skew units off peakfq 8.1.0's own reported value, tracing to "
            "emafitpr's own internally-computed as_G_mse disagreeing with a standalone "
            "call to the same mseg_all routine on identical inputs -- not a defect in "
            "flowfreq's composition of independently-verified routines. Remove this "
            "xfail if that residual is ever resolved."
        ),
    )
    def test_overall_comparison_passes(self):
        try:
            report = _wymt_report("06327450.00", "Cains Coulee")
        except FileNotFoundError as exc:
            pytest.skip(str(exc))
        assert report.comparison.passed
