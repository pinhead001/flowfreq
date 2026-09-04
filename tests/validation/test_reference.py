"""Tests for flowfreq.validation.reference.

The reference side of validation used to be a subprocess wrapper around a
PeakfqSA executable that does not exist, mock-tested only. These tests exercise
the two references that do exist: the committed golden file, and -- where the
f2py extension is built -- a live ``emafitpr`` call.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from flowfreq.validation.reference import ReferenceResult, cmoms_to_parameters

GOLDEN = (
    Path(__file__).resolve().parents[1] / "fortran_parity" / "golden" / "big_sandy_03606500.json"
)


class TestCmomsToParameters:
    """cmoms(3,3): rows are mean/variance/skew, columns are the three fits."""

    CMOMS = [
        [3.7175077114155886, 3.715335627894617, 3.7169467003909995],
        [0.08470575693954105, 0.08323178876376598, 0.08433742240041339],
        [-0.15630572656992714, 0.006601394382027989, -0.11564619077395699],
    ]

    def test_column_one_is_the_regional_fit(self):
        p = cmoms_to_parameters(self.CMOMS)
        assert p["mean_log"] == pytest.approx(3.7175077, abs=1e-6)
        assert p["skew_weighted"] == pytest.approx(-0.1563057, abs=1e-6)

    def test_column_two_is_the_at_site_fit(self):
        p = cmoms_to_parameters(self.CMOMS)
        assert p["mean_log_at_site"] == pytest.approx(3.7153356, abs=1e-6)
        assert p["skew_at_site"] == pytest.approx(0.0066014, abs=1e-6)

    def test_variance_becomes_a_standard_deviation(self):
        p = cmoms_to_parameters(self.CMOMS)
        assert p["std_log"] == pytest.approx(np.sqrt(0.08470575693954105))

    def test_diagnostics_are_included_only_when_given(self):
        assert "mse_skew" not in cmoms_to_parameters(self.CMOMS)
        assert cmoms_to_parameters(self.CMOMS, skew_mse=0.094375)["mse_skew"] == pytest.approx(
            0.094375
        )

    def test_wrong_shape_is_rejected(self):
        with pytest.raises(ValueError, match="3x3"):
            cmoms_to_parameters([[1.0, 2.0], [3.0, 4.0]])


class TestFromGolden:
    def test_loads_the_committed_golden_file(self):
        ref = ReferenceResult.from_golden(GOLDEN)
        assert ref.n_peaks == 84
        assert ref.n_historical == 3
        assert ref.n_systematic == 81

    def test_source_names_the_peakfq_version(self):
        """A benchmark that cannot say what it passed against is not a benchmark."""
        assert "peakfq 8.1.0" in ReferenceResult.from_golden(GOLDEN).source

    def test_quantiles_come_back_in_cfs_not_log10(self):
        ref = ReferenceResult.from_golden(GOLDEN)
        assert ref.quantiles[0.01] == pytest.approx(22959.4, rel=1e-4)
        lower, upper = ref.confidence_intervals[0.01]
        assert lower < ref.quantiles[0.01] < upper

    def test_variance_is_carried_through(self):
        """var_est is in log10 space and is left there; it is not a discharge."""
        assert ReferenceResult.from_golden(GOLDEN).variance[0.01] > 0.0

    def test_mgbt_sentinel_is_not_read_as_a_threshold(self):
        """gbval comes back as the -99 sentinel when nothing was flagged.

        10 ** -99 is not a PILF threshold, and reporting it as one would make
        every downstream comparison against a low-outlier cut meaningless.
        """
        ref = ReferenceResult.from_golden(GOLDEN)
        assert ref.low_outlier_count == 0
        assert ref.low_outlier_threshold == 0.0

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            ReferenceResult.from_golden(tmp_path / "nope.json")

    def test_truncated_file_raises_rather_than_returning_zeros(self, tmp_path):
        path = tmp_path / "partial.json"
        path.write_text(json.dumps({"meta": {}, "inputs": {"aeps": [0.01], "dtype": [0]}}))
        with pytest.raises(KeyError):
            ReferenceResult.from_golden(path)


@pytest.mark.requires_fortran
class TestFromEmafit:
    """The f2py bridge as a reference source. Skipped when it is not built."""

    @staticmethod
    def _live() -> ReferenceResult:
        pytest.importorskip(
            "flowfreq.peakfqr",
            reason="Fortran extension not built; run python build_fortran/build.py",
        )
        from tests.fortran_parity.cases import big_sandy_case, build_emafit_inputs

        case = big_sandy_case()
        a = build_emafit_inputs(case)
        return ReferenceResult.from_emafit(
            a["ql"],
            a["qu"],
            a["tl"],
            a["tu"],
            a["dtype"],
            case.aeps,
            case.regional_skew,
            case.regional_skew_mse,
            station_name=case.description,
            eps=case.eps,
            weight_opt=case.weight_opt,
            gbthrsh0=case.gbthrsh0,
        )

    def test_live_call_matches_the_golden_file(self):
        """If these diverge, the golden file has drifted from vendor/peakfqr."""
        live = self._live()
        golden = ReferenceResult.from_golden(GOLDEN)
        assert live.parameters == golden.parameters
        assert live.quantiles == golden.quantiles
        assert live.confidence_intervals == golden.confidence_intervals

    def test_source_distinguishes_live_from_golden(self):
        assert "live emafitpr" in self._live().source
