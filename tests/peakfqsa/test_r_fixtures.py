"""
Tests that validate fixture data extracted from peakfqr R tests.

These tests verify the fixture data is loadable and internally consistent.
Full numerical validation against Fortran output requires the f2py bridge
(Track A) and is deferred to integration tests.
"""

from __future__ import annotations

import math
import os

import pytest


class TestFortranRespecFixtures:
    """Validate fortran_respec fixture data is well-formed."""

    def test_ql_qu_same_length(self) -> None:
        from tests.peakfqsa.fixtures.fortran_respec import QL, QU

        assert len(QL) == len(QU) == 34

    def test_truth_data_length(self) -> None:
        from tests.peakfqsa.fixtures.fortran_respec import QL_TRUTH, QU_TRUTH

        assert len(QL_TRUTH) == len(QU_TRUTH) == 200

    def test_moms_p3_expected_has_three_variants(self) -> None:
        from tests.peakfqsa.fixtures.fortran_respec import MOMS_P3_EXPECTED

        assert set(MOMS_P3_EXPECTED.keys()) == {"orig", "ERL", "HWN"}
        for key, vals in MOMS_P3_EXPECTED.items():
            assert len(vals) == 3, f"{key} should have [mean, var, skew]"

    def test_p3est_ema_expected_has_three_variants(self) -> None:
        from tests.peakfqsa.fixtures.fortran_respec import P3EST_EMA_EXPECTED

        assert set(P3EST_EMA_EXPECTED.keys()) == {"orig", "ERL", "HWN"}
        for key, vals in P3EST_EMA_EXPECTED.items():
            assert len(vals) == 3

    def test_truth_moms_expected_values(self) -> None:
        from tests.peakfqsa.fixtures.fortran_respec import TRUTH_MOMS_P3_EXPECTED

        assert len(TRUTH_MOMS_P3_EXPECTED) == 3
        # mean, variance, skew should all be finite
        assert all(math.isfinite(v) for v in TRUTH_MOMS_P3_EXPECTED)

    def test_qp3_exceedance_probs(self) -> None:
        from tests.peakfqsa.fixtures.fortran_respec import QP3_EXCEEDANCE_PROBS

        assert len(QP3_EXCEEDANCE_PROBS) == 7
        # Should be in decreasing order (high to low exceedance)
        assert QP3_EXCEEDANCE_PROBS == sorted(QP3_EXCEEDANCE_PROBS, reverse=True)

    def test_censored_intervals_present(self) -> None:
        """Verify that QL has zeros (censored lower bounds) where QU > 0."""
        from tests.peakfqsa.fixtures.fortran_respec import QL, QU

        censored = [(l, u) for l, u in zip(QL, QU) if l != u]
        assert len(censored) == 3  # indices 20, 22, 23 have ql=0, qu=40


class TestSkewWeightingFixtures:
    """Validate skew_weighting fixture data is loadable."""

    def test_load_whist_cases(self) -> None:
        from tests.peakfqsa.fixtures.skew_weighting import load_whist_cases

        cases = load_whist_cases()
        assert len(cases) == 312  # 312 data rows in CSV

    def test_whist_case_values_in_range(self) -> None:
        from tests.peakfqsa.fixtures.skew_weighting import load_whist_cases

        cases = load_whist_cases()
        for case in cases:
            assert 0 < case.nu <= 1.0
            assert case.sigma2 > 0
            assert 0 < case.herl_factor < 2.0


class TestMomentsWymtFixtures:
    """Validate moments_wymt fixture metadata."""

    def test_expected_csv_files_exist(self) -> None:
        from tests.peakfqsa.fixtures.moments_wymt import (
            EXPECTED_DATA_CSV,
            EXPECTED_EMP_CSV,
            EXPECTED_INFO_CSV,
            EXPECTED_MGBT_CSV,
        )

        testdata_dir = os.path.normpath(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "..",
                "..",
                "_shared",
                "peakfqr",
                "inst",
                "testdata",
            )
        )
        for csv_name in [EXPECTED_INFO_CSV, EXPECTED_DATA_CSV, EXPECTED_EMP_CSV, EXPECTED_MGBT_CSV]:
            path = os.path.join(testdata_dir, csv_name)
            assert os.path.exists(path), f"Missing expected CSV: {path}"

    def test_column_definitions(self) -> None:
        from tests.peakfqsa.fixtures.moments_wymt import (
            MGBT_COLS,
            MOMENTS_COLS,
            QUANTILE_COLS,
            SITE_INFO_COLS,
        )

        assert "site_no" in SITE_INFO_COLS
        assert "site_no" in MOMENTS_COLS
        assert "site_no" in MGBT_COLS
        assert "site_no" in QUANTILE_COLS
