"""Arizona peak-flow regression equations (issue #102).

Sources, as NSS combines them for Arizona's peak-flow statistic group:

- Paretti, Kennedy, Turney and Veilleux, 2014, USGS SIR 2014-5211: Table 9 (equations),
  Table 10 (errors), Table 11 (ranges, n) and Table 12 (MEV, covariance; in
  SIR2014-5211_tables.xlsx). Flood regions 1-4; region 5 is a schema gap.
- Waltemeyer, 2006, USGS SIR 2006-5306 (Navajo Nation): regions 8, 11, High Elevation.

Every literal number below is copied from those reports, or from a dated live NSS
response, never from ``AZ.json`` itself.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats

from flowfreq.regression import EquationsUnavailable, OutOfRangeError, evaluate, load_state

AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]
AEPS_NAVAJO = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.002]

# In range for regions 2-4 and the Navajo regions (region 1 needs NSS's ELEV >= 7,500 ft).
BASIN = {"CONTDA": 50.0, "ELEV": 5000.0, "PRECIP": 20.0, "DRNAREA": 50.0, "BSLDEM30ff": 0.1}


@pytest.fixture(scope="module")
def lib():
    return load_state("AZ")


def _sig3(x: float) -> float:
    return float(f"{x:.3g}")


class TestLibraryShape:
    def test_status_and_regions(self, lib):
        assert lib.status == "partial"
        assert lib.issue == 102
        assert lib.regions == ["1", "2", "3", "4", "Navajo11", "Navajo8", "NavajoHighElev"]
        for region in "1234":
            assert lib.aeps(region) == AEPS
        for region in ("Navajo11", "Navajo8", "NavajoHighElev"):
            assert lib.aeps(region) == AEPS_NAVAJO
        assert len(lib.equations) == 53

    def test_region_5_is_absent(self, lib):
        """Region 5's 10^(a - b DRNAREA^-c) has no schema representation."""
        with pytest.raises(EquationsUnavailable):
            lib.equation("5", 0.01)

    @pytest.mark.parametrize("region, n", [("1", 41), ("2", 85), ("3", 68), ("4", 77)])
    def test_n_sites_from_table_11(self, lib, region, n):
        assert {e.n_sites for e in lib.equations if e.region_code == region} == {n}

    @pytest.mark.parametrize(
        "region, limits",
        [
            ("1", {"CONTDA": (1.26, 711)}),
            # Table 11 prints "1,6017"; NSS uses 16,017.
            ("2", {"CONTDA": (0.103, 16017)}),
            ("3", {"CONTDA": (0.082, 1725), "PRECIP": (3.7, 22.2), "ELEV": (283, 6404)}),
            ("4", {"CONTDA": (0.059, 18044), "PRECIP": (10.8, 33.5), "ELEV": (3274, 7451)}),
        ],
    )
    def test_limits_from_table_11(self, lib, region, limits):
        eq = lib.equation(region, 0.1)
        assert {v.code: (v.minimum, v.maximum) for v in eq.variables} == limits

    def test_variables_by_aep(self, lib):
        assert [v.code for v in lib.equation("3", 0.02).variables] == ["CONTDA", "PRECIP"]
        assert [v.code for v in lib.equation("4", 0.5).variables] == ["CONTDA"]
        v = {v.code: v for v in lib.equation("4", 0.01).variables}["ELEV"]
        assert (v.transform, v.scale) == ("identity", 0.001)

    def test_out_of_range_raises(self, lib):
        with pytest.raises(OutOfRangeError, match="ELEV=7000.0 > max 6404"):
            evaluate(lib.equation("3", 0.5), {**BASIN, "ELEV": 7000.0})


class TestPrintedForms:
    @pytest.mark.parametrize(
        "region, aep, printed",
        [
            ("1", 0.5, lambda b: 17.4 * b["CONTDA"] ** 0.655),
            ("2", 0.005, lambda b: 1028 * b["CONTDA"] ** 0.413),
            (
                "3",
                0.5,
                lambda b: 2.78
                * b["CONTDA"] ** 0.462
                * b["PRECIP"] ** 2.229
                * 10 ** (-0.351 * b["ELEV"] / 1000),
            ),
            ("3", 0.002, lambda b: 384 * b["CONTDA"] ** 0.539 * b["PRECIP"] ** 0.758),
            (
                "4",
                0.01,
                lambda b: 30.0
                * b["CONTDA"] ** 0.605
                * b["PRECIP"] ** 1.805
                * 10 ** (-0.161 * b["ELEV"] / 1000),
            ),
        ],
    )
    def test_power_form(self, lib, region, aep, printed):
        assert evaluate(lib.equation(region, aep), BASIN).flow_cfs == pytest.approx(
            printed(BASIN), rel=1e-12
        )

    @pytest.mark.parametrize(
        "region, aep, sep_pct, sem_pct, avp",
        [
            ("1", 0.5, 86.1, 83.1, 0.105),
            ("2", 0.01, 67.3, 65.6, 0.070),
            ("4", 0.01, 27.1, 24.4, 0.013),
        ],
    )
    def test_error_statistics(self, lib, region, aep, sep_pct, sem_pct, avp):
        eq = lib.equation(region, aep)
        assert 100 * math.sqrt(10 ** (math.log(10) * eq.sep_log**2) - 1) == pytest.approx(
            sep_pct, abs=1e-3
        )
        assert eq.avp == avp
        # Table 12's full-precision MEV agrees with Table 10's SEM to its rounding.
        sem_log2 = math.log(1 + (sem_pct / 100) ** 2) / math.log(10) ** 2
        assert eq.model_error_variance == pytest.approx(sem_log2, abs=2e-3)


class TestCovariance:
    def test_region_4_q1_from_table_12(self, lib):
        eq = lib.equation("4", 0.01)
        assert eq.model_error_variance == 1.0872894521135056e-2
        assert eq.covariance[0] == (8.3730441e-2, -3.3286065e-4, -6.886383e-2, 1.4090869e-3)

    def test_worked_example_09508500(self, lib):
        """p. 34: region 4 Q1 at DRNAREA 5,499, PRECIP 19.6, ELEV 5,573 is 149,714 ft3/s;
        SEp,i 0.1105; 90% interval 98,000 to 229,000 ft3/s (t = 1.666, Table 12)."""
        site = {"CONTDA": 5499.0, "PRECIP": 19.6, "ELEV": 5573.0}
        eq = lib.equation("4", 0.01)
        est = evaluate(eq, site)
        assert round(est.flow_cfs) == 149714
        x = np.array([1.0, math.log10(5499.0), math.log10(19.6), 5.573])
        sep_i = math.sqrt(eq.model_error_variance + x @ np.array(eq.covariance) @ x)
        assert sep_i == pytest.approx(0.1105, abs=2e-4)
        t = 10 ** (1.666 * sep_i)
        assert round(est.flow_cfs / t, -3) == 98000
        assert round(est.flow_cfs * t, -3) == 229000
        # evaluate() uses t(0.95, n - p) with Table 11's n = 77 and p = 4.
        assert stats.t.ppf(0.95, 77 - 4) == pytest.approx(1.666, abs=1e-3)
        assert est.interval[1] == pytest.approx(est.flow_cfs * t, rel=1e-3)


# Live NSS responses, fetched 2026-10-02 (3 significant figures).
NSS_RECORDED = [
    (BASIN, "2", 0.5, 384),
    (BASIN, "2", 0.01, 4040),
    (BASIN, "3", 0.5, 237),
    (BASIN, "3", 0.01, 15700),
    (BASIN, "4", 0.5, 735),
    (BASIN, "4", 0.01, 11200),
    (BASIN, "Navajo8", 0.01, 6840),
    (BASIN, "Navajo11", 0.01, 5770),
    ({"CONTDA": 200.0}, "1", 0.5, 559),
    ({"CONTDA": 200.0}, "1", 0.01, 4130),
]


@pytest.mark.parametrize("site, region, aep, nss", NSS_RECORDED)
def test_worked_example_matches_recorded_nss(lib, site, region, aep, nss):
    assert _sig3(evaluate(lib.equation(region, aep), site).flow_cfs) == nss


NSS_CODE_TO_REGION = {
    "GC1618": "1",
    "GC1619": "2",
    "GC1620": "3",
    "GC1621": "4",
    "GC843": "Navajo8",
    "GC844": "NavajoHighElev",
    "GC845": "Navajo11",
}
NSS_STAT_TO_AEP = {
    "PK50AEP": 0.5,
    "PK20AEP": 0.2,
    "PK10AEP": 0.1,
    "PK4AEP": 0.04,
    "PK2AEP": 0.02,
    "PK1AEP": 0.01,
    "PK0_5AEP": 0.005,
    "PK0_2AEP": 0.002,
}


@pytest.mark.requires_network
@pytest.mark.parametrize(
    "basin",
    [
        BASIN,
        {"CONTDA": 10.0, "ELEV": 4000.0, "PRECIP": 12.0, "DRNAREA": 10.0, "BSLDEM30ff": 0.05},
        {"CONTDA": 200.0, "ELEV": 8000.0, "PRECIP": 25.0, "DRNAREA": 200.0, "BSLDEM30ff": 0.15},
    ],
)
def test_matches_live_nss(lib, basin):
    """Every stored region NSS (region 5, Arizona) returns equals flowfreq at 3 significant
    figures. GC1622 (region 5) is not stored and is skipped."""
    from flowfreq.streamstats import Characteristic, estimate_flow_statistics

    chars = {k: Characteristic(k, k, "", v, "") for k, v in basin.items()}
    results, _ = estimate_flow_statistics("AZ", chars, ["PFS"])
    compared = 0
    for r in results:
        if r.region_code not in NSS_CODE_TO_REGION:
            assert r.region_code == "GC1622"
            continue
        region = NSS_CODE_TO_REGION[r.region_code]
        for code, est in r.estimates.items():
            ours = evaluate(lib.equation(region, NSS_STAT_TO_AEP[code]), basin).flow_cfs
            assert _sig3(ours) == est.value, (region, code, est.value, ours)
            compared += 1
    assert compared > 0
