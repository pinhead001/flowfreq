# Washington Region 4 covariance: a GLS refit of SIR 2016-5118

Northwest studies plan, task 4 (`docs/superpowers/plans/2026-10-11-northwest-regional-studies.md`).
This settles the question PR #65 left open: in what basis Table 7's Region 4 covariance
matrix should be read. It also validates `flowfreq.regression.gls` on a published USGS
regression.

## Outcome

**Table 7's Region 4 matrix is the regression's in the basis `[1, log10 A, P/10]`.**
`WA.json` now stores it that way:
- PRECPRIS10 carries `scale = 0.1`, and its coefficient is stored ×10, so every estimate is
  unchanged;
- Table 7 is stored as printed;
- Washington is `verified`.

Two independent lines of evidence point to P/10:

| Check | raw P | P/10 |
|---|---|---|
| √(MEV + mean x′Ux) over the 139 calibration gages vs Table 6's Sp (AEP 0.5 / 0.01 / 0.002) | 2.59 / 2.82 / 2.88 × Sp | **1.005 / 1.003 / 1.004 × Sp** |
| P variance of the refit's (X′Λ⁻¹X)⁻¹ ÷ Table 7's printed value (median over AEP 0.5, 0.01, 0.002) | 0.012 | **1.23** |

**The trade-off.** The report applies the printed matrix to raw P in its own downstream work:
- in the Flood Q Tools workbook (cell BN4);
- in Table 8's WIE weighted estimates. Inverting them gives a regression variance matching
  raw P (median ratio 0.985 over 338 gage-AEP pairs; PR #65).

flowfreq follows the regression instead. Its Region 4 prediction intervals are therefore about
3× narrower than the workbook's. Its weighted estimates also give the regression more weight
than Table 8 does. `tests/test_regression_washington.py` pins both: the stored basis, and the
report's raw-P usage.

The USGS Washington Water Science Center could confirm this with an erratum, which would
settle it either way.

## Data

All inputs are the report's own (`tools/build_wa_gls_dataset.py`):

- **Table 5** (`sir20165118_table5.xlsx`): basin characteristics, regression region and
  the "used in regression" flag. 641 gages, 463 used: Region 1 93, Region 2 89, Region 3 142,
  Region 4 139.
- **Table 2:** each gage's historical period length and record span.
- **Table 8:** station EMA quantiles (line 1).
  - Table 8 prints Washington gages only. The 51 Idaho and Oregon calibration gages (Region 1
    22, Region 2 2, Region 4 27) are flowfreq EMA estimates: peaks through WY2014, PNW skew
    -0.07 with MSE 0.18, peak codes applied, Table 2's perception thresholds not applied.
  - On 18 Washington control gages the same computation reproduces Table 8 line 1 with
    a median ratio of 1.000 and an IQR within 1 % at AEP 0.5, 0.01 and 0.002. The few
    controls 5-10 % off at the tails are where the thresholds matter.

## GLS specification

`tools/wa_region4_refit.py`, using `flowfreq.regression.gls` (WREG, Eng and others 2009):

- y = log10 of the station quantile.
- The sampling variance uses the WREG LP3 form with:
  - the station SD smoothed by OLS on the region's design;
  - the station weighted skew;
  - skew MSE 0.18;
  - n = Table 2's historical period length.
- Cross-correlation uses the Pacific Northwest concurrent-peak model (SIR 2016-5083 eq.
  B11, Z = 0.21 + exp(-0.17 - 0.0058 D)), between basin centroids. Concurrent years come
  from the Table 2 record spans.
- The designs are Table 6's: Regions 1-2 `[1, log A, P, CAN]`, Region 3 `[1, log A, log P]`,
  Region 4 `[1, log A, P]`. Region 4 is also fit with `P/10`.

## Validation of the refit

| Region (m, EMA-supplemented) | AEP | β refit | β Table 6 | MEV refit / Table 7 | diag ratio refit / Table 7 |
|---|---|---|---|---|---|
| 1 (93, 22) | 0.5 | 0.559, 0.753, 0.032, -0.0078 | 0.585, 0.745, 0.032, -0.0078 | 0.118 / 0.113 | 1.21, 1.14, 1.01, 1.02 |
| 1 | 0.01 | 1.895, 0.665, 0.021, -0.0115 | 1.966, 0.649, 0.022, -0.013 | 0.119 / 0.108 | 2.11, 1.49, 1.09, 1.15 |
| 2 (89, 2) | 0.5 | 0.177, 0.874, 0.025, -0.0021 | 0.160, 0.877, 0.0253, -0.002 | 0.087 / 0.082 | 1.27, 1.15, 1.07, 1.05 |
| 2 | 0.01 | 1.674, 0.842, 0.0153, -0.0117 | 1.836, 0.814, 0.0151, -0.0132 | 0.118 / 0.116 | 2.41, 1.30, 1.13, 1.28 |
| 3 (142, 0) | 0.5 | -1.525, 0.901, 1.744 | -1.532, 0.899, 1.75 | 0.032 / 0.031 | 0.98, 1.03, 0.97 |
| 3 | 0.01 | -0.896, 0.884, 1.669 | -0.740, 0.898, 1.587 | 0.049 / 0.043 | 0.99, 1.20, 0.95 |
| 4 (139, 27), P/10 | 0.5 | 1.179, 0.910, 0.0733 | 1.186, 0.911, 0.073 | 0.046 / 0.045 | 1.14, 1.07, 1.03 |
| 4, P/10 | 0.01 | 1.929, 0.916, 0.0408 | 1.961, 0.921, 0.040 | 0.049 / 0.047 | 1.72, 1.22, 1.23 |

Regenerated 2026-10-31 after two `flowfreq.regression.gls` fixes: concurrent years are capped
at the shorter record, and the model error variance is searched only where Σ + γ²I is positive
definite (`docs/WA_REGION2_METHOW.md`). Every value above moved by at most 0.3% (MEV 0.14%), and
the Region 4 P/10 evidence is unchanged (P-entry ratio 1.227 → 1.228).

Intercepts and slopes trade off against each other, so the useful comparison is the
predicted flow at the calibration gages. The median |Δ log10 Q| between refit and published
equations at those gages:
- **Region 3**, the only region built entirely from printed Table 8 values: ≤ 0.027 at
  every AEP;
- **Region 4:** ≤ 0.044;
- **Regions 1 and 2:** ≤ 0.03 at AEP ≥ 0.04. They drift at rare AEPs (Region 1 up to 0.12
  at AEP 0.002), where the 22 EMA-supplemented Idaho/Oregon gages lack Table 2's
  thresholds.

MEV is within 0.89-1.13 of Table 7 everywhere. The covariance diagonal agrees within
0.95-1.3 in Region 3. The intercept variance runs high (up to ~2-3×) where gages are
EMA-supplemented or the cross-correlation is approximated. The P-variance entry on which
the Region 4 decision rests is the best determined of the three (1.03-1.35).

Reproduce:

```bash
python tools/build_wa_gls_dataset.py --cache <dir>         # sites CSV from Tables 2/5/8
python tools/build_wa_gls_dataset.py --cache <dir> --ema   # ID/OR supplement (Water Data API)
python tools/wa_region4_refit.py                           # refit JSON
```
