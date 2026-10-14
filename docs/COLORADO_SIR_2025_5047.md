# Colorado SIR 2025-5047: second, opt-in equation set

Kohn, M.S., Mast, M.A., and Gross, T.A., 2026, *Peak-, mean-, and low-streamflow
regional-regression equations for natural streamflow in central and western Colorado, 2019*:
USGS Scientific Investigations Report 2025-5047, <https://doi.org/10.3133/sir20255047>.

## Where it is and how to get it

`flowfreq/data/regression/provisional/CO.json` holds the 32 peak-flow equations from Table 1:
four regions (Alpine, SubAlpine, MidElevation, Plateau), each at AEP 0.5 to 0.002. Load them with
`load_state("CO", allow_provisional=True)`. The default `load_state("CO")` still returns the
2009/2016/2006 set in `CO.json`, which is the set NSS serves.

This is a USGS report, not a flowfreq study. It is filed under `provisional/` only because that
folder is the library's one opt-in mechanism. The advisory `load_state` attaches to every
provisional set says "developed by flowfreq", and that is wrong for this file. The file's `notes`
field gives the correct provenance.

## NSS / StreamStats status (checked 2026-10-13)

- Not adopted. Public NSS (`streamstats.usgs.gov/nssservices`) still lists only citations 11, 20
  and 21 for Colorado.
- NSS has a draft of this report: citation 200, "2023, Peak-, Mean-, and Low-Streamflow ...
  2019", with regression regions 995-998 (`CO_Peakflow_*_SIR_2023_XXXX`). It has status 2,
  "Review", which means it is disabled for public users and enabled only on staging
  (`test.streamstats.usgs.gov`). The staging equation strings match Table 1 term for term, and
  the tests reproduce its returned values.
- When NSS approves it, review whether these equations should move into `CO.json`. In
  `CO.json` they would replace the Mountain, Northwest, Rio Grande and Southwest regions, and
  Foothills, Plains and Navajo would stay.

## What the schema stores and what it leaves out

- Form: `Q = 10^b0 · AREA^b1 · 10^(c1·X1 + c2·X2 [+ c3·X3])`. DRNAREA is `log10` and every other
  variable is `identity`, exactly as printed.
- Region selection is by mean basin elevation (ELEV, ft NAVD 88): Alpine > 10,490; Sub-Alpine
  9,493-10,490; Mid-Elevation 8,015-9,492; Plateau < 8,014. The library does not pick the region;
  the caller does. NSS staging uses slightly different, gap-free bounds. See the JSON notes.
- Variable codes are those printed in Table 1. NSS staging renames four of them: BSLDEM10M is
  ALTBSLDEM, ELEV is ALTELEV, I24H100Y is ALI24H100Y and PRECIP is PRECPRIS10.
- Limits are the single pooled range Table 1 prints for all regions. The report publishes no
  per-region ranges. Per Table 1.1, Alpine and Sub-Alpine gages span only about 0.8-1,450 mi2,
  against the pooled maximum of 22,629 mi2.
- Error statistics: SEP becomes `sep_log` and SME becomes `model_error_variance` (same conversion
  as `CO.json`). No covariance or per-site AVP is published, so intervals use `sep_log`. The data
  release (doi:10.5066/P9Q5AMFV), which the report says holds the final models, is "Under
  Revision" and unavailable.
- Not stored: mean, monthly, 7-day minimum and maximum and flow-duration equations (Tables 2-5),
  the drainage-area-ratio exponents (Table 6) and pseudo-R2.
- Basins in two regions: the report weights flows (not logs) by area fraction, so use
  `evaluate_weighted(..., space="linear")`.

## Checks

`tests/test_regression_colorado_2025.py` covers:

- Every coefficient, N, SEP and SME is typed from Table 1 and compared with the JSON. The JSON
  itself was parsed from the report's JATS XML (MathML) and matched against an independent parse
  of the PDF text layer, with zero disagreements.
- Evaluations at two input sets per region match NSS staging within its 3-significant-figure
  rounding. The report has no worked example or per-gage regression estimates to check against.
