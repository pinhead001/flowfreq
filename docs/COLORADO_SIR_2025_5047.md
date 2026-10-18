# Colorado SIR 2025-5047: second, opt-in equation set

Kohn, M.S., Mast, M.A., and Gross, T.A., 2026, *Peak-, mean-, and low-streamflow
regional-regression equations for natural streamflow in central and western Colorado, 2019*:
USGS Scientific Investigations Report 2025-5047, <https://doi.org/10.3133/sir20255047>.

## Where it is and how to get it

`flowfreq/data/regression/awaiting_nss/CO.json` holds the 32 peak-flow equations from Table 1:
four regions (Alpine, SubAlpine, MidElevation, Plateau), each at AEP 0.5 to 0.002. The file has
status `awaiting_nss`, which means published by the USGS but not yet served by NSS/StreamStats.

Load it with `load_state("CO", allow_awaiting_nss=True)`. The default `load_state("CO")` is
unchanged and returns the 2009/2016/2006 set in `CO.json`, which is the set NSS serves. The
advisory on this set reads "NOT YET IN STREAMSTATS (USGS-published) ...". It replaces the
"developed by flowfreq" wording that provisional sets carry.

### Why it has its own folder instead of `provisional/`

The loader binds each opt-in folder to one status. Its consistency check is "a file under
`<folder>/` must have status X". There are three reasons not to share `provisional/`:

- A provisional set (flowfreq-developed) and an awaiting-NSS set (USGS-published) are different
  claims, and each gets its own advisory.
- A state can have both at once, and each folder holds one `<CODE>.json`.
- Moving a set to `CO.json` once NSS adopts it is just a matter of moving the file.

If both opt-ins are requested, the USGS set wins, because a USGS report supersedes a flowfreq one.

## Keying: NSS staging citation 200, as is

NSS has this report only as a staging draft: citation 200, regression regions 995-998, status
2 "Review". That status means it is disabled for public NSS and StreamStats and visible only on
`test.streamstats.usgs.gov` (checked 2026-10-13). Public NSS still lists only citations 11, 20
and 21 for Colorado.

By the maintainer's decision, the file uses NSS's variable codes, units, limits and region
bounds, so StreamStats basin characteristics will feed it directly once NSS approves the
report. The coefficients are Table 1's, and NSS's equation strings match them term for term.
The report's own codes and bounds are kept for traceability in the `variable_codes` and
`region_selection` fields and in `notes`.

| Table 1 term | Report code | Stored (NSS) code | Units |
|---|---|---|---|
| AREA | DRNAREA | DRNAREA | mi2 (log10) |
| LONG | LONG_CENT | LONG_CENT | decimal degrees |
| SLOPE | BSLDEM10M | ALTBSLDEM | percent |
| MEANELEV | ELEV | ALTELEV | ft NAVD 88 |
| SWE | APR1_SWE | APR1_SWE | inches |
| P100Y24H / 100Y24H | I24H100Y | ALI24H100Y | inches |
| PPTANN | PRECIP | PRECPRIS10 | inches |

**SWE is in inches.** Table 1.2, the list of candidate characteristics, gives "snow water
equivalent (SNODAS) April 1" in meters, which raised the question of which unit the equations
use. They use inches:

- Every equation table, Tables 1-5, prints "SWE, snow-water
  equivalent on April 1 (in inches): 0–33.8" (Table 3: 0–32.3).
- A basin-mean April 1 SWE of 33.8 in (0.86 m) is ordinary for a high Colorado basin.
  33.8 m is impossible.
- The coefficients are sized for inches. Alpine's 1 % AEP coefficient, 0.0377 per unit in
  log10, spans 1.27 log units (a factor of about 19) over 0–33.8. Over a range in meters,
  0–0.86, it would move flows by under 10 %, which is too small to explain why SWE was
  selected.
- NSS staging declares `APR1_SWE` in inches with limits 0–33.8, and its evaluations match
  the stored equations.

Table 1.2's unit is the SNODAS grid's own. The report never says it converted, so this rests
on the four points above, not on a printed statement. Checked 2026-10-18;
`tests/test_regression_colorado_2025.py::test_swe_is_in_inches` pins it.

Region bounds on mean basin elevation:

- NSS (stored): Alpine `ALTELEV>=10490`; SubAlpine `[9492, 10490)`; MidElevation
  `[8014, 9492)`; Plateau `<8014`.
- Report: Alpine > 10,490; Sub-Alpine 9,493-10,490; Mid-Elevation 8,015-9,492; Plateau < 8,014.
  These leave 1-ft gaps and put exactly 10,490 ft in Sub-Alpine.

The library does not select a region; the caller does.

## What is not stored

- **Per-region limits:** the report publishes none. NSS's limits are the single pooled range
  from Table 1. Per Table 1.1, Alpine and Sub-Alpine gages span only about 0.8-1,450 mi2,
  against the pooled maximum of 22,629 mi2.
- **Covariance and per-site AVP:** neither is published, so intervals use `sep_log`. The data
  release (doi:10.5066/P9Q5AMFV) is "Under Revision" and unavailable.
- **Other tables:** Tables 2-6 (mean, monthly, low, 7-day maximum and flow-duration
  equations, and the DAR exponents) are not stored, and neither is pseudo-R2.
- **Two-region basins:** the report weights flows (not logs) by area fraction, so use
  `evaluate_weighted(..., space="linear")`.

## Checks

`tests/test_regression_colorado_2025.py` covers:

- **Table 1 values:** every coefficient, N, SEP and SME is typed from Table 1 and compared with
  the JSON. The JSON itself was parsed from the report's XML and matched against an independent
  parse of the PDF text, with zero disagreements.
- **NSS limits:** they equal Table 1's printed ranges to Table 1's rounding.
- **NSS evaluations:** at two input sets per region, results match NSS staging within its
  3-significant-figure rounding.

`tests/test_regression_awaiting_nss.py` covers the status itself: it loads only on opt-in, its
advisory text is correct, and an awaiting-NSS file in the wrong place, or with the wrong status
or state, is refused.
