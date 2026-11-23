# Monthly, summary-hydrograph and baseflow report sections, and the HYSEP windows

Approved by the owner 2026-11-22, with the HYSEP fix first.

## Intent

Reports have two kinds computed from daily flows (`low_flow`, `flow_duration`). Add three more from code flowfreq
already has (`regime`, `Hydrograph`): `monthly_statistics`, `summary_hydrograph`, `baseflow`. Each also gets a chat
tool (24 to 27). Non-goals: sharing one stored daily record between sections (each keeps its own copy, as today),
new separation methods, regulated-record detection.

## HYSEP fix (separate PR, first)

USGS WRIR 96-4040 (read 2026-11-22), "Hydrograph-separation methods": the interval 2N* is the odd integer from 3 to 11
nearest 2N, N = A^0.2. The fixed-interval bar is 2N* days wide, and the sliding and local-minimum methods look
0.5(2N*-1) days each way, so their window is 2N* too. `regime.py` computes 2N*, then doubles it (fixed 2x, the others
2x+1), and does not limit it to 3-11. Its docstring says the constant was never checked against the report.

Acceptance, from the report's French Creek example (01472157, 59.1 mi2, interval 5 days):

- text: fixed and sliding minima, April 1992: 49 cfs for April 5-9, 42 cfs on April 18;
- figure 4: local minima on April 9, 15, 21 and 24, 1992;
- Appendix 2, local minimum, 1969: monthly base flow 31.56 36.37 41.56 36.88 21.99 14.53 28.06 24.18 19.02 19.96
  19.72 30.48, annual 27.00 cfs. **The listing's dates run 11 days early:** moving today's published daily flows back
  11 days reproduces all 12 of its monthly streamflows to 0.00 cfs, and no other offset from -12 to 12 comes within
  1.8. With that offset, 5-day windows give every month within 0.24 cfs and 27.10 for the year; the present 11-day
  window is off by up to 8 cfs (25.06). The 0.4 percent left is not explained.

The report does not say how to break a tie (2N exactly even, as at 32 mi2); ours goes to the smaller interval. Port to
the 0.28 line by a PR into `release/0.28.x`; the owner decides on a 0.28.1.

## Sections

All take the station options of `low_flow` (`site_no`, `site_name`, `drainage_area_sqmi`, `daily`, `start_date`,
`end_date`), start with the record table, count negative values as missing, and refuse (`ReportError`) rather than
print an empty table.

- **`monthly_statistics`** (`annual_table`, default false): `monthly_flow_summary`, complete months only, October to
  September. Per month: years used; mean, median, standard deviation, maximum and minimum of the monthly mean flows
  (the last two with their water years); percent of annual volume. Box plot. `annual_table` adds the years by months
  table.
- **`summary_hydrograph`** (percentiles fixed at 10, 25, 50, 75, 90): the existing figure, and a table for the 1st
  and 15th of each month with years of record, minimum, the five percentiles and maximum. Percentiles are of
  non-exceedance and the text says so, since elsewhere q10 is the flow exceeded 10 percent of days. The figure's year
  count becomes water years.
- **`baseflow`** (`method`, default `ih_smoothed_minima`; `min_days`, default 350; `annual_table`, default true): the
  baseflow index pooled over complete water years, mean streamflow and baseflow (inches per year too when the area is
  known), the annual BFI table and bars, and a separated hydrograph for the complete year whose mean flow is nearest
  the median. Text: the method in words and cited, BFIs of different methods are not comparable, and regulation,
  diversion and snowmelt distort them.

## Approaches

Three kinds (chosen) against one combined `flow_regime` section (fewer downloads, but the parts could not be included
or ordered separately) and against extending `flow_duration` (one kind with unrelated outputs). HYSEP: fix it first
(chosen) against leaving it and refusing those methods in the section.

## The builder's calls

The three kinds go in a new `reporting/regime_sections.py`, not `daily.py` (510 lines), importing `daily`'s helpers.
Log axes when every plotted flow is positive, linear otherwise.
