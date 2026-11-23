# Three USGS gages with published baseflow indexes

`<site>_daily.csv.gz`: daily mean flow (cfs) from the USGS Water Data API, retrieved 2026-11-23, over the period of record
used by USGS Open-File Report 03-146 (Wolock, 2003), whose `BFI_AVE` (the mean of the annual indexes of the USGS BFI
program of Wahl and Wahl, 5-day blocks and a factor of 0.9) is the expected value in `tests/test_regime.py`:

| Site | Days | BFI_AVE |
|---|---|---|
| 04131500 | 9,953 | 0.544 |
| 12065500 | 9,137 | 0.525 |
| 13077700 | 11,050 | 0.660 |

The number of days and the mean flow agree with the report's. They are three of 23 gages with unbroken records compared
on 2026-11-23: the mean annual index was above the published one by 0.016 on average with `ih_smoothed_minima`, and below
it by 0.001 with `ih_smoothed_minima_log` (largest differences 0.061 and 0.033). These three were picked for short
records and a difference between the two methods. They include 04131500, where the logarithms do worst of the 23 (0.033
low, against 0.015 high with straight lines); why was not traced.
