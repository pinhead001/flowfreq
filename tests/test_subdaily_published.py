"""flowfreq.subdaily against a published figure for a real gage.

Source
------
Exelon Generation Company, LLC, August 2012, *Final Study Report: Downstream
Flow Ramping and Stranding Study, RSP 3.8*, Conowingo Hydroelectric Project,
FERC Project Number 405. Section 4.1.1 is Spring Surveys, Hydraulic
Conditions, and Section 4.3.1 is Fall Surveys, Hydraulic Conditions.

Retrieved 2026-09-30 from the Maryland Department of the Environment:
https://mde.maryland.gov/programs/water/wetlandsandwaterways/documents/exelonmd/ferc/conowingo-frsp-3.08.pdf

The report states that its hydraulic data were "obtained from the USGS gage
at Conowingo Dam", which is USGS 01578310, Susquehanna River at Conowingo, MD.
The report's statements checked here, verbatim:

- S4.3.1: "Conowingo Dam discharge reached 141,000 cfs on October 3, 2010."
- S4.3.1, fall Studies 9-12 (October 27 and November 3, 10 and 17, 2010):
  "The morning peak discharge prior to the onset of three studies (9, 10,
  and 12) ranged from 26,100 cfs for Study 9 to 46,200 cfs for Study 10 ...
  Morning station discharge reached 80,000 cfs prior to the start of Study
  11." Also: "The morning generation peak typically occurred from 0600-0900 h
  prevailing time." The fall surveys ran "from approximately 1100-1500 h".
- S4.1.1, spring Studies 1-4 (April 29 and May 6, 13 and 18, 2010): "Daily
  peak discharge the prior day varied from 36,500 cfs for Study 1 to 80,900
  cfs for Study 3."

What agrees, and what is pinned
-------------------------------
:func:`flowfreq.subdaily.daily_extreme_timing` was run on the gage's approved
15-minute discharge, grouped on America/New_York local days. It reproduces
141,000, 26,100, 46,200 and 80,000 cfs and the Study 1 prior-day 36,500 cfs
**exactly**, at the report's own 100-cfs precision. The 80,000 cfs morning
peak of November 10 falls at 08:45 local, inside the report's 0600-0900 h.
That date is after the 2010-11-07 fall-back, so the report's "prevailing time"
is EST there and EDT on the other fall dates. The local-day grouping is
exercised on both sides of the transition.

What does not agree, and is recorded rather than pinned
-------------------------------------------------------
- **Study 3's prior-day peak.** The report says 80,900 cfs. The approved
  record's 2010-05-12 local-day maximum is 81,100 cfs at 23:15 EDT, and
  80,900 is the 00:00 reading on May 13. That is 200 cfs (0.25%) off, twice
  the report's precision. Possible causes are record revision after 2012, or
  the report reading a different instant. The difference is not explained, so
  it is not asserted.
- **The one-hour stage declines** (S4.3.1: 2.8, 4.2, 5.6 and 3.1 ft for
  Studies 9-12) are fixed-window changes over an hour. No
  :mod:`flowfreq.subdaily` function computes that: ``ramping_rates`` reports
  per-interval (15-minute) rates. A plain rolling one-hour maximum decline on
  the same gage-height record gives 2.76, 4.20, 5.60 and 3.29 ft. Three match
  at 0.1 ft; Study 12 does not. See TODO.md.
"""

from __future__ import annotations

import pandas as pd
import pytest

from flowfreq.subdaily import daily_extreme_timing
from flowfreq.usgs import USGSgage

SITE = "01578310"  # Susquehanna River at Conowingo, MD
TZ = "America/New_York"

pytestmark = pytest.mark.requires_network


def _day(timing: pd.DataFrame, date: str) -> pd.Series:
    return timing.set_index("date").loc[pd.Timestamp(date).date()]


@pytest.fixture(scope="module")
def fall_2010() -> pd.DataFrame:
    """The whole fall hydraulic period in one download: the Oct 3 spill
    through Study 12."""
    return USGSgage(SITE).download_instantaneous_flow("2010-10-01", "2010-11-18")


class TestConowingoFall2010:
    """Exelon (2012) RSP 3.8, S4.3.1."""

    def test_october_3_spill_peak(self, fall_2010: pd.DataFrame) -> None:
        timing = daily_extreme_timing(fall_2010, tz=TZ)
        row = _day(timing, "2010-10-03")
        assert row["max_value"] == 141_000.0
        assert bool(row["complete"])

    @pytest.mark.parametrize(
        ("study", "date", "published_cfs"),
        [(9, "2010-10-27", 26_100.0), (10, "2010-11-03", 46_200.0), (11, "2010-11-10", 80_000.0)],
    )
    def test_morning_peak_before_each_survey(
        self, fall_2010: pd.DataFrame, study: int, date: str, published_cfs: float
    ) -> None:
        """The report's surveys started ~1100 h, so its "morning peak" is
        the local-day maximum of readings before 11:00. On Oct 27 and Nov 3
        the day's overall maximum is a later afternoon peak, so the full
        day's maximum would not be the reported figure."""
        local = fall_2010.index.tz_convert(TZ)
        morning = fall_2010[local.hour < 11]
        row = _day(daily_extreme_timing(morning, tz=TZ), date)
        assert row["max_value"] == published_cfs, f"Study {study}"
        # "The morning generation peak typically occurred from 0600-0900 h."
        assert 6.0 <= row["hour_of_max"] <= 9.0, f"Study {study}"

    def test_study_11_peak_is_the_whole_days_maximum_in_est(self, fall_2010: pd.DataFrame) -> None:
        """Nov 10 is after the Nov 7 fall-back. Grouped on New York local days,
        the morning peak is the whole day's maximum, at 08:45 EST."""
        row = _day(daily_extreme_timing(fall_2010, tz=TZ), "2010-11-10")
        assert row["max_value"] == 80_000.0
        assert row["hour_of_max"] == pytest.approx(8.75)


class TestConowingoSpring2010:
    """Exelon (2012) RSP 3.8, S4.1.1."""

    def test_study_1_prior_day_peak(self) -> None:
        """ "Daily peak discharge the prior day ... 36,500 cfs for Study 1"
        (Study 1 was April 29, 2010, so the prior day is April 28)."""
        iv = USGSgage(SITE).download_instantaneous_flow("2010-04-27", "2010-04-29")
        row = _day(daily_extreme_timing(iv, tz=TZ), "2010-04-28")
        assert row["max_value"] == 36_500.0
        assert bool(row["complete"])
