"""Captured-shape NWIS RDB payloads for instantaneous-value retrieval tests.

These mirror the layout of the NWIS instantaneous-values service (``/nwis/iv/``)
so the parser can be exercised without network access. Discharge values are
plausible for a snowmelt-driven Cascades stream but are not real observations.
"""

# Site 12449950, four hours of 15-minute discharge on a single day in Pacific
# Daylight Time. The rising limb through the afternoon is the ordinary
# snowmelt diel signal.
IV_BASIC = """\
# ---------------------------------- WARNING ----------------------------------
# Provisional data are subject to revision.
#
# Data for the following 1 site(s) are contained in this file
#    USGS 12449950 METHOW RIVER AT PATEROS, WA
#
# Data provided for site 12449950
#    TS   parameter     Description
#    45   00060         Discharge, cubic feet per second
#
# Data-value qualification codes included in this output:
#     A  Approved for publication
#     P  Provisional data subject to revision
#
agency_cd\tsite_no\tdatetime\ttz_cd\t45_00060\t45_00060_cd
5s\t15s\t20d\t6s\t14n\t10s
USGS\t12449950\t2022-06-15 12:00\tPDT\t1520\tA
USGS\t12449950\t2022-06-15 12:15\tPDT\t1530\tA
USGS\t12449950\t2022-06-15 12:30\tPDT\t1545\tA
USGS\t12449950\t2022-06-15 12:45\tPDT\t1560\tA
USGS\t12449950\t2022-06-15 13:00\tPDT\t1590\tA
USGS\t12449950\t2022-06-15 13:15\tPDT\t1610\tP
"""

# The autumn daylight-saving transition: 2022-11-06 02:00 PDT became 01:00 PST,
# so local wall-clock time repeats the 01:00 hour and runs BACKWARDS across the
# boundary. On a UTC axis the same records are strictly increasing.
IV_DST_FALL_BACK = """\
# Data provided for site 12449950
#    TS   parameter     Description
#    45   00060         Discharge, cubic feet per second
#
agency_cd\tsite_no\tdatetime\ttz_cd\t45_00060\t45_00060_cd
5s\t15s\t20d\t6s\t14n\t10s
USGS\t12449950\t2022-11-06 01:15\tPDT\t410\tA
USGS\t12449950\t2022-11-06 01:30\tPDT\t409\tA
USGS\t12449950\t2022-11-06 01:45\tPDT\t408\tA
USGS\t12449950\t2022-11-06 01:00\tPST\t407\tA
USGS\t12449950\t2022-11-06 01:15\tPST\t406\tA
USGS\t12449950\t2022-11-06 01:30\tPST\t405\tA
"""

# A site reporting discharge from two separate sensors (two TS/DD numbers).
IV_MULTI_SENSOR = """\
# Data provided for site 12345678
agency_cd\tsite_no\tdatetime\ttz_cd\t45_00060\t45_00060_cd\t63680_00060\t63680_00060_cd
5s\t15s\t20d\t6s\t14n\t10s\t14n\t10s
USGS\t12345678\t2022-06-15 12:00\tPDT\t1520\tA\t1490\tA
USGS\t12345678\t2022-06-15 12:15\tPDT\t1530\tA\t1495\tA
"""

# Values NWIS could not report numerically: 'Ice' affected, and an empty cell.
IV_WITH_GAPS = """\
# Data provided for site 12449950
agency_cd\tsite_no\tdatetime\ttz_cd\t45_00060\t45_00060_cd
5s\t15s\t20d\t6s\t14n\t10s
USGS\t12449950\t2022-01-15 03:00\tPST\t210\tA
USGS\t12449950\t2022-01-15 03:15\tPST\tIce\tP
USGS\t12449950\t2022-01-15 03:30\tPST\t\tP
USGS\t12449950\t2022-01-15 03:45\tPST\t205\tA
"""

# A discharge column with a time-zone abbreviation not in NWIS_TZ_OFFSETS.
IV_UNKNOWN_TZ = """\
# Data provided for site 99999999
agency_cd\tsite_no\tdatetime\ttz_cd\t45_00060\t45_00060_cd
5s\t15s\t20d\t6s\t14n\t10s
USGS\t99999999\t2022-06-15 12:00\tXYZ\t100\tA
"""

# Header and format-spec row with no records.
IV_EMPTY = """\
# Data provided for site 12449950
agency_cd\tsite_no\tdatetime\ttz_cd\t45_00060\t45_00060_cd
5s\t15s\t20d\t6s\t14n\t10s
"""

# Rows present, but the requested parameter is absent (gage height only).
IV_WRONG_PARAMETER = """\
# Data provided for site 12449950
agency_cd\tsite_no\tdatetime\ttz_cd\t45_00065\t45_00065_cd
5s\t15s\t20d\t6s\t14n\t10s
USGS\t12449950\t2022-06-15 12:00\tPDT\t3.42\tA
"""

# Body of the HTTP 400 the instantaneous service returns for a window it holds
# no records for. It is a normal outcome, not a transport failure.
IV_NO_DATA_400_BODY = "No sites/data found using the selection criteria specified\n"

# Site-service seriesCatalogOutput listing both a daily ('dv') and an
# instantaneous ('uv') discharge series. The instantaneous record starts far
# later than the daily one, which is the usual case.
SITE_SERIES_CATALOG = """\
# US Geological Survey
#
agency_cd\tsite_no\tstation_nm\tdata_type_cd\tparm_cd\tbegin_date\tend_date\tcount_nu
5s\t15s\t50s\t2s\t5s\t20d\t20d\t8n
USGS\t12449950\tMETHOW RIVER AT PATEROS, WA\tdv\t00060\t1959-10-01\t2024-09-30\t23741
USGS\t12449950\tMETHOW RIVER AT PATEROS, WA\tuv\t00060\t2007-10-01\t2024-09-30\t6210
"""

# The same listing for a site with a daily record but no instantaneous series.
SITE_SERIES_CATALOG_NO_UV = """\
# US Geological Survey
#
agency_cd\tsite_no\tstation_nm\tdata_type_cd\tparm_cd\tbegin_date\tend_date\tcount_nu
5s\t15s\t50s\t2s\t5s\t20d\t20d\t8n
USGS\t12449950\tMETHOW RIVER AT PATEROS, WA\tdv\t00060\t1959-10-01\t2024-09-30\t23741
"""

# Site-service siteOutput=expanded response, the call that carries station name,
# drainage area and decimal-degree coordinates. Trimmed to the columns flowfreq
# reads; the real response has far more.
SITE_EXPANDED = """\
# US Geological Survey
#
agency_cd\tsite_no\tstation_nm\tdec_lat_va\tdec_long_va\tdrain_area_va
5s\t15s\t50s\t16s\t16s\t8s
USGS\t03606500\tBIG SANDY RIVER AT BRUCETON, TN\t36.0389722\t-88.2450000\t205
"""

# The same response from a site NWIS has no coordinates for: the columns are
# present but empty, which is how NWIS reports a missing value.
SITE_EXPANDED_NO_COORDS = """\
# US Geological Survey
#
agency_cd\tsite_no\tstation_nm\tdec_lat_va\tdec_long_va\tdrain_area_va
5s\t15s\t50s\t16s\t16s\t8s
USGS\t03606500\tBIG SANDY RIVER AT BRUCETON, TN\t\t\t205
"""


# Daily-values payloads (``/nwis/dv/``). Site 12449500, Methow River at Twisp.
# Five consecutive days of mean daily discharge in the summer recession.
DV_BASIC = """\
# ---------------------------------- WARNING ----------------------------------
# Provisional data are subject to revision.
#
# Data provided for site 12449500
#    TS   parameter     statistic     Description
# 152418  00060         00003         Discharge, cubic feet per second (Mean)
#
agency_cd\tsite_no\tdatetime\t152418_00060_00003\t152418_00060_00003_cd
5s\t15s\t20d\t14n\t10s
USGS\t12449500\t2025-09-26\t221\tA
USGS\t12449500\t2025-09-27\t215\tA
USGS\t12449500\t2025-09-28\t210\tA
USGS\t12449500\t2025-09-29\t209\tA
USGS\t12449500\t2025-09-30\t216\tA
"""

# What NWIS actually returns for a daily-values request carrying NO date
# range: a single row, the most recent day. It is a well-formed 200 response,
# which is what makes it dangerous -- nothing about it looks like an error.
# This fixture exists so the "always send a range" behaviour has a test that
# does not need the network to fail.
DV_SINGLE_DAY = """\
# Data provided for site 12449500
#    TS   parameter     statistic     Description
# 152418  00060         00003         Discharge, cubic feet per second (Mean)
#
agency_cd\tsite_no\tdatetime\t152418_00060_00003\t152418_00060_00003_cd
5s\t15s\t20d\t14n\t10s
USGS\t12449500\t2026-09-08\t261\tP
"""


#: An instantaneous gage-height (00065) payload, same station and shape as
#: IV_BASIC so the two parameters can be compared directly. Includes an "Ice"
#: reading, which is at least as common in a stage record as in a discharge one.
IV_STAGE_BASIC = """\
# ---------------------------------- WARNING ----------------------------------
# Provisional data are subject to revision.
#
# Data provided for site 12449950
#    TS   parameter     Description
#    46   00065         Gage height, feet
#
agency_cd\tsite_no\tdatetime\ttz_cd\t46_00065\t46_00065_cd
5s\t15s\t20d\t6s\t14n\t10s
USGS\t12449950\t2022-06-15 12:00\tPDT\t4.52\tA
USGS\t12449950\t2022-06-15 12:15\tPDT\t4.55\tA
USGS\t12449950\t2022-06-15 12:30\tPDT\tIce\tP
USGS\t12449950\t2022-06-15 12:45\tPDT\t4.48\tP
"""

#: Two 00065 series at one site. Separate primary and backup stage sensors are
#: more common than duplicate discharge series, so this is the case ts_id has to
#: resolve rather than an exotic one.
IV_STAGE_MULTI_SENSOR = """\
# Data provided for site 12345678
agency_cd\tsite_no\tdatetime\ttz_cd\t46_00065\t46_00065_cd\t63681_00065\t63681_00065_cd
5s\t15s\t20d\t6s\t14n\t10s\t14n\t10s
USGS\t12345678\t2022-06-15 12:00\tPDT\t4.52\tA\t4.51\tA
USGS\t12345678\t2022-06-15 12:15\tPDT\t4.55\tA\t4.54\tA
"""

# NWIS peak-flow RDB (``/nwis/peak``) for Big Sandy (03606500), trimmed. The
# three historic peaks carry an unknown day, which NWIS writes as ``-00``;
# their values match the Water Data OGC API and the committed Big Sandy
# fixture. The last two rows are synthetic: an unknown-month zero-flow year
# (the Orestimba pattern) and an unknown day in October, which must land in
# the *next* water year.
PEAK_PARTIAL_DATES = """\
# U.S. Geological Survey
# National Water Information System
#
# Station name : Big Sandy River at Bruceton, TN
# Drainage area: 205 square miles
#
agency_cd\tsite_no\tpeak_dt\tpeak_tm\tpeak_va\tpeak_cd\tgage_ht\tgage_ht_cd\tyear_last_pk
5s\t15s\t10d\t6s\t8s\t33s\t8s\t27s\t4s
USGS\t03606500\t1897-03-00\t\t25000\t7\t\t\t
USGS\t03606500\t1919-03-00\t\t21000\t7\t\t\t
USGS\t03606500\t1926-12-00\t\t18500\t7\t\t\t
USGS\t03606500\t1930-01-09\t\t9100\t\t\t\t
USGS\t03606500\t1931-00-00\t\t0\t\t\t\t
USGS\t03606500\t1931-10-00\t\t5000\t\t\t\t
"""
