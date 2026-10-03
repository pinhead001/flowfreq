"""Turn a parsed PeakFQ ``.psf`` station into a Bulletin 17C analysis.

The other half of issue #31 (:mod:`flowfreq.psf` is the reader). Two steps,
each following the vendored peakfq 8.1.0 R code rather than reinterpreting it:

1. :func:`station_rows` builds the station's EMA records -- one
   ``(ql, qu, tl, tu, dtype)`` row per water year -- exactly as ``siteQT``
   does (``vendor/peakfqr/R/readInputs.R``). Peak qualification codes are
   applied through :func:`flowfreq.peak_codes.peak_frame_intervals`, the
   ``siteQT`` code rules, not a copy of them.
2. :func:`convert_station` maps the station's settings the way
   ``vendor/peakfqr/R/main.R`` does -- skew option and the ``r_G_mse`` sign
   encoding, low-outlier method, skew weighting, confidence level, AEP list --
   and returns a :class:`StationInputs`.

:class:`StationInputs` is not a parallel analysis structure. Its fields are
:class:`~flowfreq.bulletin17c.Bulletin17C`'s own constructor arguments, and
:meth:`StationInputs.emafit_arrays` is :class:`flowfreq.fortran_engine.EmafitArrays`,
the Fortran engine's input.

What ``Bulletin17C`` cannot express
-----------------------------------
``Bulletin17C``'s arguments -- exact systematic and historic peaks, interval
peaks (code 4/8 peaks and ``Interval`` lines), historic interval peaks (code 7
with 4 or 8), and ``(lower, upper)`` perception thresholds -- cover every row
kind ``siteQT`` builds, but not necessarily every combination a ``.psf`` can
produce (a ``MISSING`` threshold period with a peak inside it, for one).
Rather than approximate anything,
:meth:`StationInputs.bulletin17c_kwargs` rebuilds the rows each engine would
actually fit from the arguments -- ``build_emafit_arrays`` for
``engine="fortran"``, ``ExpectedMomentsAlgorithm._build_flow_intervals`` for
``engine="native"`` -- and raises :class:`UnsupportedSpecError`, naming the
years that differ, unless they match ``siteQT`` row for row.
:meth:`StationInputs.fortran_reference` always works (given the built
extension): it passes the ``siteQT`` rows to ``emafitpr`` directly, as
peakfq does.

Settings with no flowfreq equivalent (output-file and plotting options) are
logged as ignored. Nothing is dropped silently.
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from .core import PerceptionBound
from .fortran_engine import EmafitArrays, _gbthrsh0
from .peak_codes import Q_MAX, Q_MIN, peak_frame_intervals
from .psf import PsfFile, StationSpec, read_psf

if TYPE_CHECKING:
    from .bulletin17c import Bulletin17C
    from .core import FrequencyResults
    from .validation.reference import ReferenceResult

logger = logging.getLogger(__name__)

__all__ = [
    "PEAKFQ_AEPS",
    "PEAKFQ_AEPS_EXTENDED",
    "SiteRow",
    "StationInputs",
    "UnsupportedSpecError",
    "convert_peak_frame",
    "convert_psf",
    "convert_station",
    "read_psf_peaks",
    "station_rows",
]

#: AEPs peakfq 8.1.0 reports (``main.R``).
PEAKFQ_AEPS: Tuple[float, ...] = (
    0.9950,
    0.9900,
    0.9800,
    0.9750,
    0.9600,
    0.9500,
    0.9000,
    0.8000,
    0.7000,
    0.6667,
    0.6000,
    0.5704,
    0.5000,
    0.4292,
    0.4000,
    0.3000,
    0.2000,
    0.1000,
    0.0500,
    0.0400,
    0.0250,
    0.0200,
    0.0100,
    0.0050,
    0.0020,
)
#: AEPs with ``O Extended YES`` (``main.R``).
PEAKFQ_AEPS_EXTENDED: Tuple[float, ...] = PEAKFQ_AEPS + (0.0010, 0.0001)

#: ``WeightOpt`` names to ``emafitpr``'s ``wght_opt_n`` (``fortranWrappers.R``).
WEIGHT_OPTIONS: Dict[str, int] = {"HWN": 1, "ERL": 2, "INV": 3}

#: ``main.R`` initialises every skew MSE to this; for a station-skew site it is
#: what reaches ``emafitpr`` (below -98, so "station skew" in ``emafit``'s
#: encoding).
STATION_SKEW_MSE_SENTINEL = -1e99

#: Output options that shape peakfq's files and plots, not the analysis.
_OUTPUT_ONLY_OPTIONS = ("File", "Plot Style", "Plot Format", "Plot Position", "Export", "Empirical")

_DEFAULT_CONFIDENCE = 0.90


class UnsupportedSpecError(ValueError):
    """A ``.psf`` setting or record feature the requested engine cannot express."""


@dataclass(frozen=True)
class SiteRow:
    """One ``siteQT`` row, in cfs (not log10).

    ``ql == qu`` is an exactly known peak; ``dtype`` is 1 only for a peak
    carrying the historic-peak code 7.
    """

    year: int
    ql: float
    qu: float
    tl: float
    tu: float
    dtype: int

    @property
    def is_exact(self) -> bool:
        """Whether the peak is known exactly."""
        return self.ql == self.qu

    @property
    def flow(self) -> float:
        """The reported discharge of an exact row; ``Q_MIN`` reads as zero."""
        return 0.0 if self.ql <= Q_MIN else self.ql

    def log10(self) -> Tuple[int, float, float, float, float, int]:
        """``(year, log ql, log qu, log tl, log tu, dtype)``, ``emafitpr``'s units."""
        return (
            self.year,
            float(np.log10(self.ql)),
            float(np.log10(self.qu)),
            float(np.log10(self.tl)),
            float(np.log10(self.tu)),
            self.dtype,
        )


# --------------------------------------------------------------------------- #
# siteQT
# --------------------------------------------------------------------------- #


def _site_frame(peaks: pd.DataFrame, station_id: str) -> pd.DataFrame:
    """The rows of ``peaks`` for one station; all of them if there is no ``site_no``."""
    if "site_no" not in peaks.columns:
        return peaks
    ids = peaks["site_no"].astype(str).str.strip()
    return peaks.loc[ids == station_id]


def station_rows(
    spec: StationSpec,
    peaks: pd.DataFrame,
    *,
    include_urban_regulated: Optional[bool] = None,
    keep_no_info: bool = False,
) -> Tuple[SiteRow, ...]:
    """Build a station's EMA rows the way peakfq 8.1.0's ``siteQT`` does.

    In ``siteQT``'s order: the analysis period is the span of the perception
    thresholds, and peaks outside it are dropped; qualification codes set each
    peak's interval (:func:`flowfreq.peak_codes.peak_frame_intervals`); the
    ``.psf``'s ``Peak`` then ``Interval`` lines override the data for their
    years; thresholds apply in file order, later ones winning; a year with no
    peak is censored below its lower threshold, or is no information if that
    threshold is zero; every bound is floored at ``Q_MIN``; and a year whose
    interval is ``[Q_MIN, Q_MAX]`` carries no information and is removed.

    Parameters
    ----------
    spec : StationSpec
        The station's ``.psf`` block.
    peaks : pandas.DataFrame
        Peak frame (``water_year``, ``peak_flow_cfs``, optional
        ``qualification_code``). A ``site_no`` column, as
        :func:`flowfreq.watstore.read_watstore` returns, is filtered to
        ``spec.station_id``.
    include_urban_regulated : bool, optional
        Keep code 6/C peaks. Defaults to the station's ``Urb/Reg``.
    keep_no_info : bool, default False
        ``siteQT``'s ``keepNoInfo``: keep rows that carry no information.

    Returns
    -------
    tuple of SiteRow
        Sorted by water year.

    Raises
    ------
    ValueError
        No perception thresholds, no peak data, a negative discharge, a
        duplicate water year, a ``Peak``/``Interval`` year outside the
        analysis period, or a year the thresholds leave uncovered -- each a
        ``stop()`` in ``siteQT``.
    """
    sid = spec.station_id
    if not spec.thresholds:
        raise ValueError(f"station {sid} has no PCPT_Thresh lines; siteQT needs at least one")
    if include_urban_regulated is None:
        include_urban_regulated = spec.include_urban_regulated

    frame = _site_frame(peaks, sid)
    coded = peak_frame_intervals(frame, include_urban_regulated=include_urban_regulated)
    if coded.empty:
        raise ValueError(f"No peak flow data for site {sid}")
    if coded["water_year"].isna().any():
        raise ValueError(f"Some peak doesn't have associated year for site {sid}")

    start = min(t.start for t in spec.thresholds)
    end = max(t.end for t in spec.thresholds)
    years = coded["water_year"].astype(int)
    coded = coded.loc[(years >= start) & (years <= end)]
    dup = coded["water_year"].astype(int).duplicated()
    if dup.any():
        dups = sorted(coded.loc[dup, "water_year"].astype(int).unique())
        raise ValueError(f"Duplicated water year {dups} in input for site {sid}")

    ql: Dict[int, Optional[float]] = {y: None for y in range(start, end + 1)}
    qu: Dict[int, Optional[float]] = dict(ql)
    dtype: Dict[int, int] = {y: 0 for y in ql}
    for rec in coded.itertuples(index=False):
        year = int(getattr(rec, "water_year"))
        ql[year] = float(getattr(rec, "lower"))
        qu[year] = float(getattr(rec, "upper"))
        dtype[year] = 1 if bool(getattr(rec, "is_historic")) else 0

    for peak in spec.peaks:
        if peak.year not in ql:
            raise ValueError(f"Water year {peak.year} is not in analysis range for site {sid}")
        if peak.value < 0:
            hint = (
                " A legacy 'Peak YYYY -8888' excluded a year; write 'Interval YYYY 0 1E+20'."
                if peak.value == -8888
                else ""
            )
            raise ValueError(
                f"Negative discharge {peak.value} in specifications for water year "
                f"{peak.year} for site {sid}.{hint}"
            )
        ql[peak.year] = qu[peak.year] = float(peak.value)

    overridden = sorted({p.year for p in spec.peaks} & {i.year for i in spec.intervals})
    if overridden:
        logger.warning(
            "Station %s: Peak lines for %s are overwritten by Interval lines", sid, overridden
        )
    for interval in spec.intervals:
        if interval.year not in ql:
            raise ValueError(f"Water year {interval.year} is not in analysis range for site {sid}")
        ql[interval.year] = float(interval.lower)
        qu[interval.year] = float(interval.upper)

    tl: Dict[int, Optional[float]] = {y: None for y in ql}
    tu: Dict[int, Optional[float]] = dict(tl)
    for thr in spec.thresholds:
        for year in range(thr.start, thr.end + 1):
            tl[year], tu[year] = float(thr.lower), float(thr.upper)
    uncovered = [y for y in tl if tl[y] is None]
    if uncovered:
        raise ValueError(
            f"The following years in analysis period are missing perception thresholds: "
            f"{uncovered} for {sid}"
        )

    missing = [y for y in ql if ql[y] is None]
    zero_thr = [y for y in missing if (tl[y] or 0.0) <= Q_MIN]
    if zero_thr:
        logger.warning(
            "Station %s: zero perception threshold for missing years %s; "
            "these years are treated as missing data",
            sid,
            zero_thr,
        )
    rows: List[SiteRow] = []
    for year in sorted(ql):
        lo, hi, t_lo, t_hi = ql[year], qu[year], tl[year], tu[year]
        assert t_lo is not None and t_hi is not None
        if lo is None:
            lo, hi = (Q_MIN, Q_MAX) if t_lo <= Q_MIN else (Q_MIN, t_lo)
        assert hi is not None
        lo, hi, t_lo, t_hi = (max(v, Q_MIN) for v in (lo, hi, t_lo, t_hi))
        if lo == Q_MIN and hi == Q_MAX:
            t_lo = Q_MAX
        if t_lo >= Q_MAX and not keep_no_info:
            continue
        rows.append(SiteRow(year, lo, hi, t_lo, t_hi, dtype[year]))

    if not rows:
        logger.warning(
            "No data available after applying analysis specifications for site %s. "
            "Make sure the option to include urban/regulated peaks is set appropriately.",
            sid,
        )
    elif len(rows) < 10:
        logger.warning(
            "Less than 10 years of data for site %s; Bulletin 17C recommends at least 10", sid
        )
    return tuple(rows)


# --------------------------------------------------------------------------- #
# Settings
# --------------------------------------------------------------------------- #


def _skew_settings(spec: StationSpec) -> Tuple[str, float, float]:
    """``(option, r_G, r_G_mse)`` exactly as ``main.R`` passes them to ``emafit``.

    Station: ``(0, -1e99)``. Weighted: ``(GenSkew, SkewSE**2)``.
    Regional/Generalized: ``(GenSkew, -SkewSE**2)`` -- the negative sign is
    ``emafitpr``'s "generalized skew with this MSE" encoding.
    """
    sid = spec.station_id
    raw = spec.skew_option
    if raw is None:
        raise ValueError(f"station {sid} has no SkewOpt")
    option = raw.lower()
    # main.R accepts "regional"; readInputs.R's documentation and the RHOT
    # app say "Generalized". Same thing, so both are read.
    if option == "regional":
        option = "generalized"
    if option not in ("station", "weighted", "generalized"):
        raise ValueError(f"Invalid skew option {raw!r} for station {sid}")
    if option == "station":
        return option, 0.0, STATION_SKEW_MSE_SENTINEL
    gen, se = spec.generalized_skew, spec.generalized_skew_se
    if gen is None or not np.isfinite(gen):
        raise ValueError(f"Invalid regional skew coefficient at site {sid}")
    if se is None or not np.isfinite(se):
        raise ValueError(f"Invalid regional skew standard error at site {sid}")
    if se <= 0:
        raise ValueError(f"Input skew standard error is not positive at site {sid}")
    mse = se * se
    return option, float(gen), mse if option == "weighted" else -mse


def _low_outlier_settings(spec: StationSpec) -> Tuple[str, Optional[float]]:
    """``(LOType, user threshold)``, per ``main.R`` + ``fortranWrappers.R::emafit``.

    MGBT gives ``None`` (run the test). FIXED gives ``LoThresh`` in cfs. NONE
    gives ``Q_MIN``, which ``emafit``/:func:`flowfreq.fortran_engine._gbthrsh0`
    read as "no low-outlier test".
    """
    sid = spec.station_id
    lo_type = spec.fields.get("LOType")
    if lo_type is None:
        raise ValueError(f"station {sid} has no LOType (MGBT, FIXED or NONE)")
    lo_type = lo_type.upper()
    if lo_type == "MGBT":
        return lo_type, None
    if lo_type == "NONE":
        return lo_type, Q_MIN
    if lo_type == "FIXED":
        raw = spec.fields.get("LoThresh")
        try:
            value = float(raw) if raw is not None else float("nan")
        except ValueError:
            value = float("nan")
        if not np.isfinite(value):
            raise ValueError(f"Invalid low outlier threshold {raw!r} for station {sid}")
        return lo_type, value
    raise ValueError(f"Invalid low outlier option {lo_type!r} for station {sid}")


def _weight_option(spec: StationSpec) -> str:
    raw = spec.fields.get("WeightOpt")
    option = "HWN" if raw is None else raw.upper()
    if option not in WEIGHT_OPTIONS:
        raise ValueError(
            f"Invalid skew weighting option {raw!r} for station {spec.station_id}; "
            f"valid options are {sorted(WEIGHT_OPTIONS)}"
        )
    return option


@dataclass(frozen=True)
class _OutputSettings:
    confidence: float = _DEFAULT_CONFIDENCE
    aeps: Tuple[float, ...] = PEAKFQ_AEPS
    method: str = "ema"


def _output_settings(psf: Optional[PsfFile]) -> _OutputSettings:
    """Read the ``O`` lines that change the analysis; log the ones that cannot."""
    if psf is None:
        return _OutputSettings()
    options = dict(psf.output_options)

    confidence = _DEFAULT_CONFIDENCE
    raw_ci = options.pop("ConfInterval", None)
    if raw_ci is not None:
        try:
            value = float(raw_ci)
        except ValueError:
            value = float("nan")
        if 0 < value < 1:
            confidence = value
        else:
            logger.warning(
                "Invalid confidence interval option %r; defaulting to %s", raw_ci, confidence
            )

    aeps = PEAKFQ_AEPS
    for key in [k for k in options if k.upper() == "EXTENDED"]:
        raw = options.pop(key)
        if key != "EXTENDED":
            logger.warning(
                "O %s %s: honoured, but peakfq 8.1.0 itself only reads the spelling "
                "'EXTENDED' and would ignore this line",
                key,
                raw,
            )
        if raw.upper() == "YES":
            aeps = PEAKFQ_AEPS_EXTENDED
        elif raw.upper() != "NO":
            logger.warning("Invalid option O %s %s; using the standard AEPs", key, raw)

    method = "ema"
    raw_ema = options.pop("EMA", None)
    if raw_ema is not None and raw_ema.upper() == "NO":
        # peakfq 7's "EMA NO" meant a Bulletin 17B analysis. flowfreq has no
        # 17B procedure; its method of moments is the nearest thing, and not
        # the same thing, so say so.
        logger.warning(
            "O EMA NO: using flowfreq's method of moments, which is not PeakFQ 7's "
            "Bulletin 17B procedure (peakfq 8.1.0 has no non-EMA analysis at all)"
        )
        method = "mom"
    elif raw_ema is not None and raw_ema.upper() != "YES":
        raise ValueError(f"Invalid option O EMA {raw_ema!r}; expected YES or NO")

    for key, text in options.items():
        if key in _OUTPUT_ONLY_OPTIONS:
            logger.info("O %s %s: output option with no flowfreq equivalent; ignored", key, text)
        else:
            logger.warning("O %s %s: unknown output option; ignored", key, text)
    return _OutputSettings(confidence=confidence, aeps=aeps, method=method)


# --------------------------------------------------------------------------- #
# The converted station
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class StationInputs:
    """A ``.psf`` station converted to Bulletin 17C analysis arguments.

    The ``peak_flows`` ... ``user_low_outlier_threshold`` fields are
    :class:`~flowfreq.bulletin17c.Bulletin17C`'s constructor arguments. Use
    :meth:`bulletin17c_kwargs` rather than reading them directly: it checks
    that the engine you name would fit exactly :attr:`rows`.

    Attributes
    ----------
    station_id : str
        The ``.psf`` station ID, e.g. ``"06326500.00"``.
    rows : tuple of SiteRow
        The ``siteQT`` rows; the authority for everything else here.
    peak_flows, water_years : tuple
        Exactly known, non-historic peaks (zero flows as 0.0).
    historical_peaks : tuple of (int, float)
        Exactly known peaks carrying code 7.
    perception_thresholds : dict of (int, int) to float or (float, float)
        The ``PCPT_Thresh`` bounds, in ``.psf`` order (later wins): the lower
        bound alone when the upper is ``1E+20``, else ``(lower, upper)``.
    regional_skew, regional_skew_mse : float or None
        ``None`` for a station-skew analysis; otherwise ``GenSkew`` and the
        signed MSE (``SkewSE**2``, negated for a generalized skew), per
        ``main.R`` and the skew-MSE encoding in ``CLAUDE.md``.
    skew_option : str
        ``"station"``, ``"weighted"`` or ``"generalized"``.
    user_low_outlier_threshold : float or None
        ``None`` runs MGBT; a cfs value fixes the threshold; ``Q_MIN`` turns
        the test off (``LOType NONE``).
    low_outlier_method : str
        ``"MGBT"``, ``"FIXED"`` or ``"NONE"``.
    weight_option : str
        ``"HWN"``, ``"ERL"`` or ``"INV"``.
    include_urban_regulated : bool
        ``Urb/Reg``.
    method : str
        ``"ema"``, or ``"mom"`` for ``O EMA NO``.
    confidence : float
        ``O ConfInterval``, default 0.90.
    aeps : tuple of float
        peakfq's reporting AEPs (extended with ``O EXTENDED YES``).
    interval_peaks : tuple of (int, float, float)
        Non-historic years whose peak is known only as an interval
        ``(water_year, ql, qu)`` in cfs -- code 4/8 peaks and ``Interval``
        lines -- except the ``(Q_MIN, tl)`` rows ``Bulletin17C`` already
        builds for an unobserved year from its perception threshold.
    historical_interval_peaks : tuple of (int, float, float)
        Historic (code 7) years whose peak is known only as an interval --
        code 7 with 4 or 8, or an ``Interval`` line over a code 7 peak:
        ``siteQT``'s censored ``dtype = 1`` rows.
    """

    station_id: str
    rows: Tuple[SiteRow, ...]
    peak_flows: Tuple[float, ...]
    water_years: Tuple[int, ...]
    historical_peaks: Tuple[Tuple[int, float], ...]
    perception_thresholds: Dict[Tuple[int, int], PerceptionBound]
    regional_skew: Optional[float]
    regional_skew_mse: Optional[float]
    skew_option: str
    user_low_outlier_threshold: Optional[float]
    low_outlier_method: str
    weight_option: str = "HWN"
    include_urban_regulated: bool = False
    method: str = "ema"
    confidence: float = _DEFAULT_CONFIDENCE
    aeps: Tuple[float, ...] = PEAKFQ_AEPS
    interval_peaks: Tuple[Tuple[int, float, float], ...] = ()
    historical_interval_peaks: Tuple[Tuple[int, float, float], ...] = ()
    _r_g: float = field(default=0.0, repr=False)
    _r_g_mse: float = field(default=STATION_SKEW_MSE_SENTINEL, repr=False)

    # ------------------------------------------------------------------ #

    @property
    def begin_year(self) -> int:
        """First water year with a row."""
        return self.rows[0].year

    @property
    def end_year(self) -> int:
        """Last water year with a row."""
        return self.rows[-1].year

    @property
    def n_historic(self) -> int:
        """Rows carrying the historic-peak code (``HistPeaks`` in peakfq's output)."""
        return sum(r.dtype for r in self.rows)

    @property
    def n_exact(self) -> int:
        """Exactly known peaks, historic ones included."""
        return sum(1 for r in self.rows if r.is_exact)

    def emafit_arrays(self) -> EmafitArrays:
        """The ``siteQT`` rows as ``emafitpr`` input (log10).

        Returns
        -------
        EmafitArrays
        """
        logged = [r.log10() for r in self.rows]
        systematic = {r.year: r.flow for r in self.rows if r.is_exact and r.dtype == 0}
        return EmafitArrays(
            ql=np.array([r[1] for r in logged]),
            qu=np.array([r[2] for r in logged]),
            tl=np.array([r[3] for r in logged]),
            tu=np.array([r[4] for r in logged]),
            dtype=np.array([r[5] for r in logged], dtype=np.int32),
            years=np.array([r[0] for r in logged], dtype=int),
            systematic_peaks=systematic,
            n_zeros=sum(1 for v in systematic.values() if v == 0.0),
            n_censored=sum(1 for r in self.rows if not r.is_exact),
            gbthrsh0=_gbthrsh0(self.user_low_outlier_threshold),
        )

    # ------------------------------------------------------------------ #

    def _raw_kwargs(self) -> Dict[str, Any]:
        return {
            "peak_flows": np.array(self.peak_flows, dtype=float),
            "water_years": np.array(self.water_years, dtype=int),
            "regional_skew": self.regional_skew,
            "regional_skew_mse": self.regional_skew_mse,
            "historical_peaks": list(self.historical_peaks) or None,
            "perception_thresholds": dict(self.perception_thresholds) or None,
            "user_low_outlier_threshold": self.user_low_outlier_threshold,
            "interval_peaks": list(self.interval_peaks) or None,
            "historical_interval_peaks": list(self.historical_interval_peaks) or None,
        }

    def _engine_rows(self, engine: str) -> List[Tuple[int, float, float, float, float, int]]:
        """The rows ``engine`` would build from :meth:`_raw_kwargs`, in log10."""
        kw = self._raw_kwargs()
        if engine == "fortran":
            from .fortran_engine import build_emafit_arrays

            arrays = build_emafit_arrays(
                kw["peak_flows"],
                water_years=kw["water_years"],
                historical_peaks=kw["historical_peaks"],
                perception_thresholds=kw["perception_thresholds"],
                user_low_outlier_threshold=kw["user_low_outlier_threshold"],
                interval_peaks=kw["interval_peaks"],
                historical_interval_peaks=kw["historical_interval_peaks"],
            )
            return [
                (int(y), float(a), float(b), float(c), float(d), int(e))
                for y, a, b, c, d, e in zip(
                    arrays.years, arrays.ql, arrays.qu, arrays.tl, arrays.tu, arrays.dtype
                )
            ]
        from .bulletin17c import ExpectedMomentsAlgorithm

        ema = ExpectedMomentsAlgorithm(
            kw["peak_flows"],
            water_years=kw["water_years"],
            regional_skew=kw["regional_skew"],
            regional_skew_mse=kw["regional_skew_mse"],
            historical_peaks=kw["historical_peaks"],
            perception_thresholds=kw["perception_thresholds"],
            interval_peaks=kw["interval_peaks"],
            historical_interval_peaks=kw["historical_interval_peaks"],
        )
        out = []
        for iv in ema._build_flow_intervals(0.0):
            out.append(
                SiteRow(
                    year=int(iv.year),
                    ql=max(float(iv.lower), Q_MIN),
                    qu=max(float(iv.upper), Q_MIN),
                    tl=max(float(iv.perception_threshold), Q_MIN),
                    tu=max(float(iv.perception_upper), Q_MIN),
                    dtype=int(iv.is_historical),
                ).log10()
            )
        return out

    def unsupported_reasons(self, engine: str = "native") -> List[str]:
        """Why ``Bulletin17C(..., engine)`` would not reproduce peakfq here.

        Parameters
        ----------
        engine : {"native", "fortran"}
            The ``Bulletin17C.run_analysis`` engine.

        Returns
        -------
        list of str
            Empty when the analysis is exactly representable.

        Raises
        ------
        ValueError
            On an unknown engine.
        """
        if engine not in ("native", "fortran"):
            raise ValueError(f"engine must be 'native' or 'fortran', got {engine!r}")
        reasons: List[str] = []
        weighted = self.skew_option == "weighted"
        if weighted and self.weight_option != "HWN":
            reasons.append(
                f"WeightOpt {self.weight_option}: Bulletin17C implements HWN only "
                "(use fortran_reference())"
            )
        if self.method == "mom":
            if engine == "fortran":
                reasons.append("O EMA NO: the Fortran engine has no method-of-moments path")
            if self.historical_peaks or any(not r.is_exact for r in self.rows):
                reasons.append(
                    "O EMA NO: the method of moments cannot use historic or censored data"
                )
        if engine == "native" and self.skew_option == "generalized":
            reasons.append(
                "SkewOpt Generalized: the native engine has no generalized-skew path "
                "(negative r_G_mse); use engine='fortran'"
            )
        if reasons:
            return reasons

        expected = Counter(r.log10() for r in self.rows)
        actual = Counter(self._engine_rows(engine))
        if expected != actual:
            differ = sorted({row[0] for row in (expected - actual) + (actual - expected)})
            reasons.append(
                f"the {engine} engine would build different EMA rows than siteQT "
                f"for water years {differ}; use fortran_reference()"
            )
        return reasons

    def bulletin17c_kwargs(self, engine: str = "native") -> Dict[str, Any]:
        """``Bulletin17C`` constructor arguments, checked against ``siteQT``.

        Parameters
        ----------
        engine : {"native", "fortran"}
            The engine the arguments are for; the check differs because the
            two engines build their EMA rows differently.

        Returns
        -------
        dict

        Raises
        ------
        UnsupportedSpecError
            If the engine would not fit exactly :attr:`rows` with peakfq's
            settings; the message lists every reason.
        """
        reasons = self.unsupported_reasons(engine)
        if reasons:
            raise UnsupportedSpecError(f"station {self.station_id}: " + "; ".join(reasons))
        return self._raw_kwargs()

    def to_bulletin17c(self, engine: str = "native") -> "Bulletin17C":
        """A :class:`~flowfreq.bulletin17c.Bulletin17C` for this station.

        Parameters
        ----------
        engine : {"native", "fortran"}

        Returns
        -------
        Bulletin17C

        Raises
        ------
        UnsupportedSpecError
            See :meth:`bulletin17c_kwargs`.
        """
        from .bulletin17c import Bulletin17C

        return Bulletin17C(**self.bulletin17c_kwargs(engine))

    def run(self, engine: str = "native") -> "FrequencyResults":
        """Fit the station and report at peakfq's AEPs and confidence level.

        Parameters
        ----------
        engine : {"native", "fortran"}

        Returns
        -------
        FrequencyResults
            ``quantiles`` and ``confidence_limits`` evaluated at
            :attr:`aeps` and :attr:`confidence`.

        Raises
        ------
        UnsupportedSpecError
            See :meth:`bulletin17c_kwargs`.
        """
        b17c = self.to_bulletin17c(engine)
        results = b17c.run_analysis(method=self.method, engine=engine)
        aep = np.array(self.aeps, dtype=float)
        results.quantiles = b17c.compute_quantiles(aep=aep)
        results.confidence_limits = b17c.compute_confidence_limits(
            aep=aep, confidence=self.confidence
        )
        return results

    def fortran_reference(self) -> "ReferenceResult":
        """Run ``emafitpr`` on the ``siteQT`` rows with ``main.R``'s arguments.

        Handles everything a ``.psf`` can say -- intervals, upper thresholds,
        generalized skew, ERL/INV weighting -- because it bypasses
        ``Bulletin17C``, as peakfq does.

        Returns
        -------
        ReferenceResult

        Raises
        ------
        ImportError
            If the f2py extension is not built.
        UnsupportedSpecError
            For ``O EMA NO``: ``emafitpr`` is EMA only.
        """
        if self.method != "ema":
            raise UnsupportedSpecError("O EMA NO: emafitpr has no method-of-moments path")
        from .validation.reference import ReferenceResult

        arrays = self.emafit_arrays()
        return ReferenceResult.from_emafit(
            ql=arrays.ql.tolist(),
            qu=arrays.qu.tolist(),
            tl=arrays.tl.tolist(),
            tu=arrays.tu.tolist(),
            dtype=arrays.dtype.tolist(),
            aeps=list(self.aeps),
            regional_skew=self._r_g,
            regional_skew_mse=self._r_g_mse,
            station_name=self.station_id,
            eps=self.confidence,
            weight_opt=WEIGHT_OPTIONS[self.weight_option],
            gbthrsh0=arrays.gbthrsh0,
        )


def _thresholds_in_order(spec: StationSpec) -> Dict[Tuple[int, int], PerceptionBound]:
    """``PCPT_Thresh`` bounds keyed by period, preserving "later wins".

    A threshold whose upper bound is ``Q_MAX`` or more is its lower bound
    alone (a float, as ``Bulletin17C`` has always taken it); one with a finite
    upper bound is the ``(lower, upper)`` pair.

    A repeated period is moved to the end, so its last value keeps the
    priority ``siteQT`` gives it over every threshold listed before it.

    A lower bound of ``Q_MAX`` or more ("MISSING" periods) is left out: in
    ``siteQT`` such a year carries no information and is removed, which is
    what an undeclared (vacuous) year gets from ``Bulletin17C`` too, whereas
    passing it would censor the year below 1E+20. Whether that equivalence
    holds for a given station -- it does not if a peak falls inside the period
    -- is what :meth:`StationInputs.unsupported_reasons` verifies.
    """
    out: Dict[Tuple[int, int], PerceptionBound] = {}
    for thr in spec.thresholds:
        key = (thr.start, thr.end)
        out.pop(key, None)
        if thr.lower < Q_MAX:
            if thr.upper < Q_MAX:
                out[key] = (float(thr.lower), float(thr.upper))
            else:
                out[key] = float(thr.lower)
    return out


def convert_station(
    spec: StationSpec,
    peaks: pd.DataFrame,
    psf: Optional[PsfFile] = None,
) -> StationInputs:
    """Convert one ``.psf`` station and its peaks into analysis inputs.

    Parameters
    ----------
    spec : StationSpec
        The station block.
    peaks : pandas.DataFrame
        Peak frame; see :func:`station_rows`.
    psf : PsfFile, optional
        The whole file, for its ``O`` options (confidence level, extended
        AEPs, ``EMA``). Without it peakfq's defaults apply.

    Returns
    -------
    StationInputs

    Raises
    ------
    ValueError
        On an invalid or missing setting that makes ``main.R`` stop -- skew
        option, regional skew or its standard error, ``LOType``/``LoThresh``,
        ``WeightOpt`` -- or on any ``siteQT`` error (:func:`station_rows`).
    """
    option, r_g, r_g_mse = _skew_settings(spec)
    lo_method, lo_threshold = _low_outlier_settings(spec)
    weight = _weight_option(spec)
    out = _output_settings(psf)
    rows = station_rows(spec, peaks)

    n_exact = sum(1 for r in rows if r.is_exact)
    if len(rows) < 10 or n_exact < 8:
        logger.warning(
            "Station %s has %d rows and %d exact peaks; peakfq 8.1.0 skips sites with "
            "fewer than 10 rows or 8 uncensored peaks",
            spec.station_id,
            len(rows),
            n_exact,
        )

    systematic = [r for r in rows if r.is_exact and r.dtype == 0]
    historic = [r for r in rows if r.is_exact and r.dtype == 1]
    # A (Q_MIN, tl) row is what Bulletin17C builds for an unobserved year in a
    # perception period, so it needs no argument; any other censored,
    # non-historic row is an interval peak. Either way the row is the same.
    intervals = [
        r for r in rows if not r.is_exact and r.dtype == 0 and not (r.ql <= Q_MIN and r.qu == r.tl)
    ]
    # A censored historic row is always an observed peak (siteQT only flags a
    # year historic from its peak code), so it needs no such exclusion.
    historic_intervals = [r for r in rows if not r.is_exact and r.dtype == 1]
    return StationInputs(
        station_id=spec.station_id,
        rows=rows,
        peak_flows=tuple(r.flow for r in systematic),
        water_years=tuple(r.year for r in systematic),
        historical_peaks=tuple((r.year, r.flow) for r in historic),
        perception_thresholds=_thresholds_in_order(spec),
        regional_skew=None if option == "station" else r_g,
        regional_skew_mse=None if option == "station" else r_g_mse,
        skew_option=option,
        user_low_outlier_threshold=lo_threshold,
        low_outlier_method=lo_method,
        weight_option=weight,
        include_urban_regulated=spec.include_urban_regulated,
        method=out.method,
        confidence=out.confidence,
        aeps=out.aeps,
        interval_peaks=tuple((r.year, r.ql, r.qu) for r in intervals),
        historical_interval_peaks=tuple((r.year, r.ql, r.qu) for r in historic_intervals),
        _r_g=r_g,
        _r_g_mse=r_g_mse,
    )


def convert_peak_frame(
    peaks: pd.DataFrame,
    *,
    station_id: str = "site",
    regional_skew: Optional[float] = None,
    regional_skew_mse: Optional[float] = None,
    low_outlier_threshold: Optional[float] = None,
    include_urban_regulated: bool = False,
) -> StationInputs:
    """Analysis inputs for a bare peak frame, with its qualification codes applied.

    The ``.psf``-free route into :func:`convert_station`: the station gets
    the one default perception threshold PeakFQ writes for a new site
    (``PCPT_Thresh <first WY> <last WY> 0 1E+20``), so the codes are applied
    exactly as ``siteQT`` applies them -- code 7 peaks become historic, 3/O
    (and 6/C unless ``include_urban_regulated``) are removed, 4/8 become
    censored intervals -- and a water year with no peak is missing data.

    Parameters
    ----------
    peaks : pandas.DataFrame
        ``water_year``, ``peak_flow_cfs`` and (optionally)
        ``qualification_code``, as :meth:`flowfreq.usgs.USGSgage.download_peak_flow`
        returns them.
    station_id : str, optional
        Carried into messages and :attr:`StationInputs.station_id`.
    regional_skew, regional_skew_mse : float, optional
        As for ``Bulletin17C``. ``None`` is a station-skew analysis; a positive
        MSE weights the skew; a negative one is a generalized skew with
        ``MSE = -value`` (``emafitpr``'s encoding).
    low_outlier_threshold : float, optional
        A fixed PILF threshold in cfs; ``None`` runs MGBT.
    include_urban_regulated : bool, default False
        Keep code 6/C peaks.

    Returns
    -------
    StationInputs

    Raises
    ------
    ValueError
        An MSE of exactly zero (peakfq requires a positive skew standard
        error), a frame with no peaks, or any :func:`station_rows` error.
    """
    from .psf import PerceptionThreshold

    if "water_year" not in peaks.columns:
        raise ValueError("peak frame has no water_year column")
    years = pd.to_numeric(peaks["water_year"], errors="coerce").dropna()
    if years.empty:
        raise ValueError("peak frame has no water years")
    fields: Dict[str, str] = {}
    if regional_skew is None or regional_skew_mse is None:
        fields["SkewOpt"] = "Station"
    else:
        if regional_skew_mse == 0:
            raise ValueError(
                "regional_skew_mse of 0: peakfq requires a positive skew standard error"
            )
        fields["SkewOpt"] = "Weighted" if regional_skew_mse > 0 else "Generalized"
        fields["GenSkew"] = repr(float(regional_skew))
        fields["SkewSE"] = repr(float(np.sqrt(abs(regional_skew_mse))))
    if low_outlier_threshold is None:
        fields["LOType"] = "MGBT"
    else:
        fields["LOType"] = "FIXED"
        fields["LoThresh"] = repr(float(low_outlier_threshold))
    if include_urban_regulated:
        fields["Urb/Reg"] = "Yes"
    spec = StationSpec(
        station_id=station_id,
        thresholds=[PerceptionThreshold(int(years.min()), int(years.max()), 0.0, Q_MAX)],
        fields=fields,
    )
    frame = peaks.drop(columns="site_no") if "site_no" in peaks.columns else peaks
    inputs = convert_station(spec, frame)
    if regional_skew is not None and regional_skew_mse is not None:
        # Keep the caller's MSE bit-for-bit rather than sqrt()**2 of it.
        mse = float(regional_skew_mse)
        inputs = replace(inputs, regional_skew_mse=mse, _r_g_mse=mse)
    return inputs


def read_psf_peaks(psf: PsfFile, psf_path: Union[str, Path]) -> pd.DataFrame:
    """Read the peak file a ``.psf`` names, relative to the ``.psf``'s directory.

    Parameters
    ----------
    psf : PsfFile
        The parsed file.
    psf_path : str or Path
        Where the ``.psf`` lives.

    Returns
    -------
    pandas.DataFrame
        :func:`flowfreq.watstore.read_watstore`'s frame.

    Raises
    ------
    ValueError
        If the ``.psf`` has no ``I ASCI``/``I RDB`` line.
    NotImplementedError
        For ``I RDB``: read that file with your own RDB reader and pass the
        frame to :func:`convert_station`.
    FileNotFoundError
        If the named file does not exist.
    """
    spec = psf.peak_file
    if spec is None:
        raise ValueError("the .psf names no peak file (I ASCI or I RDB)")
    fmt, name = spec
    if fmt != "ASCI":
        raise NotImplementedError(
            f"I {fmt} peak files are not read here; load {name} yourself and pass "
            "the frame to convert_station"
        )
    base = Path(psf_path).parent
    path = base / name
    if not path.exists():
        # peakfq's own test data names "..._WATSTORE.txt" for a ".TXT" file;
        # on a case-sensitive filesystem only a case-insensitive match finds it.
        matches = [p for p in base.iterdir() if p.name.lower() == name.lower()]
        if not matches:
            raise FileNotFoundError(f"peak file {path} not found")
        path = matches[0]
    from .watstore import read_watstore

    return read_watstore(path)


def convert_psf(
    psf_path: Union[str, Path],
    peaks: Optional[pd.DataFrame] = None,
) -> Dict[str, StationInputs]:
    """Convert every station in a ``.psf`` file.

    Parameters
    ----------
    psf_path : str or Path
        The specification file.
    peaks : pandas.DataFrame, optional
        Peak data with a ``site_no`` column. Defaults to the file the ``.psf``
        names (:func:`read_psf_peaks`).

    Returns
    -------
    dict of str to StationInputs
        Keyed by station ID, in file order.

    Raises
    ------
    ValueError
        As :func:`convert_station`, for the first station that fails.
    """
    psf = read_psf(psf_path)
    if peaks is None:
        peaks = read_psf_peaks(psf, psf_path)
    return {sid: convert_station(spec, peaks, psf) for sid, spec in psf.stations.items()}
