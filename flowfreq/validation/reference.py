"""Reference results for validating the native EMA against peakfq 8.1.0.

This replaces ``flowfreq.peakfqsa``, which was written to drive a standalone
PeakfqSA executable. That executable does not exist and never did -- the
premise came from the original build brief (see
``AGENT_BUILD_INSTRUCTIONS_Claude.md`` §7) and was corrected later. The
subprocess wrapper, the specification-file writer and the ``.out`` parser were
therefore never run against anything; only mocks. They have been removed.

What survives is the part the validation layer actually used: a container for
"what the reference implementation says", now pointed at references that exist.
Two of them:

``from_golden``
    The committed golden file, generated from peakfq 8.1.0 by
    ``tools/gen_fortran_golden.py``. Always available, no toolchain needed.

``from_emafit``
    A live call through the f2py bridge in :mod:`flowfreq.peakfqr`. Needs the
    extension built (``python build_fortran/build.py``), and is the only way to
    get a reference for a site no golden file covers.

Both return the same shape, so :class:`~flowfreq.validation.comparisons.FrequencyComparator`
does not care which produced it.

All flow quantities from the Fortran are log10(cfs); everything this module
returns is in real space (cfs), which is what the native results use.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple, Union

import numpy as np

logger = logging.getLogger(__name__)

__all__ = ["ReferenceResult", "cmoms_to_parameters"]

# peakfqr/R/main.R: the sentinel bounds a flow interval is censored against.
QMIN = 1e-20
QMAX = 1e20


def cmoms_to_parameters(
    cmoms: Sequence[Sequence[float]],
    skew_mse: Optional[float] = None,
    skew_mse_systematic: Optional[float] = None,
    pseudo_record_length: Optional[float] = None,
    weight_factor: Optional[float] = None,
) -> Dict[str, float]:
    """Map ``emafitpr``'s ``cmoms(3,3)`` onto FlowFreq parameter names.

    Column 1 is the fit using regional information, column 2 the at-site-only
    fit, column 3 the B17B MSE formula applied at-site. Rows are mean, variance
    and skew. See ``vendor/peakfqr/R/fortranWrappers.R``.

    Parameters
    ----------
    cmoms : sequence of sequence of float
        The 3x3 central-moment matrix, row-major as ``emafitpr`` returns it.
    skew_mse, skew_mse_systematic, pseudo_record_length, weight_factor : float, optional
        ``as_G_mse_o``, ``as_G_mse_Syst_o``, ``as_G_PRL_o`` and ``Wdout``.
        Included when given; these have no native counterpart yet but are what
        the two open parity defects are argued from.

    Returns
    -------
    dict of str to float
        Parameter names matching ``Bulletin17C.to_comparison_dict()``.
    """
    m = np.asarray(cmoms, dtype=float)
    if m.shape != (3, 3):
        raise ValueError(f"cmoms must be 3x3, got {m.shape}")

    params = {
        "mean_log": float(m[0][0]),
        "std_log": float(np.sqrt(m[1][0])),
        "skew_weighted": float(m[2][0]),
        "mean_log_at_site": float(m[0][1]),
        "std_log_at_site": float(np.sqrt(m[1][1])),
        "skew_at_site": float(m[2][1]),
    }
    optional = {
        "mse_skew": skew_mse,
        "mse_skew_systematic": skew_mse_systematic,
        "pseudo_record_length": pseudo_record_length,
        "weight_factor": weight_factor,
    }
    params.update({k: float(v) for k, v in optional.items() if v is not None})
    return params


@dataclass
class ReferenceResult:
    """What the reference implementation says, in real space.

    Parameters
    ----------
    source : str
        Provenance, carried into comparison reports so a passing benchmark
        names what it passed against.
    station_name, begyear, endyear : str, int, int
        Site identification, when the source records it.
    n_peaks, n_systematic, n_historical : int
        Record composition.
    low_outlier_count : int
        PILFs detected (``gbnlow``).
    low_outlier_threshold : float
        PILF threshold in cfs (``10 ** gbval``), 0.0 when MGBT found none.
    parameters : dict of str to float
        LP3 parameters; see :func:`cmoms_to_parameters`.
    quantiles : dict of float to float
        Annual exceedance probability to discharge (cfs).
    confidence_intervals : dict of float to tuple of float
        AEP to ``(lower, upper)`` in cfs.
    variance : dict of float to float
        AEP to ``var_est``, the variance of the log10 quantile estimate.
    """

    source: str = ""
    station_name: str = ""
    begyear: int = 0
    endyear: int = 0
    n_peaks: int = 0
    n_systematic: int = 0
    n_historical: int = 0
    low_outlier_count: int = 0
    low_outlier_threshold: float = 0.0
    parameters: Dict[str, float] = field(default_factory=dict)
    quantiles: Dict[float, float] = field(default_factory=dict)
    confidence_intervals: Dict[float, Tuple[float, float]] = field(default_factory=dict)
    variance: Dict[float, float] = field(default_factory=dict)

    # ------------------------------------------------------------------ #
    # Constructors
    # ------------------------------------------------------------------ #

    @classmethod
    def from_golden(cls, path: Union[str, Path]) -> "ReferenceResult":
        """Load a golden file written by ``tools/gen_fortran_golden.py``.

        Parameters
        ----------
        path : str or Path
            Path to the golden JSON.

        Returns
        -------
        ReferenceResult

        Raises
        ------
        FileNotFoundError
            If *path* does not exist.
        KeyError
            If the file is missing a section, which means it predates the
            current generator and should be regenerated.
        """
        path = Path(path)
        golden = json.loads(path.read_text())
        meta, inputs, outputs = golden["meta"], golden["inputs"], golden["outputs"]
        version = meta.get("peakfq_version", "unknown")
        return cls._assemble(
            source=f"peakfq {version} golden file ({path.name})",
            station_name=meta.get("description", ""),
            aeps=inputs["aeps"],
            outputs=outputs,
            n_peaks=int(outputs["n"]),
            n_historical=int(np.count_nonzero(np.asarray(inputs["dtype"]) == 1)),
        )

    @classmethod
    def from_emafit(
        cls,
        ql: Sequence[float],
        qu: Sequence[float],
        tl: Sequence[float],
        tu: Sequence[float],
        dtype: Sequence[int],
        aeps: Sequence[float],
        regional_skew: float,
        regional_skew_mse: float,
        *,
        station_name: str = "",
        eps: float = 0.90,
        weight_opt: int = 1,
        gbthrsh0: float = -99.0,
        reg_m: float = 0.0,
        reg_m_mse: float = -1e99,
        reg_sd: float = 1.0,
        reg_sd_mse: float = 1e99,
    ) -> "ReferenceResult":
        """Run the vendored Fortran through the f2py bridge.

        Every flow argument is log10, exactly as ``emafitpr`` expects; the
        defaults for the non-flow arguments are the ones
        ``vendor/peakfqr/R/fortranWrappers.R`` passes.

        Parameters
        ----------
        ql, qu, tl, tu : sequence of float
            Flow interval and perception threshold bounds, log10.
        dtype : sequence of int
            0 for systematic, 1 for a peak carrying the USGS historic flag.
        aeps : sequence of float
            Annual exceedance probabilities to evaluate.
        regional_skew, regional_skew_mse : float
            ``r_G`` and ``r_G_mse``; see the MSE encoding in ``CLAUDE.md``.
        station_name : str, optional
            Carried through to the result for reporting.
        eps : float, optional
            Confidence-interval coverage, 0.90 for a 90% interval.
        weight_opt : int, optional
            1=HWN, 2=ERL, 3=INV.
        gbthrsh0 : float, optional
            ``<= -6`` computes MGBT; a larger value is used as the threshold.
        reg_m, reg_m_mse, reg_sd, reg_sd_mse : float, optional
            Regional mean and standard deviation with their MSEs.

        Returns
        -------
        ReferenceResult

        Raises
        ------
        ImportError
            If the f2py extension is not built. Build it with
            ``python build_fortran/build.py`` (needs gfortran and meson).
        """
        try:
            from flowfreq.peakfqr import emafitpr
        except ImportError as exc:  # pragma: no cover - depends on the build
            raise ImportError(
                "the Fortran bridge is not built; run python build_fortran/build.py "
                "(needs gfortran and meson), or use ReferenceResult.from_golden()"
            ) from exc

        aeps = np.asarray(aeps, dtype=float)
        out = emafitpr(
            np.asarray(ql, dtype=float),
            np.asarray(qu, dtype=float),
            np.asarray(tl, dtype=float),
            np.asarray(tu, dtype=float),
            np.asarray(dtype, dtype=np.int32),
            reg_m,
            reg_m_mse,
            reg_sd,
            reg_sd_mse,
            regional_skew,
            regional_skew_mse,
            gbthrsh0,
            1.0 - aeps,  # emafitpr takes non-exceedance probabilities
            eps,
            weight_opt,
        )
        (
            gbval,
            _gbns,
            _gbnzero,
            gbnlow,
            _gbp,
            _gbqs,
            as_g_mse,
            as_g_mse_syst,
            as_g_prl,
            cmoms,
            yp,
            ci_low,
            ci_high,
            var_est,
            wdout,
            *_rest,
        ) = out

        n_historical = int(np.count_nonzero(np.asarray(dtype) == 1))
        return cls._assemble(
            source="live emafitpr via flowfreq.peakfqr",
            station_name=station_name,
            aeps=aeps,
            outputs={
                "mgbt": {"gbval": float(gbval), "gbnlow": int(gbnlow)},
                "skew": {
                    "as_G_mse_o": float(as_g_mse),
                    "as_G_mse_Syst_o": float(as_g_mse_syst),
                    "as_G_PRL_o": float(as_g_prl),
                    "Wdout": float(wdout),
                },
                "cmoms": np.asarray(cmoms, dtype=float).tolist(),
                "quantiles": {
                    "yp": np.asarray(yp, dtype=float).tolist(),
                    "ci_low": np.asarray(ci_low, dtype=float).tolist(),
                    "ci_high": np.asarray(ci_high, dtype=float).tolist(),
                    "var_est": np.asarray(var_est, dtype=float).tolist(),
                },
            },
            n_peaks=len(ql),
            n_historical=n_historical,
        )

    # ------------------------------------------------------------------ #
    # Internals
    # ------------------------------------------------------------------ #

    @classmethod
    def _assemble(
        cls,
        *,
        source: str,
        station_name: str,
        aeps: Sequence[float],
        outputs: Mapping[str, Any],
        n_peaks: int,
        n_historical: int = 0,
    ) -> "ReferenceResult":
        """Build a result from ``emafitpr`` outputs, converting log10 to cfs."""
        aeps = [float(a) for a in aeps]
        q = outputs["quantiles"]
        mgbt = outputs.get("mgbt", {})
        skew = outputs.get("skew", {})

        gbval = float(mgbt.get("gbval", -99.0))
        # gbthrsh0 = -99 is the "run MGBT" sentinel; it comes straight back out
        # when nothing was flagged, and 10 ** -99 is not a threshold.
        threshold = 10.0**gbval if gbval > -6.0 else 0.0

        return cls(
            source=source,
            station_name=station_name,
            n_peaks=int(n_peaks),
            n_systematic=int(n_peaks) - int(n_historical),
            n_historical=int(n_historical),
            low_outlier_count=int(mgbt.get("gbnlow", 0)),
            low_outlier_threshold=threshold,
            parameters=cmoms_to_parameters(
                outputs["cmoms"],
                skew_mse=skew.get("as_G_mse_o"),
                skew_mse_systematic=skew.get("as_G_mse_Syst_o"),
                pseudo_record_length=skew.get("as_G_PRL_o"),
                weight_factor=skew.get("Wdout"),
            ),
            quantiles={aep: float(10.0 ** q["yp"][i]) for i, aep in enumerate(aeps)},
            confidence_intervals={
                aep: (float(10.0 ** q["ci_low"][i]), float(10.0 ** q["ci_high"][i]))
                for i, aep in enumerate(aeps)
            },
            variance={aep: float(q["var_est"][i]) for i, aep in enumerate(aeps)},
        )
