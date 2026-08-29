"""
hydrolib.flowio - Reading and writing flow time series

Format helpers shared by the retrieval and analysis modules. Kept apart from
:mod:`hydrolib.usgs` because they apply to any flow series, not only to one
fetched from NWIS.
"""

from __future__ import annotations

from pathlib import Path
from typing import Union

import pandas as pd


def save_flow_frame(df: pd.DataFrame, path: Union[str, Path]) -> Path:
    """Write a flow time series to disk, choosing the writer by file extension.

    Parquet is the format worth preferring for instantaneous series. It stores
    the timezone-aware index and float column with their dtypes intact, where
    CSV round-trips both through text and needs them re-inferred on read, and
    it is typically several times smaller for a series of this length. It is
    also readable from R and most GIS toolchains, which matters when the
    numbers travel further than the Python session that produced them.

    Parquet needs ``pyarrow``, which hydrolib does not require, so CSV remains
    available with no extra dependency and this function says plainly which one
    is missing rather than failing obscurely.

    Parameters
    ----------
    df : pd.DataFrame
        Frame to write, typically from
        :meth:`USGSgage.download_instantaneous_flow` or
        :meth:`USGSgage.download_daily_flow`.
    path : str or Path
        Destination. ``.parquet``/``.pq`` writes Parquet; ``.csv`` and
        ``.csv.gz`` write CSV.

    Returns
    -------
    Path
        The path written.

    Raises
    ------
    ValueError
        The extension is not one of the supported formats.
    ImportError
        Parquet was requested but no Parquet engine is installed.
    """
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix in (".parquet", ".pq"):
        try:
            df.to_parquet(path)
        except ImportError as exc:
            raise ImportError(
                f"Writing {path.name} needs a Parquet engine, which hydrolib does not "
                f"install by default. Either `pip install pyarrow` or save as .csv "
                f"(note that CSV does not preserve the tz-aware index dtype)."
            ) from exc
        return path

    if suffix == ".csv" or path.name.lower().endswith(".csv.gz"):
        df.to_csv(path)
        return path

    raise ValueError(
        f"Unsupported flow-frame format {suffix!r}; use .parquet, .pq, .csv, or .csv.gz"
    )


def load_flow_frame(path: Union[str, Path]) -> pd.DataFrame:
    """Read a flow time series written by :func:`save_flow_frame`.

    Parquet restores the frame exactly as written. CSV is re-parsed, so the
    index is reconstructed from text and a timezone-aware index comes back
    normalized to UTC rather than to the offset it was written in — the instants
    are identical, the printed offset may not be.

    Parameters
    ----------
    path : str or Path
        File to read; the format is taken from the extension.

    Returns
    -------
    pd.DataFrame
        Frame indexed by datetime.

    Raises
    ------
    ValueError
        The extension is not one of the supported formats.
    """
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix in (".parquet", ".pq"):
        return pd.read_parquet(path)

    if suffix == ".csv" or path.name.lower().endswith(".csv.gz"):
        return pd.read_csv(path, index_col=0, parse_dates=[0])

    raise ValueError(
        f"Unsupported flow-frame format {suffix!r}; use .parquet, .pq, .csv, or .csv.gz"
    )
