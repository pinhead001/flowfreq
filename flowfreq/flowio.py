"""
flowfreq.flowio - Reading and writing flow time series

Format helpers shared by the retrieval and analysis modules. Kept apart from
:mod:`flowfreq.usgs` because they apply to any flow series, not only to one
fetched from NWIS.
"""

from __future__ import annotations

from pathlib import Path
from typing import Union

import pandas as pd

#: Compression used for Parquet output unless overridden. On a year of
#: 15-minute discharge, zstd lands a few percent under the Parquet default
#: (snappy) and roughly a third the size of the equivalent CSV. Every standard
#: Parquet reader -- pyarrow, R's arrow, DuckDB, parquet-mr -- decompresses it
#: without configuration, so the default costs no portability.
DEFAULT_PARQUET_COMPRESSION: str = "zstd"


def save_flow_frame(
    df: pd.DataFrame,
    path: Union[str, Path],
    compression: str = DEFAULT_PARQUET_COMPRESSION,
) -> Path:
    """Write a flow time series to disk, choosing the writer by file extension.

    Parquet is the format to prefer, especially for instantaneous series. It
    stores the timezone-aware index and float columns with their dtypes intact,
    where CSV round-trips both through text and needs them re-inferred on read.
    It is typically several times smaller, and it is readable from R, DuckDB and
    most GIS toolchains, which matters when the numbers travel further than the
    Python session that produced them.

    CSV remains available for cases where a human has to read the file, at the
    cost of the index dtype.

    Parameters
    ----------
    df : pd.DataFrame
        Frame to write, typically from
        :meth:`flowfreq.usgs.USGSgage.download_instantaneous_flow` or
        :meth:`flowfreq.usgs.USGSgage.download_daily_flow`.
    path : str or Path
        Destination. ``.parquet``/``.pq`` writes Parquet; ``.csv`` and
        ``.csv.gz`` write CSV.
    compression : str
        Parquet compression codec, default ``"zstd"``. Accepts anything
        pyarrow supports (``"snappy"``, ``"gzip"``, ``"brotli"``, ``"lz4"``,
        ``"none"``). Ignored for CSV, whose compression follows the extension.

    Returns
    -------
    Path
        The path written.

    Raises
    ------
    ValueError
        The extension is not one of the supported formats.
    ImportError
        The Parquet engine is unavailable. ``pyarrow`` is a flowfreq
        dependency, so this indicates a damaged environment rather than a
        missing optional extra.

    Examples
    --------
    >>> save_flow_frame(iv, "methow_iv.parquet")           # doctest: +SKIP
    >>> save_flow_frame(iv, "methow_iv.csv")               # doctest: +SKIP
    """
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix in (".parquet", ".pq"):
        try:
            df.to_parquet(path, compression=compression)
        except ImportError as exc:  # pragma: no cover - requires a broken install
            raise ImportError(
                f"Writing {path.name} needs the pyarrow Parquet engine. pyarrow is a "
                f"flowfreq dependency, so this environment looks incomplete: try "
                f"`pip install --force-reinstall pyarrow`, or save as .csv to work "
                f"around it (CSV does not preserve the tz-aware index dtype)."
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

    Parquet restores the frame exactly as written, compression codec included:
    dtypes, column order and the timezone-aware index all come back unchanged.
    CSV is re-parsed, so the index is reconstructed from text and a
    timezone-aware index returns normalized to UTC rather than to the offset it
    was written in -- the instants are identical, the printed offset may not be.

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
