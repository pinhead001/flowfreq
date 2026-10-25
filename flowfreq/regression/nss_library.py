"""Build an ``nss_snapshot`` equation library from an NSS snapshot.

``tools/snapshot_nss.py`` harvests a jurisdiction's peak-flow equations from the
USGS NSS web service into ``data/nss_snapshots/<CODE>_<YYYY-MM-DD>.json``.
:func:`snapshot_to_library` turns one such snapshot into the JSON form of a
:class:`~flowfreq.regression.library.StateLibrary` with status ``nss_snapshot``,
which :func:`~flowfreq.regression.library.load_state` returns only on opt-in
(``allow_nss_snapshot=True``).

Every equation is re-parsed from its NSS string with
:func:`flowfreq.regression.nss.parse_equation` rather than read from the
snapshot's stored ``parsed`` record, so a snapshot taken before a parser
extension still gains it. A string the parser refuses, and a region NSS would
not resolve, are listed in the library's ``unparsed`` entries with the reason;
neither is dropped silently.
"""

from __future__ import annotations

import dataclasses
import logging
import math
import re
from typing import Any, Dict, List, Mapping, Optional, Tuple, Union

from flowfreq.regression.equations import Citation
from flowfreq.regression.library import SCHEMA_VERSION
from flowfreq.regression.nss import (
    NSSEquationError,
    ParsedConstant,
    ParsedEquation,
    parse_constant_equation,
    parse_equation,
    to_regression_equation,
)

logger = logging.getLogger(__name__)

_YEAR = re.compile(r"\s*(\d{4})")


def asep_to_sep_log(asep_pct: float) -> float:
    """Convert an average standard error of prediction in percent to log10 units.

    ``sqrt(ln(1 + (asep/100)**2)) / ln(10)``: the conversion WY.json documents
    for a printed SE_P in percent.

    Parameters
    ----------
    asep_pct : float
        NSS's ``ASEp`` error, in percent.

    Returns
    -------
    float
        The standard error of prediction in log10 units.

    Raises
    ------
    ValueError
        If ``asep_pct`` is not finite and positive.
    """
    if not (math.isfinite(asep_pct) and asep_pct > 0):
        raise ValueError(f"ASEp must be finite and positive, got {asep_pct}")
    return math.sqrt(math.log(1.0 + (asep_pct / 100.0) ** 2)) / math.log(10.0)


def _publication(record: Mapping[str, Any]) -> str:
    return f"{str(record.get('author') or '').strip()} {str(record.get('title') or '').strip()}"


def _citation(region: Mapping[str, Any], citations: Mapping[str, Any]) -> Citation:
    cid = region.get("citation_id")
    record = citations.get(str(cid))
    if not record:
        raise ValueError(
            f"region {region.get('region_code')}: NSS citation {cid} is not in the snapshot"
        )
    title = str(record.get("title") or "")
    year = _YEAR.match(title)
    return Citation(
        publication=_publication(record),
        table=f"NSS region {region.get('region_id')} ({region.get('region_code')}), "
        f"citation {cid}",
        year=int(year.group(1)) if year else None,
        url=str(record.get("citation_url") or ""),
    )


def _unparsed(
    region: Mapping[str, Any],
    statistic_code: str,
    aep: Optional[float],
    equation: str,
    reason: str,
) -> Dict[str, Any]:
    return {
        "region_code": str(region.get("region_code") or ""),
        "region_name": str(region.get("region_name") or ""),
        "statistic_code": statistic_code,
        "aep": aep,
        "equation": equation,
        "reason": reason,
    }


def _parse(text: str) -> Union[ParsedEquation, ParsedConstant]:
    """Parse an NSS string, reading NSS's form of a printed ``Q=0`` as a constant."""
    try:
        return parse_equation(text)
    except NSSEquationError as exc:
        try:
            return parse_constant_equation(text)
        except NSSEquationError:
            raise exc from None


def _plural(n: int, word: str) -> str:
    return f"{n} {word}" + ("" if n == 1 else "s")


def snapshot_to_library(snapshot: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """Build the JSON form of an ``nss_snapshot`` library from an NSS snapshot.

    Parameters
    ----------
    snapshot : dict
        A snapshot written by ``tools/snapshot_nss.py`` (schema
        ``flowfreq.nss_snapshot/1``).

    Returns
    -------
    dict or None
        The library (``schema_version`` 1, ``status`` ``"nss_snapshot"``), for
        :func:`~flowfreq.regression.library.parse_state` or writing to
        ``flowfreq/data/regression/nss_snapshot/<CODE>.json``; ``None`` when the
        snapshot has no regions (NSS serves no peak-flow equations there).

    Notes
    -----
    Each equation cites its NSS citation record and region, takes its limits and
    units from the region's parameters, and gets ``sep_log`` from NSS's ASEp
    (:func:`asep_to_sep_log`) where that is present and positive; otherwise
    ``sep_log`` stays unset. NSS's raw ``nss_sep`` is not stored: it is not a log10
    SEP. A statistic with no AEP (not a peak-flow AEP or recurrence-interval code)
    is skipped and counted in ``notes``.

    Raises
    ------
    ValueError
        If a region cites an NSS citation the snapshot does not carry.
    """
    regions: List[Dict[str, Any]] = list(snapshot.get("regions") or [])
    if not regions:
        return None
    state = str(snapshot["state"])
    retrieved = str(snapshot.get("retrieved_utc") or "")
    citations: Mapping[str, Any] = snapshot.get("citations") or {}

    equations: List[Dict[str, Any]] = []
    unparsed: List[Dict[str, Any]] = []
    cited: Dict[str, None] = {}
    harvested = skipped = 0
    for region in regions:
        if region.get("citation_id") is not None:
            cited[str(region["citation_id"])] = None
        if region.get("status") != "resolved":
            reason = str(region.get("unresolved_reason") or "no reason recorded")
            unparsed.append(_unparsed(region, "", None, "", f"unresolved: {reason}"))
            continue
        citation = _citation(region, citations)
        params = region.get("parameters") or []
        limits: Dict[str, Tuple[Optional[float], Optional[float]]] = {
            str(p["code"]): (p.get("min"), p.get("max")) for p in params
        }
        units = {str(p["code"]): str(p.get("units") or "") for p in params}
        for rec in region.get("equations") or []:
            harvested += 1
            code = str(rec.get("statistic_code") or "")
            aep = rec.get("aep")
            if aep is None:
                skipped += 1
                continue
            text = str(rec.get("equation") or "")
            try:
                parsed = _parse(text)
            except NSSEquationError as exc:
                unparsed.append(_unparsed(region, code, float(aep), text, str(exc)))
                continue
            eq = to_regression_equation(
                parsed,
                region_code=str(region["region_code"]),
                region_name=str(region.get("region_name") or ""),
                aep=float(aep),
                citation=citation,
                limits=limits,
                units=units,
            )
            asep = rec.get("asep_pct")
            if eq.form == "loglinear" and asep is not None and float(asep) > 0:
                eq = dataclasses.replace(eq, sep_log=asep_to_sep_log(float(asep)))
            equations.append(eq.to_dict())

    source_reports = [
        {
            "citation": _publication(citations[cid]),
            "url": str(citations[cid].get("citation_url") or ""),
            "nss": {
                "citation_id": int(cid),
                "last_year_of_data": citations[cid].get("last_year_of_data"),
            },
        }
        for cid in cited
        if cid in citations
    ]
    snapshot_file = f"data/nss_snapshots/{state}_{retrieved[:10]}.json"
    notes = (
        f"Generated by flowfreq.regression.nss_library.snapshot_to_library from "
        f"{snapshot_file} (NSS retrieved {retrieved or 'at an unrecorded time'}): "
        f"{_plural(harvested, 'equation')} harvested from {_plural(len(regions), 'region')}, "
        f"{_plural(len(equations), 'equation')} stored, {len(unparsed)} unparsed equations or regions, "
        f"{_plural(skipped, 'non-AEP statistic')} skipped."
    )
    logger.debug("%s: %s", state, notes)
    return {
        "schema_version": SCHEMA_VERSION,
        "state": state,
        "status": "nss_snapshot",
        "retrieved_utc": retrieved,
        "issue": None,
        "source_reports": source_reports,
        "notes": notes,
        "equations": equations,
        "unparsed": unparsed,
    }
