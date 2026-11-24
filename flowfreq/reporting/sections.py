"""Report section kinds: each computes its tables and figures from flowfreq, and nothing else.

A section builder takes the options a :class:`~flowfreq.reporting.spec.SectionSpec` carries
and returns the section's blocks, its references and its options *as resolved*: whatever the
builder fetched (a gage's annual peaks, a regional skew looked up by state) is written back
as plain values, so the same options build the same section again without the network.

The builders live one kind to a module (``gaged``, ``b17b``, ``regional``, ``weighted``,
``transposed``, ``engine_check``, ``comparison``, ``future``; ``daily``, ``regime_sections``
and ``summary`` hold the rest). This module is the registry, and every name that used to be
defined here is still importable from here. To replace one in a test, patch it in the module
that defines it.

Option names that look inconsistent between kinds name different things, and stored specs
use them, so they stay:

- ``target_name`` (``transposed``) is the ungaged site the flows are moved to; the donor
  gage keeps its own ``site_name`` inside the ``gage`` options.
- ``peaks_table`` (``gaged_ffa``) is the annual peaks; ``annual_table`` (``low_flow``,
  ``monthly_statistics``, ``baseflow``) is one row to a year of what the section computes.
- ``detail`` defaults on in :meth:`ReportSession.weight`, whose weighted section replaces
  the gage and regression sections, and off in :meth:`ReportSession.transpose`, whose donor
  section stays in the report.
"""

from __future__ import annotations

from typing import Any, Mapping

from ._base import (  # noqa: F401
    _AUTHORS_YEAR,
    _FIGURE_DPI,
    _INITIALS,
    _PEAK_KEYS,
    BULLETIN_17C,
    SECTION_KINDS,
    SectionContent,
    SectionKind,
    _build_text,
    _demote,
    _figure_png,
    _version,
    cfs,
    section_title,
    short_citation,
)
from .b17b import B17B_CAPTION, BULLETIN_17B, SUPERSEDED, _b17b_ffa, _build_b17b_ffa  # noqa: F401
from .comparison import (  # noqa: F401
    _COMPARED_KINDS,
    _ESTIMATE_KEYS,
    COMPARISON_CAPTION,
    _build_comparison,
    _compared,
    _comparison_label,
    _flow_table,
    _matched_aeps,
    _row_at,
)
from .engine_check import (  # noqa: F401
    _NO_EXTENSION,
    _REFERENCE_KEYS,
    _REFERENCE_PARAMETERS,
    _REFERENCE_QUANTILE_KEYS,
    _TOLERANCES,
    _build_engine_check,
    _difference,
    _extension_digest,
    _finite,
    _live_reference,
    _record_digest,
    _reference_result,
)
from .future import (  # noqa: F401
    _FACTOR_SET_KEYS,
    FUTURE_FLOW_CAPTION,
    _build_future_flow,
    _factor_set,
)
from .gaged import (  # noqa: F401
    _EXTENSION_KEYS,
    _INDEX_KEYS,
    _build_gaged_ffa,
    _extended_frame,
    _extension_blocks,
    _gaged_ffa,
    _peak_frame,
    _resolve_skew,
    _skew_sentence,
    _threshold_rows,
    download_peaks,
)
from .regional import (  # noqa: F401
    _STREAMSTATS_KEYS,
    _UNGAGED_HEADERS,
    _build_ungaged,
    _freeze_streamstats,
    _same_quantiles,
    _streamstats_notes,
    _streamstats_record,
    _ungaged,
    _used_codes,
)
from .spec import ReportError
from .transposed import (  # noqa: F401
    _EXPONENT_KEYS,
    _build_transposed,
    _positive,
    _regression_exponents,
    _site_class,
    _transposed,
)
from .weighted import _build_weighted, _nested  # noqa: F401

SECTION_KINDS.update(
    (kind.name, kind)
    for kind in (
        SectionKind(
            name="text",
            summary="A titled section of narrative only.",
            required=(),
            optional=(),
            build=_build_text,
        ),
        SectionKind(
            name="gaged_ffa",
            summary=(
                "Bulletin 17C flood frequency at a gage: record, method, fitted parameters, "
                "quantiles with confidence limits, frequency curve, stationarity screen."
            ),
            required=(),
            optional=(
                "site_no",
                "site_name",
                "drainage_area_sqmi",
                "peaks",
                "regional_skew",
                "regional_skew_se",
                "regional_skew_citation",
                "regional_skew_note",
                "regional_skew_advisory",
                "regional_skew_state",
                "regional_skew_region",
                "regional_skew_from_site",
                "allow_provisional_skew",
                "station_skew_only",
                "use_default_skew",
                "perception_thresholds",
                "low_outlier_threshold",
                "apply_peak_codes",
                "allow_regulated",
                "stationarity",
                "peaks_table",
                "extension",
            ),
            build=_build_gaged_ffa,
        ),
        SectionKind(
            name="b17b_ffa",
            summary=(
                "Bulletin 17B (superseded) flood frequency at a gage, for comparison with "
                "Bulletin 17C. Computed only when this section is requested."
            ),
            required=(),
            optional=(
                "site_no",
                "site_name",
                "drainage_area_sqmi",
                "peaks",
                "regional_skew",
                "regional_skew_se",
                "regional_skew_citation",
                "regional_skew_note",
                "regional_skew_advisory",
                "regional_skew_state",
                "regional_skew_region",
                "regional_skew_from_site",
                "allow_provisional_skew",
                "generalized_skew",
                "generalized_skew_mse",
                "station_skew_only",
                "historic_period",
                "low_outlier_threshold",
                "perception_thresholds",
                "apply_peak_codes",
                "allow_regulated",
                "hand_computation",
                "aeps",
            ),
            build=_build_b17b_ffa,
            title_label=SUPERSEDED,
        ),
        SectionKind(
            name="ungaged",
            summary=(
                "Regional regression estimate at an ungaged site: regions, inputs, quantiles "
                "with prediction limits, caveats. The basin characteristics and region are "
                "the analyst's, or come from StreamStats when none are given and are then "
                "stored in the specification."
            ),
            required=("state", "lat", "lon"),
            optional=(
                "characteristics",
                "region",
                "streamstats",
                "all_characteristics",
                "site_name",
                "aeps",
                "allow_nss_snapshot",
                "allow_provisional",
                "allow_extrapolation",
                "allow_inaccurate",
                "weighting_space",
                "min_coverage_percent",
                "supplied_characteristics",
            ),
            build=_build_ungaged,
        ),
        SectionKind(
            name="weighted",
            summary=(
                "A gage's Bulletin 17C quantiles weighted with the regional regression at the "
                "gage (Bulletin 17C Appendix 9): the weighted table, then both analyses."
            ),
            required=("gage", "regression"),
            optional=("detail",),
            build=_build_weighted,
        ),
        SectionKind(
            name="transposed",
            summary=(
                "A gage's Bulletin 17C quantiles transposed to an ungaged site on the same "
                "stream by the drainage-area ratio, with a cited regression exponent: the "
                "arithmetic, the caveats, then the analysis at the gage."
            ),
            required=("donor", "target_area_sqmi", "exponents"),
            optional=(
                "donor_area_sqmi",
                "target_name",
                "max_aep",
                "allow_out_of_range",
                "donor_site_class",
                "target_site_class",
                "allow_regulation_mismatch",
                "detail",
            ),
            build=_build_transposed,
        ),
        SectionKind(
            name="engine_check",
            summary=(
                "A gage's Bulletin 17C fit checked against a reference computation of the "
                "same record: USGS peakfq 8.1.0 through flowfreq's Fortran extension, whose "
                "numbers are then stored in the specification, or numbers the analyst gives. "
                "Parameters, quantiles and confidence limits side by side, with the "
                "differences and whether each is within tolerance."
            ),
            required=("gage",),
            optional=(
                "reference",
                "tolerance_pct",
                "parameter_tolerance_pct",
                "ci_tolerance_pct",
                "detail",
            ),
            build=_build_engine_check,
        ),
        SectionKind(
            name="future_flow",
            summary=(
                "Flood quantiles under a future condition: a gage's Bulletin 17C quantiles or "
                "a regression estimate multiplied by change factors the analyst supplies and "
                "cites, with which factors were published, interpolated or held."
            ),
            required=("factor_set",),
            optional=("gage", "regression", "aeps", "detail"),
            build=_build_future_flow,
        ),
        SectionKind(
            name="comparison",
            summary=(
                "Several estimates of the flood quantiles at one site, side by side: each a "
                "gage's Bulletin 17C fit, a transposition from a gage, a regional "
                "regression or, when requested, a Bulletin 17B (superseded) fit, with its "
                "limits, the ratio of the largest to the smallest, and each estimate's own "
                "analysis below. With a Bulletin 17B fit it adds a plot of the curves."
            ),
            required=("estimates",),
            optional=("site_name", "aeps", "detail"),
            build=_build_comparison,
            title_label=_comparison_label,
        ),
    )
)


def build_section(kind: str, options: Mapping[str, Any]) -> SectionContent:
    """Build one section.

    Parameters
    ----------
    kind : str
        A key of :data:`SECTION_KINDS`.
    options : mapping
        The kind's options (see ``SECTION_KINDS[kind].required`` and ``.optional``).

    Returns
    -------
    SectionContent

    Raises
    ------
    ReportError
        An unknown kind or option, a missing option, or an analysis that failed. Refusals
        the analyses raise themselves (a regulated record, an input outside an equation's
        range) pass through unchanged.
    """
    if kind not in SECTION_KINDS:
        raise ReportError(
            f"unknown section kind {kind!r}; choose from {', '.join(sorted(SECTION_KINDS))}"
        )
    section = SECTION_KINDS[kind]
    section.check(options)
    return section.build(dict(options))


__all__ = [
    "B17B_CAPTION",
    "BULLETIN_17B",
    "BULLETIN_17C",
    "COMPARISON_CAPTION",
    "FUTURE_FLOW_CAPTION",
    "SECTION_KINDS",
    "SectionContent",
    "SectionKind",
    "build_section",
    "cfs",
    "download_peaks",
    "short_citation",
]
