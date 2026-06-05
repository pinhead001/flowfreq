"""
Report generation for validation results.

Produces text and JSON reports from comparison and benchmark results.
"""

from __future__ import annotations

import json
import logging
from typing import Any

from hydrolib.validation.comparisons import ComparisonResult

logger = logging.getLogger(__name__)


def generate_text_report(results: dict[str, ComparisonResult]) -> str:
    """Generate a plain-text validation report.

    Parameters
    ----------
    results : dict[str, ComparisonResult]
        Benchmark name to comparison result mapping.

    Returns
    -------
    str
        Formatted text report.
    """
    lines: list[str] = []
    lines.append("HydroLib Validation Report")
    lines.append("=" * 40)

    n_pass = sum(1 for r in results.values() if r.passed)
    n_total = len(results)
    lines.append(f"Overall: {n_pass}/{n_total} passed")
    lines.append("")

    for name, result in results.items():
        status = "PASS" if result.passed else "FAIL"
        lines.append(f"[{status}] {name}")
        lines.append(f"  Max diff: {result.max_diff_pct:.3f}%")
        lines.append(f"  Tolerance: {result.tolerance_pct}%")
        lines.append(f"  {result.summary}")

        if result.parameter_diffs:
            lines.append("  Parameters:")
            for param, diff in result.parameter_diffs.items():
                lines.append(f"    {param}: {diff:.4f}%")

        if result.quantile_diffs:
            lines.append("  Quantiles:")
            for aep, diff in sorted(result.quantile_diffs.items()):
                lines.append(f"    AEP {aep}: {diff:.4f}%")

        lines.append("")

    return "\n".join(lines)


def generate_json_report(results: dict[str, ComparisonResult]) -> str:
    """Generate a JSON validation report.

    Parameters
    ----------
    results : dict[str, ComparisonResult]
        Benchmark name to comparison result mapping.

    Returns
    -------
    str
        JSON string.
    """
    report: dict[str, Any] = {}

    for name, result in results.items():
        report[name] = {
            "passed": result.passed,
            "max_diff_pct": result.max_diff_pct,
            "tolerance_pct": result.tolerance_pct,
            "summary": result.summary,
            "parameter_diffs": result.parameter_diffs,
            "quantile_diffs": {str(k): v for k, v in result.quantile_diffs.items()},
            "ci_diffs": {str(k): v for k, v in result.ci_diffs.items()},
        }

    return json.dumps(report, indent=2)
