"""Wave 1 and 2 vignette: one ungaged-site peak-flow estimate per state, offline.

Roadmap §3.2 asks each wave for "a wave vignette with one worked ungaged-site example
per state". Each example below is a real basin on or near a state line, treated as
ungaged: its basin characteristics and its regression-region area weights were captured
from StreamStats by ``tools/wave_qa.py vignette`` (``wave1_wave2_sites.json`` beside this
file). From there everything is offline:

1. each region's equation is evaluated with :func:`flowfreq.regression.evaluate`, which
   gives a 90 % prediction interval (site-specific where the report publishes a
   covariance matrix, else from the average SEP);
2. a basin that spans several regions is area-weighted with
   :func:`flowfreq.regression.evaluate_weighted` (flows, as NSS's ``areaave`` does);
   Oregon's Region 2A/2B elevation blend uses
   :func:`flowfreq.regression.oregon.estimate_region2`;
3. the gage's own Bulletin 17C 1 % AEP flow is printed alongside, as a check.

Run::

    python docs/vignettes/wave1_wave2_ungaged.py
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Tuple

from flowfreq.regression import RegressionEstimate, evaluate, evaluate_weighted, load_state
from flowfreq.regression.oregon import estimate_region2

SITES = Path(__file__).with_name("wave1_wave2_sites.json")
AEP = 0.01


def estimate(state: str, site: dict) -> Tuple[float, List[Tuple[str, float, RegressionEstimate]]]:
    """Area-weighted 1 % AEP flow and each region's estimate (with its interval)."""
    lib = load_state(state)
    chars: Dict[str, float] = site["characteristics"]
    parts: List[Tuple[str, float, RegressionEstimate]] = []
    for region, weight in site["region_weights"].items():
        if state == "OR" and region == "2":
            r2 = estimate_region2(AEP, chars, lib=lib)
            for sub, est in r2.estimates.items():
                parts.append((sub, weight * r2.weights[sub], est))
            continue
        parts.append((region, weight, evaluate(lib.equation(region, AEP), chars)))
    flow = evaluate_weighted([(w, e) for _, w, e in parts], space="linear")
    return flow, parts


def main() -> None:
    sites = json.loads(SITES.read_text(encoding="utf-8"))
    for state, site in sites.items():
        flow, parts = estimate(state, site)
        print(f"\n{state}  {site['site_no']} {site['site_name']}")
        print(
            f"    outlet {site['outlet']}, StreamStats {site['streamstats_region']} "
            f"(captured {site['captured']})"
        )
        print("    " + ", ".join(f"{k}={v:g}" for k, v in site["characteristics"].items()))
        for region, weight, est in parts:
            lo, hi = est.interval or (float("nan"), float("nan"))
            print(
                f"    region {region:>6} weight {weight:4.2f}: Q1% = {est.flow_cfs:9,.0f} cfs, "
                f"90% PI {lo:9,.0f} - {hi:9,.0f} ({est.interval_method})"
            )
        print(
            f"    area-weighted Q1% = {flow:,.0f} cfs;  gage B17C Q1% = "
            f"{site['b17c_q1pct']:,.0f} cfs ({site['b17c_n']} peaks)"
        )


if __name__ == "__main__":
    main()
