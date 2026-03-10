"""
Sounding analysis package.

Public API re-exported here so that ``from sounding import X``
continues to work exactly as before.
"""

# ── Constants ───────────────────────────────────────────────────────
from .constants import (
    STATIONS,
    STATION_WMO,
    BUFKIT_MODELS,
    PSU_MODELS,
    DATA_SOURCES,
    TORNADO_SCAN_STATIONS,
)

# ── Utilities ───────────────────────────────────────────────────────
from .utils import get_latest_sounding_time, find_nearest_station

# ── Fetchers ────────────────────────────────────────────────────────
from .fetchers import (
    fetch_iem_sounding,
    fetch_wyoming_sounding,
    fetch_rap_sounding,
    fetch_bufkit_sounding,
    fetch_psu_bufkit,
    fetch_acars_sounding,
    fetch_sounding,
    _parse_bufkit,
)

# ── VAD / VWP ───────────────────────────────────────────────────────
from .vad import fetch_vad_data, fetch_vwp_timeseries, plot_vwp

# ── Tornado risk ────────────────────────────────────────────────────
from .tornado import _quick_tornado_score, _quick_forecast_score, find_highest_tornado_risk

# ── Merge ───────────────────────────────────────────────────────────
from .merge import merge_profiles

# ── Parameters ──────────────────────────────────────────────────────
from .parameters import _mixing_ratio_from_dewpoint, compute_parameters

# ── Plotting ────────────────────────────────────────────────────────
from .plotting import plot_sounding, plot_composite_sounding

# ── CLI ─────────────────────────────────────────────────────────────
from .cli import main

__all__ = [
    "STATIONS", "STATION_WMO", "BUFKIT_MODELS", "PSU_MODELS",
    "DATA_SOURCES", "TORNADO_SCAN_STATIONS",
    "get_latest_sounding_time", "find_nearest_station",
    "fetch_iem_sounding", "fetch_wyoming_sounding", "fetch_rap_sounding",
    "fetch_bufkit_sounding", "fetch_psu_bufkit", "fetch_acars_sounding",
    "fetch_sounding", "_parse_bufkit",
    "fetch_vad_data", "fetch_vwp_timeseries", "plot_vwp",
    "_quick_tornado_score", "_quick_forecast_score", "find_highest_tornado_risk",
    "merge_profiles",
    "_mixing_ratio_from_dewpoint", "compute_parameters",
    "plot_sounding", "plot_composite_sounding",
    "main",
]
