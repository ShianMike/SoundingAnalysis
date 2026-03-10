"""
    Wind-related API endpoints for VAD profiles, VWP displays, and gridded wind fields.
"""
import base64
import io
import math
import time

import requests as _requests
from flask import Blueprint, jsonify, request

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sounding import fetch_vad_data, fetch_vwp_timeseries, plot_vwp

bp = Blueprint("wind", __name__)


@bp.route("/api/vad", methods=["GET"])
def get_vad():
    """Fetch latest NEXRAD VAD Wind Profile for a given radar."""
    radar = request.args.get("radar", "").upper()
    if not radar or len(radar) < 3:
        return jsonify({"error": "Provide a valid NEXRAD radar ID (e.g. ?radar=KTLX)"}), 400
    try:
        result = fetch_vad_data(radar)
        if not result or not result.get("winds"):
            return jsonify({"error": f"No VAD data available for {radar}"}), 404
        return jsonify(result)
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@bp.route("/api/vwp-display", methods=["GET"])
def vwp_display():
    """Generate a VWP time-height display image."""
    radar = request.args.get("radar", "").upper()
    hours = int(request.args.get("hours", 12))
    if not radar or len(radar) < 3:
        return jsonify({"error": "Provide a valid NEXRAD radar ID (e.g. ?radar=KTLX)"}), 400
    hours = max(1, min(hours, 48))
    try:
        vwp_data = fetch_vwp_timeseries(radar, hours=hours)
        if not vwp_data["snapshots"]:
            return jsonify({"error": f"No VWP data available for {radar}"}), 404
        fig = plot_vwp(vwp_data)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight",
                    facecolor=fig.get_facecolor(), edgecolor="none")
        plt.close(fig)
        buf.seek(0)
        image_b64 = base64.b64encode(buf.read()).decode("utf-8")
        return jsonify({
            "image": image_b64,
            "radar": radar,
            "snapshots": len(vwp_data["snapshots"]),
            "timeRange": {
                "start": vwp_data["snapshots"][0]["time"].strftime("%Y-%m-%d %H:%MZ"),
                "end": vwp_data["snapshots"][-1]["time"].strftime("%Y-%m-%d %H:%MZ"),
            },
        })
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


# ── Animated wind-field grid (Open-Meteo GFS) ──────────────────────
_WF_CACHE = {}
_WF_CACHE_TTL = 1800   # 30 min
_WF_MAX_HOUR_OFFSET = 12

# CONUS grid at ~2° resolution  (14 lat × 31 lon = 434 points)
_WF_LATS = list(range(24, 51, 2))     # [24, 26, 28, ..., 50]
_WF_LONS = list(range(-126, -65, 2))  # [-126, -124, ..., -66]

_VALID_LEVELS = {"surface", "500"}


@bp.route("/api/wind-field", methods=["GET"])
def wind_field():
    """Return gridded wind U/V components for CONUS.

    Query params:
      level – "surface" (10 m) or "500" (500 hPa steering flow, default).
      hourOffset – forecast hour offset from now (0-12, default 0).
    """
    level = request.args.get("level", "500")
    if level not in _VALID_LEVELS:
        level = "500"

    try:
        hour_offset = int(request.args.get("hourOffset", 0))
    except (TypeError, ValueError):
        hour_offset = 0
    hour_offset = max(0, min(hour_offset, _WF_MAX_HOUR_OFFSET))

    payload_key = f"wf_payload_{level}_{hour_offset}"
    cached = _WF_CACHE.get(payload_key)
    now = time.time()
    if cached and (now - cached["ts"]) < _WF_CACHE_TTL:
        return jsonify(cached["payload"])

    ny, nx = len(_WF_LATS), len(_WF_LONS)
    all_lats, all_lons = [], []
    for la in _WF_LATS:
        for lo in _WF_LONS:
            all_lats.append(la)
            all_lons.append(lo)

    lat_str = ",".join(f"{la:.1f}" for la in all_lats)
    lon_str = ",".join(f"{lo:.1f}" for lo in all_lons)

    source_key = f"wf_source_{level}"
    source_cached = _WF_CACHE.get(source_key)
    if source_cached and (now - source_cached["ts"]) < _WF_CACHE_TTL:
        results = source_cached["results"]
    else:
        hourly_vars = (
            "wind_speed_10m,wind_direction_10m"
            if level == "surface"
            else "wind_speed_500hPa,wind_direction_500hPa"
        )
        # Fetch the full 0..N hour block once, then slice per request.
        url = (
            f"https://api.open-meteo.com/v1/forecast?"
            f"latitude={lat_str}&longitude={lon_str}"
            f"&hourly={hourly_vars}"
            f"&forecast_hours={_WF_MAX_HOUR_OFFSET + 1}&wind_speed_unit=ms&timezone=auto"
        )
        try:
            resp = _requests.get(url, timeout=30)
            resp.raise_for_status()
            results = resp.json()
            _WF_CACHE[source_key] = {"results": results, "ts": now}
        except Exception as e:
            print(f"[WIND-FIELD] Open-Meteo error: {e}")
            # Graceful fallback: if source cache exists but stale, still use it.
            if source_cached and source_cached.get("results") is not None:
                results = source_cached["results"]
            else:
                return jsonify({"error": f"Failed to fetch wind data: {e}"}), 502

    u_data, v_data = [], []
    items = results if isinstance(results, list) else [results]
    for r in items:
        h = r.get("hourly", {})
        if level == "surface":
            spd_arr = h.get("wind_speed_10m", [0])
            dir_arr = h.get("wind_direction_10m", [0])
        else:
            spd_arr = h.get("wind_speed_500hPa", [0])
            dir_arr = h.get("wind_direction_500hPa", [0])
        idx = min(hour_offset, max(0, len(spd_arr) - 1), max(0, len(dir_arr) - 1))
        spd = (spd_arr[idx] if spd_arr else 0) or 0
        ddir = (dir_arr[idx] if dir_arr else 0) or 0
        rad = math.radians(ddir)
        u_data.append(round(-spd * math.sin(rad), 2))
        v_data.append(round(-spd * math.cos(rad), 2))

    payload = {
        "nx": nx, "ny": ny,
        "la1": _WF_LATS[0], "la2": _WF_LATS[-1],
        "lo1": _WF_LONS[0], "lo2": _WF_LONS[-1],
        "dx": _WF_LONS[1] - _WF_LONS[0],
        "dy": _WF_LATS[1] - _WF_LATS[0],
        "uData": u_data,
        "vData": v_data,
        "level": level,
        "hourOffset": hour_offset,
        "validAt": now,
    }
    _WF_CACHE[payload_key] = {"payload": payload, "ts": now}
    return jsonify(payload)
