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


# ── Animated wind-field grid (Open-Meteo GFS 10 m wind) ─────────────
_WF_CACHE = {}
_WF_CACHE_TTL = 1800   # 30 min

# CONUS grid at ~2.5° lat × 3° lon resolution
_WF_LATS = [24, 26.5, 29, 31.5, 34, 36.5, 39, 41.5, 44, 46.5, 49]
_WF_LONS = [-126, -123, -120, -117, -114, -111, -108, -105, -102,
            -99, -96, -93, -90, -87, -84, -81, -78, -75, -72, -69, -66]


@bp.route("/api/wind-field", methods=["GET"])
def wind_field():
    """Return gridded 10 m wind U/V components for CONUS."""
    cached = _WF_CACHE.get("data")
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

    url = (
        f"https://api.open-meteo.com/v1/forecast?"
        f"latitude={lat_str}&longitude={lon_str}"
        f"&current=wind_speed_10m,wind_direction_10m"
        f"&wind_speed_unit=ms&timezone=auto"
    )
    try:
        resp = _requests.get(url, timeout=30)
        resp.raise_for_status()
        results = resp.json()
    except Exception as e:
        print(f"[WIND-FIELD] Open-Meteo error: {e}")
        return jsonify({"error": f"Failed to fetch wind data: {e}"}), 502

    u_data, v_data = [], []
    items = results if isinstance(results, list) else [results]
    for r in items:
        c = r.get("current", {})
        spd = c.get("wind_speed_10m", 0) or 0
        ddir = c.get("wind_direction_10m", 0) or 0
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
    }
    _WF_CACHE["data"] = {"payload": payload, "ts": now}
    return jsonify(payload)
