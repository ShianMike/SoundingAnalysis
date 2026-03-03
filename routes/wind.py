"""
Wind-related routes: VAD and VWP display.
"""
import base64
import io

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
