"""
Flask API for the Sounding Analysis Tool.
Exposes endpoints to fetch sounding data and return analysis results
(plot image + computed parameters) to the React frontend.
"""

import base64
import io
import json
import os
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

from flask import Flask, jsonify, request
from flask_cors import CORS

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sounding import (
    STATIONS,
    STATION_WMO,
    DATA_SOURCES,
    BUFKIT_MODELS,
    TORNADO_SCAN_STATIONS,
    fetch_sounding,
    compute_parameters,
    plot_sounding,
    get_latest_sounding_time,
    find_nearest_station,
    _quick_tornado_score,
)

app = Flask(__name__)
CORS(app, resources={r"/api/*": {"origins": "*"}}, supports_credentials=False)


@app.after_request
def add_cors_headers(response):
    response.headers["Access-Control-Allow-Origin"] = "*"
    response.headers["Access-Control-Allow-Methods"] = "GET, POST, OPTIONS"
    response.headers["Access-Control-Allow-Headers"] = "Content-Type"
    return response


@app.route("/api/health", methods=["GET"])
def health():
    """Lightweight health-check endpoint."""
    return jsonify({"status": "ok"})


# ─── Helper ────────────────────────────────────────────────────────
def _fmt(val, unit_str="", decimals=0):
    """Format a pint-style value, returning None for missing data."""
    if val is None:
        return None
    try:
        v = val.magnitude if hasattr(val, "magnitude") else val
        if decimals == 0:
            return round(float(v))
        return round(float(v), decimals)
    except Exception:
        return None


# ─── Endpoints ─────────────────────────────────────────────────────

@app.route("/api/stations", methods=["GET"])
def list_stations():
    """Return all known sounding stations."""
    stations = []
    for code, (name, lat, lon) in sorted(STATIONS.items()):
        stations.append({
            "id": code,
            "name": name,
            "lat": lat,
            "lon": lon,
        })
    return jsonify(stations)


@app.route("/api/sources", methods=["GET"])
def list_sources():
    """Return available data sources and BUFKIT models."""
    sources = [{"id": k, "name": v} for k, v in DATA_SOURCES.items()]
    models = [{"id": k, "name": v} for k, v in BUFKIT_MODELS.items()]
    return jsonify({"sources": sources, "models": models})


@app.route("/api/sounding", methods=["POST", "OPTIONS"])
def get_sounding():
    if request.method == "OPTIONS":
        return "", 204
    """
    Fetch sounding data, compute parameters, generate plot.

    Expected JSON body:
      {
        "source": "obs",           // obs | rap | bufkit | era5 | acars
        "station": "OUN",          // 3-letter station ID (optional for rap/era5)
        "date": "2024061200",      // YYYYMMDDHH  (optional, defaults to latest)
        "lat": 35.22,              // required for rap/era5 when no station
        "lon": -97.46,
        "model": "hrrr",           // for bufkit only
        "fhour": 0                 // for bufkit only
      }

    Returns JSON:
      {
        "image": "<base64-encoded PNG>",
        "params": { ... },
        "meta": { ... }
      }
    """
    body = request.get_json(force=True)

    source = body.get("source", "obs").lower()
    station = body.get("station")
    date_str = body.get("date")
    lat = body.get("lat")
    lon = body.get("lon")
    model = body.get("model", "rap")
    fhour = body.get("fhour", 0)

    # Parse date
    if date_str:
        try:
            dt = datetime.strptime(str(date_str), "%Y%m%d%H").replace(
                tzinfo=timezone.utc
            )
        except ValueError:
            return jsonify({"error": f"Invalid date format '{date_str}'. Use YYYYMMDDHH."}), 400
    else:
        dt = get_latest_sounding_time()

    # Resolve station for obs if only lat/lon given
    if source == "obs" and station is None and lat is not None and lon is not None:
        station = find_nearest_station(lat, lon)

    if station:
        station = station.upper()

    # Validate
    if source == "obs" and (not station or station not in STATION_WMO):
        return jsonify({"error": f"Unknown station '{station}'. Provide a valid 3-letter ID."}), 400

    if source in ("rap", "era5") and lat is None and lon is None:
        if station and station in STATIONS:
            _, lat, lon = STATIONS[station]
        else:
            return jsonify({"error": "RAP/ERA5 source requires lat/lon or a known station."}), 400

    try:
        # 1) Fetch
        data = fetch_sounding(
            station_id=station,
            dt=dt,
            source=source,
            lat=lat,
            lon=lon,
            model=model,
            fhour=fhour,
        )

        # 2) Compute
        params = compute_parameters(data)

        # 3) Plot → base64 PNG
        plot_id = station or source.upper()
        fig = plot_sounding(data, params, plot_id, dt)

        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=180, facecolor="#0d0d0d")
        plt.close(fig)
        buf.seek(0)
        image_b64 = base64.b64encode(buf.read()).decode("utf-8")

        # 4) Serialize params
        serialized = _serialize_params(params, data, station, dt, source)

        return jsonify({
            "image": image_b64,
            "params": serialized,
            "meta": {
                "station": station,
                "stationName": STATIONS.get(station, (station or source.upper(),))[0],
                "source": source,
                "date": dt.strftime("%Y-%m-%d %HZ"),
                "levels": len(data["pressure"]),
                "sfcPressure": round(float(data["pressure"][0].magnitude)),
                "topPressure": round(float(data["pressure"][-1].magnitude)),
            },
        })

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


def _serialize_params(params, data, station, dt, source):
    """Extract key computed parameters into a JSON-friendly dict."""
    return {
        # Thermodynamic
        "sbCape": _fmt(params.get("sb_cape")),
        "sbCin": _fmt(params.get("sb_cin")),
        "sbLclM": round(params["sb_lcl_m"]) if params.get("sb_lcl_m") is not None else None,
        "sbLclP": _fmt(params.get("sb_lcl_p")),
        "muCape": _fmt(params.get("mu_cape")),
        "muCin": _fmt(params.get("mu_cin")),
        "muLclM": round(params["mu_lcl_m"]) if params.get("mu_lcl_m") is not None else None,
        "mlCape": _fmt(params.get("ml_cape")),
        "mlCin": _fmt(params.get("ml_cin")),
        "mlLclM": round(params["ml_lcl_m"]) if params.get("ml_lcl_m") is not None else None,
        "dcape": _fmt(params.get("dcape")),
        "lr03": round(params["lr_03"], 1) if params.get("lr_03") is not None else None,
        "lr36": round(params["lr_36"], 1) if params.get("lr_36") is not None else None,
        "pwat": _fmt(params.get("pwat"), decimals=1),
        "frzLevel": round(params["frz_level"]) if params.get("frz_level") is not None else None,
        "wbo": round(params["wbo"]) if params.get("wbo") is not None else None,
        "stp": round(params.get("stp", 0), 1),
        "scp": round(params.get("scp", 0), 1),
        "ship": round(params.get("ship", 0), 1),
        "dcp": round(params.get("dcp", 0), 1),
        "rh01": round(params["rh_0_1km"]) if params.get("rh_0_1km") is not None else None,
        "rh13": round(params["rh_1_3km"]) if params.get("rh_1_3km") is not None else None,
        "rh36": round(params["rh_3_6km"]) if params.get("rh_3_6km") is not None else None,
        # Kinematic
        "bwd500m": _fmt(params.get("bwd_500m")),
        "bwd1km": _fmt(params.get("bwd_1km")),
        "bwd3km": _fmt(params.get("bwd_3km")),
        "bwd6km": _fmt(params.get("bwd_6km")),
        "srh500m": _fmt(params.get("srh_500m")),
        "srh1km": _fmt(params.get("srh_1km")),
        "srh3km": _fmt(params.get("srh_3km")),
    }


@app.route("/api/risk-scan", methods=["POST", "OPTIONS"])
def risk_scan():
    """
    Scan stations and return risk scores.
    Uses ThreadPoolExecutor for parallel fetching.
    """
    if request.method == "OPTIONS":
        return "", 204

    body = request.get_json(force=True) if request.data else {}
    date_str = body.get("date")

    if date_str:
        try:
            dt = datetime.strptime(str(date_str), "%Y%m%d%H").replace(
                tzinfo=timezone.utc
            )
        except ValueError:
            return jsonify({"error": f"Invalid date format '{date_str}'."}), 400
    else:
        dt = get_latest_sounding_time()

    def _scan_one(sid):
        """Score a single station; returns dict or None."""
        try:
            score = _quick_tornado_score(sid, dt)
            if score is None:
                return None
            stp_score, raw_score, cape, srh, bwd, scp, ship, dcp = score
            return {
                "id": sid,
                "name": STATIONS.get(sid, (sid,))[0],
                "lat": STATIONS.get(sid, ("", 0, 0))[1],
                "lon": STATIONS.get(sid, ("", 0, 0))[2],
                "stp": round(stp_score, 2),
                "raw": round(raw_score, 2),
                "cape": round(cape),
                "srh": round(srh),
                "bwd": round(bwd),
                "scp": round(scp, 2),
                "ship": round(ship, 2),
                "dcp": round(dcp, 2),
            }
        except Exception:
            return None

    results = []
    station_ids = list(STATIONS.keys())

    try:
        with ThreadPoolExecutor(max_workers=3) as executor:
            futures = {executor.submit(_scan_one, sid): sid for sid in station_ids}
            for future in as_completed(futures, timeout=150):
                try:
                    r = future.result(timeout=30)
                    if r is not None:
                        results.append(r)
                except Exception:
                    pass
    except Exception as exc:
        # If timeout/error, return whatever we collected so far
        print(f"[risk-scan] partial results due to: {exc}")

    results.sort(key=lambda r: (r["stp"], r["raw"]), reverse=True)

    return jsonify({
        "date": dt.strftime("%Y-%m-%d %HZ"),
        "stations": results,
    })


# ─── Feedback ───────────────────────────────────────────────────────
FEEDBACK_FILE = os.path.join(os.path.dirname(__file__), "feedback.json")


def _load_feedback():
    if os.path.exists(FEEDBACK_FILE):
        try:
            with open(FEEDBACK_FILE, "r", encoding="utf-8") as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError):
            return []
    return []


def _save_feedback(data):
    with open(FEEDBACK_FILE, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


@app.route("/api/feedback", methods=["POST"])
def submit_feedback():
    """Store user feedback/suggestions."""
    body = request.get_json(force=True)
    text = (body.get("text") or "").strip()
    if not text:
        return jsonify({"error": "Feedback text is required"}), 400

    entry = {
        "id": int(datetime.now(timezone.utc).timestamp() * 1000),
        "type": body.get("type", "suggestion"),
        "text": text,
        "date": datetime.now(timezone.utc).isoformat(),
        "userAgent": body.get("userAgent", ""),
    }

    feedback = _load_feedback()
    feedback.append(entry)
    _save_feedback(feedback)

    return jsonify({"ok": True, "id": entry["id"]})


@app.route("/api/feedback", methods=["GET"])
def get_feedback():
    """Retrieve all feedback (for admin/developer review)."""
    return jsonify(_load_feedback())


if __name__ == "__main__":
    print("Starting Sounding API on http://localhost:5000")
    app.run(debug=True, port=5000)
