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
from datetime import datetime, timedelta, timezone
import time
import requests as _requests

from flask import Flask, jsonify, request, send_from_directory
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

FRONTEND_DIR = os.path.join(os.path.dirname(__file__), "frontend", "dist")

app = Flask(__name__, static_folder=FRONTEND_DIR, static_url_path="")
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
            "wmo": STATION_WMO.get(code, ""),
        })
    return jsonify(stations)


@app.route("/api/stations/intl", methods=["GET"])
def list_intl_stations():
    """Return IGRAv2 international stations for the global source."""
    from sounding import INTL_STATIONS
    stations = []
    for wmo_id, (name, lat, lon) in sorted(INTL_STATIONS.items()):
        stations.append({"id": wmo_id, "name": name, "lat": lat, "lon": lon})
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
    storm_motion = body.get("stormMotion")     # {direction, speed} or null
    surface_mod = body.get("surfaceMod")       # {temperature, dewpoint, wind_speed, wind_direction} or null

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

    # IGRAv2 accepts any WMO station ID; auto-convert 3-letter code if known
    if source == "igrav2":
        if station and station in STATION_WMO:
            station = STATION_WMO[station]
        if not station:
            return jsonify({"error": "IGRAv2 requires a WMO station ID (e.g. 72451, 47646)."}), 400

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
        params = compute_parameters(data, storm_motion=storm_motion, surface_mod=surface_mod)

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
        "ecape": round(params.get("ecape", 0), 1),
        "piecewiseCape": params.get("piecewise_cape", []),
        "surfaceModified": params.get("surface_modified", False),
        "customStormMotion": params.get("custom_storm_motion", False),
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
    scan_workers = int(os.environ.get("SCAN_WORKERS", "3"))

    try:
        with ThreadPoolExecutor(max_workers=scan_workers) as executor:
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


# ─── Time-Series ────────────────────────────────────────────────────

@app.route("/api/time-series", methods=["POST", "OPTIONS"])
def time_series():
    """
    Fetch sounding parameters for a station across a date range.
    Returns a time-series array of parameter snapshots.

    Expected JSON body:
      {
        "station": "OUN",
        "source": "obs",
        "startDate": "2026021200",   // YYYYMMDDHH  (optional – defaults to 7 days ago)
        "endDate": "2026022600",     // YYYYMMDDHH  (optional – defaults to latest)
        "days": 7,                   // fallback if no dates given (max 14)
        "lat": 35.22,
        "lon": -97.46,
        "model": "hrrr",
        "fhour": 0
      }
    """
    if request.method == "OPTIONS":
        return "", 204

    body = request.get_json(force=True)

    try:
        station = body.get("station")
        source = body.get("source", "obs").lower()
        lat = body.get("lat")
        lon = body.get("lon")
        model = body.get("model", "rap")
        fhour = body.get("fhour", 0)

        if station:
            station = station.upper()

        # Validate
        if source == "obs" and (not station or station not in STATION_WMO):
            return jsonify({"error": f"Unknown station '{station}'."}), 400

        if source in ("rap", "era5") and lat is None and lon is None:
            if station and station in STATIONS:
                _, lat, lon = STATIONS[station]
            else:
                return jsonify({"error": "RAP/ERA5 requires lat/lon or a known station."}), 400

        # Resolve date range
        latest = get_latest_sounding_time()
        start_str = body.get("startDate")
        end_str = body.get("endDate")

        if end_str:
            try:
                end_dt = datetime.strptime(str(end_str)[:10], "%Y%m%d%H").replace(tzinfo=timezone.utc)
            except ValueError:
                end_dt = latest
        else:
            end_dt = latest

        if start_str:
            try:
                start_dt = datetime.strptime(str(start_str)[:10], "%Y%m%d%H").replace(tzinfo=timezone.utc)
            except ValueError:
                start_dt = end_dt - timedelta(days=7)
        else:
            days = min(int(body.get("days", 7)), 14)
            start_dt = end_dt - timedelta(days=days)

        # Cap range at 14 days
        if (end_dt - start_dt).days > 14:
            start_dt = end_dt - timedelta(days=14)

        span_days = (end_dt - start_dt).days

        # Build sounding times: both 00Z+12Z for ≤7 days, 12Z only for >7
        times = []
        hours = (0, 12) if span_days <= 7 else (12,)
        d = start_dt.replace(hour=0, minute=0, second=0, microsecond=0)
        while d <= end_dt:
            for h in hours:
                t = d.replace(hour=h)
                if start_dt <= t <= end_dt:
                    times.append(t)
            d += timedelta(days=1)
        times.sort()

        def _fetch_one(dt):
            """Fetch + compute params for one time step."""
            try:
                data = fetch_sounding(
                    station_id=station, dt=dt, source=source,
                    lat=lat, lon=lon, model=model, fhour=fhour,
                )
                params = compute_parameters(data)
                return {
                    "date": dt.strftime("%Y-%m-%d %HZ"),
                    "ts": dt.isoformat(),
                    "params": _serialize_params(params, data, station, dt, source),
                }
            except Exception:
                return None

        points = []
        wall_start = time.monotonic()
        # Cloud Run allows 300s; Koyeb has 100s proxy timeout
        WALL_LIMIT = int(os.environ.get("TS_WALL_LIMIT", "85"))
        ts_workers = int(os.environ.get("TS_WORKERS", "2"))

        with ThreadPoolExecutor(max_workers=ts_workers) as executor:
            futures = {executor.submit(_fetch_one, t): t for t in times}
            for future in as_completed(futures, timeout=80):
                # Bail if we're running out of time
                if time.monotonic() - wall_start > WALL_LIMIT:
                    break
                try:
                    r = future.result(timeout=20)
                    if r is not None:
                        points.append(r)
                except Exception:
                    pass

        # Sort chronologically
        points.sort(key=lambda p: p["ts"])

        station_name = STATIONS.get(station, (station or source.upper(),))[0]

        return jsonify({
            "station": station or f"{lat},{lon}",
            "stationName": station_name,
            "source": source,
            "startDate": start_dt.strftime("%Y-%m-%d"),
            "endDate": end_dt.strftime("%Y-%m-%d %HZ"),
            "resolution": "00Z & 12Z" if span_days <= 7 else "12Z only",
            "points": points,
        })

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


# ── Compare endpoint ────────────────────────────────────────────────
@app.route("/api/compare", methods=["POST", "OPTIONS"])
def compare_soundings():
    """
    Fetch multiple soundings in parallel for side-by-side comparison.

    Expected JSON body:
      {
        "soundings": [
          { "source": "obs", "station": "OUN", "date": "2024061200" },
          { "source": "obs", "station": "FWD", "date": "2024061200" },
          ...
        ]
      }
    Max 4 soundings at a time.
    """
    if request.method == "OPTIONS":
        return "", 204

    body = request.get_json(force=True)
    items = body.get("soundings", [])

    if not items or not isinstance(items, list):
        return jsonify({"error": "Provide a 'soundings' array."}), 400
    if len(items) > 4:
        return jsonify({"error": "Maximum 4 soundings per comparison."}), 400

    compare_workers = int(os.environ.get("COMPARE_WORKERS", "2"))

    def _fetch_one_compare(item):
        """Fetch + compute + plot a single sounding item."""
        src = item.get("source", "obs").lower()
        stn = item.get("station")
        date_str = item.get("date")
        lat_v = item.get("lat")
        lon_v = item.get("lon")
        mdl = item.get("model", "rap")
        fhr = item.get("fhour", 0)

        if date_str:
            dt = datetime.strptime(str(date_str), "%Y%m%d%H").replace(
                tzinfo=timezone.utc
            )
        else:
            dt = get_latest_sounding_time()

        if src == "obs" and stn is None and lat_v is not None and lon_v is not None:
            stn = find_nearest_station(lat_v, lon_v)
        if stn:
            stn = stn.upper()

        if src in ("rap", "era5") and lat_v is None and lon_v is None:
            if stn and stn in STATIONS:
                _, lat_v, lon_v = STATIONS[stn]

        data = fetch_sounding(
            station_id=stn, dt=dt, source=src,
            lat=lat_v, lon=lon_v, model=mdl, fhour=fhr,
        )
        params = compute_parameters(data)

        plot_id = stn or src.upper()
        fig = plot_sounding(data, params, plot_id, dt)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=180, facecolor="#0d0d0d")
        plt.close(fig)
        buf.seek(0)
        image_b64 = base64.b64encode(buf.read()).decode("utf-8")

        serialized = _serialize_params(params, data, stn, dt, src)

        return {
            "image": image_b64,
            "params": serialized,
            "meta": {
                "station": stn,
                "stationName": STATIONS.get(stn, (stn or src.upper(),))[0],
                "source": src,
                "date": dt.strftime("%Y-%m-%d %HZ"),
                "levels": len(data["pressure"]),
                "sfcPressure": round(float(data["pressure"][0].magnitude)),
                "topPressure": round(float(data["pressure"][-1].magnitude)),
            },
        }

    try:
        results = []
        with ThreadPoolExecutor(max_workers=compare_workers) as executor:
            futures = {executor.submit(_fetch_one_compare, item): i for i, item in enumerate(items)}
            for future in as_completed(futures):
                idx = futures[future]
                try:
                    results.append((idx, future.result()))
                except Exception as e:
                    results.append((idx, {"error": str(e), "meta": items[idx]}))

        # Sort by original order
        results.sort(key=lambda x: x[0])
        return jsonify({"soundings": [r[1] for r in results]})

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


@app.route("/api/composite", methods=["POST", "OPTIONS"])
def composite_sounding():
    """
    Generate a composite overlay plot with multiple soundings on one Skew-T.

    Expected JSON body:
      {
        "soundings": [
          { "source": "obs", "station": "OUN", "date": "2024061200" },
          { "source": "obs", "station": "FWD", "date": "2024061200" },
          ...
        ]
      }
    Max 6 soundings.
    """
    if request.method == "OPTIONS":
        return "", 204

    body = request.get_json(force=True)
    items = body.get("soundings", [])

    if not items or not isinstance(items, list):
        return jsonify({"error": "Provide a 'soundings' array."}), 400
    if len(items) > 6:
        return jsonify({"error": "Maximum 6 soundings per composite."}), 400

    from sounding import plot_composite_sounding

    profiles = []
    errors = []
    for item in items:
        try:
            src = item.get("source", "obs").lower()
            stn = item.get("station")
            date_str = item.get("date")
            lat_v = item.get("lat")
            lon_v = item.get("lon")
            mdl = item.get("model", "rap")
            fhr = item.get("fhour", 0)

            if date_str:
                dt = datetime.strptime(str(date_str), "%Y%m%d%H").replace(tzinfo=timezone.utc)
            else:
                dt = get_latest_sounding_time()

            if src == "obs" and stn is None and lat_v is not None and lon_v is not None:
                stn = find_nearest_station(lat_v, lon_v)
            if stn:
                stn = stn.upper()

            if src in ("rap", "era5") and lat_v is None and lon_v is None:
                if stn and stn in STATIONS:
                    _, lat_v, lon_v = STATIONS[stn]

            data = fetch_sounding(
                station_id=stn, dt=dt, source=src,
                lat=lat_v, lon=lon_v, model=mdl, fhour=fhr,
            )
            params = compute_parameters(data)

            label = f"{stn or src.upper()} {dt.strftime('%d/%HZ')}"
            profiles.append({"data": data, "params": params, "label": label})
        except Exception as ex:
            errors.append({"item": item, "error": str(ex)})

    if not profiles:
        return jsonify({"error": "All soundings failed to fetch.", "details": errors}), 500

    try:
        fig = plot_composite_sounding(profiles, title="Composite Sounding Overlay")
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=180, facecolor="#0d0d0d")
        plt.close(fig)
        buf.seek(0)
        image_b64 = base64.b64encode(buf.read()).decode("utf-8")

        return jsonify({
            "image": image_b64,
            "count": len(profiles),
            "errors": errors if errors else None,
        })

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


# ─── SPC Convective Outlook ────────────────────────────────────────
_SPC_OUTLOOK_URLS = {
    "1": "https://www.spc.noaa.gov/products/outlook/day1otlk_cat.lyr.geojson",
    "2": "https://www.spc.noaa.gov/products/outlook/day2otlk_cat.lyr.geojson",
    "3": "https://www.spc.noaa.gov/products/outlook/day3otlk_cat.lyr.geojson",
}

# Simple in-memory cache: { day: { "data": ..., "ts": ... } }
_spc_cache = {}
_SPC_CACHE_TTL = 600  # 10 minutes


@app.route("/api/spc-outlook", methods=["GET"])
def spc_outlook():
    """Proxy SPC categorical outlook GeoJSON (avoids CORS)."""
    day = request.args.get("day", "1")
    if day not in _SPC_OUTLOOK_URLS:
        return jsonify({"error": f"Invalid day: {day}. Use 1, 2, or 3."}), 400

    # Return cached if fresh
    cached = _spc_cache.get(day)
    if cached and (time.time() - cached["ts"]) < _SPC_CACHE_TTL:
        return jsonify(cached["data"])

    url = _SPC_OUTLOOK_URLS[day]
    try:
        resp = _requests.get(url, timeout=15)
        resp.raise_for_status()
        geojson = resp.json()
        _spc_cache[day] = {"data": geojson, "ts": time.time()}
        return jsonify(geojson)
    except Exception as e:
        print(f"[SPC] Failed to fetch Day {day} outlook: {e}")
        return jsonify({"error": f"Failed to fetch SPC outlook: {e}"}), 502


# ─── Feedback ───────────────────────────────────────────────────────
# Use /tmp on Cloud Run (ephemeral but always writable)
_fb_dir = "/tmp" if os.environ.get("K_SERVICE") else os.path.dirname(__file__)
FEEDBACK_FILE = os.path.join(_fb_dir, "feedback.json")


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

    # Log to stdout so it appears in Cloud Run logs (file is ephemeral)
    print(f"[FEEDBACK] type={entry['type']} | {entry['text'][:200]}")

    feedback = _load_feedback()
    feedback.append(entry)
    _save_feedback(feedback)

    return jsonify({"ok": True, "id": entry["id"]})


@app.route("/api/feedback", methods=["GET"])
def get_feedback():
    """Retrieve all feedback (for admin/developer review)."""
    return jsonify(_load_feedback())


# ─── SPA catch-all ──────────────────────────────────────────────────
@app.route("/", defaults={"path": ""})
@app.route("/<path:path>")
def serve_spa(path):
    """Serve the React SPA. API routes are matched first by Flask."""
    file_path = os.path.join(FRONTEND_DIR, path)
    if path and os.path.isfile(file_path):
        return send_from_directory(FRONTEND_DIR, path)
    return send_from_directory(FRONTEND_DIR, "index.html")


if __name__ == "__main__":
    print("Starting Sounding API on http://localhost:5000")
    app.run(debug=True, port=5000)
