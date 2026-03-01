"""
Flask API for the Sounding Analysis Tool.
Exposes endpoints to fetch sounding data and return analysis results
(plot image + computed parameters) to the React frontend.
"""

import base64
import io
import json
import math
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
    fetch_vad_data,
    fetch_vwp_timeseries,
    plot_vwp,
    compute_parameters,
    plot_sounding,
    merge_profiles,
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
        "source": "obs",           // obs | rap | bufkit | acars
        "station": "OUN",          // 3-letter station ID (optional for rap)
        "date": "2024061200",      // YYYYMMDDHH  (optional, defaults to latest)
        "lat": 35.22,              // required for rap when no station
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
    vad_radar = body.get("vad")                # NEXRAD radar ID (e.g. "KTLX") or null
    smoothing = body.get("smoothing")          # Gaussian sigma (float) or null
    sr_hodograph = body.get("srHodograph", False)  # Storm-relative hodograph mode
    theme = body.get("theme", "dark")              # "dark" or "light"
    colorblind = body.get("colorblind", False)     # color-blind safe palette
    boundary_orientation = body.get("boundaryOrientation")  # degrees or null
    map_zoom = body.get("mapZoom", 1.0)                      # mini-map zoom factor (1-8)

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

    if source == "rap" and lat is None and lon is None:
        if station and station in STATIONS:
            _, lat, lon = STATIONS[station]
        else:
            return jsonify({"error": "RAP source requires lat/lon or a known station."}), 400

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
        params = compute_parameters(data, storm_motion=storm_motion, surface_mod=surface_mod, smoothing=smoothing)

        # 2b) Fetch VAD data if requested
        vad_result = None
        if vad_radar:
            try:
                vad_result = fetch_vad_data(vad_radar)
            except Exception as ve:
                print(f"[VAD] Non-fatal error fetching VAD: {ve}")

        # 3) Plot → base64 PNG
        plot_id = station or source.upper()
        fig = plot_sounding(data, params, plot_id, dt, vad_data=vad_result,
                            sr_hodograph=sr_hodograph, theme=theme,
                            colorblind=colorblind,
                            boundary_orientation=boundary_orientation,
                            map_zoom=map_zoom)

        _facecolor = "#f5f5f5" if theme == "light" else "#0d0d0d"
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=180, facecolor=_facecolor)
        plt.close(fig)
        buf.seek(0)
        image_b64 = base64.b64encode(buf.read()).decode("utf-8")

        # 4) Serialize params
        serialized = _serialize_params(params, data, station, dt, source)

        # 5) Build raw profile array for client-side export
        n = len(data["pressure"])
        profile_rows = []
        for i in range(n):
            row = {
                "p": round(float(data["pressure"][i].magnitude), 1),
                "h": round(float(data["height"][i].magnitude), 1),
                "t": round(float(data["temperature"][i].magnitude), 1),
                "td": round(float(data["dewpoint"][i].magnitude), 1),
            }
            wd_val = data["wind_direction"][i]
            ws_val = data["wind_speed"][i]
            wd_mag = float(wd_val.magnitude) if hasattr(wd_val, "magnitude") else float(wd_val)
            ws_mag = float(ws_val.magnitude) if hasattr(ws_val, "magnitude") else float(ws_val)
            row["wd"] = round(wd_mag, 1) if not math.isnan(wd_mag) else None
            row["ws"] = round(ws_mag, 1) if not math.isnan(ws_mag) else None
            profile_rows.append(row)

        # Station elevation (metres MSL) — used by CM1 input_sounding format
        stn_info = data.get("station_info", {})
        sfc_elev = stn_info.get("elev", data["height"][0].magnitude if n > 0 else 0)

        return jsonify({
            "image": image_b64,
            "params": serialized,
            "profile": profile_rows,
            "meta": {
                "station": station,
                "stationName": STATIONS.get(station, (station or source.upper(),))[0],
                "source": source,
                "date": dt.strftime("%Y-%m-%d %HZ"),
                "levels": len(data["pressure"]),
                "sfcPressure": round(float(data["pressure"][0].magnitude)),
                "topPressure": round(float(data["pressure"][-1].magnitude)),
                "sfcElevation": round(float(sfc_elev)) if sfc_elev is not None else 0,
                "vadRadar": vad_result.get("radar") if vad_result and isinstance(vad_result, dict) else None,
                "vadTime": vad_result.get("meta", {}).get("time") if vad_result and isinstance(vad_result, dict) else None,
            },
        })

    except ValueError as e:
        # Data not found / invalid input — not a server error
        return jsonify({"error": str(e)}), 400
    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


# ─── Custom Sounding Upload ───────────────────────────────────────
@app.route("/api/custom-sounding", methods=["POST", "OPTIONS"])
def custom_sounding():
    """Accept raw sounding profile data (text) and return analysis.

    Expected JSON body:
      {
        "text": "...raw sounding data lines...",
        "format": "auto" | "csv" | "sharppy" | "cm1",
        "theme": "dark" | "light",
        "colorblind": false
      }
    """
    if request.method == "OPTIONS":
        return "", 204

    import numpy as np
    from metpy.units import units as mpu

    body = request.get_json(force=True)
    raw_text = body.get("text", "").strip()
    fmt = body.get("format", "auto").lower()
    theme = body.get("theme", "dark")
    colorblind = body.get("colorblind", False)

    if not raw_text:
        return jsonify({"error": "No sounding data provided."}), 400

    try:
        lines = [l.strip() for l in raw_text.splitlines() if l.strip()]

        # ── Auto-detect format ───────────────────────────────────
        if fmt == "auto":
            first_data = [l for l in lines if not l.startswith("%") and not l.startswith("#")]
            if first_data and len(first_data[0].split()) == 3 and first_data[0].replace(".", "").replace("-", "").replace(" ", "").isdigit():
                fmt = "cm1"
            elif any("," in l for l in first_data[:5]):
                fmt = "csv"
            else:
                fmt = "sharppy"

        p_arr, h_arr, t_arr, td_arr, wd_arr, ws_arr = [], [], [], [], [], []

        if fmt == "csv":
            # Expected CSV: P(hPa), H(m), T(C), Td(C), WD(deg), WS(kt)
            for line in lines:
                if line.startswith("#") or line.startswith("P") or line.startswith("p"):
                    continue
                parts = [x.strip() for x in line.split(",")]
                if len(parts) < 6:
                    continue
                try:
                    p_arr.append(float(parts[0]))
                    h_arr.append(float(parts[1]))
                    t_arr.append(float(parts[2]))
                    td_arr.append(float(parts[3]))
                    wd_arr.append(float(parts[4]))
                    ws_arr.append(float(parts[5]))
                except ValueError:
                    continue

        elif fmt == "cm1":
            # CM1 input_sounding format:
            # Line 1: sfc_pressure(hPa) theta(K) qv(g/kg)
            # Subsequent: z(m) theta(K) qv(g/kg) u(m/s) v(m/s)
            header_parts = lines[0].split()
            sfc_p = float(header_parts[0]) * 100  # hPa -> Pa for later? No, keep hPa
            sfc_p_hPa = float(header_parts[0])
            sfc_theta = float(header_parts[1])
            sfc_qv = float(header_parts[2]) / 1000.0  # g/kg -> kg/kg

            # Compute sfc T from theta: T = theta * (P/1000)^(R/cp)
            sfc_T_K = sfc_theta * (sfc_p_hPa / 1000.0) ** 0.286
            sfc_T_C = sfc_T_K - 273.15
            # Compute sfc Td from qv: Td ≈ (243.5 * ln(e/6.112))/(17.67 - ln(e/6.112))
            e_sfc = sfc_qv * sfc_p_hPa / (0.622 + sfc_qv)  # hPa
            if e_sfc > 0:
                sfc_Td_C = (243.5 * math.log(e_sfc / 6.112)) / (17.67 - math.log(e_sfc / 6.112))
            else:
                sfc_Td_C = sfc_T_C - 20

            p_arr.append(sfc_p_hPa)
            h_arr.append(0)
            t_arr.append(sfc_T_C)
            td_arr.append(sfc_Td_C)
            wd_arr.append(0)
            ws_arr.append(0)

            for line in lines[1:]:
                parts = line.split()
                if len(parts) < 5:
                    continue
                try:
                    z_m = float(parts[0])
                    theta_K = float(parts[1])
                    qv_kgkg = float(parts[2]) / 1000.0
                    u_ms = float(parts[3])
                    v_ms = float(parts[4])
                    # Estimate pressure using hypsometric
                    p_est = sfc_p_hPa * math.exp(-z_m / 8500.0)
                    T_K = theta_K * (p_est / 1000.0) ** 0.286
                    T_C = T_K - 273.15
                    e = qv_kgkg * p_est / (0.622 + qv_kgkg)
                    if e > 0:
                        Td_C = (243.5 * math.log(e / 6.112)) / (17.67 - math.log(e / 6.112))
                    else:
                        Td_C = T_C - 20
                    spd = math.sqrt(u_ms**2 + v_ms**2) * 1.94384  # m/s -> kt
                    dirn = (math.degrees(math.atan2(-u_ms, -v_ms)) + 360) % 360
                    p_arr.append(p_est)
                    h_arr.append(z_m)
                    t_arr.append(T_C)
                    td_arr.append(Td_C)
                    wd_arr.append(dirn)
                    ws_arr.append(spd)
                except ValueError:
                    continue

        else:
            # SHARPpy format: P(hPa)  H(m)  T(C)  Td(C)  WD(deg)  WS(kt)
            # Skip header lines starting with % or non-numeric
            for line in lines:
                if line.startswith("%") or line.startswith("#") or line.startswith("SNPARM"):
                    continue
                parts = line.split()
                if len(parts) < 6:
                    continue
                try:
                    vals = [float(x) for x in parts[:6]]
                    # Skip missing data (SHARPpy uses -9999)
                    if any(v <= -9998 for v in vals):
                        continue
                    p_arr.append(vals[0])
                    h_arr.append(vals[1])
                    t_arr.append(vals[2])
                    td_arr.append(vals[3])
                    wd_arr.append(vals[4])
                    ws_arr.append(vals[5])
                except ValueError:
                    continue

        if len(p_arr) < 5:
            return jsonify({"error": f"Could not parse enough levels ({len(p_arr)} found). "
                           "Need at least 5.  Check format: P(hPa), H(m), T(°C), Td(°C), "
                           "WD(°), WS(kt)."}), 400

        # Build data dict compatible with compute_parameters / plot_sounding
        data = {
            "pressure": np.array(p_arr) * mpu.hPa,
            "height": np.array(h_arr) * mpu.meter,
            "temperature": np.array(t_arr) * mpu.degC,
            "dewpoint": np.array(td_arr) * mpu.degC,
            "wind_direction": np.array(wd_arr) * mpu.degree,
            "wind_speed": np.array(ws_arr) * mpu.knot,
            "station_info": {"name": "Custom Upload", "lat": 35.0, "lon": -97.0, "elev": h_arr[0]},
        }

        params = compute_parameters(data)

        fig = plot_sounding(data, params, "CUSTOM", datetime.now(timezone.utc),
                            theme=theme, colorblind=colorblind)

        _facecolor = "#f5f5f5" if theme == "light" else "#0d0d0d"
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=180, facecolor=_facecolor)
        plt.close(fig)
        buf.seek(0)
        image_b64 = base64.b64encode(buf.read()).decode("utf-8")

        serialized = _serialize_params(params, data, "CUSTOM",
                                       datetime.now(timezone.utc), "custom")

        return jsonify({
            "image": image_b64,
            "params": serialized,
            "meta": {
                "station": "CUSTOM",
                "stationName": "Custom Upload",
                "source": "custom",
                "date": datetime.now(timezone.utc).strftime("%Y-%m-%d %HZ"),
                "levels": len(p_arr),
                "sfcPressure": round(p_arr[0]),
                "topPressure": round(p_arr[-1]),
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
        "cape3km": params.get("cape_3km", 0),
        "cape6km": params.get("cape_6km", 0),
        "dcin": params.get("dcin", 0),
        "ncape": params.get("ncape", 0),
        "surfaceModified": params.get("surface_modified", False),
        "customStormMotion": params.get("custom_storm_motion", False),
        "smoothingApplied": params.get("smoothing_applied", False),
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
        "esrh": _fmt(params.get("esrh")),
        "ebwd": _fmt(params.get("ebwd")),
        "eilBot": round(params["eil_bot_h"]) if params.get("eil_bot_h") is not None else None,
        "eilTop": round(params["eil_top_h"]) if params.get("eil_top_h") is not None else None,
    }


@app.route("/api/vad", methods=["GET"])
def get_vad():
    """
    Fetch latest NEXRAD VAD Wind Profile for a given radar.
    Query params: ?radar=KTLX
    Returns JSON with wind data at various altitudes.
    """
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


@app.route("/api/vwp-display", methods=["GET"])
def vwp_display():
    """
    Generate a VWP (VAD Wind Profile) time-height display image.
    Query params: ?radar=KTLX&hours=12
    Returns JSON with base64 PNG image.
    """
    radar = request.args.get("radar", "").upper()
    hours = int(request.args.get("hours", 12))
    if not radar or len(radar) < 3:
        return jsonify({"error": "Provide a valid NEXRAD radar ID (e.g. ?radar=KTLX)"}), 400
    hours = max(1, min(hours, 48))  # clamp to 1-48
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
            }
        })
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


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


# ── Merge Profiles endpoint ─────────────────────────────────────────
@app.route("/api/merge-profiles", methods=["POST", "OPTIONS"])
def merge_profiles_endpoint():
    """
    Merge two soundings into a weighted-average blended profile.

    Expected JSON body:
      {
        "soundings": [
          { "source": "obs", "station": "OUN", "date": "2024061200" },
          { "source": "obs", "station": "FWD", "date": "2024061200" }
        ],
        "weight": 0.5,
        "theme": "dark",
        "colorblind": false
      }
    weight = proportion of profile A (0–1).  Default 0.5.
    """
    if request.method == "OPTIONS":
        return "", 204

    body = request.get_json(force=True)
    items = body.get("soundings", [])
    weight = float(body.get("weight", 0.5))
    theme = body.get("theme", "dark")
    cb = body.get("colorblind", False)

    if len(items) != 2:
        return jsonify({"error": "Exactly 2 soundings required for merging."}), 400

    try:
        datasets = []
        labels = []
        for item in items:
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
            datasets.append(data)
            labels.append(f"{stn or src.upper()} {dt.strftime('%d/%HZ')}")

        # Merge
        merged_data = merge_profiles(datasets[0], datasets[1], weight_a=weight)
        merged_params = compute_parameters(merged_data)

        merge_label = f"MERGE: {int(weight*100)}% {labels[0]} + {int((1-weight)*100)}% {labels[1]}"
        plot_id = "MERGED"
        fig = plot_sounding(merged_data, merged_params, plot_id, dt, theme=theme,
                            colorblind=cb)

        _facecolor = "#f5f5f5" if theme == "light" else "#0d0d0d"
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=180, facecolor=_facecolor)
        plt.close(fig)
        buf.seek(0)
        image_b64 = base64.b64encode(buf.read()).decode("utf-8")

        serialized = _serialize_params(merged_params, merged_data, "MERGED", dt, "merge")

        return jsonify({
            "image": image_b64,
            "params": serialized,
            "label": merge_label,
            "meta": {
                "station": "MERGED",
                "stationName": merge_label,
                "source": "merge",
                "date": dt.strftime("%Y-%m-%d %HZ"),
                "levels": len(merged_data["pressure"]),
                "weight": weight,
                "profileA": labels[0],
                "profileB": labels[1],
            },
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
