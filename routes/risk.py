"""
Risk-scan and time-series routes.
"""
import os
import time
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta, timezone

from flask import Blueprint, jsonify, request

from sounding import (
    STATIONS, STATION_WMO, BUFKIT_MODELS, TORNADO_SCAN_STATIONS,
    fetch_sounding, compute_parameters,
    get_latest_sounding_time, _quick_tornado_score, _quick_forecast_score,
)
from .helpers import _serialize_params, _nan_safe

bp = Blueprint("risk", __name__)


@bp.route("/api/risk-scan", methods=["POST", "OPTIONS"])
def risk_scan():
    """Scan stations and return risk scores."""
    if request.method == "OPTIONS":
        return "", 204

    body = request.get_json(force=True) if request.data else {}
    date_str = body.get("date")

    if date_str:
        try:
            dt = datetime.strptime(str(date_str), "%Y%m%d%H").replace(tzinfo=timezone.utc)
        except ValueError:
            return jsonify({"error": f"Invalid date format '{date_str}'."}), 400
    else:
        dt = get_latest_sounding_time()

    def _scan_one(sid):
        try:
            score = _quick_tornado_score(sid, dt)
            if score is None:
                return None
            stp_score, raw_score, cape, srh, bwd, scp, ship, dcp, bwd_dir = score
            result = {
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
            if bwd_dir is not None:
                result["bwdDir"] = round(bwd_dir)
            return result
        except Exception:
            return None

    results = []
    station_ids = list(STATIONS.keys())
    scan_workers = int(os.environ.get("SCAN_WORKERS", "20"))

    try:
        with ThreadPoolExecutor(max_workers=scan_workers) as executor:
            futures = {executor.submit(_scan_one, sid): sid for sid in station_ids}
            for future in as_completed(futures, timeout=90):
                try:
                    r = future.result(timeout=15)
                    if r is not None:
                        results.append(r)
                except Exception:
                    pass
    except Exception as exc:
        print(f"[risk-scan] partial results due to: {exc}")

    results.sort(key=lambda r: (r["stp"], r["raw"]), reverse=True)

    return jsonify({
        "date": dt.strftime("%Y-%m-%d %HZ"),
        "stations": results,
    })


# ─── Forecast Risk Scan ────────────────────────────────────────────
@bp.route("/api/forecast-risk-scan", methods=["POST", "OPTIONS"])
def forecast_risk_scan():
    """Scan stations using model/forecast data and return risk scores."""
    if request.method == "OPTIONS":
        return "", 204

    body = request.get_json(force=True) if request.data else {}
    model = body.get("model", "hrrr").lower()
    fhour = int(body.get("fhour", 0))

    if model not in BUFKIT_MODELS:
        return jsonify({"error": f"Unknown model '{model}'."}), 400

    # Model init cycles and data availability lag
    # HRRR/RAP run hourly but archive has ~2-3h delay
    # NAM/GFS run at 00/06/12/18Z with ~4-5h delay
    MODEL_CYCLES = {
        "hrrr":    {"interval": 1,  "lag": 3},
        "rap":     {"interval": 1,  "lag": 3},
        "nam":     {"interval": 6,  "lag": 5},
        "namnest": {"interval": 6,  "lag": 5},
        "gfs":     {"interval": 6,  "lag": 6},
        "sref":    {"interval": 6,  "lag": 6},
    }
    cycle_info = MODEL_CYCLES.get(model, {"interval": 6, "lag": 5})
    now = datetime.now(timezone.utc)
    # Find most recent init cycle that should be available
    candidate = now - timedelta(hours=cycle_info["lag"])
    if cycle_info["interval"] == 1:
        dt = candidate.replace(minute=0, second=0, microsecond=0)
    else:
        cycle_hour = (candidate.hour // cycle_info["interval"]) * cycle_info["interval"]
        dt = candidate.replace(hour=cycle_hour, minute=0, second=0, microsecond=0)

    # Valid time = init time + forecast hour
    valid_time = dt + timedelta(hours=fhour)

    station_ids = list(TORNADO_SCAN_STATIONS)

    def _scan_one(sid):
        try:
            score = _quick_forecast_score(sid, dt, model=model, fhour=fhour)
            if score is None:
                return None
            stp_score, raw_score, cape, srh, bwd, scp, ship, dcp, bwd_dir = score
            result = {
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
            if bwd_dir is not None:
                result["bwdDir"] = round(bwd_dir)
            return result
        except Exception as e:
            print(f"[forecast-scan] {sid} failed: {e}")
            return None

    results = []
    scan_workers = int(os.environ.get("SCAN_WORKERS", "20"))
    print(f"[forecast-risk-scan] model={model} init={dt:%Y-%m-%d %H}Z fhour={fhour} "
          f"valid={valid_time:%Y-%m-%d %H}Z stations={len(station_ids)}")

    try:
        with ThreadPoolExecutor(max_workers=scan_workers) as executor:
            futures = {executor.submit(_scan_one, sid): sid for sid in station_ids}
            for future in as_completed(futures, timeout=90):
                try:
                    r = future.result(timeout=15)
                    if r is not None:
                        results.append(r)
                except Exception:
                    pass
    except Exception as exc:
        print(f"[forecast-risk-scan] partial results due to: {exc}")

    results.sort(key=lambda r: (r["stp"], r["raw"]), reverse=True)

    return jsonify({
        "date": valid_time.strftime("%Y-%m-%d %HZ"),
        "model": model.upper(),
        "fhour": fhour,
        "initTime": dt.strftime("%Y-%m-%d %HZ"),
        "stations": results,
    })


# ─── Time-Series ────────────────────────────────────────────────────
@bp.route("/api/time-series", methods=["POST", "OPTIONS"])
def time_series():
    """Fetch sounding parameters for a station across a date range."""
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

        if source == "obs" and (not station or station not in STATION_WMO):
            return jsonify({"error": f"Unknown station '{station}'."}), 400

        if source in ("rap", "era5") and lat is None and lon is None:
            if station and station in STATIONS:
                _, lat, lon = STATIONS[station]
            else:
                return jsonify({"error": "RAP/ERA5 requires lat/lon or a known station."}), 400

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

        if (end_dt - start_dt).days > 14:
            start_dt = end_dt - timedelta(days=14)

        span_days = (end_dt - start_dt).days

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

        def _fetch_one(dt_val):
            try:
                data = fetch_sounding(
                    station_id=station, dt=dt_val, source=source,
                    lat=lat, lon=lon, model=model, fhour=fhour,
                )
                params = compute_parameters(data)
                return {
                    "date": dt_val.strftime("%Y-%m-%d %HZ"),
                    "ts": dt_val.isoformat(),
                    "params": _serialize_params(params, data, station, dt_val, source),
                }
            except Exception:
                return None

        points = []
        wall_start = time.monotonic()
        WALL_LIMIT = int(os.environ.get("TS_WALL_LIMIT", "85"))
        ts_workers = int(os.environ.get("TS_WORKERS", "2"))

        with ThreadPoolExecutor(max_workers=ts_workers) as executor:
            futures = {executor.submit(_fetch_one, t): t for t in times}
            for future in as_completed(futures, timeout=80):
                if time.monotonic() - wall_start > WALL_LIMIT:
                    break
                try:
                    r = future.result(timeout=20)
                    if r is not None:
                        points.append(r)
                except Exception:
                    pass

        points.sort(key=lambda p: p["ts"])

        station_name = STATIONS.get(station, (station or source.upper(),))[0]

        return jsonify(_nan_safe({
            "station": station or f"{lat},{lon}",
            "stationName": station_name,
            "source": source,
            "startDate": start_dt.strftime("%Y-%m-%d"),
            "endDate": end_dt.strftime("%Y-%m-%d %HZ"),
            "resolution": "00Z & 12Z" if span_days <= 7 else "12Z only",
            "points": points,
        }))

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500
