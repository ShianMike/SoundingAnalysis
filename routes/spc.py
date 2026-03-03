"""
SPC Convective Outlook proxy routes.
"""
import time

import requests as _requests
from flask import Blueprint, jsonify, request

from sounding import STATIONS

bp = Blueprint("spc", __name__)

# Keyed by "{day}_{type}" — type: cat, ci_torn, ci_wind, ci_hail
_SPC_OUTLOOK_URLS = {
    "1_cat": "https://www.spc.noaa.gov/products/outlook/day1otlk_cat.lyr.geojson",
    "2_cat": "https://www.spc.noaa.gov/products/outlook/day2otlk_cat.lyr.geojson",
    "3_cat": "https://www.spc.noaa.gov/products/outlook/day3otlk_cat.lyr.geojson",
    "1_ci_torn": "https://www.spc.noaa.gov/products/outlook/day1otlk_CI_torn.lyr.geojson",
    "2_ci_torn": "https://www.spc.noaa.gov/products/outlook/day2otlk_CI_torn.lyr.geojson",
    "3_ci_torn": "https://www.spc.noaa.gov/products/outlook/day3otlk_CI_torn.lyr.geojson",
    "1_ci_wind": "https://www.spc.noaa.gov/products/outlook/day1otlk_CI_wind.lyr.geojson",
    "2_ci_wind": "https://www.spc.noaa.gov/products/outlook/day2otlk_CI_wind.lyr.geojson",
    "3_ci_wind": "https://www.spc.noaa.gov/products/outlook/day3otlk_CI_wind.lyr.geojson",
    "1_ci_hail": "https://www.spc.noaa.gov/products/outlook/day1otlk_CI_hail.lyr.geojson",
    "2_ci_hail": "https://www.spc.noaa.gov/products/outlook/day2otlk_CI_hail.lyr.geojson",
    "3_ci_hail": "https://www.spc.noaa.gov/products/outlook/day3otlk_CI_hail.lyr.geojson",
}

_spc_cache = {}
_SPC_CACHE_TTL = 600  # 10 minutes


@bp.route("/api/spc-outlook", methods=["GET"])
def spc_outlook():
    """Proxy SPC outlook GeoJSON."""
    day = request.args.get("day", "1")
    otype = request.args.get("type", "cat")
    if day not in ("1", "2", "3"):
        return jsonify({"error": f"Invalid day: {day}. Use 1, 2, or 3."}), 400
    valid_types = ("cat", "ci_torn", "ci_wind", "ci_hail")
    if otype not in valid_types:
        return jsonify({"error": f"Invalid type: {otype}. Use one of: {', '.join(valid_types)}"}), 400

    cache_key = f"{day}_{otype}"
    if cache_key not in _SPC_OUTLOOK_URLS:
        return jsonify({"error": f"No SPC product for day={day} type={otype}"}), 400

    cached = _spc_cache.get(cache_key)
    if cached and (time.time() - cached["ts"]) < _SPC_CACHE_TTL:
        return jsonify(cached["data"])

    url = _SPC_OUTLOOK_URLS[cache_key]
    try:
        resp = _requests.get(url, timeout=15)
        resp.raise_for_status()
        geojson = resp.json()
        _spc_cache[cache_key] = {"data": geojson, "ts": time.time()}
        return jsonify(geojson)
    except Exception as e:
        print(f"[SPC] Failed to fetch Day {day} {otype} outlook: {e}")
        return jsonify({"error": f"Failed to fetch SPC outlook: {e}"}), 502


@bp.route("/api/spc-outlook-stations", methods=["GET"])
def spc_outlook_stations():
    """Return sounding stations within SPC outlook polygons."""
    day = request.args.get("day", "1")
    otype = request.args.get("type", "cat")
    min_risk = request.args.get("minRisk", "")

    valid_types = ("cat", "ci_torn", "ci_wind", "ci_hail")
    if day not in ("1", "2", "3"):
        return jsonify({"error": f"Invalid day: {day}"}), 400
    if otype not in valid_types:
        return jsonify({"error": f"Invalid type: {otype}"}), 400

    cache_key = f"{day}_{otype}"
    if cache_key not in _SPC_OUTLOOK_URLS:
        return jsonify({"error": f"No SPC product for day={day} type={otype}"}), 400

    cached = _spc_cache.get(cache_key)
    if cached and (time.time() - cached["ts"]) < _SPC_CACHE_TTL:
        geojson = cached["data"]
    else:
        url = _SPC_OUTLOOK_URLS[cache_key]
        try:
            resp = _requests.get(url, timeout=15)
            resp.raise_for_status()
            geojson = resp.json()
            _spc_cache[cache_key] = {"data": geojson, "ts": time.time()}
        except Exception as e:
            return jsonify({"error": f"Failed to fetch outlook: {e}"}), 502

    try:
        from shapely.geometry import shape, Point
    except ImportError:
        return _spc_stations_bbox_fallback(geojson, min_risk)

    CAT_ORDER = ["TSTM", "MRGL", "SLGT", "ENH", "MDT", "HIGH"]
    min_idx = CAT_ORDER.index(min_risk) if min_risk in CAT_ORDER else 0

    results = []
    seen = set()
    features = geojson.get("features", [])

    for feat in features:
        label = feat.get("properties", {}).get("LABEL", feat.get("properties", {}).get("LABEL2", ""))
        if otype == "cat" and min_risk:
            if label in CAT_ORDER and CAT_ORDER.index(label) < min_idx:
                continue

        try:
            poly = shape(feat["geometry"])
        except Exception:
            continue

        for sid, (name, lat, lon) in STATIONS.items():
            if sid in seen:
                continue
            pt = Point(lon, lat)
            try:
                if poly.contains(pt) or poly.touches(pt):
                    seen.add(sid)
                    results.append({
                        "id": sid, "name": name, "lat": lat, "lon": lon,
                        "riskLabel": label,
                    })
            except Exception:
                pass

    results.sort(key=lambda r: r["id"])
    return jsonify({
        "day": int(day),
        "type": otype,
        "stations": results,
        "count": len(results),
    })


def _spc_stations_bbox_fallback(geojson, min_risk):
    """Fallback when shapely is not installed — uses bounding-box approximation."""
    results = []
    seen = set()
    for feat in geojson.get("features", []):
        label = feat.get("properties", {}).get("LABEL", "")
        geom = feat.get("geometry", {})
        coords = geom.get("coordinates", [])
        if not coords:
            continue
        all_lons, all_lats = [], []

        def _flatten(c):
            if isinstance(c, (list, tuple)) and len(c) >= 2 and isinstance(c[0], (int, float)):
                all_lons.append(c[0])
                all_lats.append(c[1])
            elif isinstance(c, (list, tuple)):
                for sub in c:
                    _flatten(sub)
        _flatten(coords)
        if not all_lons:
            continue
        min_lon, max_lon = min(all_lons), max(all_lons)
        min_lat, max_lat = min(all_lats), max(all_lats)

        for sid, (name, lat, lon) in STATIONS.items():
            if sid in seen:
                continue
            if min_lat <= lat <= max_lat and min_lon <= lon <= max_lon:
                seen.add(sid)
                results.append({
                    "id": sid, "name": name, "lat": lat, "lon": lon,
                    "riskLabel": label,
                })

    results.sort(key=lambda r: r["id"])
    return jsonify({"day": 0, "type": "cat", "stations": results, "count": len(results)})
