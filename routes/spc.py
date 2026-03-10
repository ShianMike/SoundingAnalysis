"""
SPC Convective Outlook proxy routes.
"""
import re
import time

import requests as _requests
from flask import Blueprint, jsonify, request

from sounding import STATIONS

bp = Blueprint("spc", __name__)

# Keyed by "{day}_{type}"
# Types: cat, torn, wind, hail (probabilistic — includes embedded CI groups)
_SPC_OUTLOOK_URLS = {
    # Categorical
    "1_cat": "https://www.spc.noaa.gov/products/outlook/day1otlk_cat.lyr.geojson",
    "2_cat": "https://www.spc.noaa.gov/products/outlook/day2otlk_cat.lyr.geojson",
    "3_cat": "https://www.spc.noaa.gov/products/outlook/day3otlk_cat.lyr.geojson",
    # Probabilistic (Day 1-2 only; Day 3 has no prob outlooks)
    "1_torn": "https://www.spc.noaa.gov/products/outlook/day1otlk_torn.lyr.geojson",
    "2_torn": "https://www.spc.noaa.gov/products/outlook/day2otlk_torn.lyr.geojson",
    "1_wind": "https://www.spc.noaa.gov/products/outlook/day1otlk_wind.lyr.geojson",
    "2_wind": "https://www.spc.noaa.gov/products/outlook/day2otlk_wind.lyr.geojson",
    "1_hail": "https://www.spc.noaa.gov/products/outlook/day1otlk_hail.lyr.geojson",
    "2_hail": "https://www.spc.noaa.gov/products/outlook/day2otlk_hail.lyr.geojson",
    # Day 4-8 extended-range (probability-based, single product per day)
    "4_prob": "https://www.spc.noaa.gov/products/exper/day4-8/day4prob.lyr.geojson",
    "5_prob": "https://www.spc.noaa.gov/products/exper/day4-8/day5prob.lyr.geojson",
    "6_prob": "https://www.spc.noaa.gov/products/exper/day4-8/day6prob.lyr.geojson",
    "7_prob": "https://www.spc.noaa.gov/products/exper/day4-8/day7prob.lyr.geojson",
    "8_prob": "https://www.spc.noaa.gov/products/exper/day4-8/day8prob.lyr.geojson",
}

_VALID_DAYS = {"1", "2", "3", "4", "5", "6", "7", "8"}

# SPC discussion text URLs (plain-text ACUS products)
_SPC_DISCUSSION_URLS = {
    "1": "https://www.spc.noaa.gov/products/outlook/day1otlk.txt",
    "2": "https://www.spc.noaa.gov/products/outlook/day2otlk.txt",
    "3": "https://www.spc.noaa.gov/products/outlook/day3otlk.txt",
    "48": "https://www.spc.noaa.gov/products/exper/day4-8/day48prob.txt",
}

_spc_cache = {}
_SPC_CACHE_TTL = 600  # 10 minutes


@bp.route("/api/spc-outlook", methods=["GET"])
def spc_outlook():
    """Proxy SPC outlook GeoJSON."""
    day = request.args.get("day", "1")
    otype = request.args.get("type", "cat")
    if day not in _VALID_DAYS:
        return jsonify({"error": f"Invalid day: {day}. Use 1-8."}), 400
    # Day 4-8 only supports "prob" type
    if int(day) >= 4:
        otype = "prob"
    valid_types = ("cat", "torn", "wind", "hail", "prob")
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


def _parse_discussion(raw_text, day):
    """Extract the narrative discussion from an SPC outlook text product."""
    # The discussion block sits between ...SUMMARY... or ...DISCUSSION...
    # and the next section marker or the end of the narrative area.
    text = raw_text.replace("\r\n", "\n")

    # Try to extract between "...DISCUSSION..." and the end marker ($$, .PREV, or similar)
    disc_match = re.search(
        r"\.{3}(?:SUMMARY|DISCUSSION)\.{3}\s*\n(.*?)(?:\n\$\$|\n\.PREV|\n&&|\Z)",
        text, re.DOTALL | re.IGNORECASE,
    )
    if disc_match:
        body = disc_match.group(1).strip()
    else:
        # Fallback: grab everything after the header lines
        body = text.strip()

    return body


@bp.route("/api/spc-discussion", methods=["GET"])
def spc_discussion():
    """Proxy SPC convective outlook discussion text."""
    day = request.args.get("day", "1")
    if day not in _VALID_DAYS:
        return jsonify({"error": f"Invalid day: {day}. Use 1-8."}), 400

    # Day 4-8 all share the same product
    url_key = "48" if int(day) >= 4 else day
    cache_key = f"disc_{url_key}"

    cached = _spc_cache.get(cache_key)
    if cached and (time.time() - cached["ts"]) < _SPC_CACHE_TTL:
        return jsonify(cached["data"])

    url = _SPC_DISCUSSION_URLS.get(url_key)
    if not url:
        return jsonify({"error": f"No discussion product for day {day}"}), 400

    try:
        resp = _requests.get(url, timeout=15)
        resp.raise_for_status()
        raw = resp.text
        body = _parse_discussion(raw, int(day))
        result = {"day": int(day), "text": body}
        _spc_cache[cache_key] = {"data": result, "ts": time.time()}
        return jsonify(result)
    except Exception as e:
        print(f"[SPC] Failed to fetch Day {day} discussion: {e}")
        return jsonify({"error": f"Failed to fetch discussion: {e}"}), 502


@bp.route("/api/spc-outlook-stations", methods=["GET"])
def spc_outlook_stations():
    """Return sounding stations within SPC outlook polygons."""
    day = request.args.get("day", "1")
    otype = request.args.get("type", "cat")
    min_risk = request.args.get("minRisk", "")

    valid_types = ("cat", "torn", "wind", "hail", "prob")
    if day not in _VALID_DAYS:
        return jsonify({"error": f"Invalid day: {day}. Use 1-8."}), 400
    # Day 4-8 only supports "prob" type
    if int(day) >= 4:
        otype = "prob"
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
