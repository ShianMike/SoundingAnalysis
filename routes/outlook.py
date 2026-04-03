"""
Custom Severe Weather Outlook — formula-based threat classification.

Uses latest observed radiosonde data to classify stations into
SPC-style categorical and hazard-specific risk levels using our own formulas.
Returns GeoJSON polygons for map rendering.
"""
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta, timezone
from functools import lru_cache

from flask import Blueprint, jsonify, request
from shapely.geometry import MultiPoint, mapping
from shapely.ops import unary_union

from sounding import (
    STATIONS, TORNADO_SCAN_STATIONS,
    _quick_tornado_score,
    get_latest_sounding_time,
)

bp = Blueprint("outlook", __name__)

# ── Result cache keyed on sounding cycle (e.g. "2026-04-03T00Z") ────
# Avoids re-scanning when user refreshes within the same cycle.
_outlook_cache = {}  # { cycle_key: { "result": ..., "ts": float } }
_CACHE_TTL = 600     # 10 minutes

# ── Category constants ──────────────────────────────────────────
CAT_LEVELS = ["NONE", "TSTM", "MRGL", "SLGT", "ENH", "MDT", "HIGH"]
CAT_ORDER = {k: i for i, k in enumerate(CAT_LEVELS)}


# ── Classification helpers ──────────────────────────────────────

def _classify_categorical(stp, scp, ship, dcp, cape, srh, bwd):
    """Overall categorical risk — mirrors SPC Day-1 style tiers."""
    if stp >= 8 and scp >= 8 and cape >= 2500:
        return "HIGH"
    if (stp >= 4 and scp >= 4) or \
       (stp >= 3 and cape >= 2000 and srh >= 200 and bwd >= 40):
        return "MDT"
    if stp >= 2 or scp >= 4 or ship >= 2 or \
       (cape >= 1500 and srh >= 150 and bwd >= 35):
        return "ENH"
    if stp >= 1 or scp >= 2 or ship >= 1 or \
       (cape >= 1000 and srh >= 75 and bwd >= 30):
        return "SLGT"
    if stp >= 0.5 or scp >= 1 or ship >= 0.5 or \
       (cape >= 500 and bwd >= 25):
        return "MRGL"
    if cape >= 250:
        return "TSTM"
    return "NONE"


def _classify_tornado(stp, srh, bwd, cape):
    if stp >= 6 and srh >= 300:
        return "EXTREME"
    if stp >= 3 and srh >= 200 and bwd >= 35:
        return "HIGH"
    if stp >= 1 and srh >= 100:
        return "MOD"
    if stp >= 0.5:
        return "LOW"
    return "NONE"


def _classify_wind(dcp, bwd, cape):
    if dcp >= 3 and bwd >= 40 and cape >= 1000:
        return "HIGH"
    if (dcp >= 1.5 and bwd >= 30) or (bwd >= 40 and cape >= 1000):
        return "MOD"
    if (dcp >= 0.5 and bwd >= 30) or (bwd >= 35 and cape >= 500):
        return "LOW"
    return "NONE"


def _classify_hail(ship, bwd, cape):
    if ship >= 2 and bwd >= 40 and cape >= 1500:
        return "HIGH"
    if ship >= 1 and bwd >= 30:
        return "MOD"
    if ship >= 0.5:
        return "LOW"
    return "NONE"


def _max_cat(*cats):
    """Return the highest categorical level from a list."""
    best = "NONE"
    for c in cats:
        if CAT_ORDER.get(c, 0) > CAT_ORDER.get(best, 0):
            best = c
    return best


# ── Polygon generation ──────────────────────────────────────────

# Graduated buffer per category (degrees); lower categories get wider radii
# so visual nesting is guaranteed even when the same stations span two levels.
_CAT_BUFFER = {
    "TSTM": 2.0,   # ~220 km — broad thunderstorm area
    "MRGL": 1.7,   # ~190 km
    "SLGT": 1.4,   # ~155 km
    "ENH":  1.2,   # ~130 km
    "MDT":  1.0,   # ~110 km
    "HIGH": 0.8,   # ~90 km
}

_HAZARD_BUFFER = {
    "LOW":     1.8,
    "MOD":     1.4,
    "HIGH":    1.1,
    "EXTREME": 0.9,
}

# Map internal hazard level names → SPC-style probability labels
_TORNADO_LEVELS = [
    ("LOW",     "0.02"),
    ("MOD",     "0.05"),
    ("HIGH",    "0.10"),
    ("EXTREME", "0.30"),
]

_WIND_LEVELS = [
    ("LOW",  "0.05"),
    ("MOD",  "0.15"),
    ("HIGH", "0.30"),
]

_HAIL_LEVELS = [
    ("LOW",  "0.05"),
    ("MOD",  "0.15"),
    ("HIGH", "0.30"),
]


def _stations_to_polygon(pts, buffer_deg):
    """Convex hull of station cluster, then buffer for smooth edges."""
    mp = MultiPoint(pts)
    if len(pts) == 1:
        geom = mp.buffer(buffer_deg, resolution=32)
    elif len(pts) == 2:
        geom = mp.convex_hull.buffer(buffer_deg, resolution=32)
    else:
        geom = mp.convex_hull.buffer(buffer_deg, resolution=32)
    return geom.simplify(0.08, preserve_topology=True)


def _build_geojson(stations, valid_iso):
    """
    Build a GeoJSON FeatureCollection with blob polygons for each
    categorical risk level.  Each level encompasses all stations at
    that level *or above* so polygons nest naturally.
    """
    render_levels = [lv for lv in CAT_LEVELS if lv != "NONE"]
    features = []
    for level in render_levels:
        min_order = CAT_ORDER[level]
        pts = [
            (s["lon"], s["lat"])
            for s in stations
            if CAT_ORDER.get(s["categorical"], 0) >= min_order
        ]
        if not pts:
            continue

        buf = _CAT_BUFFER.get(level, 1.5)
        geom = _stations_to_polygon(pts, buf)

        features.append({
            "type": "Feature",
            "properties": {
                "LABEL": level,
                "VALID_ISO": valid_iso,
                "station_count": len(pts),
            },
            "geometry": mapping(geom),
        })

    return {"type": "FeatureCollection", "features": features}


def _build_hazard_geojson(stations, hazard_key, level_map, valid_iso):
    """
    Build GeoJSON for a probabilistic hazard type (tornado / wind / hail).
    level_map is a list of (internal_level, spc_prob_label).
    """
    level_order = {name: i for i, (name, _) in enumerate(level_map)}
    level_order["NONE"] = -1

    features = []
    for level_name, prob_label in level_map:
        min_order = level_order[level_name]
        pts = [
            (s["lon"], s["lat"])
            for s in stations
            if level_order.get(s[hazard_key], -1) >= min_order
        ]
        if not pts:
            continue

        buf = _HAZARD_BUFFER.get(level_name, 1.5)
        geom = _stations_to_polygon(pts, buf)

        features.append({
            "type": "Feature",
            "properties": {
                "LABEL": prob_label,
                "VALID_ISO": valid_iso,
                "station_count": len(pts),
            },
            "geometry": mapping(geom),
        })

    return {"type": "FeatureCollection", "features": features}


# ── Endpoint ────────────────────────────────────────────────────

@bp.route("/api/outlook", methods=["POST", "OPTIONS"])
def generate_outlook():
    """Generate a custom severe-weather outlook from latest observed soundings."""
    if request.method == "OPTIONS":
        return "", 204

    # Always use latest observed radiosonde data
    dt = get_latest_sounding_time()
    valid_time = dt
    cycle_key = dt.strftime("%Y%m%dT%HZ")

    # Return cached result if still fresh
    cached = _outlook_cache.get(cycle_key)
    if cached and (time.time() - cached["ts"]) < _CACHE_TTL:
        print(f"[outlook] cache hit for {cycle_key}")
        return jsonify(cached["result"])

    station_ids = list(STATIONS.keys())

    # High parallelism — these are I/O-bound HTTP fetches
    scan_workers = int(os.environ.get("SCAN_WORKERS", "20"))
    raw_results = []

    t0 = time.time()
    try:
        with ThreadPoolExecutor(max_workers=scan_workers) as pool:
            futs = {pool.submit(_quick_tornado_score, sid, dt): sid
                    for sid in station_ids}
            for fut in as_completed(futs, timeout=90):
                try:
                    score = fut.result(timeout=12)
                    if score is None:
                        continue
                    sid = futs[fut]
                    stp, raw, cape, srh, bwd, scp, ship, dcp, bwd_dir = score
                    raw_results.append({
                        "id": sid,
                        "name": STATIONS.get(sid, (sid,))[0],
                        "lat": STATIONS.get(sid, ("", 0, 0))[1],
                        "lon": STATIONS.get(sid, ("", 0, 0))[2],
                        "stp": round(stp, 2),
                        "scp": round(scp, 2),
                        "ship": round(ship, 2),
                        "dcp": round(dcp, 2),
                        "cape": round(cape),
                        "srh": round(srh),
                        "bwd": round(bwd),
                    })
                except Exception:
                    pass
    except Exception as exc:
        print(f"[outlook] partial results: {exc}")

    elapsed = time.time() - t0
    print(f"[outlook] scanned {len(raw_results)}/{len(station_ids)} "
          f"stations in {elapsed:.1f}s")

    # Classify each station
    for s in raw_results:
        s["categorical"] = _classify_categorical(
            s["stp"], s["scp"], s["ship"], s["dcp"],
            s["cape"], s["srh"], s["bwd"],
        )
        s["tornado"] = _classify_tornado(s["stp"], s["srh"], s["bwd"], s["cape"])
        s["wind"] = _classify_wind(s["dcp"], s["bwd"], s["cape"])
        s["hail"] = _classify_hail(s["ship"], s["bwd"], s["cape"])

    # Sort by categorical severity descending
    raw_results.sort(key=lambda r: (
        CAT_ORDER.get(r["categorical"], 0),
        r["stp"], r["scp"],
    ), reverse=True)

    # Build GeoJSON polygons for map rendering
    valid_iso = valid_time.strftime("%Y-%m-%dT%H:%M:%SZ")
    geojson = _build_geojson(raw_results, valid_iso)
    tornado_geo = _build_hazard_geojson(raw_results, "tornado",
                                        _TORNADO_LEVELS, valid_iso)
    wind_geo = _build_hazard_geojson(raw_results, "wind",
                                     _WIND_LEVELS, valid_iso)
    hail_geo = _build_hazard_geojson(raw_results, "hail",
                                     _HAIL_LEVELS, valid_iso)

    result = {
        "geojson": geojson,
        "tornadoGeo": tornado_geo,
        "windGeo": wind_geo,
        "hailGeo": hail_geo,
        "validTime": valid_time.strftime("%Y-%m-%d %HZ"),
        "stationCount": len(raw_results),
    }

    # Cache for this sounding cycle
    _outlook_cache[cycle_key] = {"result": result, "ts": time.time()}
    # Evict old cycles
    for k in list(_outlook_cache):
        if k != cycle_key:
            del _outlook_cache[k]

    return jsonify(result)
