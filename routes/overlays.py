"""
Backend proxy routes for map overlays.

These proxy external APIs (NWS, SPC, IEM, Spotter Network) so the
browser never makes cross-origin requests that would be blocked by
CORS or Content Security Policy.  Each endpoint is cached briefly to
avoid hammering upstream services.
"""

import time
import re

import requests as _requests
from flask import Blueprint, jsonify

bp = Blueprint("overlays", __name__)

# ── Cache helpers ─────────────────────────────────────────────────
_cache = {}


def _cached(key, ttl):
    """Return cached (data, True) if fresh, else (None, False)."""
    entry = _cache.get(key)
    if entry and (time.time() - entry["ts"]) < ttl:
        return entry["data"], True
    return None, False


def _store(key, data):
    _cache[key] = {"data": data, "ts": time.time()}


_HEADERS = {"User-Agent": "SoundingAnalysis/1.0 (proxy)"}

# ── 1. NWS Active Warnings ───────────────────────────────────────
NWS_ALERTS_URL = "https://api.weather.gov/alerts/active?status=actual"
NWS_CACHE_TTL = 60  # 1 min


@bp.route("/api/overlays/warnings", methods=["GET"])
def nws_warnings():
    """Proxy NWS active warnings (GeoJSON)."""
    cached, hit = _cached("warnings", NWS_CACHE_TTL)
    if hit:
        return jsonify(cached)

    try:
        resp = _requests.get(
            NWS_ALERTS_URL,
            headers={**_HEADERS, "Accept": "application/geo+json"},
            timeout=15,
        )
        resp.raise_for_status()
        data = resp.json()
        _store("warnings", data)
        return jsonify(data)
    except Exception as e:
        print(f"[overlays] NWS warnings fetch error: {e}")
        return jsonify({"type": "FeatureCollection", "features": []}), 502


# ── 2. TVS / Mesocyclone Storm Attributes (IEM) ──────────────────
IEM_STORM_ATTR_URL = (
    "https://mesonet.agron.iastate.edu/geojson/nexrad_attr.geojson"
)
STORM_CACHE_TTL = 30  # 30s — data changes every volume scan


@bp.route("/api/overlays/storm-attributes", methods=["GET"])
def storm_attributes():
    """Proxy IEM NEXRAD storm attributes (TVS / Mesocyclone) GeoJSON."""
    cached, hit = _cached("storm_attr", STORM_CACHE_TTL)
    if hit:
        return jsonify(cached)

    try:
        resp = _requests.get(
            IEM_STORM_ATTR_URL, headers=_HEADERS, timeout=15
        )
        resp.raise_for_status()
        data = resp.json()
        _store("storm_attr", data)
        return jsonify(data)
    except Exception as e:
        print(f"[overlays] Storm attr fetch error: {e}")
        return jsonify({"type": "FeatureCollection", "features": []}), 502


# ── 3. SPC Mesoscale Discussions ──────────────────────────────────
# The old IEM spc_mcd.geojson endpoint was removed.
# We now use the NWS MapServer with a phenom query for Mesoscale Discussions.
# Layer 1 (WatchesWarnings) contains MDs when queried with phenom='SV' sig='Y' (SPS)
# but MDs are actually issued by SPC separately.  The most reliable public
# source is the SPC's own GeoJSON served from their products page.

SPC_MD_URL = (
    "https://www.spc.noaa.gov/products/md/md.geojson"
)

# Fallback: scrape the MDs page and build minimal GeoJSON from the MD
# polygon coordinates embedded in each individual MD product page.
SPC_MD_FALLBACK_URL = (
    "https://mapservices.weather.noaa.gov/eventdriven/rest/services/WWA/"
    "watch_warn_adv/MapServer/1/query"
    "?where=prod_type+LIKE+%27%25Mesoscale%25%27"
    "&f=geojson&outFields=prod_type,event,wfo,expiration,issuance"
    "&returnGeometry=true&outSR=4326"
)

MD_CACHE_TTL = 120  # 2 min


@bp.route("/api/overlays/spc-mds", methods=["GET"])
def spc_mds():
    """Proxy SPC Mesoscale Discussions (GeoJSON)."""
    cached, hit = _cached("spc_mds", MD_CACHE_TTL)
    if hit:
        return jsonify(cached)

    # Try the primary SPC GeoJSON feed
    data = _try_fetch_geojson(SPC_MD_URL)

    # Fallback: NWS MapServer query for MDs
    if data is None:
        data = _try_fetch_geojson(SPC_MD_FALLBACK_URL)

    if data is None:
        data = {"type": "FeatureCollection", "features": []}

    _store("spc_mds", data)
    return jsonify(data)


def _try_fetch_geojson(url):
    """Attempt to fetch GeoJSON from url, return dict or None."""
    try:
        resp = _requests.get(url, headers=_HEADERS, timeout=15)
        resp.raise_for_status()
        data = resp.json()
        if "features" in data:
            return data
    except Exception as e:
        print(f"[overlays] GeoJSON fetch failed ({url[:60]}…): {e}")
    return None


# ── 4. SPC Watches (NWS MapServer) ───────────────────────────────
SPC_WATCH_URL = (
    "https://mapservices.weather.noaa.gov/eventdriven/rest/services/WWA/"
    "watch_warn_adv/MapServer/1/query"
    "?where=(phenom%3D%27TO%27+OR+phenom%3D%27SV%27)+AND+sig%3D%27A%27"
    "&f=geojson&outFields=prod_type,event,sig,phenom,wfo,expiration,issuance"
    "&returnGeometry=true&outSR=4326"
)
WATCH_CACHE_TTL = 120  # 2 min


@bp.route("/api/overlays/spc-watches", methods=["GET"])
def spc_watches():
    """Proxy SPC Watches from NWS MapServer (GeoJSON)."""
    cached, hit = _cached("spc_watches", WATCH_CACHE_TTL)
    if hit:
        return jsonify(cached)

    try:
        resp = _requests.get(SPC_WATCH_URL, headers=_HEADERS, timeout=15)
        resp.raise_for_status()
        data = resp.json()
        # Normalise property names for the frontend MdWatchLayer
        for f in data.get("features", []):
            p = f.get("properties", {})
            p["type"] = p.get("prod_type", p.get("type", ""))
            p["number"] = p.get("event", p.get("number", ""))
        _store("spc_watches", data)
        return jsonify(data)
    except Exception as e:
        print(f"[overlays] SPC watch fetch error: {e}")
        return jsonify({"type": "FeatureCollection", "features": []}), 502


# ── 5. Spotter Network ───────────────────────────────────────────
SN_POSITIONS_URL = "https://www.spotternetwork.org/feeds/gr.txt"
SN_REPORTS_URL = "https://www.spotternetwork.org/feeds/reports.txt"
SN_CACHE_TTL = 30  # 30s — positions change rapidly


@bp.route("/api/overlays/spotters", methods=["GET"])
def spotter_network():
    """Proxy Spotter Network position + report feeds."""
    cached, hit = _cached("spotters", SN_CACHE_TTL)
    if hit:
        return jsonify(cached)

    positions = []
    reports = []

    try:
        pos_resp = _requests.get(SN_POSITIONS_URL, headers=_HEADERS, timeout=10)
        if pos_resp.ok:
            positions = _parse_spotter_positions(pos_resp.text)
    except Exception as e:
        print(f"[overlays] Spotter positions error: {e}")

    try:
        rep_resp = _requests.get(SN_REPORTS_URL, headers=_HEADERS, timeout=10)
        if rep_resp.ok:
            reports = _parse_spotter_reports(rep_resp.text)
    except Exception as e:
        print(f"[overlays] Spotter reports error: {e}")

    data = {"positions": positions, "reports": reports}
    _store("spotters", data)
    return jsonify(data)


# ── Spotter Network parsers (mirror the frontend logic) ──────────

SN_REPORT_TYPES = {
    1: "Tornado", 2: "Funnel Cloud", 3: "Rotating Wall Cloud",
    4: "Hail", 5: "Wind Damage", 6: "Flooding", 7: "Other",
}


def _parse_spotter_positions(text):
    """Parse the GR-format position feed."""
    lines = text.splitlines()
    spotters = []
    for i, raw in enumerate(lines):
        line = raw.strip()
        if not line.startswith("Object:"):
            continue
        coords = line[7:].strip().split(",")
        try:
            lat, lon = float(coords[0]), float(coords[1])
        except (ValueError, IndexError):
            continue
        name = time_str = heading = ""
        for j in range(i + 1, min(i + 5, len(lines))):
            l = lines[j].strip()
            if l.startswith("Icon:") and '"' in l:
                m = re.search(r'"([^"]*)"', l)
                if m:
                    parts = m.group(1).split("\\n")
                    name = parts[0] if parts else ""
                    time_str = parts[1] if len(parts) > 1 else ""
                    heading = parts[2] if len(parts) > 2 else ""
                break
        spotters.append({
            "lat": lat, "lon": lon,
            "name": name, "time": time_str, "heading": heading,
        })
    return spotters


def _parse_spotter_reports(text):
    """Parse the report feed."""
    reports = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line.startswith("Icon:"):
            continue
        m = re.match(
            r'^Icon:\s*([\d.-]+),([\d.-]+),\d+,(\d+),(\d+),"(.+)"$', line
        )
        if not m:
            continue
        lat, lon = float(m.group(1)), float(m.group(2))
        sheet = int(m.group(3))
        idx = int(m.group(4))
        tooltip = m.group(5)
        parts = tooltip.split("\\n")
        reporter = re.sub(r"^Reported By:\s*", "", parts[0]) if parts else ""
        rtype = parts[1] if len(parts) > 1 else "Other"
        time_part = re.sub(r"^Time:\s*", "", parts[2]) if len(parts) > 2 else ""
        rest = " · ".join(parts[3:]).lstrip("Size: ").lstrip("Notes: ") if len(parts) > 3 else ""
        reports.append({
            "lat": lat, "lon": lon,
            "reporter": reporter, "type": rtype,
            "time": time_part, "notes": rest,
            "sheet": sheet, "idx": idx,
        })
    return reports
