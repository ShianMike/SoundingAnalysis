"""
Meta routes: health check, station list, data-source catalogue, ACARS airports.
"""
from flask import Blueprint, jsonify, request

from sounding import STATIONS, STATION_WMO, DATA_SOURCES, BUFKIT_MODELS, PSU_MODELS

bp = Blueprint("meta", __name__)

# ── Major ACARS/AMDAR airports (high commercial traffic = reliable profiles) ──
ACARS_AIRPORTS = [
    ("KATL", "Atlanta Hartsfield", 33.64, -84.43),
    ("KORD", "Chicago O'Hare", 41.97, -87.91),
    ("KDFW", "Dallas/Fort Worth", 32.90, -97.04),
    ("KDEN", "Denver International", 39.86, -104.67),
    ("KJFK", "New York JFK", 40.64, -73.78),
    ("KLAX", "Los Angeles Intl", 33.94, -118.41),
    ("KSFO", "San Francisco Intl", 37.62, -122.38),
    ("KIAH", "Houston Intercontinental", 29.98, -95.34),
    ("KMSP", "Minneapolis-St Paul", 44.88, -93.22),
    ("KDTW", "Detroit Metro", 42.21, -83.35),
    ("KBOS", "Boston Logan", 42.36, -71.01),
    ("KSEA", "Seattle-Tacoma", 47.45, -122.31),
    ("KPHL", "Philadelphia Intl", 39.87, -75.24),
    ("KPHX", "Phoenix Sky Harbor", 33.44, -112.01),
    ("KMCO", "Orlando Intl", 28.43, -81.31),
    ("KEWR", "Newark Liberty", 40.69, -74.17),
    ("KMIA", "Miami Intl", 25.79, -80.29),
    ("KFLL", "Fort Lauderdale", 26.07, -80.15),
    ("KCLT", "Charlotte Douglas", 35.21, -80.94),
    ("KLAS", "Las Vegas McCarran", 36.08, -115.15),
    ("KBWI", "Baltimore-Washington", 39.18, -76.67),
    ("KSLC", "Salt Lake City", 40.79, -111.98),
    ("KSTL", "St. Louis Lambert", 38.75, -90.37),
    ("KMCI", "Kansas City Intl", 39.30, -94.71),
    ("KMKE", "Milwaukee Mitchell", 42.95, -87.90),
    ("KPIT", "Pittsburgh Intl", 40.49, -80.23),
    ("KIND", "Indianapolis Intl", 39.72, -86.29),
    ("KCMH", "Columbus OH", 39.99, -82.89),
    ("KCVG", "Cincinnati N. KY", 39.05, -84.66),
    ("KBNA", "Nashville Intl", 36.13, -86.68),
    ("KMEM", "Memphis Intl", 35.04, -89.98),
    ("KPBI", "West Palm Beach", 26.68, -80.10),
    ("KRDU", "Raleigh-Durham", 35.88, -78.79),
    ("KSMF", "Sacramento Intl", 38.70, -121.59),
    ("KPDX", "Portland Intl", 45.59, -122.59),
    ("KSAN", "San Diego Lindbergh", 32.73, -117.19),
    ("KTPA", "Tampa Intl", 27.98, -82.53),
    ("KHOU", "Houston Hobby", 29.65, -95.28),
    ("KAUS", "Austin-Bergstrom", 30.19, -97.67),
    ("KSAT", "San Antonio Intl", 29.53, -98.47),
    ("KOAK", "Oakland Intl", 37.72, -122.22),
    ("KSJC", "San Jose Intl", 37.36, -121.93),
    ("KDCA", "Washington Reagan", 38.85, -77.04),
    ("KLGA", "New York LaGuardia", 40.78, -73.87),
    ("KMDW", "Chicago Midway", 41.79, -87.75),
    ("KDAL", "Dallas Love Field", 32.85, -96.85),
    ("KOMA", "Omaha Eppley", 41.30, -95.89),
    ("KOKC", "Oklahoma City", 35.39, -97.60),
    ("KTUL", "Tulsa Intl", 36.20, -95.89),
    ("KLIT", "Little Rock", 34.73, -92.22),
    ("KMSN", "Madison WI", 43.14, -89.34),
    ("KDSM", "Des Moines", 41.53, -93.66),
    ("KSGF", "Springfield MO", 37.25, -93.39),
    ("KICT", "Wichita Eisenhower", 37.65, -97.43),
    ("KABQ", "Albuquerque Sunport", 35.04, -106.61),
    ("KELP", "El Paso Intl", 31.81, -106.38),
    ("KJAN", "Jackson MS", 32.31, -90.08),
    ("KBHM", "Birmingham AL", 33.56, -86.75),
    ("KLEX", "Lexington KY", 38.04, -84.61),
    ("KSDF", "Louisville KY", 38.17, -85.74),
]


@bp.route("/api/health", methods=["GET"])
def health():
    """Lightweight health-check endpoint."""
    return jsonify({"status": "ok"})


@bp.route("/api/stations", methods=["GET"])
def list_stations():
    """Return all known sounding stations with available data sources."""
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


@bp.route("/api/sources", methods=["GET"])
def list_sources():
    """Return available data sources and BUFKIT models."""
    sources = [{"id": k, "name": v} for k, v in DATA_SOURCES.items()]
    models = [{"id": k, "name": v} for k, v in BUFKIT_MODELS.items()]
    psu_models = [{"id": k, "name": v} for k, v in PSU_MODELS.items()]
    return jsonify({"sources": sources, "models": models, "psuModels": psu_models})


@bp.route("/api/acars-airports", methods=["GET"])
def acars_airports():
    """Return major ACARS/AMDAR airport locations for the map overlay."""
    airports = [
        {"id": code, "name": name, "lat": lat, "lon": lon}
        for code, name, lat, lon in ACARS_AIRPORTS
    ]
    return jsonify(airports)
