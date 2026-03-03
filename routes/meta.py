"""
Meta routes: health check, station list, data-source catalogue.
"""
from flask import Blueprint, jsonify

from sounding import STATIONS, STATION_WMO, DATA_SOURCES, BUFKIT_MODELS, PSU_MODELS

bp = Blueprint("meta", __name__)


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
