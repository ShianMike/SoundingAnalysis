"""
Small utility helpers (time resolution, nearest-station lookup).
"""
from datetime import datetime, timedelta, timezone

from .constants import STATIONS


def get_latest_sounding_time():
    """Return the most recent standard sounding time (00Z or 12Z)."""
    now = datetime.now(timezone.utc)
    if now.hour >= 12:
        sounding_hour = 12
    else:
        sounding_hour = 0
    t = now.replace(hour=sounding_hour, minute=0, second=0, microsecond=0)
    # If it's too early for this cycle, go back
    if (now - t).total_seconds() < 5400:  # within 1.5 hours, data may not be up
        t -= timedelta(hours=12)
    return t


def find_nearest_station(lat, lon):
    """Find the nearest sounding station to a lat/lon."""
    best = None
    best_dist = 1e9
    for code, (name, slat, slon) in STATIONS.items():
        dist = ((slat - lat)**2 + (slon - lon)**2)**0.5
        if dist < best_dist:
            best_dist = dist
            best = code
