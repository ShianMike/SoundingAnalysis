"""
VERTICAL PROFILE ANALYSIS TOOL
===============================
Fetches real upper-air sounding data from multiple sources
and produces a comprehensive analysis plot similar to SounderPy.

Data Sources:
  - Observed radiosondes (IEM / University of Wyoming)
  - RAP model analysis for any lat/lon (NCEI THREDDS, requires siphon)
  - BUFKIT forecast soundings: HRRR, RAP, NAM, GFS, etc. (Iowa State)
  - ACARS/AMDAR aircraft observations at airports (IEM)

Features:
  - Skew-T Log-P diagram with temperature, dewpoint, wet-bulb, and parcel traces
  - Hodograph with storm-motion vectors
  - Thermodynamic indices (CAPE, CIN, LCL, LFC, EL, etc.)
  - Kinematic indices (SRH, bulk wind difference, etc.)
  - Wind barbs
  - Storm-relative wind & streamwise vorticity profiles

Usage:
  python sounding.py                                               # Interactive menu
  python sounding.py --station OUN                                 # Latest observed
  python sounding.py --station OUN --date 2024061200               # Specific date
  python sounding.py --source rap --lat 36.4 --lon -99.4           # RAP at any point
  python sounding.py --source bufkit --model hrrr --station OUN    # HRRR forecast
  python sounding.py --source acars --station KDFW                 # Aircraft obs
  python sounding.py --list-sources                                # Show all sources
"""

import argparse
import sys
import warnings
from datetime import datetime, timedelta, timezone
from io import StringIO

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for saving

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as path_effects
import numpy as np
import requests
from scipy.ndimage import gaussian_filter1d

import metpy.calc as mpcalc
from metpy.plots import SkewT, Hodograph
from metpy.units import units, pandas_dataframe_to_unit_arrays
import metpy.constants as mpconst

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────────────
# STATION LIST (common US upper-air sites)
# ─────────────────────────────────────────────────────────────────────
STATIONS = {
    "OUN": ("Norman, OK", 35.22, -97.46),
    "FWD": ("Fort Worth, TX", 32.83, -97.30),
    "AMA": ("Amarillo, TX", 35.23, -101.71),
    "DDC": ("Dodge City, KS", 37.77, -99.97),
    "TOP": ("Topeka, KS", 39.07, -95.62),
    "SGF": ("Springfield, MO", 37.23, -93.40),
    "LZK": ("Little Rock, AR", 34.83, -92.26),
    "SHV": ("Shreveport, LA", 32.45, -93.84),
    "BMX": ("Birmingham, AL", 33.18, -86.77),
    "JAN": ("Jackson, MS", 32.32, -90.08),
    "LCH": ("Lake Charles, LA", 30.12, -93.23),
    "BNA": ("Nashville, TN", 36.25, -86.56),
    "ILN": ("Wilmington, OH", 39.42, -83.82),
    "DTX": ("Detroit, MI", 42.70, -83.47),
    "DVN": ("Davenport, IA", 41.61, -90.58),
    "OAX": ("Omaha, NE", 41.32, -96.37),
    "LBF": ("North Platte, NE", 41.13, -100.68),
    "ABR": ("Aberdeen, SD", 45.45, -98.41),
    "UNR": ("Rapid City, SD", 44.07, -103.21),
    "BIS": ("Bismarck, ND", 46.77, -100.75),
    "GJT": ("Grand Junction, CO", 39.12, -108.53),
    "DNR": ("Denver, CO", 39.77, -104.88),
    "ABQ": ("Albuquerque, NM", 35.04, -106.62),
    "EPZ": ("El Paso, TX", 31.87, -106.70),
    "TFX": ("Great Falls, MT", 47.46, -111.38),
    "SLC": ("Salt Lake City, UT", 40.77, -111.97),
    "BOI": ("Boise, ID", 43.57, -116.22),
    "MFR": ("Medford, OR", 42.37, -122.87),
    "OTX": ("Spokane, WA", 47.68, -117.63),
    "UIL": ("Quillayute, WA", 47.95, -124.55),
    "REV": ("Reno, NV", 39.57, -119.80),
    "VBG": ("Vandenberg, CA", 34.75, -120.57),
    "NKX": ("San Diego, CA", 32.87, -117.15),
    "TUS": ("Tucson, AZ", 32.23, -110.95),
    "FGZ": ("Flagstaff, AZ", 35.23, -111.82),
    "IAD": ("Sterling, VA", 38.98, -77.48),
    "WAL": ("Wallops Island, VA", 37.94, -75.47),
    "MHX": ("Morehead City, NC", 34.78, -76.88),
    "GSO": ("Greensboro, NC", 36.10, -79.95),
    "CHS": ("Charleston, SC", 32.90, -80.03),
    "JAX": ("Jacksonville, FL", 30.50, -81.70),
    "TBW": ("Tampa Bay, FL", 27.70, -82.40),
    "MFL": ("Miami, FL", 25.75, -80.38),
    "TLH": ("Tallahassee, FL", 30.40, -84.35),
    "BUF": ("Buffalo, NY", 42.93, -78.73),
    "ALB": ("Albany, NY", 42.75, -73.80),
    "OKX": ("Upton, NY", 40.87, -72.87),
    "GYX": ("Gray, ME", 43.90, -70.25),
    "CHH": ("Chatham, MA", 41.67, -69.97),
    "CAR": ("Caribou, ME", 46.87, -68.02),
    "PIT": ("Pittsburgh, PA", 40.53, -80.23),
    "RNK": ("Blacksburg, VA", 37.20, -80.40),
    "MPX": ("Minneapolis, MN", 44.85, -93.57),
    "GRB": ("Green Bay, WI", 44.48, -88.13),
    "ILX": ("Lincoln, IL", 40.15, -89.34),
    "APX": ("Gaylord, MI", 44.90, -84.72),
    "INL": ("International Falls, MN", 48.57, -93.38),
    "CRP": ("Corpus Christi, TX", 27.77, -97.50),
    "MAF": ("Midland, TX", 31.95, -102.18),
    "DRT": ("Del Rio, TX", 29.37, -100.92),
    "BRO": ("Brownsville, TX", 25.92, -97.42),
    "KEY": ("Key West, FL", 24.55, -81.75),
    "RIW": ("Riverton, WY", 43.07, -108.48),
    "LMN": ("Lamont, OK", 36.60, -97.48),
}

# University of Wyoming station IDs (3-letter -> WMO number lookup)
STATION_WMO = {
    "OUN": "72357", "FWD": "72249", "AMA": "72363", "DDC": "72451",
    "TOP": "72456", "SGF": "72440", "LZK": "72340", "SHV": "72248",
    "BMX": "72230", "JAN": "72235", "LCH": "72240", "BNA": "72327",
    "ILN": "72426", "DTX": "72632", "DVN": "72558", "OAX": "72558",
    "LBF": "72562", "ABR": "72659", "UNR": "72662", "BIS": "72764",
    "GJT": "72476", "DNR": "72469", "ABQ": "72365", "EPZ": "72270",
    "TFX": "72776", "SLC": "72572", "BOI": "72681", "MFR": "72597",
    "OTX": "72786", "UIL": "72797", "REV": "72489", "VBG": "72393",
    "NKX": "72293", "TUS": "72274", "FGZ": "72376", "IAD": "72403",
    "WAL": "72402", "MHX": "72305", "GSO": "72317", "CHS": "72208",
    "JAX": "72206", "TBW": "72210", "MFL": "72202", "TLH": "72214",
    "BUF": "72528", "ALB": "72518", "OKX": "72501", "GYX": "74389",
    "CHH": "74494", "CAR": "72712", "PIT": "72520", "RNK": "72318",
    "MPX": "72649", "GRB": "72645", "ILX": "74560", "APX": "72634",
    "INL": "72747", "CRP": "72251", "MAF": "72265", "DRT": "72261",
    "BRO": "72250", "KEY": "72201", "RIW": "72672", "DVN": "74455",
    "LMN": "74646", "OAX": "72558",
}


# ─────────────────────────────────────────────────────────────────────
# DATA FETCHING
# ─────────────────────────────────────────────────────────────────────
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
    return best


def fetch_iem_sounding(station_id, dt, quiet=False):
    """
    Fetch sounding data from the Iowa Environmental Mesonet (IEM).
    This is generally more reliable than the UWyo server.
    Returns parsed arrays with units.
    """
    # IEM API works with 3-letter station IDs directly
    url = (
        f"https://mesonet.agron.iastate.edu/json/raob.py?"
        f"ts={dt.strftime('%Y%m%d%H%M')}&station={station_id}"
    )
    
    if not quiet:
        print(f"  Fetching from IEM: {url}")
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    result = resp.json()
    
    profiles = result.get("profiles", [])
    if not profiles:
        raise ValueError(f"No IEM sounding data for {station_id} at {dt}")
    
    profile = profiles[0]  # first match
    levels = profile.get("profile", [])
    
    # Station info - the "station" field may be a string or dict
    station_field = profile.get("station", {})
    if isinstance(station_field, dict):
        station_info = {
            "lat": station_field.get("latitude"),
            "lon": station_field.get("longitude"),
            "elev": station_field.get("elevation"),
            "name": station_field.get("name", ""),
        }
    else:
        # station is just a string like "KOUN"
        stn_entry = STATIONS.get(station_id, (str(station_field), 0, 0))
        station_info = {
            "lat": stn_entry[1],
            "lon": stn_entry[2],
            "name": stn_entry[0],
        }
    
    pressure, height, temp, dewpoint, wind_dir, wind_spd = [], [], [], [], [], []
    
    for lvl in levels:
        p = lvl.get("pres")
        h = lvl.get("hght")
        t = lvl.get("tmpc")
        td = lvl.get("dwpc")
        wd = lvl.get("drct")
        ws = lvl.get("sknt")
        
        # Skip levels missing critical thermo data (p, h, t, td)
        if any(v is None for v in [p, h, t, td]):
            continue
        # Skip clearly bogus values
        if p < 10 or t < -999 or td < -999:
            continue
        
        # Keep level even if wind is missing — interpolate wind later
        pressure.append(float(p))
        height.append(float(h))
        temp.append(float(t))
        dewpoint.append(float(td))
        wind_dir.append(float(wd) if wd is not None else np.nan)
        wind_spd.append(float(ws) if ws is not None else np.nan)
    
    if len(pressure) < 5:
        raise ValueError(
            f"Insufficient data points ({len(pressure)}) from IEM for {station_id}"
        )
    
    # Interpolate missing wind data.
    # Interpolating direction directly can create large vector errors
    # across 0/360-degree wrap, so interpolate u/v components instead.
    wd_arr = np.array(wind_dir)
    ws_arr = np.array(wind_spd)
    valid_wind = ~np.isnan(wd_arr) & ~np.isnan(ws_arr)
    has_wind = valid_wind.copy()  # True for levels with real wind obs
    if np.any(~valid_wind) and np.sum(valid_wind) >= 2:
        h_arr = np.array(height)
        u_valid, v_valid = mpcalc.wind_components(
            ws_arr[valid_wind] * units.knot,
            wd_arr[valid_wind] * units.degree,
        )
        u_all = np.interp(h_arr, h_arr[valid_wind], u_valid.to("knot").magnitude)
        v_all = np.interp(h_arr, h_arr[valid_wind], v_valid.to("knot").magnitude)
        wd_arr = mpcalc.wind_direction(
            u_all * units.knot, v_all * units.knot
        ).to("degree").magnitude
        ws_arr = mpcalc.wind_speed(
            u_all * units.knot, v_all * units.knot
        ).to("knot").magnitude
    
    data = {
        "pressure": np.array(pressure) * units.hPa,
        "height": np.array(height) * units.meter,
        "temperature": np.array(temp) * units.degC,
        "dewpoint": np.array(dewpoint) * units.degC,
        "wind_direction": wd_arr * units.degree,
        "wind_speed": ws_arr * units.knot,
        "has_wind": has_wind,
        "station_info": station_info,
    }
    
    return data


def fetch_wyoming_sounding(station_id, dt):
    """
    Fetch sounding data from the University of Wyoming (fallback).
    Returns parsed arrays with units.
    """
    wmo = STATION_WMO.get(station_id, station_id)

    dt_str = f"{dt.year}-{dt.month:02d}-{dt.day:02d} {dt.hour:02d}:00:00"
    url = "https://weather.uwyo.edu/wsgi/sounding"
    params = {"type": "TEXT:LIST", "datetime": dt_str, "id": wmo}

    print(f"  Fetching from UWyo: {url}  params={params}")
    resp = requests.get(url, params=params, timeout=20)
    resp.raise_for_status()
    html = resp.text

    if "Can't get" in html or "No data" in html:
        raise ValueError(f"No sounding data available for {station_id} at {dt}")

    html_lower = html.lower()
    pre_start = html_lower.find("<pre>")
    pre_end = html_lower.find("</pre>")
    if pre_start == -1 or pre_end == -1:
        raise ValueError("Could not parse sounding data from response")
    
    raw = html[pre_start+5:pre_end].strip()
    lines = raw.split("\n")
    
    header_idx = None
    data_start = None
    for i, line in enumerate(lines):
        if "PRES" in line and "HGHT" in line and "TEMP" in line:
            header_idx = i
        if header_idx is not None and line.strip().startswith("---"):
            data_start = i + 1
            break
    
    if data_start is None:
        raise ValueError("Could not find data in sounding response")
    
    pressure, height, temp, dewpoint, wind_dir, wind_spd = [], [], [], [], [], []
    
    for line in lines[data_start:]:
        line = line.strip()
        if not line or line.lower().startswith("</pre>"):
            break
        parts = line.split()
        if len(parts) < 8:
            continue
        try:
            p = float(parts[0])
            h = float(parts[1])
            t = float(parts[2])
            td = float(parts[3])
            wd = float(parts[6])
            ws = float(parts[7])
            pressure.append(p)
            height.append(h)
            temp.append(t)
            dewpoint.append(td)
            wind_dir.append(wd)
            wind_spd.append(ws)
        except (ValueError, IndexError):
            continue

    if len(pressure) < 5:
        raise ValueError(f"Insufficient data ({len(pressure)} levels) from UWyo")
    
    station_info = {}
    info_start = html.find("Station number:")
    if info_start > -1:
        info_block = html[info_start:info_start+500]
        for info_line in info_block.split("\n"):
            if "Station latitude:" in info_line:
                try: station_info["lat"] = float(info_line.split(":")[-1].strip())
                except: pass
            if "Station longitude:" in info_line:
                try: station_info["lon"] = float(info_line.split(":")[-1].strip())
                except: pass
            if "Station elevation:" in info_line:
                try: station_info["elev"] = float(info_line.split(":")[-1].strip().split()[0])
                except: pass
    
    data = {
        "pressure": np.array(pressure) * units.hPa,
        "height": np.array(height) * units.meter,
        "temperature": np.array(temp) * units.degC,
        "dewpoint": np.array(dewpoint) * units.degC,
        "wind_direction": np.array(wind_dir) * units.degree,
        "wind_speed": np.array(wind_spd) * units('m/s'),
        "has_wind": np.ones(len(pressure), dtype=bool),
        "station_info": station_info,
    }
    
    return data


# ─────────────────────────────────────────────────────────────────────
# ADDITIONAL DATA SOURCES
# ─────────────────────────────────────────────────────────────────────
# Available BUFKIT models (for --model flag)
BUFKIT_MODELS = {
    "rap":     "Rapid Refresh (RAP) - hourly, CONUS 13 km",
    "hrrr":    "High-Res Rapid Refresh (HRRR) - hourly, CONUS 3 km",
    "nam":     "North American Mesoscale (NAM) - hourly, 12 km",
    "namnest": "NAM Nest - hourly, 3 km CONUS",
    "gfs":     "Global Forecast System (GFS) - 3-hourly, global",
    "sref":    "Short-Range Ensemble Forecast (SREF) - 3-hourly",
}


def fetch_rap_sounding(lat, lon, dt):
    """
    Fetch RAP model analysis sounding for any lat/lon via NCEI THREDDS.
    Provides analysis profiles at ~13 km resolution over CONUS.
    Available from ~2012 to near-present (with ~1-2 day archival lag).
    """
    try:
        from siphon.catalog import TDSCatalog
    except ImportError:
        raise ImportError(
            "siphon is required for RAP data. Install with: pip install siphon"
        )

    print(f"  Fetching RAP analysis for ({lat:.2f}, {lon:.2f}) at {dt}...")

    # NCEI archive of RAP 130-km analysis GRIB2 files
    cat_url = (
        f"https://www.ncei.noaa.gov/thredds/catalog/model-rap130anl/"
        f"{dt:%Y%m}/{dt:%Y%m%d}/catalog.xml"
    )

    try:
        cat = TDSCatalog(cat_url)
    except Exception as e:
        raise ValueError(f"Cannot access RAP archive for {dt:%Y-%m-%d}: {e}")

    # Dataset naming convention: rap_130_YYYYMMDD_HH00_000.grb2
    target = f"rap_130_{dt:%Y%m%d}_{dt:%H}00_000.grb2"
    if target not in cat.datasets:
        available = sorted(cat.datasets.keys())
        match = next((n for n in available if f"_{dt:%H}00_000" in n), None)
        if match is None:
            raise ValueError(
                f"RAP file not found for {dt:%Y-%m-%d %H}Z. "
                f"Available files: {available[:5]}…"
            )
        target = match

    ds = cat.datasets[target]
    ncss = ds.subset()

    query = ncss.query()
    query.lonlat_point(lon, lat)
    query.time(dt)
    query.accept("netcdf4")
    query.variables(
        "Temperature_isobaric",
        "Relative_humidity_isobaric",
        "Geopotential_height_isobaric",
        "u-component_of_wind_isobaric",
        "v-component_of_wind_isobaric",
    )

    print(f"  Querying NCSS for point ({lat:.2f}, {lon:.2f})...")
    nc = ncss.get_data(query)

    # ── Helper to find a variable by several possible names ──
    def _var(nc, *names):
        for n in names:
            if n in nc.variables:
                return nc.variables[n][:]
        raise KeyError(
            f"None of {names} found.  Available: {list(nc.variables.keys())}"
        )

    pres_pa = _var(nc, "isobaric", "isobaric1", "isobaric3", "pressure")
    temp_k  = _var(nc, "Temperature_isobaric")
    rh_pct  = _var(nc, "Relative_humidity_isobaric")
    hgt_m   = _var(nc, "Geopotential_height_isobaric")
    u_ms    = _var(nc, "u-component_of_wind_isobaric")
    v_ms    = _var(nc, "v-component_of_wind_isobaric")
    nc.close()

    # Flatten leading (time, …) dimensions
    if temp_k.ndim > 1:
        temp_k = temp_k.squeeze()
        rh_pct = rh_pct.squeeze()
        hgt_m  = hgt_m.squeeze()
        u_ms   = u_ms.squeeze()
        v_ms   = v_ms.squeeze()

    pres_hpa = pres_pa / 100.0
    temp_c   = temp_k - 273.15

    # RH → dewpoint (Magnus-Tetens)
    a, b = 17.625, 243.04
    alpha = (a * temp_c) / (b + temp_c) + np.log(np.clip(rh_pct, 1, 100) / 100.0)
    td_c = (b * alpha) / (a - alpha)

    # Wind components → speed / direction
    wspd_kt = np.sqrt(u_ms**2 + v_ms**2) * 1.94384
    wdir_deg = (np.degrees(np.arctan2(-u_ms, -v_ms)) + 360) % 360

    # Sort surface-first (descending pressure), keep ≥ 50 hPa
    order = np.argsort(pres_hpa)[::-1]
    mask = pres_hpa[order] >= 50
    idx = order[mask]

    data = {
        "pressure":       np.array(pres_hpa[idx])  * units.hPa,
        "height":         np.array(hgt_m[idx])     * units.meter,
        "temperature":    np.array(temp_c[idx])    * units.degC,
        "dewpoint":       np.array(td_c[idx])      * units.degC,
        "wind_direction": np.array(wdir_deg[idx])  * units.degree,
        "wind_speed":     np.array(wspd_kt[idx])   * units.knot,
        "has_wind":       np.ones(len(idx), dtype=bool),
        "station_info": {
            "lat": lat, "lon": lon,
            "name": f"RAP Analysis ({lat:.2f}, {lon:.2f})",
        },
    }
    return data


# ─── BUFKIT FORECAST SOUNDINGS ─────────────────────────────────────
def _parse_bufkit(text, fhour=0):
    """
    Parse a BUFKIT text file and return sounding data for a given
    forecast hour (fhour=0 → analysis).
    """
    lines = text.split("\n")

    # ── Locate SNPARM / STNPRM column definitions ──
    snparm_cols, stnprm_cols = None, None
    for line in lines:
        if "SNPARM" in line and "=" in line:
            raw = line.split("=", 1)[1].strip().rstrip(";")
            snparm_cols = [c.strip() for c in raw.split(";") if c.strip()]
        if "STNPRM" in line and "=" in line:
            raw = line.split("=", 1)[1].strip().rstrip(";")
            stnprm_cols = [c.strip() for c in raw.split(";") if c.strip()]

    if snparm_cols is None:
        raise ValueError("Could not find SNPARM header in BUFKIT data")
    ncols = len(snparm_cols)

    col = {name: i for i, name in enumerate(snparm_cols)}

    # ── Extract station metadata ──
    slat = slon = selv = 0.0
    stn_name = ""
    for line in lines:
        s = line.strip()
        if "STID" in s and "=" in s:
            try:
                stn_name = s.split("STID")[1].split("=")[1].strip().split()[0]
            except (IndexError, ValueError):
                pass
        if "SLAT" in s and "=" in s:
            try:
                slat = float(s.split("SLAT")[1].split("=")[1].strip().split()[0])
            except (IndexError, ValueError):
                pass
        if "SLON" in s and "=" in s:
            try:
                slon = float(s.split("SLON")[1].split("=")[1].strip().split()[0])
            except (IndexError, ValueError):
                pass
        if "SELV" in s and "=" in s:
            try:
                selv = float(s.split("SELV")[1].split("=")[1].strip().split()[0])
            except (IndexError, ValueError):
                pass

    # ── Split file into per-forecast-hour blocks at STIM lines ──
    stim_idx = [i for i, l in enumerate(lines) if l.strip().startswith("STIM")]
    if fhour >= len(stim_idx):
        raise ValueError(
            f"Forecast hour {fhour} not available (max f{len(stim_idx)-1:03d})"
        )

    blk_start = stim_idx[fhour] + 1
    blk_end = stim_idx[fhour + 1] if fhour + 1 < len(stim_idx) else len(lines)

    # Collect all numeric values in the block, then chunk by ncols
    all_vals = []
    for line in lines[blk_start:blk_end]:
        s = line.strip()
        if not s:
            continue
        try:
            vals = [float(v) for v in s.split()]
            all_vals.extend(vals)
        except ValueError:
            continue

    # Sounding levels: first N*ncols values, then remaining are surface params
    n_levels = len(all_vals) // ncols
    if n_levels < 3:
        raise ValueError(f"Only {n_levels} levels found in BUFKIT block")

    sounding_vals = all_vals[: n_levels * ncols]

    pressure, height, temp, dewpoint, wind_dir, wind_spd = [], [], [], [], [], []
    for lev in range(n_levels):
        row = sounding_vals[lev * ncols : (lev + 1) * ncols]
        try:
            p  = row[col.get("PRES", 0)]
            t  = row[col.get("TMPC", 1)]
            td = row[col.get("DWPC", 3)]
            wd = row[col.get("DRCT", 5)]
            ws = row[col.get("SKNT", 6)]
            h  = row[col.get("HGHT", 9)]
        except IndexError:
            continue

        if p < 10 or t < -999 or td < -999:
            continue

        pressure.append(p)
        temp.append(t)
        dewpoint.append(td)
        wind_dir.append(wd)
        wind_spd.append(ws)
        height.append(h)

    if len(pressure) < 5:
        raise ValueError(f"Insufficient BUFKIT data ({len(pressure)} levels)")

    return {
        "pressure":       np.array(pressure)  * units.hPa,
        "height":         np.array(height)    * units.meter,
        "temperature":    np.array(temp)      * units.degC,
        "dewpoint":       np.array(dewpoint)  * units.degC,
        "wind_direction": np.array(wind_dir)  * units.degree,
        "wind_speed":     np.array(wind_spd)  * units.knot,
        "has_wind":       np.ones(len(pressure), dtype=bool),
        "station_info": {
            "lat": slat, "lon": slon, "elev": selv,
            "name": stn_name,
        },
    }


def fetch_bufkit_sounding(station_id, dt, model="rap", fhour=0):
    """
    Fetch BUFKIT model-forecast sounding from Iowa State University archive.
    
    Parameters
    ----------
    station_id : str   3-letter station (e.g. 'OUN')
    dt         : datetime  Model initialisation time (00/06/12/18 Z)
    model      : str   One of: rap, hrrr, nam, namnest, gfs, sref
    fhour      : int   Forecast hour (0 = analysis, 1, 2, …)
    """
    model = model.lower()
    if model not in BUFKIT_MODELS:
        raise ValueError(
            f"Unknown model '{model}'. Options: {list(BUFKIT_MODELS.keys())}"
        )

    stn = station_id.lower()
    stn_icao = f"k{stn}" if len(stn) == 3 else stn

    # Iowa State BUFKIT archive
    url = (
        f"https://mtarchive.geol.iastate.edu/"
        f"{dt:%Y}/{dt:%m}/{dt:%d}/bufkit/{dt:%H}/{model}/{stn_icao}.buf"
    )

    print(f"  Fetching BUFKIT {model.upper()} f{fhour:03d} for {station_id} "
          f"init {dt:%Y-%m-%d %H}Z...")
    print(f"  URL: {url}")

    resp = requests.get(url, timeout=45)
    if resp.status_code == 404:
        # Try without K prefix (some stations use 3-letter IDs)
        url2 = url.replace(f"/{stn_icao}.buf", f"/{stn}.buf")
        resp = requests.get(url2, timeout=45)
        if resp.status_code == 404:
            raise ValueError(
                f"BUFKIT data not found for {station_id} "
                f"({model.upper()}) at {dt:%Y-%m-%d %H}Z. "
                f"Archive may not cover this station/time."
            )
    resp.raise_for_status()

    data = _parse_bufkit(resp.text, fhour=fhour)

    # Enrich station_info with our local metadata when available
    if station_id.upper() in STATIONS:
        name, lat, lon = STATIONS[station_id.upper()]
        data["station_info"]["name"] = (
            f"{name} ({model.upper()} f{fhour:03d})"
        )
        data["station_info"]["lat"] = lat
        data["station_info"]["lon"] = lon

    return data


# ─── ACARS / AMDAR AIRCRAFT OBSERVATIONS ──────────────────────────
def fetch_acars_sounding(airport, dt):
    """
    Fetch ACARS/AMDAR aircraft-observation profiles from the IEM.
    
    Parameters
    ----------
    airport : str  ICAO airport code (e.g. 'KORD', 'KDFW') or 3-letter ID.
    dt      : datetime  Target time (nearest available profile is returned).
    """
    if len(airport) == 3:
        airport = f"K{airport.upper()}"
    airport = airport.upper()

    # IEM ACARS endpoint returns profiles within ±3 h of the given date
    url = (
        f"https://mesonet.agron.iastate.edu/json/acars_profiles.py?"
        f"airport={airport}&date={dt:%Y-%m-%d}"
    )
    print(f"  Fetching ACARS profiles for {airport} on {dt:%Y-%m-%d}...")
    print(f"  URL: {url}")

    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    result = resp.json()

    profiles = result.get("profiles", [])
    if not profiles:
        raise ValueError(
            f"No ACARS profiles found for {airport} on {dt:%Y-%m-%d}. "
            f"ACARS coverage depends on flight traffic at this airport."
        )

    # Pick the profile whose timestamp is closest to dt
    best_profile = None
    best_delta = 1e18
    for prof in profiles:
        ts_str = prof.get("utcvalid", prof.get("ts", ""))
        try:
            ts = datetime.strptime(ts_str[:16], "%Y-%m-%dT%H:%M").replace(
                tzinfo=timezone.utc
            )
        except (ValueError, TypeError):
            continue
        delta = abs((ts - dt).total_seconds())
        if delta < best_delta:
            best_delta = delta
            best_profile = prof

    if best_profile is None:
        raise ValueError(f"No usable ACARS profiles found for {airport}")

    levels = best_profile.get("profile", [])
    if not levels:
        raise ValueError(f"Selected ACARS profile is empty for {airport}")

    pressure, height, temp, dewpoint, wind_dir, wind_spd = (
        [], [], [], [], [], [],
    )

    for lvl in levels:
        p  = lvl.get("pressure")
        h  = lvl.get("altitude")
        t  = lvl.get("temperature")
        td = lvl.get("dewpoint")
        wd = lvl.get("direction")
        ws = lvl.get("speed")

        if p is None or h is None or t is None:
            continue

        pressure.append(float(p))
        height.append(float(h))
        temp.append(float(t))
        dewpoint.append(float(td) if td is not None else float(t) - 30)
        wind_dir.append(float(wd) if wd is not None else np.nan)
        wind_spd.append(float(ws) if ws is not None else np.nan)

    if len(pressure) < 5:
        raise ValueError(
            f"Insufficient ACARS data ({len(pressure)} levels) for {airport}"
        )

    # Sort by pressure descending (surface first)
    order = np.argsort(pressure)[::-1]
    pressure  = np.array(pressure)[order]
    height    = np.array(height)[order]
    temp      = np.array(temp)[order]
    dewpoint  = np.array(dewpoint)[order]
    wind_dir  = np.array(wind_dir)[order]
    wind_spd  = np.array(wind_spd)[order]

    # Interpolate missing winds (same method as IEM RAOB fetcher)
    valid = ~np.isnan(wind_dir) & ~np.isnan(wind_spd)
    has_wind = valid.copy()
    if np.any(~valid) and np.sum(valid) >= 2:
        u_v, v_v = mpcalc.wind_components(
            wind_spd[valid] * units.knot, wind_dir[valid] * units.degree
        )
        u_all = np.interp(height, height[valid], u_v.to("knot").m)
        v_all = np.interp(height, height[valid], v_v.to("knot").m)
        wind_dir = mpcalc.wind_direction(
            u_all * units.knot, v_all * units.knot
        ).to("degree").m
        wind_spd = mpcalc.wind_speed(
            u_all * units.knot, v_all * units.knot
        ).to("knot").m

    data = {
        "pressure":       pressure  * units.hPa,
        "height":         height    * units.meter,
        "temperature":    temp      * units.degC,
        "dewpoint":       dewpoint  * units.degC,
        "wind_direction": wind_dir  * units.degree,
        "wind_speed":     wind_spd  * units.knot,
        "has_wind":       has_wind,
        "station_info": {
            "lat": best_profile.get("latitude", 0),
            "lon": best_profile.get("longitude", 0),
            "name": f"ACARS {airport}",
        },
    }
    return data


# ─────────────────────────────────────────────────────────────────────
# VAD WIND PROFILE (NEXRAD Level-III product 48)
# ─────────────────────────────────────────────────────────────────────
def fetch_vad_data(radar_id):
    """
    Fetch the latest NEXRAD VAD Wind Profile (VWP) for a given radar.

    Uses the NWS TGFTP server to get the latest Level-III product 48 file,
    then parses it with MetPy's Level3File.

    Returns a list of dicts:
      [{"alt_ft": 2000, "alt_m": 610, "dir": 230, "spd_kt": 35, "u_kt": ..., "v_kt": ...}, ...]
    or an empty list on failure.
    """
    from metpy.io import Level3File
    import io as _io

    radar = radar_id.upper()
    # Remove leading 'K' for the TGFTP path if it's a 4-letter ICAO
    radar_lower = radar.lower()
    url = (
        f"https://tgftp.nws.noaa.gov/SL.us008001/DF.of/DC.radar/"
        f"DS.48vwp/SI.{radar_lower}/sn.last"
    )
    try:
        resp = requests.get(url, timeout=12, headers={"User-Agent": "SoundingAnalysis/1.0"})
        resp.raise_for_status()
    except Exception as e:
        print(f"[VAD] Failed to fetch VWP for {radar}: {e}")
        return []

    try:
        f = Level3File(_io.BytesIO(resp.content))
    except Exception as e:
        print(f"[VAD] Failed to parse Level3 VWP for {radar}: {e}")
        return []

    # Extract tabular data from tab_pages[0] — the VAD Algorithm Output table
    if not f.tab_pages or len(f.tab_pages) < 1:
        print(f"[VAD] No tab_pages in VWP for {radar}")
        return []

    full_text = "".join(f.tab_pages[0])
    import re as _re

    winds = []
    # Each data line: altitude(100ft)  U(m/s)  V(m/s)  W(cm/s)  DIR(deg)  SPD(kts) ...
    for match in _re.finditer(
        r"^\s*(\d{3})\s+([\d.\-]+|NA)\s+([\d.\-]+|NA)\s+\S+\s+(\d+|NA)\s+(\d+|NA)",
        full_text,
        _re.MULTILINE,
    ):
        alt_100ft = int(match.group(1))
        alt_ft = alt_100ft * 100
        dir_str = match.group(4)
        spd_str = match.group(5)
        if dir_str == "NA" or spd_str == "NA":
            continue
        wdir = int(dir_str)
        wspd = int(spd_str)
        alt_m = alt_ft * 0.3048
        # Compute u, v in knots from direction/speed
        u_kt = -wspd * np.sin(np.radians(wdir))
        v_kt = -wspd * np.cos(np.radians(wdir))
        winds.append({
            "alt_ft": alt_ft,
            "alt_m": round(alt_m),
            "dir": wdir,
            "spd_kt": wspd,
            "u_kt": round(float(u_kt), 1),
            "v_kt": round(float(v_kt), 1),
        })

    # Also attach metadata
    meta = {}
    if hasattr(f, "metadata"):
        m = f.metadata
        if "vol_time" in m:
            meta["time"] = m["vol_time"].strftime("%Y-%m-%d %H:%MZ")
        if "max" in m:
            meta["max_wind_kt"] = m["max"]
        if "dir_max" in m:
            meta["max_wind_dir"] = m["dir_max"]
        if "alt_max" in m:
            meta["max_wind_alt_ft"] = m["alt_max"]

    return {"winds": winds, "radar": radar, "meta": meta}


# ─────────────────────────────────────────────────────────────────────
# VWP TIME-HEIGHT DISPLAY (multiple VWP snapshots over time)
# ─────────────────────────────────────────────────────────────────────
def _parse_vwp_file(raw_bytes):
    """Parse a single Level-III product 48 file into (datetime, winds_list)."""
    from metpy.io import Level3File
    import io as _io
    import re as _re

    f = Level3File(_io.BytesIO(raw_bytes))
    vol_time = f.metadata.get("vol_time")
    if vol_time is None:
        return None, []

    if not f.tab_pages or len(f.tab_pages) < 1:
        return vol_time, []

    full_text = "".join(f.tab_pages[0])
    winds = []
    for match in _re.finditer(
        r"^\s*(\d{3})\s+([\d.\-]+|NA)\s+([\d.\-]+|NA)\s+\S+\s+(\d+|NA)\s+(\d+|NA)",
        full_text, _re.MULTILINE,
    ):
        alt_100ft = int(match.group(1))
        alt_ft = alt_100ft * 100
        dir_str, spd_str = match.group(4), match.group(5)
        if dir_str == "NA" or spd_str == "NA":
            continue
        wdir, wspd = int(dir_str), int(spd_str)
        winds.append({"alt_ft": alt_ft, "dir": wdir, "spd_kt": wspd})
    return vol_time, winds


def fetch_vwp_timeseries(radar_id, hours=12, max_files=250):
    """
    Fetch multiple NEXRAD VWP snapshots from the NWS TGFTP archive.

    Returns:
      {"radar": str, "snapshots": [{"time": datetime, "winds": [...]}, ...]}
    Snapshots are sorted chronologically (oldest first).
    """
    from datetime import datetime, timezone, timedelta

    radar = radar_id.upper()
    radar_lower = radar.lower()
    base_url = (
        f"https://tgftp.nws.noaa.gov/SL.us008001/DF.of/DC.radar/"
        f"DS.48vwp/SI.{radar_lower}"
    )
    cutoff = datetime.now(timezone.utc) - timedelta(hours=hours)
    snapshots = []

    for i in range(1, max_files + 1):
        fn = f"sn.{i:04d}"
        url = f"{base_url}/{fn}"
        try:
            resp = requests.get(url, timeout=8, headers={"User-Agent": "SoundingAnalysis/1.0"})
            if resp.status_code != 200:
                break
            vol_time, winds = _parse_vwp_file(resp.content)
            if vol_time is None:
                continue
            # Make timezone-aware if needed
            if vol_time.tzinfo is None:
                vol_time = vol_time.replace(tzinfo=timezone.utc)
            if vol_time < cutoff:
                break
            if winds:
                snapshots.append({"time": vol_time, "winds": winds})
        except Exception as e:
            print(f"[VWP-TS] Error fetching {fn} for {radar}: {e}")
            continue

    # Sort chronologically (oldest first)
    snapshots.sort(key=lambda s: s["time"])
    return {"radar": radar, "snapshots": snapshots}


def plot_vwp(vwp_data):
    """
    Generate a VWP (VAD Wind Profile) time-height display.

    Classic radar product: wind barbs stacked by height, spread across time.
    Returns a matplotlib Figure.
    """
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from matplotlib import patheffects as path_effects

    radar = vwp_data["radar"]
    snapshots = vwp_data["snapshots"]

    if not snapshots:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, f"No VWP data available for {radar}",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=16, color="#aaa")
        ax.set_facecolor("#1a1a2e")
        fig.set_facecolor("#1a1a2e")
        return fig

    BG = "#1a1a2e"
    TEXT = "#e0e0e0"
    GRID = "#333355"

    # Collect all unique altitudes
    all_alts = set()
    for snap in snapshots:
        for w in snap["winds"]:
            all_alts.add(w["alt_ft"])
    alt_levels = sorted(all_alts)

    # Filter to reasonable range (surface to 50,000 ft)
    alt_levels = [a for a in alt_levels if a <= 50000]

    fig, ax = plt.subplots(figsize=(14, 8))
    fig.set_facecolor(BG)
    ax.set_facecolor(BG)

    # Speed color scale
    def spd_color(spd_kt):
        if spd_kt >= 100:
            return "#ff2266"
        elif spd_kt >= 70:
            return "#ff6633"
        elif spd_kt >= 50:
            return "#ffaa00"
        elif spd_kt >= 30:
            return "#00cc88"
        elif spd_kt >= 15:
            return "#44aaff"
        else:
            return "#8888aa"

    # Plot wind barbs at each time/height intersection
    times_plotted = []
    for snap in snapshots:
        t = mdates.date2num(snap["time"])
        times_plotted.append(t)
        for w in snap["winds"]:
            alt_ft = w["alt_ft"]
            if alt_ft > 50000:
                continue
            alt_kft = alt_ft / 1000.0
            spd = w["spd_kt"]
            wdir = w["dir"]
            # Convert to u, v for barbs (meteorological convention)
            u = -spd * np.sin(np.radians(wdir))
            v = -spd * np.cos(np.radians(wdir))
            color = spd_color(spd)
            ax.barbs(t, alt_kft, u, v,
                     length=5.5, linewidth=0.6,
                     barbcolor=color, flagcolor=color,
                     sizes=dict(emptybarb=0.04),
                     zorder=5, clip_on=True)

    # Axes formatting
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%MZ"))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_minor_locator(mdates.MinuteLocator(byminute=[0, 15, 30, 45]))

    # Y-axis: altitude in kft
    ax.set_ylabel("Altitude (kft AGL)", color=TEXT, fontsize=11, fontweight="bold")
    ax.set_xlabel("Time (UTC)", color=TEXT, fontsize=11, fontweight="bold")

    # Set y range based on data
    if alt_levels:
        max_alt_kft = min(alt_levels[-1] / 1000.0 + 2, 55)
    else:
        max_alt_kft = 40
    ax.set_ylim(0, max_alt_kft)

    # Grid
    ax.grid(True, axis="both", color=GRID, linewidth=0.5, alpha=0.6)
    ax.tick_params(colors=TEXT, labelsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.spines["left"].set_color(GRID)

    # Add secondary y-axis in km
    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim()[0] * 0.3048, ax.get_ylim()[1] * 0.3048)
    ax2.set_ylabel("Altitude (km AGL)", color=TEXT, fontsize=10, fontweight="bold")
    ax2.tick_params(colors=TEXT, labelsize=9)
    ax2.spines["top"].set_visible(False)
    ax2.spines["left"].set_visible(False)
    ax2.spines["right"].set_color(GRID)
    ax2.spines["bottom"].set_color(GRID)

    # Title
    if snapshots:
        t0 = snapshots[0]["time"].strftime("%Y-%m-%d %H:%MZ")
        t1 = snapshots[-1]["time"].strftime("%H:%MZ")
        title = f"VWP — {radar}  |  {t0} → {t1}  |  {len(snapshots)} scans"
    else:
        title = f"VWP — {radar}"
    ax.set_title(title, color=TEXT, fontsize=13, fontweight="bold", pad=12)

    # Speed legend
    legend_items = [
        (100, "≥100 kt", "#ff2266"),
        (70, "70-99 kt", "#ff6633"),
        (50, "50-69 kt", "#ffaa00"),
        (30, "30-49 kt", "#00cc88"),
        (15, "15-29 kt", "#44aaff"),
        (0, "<15 kt", "#8888aa"),
    ]
    # Build legend with colored text
    legend_y = 0.97
    for _, label, color in legend_items:
        ax.text(1.12, legend_y, label, transform=ax.transAxes,
                fontsize=8, color=color, fontfamily="monospace", fontweight="bold",
                va="top", ha="left",
                path_effects=[path_effects.withStroke(linewidth=2, foreground=BG)])
        legend_y -= 0.045

    fig.tight_layout()
    fig.subplots_adjust(right=0.85)  # Make room for right labels
    return fig


# ─── PSU BUFKIT FEED ────────────────────────────────────────────────
PSU_MODELS = {
    "rap":     "Rapid Refresh (RAP)",
    "nam":     "NAM 12 km",
    "namnest": "NAM Nest 3 km",
    "gfs":     "GFS global",
    "hrrr":    "HRRR 3 km",
    "nam4km":  "NAM 4-km CONUS",
    "hiresw":  "HiResW NMMB / ARW",
    "sref":    "SREF ensemble mean",
}


def fetch_psu_bufkit(station_id, model="rap", fhour=0):
    """
    Fetch latest-run BUFKIT sounding from Penn State's e-wall server.

    PSU provides real-time (latest model cycle) BUFKIT profiles for many
    NWS models.  Unlike the Iowa State archive, only the most recent
    model run is available — no date-based archive.

    URL pattern:
        https://www.meteo.psu.edu/bufkit/data/{MODEL}/{MODEL}_{station}.buf

    Parameters
    ----------
    station_id : str   3-letter station (e.g. 'OUN') or 4-letter ICAO
    model      : str   One of: rap, nam, namnest, gfs, hrrr, nam4km, hiresw, sref
    fhour      : int   Forecast hour offset within the file (0 = analysis)
    """
    model = model.lower()
    if model not in PSU_MODELS:
        raise ValueError(
            f"Unknown PSU model '{model}'. Options: {list(PSU_MODELS.keys())}"
        )

    stn = station_id.upper()
    stn_lower = station_id.lower()
    stn_k = f"k{stn_lower}" if len(stn_lower) == 3 else stn_lower

    # Model name mapping for PSU URL
    model_url = model.upper()

    # PSU uses uppercase model dir and MODEL_station.buf naming
    url = (
        f"https://www.meteo.psu.edu/bufkit/data/"
        f"{model_url}/{model_url}_{stn_k}.buf"
    )

    print(f"  Fetching PSU BUFKIT {model.upper()} f{fhour:03d} for {station_id}...")
    print(f"  URL: {url}")

    resp = requests.get(url, timeout=45)
    if resp.status_code == 404:
        # Try alternate naming — 3-letter without K prefix
        url2 = url.replace(f"_{stn_k}.buf", f"_{stn_lower}.buf")
        resp = requests.get(url2, timeout=45)
        if resp.status_code == 404:
            raise ValueError(
                f"PSU BUFKIT data not found for {station_id} ({model.upper()}). "
                f"The station may not be available in Penn State's feed."
            )
    resp.raise_for_status()

    data = _parse_bufkit(resp.text, fhour=fhour)

    # Enrich station_info
    if stn in STATIONS:
        name, lat, lon = STATIONS[stn]
        data["station_info"]["name"] = (
            f"{name} (PSU {model.upper()} f{fhour:03d})"
        )
        data["station_info"]["lat"] = lat
        data["station_info"]["lon"] = lon
    else:
        sinfo = data["station_info"]
        sinfo["name"] = f"{sinfo.get('name', station_id)} (PSU {model.upper()} f{fhour:03d})"

    return data


# ─────────────────────────────────────────────────────────────────────
# UNIFIED FETCH DISPATCHER
# ─────────────────────────────────────────────────────────────────────
# Recognised source keywords (for --source flag / interactive menu)
DATA_SOURCES = {
    "obs":    "Observed radiosonde (IEM / UWyo)",
    "rap":    "RAP model analysis - any lat/lon, CONUS (NCEI THREDDS)",
    "bufkit": "BUFKIT forecast models - station-based (Iowa State)",
    "psu":    "PSU BUFKIT feed - latest run, station-based (Penn State)",
    "acars":  "ACARS/AMDAR aircraft obs - airport-based (IEM)",
}


def fetch_sounding(station_id, dt, source="obs", lat=None, lon=None,
                   model="rap", fhour=0):
    """
    Unified data-fetching dispatcher.

    Parameters
    ----------
    station_id : str or None
        3-letter station / ICAO airport code.
    dt : datetime
        Target date/time (UTC).
    source : str
        One of: obs, rap, bufkit, acars.
    lat, lon : float or None
        Required for rap (point-based source).
    model : str
        BUFKIT model name (rap, hrrr, nam, …).  Ignored for other sources.
    fhour : int
        BUFKIT forecast hour (0 = analysis).     Ignored for other sources.
    """
    source = source.lower()
    errors = []

    # ── Observed (default) ──────────────────────────────────────────
    if source == "obs":
        try:
            return fetch_iem_sounding(station_id, dt)
        except Exception as e:
            errors.append(f"IEM: {e}")
            print(f"  IEM fetch failed: {e}")
        try:
            return fetch_wyoming_sounding(station_id, dt)
        except Exception as e:
            errors.append(f"UWyo: {e}")
            print(f"  UWyo fetch failed: {e}")
        raise ValueError(
            "Could not fetch observed sounding:\n    "
            + "\n    ".join(errors)
        )

    # ── RAP Analysis ────────────────────────────────────────────────
    if source == "rap":
        if lat is None or lon is None:
            if station_id and station_id.upper() in STATIONS:
                _, lat, lon = STATIONS[station_id.upper()]
            else:
                raise ValueError(
                    "RAP source requires --lat / --lon (or a known --station)."
                )
        return fetch_rap_sounding(lat, lon, dt)

    # ── BUFKIT ──────────────────────────────────────────────────────
    if source == "bufkit":
        if station_id is None:
            raise ValueError("BUFKIT source requires --station.")
        return fetch_bufkit_sounding(station_id, dt, model=model, fhour=fhour)

    # ── PSU BUFKIT (latest run) ─────────────────────────────────────
    if source == "psu":
        if station_id is None:
            raise ValueError("PSU source requires --station.")
        return fetch_psu_bufkit(station_id, model=model, fhour=fhour)

    # ── ACARS ───────────────────────────────────────────────────────
    if source == "acars":
        airport = station_id or ""
        if not airport:
            raise ValueError("ACARS source requires --station (airport code).")
        return fetch_acars_sounding(airport, dt)


    raise ValueError(
        f"Unknown source '{source}'. Options: {list(DATA_SOURCES.keys())}"
    )


# ─────────────────────────────────────────────────────────────────────
# AUTO TORNADO RISK STATION SELECTION
# ─────────────────────────────────────────────────────────────────────
# Key stations spread across tornado-prone regions of the US
TORNADO_SCAN_STATIONS = [
    "OUN", "FWD", "DDC", "AMA", "TOP", "LMN", "OAX", "SGF", "LZK",
    "SHV", "JAN", "BMX", "ILX", "DVN", "MPX", "BNA", "LBF", "MAF",
    "CRP", "DRT", "BRO", "ABQ", "DNR", "GRB", "APX", "ILN", "RNK",
    "IAD", "MHX", "TLH", "TBW", "JAX",
]


def _quick_tornado_score(station_id, dt):
    """
    Fetch sounding for a station and compute quick severe weather composite scores.
    Returns (stp_score, raw_score, cape, srh, bwd, scp, ship, dcp) or None on failure.
    """
    try:
        data = fetch_iem_sounding(station_id, dt, quiet=True)
    except Exception:
        return None

    p = data["pressure"]
    T = data["temperature"]
    Td = data["dewpoint"]
    h = data["height"]
    wdir = data["wind_direction"]
    wspd = data["wind_speed"]

    try:
        u_wind, v_wind = mpcalc.wind_components(wspd, wdir)

        # Surface-based CAPE
        sb_cape, sb_cin = mpcalc.surface_based_cape_cin(p, T, Td)
        cape_val = sb_cape.magnitude if hasattr(sb_cape, "magnitude") else 0
        cin_val = sb_cin.magnitude if hasattr(sb_cin, "magnitude") else 0

        # Most-Unstable CAPE (search lowest 300 hPa)
        try:
            mu_search_top = p[0] - 300 * units.hPa
            mu_mask = p >= mu_search_top
            mu_n = int(np.sum(mu_mask))
            if mu_n < 2:
                mu_n = min(30, len(p))
            mu_idx = np.argmax(
                np.array(mpcalc.equivalent_potential_temperature(
                    p[:mu_n], T[:mu_n], Td[:mu_n]
                ).magnitude)
            )
            mu_prof = mpcalc.parcel_profile(p[mu_idx:], T[mu_idx], Td[mu_idx]).to("degC")
            mu_cape_v, mu_cin_v = mpcalc.cape_cin(p[mu_idx:], T[mu_idx:], Td[mu_idx:], mu_prof)
            mu_cape_val = float(mu_cape_v.magnitude)
            mu_cin_val = float(mu_cin_v.magnitude)
        except Exception:
            mu_cape_val = cape_val
            mu_cin_val = cin_val

        # Height AGL
        h_agl = (h - h[0]).to("meter")

        # Bunkers & SRH
        try:
            rm_u, rm_v, _, _ = mpcalc.bunkers_storm_motion(p, u_wind, v_wind, h)
            _, _, total_srh_1km = mpcalc.storm_relative_helicity(
                h_agl, u_wind, v_wind, 1000 * units.meter,
                storm_u=rm_u, storm_v=rm_v
            )
            _, _, total_srh_3km = mpcalc.storm_relative_helicity(
                h_agl, u_wind, v_wind, 3000 * units.meter,
                storm_u=rm_u, storm_v=rm_v
            )
            srh_val = total_srh_1km.magnitude if hasattr(total_srh_1km, "magnitude") else 0
            srh3_val = total_srh_3km.magnitude if hasattr(total_srh_3km, "magnitude") else 0
        except Exception:
            srh_val = 0
            srh3_val = 0

        # 0-6 km BWD
        try:
            bwd_u, bwd_v = mpcalc.bulk_shear(p, u_wind, v_wind,
                                             height=h_agl,
                                             depth=6000 * units.meter)
            bwd_val = np.sqrt(bwd_u**2 + bwd_v**2).to("knot").magnitude
            bwd_ms = np.sqrt(bwd_u**2 + bwd_v**2).to("m/s").magnitude
        except Exception:
            bwd_val = 0
            bwd_ms = 0

        # STP-like composite
        cape_term = min(cape_val / 1500.0, 3.0)
        srh_term = min(srh_val / 150.0, 3.0)
        bwd_term = min(bwd_val / 20.0, 3.0) if bwd_val >= 12 else 0
        cin_term = min((200.0 + cin_val) / 150.0, 1.0) if cin_val > -250 else 0
        stp_score = cape_term * srh_term * bwd_term * cin_term

        # Additive raw score
        raw_score = (max(cape_val, 0) / 1500.0
                     + max(srh_val, 0) / 150.0
                     + bwd_val / 40.0)

        # SCP = (muCAPE/1000) * (SRH_3km/50) * (BWD_6km_ms/20)
        _scp_bwd = bwd_ms / 20.0 if bwd_ms >= 10.0 else 0.0
        scp_score = (mu_cape_val / 1000.0) * (srh3_val / 50.0) * _scp_bwd

        # SHIP = (muCAPE * mixRatio * LR_7-5 * (-T500) * BWD_6km) / 42M
        try:
            _sfc_mr = float(_mixing_ratio_from_dewpoint(p[0], Td[0]).to("g/kg").magnitude)
            _p700_i = int(np.argmin(np.abs(p.magnitude - 700.0)))
            _p500_i = int(np.argmin(np.abs(p.magnitude - 500.0)))
            _t700 = float(T.magnitude[_p700_i])
            _t500 = float(T.magnitude[_p500_i])
            _lr75 = (_t700 - _t500) / ((h.magnitude[_p500_i] - h.magnitude[_p700_i]) / 1000.0)
            _neg500 = max(-_t500, 0)
            ship_raw = (mu_cape_val * _sfc_mr * _lr75 * _neg500 * bwd_ms) / 42_000_000.0
            if mu_cape_val < 1300 or mu_cin_val < -200:
                ship_raw = 0.0
            ship_score = max(ship_raw, 0)
        except Exception:
            ship_score = 0

        # DCP = (DCAPE/980) * (muCAPE/2000) * (BWD/20) * (meanWind_0-6/16)
        try:
            dcape_v, _dp, _dt = mpcalc.downdraft_cape(p, T, Td)
            _dcape_f = float(dcape_v.magnitude)
            # Mean wind 0-6 km
            _mask06 = h_agl.magnitude <= 6000.0
            _u06 = u_wind[_mask06].to("m/s").magnitude
            _v06 = v_wind[_mask06].to("m/s").magnitude
            _mw06 = np.sqrt(np.mean(_u06)**2 + np.mean(_v06)**2)
            dcp_score = (_dcape_f / 980.0) * (mu_cape_val / 2000.0) * (bwd_ms / 20.0) * (_mw06 / 16.0)
        except Exception:
            dcp_score = 0

        return (float(stp_score), float(raw_score), cape_val, srh_val, bwd_val,
                float(scp_score), float(ship_score), float(dcp_score))

    except Exception:
        return None


def find_highest_tornado_risk(dt, stations=None):
    """
    Scan stations, display a ranked table with thermodynamic values,
    and let the user pick one by number or station ID.
    Pressing Enter (blank) auto-selects the highest-risk station.
    """
    if stations is None:
        stations = TORNADO_SCAN_STATIONS

    print(f"\n  Scanning {len(stations)} stations for tornado risk parameters...")
    print(f"  (this fetches a quick sounding from each site - may take a moment)\n")

    # Collect results: (station_id, stp, raw, cape, srh, bwd, name, scp, ship, dcp)
    results = []
    for sid in stations:
        result = _quick_tornado_score(sid, dt)
        if result is None:
            continue
        stp_score, raw_score, cape, srh, bwd, scp, ship, dcp = result
        name = STATIONS.get(sid, (sid,))[0]
        results.append((sid, stp_score, raw_score, cape, srh, bwd, name, scp, ship, dcp))

    if not results:
        raise ValueError("Could not fetch data from any station")

    # Sort by (STP desc, Raw desc) so best is at the top
    results.sort(key=lambda r: (r[1], r[2]), reverse=True)

    # ── Print ranked table ──────────────────────────────────────────
    print(f"  {'#':>3s}  {'ID':5s}  {'Station':22s}  {'STP':>7s}  {'SCP':>7s}  {'SHIP':>7s}  {'DCP':>7s}  {'CAPE':>8s}  "
          f"{'0-1 SRH':>8s}  {'0-6 BWD':>7s}")
    print(f"  {'-'*92}")

    for idx, (sid, stp, raw, cape, srh, bwd, name, scp, ship, dcp) in enumerate(results, 1):
        highlight = " <--" if idx == 1 else ""
        print(f"  {idx:3d}  {sid:5s}  {name:22s}  {stp:7.2f}  {scp:7.2f}  {ship:7.2f}  {dcp:7.2f}  "
              f"{cape:8.0f}  {srh:8.0f}  {bwd:7.0f}{highlight}")

    best_sid = results[0][0]
    best_name = results[0][6]

    # ── Interactive selection ────────────────────────────────────────
    print(f"\n  >> Highest risk: #{1} {best_sid} ({best_name})  "
          f"STP={results[0][1]:.2f}  CAPE={results[0][3]:.0f}  "
          f"SRH={results[0][4]:.0f}  BWD={results[0][5]:.0f}kt")
    print(f"    Enter a number (1-{len(results)}), station ID, "
          f"or press Enter to auto-select #{1}:")

    try:
        choice = input("\n  Select station: ").strip()
    except (EOFError, KeyboardInterrupt):
        choice = ""

    if choice == "":
        # Auto-select highest
        print(f"  => Auto-selected: {best_sid} ({best_name})")
        return best_sid

    # Try as a number
    try:
        num = int(choice)
        if 1 <= num <= len(results):
            picked = results[num - 1]
            print(f"  => Selected: {picked[0]} ({picked[6]})")
            return picked[0]
        else:
            print(f"  Invalid number. Auto-selecting {best_sid}.")
            return best_sid
    except ValueError:
        pass

    # Try as a station ID (case-insensitive)
    choice_upper = choice.upper()
    for sid, stp, raw, cape, srh, bwd, name in results:
        if sid == choice_upper:
            print(f"  => Selected: {sid} ({name})")
            return sid

    # Also allow matching by partial name
    for sid, stp, raw, cape, srh, bwd, name in results:
        if choice_upper in name.upper():
            print(f"  => Matched: {sid} ({name})")
            return sid

    print(f"  Station '{choice}' not found in results. Auto-selecting {best_sid}.")
    return best_sid


# ─────────────────────────────────────────────────────────────────────
# PROFILE MERGING
# ─────────────────────────────────────────────────────────────────────

def merge_profiles(data_a, data_b, weight_a=0.5):
    """Merge two sounding profiles into a single blended profile.

    Both profiles are interpolated onto a common pressure grid, then T, Td,
    u, and v are combined via weighted average.

    Parameters
    ----------
    data_a, data_b : dict
        Sounding dicts as returned by ``fetch_sounding``.
    weight_a : float
        Weight for profile A (0–1).  Profile B gets ``1 - weight_a``.

    Returns
    -------
    dict
        A merged sounding dict with the same keys, ready for
        ``compute_parameters`` / ``plot_sounding``.
    """
    weight_b = 1.0 - weight_a

    # Build common pressure grid: surface = max of the two, top = min
    p_a = data_a["pressure"].to("hPa").magnitude
    p_b = data_b["pressure"].to("hPa").magnitude

    p_top = max(p_a.min(), p_b.min())
    p_bot = min(p_a.max(), p_b.max())

    # 5-hPa spacing, descending (surface → top)
    common_p = np.arange(p_bot, p_top - 1, -5.0)
    common_p = common_p[common_p >= p_top]

    def _interp(p_orig, vals):
        """Log-pressure linear interpolation (p decreasing)."""
        lp_orig = np.log(p_orig)
        lp_target = np.log(common_p)
        return np.interp(lp_target, lp_orig[::-1], vals[::-1])

    # Interpolate thermodynamic fields
    T_a  = _interp(p_a, data_a["temperature"].to("degC").magnitude)
    T_b  = _interp(p_b, data_b["temperature"].to("degC").magnitude)
    Td_a = _interp(p_a, data_a["dewpoint"].to("degC").magnitude)
    Td_b = _interp(p_b, data_b["dewpoint"].to("degC").magnitude)

    u_a, v_a = mpcalc.wind_components(data_a["wind_speed"], data_a["wind_direction"])
    u_b, v_b = mpcalc.wind_components(data_b["wind_speed"], data_b["wind_direction"])

    u_a_i = _interp(p_a, u_a.to("knot").magnitude)
    u_b_i = _interp(p_b, u_b.to("knot").magnitude)
    v_a_i = _interp(p_a, v_a.to("knot").magnitude)
    v_b_i = _interp(p_b, v_b.to("knot").magnitude)

    # Height (use hypsometric from merged temperature)
    h_a = _interp(p_a, data_a["height"].to("meter").magnitude)
    h_b = _interp(p_b, data_b["height"].to("meter").magnitude)

    # Weighted blend
    T_m  = weight_a * T_a  + weight_b * T_b
    Td_m = weight_a * Td_a + weight_b * Td_b
    u_m  = weight_a * u_a_i + weight_b * u_b_i
    v_m  = weight_a * v_a_i + weight_b * v_b_i
    h_m  = weight_a * h_a   + weight_b * h_b

    # Ensure Td <= T
    Td_m = np.minimum(Td_m, T_m)

    # Convert wind components back to direction / speed
    wspd_m = np.sqrt(u_m**2 + v_m**2) * units.knot
    wdir_m = (np.degrees(np.arctan2(-u_m, -v_m)) % 360) * units.degree

    # Build station_info from profile A
    info_a = data_a.get("station_info", {})
    info_b = data_b.get("station_info", {})

    merged = {
        "pressure":       common_p * units.hPa,
        "temperature":    T_m * units("degC"),
        "dewpoint":       Td_m * units("degC"),
        "height":         h_m * units.meter,
        "wind_direction": wdir_m,
        "wind_speed":     wspd_m,
        "station_info": {
            "lat":  info_a.get("lat", info_b.get("lat")),
            "lon":  info_a.get("lon", info_b.get("lon")),
            "elev": (weight_a * info_a.get("elev", h_a[0])
                     + weight_b * info_b.get("elev", h_b[0])),
        },
    }
    return merged


# ─────────────────────────────────────────────────────────────────────
# CALCULATIONS
# ─────────────────────────────────────────────────────────────────────

def _mixing_ratio_from_dewpoint(pressure, dewpoint):
    """Compute mixing ratio from dewpoint (MetPy 1.7+ compatible).

    Replaces deprecated mpcalc.mixing_ratio_from_dewpoint by using:
        e = saturation_vapor_pressure(Td)
        r = mixing_ratio(e, p)
    """
    e = mpcalc.saturation_vapor_pressure(dewpoint)
    return mpcalc.mixing_ratio(e, pressure)


def compute_parameters(data, storm_motion=None, surface_mod=None, smoothing=None):
    """Compute a comprehensive set of thermodynamic & kinematic parameters.

    Parameters
    ----------
    data : dict
        Sounding data with pressure, temperature, dewpoint, height, wind_direction, wind_speed.
    storm_motion : dict or None
        Custom storm motion override: {"direction": degrees, "speed": knots}.
        When provided, overrides the Bunkers right-mover storm motion.
    surface_mod : dict or None
        Surface modification: {"temperature": °C, "dewpoint": °C,
        "wind_speed": knots, "wind_direction": degrees}.
        When provided, replaces surface-level values before computation.
    smoothing : float or None
        Gaussian smoothing sigma (in number of data levels).
        Typical values: 2-5. Applied to T, Td, u, v before computation.
        Useful for noisy profiles (e.g. ACARS). Preserves surface values.
    """
    p = data["pressure"]
    T = data["temperature"].copy()
    Td = data["dewpoint"].copy()
    h = data["height"]
    wdir = data["wind_direction"].copy()
    wspd = data["wind_speed"].copy()

    # ── Surface modification ─────────────────────────────────────────
    if surface_mod:
        if surface_mod.get("temperature") is not None:
            T[0] = units.Quantity(float(surface_mod["temperature"]), "degC")
        if surface_mod.get("dewpoint") is not None:
            Td[0] = units.Quantity(float(surface_mod["dewpoint"]), "degC")
        if surface_mod.get("wind_speed") is not None:
            wspd[0] = units.Quantity(float(surface_mod["wind_speed"]), "knot")
        if surface_mod.get("wind_direction") is not None:
            wdir[0] = units.Quantity(float(surface_mod["wind_direction"]), "degree")
        params_mod_applied = True
    else:
        params_mod_applied = False

    params = {}
    params["surface_modified"] = params_mod_applied
    params["smoothing_applied"] = smoothing is not None and smoothing > 0

    # ── Profile smoothing (Gaussian) ─────────────────────────────────
    if smoothing and smoothing > 0 and len(T) > 5:
        sigma = float(smoothing)
        # Smooth T and Td (preserving surface values)
        T_raw = T.magnitude.copy()
        Td_raw = Td.magnitude.copy()
        T_smooth = gaussian_filter1d(T_raw, sigma=sigma, mode='nearest')
        Td_smooth = gaussian_filter1d(Td_raw, sigma=sigma, mode='nearest')
        # Keep surface value unchanged
        T_smooth[0] = T_raw[0]
        Td_smooth[0] = Td_raw[0]
        # Ensure Td <= T at every level
        Td_smooth = np.minimum(Td_smooth, T_smooth)
        T = T_smooth * T.units
        Td = Td_smooth * Td.units

        # Smooth wind components (convert to u/v, smooth, convert back)
        u_raw, v_raw = mpcalc.wind_components(wspd, wdir)
        u_smooth = gaussian_filter1d(u_raw.to("knot").magnitude, sigma=sigma, mode='nearest')
        v_smooth = gaussian_filter1d(v_raw.to("knot").magnitude, sigma=sigma, mode='nearest')
        # Preserve surface wind
        u_smooth[0] = u_raw.to("knot").magnitude[0]
        v_smooth[0] = v_raw.to("knot").magnitude[0]
        wspd = mpcalc.wind_speed(u_smooth * units.knot, v_smooth * units.knot)
        wdir = mpcalc.wind_direction(u_smooth * units.knot, v_smooth * units.knot)
        print(f"  Smoothing applied: sigma={sigma}, {len(T)} levels")

    # Wind components
    u, v = mpcalc.wind_components(wspd, wdir)
    params["u"] = u
    params["v"] = v
    
    # Wet-bulb temperature
    try:
        params["wetbulb"] = mpcalc.wet_bulb_temperature(p, T, Td)
    except:
        params["wetbulb"] = None
    
    # Virtual temperature
    try:
        mr = mpcalc.mixing_ratio_from_relative_humidity(p, T, mpcalc.relative_humidity_from_dewpoint(T, Td))
        params["virtual_temp"] = mpcalc.virtual_temperature(T, mr)
    except:
        params["virtual_temp"] = None
    
    # Surface-based parcel
    try:
        sb_prof = mpcalc.parcel_profile(p, T[0], Td[0]).to("degC")
        params["sb_profile"] = sb_prof
        sb_cape, sb_cin = mpcalc.cape_cin(p, T, Td, sb_prof)
        params["sb_cape"] = sb_cape
        params["sb_cin"] = sb_cin
        params["sb_lcl_p"], params["sb_lcl_t"] = mpcalc.lcl(p[0], T[0], Td[0])
        # Compute LCL height in meters AGL using hypsometric approximation
        try:
            lcl_idx = np.argmin(np.abs(p.magnitude - params["sb_lcl_p"].magnitude))
            sb_lcl_h_msl = h.magnitude[lcl_idx]
            if lcl_idx > 0:
                frac = (p.magnitude[lcl_idx-1] - params["sb_lcl_p"].magnitude) / \
                       (p.magnitude[lcl_idx-1] - p.magnitude[lcl_idx])
                sb_lcl_h_msl = h.magnitude[lcl_idx-1] + frac * (h.magnitude[lcl_idx] - h.magnitude[lcl_idx-1])
            params["sb_lcl_m"] = sb_lcl_h_msl - h[0].magnitude
        except:
            params["sb_lcl_m"] = None
        try:
            params["sb_lfc_p"], params["sb_lfc_t"] = mpcalc.lfc(p, T, Td)
        except:
            params["sb_lfc_p"] = None
        try:
            params["sb_el_p"], params["sb_el_t"] = mpcalc.el(p, T, Td)
        except:
            params["sb_el_p"] = None
    except Exception as e:
        print(f"  Warning: SB parcel calc failed: {e}")
        params["sb_cape"] = 0 * units("J/kg")
        params["sb_cin"] = 0 * units("J/kg")
    
    # Most-Unstable parcel (search within bottom 300 hPa)
    try:
        # Limit search to the lowest 300 hPa above surface
        mu_search_top = p[0] - 300 * units.hPa
        mu_mask = p >= mu_search_top
        mu_search_n = int(np.sum(mu_mask))
        if mu_search_n < 2:
            mu_search_n = min(30, len(p))
        mu_idx = np.argmax(
            np.array(mpcalc.equivalent_potential_temperature(
                p[:mu_search_n], T[:mu_search_n], Td[:mu_search_n]
            ).magnitude)
        )
        mu_p, mu_t, mu_td = p[mu_idx], T[mu_idx], Td[mu_idx]
        mu_prof = mpcalc.parcel_profile(p[mu_idx:], mu_t, mu_td).to("degC")
        params["mu_profile"] = mu_prof
        params["mu_start_idx"] = mu_idx
        mu_cape, mu_cin = mpcalc.cape_cin(p[mu_idx:], T[mu_idx:], Td[mu_idx:], mu_prof)
        params["mu_cape"] = mu_cape
        params["mu_cin"] = mu_cin
        params["mu_lcl_p"], params["mu_lcl_t"] = mpcalc.lcl(mu_p, mu_t, mu_td)
        try:
            mu_lcl_idx = np.argmin(np.abs(p.magnitude - params["mu_lcl_p"].magnitude))
            mu_lcl_h_msl = h.magnitude[mu_lcl_idx]
            if mu_lcl_idx > 0:
                frac = (p.magnitude[mu_lcl_idx-1] - params["mu_lcl_p"].magnitude) / \
                       (p.magnitude[mu_lcl_idx-1] - p.magnitude[mu_lcl_idx])
                mu_lcl_h_msl = h.magnitude[mu_lcl_idx-1] + frac * (h.magnitude[mu_lcl_idx] - h.magnitude[mu_lcl_idx-1])
            params["mu_lcl_m"] = mu_lcl_h_msl - h[0].magnitude
        except:
            params["mu_lcl_m"] = None
        # MU LFC and EL (needed for effective inflow layer & NCAPE)
        try:
            params["mu_lfc_p"], params["mu_lfc_t"] = mpcalc.lfc(p[mu_idx:], T[mu_idx:], Td[mu_idx:], mu_prof)
        except:
            params["mu_lfc_p"] = None
        try:
            params["mu_el_p"], params["mu_el_t"] = mpcalc.el(p[mu_idx:], T[mu_idx:], Td[mu_idx:], mu_prof)
        except:
            params["mu_el_p"] = None
    except Exception as e:
        print(f"  Warning: MU parcel calc failed: {e}")
        params["mu_cape"] = 0 * units("J/kg")
        params["mu_cin"] = 0 * units("J/kg")
    
    # Mixed-layer parcel (100 hPa deep)
    try:
        ml_p, ml_t, ml_td = mpcalc.mixed_parcel(p, T, Td, depth=100 * units.hPa)
        ml_prof = mpcalc.parcel_profile(p, ml_t, ml_td).to("degC")
        params["ml_profile"] = ml_prof
        ml_cape, ml_cin = mpcalc.cape_cin(p, T, Td, ml_prof)
        params["ml_cape"] = ml_cape
        params["ml_cin"] = ml_cin
        params["ml_lcl_p"], params["ml_lcl_t"] = mpcalc.lcl(p[0], ml_t, ml_td)
        try:
            ml_lcl_idx = np.argmin(np.abs(p.magnitude - params["ml_lcl_p"].magnitude))
            ml_lcl_h_msl = h.magnitude[ml_lcl_idx]
            if ml_lcl_idx > 0:
                frac = (p.magnitude[ml_lcl_idx-1] - params["ml_lcl_p"].magnitude) / \
                       (p.magnitude[ml_lcl_idx-1] - p.magnitude[ml_lcl_idx])
                ml_lcl_h_msl = h.magnitude[ml_lcl_idx-1] + frac * (h.magnitude[ml_lcl_idx] - h.magnitude[ml_lcl_idx-1])
            params["ml_lcl_m"] = ml_lcl_h_msl - h[0].magnitude
        except:
            params["ml_lcl_m"] = None
        # ML LFC and EL
        try:
            params["ml_lfc_p"], params["ml_lfc_t"] = mpcalc.lfc(p, T, Td, ml_prof)
            _ml_lfc_idx = np.argmin(np.abs(p.magnitude - params["ml_lfc_p"].magnitude))
            _ml_lfc_h_msl = h.magnitude[_ml_lfc_idx]
            if _ml_lfc_idx > 0:
                _frac = (p.magnitude[_ml_lfc_idx-1] - params["ml_lfc_p"].magnitude) / \
                        (p.magnitude[_ml_lfc_idx-1] - p.magnitude[_ml_lfc_idx])
                _ml_lfc_h_msl = h.magnitude[_ml_lfc_idx-1] + _frac * (h.magnitude[_ml_lfc_idx] - h.magnitude[_ml_lfc_idx-1])
            params["ml_lfc_m"] = _ml_lfc_h_msl - h[0].magnitude
        except:
            params["ml_lfc_p"] = None
            params["ml_lfc_m"] = None
        try:
            params["ml_el_p"], params["ml_el_t"] = mpcalc.el(p, T, Td, ml_prof)
            _ml_el_idx = np.argmin(np.abs(p.magnitude - params["ml_el_p"].magnitude))
            _ml_el_h_msl = h.magnitude[_ml_el_idx]
            if _ml_el_idx > 0:
                _frac = (p.magnitude[_ml_el_idx-1] - params["ml_el_p"].magnitude) / \
                        (p.magnitude[_ml_el_idx-1] - p.magnitude[_ml_el_idx])
                _ml_el_h_msl = h.magnitude[_ml_el_idx-1] + _frac * (h.magnitude[_ml_el_idx] - h.magnitude[_ml_el_idx-1])
            params["ml_el_m"] = _ml_el_h_msl - h[0].magnitude
        except:
            params["ml_el_p"] = None
            params["ml_el_m"] = None
    except Exception as e:
        print(f"  Warning: ML parcel calc failed: {e}")
        params["ml_cape"] = 0 * units("J/kg")
        params["ml_cin"] = 0 * units("J/kg")
    
    # DCAPE (Downdraft CAPE) and downdraft parcel profile
    # Computed early so DCP and DCIN can reference it
    # MetPy returns (dcape, down_pressure, down_parcel_trace)
    try:
        dcape_val, _down_p, dtemp = mpcalc.downdraft_cape(p, T, Td)
        params["dcape"] = dcape_val
        params["dcape_profile"] = dtemp.to("degC")
        params["dcape_pressure"] = _down_p  # pressure levels for the downdraft parcel
    except Exception as e:
        print(f"  Warning: DCAPE calc failed: {e}")
        params["dcape"] = None
        params["dcape_profile"] = None
        params["dcape_pressure"] = None

    # Significant tornado parameter (STP) — computed after Bunkers & SRH below
    
    # Height AGL for SRH and BWD calculations
    h_agl_calc = (h - h[0]).to("meter")
    
    # ── Interpolate to uniform 100 m vertical spacing for kinematic calcs ──
    # This matches SounderPy's approach and produces more stable SRH/BWD
    # integrals than using raw (unevenly-spaced) RAOB levels.
    max_agl = min(float(h_agl_calc[-1].magnitude), 12000.0)
    h_interp = np.arange(0, max_agl + 100, 100) * units.meter  # 0, 100, 200, ... m AGL
    u_interp = np.interp(h_interp.magnitude, h_agl_calc.magnitude, u.to("knot").magnitude) * units.knot
    v_interp = np.interp(h_interp.magnitude, h_agl_calc.magnitude, v.to("knot").magnitude) * units.knot
    p_interp = np.interp(h_interp.magnitude, h_agl_calc.magnitude, p.magnitude) * p.units
    # Also store MSL heights for Bunkers
    h_msl_interp = h_interp + h[0]
    
    params["h_interp"] = h_interp        # 100 m AGL grid
    params["u_interp"] = u_interp
    params["v_interp"] = v_interp
    
    # Bunkers storm motion (uses interpolated profiles for consistency)
    try:
        rm, lm, mw = mpcalc.bunkers_storm_motion(p_interp, u_interp, v_interp, h_msl_interp)
        params["rm_u"], params["rm_v"] = rm
        params["lm_u"], params["lm_v"] = lm
        params["mw_u"], params["mw_v"] = mw
    except Exception as e:
        print(f"  Warning: Bunkers calc failed: {e}")
        params["rm_u"] = params["rm_v"] = 0 * units("m/s")
        params["lm_u"] = params["lm_v"] = 0 * units("m/s")
        params["mw_u"] = params["mw_v"] = 0 * units("m/s")

    # ── Custom storm motion override ─────────────────────────────────
    if storm_motion and storm_motion.get("direction") is not None and storm_motion.get("speed") is not None:
        try:
            sm_dir = float(storm_motion["direction"]) * units.degree
            sm_spd = float(storm_motion["speed"]) * units.knot
            sm_u, sm_v = mpcalc.wind_components(sm_spd, sm_dir)
            params["rm_u"] = sm_u.to("m/s")
            params["rm_v"] = sm_v.to("m/s")
            params["custom_storm_motion"] = True
            print(f"  Using custom storm motion: {sm_dir.magnitude}° @ {sm_spd.magnitude} kt")
        except Exception as e:
            print(f"  Warning: Custom storm motion failed: {e}")
            params["custom_storm_motion"] = False
    else:
        params["custom_storm_motion"] = False
    
    # Storm-relative helicity (0-500m, 0-1km, 0-3km) on interpolated grid
    try:
        _, _, params["srh_500m"] = mpcalc.storm_relative_helicity(
            h_interp, u_interp, v_interp, depth=500 * units.meter,
            storm_u=params["rm_u"], storm_v=params["rm_v"]
        )
        _, _, params["srh_1km"] = mpcalc.storm_relative_helicity(
            h_interp, u_interp, v_interp, depth=1000 * units.meter,
            storm_u=params["rm_u"], storm_v=params["rm_v"]
        )
        _, _, params["srh_3km"] = mpcalc.storm_relative_helicity(
            h_interp, u_interp, v_interp, depth=3000 * units.meter,
            storm_u=params["rm_u"], storm_v=params["rm_v"]
        )
    except Exception as e:
        print(f"  Warning: SRH calc failed: {e}")
        params["srh_500m"] = 0 * units("m^2/s^2")
        params["srh_1km"] = 0 * units("m^2/s^2")
        params["srh_3km"] = 0 * units("m^2/s^2")
    
    # Streamwiseness profile (fraction of horizontal vorticity that is streamwise)
    try:
        rm_u_ms = params["rm_u"].to("m/s").magnitude
        rm_v_ms = params["rm_v"].to("m/s").magnitude
        u_ms = u_interp.to("m/s").magnitude
        v_ms = v_interp.to("m/s").magnitude
        dz = 100.0  # 100 m spacing
        
        # Storm-relative wind components
        u_sr = u_ms - rm_u_ms
        v_sr = v_ms - rm_v_ms
        
        # Horizontal vorticity from vertical wind shear: ωh = (dv/dz, -du/dz)
        dudz = np.gradient(u_ms, dz)
        dvdz = np.gradient(v_ms, dz)
        omega_x = dvdz     # crosswise component of horiz vorticity
        omega_y = -dudz    # other component
        omega_h_mag = np.sqrt(omega_x**2 + omega_y**2)
        
        # Storm-relative wind unit vector
        sr_spd = np.sqrt(u_sr**2 + v_sr**2)
        sr_spd_safe = np.where(sr_spd > 0.1, sr_spd, 0.1)  # avoid div by zero
        sr_hat_u = u_sr / sr_spd_safe
        sr_hat_v = v_sr / sr_spd_safe
        
        # Streamwise vorticity = dot(ωh, SR_hat)
        omega_s = omega_x * sr_hat_u + omega_y * sr_hat_v
        
        # Streamwiseness = |ω_s| / |ω_h|  (0 = fully crosswise, 1 = fully streamwise)
        omega_h_safe = np.where(omega_h_mag > 1e-6, omega_h_mag, 1e-6)
        streamwiseness = np.abs(omega_s) / omega_h_safe
        streamwiseness = np.clip(streamwiseness, 0, 1)
        
        # Also store the sign: positive = cyclonic (streamwise), negative = anticyclonic
        sw_signed = np.sign(omega_s) * streamwiseness
        
        params["streamwiseness"] = streamwiseness
        params["streamwiseness_signed"] = sw_signed
        params["streamwiseness_height"] = h_interp.magnitude / 1000.0  # km AGL
    except Exception as e:
        print(f"  Warning: Streamwiseness calc failed: {e}")
        params["streamwiseness"] = None
    
    # Bulk wind difference (shear) on interpolated grid
    try:
        bwd_u_05, bwd_v_05 = mpcalc.bulk_shear(p_interp, u_interp, v_interp,
                                               height=h_interp, depth=500 * units.meter)
        bwd_u_1, bwd_v_1 = mpcalc.bulk_shear(p_interp, u_interp, v_interp,
                                               height=h_interp, depth=1000 * units.meter)
        bwd_u_3, bwd_v_3 = mpcalc.bulk_shear(p_interp, u_interp, v_interp,
                                               height=h_interp, depth=3000 * units.meter)
        bwd_u_6, bwd_v_6 = mpcalc.bulk_shear(p_interp, u_interp, v_interp,
                                               height=h_interp, depth=6000 * units.meter)
        params["bwd_500m"] = np.sqrt(bwd_u_05**2 + bwd_v_05**2).to("knot")
        params["bwd_1km"] = np.sqrt(bwd_u_1**2 + bwd_v_1**2).to("knot")
        params["bwd_3km"] = np.sqrt(bwd_u_3**2 + bwd_v_3**2).to("knot")
        params["bwd_6km"] = np.sqrt(bwd_u_6**2 + bwd_v_6**2).to("knot")
    except Exception as e:
        print(f"  Warning: BWD calc failed: {e}")
        params["bwd_500m"] = 0 * units.knot
        params["bwd_1km"] = 0 * units.knot
        params["bwd_3km"] = 0 * units.knot
        params["bwd_6km"] = 0 * units.knot
    
    # STP — fixed-layer (SB CAPE, 0-1km SRH, 0-6km BWD, SB LCL)
    # MetPy's significant_tornado expects LCL HEIGHT in meters, not pressure
    try:
        _stp_lcl_h = params.get("sb_lcl_m")
        if _stp_lcl_h is not None:
            _stp_result = mpcalc.significant_tornado(
                params["sb_cape"],
                _stp_lcl_h * units.meter,
                params["srh_1km"],
                params["bwd_6km"].to("m/s")
            )
            params["stp"] = round(float(np.asarray(_stp_result.magnitude).flat[0]), 2)
        else:
            params["stp"] = 0
    except Exception as e:
        print(f"  Warning: STP calc failed: {e}")
        params["stp"] = 0
    
    # (SCP and STP-Eff moved after effective layer computations below)
    
    # ── Significant Hail Parameter (SHIP) ──
    # SHIP = (muCAPE × mixRatio × LR_7-5 × (-T500) × BWD_6km) / 42_000_000
    # Capped at 0 when muCAPE < 1300 or muCIN > -50 (per SPC guidelines)
    try:
        _mu_cape_ship = float(params.get("mu_cape", 0 * units("J/kg")).magnitude)
        _mu_cin_ship = float(params.get("mu_cin", 0 * units("J/kg")).magnitude)
        # Surface mixing ratio (g/kg)
        _sfc_mr = float(_mixing_ratio_from_dewpoint(p[0], Td[0]).to("g/kg").magnitude)
        # 700-500 hPa lapse rate
        _p700_idx = int(np.argmin(np.abs(p.magnitude - 700.0)))
        _p500_idx = int(np.argmin(np.abs(p.magnitude - 500.0)))
        _t700 = float(T.magnitude[_p700_idx])
        _t500 = float(T.magnitude[_p500_idx])
        _lr_75 = (_t700 - _t500) / ((h.magnitude[_p500_idx] - h.magnitude[_p700_idx]) / 1000.0)
        _neg_t500 = max(-_t500, 0)  # magnitude of T500 below freezing
        _bwd6_ship = float(params.get("bwd_6km", 0 * units.knot).to("m/s").magnitude)
        _ship_raw = (_mu_cape_ship * _sfc_mr * _lr_75 * _neg_t500 * _bwd6_ship) / 42_000_000.0
        # Zero out if CAPE is too low or CIN is too weak (too much inhibition)
        if _mu_cape_ship < 1300 or _mu_cin_ship < -200:
            _ship_raw = 0.0
        params["ship"] = round(max(_ship_raw, 0), 2)
    except:
        params["ship"] = 0
    
    # ── Derecho Composite Parameter (DCP) ──
    # DCP = (DCAPE/980) × (muCAPE/2000) × (BWD_6km_ms/20) × (meanWind_0-6km/16 m/s)
    try:
        _dcape_val = float(params.get("dcape", 0 * units("J/kg")).magnitude) if params.get("dcape") is not None else 0.0
        _mu_cape_dcp = float(params.get("mu_cape", 0 * units("J/kg")).magnitude)
        _bwd6_dcp = float(params.get("bwd_6km", 0 * units.knot).to("m/s").magnitude)
        # Mean wind 0-6 km from interpolated grid
        _mask_06 = h_interp.magnitude <= 6000.0
        _u06 = u_interp[_mask_06].to("m/s").magnitude
        _v06 = v_interp[_mask_06].to("m/s").magnitude
        _mean_u06 = np.mean(_u06)
        _mean_v06 = np.mean(_v06)
        _mean_wind_06 = np.sqrt(_mean_u06**2 + _mean_v06**2)
        params["dcp"] = round(
            (_dcape_val / 980.0) * (_mu_cape_dcp / 2000.0) *
            (_bwd6_dcp / 20.0) * (_mean_wind_06 / 16.0),
            2
        )
    except:
        params["dcp"] = 0

    # ── ECAPE (Entraining CAPE) — Peters et al. 2023 ─────────────────
    # Simplified analytic approximation:
    #   ECAPE ≈ MUCAPE² / (MUCAPE + 2σ² V̄²_sr)
    # where σ ≈ 1.6, V̄_sr = mean 0-6 km storm-relative wind speed
    try:
        _mu_cape_ecape = float(params.get("mu_cape", 0 * units("J/kg")).magnitude)
        if _mu_cape_ecape > 0 and params.get("rm_u") is not None:
            _rm_u_ms = params["rm_u"].to("m/s").magnitude
            _rm_v_ms = params["rm_v"].to("m/s").magnitude
            _u_interp_ms = u_interp.to("m/s").magnitude
            _v_interp_ms = v_interp.to("m/s").magnitude
            _mask_06_ecape = h_interp.magnitude <= 6000.0
            _sr_u = _u_interp_ms[_mask_06_ecape] - _rm_u_ms
            _sr_v = _v_interp_ms[_mask_06_ecape] - _rm_v_ms
            _vsr_mean_sq = np.mean(_sr_u**2 + _sr_v**2)
            _sigma = 1.6
            _ecape = _mu_cape_ecape**2 / (_mu_cape_ecape + 2.0 * _sigma**2 * _vsr_mean_sq)
            params["ecape"] = round(max(_ecape, 0), 1)
        else:
            params["ecape"] = 0
    except Exception as e:
        print(f"  Warning: ECAPE calc failed: {e}")
        params["ecape"] = 0

    # ── Effective Inflow Layer (Thompson et al. 2007) ────────────────
    # Bottom: lowest level where CAPE ≥ 100 J/kg AND CIN > -250 J/kg
    # Top:    highest level below EL meeting same criteria
    try:
        _eil_bot_p = None
        _eil_top_p = None
        _eil_bot_h = None
        _eil_top_h = None
        _sfc_h_m = h[0].magnitude
        _n_levels = len(p)
        _cape_thresh = 100.0   # J/kg
        _cin_thresh = -250.0   # J/kg
        for _i in range(_n_levels):
            if p.magnitude[_i] < 500:  # don't search above 500 hPa
                break
            try:
                _lp = mpcalc.parcel_profile(p[_i:], T[_i], Td[_i]).to("degC")
                _lcape, _lcin = mpcalc.cape_cin(p[_i:], T[_i:], Td[_i:], _lp)
                _lcape_v = float(_lcape.magnitude)
                _lcin_v = float(_lcin.magnitude)
                if _lcape_v >= _cape_thresh and _lcin_v >= _cin_thresh:
                    if _eil_bot_p is None:
                        _eil_bot_p = p.magnitude[_i]
                        _eil_bot_h = h.magnitude[_i] - _sfc_h_m
                    _eil_top_p = p.magnitude[_i]
                    _eil_top_h = h.magnitude[_i] - _sfc_h_m
                elif _eil_bot_p is not None:
                    # Exited the effective layer
                    break
            except Exception:
                continue
        params["eil_bot_p"] = _eil_bot_p
        params["eil_top_p"] = _eil_top_p
        params["eil_bot_h"] = _eil_bot_h
        params["eil_top_h"] = _eil_top_h
    except Exception as e:
        print(f"  Warning: Effective inflow layer calc failed: {e}")
        params["eil_bot_p"] = None
        params["eil_top_p"] = None
        params["eil_bot_h"] = None
        params["eil_top_h"] = None

    # ── Effective SRH (within effective inflow layer) ────────────────
    try:
        if params["eil_bot_h"] is not None and params["eil_top_h"] is not None:
            _eil_depth = params["eil_top_h"] - params["eil_bot_h"]
            if _eil_depth > 0:
                _eil_bot_agl = params["eil_bot_h"]
                _eil_top_agl = params["eil_top_h"]
                # Subset the interp grid to the effective inflow layer
                _eil_mask = (h_interp.magnitude >= _eil_bot_agl) & (h_interp.magnitude <= _eil_top_agl)
                if np.sum(_eil_mask) >= 2:
                    _h_eil = h_interp[_eil_mask]
                    _u_eil = u_interp[_eil_mask]
                    _v_eil = v_interp[_eil_mask]
                    _, _, _esrh = mpcalc.storm_relative_helicity(
                        _h_eil, _u_eil, _v_eil,
                        depth=(_eil_top_agl - _eil_bot_agl) * units.meter,
                        storm_u=params["rm_u"], storm_v=params["rm_v"]
                    )
                    params["esrh"] = _esrh
                else:
                    params["esrh"] = 0 * units("m^2/s^2")
            else:
                params["esrh"] = 0 * units("m^2/s^2")
        else:
            params["esrh"] = 0 * units("m^2/s^2")
    except Exception as e:
        print(f"  Warning: Effective SRH calc failed: {e}")
        params["esrh"] = 0 * units("m^2/s^2")

    # ── Effective BWD (shear across effective inflow layer) ──────────
    try:
        if params["eil_bot_h"] is not None and params["eil_top_h"] is not None:
            _eil_bot_agl = params["eil_bot_h"]
            # Effective BWD uses half the depth of the effective layer as the
            # "effective shear" top, but at least 1500m and no more than
            # half the EL height. SPC convention: from eil_bot to 50% of EL height (capped at ~half EL)
            _el_h = None
            if params.get("mu_el_p") is not None:
                _el_idx = np.argmin(np.abs(p.magnitude - params["mu_el_p"].magnitude))
                _el_h = h.magnitude[_el_idx] - h[0].magnitude
            if _el_h is not None and _el_h > 0:
                _ebwd_top = max(min(_el_h * 0.5, 10000.0), 1500.0)
            else:
                _ebwd_top = 6000.0  # fallback to 0-6 km
            # Get winds at EIL bottom and EBWD top
            _eb_bot_idx = np.argmin(np.abs(h_interp.magnitude - _eil_bot_agl))
            _eb_top_idx = np.argmin(np.abs(h_interp.magnitude - _ebwd_top))
            _ebwd_u = u_interp[_eb_top_idx].to("knot").magnitude - u_interp[_eb_bot_idx].to("knot").magnitude
            _ebwd_v = v_interp[_eb_top_idx].to("knot").magnitude - v_interp[_eb_bot_idx].to("knot").magnitude
            params["ebwd"] = np.sqrt(_ebwd_u**2 + _ebwd_v**2) * units.knot
        else:
            params["ebwd"] = 0 * units.knot
    except Exception as e:
        print(f"  Warning: Effective BWD calc failed: {e}")
        params["ebwd"] = 0 * units.knot

    # ── Effective-Layer STP (Thompson et al. 2012) ──────────────
    # STP_eff = (mlCAPE/1500) × (ESRH/150) × (EBWD/20 m/s) × ((2000-mlLCL)/1000) × ((mlCIN+200)/150)
    # Uses ML CAPE, effective SRH, effective BWD, ML LCL, ML CIN
    try:
        _ml_cape_v = float(params.get("ml_cape", 0 * units("J/kg")).magnitude)
        _ml_cin_v = float(params.get("ml_cin", 0 * units("J/kg")).magnitude)
        _ml_lcl_m = params.get("ml_lcl_m")
        _esrh_v = float(params.get("esrh", 0 * units("m^2/s^2")).magnitude)
        _ebwd_v = float(params.get("ebwd", 0 * units.knot).to("m/s").magnitude)

        if _ml_lcl_m is not None and _ebwd_v >= 12.5 and _ml_cape_v > 0:
            # CAPE term
            _stp_cape = _ml_cape_v / 1500.0
            # ESRH term
            _stp_esrh = _esrh_v / 150.0
            # EBWD term (capped at 1.5 for EBWD > 30 m/s, zeroed below 12.5 m/s)
            _stp_ebwd = min(_ebwd_v / 20.0, 1.5)
            # LCL term (capped at 1.0 for LCL < 1000m, zeroed for LCL > 2000m)
            if _ml_lcl_m < 1000.0:
                _stp_lcl = 1.0
            elif _ml_lcl_m > 2000.0:
                _stp_lcl = 0.0
            else:
                _stp_lcl = (2000.0 - _ml_lcl_m) / 1000.0
            # CIN term (capped at 1.0 for CIN > -50, zeroed for CIN < -200)
            if _ml_cin_v >= -50.0:
                _stp_cin = 1.0
            elif _ml_cin_v <= -200.0:
                _stp_cin = 0.0
            else:
                _stp_cin = (_ml_cin_v + 200.0) / 150.0
            params["stp_eff"] = round(_stp_cape * _stp_esrh * _stp_ebwd * _stp_lcl * _stp_cin, 2)
        else:
            params["stp_eff"] = 0
    except Exception as e:
        print(f"  Warning: STP-Eff calc failed: {e}")
        params["stp_eff"] = 0

    # ── Warm Cloud Depth (WCD) ─────────────────────────────────────
    # WCD = Freezing Level (AGL) - SB LCL height (AGL)
    # Critical for hail melting assessment and precipitation efficiency
    try:
        _frz = params.get("frz_level")
        _sb_lcl = params.get("sb_lcl_m")
        if _frz is not None and _sb_lcl is not None and _frz > _sb_lcl:
            params["wcd"] = round(_frz - _sb_lcl)
        else:
            params["wcd"] = None
    except Exception:
        params["wcd"] = None

    # ── Supercell Composite Parameter (SCP) ── Thompson et al. 2004
    # Uses MetPy built-in: SCP = (muCAPE/1000) × (ESRH/50) × (EBW/20 m/s)
    # The BWD term is capped at 1.0 (effective shear > 20 m/s) and zeroed < 10 m/s
    try:
        _scp_esrh = params.get("esrh", 0 * units("m^2/s^2"))
        _scp_ebwd = params.get("ebwd", 0 * units.knot).to("m/s")
        _scp_result = mpcalc.supercell_composite(
            params["mu_cape"], _scp_esrh, _scp_ebwd
        )
        params["scp"] = round(float(np.asarray(_scp_result.magnitude).flat[0]), 2)
    except Exception as e:
        print(f"  Warning: SCP calc failed: {e}")
        params["scp"] = 0

    # ── 3CAPE and 6CAPE (0-3 km and 0-6 km CAPE) ────────────────────
    # Computed using MU parcel profile, truncated at 3 km and 6 km AGL.
    # Uses virtual temperature correction for consistency with MetPy CAPE:
    #   Tv = T_K × (1 + 0.61r)
    # Parcel: below LCL uses surface mixing ratio; above LCL uses saturation MR
    # Environment: uses actual mixing ratio from dewpoint
    try:
        _mu_cape_3 = 0
        _mu_cape_6 = 0
        if params.get("mu_profile") is not None and params.get("mu_start_idx") is not None:
            _mu_si = params["mu_start_idx"]
            _T_env_K = (T[_mu_si:].to("degC").magnitude + 273.15)
            _T_par_K = (params["mu_profile"].to("degC").magnitude + 273.15)
            _h_mu = h[_mu_si:].magnitude
            _p_mu = p[_mu_si:]
            _sfc_h = h[0].magnitude

            # Environment mixing ratio from dewpoint
            try:
                _mr_env = _mixing_ratio_from_dewpoint(
                    _p_mu, Td[_mu_si:]
                ).magnitude  # kg/kg
            except:
                _mr_env = np.zeros(len(_T_env_K))

            # Parcel mixing ratio: surface value below LCL, saturation above
            _lcl_p_val = None
            if params.get("mu_lcl_p") is not None:
                _lcl_p_val = params["mu_lcl_p"].magnitude
            try:
                _r_sfc = _mixing_ratio_from_dewpoint(
                    p[_mu_si], Td[_mu_si]
                ).magnitude  # kg/kg (surface value)
            except:
                _r_sfc = 0.0
            try:
                _r_sat = mpcalc.saturation_mixing_ratio(
                    _p_mu, params["mu_profile"]
                ).magnitude  # kg/kg
            except:
                _r_sat = np.full(len(_T_par_K), _r_sfc)
            # Below LCL: use surface r; above LCL: use saturated r
            _mr_par = np.where(
                _p_mu.magnitude >= (_lcl_p_val if _lcl_p_val else 0),
                _r_sfc,  # below LCL (higher pressure)
                _r_sat   # above LCL (lower pressure)
            )

            # Virtual temperatures
            _Tv_env = _T_env_K * (1 + 0.61 * _mr_env)
            _Tv_par = _T_par_K * (1 + 0.61 * _mr_par)

            _g = 9.81
            for _depth_name, _depth_m in [("3", 3000), ("6", 6000)]:
                _h_top = _sfc_h + _depth_m
                _mask_d = _h_mu <= _h_top
                if np.sum(_mask_d) >= 2:
                    _buoy_d = (_Tv_par[_mask_d] - _Tv_env[_mask_d]) / _Tv_env[_mask_d]
                    _buoy_d[_buoy_d < 0] = 0  # only positive buoyancy
                    _h_d = _h_mu[_mask_d]
                    _cape_d = np.trapezoid(_buoy_d * _g, _h_d)
                    if _depth_name == "3":
                        _mu_cape_3 = max(round(float(_cape_d), 1), 0)
                    else:
                        _mu_cape_6 = max(round(float(_cape_d), 1), 0)
        params["cape_3km"] = _mu_cape_3
        params["cape_6km"] = _mu_cape_6
    except Exception as e:
        print(f"  Warning: 3CAPE/6CAPE calc failed: {e}")
        params["cape_3km"] = 0
        params["cape_6km"] = 0

    # ── DCIN (Downdraft CIN) ─────────────────────────────────────────
    # DCIN measures the inhibition of downdrafts reaching the surface.
    # Uses the downdraft parcel profile from MetPy's downdraft_cape.
    # DCIN = ∫(Tv_parcel - Tv_env) * g / Tv_env dz  (only where parcel > env, indicating
    # the downdraft is warmer than environment = inhibition for downdraft reaching surface)
    try:
        _dcin = 0
        if params.get("dcape_profile") is not None and params.get("dcape_pressure") is not None:
            _dd_p = params["dcape_pressure"].magnitude  # hPa
            _dd_T = params["dcape_profile"].to("degC").magnitude
            # Interpolate environment temperature and heights to the downdraft pressure levels
            _env_T_interp = np.interp(_dd_p, p.magnitude[::-1], T.to("degC").magnitude[::-1])
            _h_interp_dd = np.interp(_dd_p, p.magnitude[::-1], h.magnitude[::-1])
            _sfc_h = h[0].magnitude
            # Focus on the sub-cloud layer (surface to ~3 km AGL)
            _h_agl_dd = _h_interp_dd - _sfc_h
            _mask_sc = _h_agl_dd <= 3000
            if np.sum(_mask_sc) >= 2:
                _dd_sub = _dd_T[_mask_sc]
                _env_sub = _env_T_interp[_mask_sc]
                _h_sub = _h_interp_dd[_mask_sc]
                # Positive buoyancy for downdraft = inhibition (parcel warmer than env)
                _buoy_dcin = _dd_sub - _env_sub
                _pos_buoy = _buoy_dcin.copy()
                _pos_buoy[_pos_buoy < 0] = 0
                _dcin = -np.trapezoid(_pos_buoy * 9.81 / (273.15 + _env_sub), _h_sub)
                _dcin = round(float(min(_dcin, 0)), 1)
        params["dcin"] = _dcin
    except Exception as e:
        print(f"  Warning: DCIN calc failed: {e}")
        params["dcin"] = 0

    # ── MU NCAPE (Normalized CAPE) ───────────────────────────────────
    # NCAPE = MUCAPE / (EL_height - LFC_height)  [J/kg/m]
    # Measures buoyancy intensity per unit depth
    try:
        _mu_cape_ncape = float(params.get("mu_cape", 0 * units("J/kg")).magnitude)
        _ncape = 0
        if _mu_cape_ncape > 0 and params.get("mu_lfc_p") is not None and params.get("mu_el_p") is not None:
            _lfc_idx = np.argmin(np.abs(p.magnitude - params["mu_lfc_p"].magnitude))
            _el_idx = np.argmin(np.abs(p.magnitude - params["mu_el_p"].magnitude))
            _lfc_h = h.magnitude[_lfc_idx]
            _el_h = h.magnitude[_el_idx]
            _depth_m = _el_h - _lfc_h
            if _depth_m > 100:
                _ncape = round(_mu_cape_ncape / _depth_m, 3)
        params["ncape"] = _ncape
    except Exception as e:
        print(f"  Warning: NCAPE calc failed: {e}")
        params["ncape"] = 0

    # ── Piecewise CAPE (500 hPa layers from LFC to EL) ──────────────
    try:
        _pw_layers = []
        _mu_cape_pw = float(params.get("mu_cape", 0 * units("J/kg")).magnitude)
        if _mu_cape_pw > 0 and params.get("mu_profile") is not None:
            _T_env = T.magnitude
            _T_parcel = params["mu_profile"].magnitude
            _p_arr = p.magnitude
            _h_arr = h.magnitude
            # Define layers by pressure: each ~50 hPa thick from 900 to 200 hPa
            _layer_edges = list(range(900, 150, -50))
            for i in range(len(_layer_edges) - 1):
                _p_top = _layer_edges[i + 1]
                _p_bot = _layer_edges[i]
                _mask_pw = (_p_arr <= _p_bot) & (_p_arr >= _p_top)
                if np.sum(_mask_pw) < 2:
                    continue
                _buoy = _T_parcel[_mask_pw] - _T_env[_mask_pw]
                _pos = np.sum(_buoy[_buoy > 0]) * 9.81 / 273.15  # crude integration
                _neg = np.sum(_buoy[_buoy < 0]) * 9.81 / 273.15
                if abs(_pos) > 0.1 or abs(_neg) > 0.1:
                    _pw_layers.append({
                        "p_bot": _p_bot,
                        "p_top": _p_top,
                        "cape": round(float(_pos * 50.0), 1),  # scale by layer thickness
                        "cin": round(float(_neg * 50.0), 1),
                    })
        params["piecewise_cape"] = _pw_layers
    except Exception as e:
        print(f"  Warning: Piecewise CAPE calc failed: {e}")
        params["piecewise_cape"] = []

    # Freezing level (AGL)
    try:
        zero_crossings = np.where(np.diff(np.sign(T.magnitude)))[0]
        if len(zero_crossings) > 0:
            idx = zero_crossings[0]
            frac = -T.magnitude[idx] / (T.magnitude[idx+1] - T.magnitude[idx])
            frz_h_msl = h.magnitude[idx] + frac * (h.magnitude[idx+1] - h.magnitude[idx])
            params["frz_level"] = frz_h_msl - h[0].magnitude  # AGL
        else:
            params["frz_level"] = None
    except:
        params["frz_level"] = None
    
    # Precipitable water
    try:
        params["pwat"] = mpcalc.precipitable_water(p, Td)
    except:
        params["pwat"] = None
    
    # (DCAPE already computed earlier — before DCP/DCIN)
    
    # Lapse rates (Γ0-3, Γ3-6) in °C/km
    try:
        sfc_h_m = h[0].magnitude
        h_m = h.magnitude
        T_c = T.magnitude
        
        # 0-3 km lapse rate
        mask_03 = (h_m >= sfc_h_m) & (h_m <= sfc_h_m + 3000)
        if np.sum(mask_03) >= 2:
            idx_03 = np.where(mask_03)[0]
            dT_03 = T_c[idx_03[0]] - T_c[idx_03[-1]]
            dZ_03 = (h_m[idx_03[-1]] - h_m[idx_03[0]]) / 1000.0
            params["lr_03"] = dT_03 / dZ_03 if dZ_03 > 0 else None
        else:
            params["lr_03"] = None
        
        # 3-6 km lapse rate
        mask_36 = (h_m >= sfc_h_m + 3000) & (h_m <= sfc_h_m + 6000)
        if np.sum(mask_36) >= 2:
            idx_36 = np.where(mask_36)[0]
            dT_36 = T_c[idx_36[0]] - T_c[idx_36[-1]]
            dZ_36 = (h_m[idx_36[-1]] - h_m[idx_36[0]]) / 1000.0
            params["lr_36"] = dT_36 / dZ_36 if dZ_36 > 0 else None
        else:
            params["lr_36"] = None
    except:
        params["lr_03"] = None
        params["lr_36"] = None
    
    # Wet-bulb zero height (WBO) — height AGL where wet-bulb = 0°C
    try:
        if params.get("wetbulb") is not None:
            wb = params["wetbulb"].magnitude
            wb_crossings = np.where(np.diff(np.sign(wb)))[0]
            if len(wb_crossings) > 0:
                idx_wb = wb_crossings[0]
                frac_wb = -wb[idx_wb] / (wb[idx_wb+1] - wb[idx_wb])
                wbo_msl = h_m[idx_wb] + frac_wb * (h_m[idx_wb+1] - h_m[idx_wb])
                params["wbo"] = wbo_msl - h[0].magnitude
            else:
                params["wbo"] = None
        else:
            params["wbo"] = None
    except:
        params["wbo"] = None
    
    # Relative humidity layers
    try:
        rh = mpcalc.relative_humidity_from_dewpoint(T, Td) * 100
        params["rh"] = rh
        
        # Layer averages
        sfc_h = h[0].magnitude
        for layer_name, bot, top in [
            ("rh_0_1km", 0, 1000),
            ("rh_1_3km", 1000, 3000),
            ("rh_3_6km", 3000, 6000),
        ]:
            mask = (h.magnitude >= sfc_h + bot) & (h.magnitude <= sfc_h + top)
            if np.any(mask):
                params[layer_name] = np.mean(rh.magnitude[mask])
            else:
                params[layer_name] = None
    except:
        params["rh"] = None

    # ── Winter Weather / Precip Type (Bourgouin Method) ──────────────
    # Uses warm-nose and cold-layer energy areas above/below 0°C to
    # classify precipitation type: Rain / Snow / Ice Pellets / Freezing Rain
    try:
        _T_c_pt = T.to("degC").magnitude
        _h_m_pt = h.magnitude
        _sfc_T_pt = _T_c_pt[0]

        # Find all freezing level crossings (sign changes)
        _sign_pt = np.sign(_T_c_pt)
        _crossings_pt = np.where(np.diff(_sign_pt))[0]

        _warm_area = 0.0  # J/kg (warm nose positive area)
        _cold_area = 0.0  # J/kg (cold layer negative area)
        _precip_type = "N/A"

        if len(_crossings_pt) >= 2:
            # Warm nose: integrate positive T between 1st and 2nd crossing
            _wn_bot = _crossings_pt[0]
            _wn_top = _crossings_pt[1]
            for _k in range(_wn_bot, _wn_top):
                _dz = _h_m_pt[_k + 1] - _h_m_pt[_k]
                _T_avg = (_T_c_pt[_k] + _T_c_pt[_k + 1]) / 2.0
                if _T_avg > 0:
                    _warm_area += 9.81 * _T_avg / 273.15 * _dz

            # Cold layer: from surface to 1st crossing
            for _k in range(0, _crossings_pt[0]):
                _dz = _h_m_pt[_k + 1] - _h_m_pt[_k]
                _T_avg = (_T_c_pt[_k] + _T_c_pt[_k + 1]) / 2.0
                if _T_avg < 0:
                    _cold_area += 9.81 * abs(_T_avg) / 273.15 * _dz

            # Bourgouin (2000) classification thresholds
            if _warm_area < 5.6:
                _precip_type = "Snow"
            elif _warm_area >= 13.2:
                if _cold_area < 5.6:
                    _precip_type = "Rain"
                elif _cold_area >= 13.2:
                    _precip_type = "Ice Pellets"
                else:
                    _precip_type = "Freezing Rain"
            else:
                if _cold_area < 5.6:
                    _precip_type = "Freezing Rain"
                elif _cold_area >= 66:
                    _precip_type = "Ice Pellets"
                else:
                    _precip_type = "Frzg Rain/Sleet"
        elif len(_crossings_pt) == 1:
            _precip_type = "Rain" if _sfc_T_pt > 0 else "Snow"
        else:
            _precip_type = "Rain" if _sfc_T_pt > 0 else "Snow"

        params["precip_type"] = _precip_type
        params["warm_layer_energy"] = round(_warm_area, 1)
        params["cold_layer_energy"] = round(_cold_area, 1)
    except Exception as e:
        print(f"  Warning: Precip type calc failed: {e}")
        params["precip_type"] = "N/A"
        params["warm_layer_energy"] = 0
        params["cold_layer_energy"] = 0

    # ── Microburst / Downburst Composites ────────────────────────────
    # WMSI: Wet Microburst Severity Index ≈ CAPE × Γ0-3 / 1000
    # MDPI: Microburst Day Potential Index = (θe_sfc − θe_min_0-6km) / 20
    # Max Gust Potential: V = √(2 × DCAPE) converted to knots
    try:
        _mu_cape_mb = float(params.get("mu_cape", 0 * units("J/kg")).magnitude)
        _dcape_mb = float(params.get("dcape", 0 * units("J/kg")).magnitude) if params.get("dcape") is not None else 0

        # WMSI
        _lr03_mb = params.get("lr_03", 0)
        if _lr03_mb is None:
            _lr03_mb = 0
        _wmsi = (_mu_cape_mb * max(float(_lr03_mb), 0) / 1000.0) if _mu_cape_mb > 0 else 0
        params["wmsi"] = round(_wmsi, 1)

        # MDPI: θe deficit surface to min in 0-6 km
        _sfc_h_mb = h[0].magnitude
        _mask_06_mb = (h.magnitude >= _sfc_h_mb) & (h.magnitude <= _sfc_h_mb + 6000)
        if np.sum(_mask_06_mb) >= 2:
            _theta_e_all = mpcalc.equivalent_potential_temperature(p, T, Td)
            _theta_e_06 = _theta_e_all[_mask_06_mb].magnitude
            _theta_e_sfc = _theta_e_all[0].magnitude
            _theta_e_min = np.min(_theta_e_06)
            _mdpi = (_theta_e_sfc - _theta_e_min) / 20.0
            params["mdpi"] = round(max(float(_mdpi), 0), 2)
        else:
            params["mdpi"] = 0

        # Max downburst gust potential
        if _dcape_mb > 0:
            _max_gust_ms = np.sqrt(2.0 * _dcape_mb)
            params["max_gust"] = round(float(_max_gust_ms * 1.94384), 1)
        else:
            params["max_gust"] = 0
    except Exception as e:
        print(f"  Warning: Microburst composites calc failed: {e}")
        params["wmsi"] = 0
        params["mdpi"] = 0
        params["max_gust"] = 0

    # ── Corfidi Vectors (Corfidi 2003) ───────────────────────────────
    # V_CL  = mean wind 850-300 hPa (cloud-layer mean)
    # V_LLJ = max wind in 0-1.5 km AGL (low-level jet)
    # Corfidi Upwind  = V_CL − V_LLJ  (back-building / upwind propagation)
    # Corfidi Downwind = 2·V_CL − V_LLJ (forward / downwind propagation)
    try:
        # Cloud-layer mean wind (850-300 hPa)
        _cl_mask = (p_interp.magnitude >= 300) & (p_interp.magnitude <= 850)
        if np.sum(_cl_mask) >= 2:
            _cl_u = float(np.mean(u_interp[_cl_mask].to("knot").magnitude))
            _cl_v = float(np.mean(v_interp[_cl_mask].to("knot").magnitude))
        else:
            _cl_u = float(params["mw_u"].to("knot").magnitude)
            _cl_v = float(params["mw_v"].to("knot").magnitude)

        # LLJ: max wind speed in 0-1.5 km AGL
        _llj_mask = h_interp.magnitude <= 1500
        if np.sum(_llj_mask) >= 2:
            _u_llj_all = u_interp[_llj_mask].to("knot").magnitude
            _v_llj_all = v_interp[_llj_mask].to("knot").magnitude
            _spd_llj = np.sqrt(_u_llj_all**2 + _v_llj_all**2)
            _max_idx = int(np.argmax(_spd_llj))
            _llj_u = float(_u_llj_all[_max_idx])
            _llj_v = float(_v_llj_all[_max_idx])
        else:
            _llj_u = float(u_interp[0].to("knot").magnitude)
            _llj_v = float(v_interp[0].to("knot").magnitude)

        params["corfidi_up_u"] = round(_cl_u - _llj_u, 1)
        params["corfidi_up_v"] = round(_cl_v - _llj_v, 1)
        params["corfidi_dn_u"] = round(2 * _cl_u - _llj_u, 1)
        params["corfidi_dn_v"] = round(2 * _cl_v - _llj_v, 1)
        params["corfidi_up_spd"] = round(np.sqrt(params["corfidi_up_u"]**2 + params["corfidi_up_v"]**2), 1)
        params["corfidi_dn_spd"] = round(np.sqrt(params["corfidi_dn_u"]**2 + params["corfidi_dn_v"]**2), 1)
    except Exception as e:
        print(f"  Warning: Corfidi vectors calc failed: {e}")
        params["corfidi_up_u"] = None
        params["corfidi_up_v"] = None
        params["corfidi_dn_u"] = None
        params["corfidi_dn_v"] = None
        params["corfidi_up_spd"] = None
        params["corfidi_dn_spd"] = None

    # ── Fire Weather Indices ─────────────────────────────────────────
    # Fosberg FWI: surface T, RH, and wind speed
    # Haines Index: mid-level stability + moisture (850-700 hPa)
    # Hot-Dry-Windy (HDW): max(VPD × wind) in lowest 500 m
    try:
        _sfc_T_fw = float(T[0].to("degC").magnitude)
        _sfc_rh_fw = float(mpcalc.relative_humidity_from_dewpoint(T[0], Td[0]).magnitude * 100)
        _sfc_wspd_mph = float(np.sqrt(u[0]**2 + v[0]**2).to("mph").magnitude)

        # Fosberg FWI — equilibrium moisture content
        if _sfc_rh_fw <= 10:
            _m_fw = 0.03229 + 0.281073 * _sfc_rh_fw - 0.000578 * _sfc_rh_fw * _sfc_T_fw
        elif _sfc_rh_fw <= 50:
            _m_fw = 2.22749 + 0.160107 * _sfc_rh_fw - 0.01478 * _sfc_T_fw
        else:
            _m_fw = 21.0606 + 0.005565 * _sfc_rh_fw**2 - 0.00035 * _sfc_rh_fw * _sfc_T_fw - 0.483199 * _sfc_rh_fw
        _m_fw = max(_m_fw, 0.1)
        _eta_fw = 1.0 - 2.0 * (_m_fw / 30.0) + 1.5 * (_m_fw / 30.0)**2 - 0.5 * (_m_fw / 30.0)**3
        _fwi_val = _eta_fw * np.sqrt(1 + _sfc_wspd_mph**2) / 0.3002
        params["fosberg_fwi"] = round(float(min(max(_fwi_val, 0), 100)), 1)

        # Haines Index (mid-level variant: 850-700 hPa)
        _p850i = int(np.argmin(np.abs(p.magnitude - 850)))
        _p700i = int(np.argmin(np.abs(p.magnitude - 700)))
        _T850_fw = float(T[_p850i].to("degC").magnitude)
        _T700_fw = float(T[_p700i].to("degC").magnitude)
        _Td850_fw = float(Td[_p850i].to("degC").magnitude)
        _stab_fw = _T850_fw - _T700_fw
        _A_fw = 1 if _stab_fw < 6 else (2 if _stab_fw <= 10 else 3)
        _depr_fw = _T850_fw - _Td850_fw
        _B_fw = 1 if _depr_fw < 6 else (2 if _depr_fw <= 12 else 3)
        params["haines"] = _A_fw + _B_fw  # Range: 2-6

        # Hot-Dry-Windy Index (Srock et al. 2018)
        _sfc_h_fw = h[0].magnitude
        _mask_500m_fw = h.magnitude <= _sfc_h_fw + 500
        if np.sum(_mask_500m_fw) >= 2:
            _T_fw500 = T[_mask_500m_fw].to("degC").magnitude
            _Td_fw500 = Td[_mask_500m_fw].to("degC").magnitude
            _wspd_fw500 = np.sqrt(u[_mask_500m_fw]**2 + v[_mask_500m_fw]**2).to("m/s").magnitude
            _es_T_fw = 6.112 * np.exp(17.67 * _T_fw500 / (_T_fw500 + 243.5))
            _es_Td_fw = 6.112 * np.exp(17.67 * _Td_fw500 / (_Td_fw500 + 243.5))
            _vpd_fw = _es_T_fw - _es_Td_fw
            params["hdw"] = round(float(np.max(_vpd_fw * _wspd_fw500)), 1)
        else:
            params["hdw"] = 0
    except Exception as e:
        print(f"  Warning: Fire weather calc failed: {e}")
        params["fosberg_fwi"] = None
        params["haines"] = None
        params["hdw"] = None

    # ── Sounding-Derived Hazard Classification ───────────────────────
    # Auto-scores: TORNADO, HAIL, WIND, FLOOD with LOW / MOD / HIGH
    try:
        _hazards = []
        _mu_cape_hz = float(params.get("mu_cape", 0 * units("J/kg")).magnitude)
        _stp_hz = params.get("stp", 0)
        _stp_eff_hz = params.get("stp_eff", 0)
        _scp_hz = params.get("scp", 0)
        _ship_hz = params.get("ship", 0)
        _dcp_hz = params.get("dcp", 0)
        _dcape_hz = float(params.get("dcape", 0 * units("J/kg")).magnitude) if params.get("dcape") is not None else 0
        _bwd6_hz = float(params.get("bwd_6km", 0 * units.knot).magnitude)
        _srh1_hz = float(params.get("srh_1km", 0 * units("m^2/s^2")).magnitude)
        _ml_lcl_hz = params.get("ml_lcl_m", 9999)
        if _ml_lcl_hz is None:
            _ml_lcl_hz = 9999
        _pw_hz = float(params.get("pwat", 0 * units.mm).magnitude) if params.get("pwat") is not None else 0

        # Tornado
        _tor = 0
        if _stp_eff_hz >= 4:
            _tor = 3
        elif _stp_eff_hz >= 1:
            _tor = 2
        elif _stp_hz >= 1 or (_srh1_hz >= 150 and _ml_lcl_hz < 1500 and _mu_cape_hz >= 500):
            _tor = 1
        if _tor > 0:
            _hazards.append({"type": "TORNADO", "level": ["LOW", "MOD", "HIGH"][_tor - 1]})

        # Hail
        _hal = 0
        if _ship_hz >= 3:
            _hal = 3
        elif _ship_hz >= 1.5:
            _hal = 2
        elif _ship_hz >= 0.5 or (_scp_hz >= 2 and _mu_cape_hz >= 1500):
            _hal = 1
        if _hal > 0:
            _hazards.append({"type": "HAIL", "level": ["LOW", "MOD", "HIGH"][_hal - 1]})

        # Wind
        _wnd = 0
        if _dcp_hz >= 6 or (_dcape_hz >= 1200 and _bwd6_hz >= 40):
            _wnd = 3
        elif _dcp_hz >= 3 or _dcape_hz >= 800:
            _wnd = 2
        elif _dcp_hz >= 1 or _dcape_hz >= 400:
            _wnd = 1
        if _wnd > 0:
            _hazards.append({"type": "WIND", "level": ["LOW", "MOD", "HIGH"][_wnd - 1]})

        # Flood
        _fld = 0
        if _pw_hz >= 50 and _mu_cape_hz >= 1000:
            _fld = 3
        elif _pw_hz >= 40 and _mu_cape_hz >= 500:
            _fld = 2
        elif _pw_hz >= 30:
            _fld = 1
        if _fld > 0:
            _hazards.append({"type": "FLOOD", "level": ["LOW", "MOD", "HIGH"][_fld - 1]})

        params["hazards"] = _hazards
    except Exception as e:
        print(f"  Warning: Hazard classification failed: {e}")
        params["hazards"] = []

    # ── Temperature Advection Profile ────────────────────────────────
    # Wind veering with height = warm-air advection (WAA)
    # Wind backing with height = cold-air advection (CAA)
    # Assessed in 1 km layers from 0-6 km AGL
    try:
        _adv_layers = []
        _u_ms_adv = u_interp.to("m/s").magnitude
        _v_ms_adv = v_interp.to("m/s").magnitude
        _h_agl_adv = h_interp.magnitude

        for _bot_km in range(0, 6):
            _top_km = _bot_km + 1
            _bot_m = _bot_km * 1000
            _top_m = _top_km * 1000
            _bot_idx = int(np.argmin(np.abs(_h_agl_adv - _bot_m)))
            _top_idx = int(np.argmin(np.abs(_h_agl_adv - _top_m)))

            if _bot_idx == _top_idx:
                continue

            _dir_bot = float(np.degrees(np.arctan2(-_u_ms_adv[_bot_idx], -_v_ms_adv[_bot_idx])) % 360)
            _dir_top = float(np.degrees(np.arctan2(-_u_ms_adv[_top_idx], -_v_ms_adv[_top_idx])) % 360)

            _delta = (_dir_top - _dir_bot + 540) % 360 - 180

            if abs(_delta) < 5:
                _adv_type = "NEUTRAL"
            elif _delta > 0:
                _adv_type = "WAA"
            else:
                _adv_type = "CAA"

            _adv_layers.append({
                "layer": f"{_bot_km}-{_top_km} km",
                "type": _adv_type,
                "turn": round(float(_delta), 1)
            })

        params["temp_advection"] = _adv_layers
    except Exception as e:
        print(f"  Warning: Temperature advection calc failed: {e}")
        params["temp_advection"] = []

    return params


# ─────────────────────────────────────────────────────────────────────
# PLOTTING
# ─────────────────────────────────────────────────────────────────────
def plot_sounding(data, params, station_id, dt, vad_data=None, sr_hodograph=False,
                  theme="dark", colorblind=False, boundary_orientation=None,
                  map_zoom=1.0):
    """Create a comprehensive sounding analysis figure.

    Parameters
    ----------
    sr_hodograph : bool
        If True, plot the hodograph in storm-relative frame (subtract storm
        motion from all winds so SM is at the origin).
    theme : str
        'dark' (default) or 'light'.
    colorblind : bool
        If True, use color-blind safe palette.
    boundary_orientation : float or None
        Meteorological orientation of a surface boundary in degrees (0-360).
        This draws a dashed line across the hodograph at the given angle,
        representing the boundary orientation (e.g. an outflow boundary,
        front, or dryline).  The line is drawn perpendicular to the
        boundary-normal vector.
    map_zoom : float
        Zoom factor for CONUS mini-map inset (1.0 = full CONUS, >1 = zoomed in
        on station, max ~8).  Higher values narrow the map extent around the
        station location.
    """
    p = data["pressure"]
    T = data["temperature"]
    Td = data["dewpoint"]
    h = data["height"]
    wdir = data["wind_direction"]
    wspd = data["wind_speed"]
    u, v = params["u"], params["v"]
    info = data.get("station_info", {})
    
    station_name = STATIONS.get(station_id, (station_id, 0, 0))[0]
    lat = info.get("lat", STATIONS.get(station_id, ("", 0, 0))[1])
    lon = info.get("lon", STATIONS.get(station_id, ("", 0, 0))[2])
    
    # ── Theme colors ─────────────────────────────────────────────────
    if theme == "light":
        BG        = "#f5f5f5"
        BG_PANEL  = "#e8e8e8"
        FG        = "#1a1a1a"
        FG_DIM    = "#444444"
        FG_FAINT  = "#888888"
        GRID_CLR  = "#cccccc"
        BORDER    = "#aaaaaa"
        ACCENT    = "#2563eb"
    else:
        BG        = "#0d0d0d"
        BG_PANEL  = "#141414"
        FG        = "#e8e8e8"
        FG_DIM    = "#b0b0b0"
        FG_FAINT  = "#707070"
        GRID_CLR  = "#333333"
        BORDER    = "#444444"
        ACCENT    = "#55bbee"

    # ── Trace colors (standard vs color-blind safe) ──────────────────
    if colorblind:
        # Wong 2011 / Okabe-Ito palette
        CLR_TEMP      = "#D55E00"   # vermillion (temperature)
        CLR_DEW       = "#0072B2"   # blue (dewpoint)
        CLR_WETBULB   = "#009E73"   # bluish green (wetbulb)
        CLR_VTEMP     = "#E69F00"   # orange (virtual temp)
        CLR_SB_PARCEL = "#CC79A7"   # reddish purple (SB parcel)
        CLR_MU_PARCEL = "#F0E442"   # yellow (MU parcel)
        CLR_ML_PARCEL = "#56B4E9"   # sky blue (ML parcel)
        CLR_CAPE_FILL = "#D55E00"   # vermillion
        CLR_CIN_FILL  = "#0072B2"   # blue
        CLR_LCL       = "#009E73"   # bluish green
        CLR_LFC       = "#E69F00"   # orange
        CLR_EL        = "#56B4E9"   # sky blue
        CLR_VAD       = "#009E73"   # green
    else:
        CLR_TEMP      = "red"
        CLR_DEW       = "blue"
        CLR_WETBULB   = "cyan"
        CLR_VTEMP     = "red"
        CLR_SB_PARCEL = "#ff8800"
        CLR_MU_PARCEL = FG
        CLR_ML_PARCEL = "#dd44dd"
        CLR_CAPE_FILL = "#ff3333"
        CLR_CIN_FILL  = "#4488ff"
        CLR_LCL       = "#44cc44"
        CLR_LFC       = "#ddaa22"
        CLR_EL        = "#4499ee"
        CLR_VAD       = "#00ff88"
    
    # Helper
    def fv(val, unit_str="", decimals=0):
        if val is None:
            return "---"
        try:
            v = val.magnitude if hasattr(val, "magnitude") else val
            if decimals == 0:
                return f"{v:.0f}{' ' + unit_str if unit_str else ''}"
            else:
                return f"{v:.{decimals}f}{' ' + unit_str if unit_str else ''}"
        except:
            return "---"
    
    # ── Create figure ────────────────────────────────────────────────
    fig = plt.figure(figsize=(28, 12), facecolor=BG)
    fig.patch.set_facecolor(BG)
    
    # Define grid for layout — 5 columns: SkewT | Hodograph | SRW | SRH | Theta
    gs = gridspec.GridSpec(
        2, 5, figure=fig,
        width_ratios=[1.8, 1.6, 0.36, 0.38, 0.42],
        height_ratios=[3.0, 1.0],
        hspace=0.12, wspace=0.12,
        left=0.05, right=0.97, top=0.94, bottom=0.04
    )
    
    # ── TITLE BAR ────────────────────────────────────────────────────
    title_tags = []
    if params.get("smoothing_applied"):
        title_tags.append("SMOOTHED")
    if sr_hodograph:
        title_tags.append("SR HODO")
    tag_str = f" [{', '.join(title_tags)}]" if title_tags else ""
    title_str = (
        f"OBSERVED UPPER-AIR SOUNDING | {station_id} | "
        f"VALID: {dt.strftime('%m/%d/%Y %HZ')}{tag_str}"
    )
    fig.suptitle(title_str, fontsize=20, fontweight="bold",
                 color=FG, y=0.985, x=0.36, ha="center",
                 fontfamily="monospace")
    
    # (Station info placed in bottom-right footer)
    
    # ════════════════════════════════════════════════════════════════
    # SKEW-T LOG-P DIAGRAM
    # ════════════════════════════════════════════════════════════════
    skew = SkewT(fig, rotation=40, subplot=gs[0, 0])
    skew.ax.set_facecolor(BG)
    skew.ax.set_aspect('auto')  # override MetPy's fixed aspect so SkewT fills the grid cell
    
    for spine in skew.ax.spines.values():
        spine.set_color(BORDER)
    skew.ax.tick_params(colors=FG_DIM, labelsize=11, width=1.2)
    skew.ax.set_xlabel("Temperature (°C)", color=FG_DIM, fontsize=12,
                       fontfamily="monospace", fontweight="bold")
    skew.ax.set_ylabel("Pressure (hPa)", color=FG_DIM, fontsize=12,
                       fontfamily="monospace", fontweight="bold")
    
    # Reference lines (subtle)
    skew.plot_dry_adiabats(colors=GRID_CLR, alpha=0.4, linewidth=0.5)
    skew.plot_moist_adiabats(colors=GRID_CLR, alpha=0.4, linewidth=0.5)
    skew.plot_mixing_lines(colors=GRID_CLR, alpha=0.25, linewidth=0.4)
    
    # Panel label
    skew.ax.text(0.02, 0.98, "SKEW-T LOG-P", transform=skew.ax.transAxes,
                 fontsize=13, color=FG, fontfamily="monospace",
                 fontweight="bold", va="top", ha="left", alpha=0)  # hidden, merged into legend title
    
    # Temperature (red, solid thick)
    skew.plot(p, T, color="red", linewidth=4.0, zorder=6, label="TEMPERATURE")
    # Dewpoint (blue, solid thick)
    skew.plot(p, Td, color="blue", linewidth=4.0, zorder=6, label="DEWPOINT")
    
    # Wet-bulb temperature (cyan, solid thin)
    if params.get("wetbulb") is not None:
        skew.plot(p, params["wetbulb"], color=CLR_WETBULB, linewidth=1.5,
                  linestyle="-", alpha=0.85, zorder=5, label="WETBULB TEMP")
    
    # Virtual temperature (red, dotted)
    if params.get("virtual_temp") is not None:
        skew.plot(p, params["virtual_temp"], color=CLR_VTEMP, linewidth=1.5,
                  linestyle=":", alpha=0.7, zorder=4, label="VIRTUAL TEMP")
    
    # Downdraft parcel trace (gray, dashed)
    if params.get("dcape_profile") is not None:
        _dcape_p = params.get("dcape_pressure", p)
        skew.plot(_dcape_p, params["dcape_profile"], color="gray", linewidth=1.5,
                  linestyle="--", alpha=0.85, zorder=7, label="DWNDRFT PARCEL")
    
    # SB parcel trace (orange, dashed) + CAPE/CIN shading
    if params.get("sb_profile") is not None:
        sb_prof = params["sb_profile"]
        skew.plot(p, sb_prof, color=CLR_SB_PARCEL, linewidth=2.0,
                  linestyle="--", alpha=0.85, zorder=5, label="SB PARCEL")
        # --- CAPE shading (red) — parcel warmer than environment ---
        try:
            sb_T = sb_prof.to("degC").magnitude
            env_T = T.to("degC").magnitude
            cape_mask = sb_T > env_T
            skew.ax.fill_betweenx(
                p.magnitude, sb_T, env_T,
                where=cape_mask, facecolor=CLR_CAPE_FILL, alpha=0.18,
                interpolate=True, zorder=3
            )
            # --- CIN shading (blue) — parcel cooler than environment ---
            cin_mask = sb_T < env_T
            skew.ax.fill_betweenx(
                p.magnitude, sb_T, env_T,
                where=cin_mask, facecolor=CLR_CIN_FILL, alpha=0.12,
                interpolate=True, zorder=3
            )
        except Exception:
            pass
    
    # MU parcel trace (white/light, dashed)
    if params.get("mu_profile") is not None and params.get("mu_start_idx") is not None:
        mu_si = params["mu_start_idx"]
        skew.plot(p[mu_si:], params["mu_profile"], color=CLR_MU_PARCEL, linewidth=2.0,
                  linestyle="--", alpha=0.75, zorder=5, label="MU PARCEL")
    
    # ML parcel trace (magenta, dashed)
    if params.get("ml_profile") is not None:
        skew.plot(p, params["ml_profile"], color=CLR_ML_PARCEL, linewidth=2.0,
                  linestyle="--", alpha=0.75, zorder=5, label="ML PARCEL")
    
    # --- Highlighted isotherms: 0°C and -20°C ---
    skew.ax.axvline(0, color="#44bbee", linestyle="--", alpha=0.5,
                    linewidth=1.2, zorder=2)
    skew.ax.axvline(-20, color="#8888ff", linestyle="--", alpha=0.4,
                    linewidth=1.0, zorder=2)
    
    # --- Surface T and Td labels in °F ---
    try:
        sfc_T_F = T[0].to("degF").magnitude
        sfc_Td_F = Td[0].to("degF").magnitude
        skew.ax.annotate(
            f"{sfc_T_F:.0f}°F", xy=(T[0].magnitude, p[0].magnitude),
            xytext=(8, -15), textcoords="offset points",
            fontsize=11, color=CLR_TEMP, fontweight="bold",
            fontfamily="monospace", zorder=10,
            path_effects=[path_effects.withStroke(linewidth=3, foreground=BG)]
        )
        skew.ax.annotate(
            f"{sfc_Td_F:.0f}°F", xy=(Td[0].magnitude, p[0].magnitude),
            xytext=(-30, -15), textcoords="offset points",
            fontsize=11, color=CLR_DEW, fontweight="bold",
            fontfamily="monospace", zorder=10,
            path_effects=[path_effects.withStroke(linewidth=3, foreground=BG)]
        )
    except Exception:
        pass
    
    # --- Wet-Bulb Zero (WBZ) level annotation ---
    if params.get("wetbulb") is not None:
        try:
            wb = params["wetbulb"].magnitude
            # Find where wetbulb crosses 0°C
            for i in range(len(wb) - 1):
                if wb[i] >= 0 and wb[i+1] < 0:
                    # Linear interpolation for exact crossing pressure
                    frac = wb[i] / (wb[i] - wb[i+1])
                    wbz_p = p.magnitude[i] + frac * (p.magnitude[i+1] - p.magnitude[i])
                    wbz_h_msl = h.magnitude[i] + frac * (h.magnitude[i+1] - h.magnitude[i])
                    wbz_h_agl = wbz_h_msl - h[0].magnitude
                    skew.ax.axhline(y=wbz_p, color=CLR_WETBULB, linestyle=":",
                                    alpha=0.4, linewidth=0.8)
                    skew.ax.text(
                        skew.ax.get_xlim()[0] + 5, wbz_p,
                        f"WBZ ({wbz_h_agl:.0f}m)",
                        color=CLR_WETBULB, fontsize=10, va="center",
                        fontfamily="monospace", fontweight="bold",
                        path_effects=[path_effects.withStroke(linewidth=3, foreground=BG)]
                    )
                    break
        except Exception:
            pass
    
    # Wind barbs — sounding obs column + optional VAD column, both inside Skew-T
    has_vad = vad_data and isinstance(vad_data, dict) and vad_data.get("winds")
    barb_interval = max(1, len(p) // 40)
    obs_xloc = 0.98 if has_vad else 0.98
    try:
        obs_barbs = skew.plot_barbs(
            p[::barb_interval], u[::barb_interval], v[::barb_interval],
            color=FG, length=6, linewidth=0.8,
            xloc=obs_xloc, x_clip_radius=0.08
        )
        if obs_barbs:
            obs_barbs.set_clip_on(True)
            obs_barbs.set_clip_box(skew.ax.bbox)
    except:
        pass

    # VAD wind barbs on Skew-T (left column, green)
    if has_vad:
        try:
            vad_winds = vad_data["winds"]
            h_msl = h.to("meter").magnitude
            p_mag = p.magnitude
            sfc_h = h_msl[0]

            vad_p_list, vad_u_list, vad_v_list = [], [], []
            for w in vad_winds:
                alt_agl_m = w["alt_m"]
                alt_msl_m = alt_agl_m + sfc_h
                if alt_msl_m < h_msl[0] or alt_msl_m > h_msl[-1]:
                    continue
                p_interp_val = np.interp(alt_msl_m, h_msl, p_mag)
                if p_interp_val < 100:
                    continue
                vad_p_list.append(p_interp_val)
                vad_u_list.append(w["u_kt"])
                vad_v_list.append(w["v_kt"])

            if vad_p_list:
                VAD_COLOR = CLR_VAD
                vad_p_arr = np.array(vad_p_list) * units.hPa
                vad_u_arr = np.array(vad_u_list) * units.knot
                vad_v_arr = np.array(vad_v_list) * units.knot
                vad_barbs = skew.plot_barbs(
                    vad_p_arr, vad_u_arr, vad_v_arr,
                    color=VAD_COLOR, length=6, linewidth=0.8,
                    xloc=0.91, x_clip_radius=0.08
                )
                if vad_barbs:
                    vad_barbs.set_clip_on(True)
                    vad_barbs.set_clip_box(skew.ax.bbox)
        except Exception as e:
            print(f"[PLOT] VAD barbs on Skew-T failed: {e}")

    # Annotate key levels
    pe = path_effects.withStroke(linewidth=3, foreground=BG)
    
    if params.get("sb_lcl_p") is not None:
        skew.ax.axhline(y=params["sb_lcl_p"].magnitude, color=CLR_LCL,
                        linestyle="--", alpha=0.6, linewidth=1.0)
        skew.ax.text(
            skew.ax.get_xlim()[1] - 2, params["sb_lcl_p"].magnitude,
            f"←SBLCL ({fv(params['sb_lcl_p'],'hPa')})",
            color=CLR_LCL, fontsize=12, va="center",
            fontfamily="monospace", fontweight="bold", path_effects=[pe]
        )
    
    if params.get("sb_lfc_p") is not None:
        skew.ax.axhline(y=params["sb_lfc_p"].magnitude, color=CLR_LFC,
                        linestyle="--", alpha=0.5, linewidth=1.0)
        skew.ax.text(
            skew.ax.get_xlim()[1] - 2, params["sb_lfc_p"].magnitude,
            f"←LFC", color=CLR_LFC, fontsize=12, va="center",
            fontfamily="monospace", fontweight="bold", path_effects=[pe]
        )
    
    if params.get("sb_el_p") is not None:
        skew.ax.axhline(y=params["sb_el_p"].magnitude, color=CLR_EL,
                        linestyle="--", alpha=0.5, linewidth=1.0)
        skew.ax.text(
            skew.ax.get_xlim()[1] - 2, params["sb_el_p"].magnitude,
            f"←EL", color=CLR_EL, fontsize=12, va="center",
            fontfamily="monospace", fontweight="bold", path_effects=[pe]
        )
    
    if params.get("frz_level") is not None:
        # Find pressure at freezing level (frz_level is AGL, h is MSL)
        frz_h_agl = params["frz_level"]
        frz_h_msl = frz_h_agl + h[0].magnitude
        frz_idx = np.argmin(np.abs(h.magnitude - frz_h_msl))
        frz_p = p.magnitude[frz_idx]
        skew.ax.axhline(y=frz_p, color="#44bbee", linestyle=":",
                        alpha=0.5, linewidth=1.0)
        skew.ax.text(
            skew.ax.get_xlim()[1] - 2, frz_p,
            f"←FRZ ({frz_h_agl:.0f}m)", color="#44bbee",
            fontsize=12, va="center", fontfamily="monospace",
            fontweight="bold", path_effects=[pe]
        )
    
    # --- Dendritic Growth Zone (DGZ) shading: -12°C to -17°C ---
    # --- Piecewise CAPE colored bars on left side of Skew-T ---
    try:
        _pw_data = params.get("piecewise_cape", [])
        if _pw_data and len(_pw_data) > 0:
            _max_pw_cape = max((lyr["cape"] for lyr in _pw_data), default=1)
            if _max_pw_cape > 0:
                _xlim = skew.ax.get_xlim()
                _bar_x_start = _xlim[0] + 1.0
                _bar_max_width = 10.0  # max width in °C units
                for _lyr in _pw_data:
                    if _lyr["cape"] > 0:
                        _bar_width = (_lyr["cape"] / _max_pw_cape) * _bar_max_width
                        # Map CAPE magnitude to color: green → yellow → red
                        _frac = min(_lyr["cape"] / max(_max_pw_cape, 500), 1.0)
                        if _frac < 0.5:
                            _r = int(255 * _frac * 2)
                            _g = 255
                        else:
                            _r = 255
                            _g = int(255 * (1 - (_frac - 0.5) * 2))
                        _bar_color = f"#{_r:02x}{_g:02x}00"
                        skew.ax.barh(
                            (_lyr["p_bot"] + _lyr["p_top"]) / 2.0,
                            _bar_width,
                            height=(_lyr["p_bot"] - _lyr["p_top"]),
                            left=_bar_x_start,
                            color=_bar_color, alpha=0.35, zorder=2,
                            edgecolor=_bar_color, linewidth=0.5
                        )
    except Exception:
        pass

    # --- Color-coded CAPE badge (upper-right of Skew-T) ---
    try:
        _sb_cape_val = float(params.get("sb_cape", 0 * units("J/kg")).magnitude)
        if _sb_cape_val > 0:
            if _sb_cape_val >= 4000:
                _badge_color = "#ff0000"
                _badge_label = "EXTREME"
            elif _sb_cape_val >= 3000:
                _badge_color = "#ff4400"
                _badge_label = "HIGH"
            elif _sb_cape_val >= 2000:
                _badge_color = "#ff8800"
                _badge_label = "MODERATE"
            elif _sb_cape_val >= 1000:
                _badge_color = "#ffcc00"
                _badge_label = "MARGINAL"
            else:
                _badge_color = "#88cc44"
                _badge_label = "LOW"
            skew.ax.text(
                0.98, 0.98,
                f"SBCAPE\n{_sb_cape_val:.0f}\n{_badge_label}",
                transform=skew.ax.transAxes,
                fontsize=12, fontweight="bold", fontfamily="monospace",
                color=_badge_color, ha="right", va="top",
                bbox=dict(boxstyle="round,pad=0.4", facecolor=BG,
                         edgecolor=_badge_color, alpha=0.9, linewidth=2),
                zorder=20,
                path_effects=[path_effects.withStroke(linewidth=1, foreground=BG)]
            )
    except Exception:
        pass

    # Find pressure levels where T crosses -12°C and -17°C for the band
    try:
        T_mag = T.to("degC").magnitude
        p_mag = p.magnitude
        # Find pressure at -12°C and -17°C by interpolation
        dgz_top_p, dgz_bot_p = None, None
        for i in range(len(T_mag) - 1):
            if dgz_bot_p is None and ((T_mag[i] >= -12 and T_mag[i+1] < -12) or
                                       (T_mag[i] <= -12 and T_mag[i+1] > -12)):
                frac = (T_mag[i] - (-12)) / (T_mag[i] - T_mag[i+1])
                dgz_bot_p = p_mag[i] + frac * (p_mag[i+1] - p_mag[i])
            if dgz_top_p is None and ((T_mag[i] >= -17 and T_mag[i+1] < -17) or
                                       (T_mag[i] <= -17 and T_mag[i+1] > -17)):
                frac = (T_mag[i] - (-17)) / (T_mag[i] - T_mag[i+1])
                dgz_top_p = p_mag[i] + frac * (p_mag[i+1] - p_mag[i])
        if dgz_bot_p is not None and dgz_top_p is not None:
            skew.ax.axhspan(dgz_top_p, dgz_bot_p, color="#00ccff", alpha=0.06, zorder=1)
            skew.ax.text(
                skew.ax.get_xlim()[0] + 3, (dgz_bot_p + dgz_top_p) / 2,
                "DGZ", color="#00ccff", fontsize=9, fontweight="bold",
                fontfamily="monospace", alpha=0.6, va="center",
                path_effects=[path_effects.withStroke(linewidth=2, foreground=BG)]
            )
    except Exception:
        pass
    
    # --- Hail Growth Zone (HGZ) shading: -10°C to -30°C ---
    try:
        hgz_top_p, hgz_bot_p = None, None
        for i in range(len(T_mag) - 1):
            if hgz_bot_p is None and ((T_mag[i] >= -10 and T_mag[i+1] < -10) or
                                       (T_mag[i] <= -10 and T_mag[i+1] > -10)):
                frac = (T_mag[i] - (-10)) / (T_mag[i] - T_mag[i+1])
                hgz_bot_p = p_mag[i] + frac * (p_mag[i+1] - p_mag[i])
            if hgz_top_p is None and ((T_mag[i] >= -30 and T_mag[i+1] < -30) or
                                       (T_mag[i] <= -30 and T_mag[i+1] > -30)):
                frac = (T_mag[i] - (-30)) / (T_mag[i] - T_mag[i+1])
                hgz_top_p = p_mag[i] + frac * (p_mag[i+1] - p_mag[i])
        if hgz_bot_p is not None and hgz_top_p is not None:
            skew.ax.axhspan(hgz_top_p, hgz_bot_p, color="#22ff88", alpha=0.04, zorder=1)
            skew.ax.text(
                skew.ax.get_xlim()[0] + 3, hgz_top_p,
                "HGZ", color="#22ff88", fontsize=9, fontweight="bold",
                fontfamily="monospace", alpha=0.5, va="top",
                path_effects=[path_effects.withStroke(linewidth=2, foreground=BG)]
            )
    except Exception:
        pass
    
    # --- Effective Inflow Layer (EIL) shading on Skew-T ---
    try:
        _eil_bot_p_skew = params.get("eil_bot_p")
        _eil_top_p_skew = params.get("eil_top_p")
        if _eil_bot_p_skew is not None and _eil_top_p_skew is not None:
            # Shade the effective inflow layer as a vertical band on the left edge
            skew.ax.axhspan(_eil_top_p_skew, _eil_bot_p_skew,
                           xmin=0, xmax=0.04,
                           color="#44ddaa", alpha=0.45, zorder=5)
            # Label
            _eil_mid_p = (_eil_bot_p_skew + _eil_top_p_skew) / 2
            _eil_bot_h_v = params.get("eil_bot_h", 0)
            _eil_top_h_v = params.get("eil_top_h", 0)
            skew.ax.text(
                skew.ax.get_xlim()[0] + 1, _eil_mid_p,
                f"EIL", color="#44ddaa",
                fontsize=8, fontweight="bold", fontfamily="monospace",
                alpha=0.9, va="center", ha="left",
                path_effects=[path_effects.withStroke(linewidth=2, foreground=BG)]
            )
    except Exception:
        pass

    # --- PBL (Planetary Boundary Layer) top marker ---
    # Estimate PBL using the ML (mixed-layer) LCL height as proxy
    try:
        ml_lcl_m = params.get("ml_lcl_m")
        if ml_lcl_m is not None and ml_lcl_m > 0:
            pbl_h_msl = ml_lcl_m + h[0].magnitude
            pbl_idx = np.argmin(np.abs(h.magnitude - pbl_h_msl))
            pbl_p = p.magnitude[pbl_idx]
            skew.ax.axhline(y=pbl_p, color="#bb88ff", linestyle="-.",
                            alpha=0.5, linewidth=1.0)
            skew.ax.text(
                skew.ax.get_xlim()[1] - 2, pbl_p,
                f"←PBL ({ml_lcl_m:.0f}m)", color="#bb88ff",
                fontsize=10, va="center", fontfamily="monospace",
                fontweight="bold",
                path_effects=[path_effects.withStroke(linewidth=3, foreground=BG)]
            )
    except Exception:
        pass
    
    # Set axis limits
    skew.ax.set_xlim(-100, 50)
    skew.ax.set_ylim(1050, 100)
    
    # Height labels — positioned inside the plot, offset from the y-axis
    h_agl = (h - h[0]).to("meter").magnitude
    for target_km in [1, 2, 3, 5, 7, 9, 12]:
        target_m = target_km * 1000
        idx = np.argmin(np.abs(h_agl - target_m))
        if abs(h_agl[idx] - target_m) < 500 and p.magnitude[idx] > 100:
            skew.ax.annotate(
                f"{target_km} km",
                xy=(0.01, p.magnitude[idx]),
                xycoords=("axes fraction", "data"),
                fontsize=10, color=ACCENT, fontfamily="monospace",
                fontweight="bold", va="center", ha="left",
                bbox=dict(boxstyle="round,pad=0.15", facecolor=BG,
                         edgecolor="none", alpha=0.85)
            )
    
    # Surface elevation label
    sfc_elev = h[0].magnitude
    skew.ax.annotate(
        f"SFC {sfc_elev:.0f}m",
        xy=(0.01, p.magnitude[0]),
        xycoords=("axes fraction", "data"),
        fontsize=10, color=FG_DIM, fontfamily="monospace",
        fontweight="bold", va="top", ha="left",
        bbox=dict(boxstyle="round,pad=0.15", facecolor=BG,
                 edgecolor="none", alpha=0.85)
    )
    
    # Legend — add OBS/VAD barb entries if VAD is present
    from matplotlib.lines import Line2D
    extra_handles, extra_labels = [], []
    if has_vad:
        extra_handles.append(Line2D([], [], color=FG, marker=r'$\rightarrow$',
                                    markersize=8, linestyle='None'))
        extra_labels.append('OBS BARBS')
        extra_handles.append(Line2D([], [], color=CLR_VAD, marker=r'$\rightarrow$',
                                    markersize=8, linestyle='None'))
        extra_labels.append('VAD BARBS')
    existing_handles, existing_labels = skew.ax.get_legend_handles_labels()
    all_handles = existing_handles + extra_handles
    all_labels = existing_labels + extra_labels
    legend = skew.ax.legend(
        all_handles, all_labels,
        loc="upper left", fontsize=10, facecolor=BG,
        edgecolor=BORDER, labelcolor=FG,
        framealpha=0.9, borderpad=0.5,
        title="SKEW-T LOG-P", title_fontproperties={"size": 12, "weight": "bold",
        "family": "monospace"}
    )
    legend.get_title().set_color(FG)
    # Color the OBS/VAD legend text to match their barb colors
    if has_vad:
        texts = legend.get_texts()
        for t in texts:
            if t.get_text() == 'OBS BARBS':
                t.set_color(FG)
            elif t.get_text() == 'VAD BARBS':
                t.set_color(CLR_VAD)
    
    # ════════════════════════════════════════════════════════════════
    # HODOGRAPH
    # ════════════════════════════════════════════════════════════════
    ax_hodo = fig.add_subplot(gs[0, 1])
    ax_hodo.set_facecolor(BG_PANEL)

    # --- Use the pre-computed 100-m interpolated winds from params ---
    hodo_u = params["u_interp"].to("knot").magnitude.copy()
    hodo_v = params["v_interp"].to("knot").magnitude.copy()
    hodo_z = params["h_interp"].to("meter").magnitude.copy()

    # --- Storm-relative hodograph mode ---
    # Subtract storm motion from all winds so SM sits at origin (0,0)
    sr_offset_u = 0.0
    sr_offset_v = 0.0
    if sr_hodograph and params.get("rm_u") is not None:
        sr_offset_u = params["rm_u"].to("knot").magnitude
        sr_offset_v = params["rm_v"].to("knot").magnitude
        hodo_u = hodo_u - sr_offset_u
        hodo_v = hodo_v - sr_offset_v

    # Cap at 9 km for hodograph bounds
    if hodo_z.max() > 9001:
        hodo_bound_idx = np.argmin(np.abs(hodo_z - 9000))
    else:
        hodo_bound_idx = len(hodo_u) - 1
    u_hodo = hodo_u[:hodo_bound_idx + 1]
    v_hodo = hodo_v[:hodo_bound_idx + 1]

    # --- Dynamic hodograph bounds ---
    x_min = u_hodo.min()
    y_min = v_hodo.min()
    x_max = u_hodo.max()
    y_max = v_hodo.max()

    y_Maxlimit = y_max + 30
    x_Maxlimit = x_max + 30
    y_Minlimit = y_min - 45
    x_Minlimit = x_min - 45

    # --- Create hodograph object ---
    hodo = Hodograph(ax_hodo, component_range=160.)
    try:
        hodo.ax.set_xlim(x_Minlimit, x_Maxlimit)
        hodo.ax.set_ylim(y_Minlimit, y_Maxlimit)
    except Exception:
        hodo.ax.set_xlim(-65, 65)
        hodo.ax.set_ylim(-65, 65)

    # Two grid layers — solid 20-kt and dashed 10-kt
    hodo.add_grid(increment=20, color=FG, linestyle='-', linewidth=1.5, alpha=0.2)
    hodo.add_grid(increment=10, color=FG, linewidth=1, linestyle='--', alpha=0.2)

    hodo.ax.set_facecolor(BG_PANEL)
    for spine in ax_hodo.spines.values():
        spine.set_color(BORDER)
    # Override MetPy's internal aspect='equal' so the hodograph fills the grid cell
    hodo.ax.set_aspect('auto')

    # Remove tick labels / axis labels
    hodo.ax.set_yticklabels([])
    hodo.ax.set_xticklabels([])
    hodo.ax.set_xticks([])
    hodo.ax.set_yticks([])
    hodo.ax.set_xlabel(' ')
    hodo.ax.set_ylabel(' ')

    # --- Velocity ring annotations ---
    for i in range(10, 130, 20):
        hodo.ax.annotate(str(i), (i, 0), xytext=(0, 2),
                         textcoords='offset pixels', clip_on=True,
                         fontsize=12, weight='bold', alpha=0.2, zorder=0, color=FG)
    for i in range(10, 130, 20):
        hodo.ax.annotate(str(i), (0, i), xytext=(0, 2),
                         textcoords='offset pixels', clip_on=True,
                         fontsize=12, weight='bold', alpha=0.2, zorder=0, color=FG)
    for i in range(10, 130, 20):
        hodo.ax.annotate(str(i), (-i, 0), xytext=(0, 2),
                         textcoords='offset pixels', clip_on=True,
                         fontsize=12, weight='bold', alpha=0.2, zorder=0, color=FG)
    for i in range(10, 130, 20):
        hodo.ax.annotate(str(i), (0, -i), xytext=(0, 2),
                         textcoords='offset pixels', clip_on=True,
                         fontsize=12, weight='bold', alpha=0.2, zorder=0, color=FG)

    # --- Height markers ---
    n_hodo = len(hodo_u)
    # 0.5 km marker
    idx_05 = 5  # index 5 = 500 m
    if idx_05 < n_hodo:
        hodo.ax.plot(hodo_u[idx_05], hodo_v[idx_05], '.',
                     color='white', markeredgecolor='black',
                     alpha=1, markersize=30, zorder=5, clip_on=True)
        hodo.ax.annotate('.5', (hodo_u[idx_05], hodo_v[idx_05]),
                         weight='bold', fontsize=13, color='black',
                         xytext=(0.02, -5), textcoords='offset pixels',
                         horizontalalignment='center', clip_on=True, zorder=6)

    # km markers at 1.5, 2.5, 3.5, ... (every other 500m level)
    hgt_lvl_indices = list(range(5, min(91, n_hodo), 5))  # every 500m: 5,10,15,...,90
    if len(hgt_lvl_indices) > 0:
        hgt_lvl_indices.pop(0)  # remove h05 (already plotted above)
    for lvl in hgt_lvl_indices[1::2]:  # every other level starting from index 1
        if lvl < n_hodo and lvl < 130:
            km_label = str(int(round(hodo_z[lvl] / 1000, 0)))
            hodo.ax.plot(hodo_u[lvl], hodo_v[lvl], '.',
                         color='white', markeredgecolor='black',
                         alpha=1, markersize=30, zorder=5, clip_on=True)
            hodo.ax.annotate(km_label, (hodo_u[lvl], hodo_v[lvl]),
                             weight='bold', fontsize=13, color='black',
                             xytext=(0.02, -5), textcoords='offset pixels',
                             horizontalalignment='center', clip_on=True, zorder=5.1)

    # --- Plot hodograph line in SounderPy colors ---
    hodo_color = ['purple', 'red', 'darkorange', 'gold', '#fff09f']
    n_full = len(hodo_u)
    # 0-1km: idx 0-10, 1-3km: 10-30, 3-6km: 30-60, 6-9km: 60-90, 9-12km: 90-120
    hodo.ax.plot(hodo_u[0:min(10+1, n_full)], hodo_v[0:min(10+1, n_full)],
                 color=hodo_color[0], linewidth=7, zorder=4, clip_on=True)
    hodo.ax.plot(hodo_u[10:min(30+1, n_full)], hodo_v[10:min(30+1, n_full)],
                 color=hodo_color[1], linewidth=7, zorder=4, clip_on=True)
    hodo.ax.plot(hodo_u[30:min(60+1, n_full)], hodo_v[30:min(60+1, n_full)],
                 color=hodo_color[2], linewidth=7, zorder=4, clip_on=True)
    hodo.ax.plot(hodo_u[60:min(90+1, n_full)], hodo_v[60:min(90+1, n_full)],
                 color=hodo_color[3], linewidth=7, zorder=4, clip_on=True)
    hodo.ax.plot(hodo_u[90:min(120+1, n_full)], hodo_v[90:min(120+1, n_full)],
                 color=hodo_color[4], linewidth=7, zorder=4, clip_on=True)

    # --- Storm motion annotations ---
    if params.get("rm_u") is not None:
        rm_u_kt = params["rm_u"].to("knot").magnitude - sr_offset_u
        rm_v_kt = params["rm_v"].to("knot").magnitude - sr_offset_v
        lm_u_kt = params["lm_u"].to("knot").magnitude - sr_offset_u
        lm_v_kt = params["lm_v"].to("knot").magnitude - sr_offset_v
        mw_u_kt = params["mw_u"].to("knot").magnitude - sr_offset_u
        mw_v_kt = params["mw_v"].to("knot").magnitude - sr_offset_v

        sm_u_kt = rm_u_kt
        sm_v_kt = rm_v_kt

        # RM / LM / MW text labels
        hodo.ax.text(rm_u_kt + 0.5, rm_v_kt - 0.5, 'RM',
                     weight='bold', ha='left', fontsize=15, zorder=7,
                     alpha=0.9, color=FG)
        hodo.ax.text(lm_u_kt + 0.5, lm_v_kt - 0.5, 'LM',
                     weight='bold', ha='left', fontsize=15, zorder=7,
                     alpha=0.9, color=FG)
        hodo.ax.text(mw_u_kt + 0.5, mw_v_kt - 0.5, 'MW',
                     weight='bold', ha='left', fontsize=15, zorder=7,
                     alpha=0.9, color=FG)

        # DTM (Deviant Tornado Motion)
        dtm_u_kt = mw_u_kt + (rm_v_kt - mw_v_kt)
        dtm_v_kt = mw_v_kt - (rm_u_kt - mw_u_kt)
        hodo.ax.text(dtm_u_kt, dtm_v_kt + 2, 'DTM',
                     weight='bold', fontsize=12, color='brown',
                     ha='center', zorder=7)
        hodo.ax.plot(dtm_u_kt, dtm_v_kt, marker='v', color='brown',
                     markersize=8, zorder=7, alpha=0.8, ls='')

        # SM arrow / origin marker
        if sr_hodograph:
            # In SR mode, SM is at origin — mark it with a crosshair
            hodo.ax.plot(0, 0, '+', color=ACCENT, markersize=15, markeredgewidth=2.5,
                         zorder=7, alpha=0.9)
            hodo.ax.text(2, -3, 'SM', weight='bold', fontsize=12, color=ACCENT,
                         ha='left', zorder=7, alpha=0.9)
        else:
            # Ground-relative: arrow from origin to storm motion
            hodo.ax.arrow(0, 0, sm_u_kt - 0.3, sm_v_kt - 0.3,
                          linewidth=3, color=FG, alpha=0.2,
                          label='SM Vector', length_includes_head=True,
                          head_width=0.6)

        # --- Effective inflow layer SRH fill (using actual EIL bounds) ---
        _eil_bot_h = params.get("eil_bot_h")
        _eil_top_h = params.get("eil_top_h")
        if _eil_bot_h is not None and _eil_top_h is not None:
            eil_bot_idx = max(0, min(int(round(_eil_bot_h / 100)), n_full - 1))
            eil_top_idx = max(0, min(int(round(_eil_top_h / 100)), n_full - 1))
        else:
            # Fallback to 0-3 km if no effective layer
            eil_bot_idx = 0
            eil_top_idx = min(30, n_full - 1)
        _eil_label = (f'EIL {int(_eil_bot_h)}-{int(_eil_top_h)}m SRH'
                      if _eil_bot_h is not None else '0-3 SRH')
        hodo.ax.plot(
            (sm_u_kt, hodo_u[eil_bot_idx]),
            (sm_v_kt, hodo_v[eil_bot_idx]),
            linestyle='-', linewidth=2.3, alpha=0.5, zorder=3,
            color='lightblue', label='Effective Inflow Layer')
        hodo.ax.plot(
            (sm_u_kt, hodo_u[eil_top_idx]),
            (sm_v_kt, hodo_v[eil_top_idx]),
            linestyle='-', linewidth=2.3, alpha=0.5, zorder=3,
            color='lightblue')
        hodo.ax.fill(
            np.append(hodo_u[eil_bot_idx:eil_top_idx + 1], sm_u_kt),
            np.append(hodo_v[eil_bot_idx:eil_top_idx + 1], sm_v_kt),
            'lightblue', alpha=0.3, zorder=2, label=_eil_label)

        # --- Corfidi MCS motion vectors (Corfidi 2003) ---
        if params.get("corfidi_up_u") is not None and params.get("corfidi_dn_u") is not None:
            cup_u = params["corfidi_up_u"]
            cup_v = params["corfidi_up_v"]
            cdn_u = params["corfidi_dn_u"]
            cdn_v = params["corfidi_dn_v"]
            hodo.ax.plot(cup_u, cup_v, 's', markersize=9, color='orange',
                         markeredgecolor='white', markeredgewidth=0.8,
                         zorder=15, alpha=0.85, clip_on=True)
            hodo.ax.text(cup_u, cup_v + 3.5, 'CU', weight='bold', fontsize=10,
                         color='orange', ha='center', alpha=0.8, clip_on=True)
            hodo.ax.plot(cdn_u, cdn_v, 's', markersize=9, color='#ff4444',
                         markeredgecolor='white', markeredgewidth=0.8,
                         zorder=15, alpha=0.85, clip_on=True)
            hodo.ax.text(cdn_u, cdn_v + 3.5, 'CD', weight='bold', fontsize=10,
                         color='#ff4444', ha='center', alpha=0.8, clip_on=True)
        else:
            # Fallback to simplified MCS markers
            bwd_u = (hodo_u[min(60, n_full-1)] - hodo_u[0])
            bwd_v = (hodo_v[min(60, n_full-1)] - hodo_v[0])
            bwd_mag = np.sqrt(bwd_u**2 + bwd_v**2)
            if bwd_mag > 0.1:
                us_u = mw_u_kt + 7.5 * 1.94384 * (-bwd_u / bwd_mag)
                us_v = mw_v_kt + 7.5 * 1.94384 * (-bwd_v / bwd_mag)
                ds_u = mw_u_kt + 7.5 * 1.94384 * (bwd_u / bwd_mag)
                ds_v = mw_v_kt + 7.5 * 1.94384 * (bwd_v / bwd_mag)
            else:
                us_u, us_v = mw_u_kt, mw_v_kt
                ds_u, ds_v = mw_u_kt, mw_v_kt
            hodo.ax.text(us_u, us_v, 'UP', weight='bold', fontsize=12,
                         color='orange', ha='center', alpha=0.5, clip_on=True)
            hodo.ax.text(ds_u, ds_v, 'DN', weight='bold', fontsize=12,
                         color='orange', ha='center', alpha=0.5, clip_on=True)

    # --- Storm motion info text box ---
    sm_lines = []
    if params.get("rm_u") is not None:
        rm_u_kt = params["rm_u"].to("knot").magnitude
        rm_v_kt = params["rm_v"].to("knot").magnitude
        lm_u_kt = params["lm_u"].to("knot").magnitude
        lm_v_kt = params["lm_v"].to("knot").magnitude
        mw_u_kt = params["mw_u"].to("knot").magnitude
        mw_v_kt = params["mw_v"].to("knot").magnitude

        def _wind_to_dir_str(u_val, v_val):
            deg = int(np.degrees(np.arctan2(-u_val, -v_val)) % 360)
            dirs = ['N','NNE','NE','ENE','E','ESE','SE','SSE',
                    'S','SSW','SW','WSW','W','WNW','NW','NNW']
            idx = int((deg + 11.25) / 22.5) % 16
            return dirs[idx]
        def _wind_spd(u_val, v_val):
            return int(np.sqrt(u_val**2 + v_val**2))

        dtm_u_kt = mw_u_kt + (rm_v_kt - mw_v_kt)
        dtm_v_kt = mw_v_kt - (rm_u_kt - mw_u_kt)

        bwd_u = (hodo_u[min(60, n_full-1)] - hodo_u[0])
        bwd_v = (hodo_v[min(60, n_full-1)] - hodo_v[0])
        bwd_mag = np.sqrt(bwd_u**2 + bwd_v**2)
        if bwd_mag > 0.1:
            us_u = mw_u_kt + 7.5 * 1.94384 * (-bwd_u / bwd_mag)
            us_v = mw_v_kt + 7.5 * 1.94384 * (-bwd_v / bwd_mag)
            ds_u = mw_u_kt + 7.5 * 1.94384 * (bwd_u / bwd_mag)
            ds_v = mw_v_kt + 7.5 * 1.94384 * (bwd_v / bwd_mag)
        else:
            us_u, us_v = mw_u_kt, mw_v_kt
            ds_u, ds_v = mw_u_kt, mw_v_kt

        # --- Critical Angle ---
        # Angle between the 0-500m shear vector and the storm-relative inflow vector
        try:
            # 0-500m shear vector (surface to 500m AGL)
            shr_u = hodo_u[min(5, n_full-1)] - hodo_u[0]
            shr_v = hodo_v[min(5, n_full-1)] - hodo_v[0]
            # Storm-relative inflow vector (surface wind - storm motion)
            sri_u = hodo_u[0] - sm_u_kt
            sri_v = hodo_v[0] - sm_v_kt
            shr_mag = np.sqrt(shr_u**2 + shr_v**2)
            sri_mag = np.sqrt(sri_u**2 + sri_v**2)
            if shr_mag > 0.5 and sri_mag > 0.5:
                cos_angle = (shr_u * sri_u + shr_v * sri_v) / (shr_mag * sri_mag)
                cos_angle = np.clip(cos_angle, -1, 1)
                crit_angle = np.degrees(np.arccos(cos_angle))
                crit_angle_str = f"CRIT∠: {crit_angle:.0f}°"
            else:
                crit_angle_str = "CRIT∠: N/A"
        except Exception:
            crit_angle_str = "CRIT∠: N/A"
        
        # For the info text box, use original (ground-relative) values
        _info_rm_u = params["rm_u"].to("knot").magnitude
        _info_rm_v = params["rm_v"].to("knot").magnitude
        _info_lm_u = params["lm_u"].to("knot").magnitude
        _info_lm_v = params["lm_v"].to("knot").magnitude
        _info_mw_u = params["mw_u"].to("knot").magnitude
        _info_mw_v = params["mw_v"].to("knot").magnitude
        _info_dtm_u = _info_mw_u + (_info_rm_v - _info_mw_v)
        _info_dtm_v = _info_mw_v - (_info_rm_u - _info_mw_u)

        # US/DS use BWD from hodo winds (which is the same in both frames)
        bwd_u = (hodo_u[min(60, n_full-1)] - hodo_u[0])
        bwd_v = (hodo_v[min(60, n_full-1)] - hodo_v[0])
        bwd_mag = np.sqrt(bwd_u**2 + bwd_v**2)
        if bwd_mag > 0.1:
            _info_us_u = _info_mw_u + 7.5 * 1.94384 * (-bwd_u / bwd_mag)
            _info_us_v = _info_mw_v + 7.5 * 1.94384 * (-bwd_v / bwd_mag)
            _info_ds_u = _info_mw_u + 7.5 * 1.94384 * (bwd_u / bwd_mag)
            _info_ds_v = _info_mw_v + 7.5 * 1.94384 * (bwd_v / bwd_mag)
        else:
            _info_us_u, _info_us_v = _info_mw_u, _info_mw_v
            _info_ds_u, _info_ds_v = _info_mw_u, _info_mw_v

        sm_lines = [
            f"{'STORM-RELATIVE | ' if sr_hodograph else ''}SM: RIGHT MOVING",
            f"RM: {_wind_to_dir_str(_info_rm_u, _info_rm_v)} @ {_wind_spd(_info_rm_u, _info_rm_v)} kts",
            f"LM: {_wind_to_dir_str(_info_lm_u, _info_lm_v)} @ {_wind_spd(_info_lm_u, _info_lm_v)} kts",
            f"MW: {_wind_to_dir_str(_info_mw_u, _info_mw_v)} @ {_wind_spd(_info_mw_u, _info_mw_v)} kts",
            f"DTM: {_wind_to_dir_str(_info_dtm_u, _info_dtm_v)} @ {_wind_spd(_info_dtm_u, _info_dtm_v)} kts",
            f"US: {_wind_to_dir_str(_info_us_u, _info_us_v)} @ {_wind_spd(_info_us_u, _info_us_v)} kts",
            f"DS: {_wind_to_dir_str(_info_ds_u, _info_ds_v)} @ {_wind_spd(_info_ds_u, _info_ds_v)} kts",
            crit_angle_str,
        ]

    sm_text = "\n".join(sm_lines) if sm_lines else ""
    if sm_text:
        ax_hodo.text(
            0.02, 0.98, sm_text,
            transform=ax_hodo.transAxes, fontsize=11,
            color=FG, fontfamily="monospace", fontweight="bold",
            va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.3", facecolor=BG,
                     edgecolor=BORDER, alpha=0.92)
        )

    # --- VAD Wind Profile overlay ---
    if vad_data and isinstance(vad_data, dict) and vad_data.get("winds"):
        vad_winds = vad_data["winds"]
        vad_u = np.array([w["u_kt"] for w in vad_winds]) - sr_offset_u
        vad_v = np.array([w["v_kt"] for w in vad_winds]) - sr_offset_v
        vad_alt_m = np.array([w["alt_m"] for w in vad_winds])
        vad_alt_agl = vad_alt_m  # VWP altitudes are already AGL (above radar level ≈ AGL)

        VAD_COLOR = CLR_VAD
        # Draw VAD hodograph line
        hodo.ax.plot(vad_u, vad_v, color=VAD_COLOR, linewidth=2.5,
                     linestyle='-', alpha=0.85, zorder=8, clip_on=True)
        # Draw dots at each VAD level
        hodo.ax.scatter(vad_u, vad_v, c=VAD_COLOR, s=30, zorder=9,
                       edgecolors='black', linewidths=0.5, alpha=0.9, clip_on=True)

        # Height labels at select VAD levels (every ~1 km)
        labeled_km = set()
        for i, alt in enumerate(vad_alt_agl):
            km_val = round(alt / 1000.0)
            if km_val >= 1 and km_val not in labeled_km and alt >= 800:
                labeled_km.add(km_val)
                hodo.ax.annotate(
                    f"{km_val}k", (vad_u[i], vad_v[i]),
                    xytext=(5, 5), textcoords='offset pixels',
                    fontsize=7, color=VAD_COLOR, fontweight='bold',
                    alpha=0.8, clip_on=True,
                    path_effects=[path_effects.withStroke(linewidth=2, foreground=BG)]
                )

        # VAD label
        vad_label = f"VAD: {vad_data.get('radar', '???')}"
        vad_meta = vad_data.get("meta", {})
        if vad_meta.get("time"):
            vad_label += f" | {vad_meta['time']}"
        hodo.ax.text(
            0.98, 0.02, vad_label,
            transform=ax_hodo.transAxes, fontsize=8,
            color=VAD_COLOR, fontfamily="monospace", fontweight="bold",
            va="bottom", ha="right", alpha=0.9,
            bbox=dict(boxstyle="round,pad=0.2", facecolor=BG,
                     edgecolor=VAD_COLOR, alpha=0.8, linewidth=0.8)
        )

    # --- Boundary orientation line ---
    if boundary_orientation is not None:
        try:
            bdry_deg = float(boundary_orientation) % 360
            # The orientation defines the direction the boundary *runs* along
            # (like a front oriented NE-SW = 045°).  We draw a line in that
            # direction through the origin (or center of hodograph).
            bdry_rad = np.radians(bdry_deg)
            # Meteorological convention: 0° = N, 90° = E
            # Convert to math coordinates: x = sin(θ), y = cos(θ)
            dx = np.sin(bdry_rad)
            dy = np.cos(bdry_rad)
            # Line length: use current axis limits
            xlim = hodo.ax.get_xlim()
            ylim = hodo.ax.get_ylim()
            span = max(abs(xlim[1] - xlim[0]), abs(ylim[1] - ylim[0]))
            L = span * 0.8

            # Center the line on the surface wind
            cx, cy = hodo_u[0], hodo_v[0]

            bdry_color = "#ff44ff" if not colorblind else "#ff8800"
            hodo.ax.plot(
                [cx - L * dx, cx + L * dx],
                [cy - L * dy, cy + L * dy],
                color=bdry_color, linewidth=2.5, linestyle="--",
                alpha=0.85, zorder=10, clip_on=True,
            )
            # Small label at one end
            label_x = cx + L * 0.65 * dx
            label_y = cy + L * 0.65 * dy
            hodo.ax.text(
                label_x, label_y,
                f"BDRY {int(bdry_deg)}°",
                color=bdry_color, fontsize=10, fontweight="bold",
                fontfamily="monospace", ha="center", va="bottom",
                alpha=0.9, zorder=11, clip_on=True,
                path_effects=[path_effects.withStroke(linewidth=3, foreground=BG)],
            )
        except Exception:
            pass
    
    # ════════════════════════════════════════════════════════════════
    # STORM-RELATIVE WIND PROFILE (right side panel)
    # ════════════════════════════════════════════════════════════════
    ax_srw = fig.add_subplot(gs[0, 2])
    ax_srw.set_facecolor(BG_PANEL)
    
    if params.get("rm_u") is not None:
        # Compute storm-relative wind
        rm_u_ms = params["rm_u"].to("m/s").magnitude
        rm_v_ms = params["rm_v"].to("m/s").magnitude
        sru = u.to("m/s").magnitude - rm_u_ms
        srv = v.to("m/s").magnitude - rm_v_ms
        sr_spd = np.sqrt(sru**2 + srv**2) * 1.94384  # to knots
        
        # Plot SR wind speed vs height AGL
        h_km = h_agl / 1000.0
        ax_srw.plot(sr_spd, h_km, color="#ff8800", linewidth=2.5)
        ax_srw.fill_betweenx(h_km, 0, sr_spd, alpha=0.2, color="#ff8800")
        ax_srw.set_xlim(0, max(sr_spd.max() * 1.1, 20))
        ax_srw.set_ylim(0, min(h_km.max(), 12))
        ax_srw.set_xlabel("SRW (kt)", color=FG_DIM, fontsize=10,
                         fontfamily="monospace", fontweight="bold")
        ax_srw.set_ylabel("Height AGL (km)", color=FG_DIM, fontsize=10,
                         fontfamily="monospace", fontweight="bold")
        ax_srw.set_title("STORM-REL\nWIND", color=FG, fontsize=10,
                        fontfamily="monospace", fontweight="bold", pad=3)
    else:
        ax_srw.text(0.5, 0.5, "N/A", transform=ax_srw.transAxes,
                   color=FG_FAINT, ha="center", fontsize=12)
    
    ax_srw.tick_params(colors=FG_DIM, labelsize=9, width=1.2)
    for spine in ax_srw.spines.values():
        spine.set_color(BORDER)
    ax_srw.grid(True, alpha=0.25, color=GRID_CLR)
    
    # ════════════════════════════════════════════════════════════════
    # STREAMWISENESS PROFILE
    # ════════════════════════════════════════════════════════════════
    ax_sw = fig.add_subplot(gs[0, 3])
    ax_sw.set_facecolor(BG_PANEL)
    
    if params.get("streamwiseness") is not None:
        sw_vals = params["streamwiseness"] * 100  # percent
        sw_h = params["streamwiseness_height"]
        
        # Color by streamwise (positive=cyclonic) vs anticyclonic
        sw_signed = params["streamwiseness_signed"]
        
        # Plot streamwiseness percentage vs height
        ax_sw.plot(sw_vals, sw_h, color="#44ddaa", linewidth=2.5, zorder=5)
        ax_sw.fill_betweenx(sw_h, 0, sw_vals,
                            where=(sw_signed >= 0),
                            alpha=0.2, color="#ff3333", interpolate=True,
                            label="Cyclonic")
        ax_sw.fill_betweenx(sw_h, 0, sw_vals,
                            where=(sw_signed < 0),
                            alpha=0.2, color="#4488ff", interpolate=True,
                            label="Anticyclonic")
        
        # Mark key depth levels
        for depth_km, label, color in [
            (0.5, "500m", "#aaaaaa"),
            (1.0, "1km", "#ff8800"),
            (3.0, "3km", "#ffcc00"),
        ]:
            idx_depth = np.argmin(np.abs(sw_h - depth_km))
            if idx_depth < len(sw_vals):
                val = sw_vals[idx_depth]
                ax_sw.axhline(y=depth_km, color=color, linestyle="--",
                              alpha=0.4, linewidth=0.8)
                ax_sw.plot(val, depth_km, "o", color=color, markersize=6,
                           markeredgecolor="white", markeredgewidth=0.5, zorder=6)
                ax_sw.annotate(f"{val:.0f}%", (val, depth_km),
                              xytext=(5, 3), textcoords="offset points",
                              fontsize=9, color=color, fontweight="bold",
                              fontfamily="monospace")
        
        # Axis limits
        ax_sw.set_xlim(0, 105)
        ax_sw.set_ylim(0, min(sw_h.max(), 6))
        ax_sw.set_xlabel("Streamwiseness (%)", color=FG_DIM, fontsize=10,
                         fontfamily="monospace", fontweight="bold")
        ax_sw.set_ylabel("Height AGL (km)", color=FG_DIM, fontsize=10,
                         fontfamily="monospace", fontweight="bold")
        ax_sw.set_title("STREAM-\nWISENESS", color=FG, fontsize=10,
                        fontfamily="monospace", fontweight="bold", pad=3)
        
        # Legend
        leg = ax_sw.legend(loc="lower right", fontsize=7, facecolor=BG,
                           edgecolor=BORDER, labelcolor=FG, framealpha=0.9)
    else:
        ax_sw.text(0.5, 0.5, "N/A", transform=ax_sw.transAxes,
                   color=FG_FAINT, ha="center", fontsize=12)
        ax_sw.set_title("STREAM-\nWISENESS", color=FG, fontsize=10,
                        fontfamily="monospace", fontweight="bold", pad=3)
    
    ax_sw.tick_params(colors=FG_DIM, labelsize=9, width=1.2)
    for spine in ax_sw.spines.values():
        spine.set_color(BORDER)
    ax_sw.grid(True, alpha=0.25, color=GRID_CLR)
    
    # ════════════════════════════════════════════════════════════════
    # THETA / THETA-E PROFILE
    # ════════════════════════════════════════════════════════════════
    ax_th = fig.add_subplot(gs[0, 4])
    ax_th.set_facecolor(BG_PANEL)
    
    try:
        # Compute potential temperature (θ) and equivalent potential temperature (θe)
        theta = mpcalc.potential_temperature(p, T).to("K").magnitude
        theta_e = mpcalc.equivalent_potential_temperature(p, T, Td).to("K").magnitude
        
        h_km = h_agl / 1000.0
        
        # Plot θ (orange) and θe (cyan)
        ax_th.plot(theta, h_km, color="#ff8800", linewidth=2.2, label="θ", zorder=5)
        ax_th.plot(theta_e, h_km, color="#00ccff", linewidth=2.2, label="θe", zorder=5)
        
        # Fill between θ and θe (shows moisture)
        ax_th.fill_betweenx(h_km, theta, theta_e, alpha=0.12, color="#00ccff")
        
        # Mark key heights
        for depth_km, lbl, clr in [
            (1.0, "1km", "#ff8800"),
            (3.0, "3km", "#ffcc00"),
            (6.0, "6km", "#44ddaa"),
        ]:
            idx = np.argmin(np.abs(h_km - depth_km))
            if idx < len(theta) and depth_km <= h_km.max():
                ax_th.axhline(y=depth_km, color=clr, linestyle="--",
                              alpha=0.35, linewidth=0.7)
        
        # Mark lapse rate stability regions
        # dθe/dz < 0 → conditionally unstable (red shading)
        if len(theta_e) > 2:
            dthe_dz = np.gradient(theta_e, h_km)
            unstable = dthe_dz < 0
            ax_th.fill_betweenx(h_km, ax_th.get_xlim()[0] if ax_th.get_xlim()[0] > 0 else theta.min() - 5,
                                theta.min() - 5, where=unstable,
                                alpha=0.0)  # placeholder, actual shading below
        
        # Auto x-limits based on data range
        th_min = min(theta.min(), theta_e.min()) - 5
        th_max = max(theta.max(), theta_e.max()) + 5
        ax_th.set_xlim(th_min, th_max)
        ax_th.set_ylim(0, min(h_km.max(), 12))
        
        ax_th.set_xlabel("θ / θe (K)", color=FG_DIM, fontsize=10,
                         fontfamily="monospace", fontweight="bold")
        ax_th.set_ylabel("Height AGL (km)", color=FG_DIM, fontsize=10,
                         fontfamily="monospace", fontweight="bold")
        ax_th.set_title("θ / θe\nPROFILE", color=FG, fontsize=10,
                        fontfamily="monospace", fontweight="bold", pad=3)
        
        # Legend
        leg_th = ax_th.legend(loc="lower right", fontsize=8, facecolor=BG,
                              edgecolor=BORDER, labelcolor=FG, framealpha=0.9)
    except Exception as ex:
        print(f"[Theta panel] Error: {ex}")
        ax_th.text(0.5, 0.5, "N/A", transform=ax_th.transAxes,
                   color=FG_FAINT, ha="center", fontsize=12)
        ax_th.set_title("θ / θe\nPROFILE", color=FG, fontsize=10,
                        fontfamily="monospace", fontweight="bold", pad=3)
    
    ax_th.tick_params(colors=FG_DIM, labelsize=9, width=1.2)
    for spine in ax_th.spines.values():
        spine.set_color(BORDER)
    ax_th.grid(True, alpha=0.25, color=GRID_CLR)
    
    # ════════════════════════════════════════════════════════════════
    # MINI MAP INSET (station location)
    # ════════════════════════════════════════════════════════════════
    # Detailed CONUS outline — traced clockwise from Pacific NW
    _conus_lon = [
        # Pacific NW coast (WA)
        -124.7, -124.6, -124.1, -123.2, -122.9, -122.8, -123.0, -123.5,
        -124.1, -124.6, -124.4, -124.2, -124.0, -123.9, -124.4,
        # OR coast
        -124.6, -124.5, -124.2, -124.3, -124.6, -124.5, -124.2, -124.4,
        -124.3, -124.1,
        # CA coast
        -124.2, -124.0, -123.7, -122.4, -122.0, -121.8, -122.0, -122.4,
        -122.5, -121.9, -121.3, -120.9, -120.6, -120.6, -120.2, -119.5,
        -118.5, -117.9, -117.6, -117.2, -117.1, -117.1,
        # US-Mexico border
        -114.7, -111.1, -109.0, -108.2, -106.6, -104.9, -103.3, -101.4,
        -100.0, -99.2, -97.8, -97.1,
        # TX Gulf coast
        -97.2, -97.0, -96.8, -96.4, -95.5, -94.9, -94.7, -93.8, -93.5,
        -93.2, -93.0, -91.6, -90.6,
        # LA coast
        -90.1, -89.7, -89.5, -89.2, -89.4, -89.0, -88.8, -88.6, -88.9,
        -88.5, -88.4,
        # MS-AL-FL Gulf coast
        -88.3, -87.6, -87.2, -86.5, -85.8, -85.5, -85.0, -84.0, -83.5,
        -82.8, -82.2, -81.8, -81.1, -80.4, -80.0, -80.1,
        # FL Atlantic coast
        -80.4, -80.6, -80.3, -80.0, -80.1, -80.5, -81.0, -81.3, -81.3,
        -81.0, -80.7, -80.5,
        # GA-SC-NC coast
        -80.8, -81.1, -80.8, -79.9, -79.1, -78.5, -77.9, -77.7, -76.5,
        -75.8, -75.5, -75.5, -76.0,
        # VA-MD-DE-NJ coast
        -75.6, -75.7, -76.0, -75.5, -74.9, -74.5, -74.0, -73.9,
        # NY-CT-RI-MA coast
        -74.0, -73.6, -72.8, -72.0, -71.2, -70.2, -70.0, -70.8, -71.4,
        -71.0, -70.5, -70.0, -69.9,
        # ME coast
        -69.8, -69.0, -68.5, -67.8, -67.0, -67.0,
        # Northern border — ME to MN (US-Canada)
        -67.1, -67.8, -69.0, -69.2, -71.1, -71.5, -73.4, -74.7, -75.0,
        -76.8, -79.0, -79.5, -82.4, -83.5, -84.1, -84.8, -85.0, -88.4,
        -89.6, -90.0, -89.5, -90.8, -92.0, -92.2, -94.6, -94.6, -95.1,
        # ND-MT-WA border (49°N)
        -95.2, -97.0, -99.0, -100.0, -102.0, -104.0, -106.0, -109.0,
        -111.0, -113.0, -116.0, -117.0, -120.0, -122.0, -122.8, -123.0,
        # WA coast back to start
        -123.2, -124.1, -124.7,
    ]
    _conus_lat = [
        # Pacific NW coast (WA)
        40.3, 42.0, 43.4, 46.2, 46.7, 47.1, 47.5, 48.0,
        48.3, 48.1, 47.5, 47.0, 46.6, 46.2, 44.7,
        # OR coast
        44.0, 43.5, 43.3, 42.8, 42.4, 42.0, 41.8, 41.0,
        40.5, 40.0,
        # CA coast
        39.5, 39.0, 38.5, 37.8, 37.5, 37.0, 36.8, 36.5,
        36.2, 36.0, 35.7, 35.3, 35.0, 34.8, 34.1, 33.4,
        33.0, 32.8, 32.6, 32.7, 33.0, 32.5,
        # US-Mexico border
        32.5, 31.3, 31.3, 31.4, 31.9, 30.6, 29.1, 29.8,
        28.2, 26.4, 26.0, 25.8,
        # TX Gulf coast
        26.3, 27.5, 28.2, 28.6, 28.7, 29.2, 29.6, 29.8, 29.8,
        29.5, 29.7, 30.0, 29.4,
        # LA coast
        29.1, 29.1, 29.3, 29.1, 28.9, 29.1, 28.9, 29.2, 29.3,
        30.2, 30.4,
        # MS-AL-FL Gulf coast
        30.4, 30.3, 30.3, 30.4, 30.2, 29.9, 29.5, 29.9, 30.0,
        29.9, 29.5, 29.1, 28.5, 27.0, 26.0, 25.3,
        # FL Atlantic coast
        25.5, 25.8, 26.5, 27.2, 28.0, 28.5, 29.0, 29.8, 30.5,
        30.8, 31.2, 32.0,
        # GA-SC-NC coast
        32.1, 32.1, 32.6, 33.1, 33.2, 33.7, 33.9, 34.3, 34.7,
        35.2, 35.8, 36.5, 36.9,
        # VA-MD-DE-NJ coast
        37.2, 37.5, 38.0, 38.5, 39.0, 39.3, 39.5, 40.5,
        # NY-CT-RI-MA coast
        40.7, 40.6, 40.8, 41.0, 41.5, 41.7, 42.0, 42.3, 42.0,
        41.5, 41.5, 41.8, 43.4,
        # ME coast
        43.9, 44.2, 44.4, 44.5, 44.9, 47.3,
        # Northern border — ME to MN
        47.4, 47.1, 47.2, 47.5, 45.0, 45.0, 45.0, 45.0, 44.8,
        43.6, 43.5, 43.2, 46.0, 46.1, 46.4, 46.0, 46.5, 48.2,
        48.0, 48.1, 47.1, 47.0, 46.7, 46.1, 46.8, 49.0, 49.0,
        # ND-MT-WA border (49°N)
        49.0, 49.0, 49.0, 49.0, 49.0, 49.0, 49.0, 49.0,
        49.0, 49.0, 49.0, 49.0, 49.0, 48.8, 48.7, 48.4,
        # WA coast back to start
        47.0, 44.0, 40.3,
    ]
    
    # ════════════════════════════════════════════════════════════════
    # PARAMETER TABLES (bottom)
    # ════════════════════════════════════════════════════════════════
    ax_params = fig.add_subplot(gs[1, 0])
    ax_params.set_facecolor(BG)
    ax_params.axis("off")
    
    # ── THERMODYNAMIC TABLE ──
    header = f"{'':3s}{'':8s}  {'CAPE':>10s}  {'CIN':>10s}  {'LCL':>8s}"
    
    sb_cape = fv(params.get("sb_cape"), "J/kg")
    sb_cin = fv(params.get("sb_cin"), "J/kg")
    sb_lcl = f"{params['sb_lcl_m']:.0f} m" if params.get('sb_lcl_m') is not None else fv(params.get("sb_lcl_p"), "hPa")
    
    mu_cape = fv(params.get("mu_cape"), "J/kg")
    mu_cin = fv(params.get("mu_cin"), "J/kg")
    mu_lcl = f"{params['mu_lcl_m']:.0f} m" if params.get('mu_lcl_m') is not None else fv(params.get("mu_lcl_p"), "hPa")
    
    ml_cape = fv(params.get("ml_cape"), "J/kg")
    ml_cin = fv(params.get("ml_cin"), "J/kg")
    ml_lcl = f"{params['ml_lcl_m']:.0f} m" if params.get('ml_lcl_m') is not None else fv(params.get("ml_lcl_p"), "hPa")
    
    dcape_v = fv(params.get("dcape"), "J/kg") if params.get("dcape") else "---"
    dcin_v = f"{params.get('dcin', 0)} J/kg" if params.get("dcin") else "0 J/kg"
    ncape_v = f"{params.get('ncape', 0):.3f}" if params.get("ncape") else "0"
    cape3_v = f"{params.get('cape_3km', 0)} J/kg"
    cape6_v = f"{params.get('cape_6km', 0)} J/kg"
    lr03 = f"{params['lr_03']:.1f}" if params.get('lr_03') is not None else "---"
    lr36 = f"{params['lr_36']:.1f}" if params.get('lr_36') is not None else "---"
    
    pwat = fv(params.get("pwat"), "mm", 1) if params.get("pwat") else "---"
    frz = f"{params['frz_level']:.0f}m" if params.get("frz_level") else "---"
    wbo_v = f"{params['wbo']:.0f}m" if params.get('wbo') is not None else "---"
    wcd_v = f"{params['wcd']:.0f}m" if params.get('wcd') is not None else "---"
    rh_01 = f"{params.get('rh_0_1km', 0):.0f}%" if params.get('rh_0_1km') else "---"
    rh_13 = f"{params.get('rh_1_3km', 0):.0f}%" if params.get('rh_1_3km') else "---"
    rh_36 = f"{params.get('rh_3_6km', 0):.0f}%" if params.get('rh_3_6km') else "---"
    stp_v = fv(params.get("stp"), "", 1)
    stp_eff_v = fv(params.get("stp_eff"), "", 1)
    ml_lfc_v = f"{params['ml_lfc_m']:.0f}m" if params.get('ml_lfc_m') is not None else "---"
    ml_el_v = f"{params['ml_el_m']:.0f}m" if params.get('ml_el_m') is not None else "---"
    ecape_v = f"{params.get('ecape', 0)} J/kg" if params.get('ecape') else "0 J/kg"
    precip_type_v = params.get("precip_type", "N/A")
    fosberg_v = f"{params['fosberg_fwi']:.0f}" if params.get('fosberg_fwi') is not None else "---"
    haines_v = f"{params['haines']}" if params.get('haines') is not None else "---"
    hdw_v = f"{params['hdw']:.0f}" if params.get('hdw') is not None else "---"
    
    # Rows: (text, color, y_position) — tightly spaced
    thermo_rows = [
        ("THERMODYNAMIC", ACCENT, 0.97),
        (header, FG_FAINT, 0.91),
        (f"   {'SB:':8s}  {sb_cape:>10s}  {sb_cin:>10s}  {sb_lcl:>8s}", CLR_SB_PARCEL, 0.85),
        (f"   {'MU:':8s}  {mu_cape:>10s}  {mu_cin:>10s}  {mu_lcl:>8s}", CLR_MU_PARCEL, 0.79),
        (f"   {'ML:':8s}  {ml_cape:>10s}  {ml_cin:>10s}  {ml_lcl:>8s}", CLR_ML_PARCEL, 0.73),
        (f"   ML LFC: {ml_lfc_v}  |  ML EL: {ml_el_v}  |  WCD: {wcd_v}", FG_DIM, 0.67),
        (f"   DCAPE: {dcape_v}  |  ECAPE: {ecape_v}", FG_DIM, 0.61),
        (f"   3CAPE: {cape3_v}  |  6CAPE: {cape6_v}  |  NCAPE: {ncape_v}", FG_DIM, 0.55),
        (f"   DCIN: {dcin_v}", FG_DIM, 0.49),
        (f"   \u03930-3: {lr03} \u00b0C/km   \u03933-6: {lr36} \u00b0C/km", FG_DIM, 0.43),
        (f"   PWAT: {pwat}  |  FRZ: {frz}  |  WBO: {wbo_v}", FG_DIM, 0.37),
        (f"   RH  0-1km: {rh_01}  1-3km: {rh_13}  3-6km: {rh_36}", FG_DIM, 0.31),
        (f"   Precip Type: {precip_type_v}", "#22d3ee", 0.25),
        (f"   STP: {stp_v}  STPeff: {stp_eff_v}  |  SCP: {fv(params.get('scp'), '', 1)}  |  SHIP: {fv(params.get('ship'), '', 1)}", ACCENT, 0.18),
        (f"   DCP: {fv(params.get('dcp'), '', 1)}  |  FireWx: FWI {fosberg_v}  Haines {haines_v}", ACCENT, 0.11),
        (f"   {'[SURFACE MODIFIED]' if params.get('surface_modified') else ''}{'  [CUSTOM SM]' if params.get('custom_storm_motion') else ''}", "#ff5555" if params.get('surface_modified') or params.get('custom_storm_motion') else BG, 0.04),
    ]
    
    for text, color, y_pos in thermo_rows:
        ax_params.text(0.01, y_pos, text,
                      transform=ax_params.transAxes,
                      fontsize=12, color=color, fontfamily="monospace",
                      fontweight="bold", va="top")
    
    # ── KINEMATIC TABLE ──
    ax_kin = fig.add_subplot(gs[1, 1:5])
    ax_kin.set_facecolor(BG)
    ax_kin.axis("off")
    
    kin_header = f"{'':3s}{'':8s} {'BWD':>8s}   {'SRH':>12s}   {'SRW':>8s}"
    
    bwd05 = fv(params.get("bwd_500m"), "kt")
    bwd1 = fv(params.get("bwd_1km"), "kt")
    bwd3 = fv(params.get("bwd_3km"), "kt")
    bwd6 = fv(params.get("bwd_6km"), "kt")
    srh05 = fv(params.get("srh_500m"), "m²/s²")
    srh1 = fv(params.get("srh_1km"), "m²/s²")
    srh3 = fv(params.get("srh_3km"), "m²/s²")
    esrh_v = fv(params.get("esrh"), "m²/s²")
    ebwd_v = fv(params.get("ebwd"), "kt")
    eil_bot_v = f"{params['eil_bot_h']:.0f}" if params.get("eil_bot_h") is not None else "---"
    eil_top_v = f"{params['eil_top_h']:.0f}" if params.get("eil_top_h") is not None else "---"
    
    # Storm-relative wind by layer
    def _srw_layer(depth_m):
        """Mean SR wind speed (kt) in 0-depth_m layer."""
        try:
            if params.get("rm_u") is None:
                return "---"
            rm_u_ms = params["rm_u"].to("m/s").magnitude
            rm_v_ms = params["rm_v"].to("m/s").magnitude
            mask = h_agl <= depth_m
            if np.sum(mask) < 2:
                return "---"
            sru = u.to("m/s").magnitude[mask] - rm_u_ms
            srv = v.to("m/s").magnitude[mask] - rm_v_ms
            sr_mean = np.mean(np.sqrt(sru**2 + srv**2)) * 1.94384
            return f"{sr_mean:.0f} kt"
        except Exception:
            return "---"
    
    srw05 = _srw_layer(500)
    srw1 = _srw_layer(1000)
    srw3 = _srw_layer(3000)
    srw6 = _srw_layer(6000)
    
    kin_rows = [
        ("KINEMATIC", ACCENT, 0.97),
        (kin_header, FG_FAINT, 0.89),
        (f"   {'0-500m:':8s} {bwd05:>8s}   {srh05:>12s}   {srw05:>8s}", "#ff5555", 0.81),
        (f"   {'0-1km:':8s} {bwd1:>8s}   {srh1:>12s}   {srw1:>8s}", "#ff3333", 0.73),
        (f"   {'0-3km:':8s} {bwd3:>8s}   {srh3:>12s}   {srw3:>8s}", "#ff8800", 0.65),
        (f"   {'0-6km:':8s} {bwd6:>8s}   {'':12s}   {srw6:>8s}", "#ffcc00", 0.57),
        (f"   EFFECTIVE: BWD {ebwd_v}  ESRH {esrh_v}", "#44ddaa", 0.49),
        (f"   EIL: {eil_bot_v}–{eil_top_v} m AGL", "#44ddaa", 0.41),
    ]
    
    # Bunkers
    if params.get("rm_u") is not None:
        rm_spd = np.sqrt(params["rm_u"]**2 + params["rm_v"]**2).to("knot")
        rm_toward = np.degrees(np.arctan2(params["rm_u"].magnitude,
                     params["rm_v"].magnitude)) % 360
        lm_spd = np.sqrt(params["lm_u"]**2 + params["lm_v"]**2).to("knot")
        lm_toward = np.degrees(np.arctan2(params["lm_u"].magnitude,
                     params["lm_v"].magnitude)) % 360
        
        kin_rows.append(("   BUNKERS STORM MOTION:", ACCENT, 0.33))
        kin_rows.append((
            f"    RM: {rm_toward:.0f}° @ {rm_spd.magnitude:.0f} kt  |  "
            f"LM: {lm_toward:.0f}° @ {lm_spd.magnitude:.0f} kt",
            FG_DIM, 0.25
        ))

    # Corfidi MCS motion vectors
    _cup_spd = params.get("corfidi_up_spd")
    _cdn_spd = params.get("corfidi_dn_spd")
    if _cup_spd is not None and _cdn_spd is not None:
        kin_rows.append(("   CORFIDI MCS MOTION:", "orange", 0.17))
        kin_rows.append((
            f"    UPW: {_cup_spd:.0f} kt  |  DNW: {_cdn_spd:.0f} kt",
            FG_DIM, 0.09
        ))
    
    for text, color, y_pos in kin_rows:
        ax_kin.text(0.01, y_pos, text,
                   transform=ax_kin.transAxes,
                   fontsize=12, color=color, fontfamily="monospace",
                   fontweight="bold", va="top")
    
    # ════════════════════════════════════════════════════════════════
    # MINI MAP INSET (station location) — inside kinematic panel
    # ════════════════════════════════════════════════════════════════
    ax_map = ax_kin.inset_axes([0.55, 0.05, 0.44, 0.90])  # right side of kinematic panel
    ax_map.set_facecolor(BG_PANEL)
    ax_map.plot(_conus_lon, _conus_lat, color=FG_FAINT, linewidth=0.8, zorder=2)
    ax_map.plot(lon, lat, marker="o", color="#ff3333", markersize=6,
                markeredgecolor=FG, markeredgewidth=0.8, zorder=10)
    # Apply zoom: 1.0 = full CONUS, higher = zoomed in on station
    _mz = max(1.0, min(float(map_zoom), 8.0))
    _full_lon_range = 63.0   # -128 to -65
    _full_lat_range = 28.0   # 23 to 51
    _half_lon = (_full_lon_range / _mz) / 2
    _half_lat = (_full_lat_range / _mz) / 2
    _cx = max(-128 + _half_lon, min(lon, -65 - _half_lon))
    _cy = max(23 + _half_lat, min(lat, 51 - _half_lat))
    ax_map.set_xlim(_cx - _half_lon, _cx + _half_lon)
    ax_map.set_ylim(_cy - _half_lat, _cy + _half_lat)
    ax_map.set_aspect(1.3)
    ax_map.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    for spine in ax_map.spines.values():
        spine.set_color(BORDER)
        spine.set_linewidth(0.5)

    # ── FOOTER ───────────────────────────────────────────────────────
    # Determine data source label from station_info
    stn_name_lower = info.get("name", "").lower()
    if "rap analysis" in stn_name_lower:
        source_label = "RAP Model Analysis (NCEI)"
    elif "acars" in stn_name_lower:
        source_label = "ACARS/AMDAR Aircraft Obs (IEM)"
    elif any(m in stn_name_lower for m in ("rap f", "hrrr f", "nam f", "gfs f", "sref f")):
        source_label = "BUFKIT Forecast (Iowa State)"
    else:
        source_label = "Iowa Environmental Mesonet"

    fig.text(0.04, 0.012,
             "VERTICAL PROFILE ANALYSIS TOOL | "
             f"Data: {source_label} | "
             f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%MZ')}",
             fontsize=10, color=FG_FAINT, fontfamily="monospace",
             fontweight="bold")
    fig.text(0.96, 0.012,
             f"Station: {station_name}  |  {lat:.4f}, {lon:.4f}",
             fontsize=10, color=ACCENT, ha="right", va="bottom",
             fontfamily="monospace", fontweight="bold")
    
    return fig


def plot_composite_sounding(profiles, title="Composite Sounding Overlay"):
    """
    Plot multiple sounding profiles overlaid on a single Skew-T.

    Parameters
    ----------
    profiles : list of dict
        Each dict has: {"data": sounding_data, "params": computed_params,
                        "label": str, "color": str (optional)}
    title : str
        Plot title.

    Returns
    -------
    matplotlib.figure.Figure
    """
    BG        = "#0d0d0d"
    FG        = "#e8e8e8"
    FG_DIM    = "#b0b0b0"
    GRID_CLR  = "#333333"
    BORDER    = "#444444"
    PALETTE   = ["#ef4444", "#3b82f6", "#22c55e", "#f59e0b", "#a855f7",
                 "#ec4899", "#06b6d4", "#84cc16"]

    fig = plt.figure(figsize=(14, 12), facecolor=BG)
    fig.patch.set_facecolor(BG)

    gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[2.2, 1],
                           hspace=0.1, wspace=0.12,
                           left=0.06, right=0.96, top=0.93, bottom=0.06)

    fig.suptitle(title, fontsize=15, fontweight="bold",
                 color=FG, y=0.98, fontfamily="monospace")

    # ── Skew-T ──
    skew = SkewT(fig, rotation=40, subplot=gs[0])
    skew.ax.set_facecolor(BG)
    skew.ax.set_aspect('auto')
    for spine in skew.ax.spines.values():
        spine.set_color(BORDER)
    skew.ax.tick_params(colors=FG_DIM, labelsize=8)
    skew.ax.set_xlabel("Temperature (°C)", color=FG_DIM, fontsize=9, fontfamily="monospace")
    skew.ax.set_ylabel("Pressure (hPa)", color=FG_DIM, fontsize=9, fontfamily="monospace")
    skew.ax.set_xlim(-40, 50)
    skew.ax.set_ylim(1050, 100)

    # Reference lines
    for t in range(-30, 40, 10):
        skew.ax.axvline(t, color=GRID_CLR, linestyle=":", linewidth=0.5, alpha=0.5)
    try:
        skew.plot_dry_adiabats(colors=GRID_CLR, linewidths=0.4, alpha=0.3)
        skew.plot_moist_adiabats(colors=GRID_CLR, linewidths=0.4, alpha=0.3)
        skew.plot_mixing_lines(colors=GRID_CLR, linewidths=0.4, alpha=0.3)
    except:
        pass

    legend_entries = []
    for i, prof in enumerate(profiles):
        color = prof.get("color", PALETTE[i % len(PALETTE)])
        label = prof.get("label", f"Profile {i+1}")
        data = prof["data"]
        p = data["pressure"]
        T = data["temperature"]
        Td = data["dewpoint"]

        skew.plot(p, T, color=color, linewidth=1.8, alpha=0.9)
        skew.plot(p, Td, color=color, linewidth=1.2, alpha=0.6, linestyle="--")
        legend_entries.append((color, label))

    # ── Hodograph ──
    ax_hodo = fig.add_subplot(gs[1], projection="polar")
    ax_hodo.set_facecolor(BG)
    # Convert to non-polar for hodograph
    fig.delaxes(ax_hodo)
    ax_hodo = fig.add_subplot(gs[1])
    ax_hodo.set_facecolor(BG)
    ax_hodo.set_aspect("equal")
    for spine in ax_hodo.spines.values():
        spine.set_color(BORDER)
    ax_hodo.tick_params(colors=FG_DIM, labelsize=7)
    ax_hodo.set_title("Hodograph Overlay", color=FG, fontsize=10,
                       fontfamily="monospace", fontweight="bold")

    for i, prof in enumerate(profiles):
        color = prof.get("color", PALETTE[i % len(PALETTE)])
        data = prof["data"]
        params = prof.get("params", {})
        if params.get("u") is not None and params.get("v") is not None:
            u_kt = params["u"].to("knot").magnitude
            v_kt = params["v"].to("knot").magnitude
            ax_hodo.plot(u_kt, v_kt, color=color, linewidth=1.5, alpha=0.8)

    ax_hodo.axhline(0, color=GRID_CLR, linewidth=0.5)
    ax_hodo.axvline(0, color=GRID_CLR, linewidth=0.5)
    ax_hodo.set_xlabel("U (kt)", color=FG_DIM, fontsize=8, fontfamily="monospace")
    ax_hodo.set_ylabel("V (kt)", color=FG_DIM, fontsize=8, fontfamily="monospace")

    # ── Legend ──
    for i, (color, label) in enumerate(legend_entries):
        fig.text(0.06, 0.03 - i * 0.018, f"━━ {label}",
                 fontsize=9, color=color, fontfamily="monospace", fontweight="bold")

    return fig


# ─────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Sounding Analysis Tool — fetch & plot real upper-air data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Data sources (--source):
  obs      Observed radiosonde (IEM / UWyo) — default
  rap      RAP model analysis — any lat/lon, CONUS (needs siphon)
  bufkit   BUFKIT forecast models — station-based (Iowa State)
  acars    ACARS/AMDAR aircraft obs — airport-based (IEM)

Examples:
  python sounding.py --station OUN                                 # Observed
  python sounding.py --source rap --lat 35.2 --lon -97.4           # RAP at point
  python sounding.py --source bufkit --model hrrr --station OUN    # HRRR forecast
  python sounding.py --source acars --station KDFW                 # Aircraft obs
""",
    )
    parser.add_argument("--station", type=str, default=None,
                        help="3-letter station ID (e.g. OUN) or ICAO airport "
                             "code for ACARS (e.g. KDFW)")
    parser.add_argument("--date", type=str, default=None,
                        help="Date/time as YYYYMMDDHH (e.g. 2024061200). "
                             "Defaults to most recent sounding time.")
    parser.add_argument("--lat", type=float, default=None,
                        help="Latitude — required for rap if no --station")
    parser.add_argument("--lon", type=float, default=None,
                        help="Longitude — required for rap if no --station")
    parser.add_argument("--source", type=str, default=None,
                        choices=list(DATA_SOURCES.keys()),
                        help="Data source (default: obs). See below for details.")
    parser.add_argument("--model", type=str, default="rap",
                        choices=list(BUFKIT_MODELS.keys()),
                        help="BUFKIT model name (default: rap). "
                             "Only used with --source bufkit.")
    parser.add_argument("--fhour", type=int, default=0,
                        help="BUFKIT forecast hour, 0 = analysis (default: 0). "
                             "Only used with --source bufkit.")
    parser.add_argument("--output", type=str, default=None,
                        help="Output filename (default: auto-generated)")
    parser.add_argument("--list-stations", action="store_true",
                        help="List all available stations")
    parser.add_argument("--list-sources", action="store_true",
                        help="List all available data sources")

    args = parser.parse_args()

    # ── Informational listings ──────────────────────────────────────
    if args.list_stations:
        print("\nAvailable Sounding Stations:")
        print(f"  {'ID':5s} {'Name':30s} {'Lat':>8s} {'Lon':>10s}")
        print("  " + "-" * 55)
        for code in sorted(STATIONS):
            name, lat, lon = STATIONS[code]
            print(f"  {code:5s} {name:30s} {lat:8.2f} {lon:10.2f}")
        return

    if args.list_sources:
        print("\nAvailable Data Sources (--source):")
        for key, desc in DATA_SOURCES.items():
            print(f"  {key:8s}  {desc}")
        print("\nBUFKIT Models (--model, used with --source bufkit):")
        for key, desc in BUFKIT_MODELS.items():
            print(f"  {key:8s}  {desc}")
        return

    # ── Interactive source selection when no CLI flags given ────────
    source = args.source
    if source is None and args.station is None and args.lat is None:
        print("\n  +====================================================+")
        print("  |       VERTICAL PROFILE ANALYSIS TOOL              |")
        print("  +====================================================+")
        print("  |  Select a data source:                            |")
        print("  |                                                   |")
        print("  |  [1] Observed radiosonde   (IEM / UWyo)           |")
        print("  |  [2] RAP model analysis    (any lat/lon, CONUS)   |")
        print("  |  [3] BUFKIT forecast model (station-based)        |")
        print("  |  [5] ACARS aircraft obs    (airport-based)        |")
        print("  |                                                   |")
        print("  |  [0] Auto-select (tornado risk scan, observed)    |")
        print("  +====================================================+")
        try:
            choice = input("\n  Enter choice [0-5]: ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            return

        source_map = {
            "0": None, "1": "obs", "2": "rap", "3": "bufkit",
            "4": "acars",
        }
        source = source_map.get(choice)
        if choice not in source_map:
            print(f"  Invalid choice '{choice}'. Defaulting to observed.")
            source = "obs"

        # Additional interactive prompts depending on source
        if source == "rap":
            try:
                lat_str = input("  Latitude  (e.g. 35.22): ").strip()
                lon_str = input("  Longitude (e.g. -97.46): ").strip()
                args.lat = float(lat_str)
                args.lon = float(lon_str)
            except (ValueError, EOFError, KeyboardInterrupt):
                print("  Invalid coordinates. Aborting.")
                return

        if source in ("obs", "bufkit"):
            if args.station is None:
                try:
                    stn = input(
                        "  Station ID (e.g. OUN) or press Enter for auto-select: "
                    ).strip()
                    if stn:
                        args.station = stn
                except (EOFError, KeyboardInterrupt):
                    print()
                    return

        if source == "bufkit":
            print(f"\n  Available BUFKIT models:")
            for i, (key, desc) in enumerate(BUFKIT_MODELS.items(), 1):
                print(f"    [{i}] {key:8s} - {desc}")
            try:
                m_choice = input(
                    "  Model [1-6, default=1 (rap)]: "
                ).strip()
                if m_choice:
                    model_keys = list(BUFKIT_MODELS.keys())
                    idx = int(m_choice) - 1
                    if 0 <= idx < len(model_keys):
                        args.model = model_keys[idx]
            except (ValueError, EOFError, KeyboardInterrupt):
                pass
            try:
                fh = input("  Forecast hour [default=0]: ").strip()
                if fh:
                    args.fhour = int(fh)
            except (ValueError, EOFError, KeyboardInterrupt):
                pass

        if source == "acars":
            if args.station is None:
                try:
                    args.station = input(
                        "  Airport ICAO code (e.g. KDFW, KORD): "
                    ).strip()
                except (EOFError, KeyboardInterrupt):
                    print()
                    return

    if source is None:
        source = "obs"

    # ── Determine station / coordinates ─────────────────────────────
    station = args.station
    lat = args.lat
    lon = args.lon

    if station is None and lat is not None and lon is not None:
        if source == "obs":
            station = find_nearest_station(lat, lon)
            print(f"  Nearest station to ({lat}, {lon}): {station} "
                  f"({STATIONS[station][0]})")

    # ── Determine time ──────────────────────────────────────────────
    if args.date:
        try:
            dt = datetime.strptime(args.date, "%Y%m%d%H").replace(
                tzinfo=timezone.utc
            )
        except ValueError:
            print(f"  ERROR: Invalid date format '{args.date}'. Use YYYYMMDDHH.")
            sys.exit(1)
    else:
        dt = get_latest_sounding_time()
        print(f"  Using most recent sounding time: {dt.strftime('%Y-%m-%d %HZ')}")

    # ── Auto-select station (obs only, no station given) ────────────
    if source == "obs" and station is None:
        station = find_highest_tornado_risk(dt)

    if station:
        station = station.upper()

    # Validate station for obs source
    if source == "obs" and station not in STATION_WMO:
        print(f"  ERROR: Unknown station '{station}'. "
              f"Use --list-stations to see options.")
        sys.exit(1)

    # ── Build descriptive label ─────────────────────────────────────
    if source == "rap" and lat is not None:
        label = f"{source.upper()} ({lat:.2f}, {lon:.2f})"
    elif station:
        label = f"{station} ({STATIONS.get(station, (station,))[0]})"
    else:
        label = source.upper()

    print(f"\n{'='*60}")
    print(f"  SOUNDING ANALYSIS: {label}")
    print(f"  Source: {DATA_SOURCES.get(source, source)}")
    print(f"  Valid: {dt.strftime('%B %d, %Y at %HZ')}")
    if source == "bufkit":
        print(f"  Model: {args.model.upper()}  |  Forecast hour: f{args.fhour:03d}")
    print(f"{'='*60}\n")

    # ── Fetch data ──────────────────────────────────────────────────
    print("  [1/3] Fetching sounding data...")
    try:
        data = fetch_sounding(
            station_id=station,
            dt=dt,
            source=source,
            lat=lat,
            lon=lon,
            model=args.model,
            fhour=args.fhour,
        )
    except Exception as e:
        print(f"\n  ERROR: Could not fetch data: {e}")
        print(f"  Try a different station/location or time.")
        if source == "obs":
            print(f"  Tip: Use --date to specify a past sounding "
                  f"(e.g. --date "
                  f"{(dt - timedelta(hours=12)).strftime('%Y%m%d%H')})")
        sys.exit(1)

    n_levels = len(data["pressure"])
    sfc_p = data["pressure"][0].magnitude
    top_p = data["pressure"][-1].magnitude
    print(f"  [OK] Retrieved {n_levels} levels "
          f"({sfc_p:.0f} hPa to {top_p:.0f} hPa)")

    # ── Compute parameters ──────────────────────────────────────────
    print("  [2/3] Computing thermodynamic & kinematic parameters...")
    params = compute_parameters(data)

    sb_cape_val = params.get("sb_cape")
    if sb_cape_val is not None and hasattr(sb_cape_val, "magnitude"):
        print(f"  [OK] SB CAPE: {sb_cape_val.magnitude:.0f} J/kg")
    mu_cape_val = params.get("mu_cape")
    if mu_cape_val is not None and hasattr(mu_cape_val, "magnitude"):
        print(f"    MU CAPE: {mu_cape_val.magnitude:.0f} J/kg")

    # ── Plot ────────────────────────────────────────────────────────
    print("  [3/3] Generating analysis plot...")
    # Use station id for plot; for lat/lon sources, create a synthetic id
    plot_id = station if station else f"{source.upper()}"
    fig = plot_sounding(data, params, plot_id, dt)

    # ── Save ────────────────────────────────────────────────────────
    if args.output:
        outfile = args.output
    else:
        tag = station or f"{source}_{lat:.1f}_{lon:.1f}".replace(".", "p").replace("-", "m")
        outfile = f"sounding_{tag}_{dt.strftime('%Y%m%d_%HZ')}.png"

    fig.savefig(outfile, dpi=180, facecolor="#0d0d0d")
    plt.close(fig)
    print(f"\n  [OK] Saved to: {outfile}")
    print(f"    Resolution: 180 DPI")
    print(f"\n{'='*60}")
    print("  Done!")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
