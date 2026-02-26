"""
VERTICAL PROFILE ANALYSIS TOOL
===============================
Fetches real upper-air sounding data from multiple sources
and produces a comprehensive analysis plot similar to SounderPy.

Data Sources:
  - Observed radiosondes (IEM / University of Wyoming)
  - RAP model analysis for any lat/lon (NCEI THREDDS, requires siphon)
  - BUFKIT forecast soundings: HRRR, RAP, NAM, GFS, etc. (Iowa State)
  - ERA5 global reanalysis for any lat/lon (CDS API, requires cdsapi + key)
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
  python sounding.py --source era5 --lat 36.4 --lon -99.4 --date 2020050100
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
    
    url = (
        f"https://weather.uwyo.edu/cgi-bin/sounding.py?"
        f"region=naconf&TYPE=TEXT%3ALIST&YEAR={dt.year}&MONTH={dt.month:02d}"
        f"&FROM={dt.day:02d}{dt.hour:02d}&TO={dt.day:02d}{dt.hour:02d}"
        f"&STNM={wmo}"
    )
    
    print(f"  Fetching from UWyo: {url}")
    resp = requests.get(url, timeout=20)
    resp.raise_for_status()
    html = resp.text
    
    if "Can't get" in html or "No data" in html:
        raise ValueError(f"No sounding data available for {station_id} at {dt}")
    
    pre_start = html.find("<pre>")
    pre_end = html.find("</pre>")
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
        if not line or line.startswith("</pre>"):
            break
        parts = line.split()
        if len(parts) < 7:
            continue
        try:
            p = float(parts[0])
            h = float(parts[1])
            t = float(parts[2])
            td = float(parts[3])
            wd = float(parts[5])
            ws = float(parts[6])
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
        "wind_speed": np.array(wind_spd) * units.knot,
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


# ─── ERA5 REANALYSIS ──────────────────────────────────────────────
# ERA5 pressure levels (hPa) for retrieval
_ERA5_LEVELS = [
    "1000", "975", "950", "925", "900", "875", "850", "825", "800",
    "775", "750", "700", "650", "600", "550", "500", "450", "400",
    "350", "300", "250", "225", "200", "175", "150", "125", "100",
    "70", "50", "30", "20", "10",
]


def fetch_era5_sounding(lat, lon, dt):
    """
    Fetch ERA5 global reanalysis sounding for any lat/lon.
    Requires a free CDS API key — see https://cds.climate.copernicus.eu

    Setup:
      1. Create a free account at https://cds.climate.copernicus.eu
      2. Accept the ERA5 licence terms.
      3. Copy your API key from your profile page.
      4. Create file  ~/.cdsapirc  with:
            url: https://cds.climate.copernicus.eu/api
            key: <YOUR-API-KEY>
    """
    try:
        import cdsapi
    except ImportError:
        raise ImportError(
            "cdsapi is required for ERA5 data. Install with: pip install cdsapi"
        )
    try:
        import netCDF4
    except ImportError:
        raise ImportError(
            "netCDF4 is required for ERA5 data. Install with: pip install netCDF4"
        )

    import tempfile, os

    print(f"  Fetching ERA5 reanalysis for ({lat:.2f}, {lon:.2f}) "
          f"at {dt:%Y-%m-%d %H}Z...")
    print("  (ERA5 requests may take 1-5 minutes to process)")

    c = cdsapi.Client(progress=False)

    # Small bounding box around the point (0.5° pad)
    area = [lat + 0.5, lon - 0.5, lat - 0.5, lon + 0.5]

    tmp_file = os.path.join(tempfile.gettempdir(), "era5_sounding_tmp.nc")

    c.retrieve(
        "reanalysis-era5-pressure-levels",
        {
            "product_type": "reanalysis",
            "variable": [
                "temperature",
                "u_component_of_wind",
                "v_component_of_wind",
                "specific_humidity",
                "geopotential",
            ],
            "pressure_level": _ERA5_LEVELS,
            "year":  f"{dt.year}",
            "month": f"{dt.month:02d}",
            "day":   f"{dt.day:02d}",
            "time":  f"{dt.hour:02d}:00",
            "area":  area,
            "data_format": "netcdf",
        },
        tmp_file,
    )

    ds = netCDF4.Dataset(tmp_file, "r")

    # Variable names in ERA5 netCDF: t, u, v, q, z  (or long names)
    def _v(ds, *names):
        for n in names:
            if n in ds.variables:
                return ds.variables[n][:]
        raise KeyError(f"None of {names} in ERA5 file. "
                       f"Available: {list(ds.variables.keys())}")

    pres_pa  = _v(ds, "level", "pressure_level", "isobaric")
    temp_k   = _v(ds, "t", "temperature")
    u_ms     = _v(ds, "u", "u_component_of_wind")
    v_ms     = _v(ds, "v", "v_component_of_wind")
    q_kgkg   = _v(ds, "q", "specific_humidity")
    z_m2s2   = _v(ds, "z", "geopotential")
    ds.close()

    # Remove temp file
    try:
        os.remove(tmp_file)
    except OSError:
        pass

    # ERA5 arrays shape: (time, level, lat, lon)   — squeeze to (level,)
    temp_k  = temp_k.squeeze()
    u_ms    = u_ms.squeeze()
    v_ms    = v_ms.squeeze()
    q_kgkg  = q_kgkg.squeeze()
    z_m2s2  = z_m2s2.squeeze()

    # If still >1-D after squeeze (multiple lat/lon in box), pick nearest
    if temp_k.ndim > 1:
        # Take centre point of the small box
        mid_lat = temp_k.shape[-2] // 2
        mid_lon = temp_k.shape[-1] // 2
        temp_k  = temp_k[:, mid_lat, mid_lon]
        u_ms    = u_ms[:, mid_lat, mid_lon]
        v_ms    = v_ms[:, mid_lat, mid_lon]
        q_kgkg  = q_kgkg[:, mid_lat, mid_lon]
        z_m2s2  = z_m2s2[:, mid_lat, mid_lon]

    pres_hpa = np.array(pres_pa, dtype=float)        # levels are already in hPa
    temp_c   = np.array(temp_k, dtype=float) - 273.15
    hgt_m    = np.array(z_m2s2, dtype=float) / 9.80665   # geopotential → metres

    # Specific humidity → dewpoint
    # e = q * p / (0.622 + 0.378*q)     (vapour pressure in hPa)
    q = np.clip(np.array(q_kgkg, dtype=float), 1e-10, None)
    e_hpa = q * pres_hpa / (0.622 + 0.378 * q)
    # Magnus invert: Td = b * ln(e/6.112) / (a - ln(e/6.112))
    a, b = 17.625, 243.04
    ln_e = np.log(np.clip(e_hpa / 6.112, 1e-10, None))
    td_c = (b * ln_e) / (a - ln_e)

    # Wind → speed & direction
    wspd_kt  = np.sqrt(np.array(u_ms, dtype=float)**2
                       + np.array(v_ms, dtype=float)**2) * 1.94384
    wdir_deg = (np.degrees(np.arctan2(-np.array(u_ms, dtype=float),
                                       -np.array(v_ms, dtype=float))) + 360) % 360

    # Sort descending pressure (surface first), keep ≥ 10 hPa
    order = np.argsort(pres_hpa)[::-1]
    mask = pres_hpa[order] >= 10
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
            "name": f"ERA5 Reanalysis ({lat:.2f}, {lon:.2f})",
        },
    }
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
# UNIFIED FETCH DISPATCHER
# ─────────────────────────────────────────────────────────────────────
# Recognised source keywords (for --source flag / interactive menu)
DATA_SOURCES = {
    "obs":    "Observed radiosonde (IEM / UWyo)",
    "rap":    "RAP model analysis - any lat/lon, CONUS (NCEI THREDDS)",
    "bufkit": "BUFKIT forecast models - station-based (Iowa State)",
    "era5":   "ERA5 global reanalysis - any lat/lon, 1940-present (CDS API)",
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
        One of: obs, rap, bufkit, era5, acars.
    lat, lon : float or None
        Required for rap / era5 (point-based sources).
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

    # ── ERA5 ────────────────────────────────────────────────────────
    if source == "era5":
        if lat is None or lon is None:
            if station_id and station_id.upper() in STATIONS:
                _, lat, lon = STATIONS[station_id.upper()]
            else:
                raise ValueError(
                    "ERA5 source requires --lat / --lon (or a known --station)."
                )
        return fetch_era5_sounding(lat, lon, dt)

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
    Fetch sounding for a station and compute a quick tornado composite score.
    Returns (score, cape, srh, bwd) or None on failure.
    Score = STP if computable, else composite of CAPE * SRH * BWD.
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

        # 0-1 km SRH
        h_agl = (h - h[0]).to("meter")
        try:
            rm_u, rm_v, _, _ = mpcalc.bunkers_storm_motion(p, u_wind, v_wind, h)
            pos_srh, neg_srh, total_srh = mpcalc.storm_relative_helicity(
                h_agl, u_wind, v_wind, 1000 * units.meter,
                storm_u=rm_u, storm_v=rm_v
            )
            srh_val = total_srh.magnitude if hasattr(total_srh, "magnitude") else 0
        except Exception:
            srh_val = 0

        # 0-6 km BWD
        try:
            bwd_u, bwd_v = mpcalc.bulk_shear(p, u_wind, v_wind,
                                             height=h_agl,
                                             depth=6000 * units.meter)
            bwd_val = np.sqrt(bwd_u**2 + bwd_v**2).to("knot").magnitude
        except Exception:
            bwd_val = 0

        # STP-like composite
        # STP = (CAPE/1500) * (SRH/150) * (BWD/20 kt) * ((200 + CIN)/150)
        cape_term = min(cape_val / 1500.0, 3.0)
        srh_term = min(srh_val / 150.0, 3.0)
        bwd_term = min(bwd_val / 20.0, 3.0) if bwd_val >= 12 else 0
        cin_term = min((200.0 + cin_val) / 150.0, 1.0) if cin_val > -250 else 0

        # Primary: STP-like multiplicative score (may be 0 in benign wx)
        stp_score = cape_term * srh_term * bwd_term * cin_term

        # Secondary: additive composite so individual parameters still
        # differentiate stations when the STP score ties at 0.
        raw_score = (max(cape_val, 0) / 1500.0
                     + max(srh_val, 0) / 150.0
                     + bwd_val / 40.0)

        return (float(stp_score), float(raw_score), cape_val, srh_val, bwd_val)

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

    # Collect results: (station_id, stp, raw, cape, srh, bwd, name)
    results = []
    for sid in stations:
        result = _quick_tornado_score(sid, dt)
        if result is None:
            continue
        stp_score, raw_score, cape, srh, bwd = result
        name = STATIONS.get(sid, (sid,))[0]
        results.append((sid, stp_score, raw_score, cape, srh, bwd, name))

    if not results:
        raise ValueError("Could not fetch data from any station")

    # Sort by (STP desc, Raw desc) so best is at the top
    results.sort(key=lambda r: (r[1], r[2]), reverse=True)

    # ── Print ranked table ──────────────────────────────────────────
    print(f"  {'#':>3s}  {'ID':5s}  {'Station':22s}  {'STP':>7s}  {'CAPE':>8s}  "
          f"{'0-1 SRH':>8s}  {'0-6 BWD':>7s}")
    print(f"  {'-'*68}")

    for idx, (sid, stp, raw, cape, srh, bwd, name) in enumerate(results, 1):
        highlight = " <--" if idx == 1 else ""
        print(f"  {idx:3d}  {sid:5s}  {name:22s}  {stp:7.2f}  "
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
# CALCULATIONS
# ─────────────────────────────────────────────────────────────────────
def compute_parameters(data):
    """Compute a comprehensive set of thermodynamic & kinematic parameters."""
    p = data["pressure"]
    T = data["temperature"]
    Td = data["dewpoint"]
    h = data["height"]
    wdir = data["wind_direction"]
    wspd = data["wind_speed"]
    
    params = {}
    
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
    except Exception as e:
        print(f"  Warning: ML parcel calc failed: {e}")
        params["ml_cape"] = 0 * units("J/kg")
        params["ml_cin"] = 0 * units("J/kg")
    
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
    
    # STP (now that we have SRH and BWD)
    try:
        params["stp"] = float(mpcalc.significant_tornado(
            params["sb_cape"],
            params["sb_lcl_p"] if params.get("sb_lcl_p") is not None else p[0],
            params["srh_1km"],
            params["bwd_6km"].to("m/s")
        ).magnitude)
    except:
        params["stp"] = 0
    
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
    
    # DCAPE (Downdraft CAPE)
    try:
        dcape_val, dtemp = mpcalc.downdraft_cape(p, T, Td)
        params["dcape"] = dcape_val
    except:
        params["dcape"] = None
    
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
    
    return params


# ─────────────────────────────────────────────────────────────────────
# PLOTTING
# ─────────────────────────────────────────────────────────────────────
def plot_sounding(data, params, station_id, dt):
    """Create a comprehensive sounding analysis figure (dark theme)."""
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
    
    # ── Theme colors (dark) ──────────────────────────────────────────
    BG        = "#0d0d0d"
    BG_PANEL  = "#141414"
    FG        = "#e8e8e8"
    FG_DIM    = "#b0b0b0"
    FG_FAINT  = "#707070"
    GRID_CLR  = "#333333"
    BORDER    = "#444444"
    ACCENT    = "#55bbee"
    
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
    fig = plt.figure(figsize=(26, 12), facecolor=BG)
    fig.patch.set_facecolor(BG)
    
    # Define grid for layout — 4 columns: SkewT | Hodograph | SRW | SRH
    gs = gridspec.GridSpec(
        2, 4, figure=fig,
        width_ratios=[1.8, 1.6, 0.38, 0.42],
        height_ratios=[3.0, 1.0],
        hspace=0.12, wspace=0.12,
        left=0.05, right=0.97, top=0.94, bottom=0.04
    )
    
    # ── TITLE BAR ────────────────────────────────────────────────────
    title_str = (
        f"OBSERVED UPPER-AIR SOUNDING | {station_id} | "
        f"VALID: {dt.strftime('%m/%d/%Y %HZ')}"
    )
    fig.suptitle(title_str, fontsize=17, fontweight="bold",
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
    skew.ax.tick_params(colors=FG_DIM, labelsize=8, width=1.2)
    skew.ax.set_xlabel("Temperature (°C)", color=FG_DIM, fontsize=9,
                       fontfamily="monospace", fontweight="bold")
    skew.ax.set_ylabel("Pressure (hPa)", color=FG_DIM, fontsize=9,
                       fontfamily="monospace", fontweight="bold")
    
    # Reference lines (subtle)
    skew.plot_dry_adiabats(colors=GRID_CLR, alpha=0.4, linewidth=0.5)
    skew.plot_moist_adiabats(colors=GRID_CLR, alpha=0.4, linewidth=0.5)
    skew.plot_mixing_lines(colors=GRID_CLR, alpha=0.25, linewidth=0.4)
    
    # Panel label
    skew.ax.text(0.02, 0.98, "SKEW-T LOG-P", transform=skew.ax.transAxes,
                 fontsize=11, color=FG, fontfamily="monospace",
                 fontweight="bold", va="top", ha="left", alpha=0)  # hidden, merged into legend title
    
    # Temperature & Dewpoint (main traces)
    skew.plot(p, T, "r", linewidth=3.5, zorder=6)
    skew.plot(p, Td, color="#4488ff", linewidth=3.5, zorder=6)
    
    # Wet-bulb temperature
    if params.get("wetbulb") is not None:
        skew.plot(p, params["wetbulb"], color="#44bbdd", linewidth=1.5,
                  linestyle="--", alpha=0.7, zorder=5)
    
    # Virtual temperature
    if params.get("virtual_temp") is not None:
        skew.plot(p, params["virtual_temp"], color="#ff8888", linewidth=1.0,
                  linestyle="--", alpha=0.4, zorder=4)
    
    # Parcel traces
    if "sb_profile" in params:
        skew.plot(p, params["sb_profile"], color="#ff8800", linewidth=2.5,
                  linestyle="--", zorder=7, label="SB Parcel")
        try:
            skew.shade_cape(p, T, params["sb_profile"], alpha=0.25)
            skew.shade_cin(p, T, params["sb_profile"], alpha=0.18)
        except:
            pass
    
    if "ml_profile" in params:
        skew.plot(p, params["ml_profile"], color="#dd44dd", linewidth=2.0,
                  linestyle="--", alpha=0.8, zorder=7, label="ML Parcel")
    
    if "mu_profile" in params:
        idx = params.get("mu_start_idx", 0)
        skew.plot(p[idx:], params["mu_profile"], color="#aaaaaa",
                  linewidth=1.5, linestyle=":", alpha=0.6, zorder=7,
                  label="MU Parcel")
    
    # Wind barbs
    barb_interval = max(1, len(p) // 40)
    try:
        skew.plot_barbs(
            p[::barb_interval], u[::barb_interval], v[::barb_interval],
            color=FG, length=6, linewidth=0.9,
            xloc=1.05, zorder=10
        )
    except:
        pass
    
    # Annotate key levels
    pe = path_effects.withStroke(linewidth=3, foreground=BG)
    
    if params.get("sb_lcl_p") is not None:
        skew.ax.axhline(y=params["sb_lcl_p"].magnitude, color="#44cc44",
                        linestyle="--", alpha=0.6, linewidth=1.0)
        skew.ax.text(
            skew.ax.get_xlim()[1] - 2, params["sb_lcl_p"].magnitude,
            f"←SBLCL ({fv(params['sb_lcl_p'],'hPa')})",
            color="#44cc44", fontsize=8, va="center",
            fontfamily="monospace", fontweight="bold", path_effects=[pe]
        )
    
    if params.get("sb_lfc_p") is not None:
        skew.ax.axhline(y=params["sb_lfc_p"].magnitude, color="#ddaa22",
                        linestyle="--", alpha=0.5, linewidth=1.0)
        skew.ax.text(
            skew.ax.get_xlim()[1] - 2, params["sb_lfc_p"].magnitude,
            f"←LFC", color="#ddaa22", fontsize=8, va="center",
            fontfamily="monospace", fontweight="bold", path_effects=[pe]
        )
    
    if params.get("sb_el_p") is not None:
        skew.ax.axhline(y=params["sb_el_p"].magnitude, color="#4499ee",
                        linestyle="--", alpha=0.5, linewidth=1.0)
        skew.ax.text(
            skew.ax.get_xlim()[1] - 2, params["sb_el_p"].magnitude,
            f"←EL", color="#4499ee", fontsize=8, va="center",
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
            fontsize=8, va="center", fontfamily="monospace",
            fontweight="bold", path_effects=[pe]
        )
    
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
                fontsize=7.5, color=ACCENT, fontfamily="monospace",
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
        fontsize=7.5, color=FG_DIM, fontfamily="monospace",
        fontweight="bold", va="top", ha="left",
        bbox=dict(boxstyle="round,pad=0.15", facecolor=BG,
                 edgecolor="none", alpha=0.85)
    )
    
    # Legend
    legend = skew.ax.legend(
        loc="upper left", fontsize=8, facecolor=BG,
        edgecolor=BORDER, labelcolor=FG,
        framealpha=0.9, borderpad=0.5,
        title="SKEW-T LOG-P", title_fontproperties={"size": 10, "weight": "bold",
        "family": "monospace"}
    )
    legend.get_title().set_color(FG)
    
    # ════════════════════════════════════════════════════════════════
    # HODOGRAPH
    # ════════════════════════════════════════════════════════════════
    ax_hodo = fig.add_subplot(gs[0, 1])
    ax_hodo.set_facecolor(BG_PANEL)

    # --- Use the pre-computed 100-m interpolated winds from params ---
    hodo_u = params["u_interp"].to("knot").magnitude.copy()
    hodo_v = params["v_interp"].to("knot").magnitude.copy()
    hodo_z = params["h_interp"].to("meter").magnitude.copy()

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
                         weight='bold', fontsize=11, color='black',
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
                             weight='bold', fontsize=11, color='black',
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
        rm_u_kt = params["rm_u"].to("knot").magnitude
        rm_v_kt = params["rm_v"].to("knot").magnitude
        lm_u_kt = params["lm_u"].to("knot").magnitude
        lm_v_kt = params["lm_v"].to("knot").magnitude
        mw_u_kt = params["mw_u"].to("knot").magnitude
        mw_v_kt = params["mw_v"].to("knot").magnitude

        sm_u_kt = rm_u_kt
        sm_v_kt = rm_v_kt

        # RM / LM / MW text labels
        hodo.ax.text(rm_u_kt + 0.5, rm_v_kt - 0.5, 'RM',
                     weight='bold', ha='left', fontsize=14, zorder=7,
                     alpha=0.9, color=FG)
        hodo.ax.text(lm_u_kt + 0.5, lm_v_kt - 0.5, 'LM',
                     weight='bold', ha='left', fontsize=14, zorder=7,
                     alpha=0.9, color=FG)
        hodo.ax.text(mw_u_kt + 0.5, mw_v_kt - 0.5, 'MW',
                     weight='bold', ha='left', fontsize=14, zorder=7,
                     alpha=0.9, color=FG)

        # DTM (Deviant Tornado Motion)
        dtm_u_kt = mw_u_kt + (rm_v_kt - mw_v_kt)
        dtm_v_kt = mw_v_kt - (rm_u_kt - mw_u_kt)
        hodo.ax.text(dtm_u_kt, dtm_v_kt + 2, 'DTM',
                     weight='bold', fontsize=10, color='brown',
                     ha='center', zorder=7)
        hodo.ax.plot(dtm_u_kt, dtm_v_kt, marker='v', color='brown',
                     markersize=8, zorder=7, alpha=0.8, ls='')

        # SM arrow from origin to storm motion
        hodo.ax.arrow(0, 0, sm_u_kt - 0.3, sm_v_kt - 0.3,
                      linewidth=3, color=FG, alpha=0.2,
                      label='SM Vector', length_includes_head=True,
                      head_width=0.6)

        # --- Effective inflow layer SRH fill (0-3 km, using RM) ---
        eil_bot_idx = 0
        eil_top_idx = min(30, n_full - 1)  # 3 km
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
            'lightblue', alpha=0.3, zorder=2, label='0-3 SRH')

        # --- MCS motion markers (upshear / downshear) ---
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

        sm_lines = [
            f"SM: RIGHT MOVING | {_wind_to_dir_str(rm_u_kt, rm_v_kt)} @ {_wind_spd(rm_u_kt, rm_v_kt)} kts",
            f"RM: {_wind_to_dir_str(rm_u_kt, rm_v_kt)} @ {_wind_spd(rm_u_kt, rm_v_kt)} kts",
            f"LM: {_wind_to_dir_str(lm_u_kt, lm_v_kt)} @ {_wind_spd(lm_u_kt, lm_v_kt)} kts",
            f"MW: {_wind_to_dir_str(mw_u_kt, mw_v_kt)} @ {_wind_spd(mw_u_kt, mw_v_kt)} kts",
            f"DTM: {_wind_to_dir_str(dtm_u_kt, dtm_v_kt)} @ {_wind_spd(dtm_u_kt, dtm_v_kt)} kts",
            f"US: {_wind_to_dir_str(us_u, us_v)} @ {_wind_spd(us_u, us_v)} kts",
            f"DS: {_wind_to_dir_str(ds_u, ds_v)} @ {_wind_spd(ds_u, ds_v)} kts",
        ]

    sm_text = "\n".join(sm_lines) if sm_lines else ""
    if sm_text:
        ax_hodo.text(
            0.02, 0.98, sm_text,
            transform=ax_hodo.transAxes, fontsize=9,
            color=FG, fontfamily="monospace", fontweight="bold",
            va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.3", facecolor=BG,
                     edgecolor=BORDER, alpha=0.92)
        )
    
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
        ax_srw.set_xlabel("SRW (kt)", color=FG_DIM, fontsize=8,
                         fontfamily="monospace", fontweight="bold")
        ax_srw.set_ylabel("Height AGL (km)", color=FG_DIM, fontsize=8,
                         fontfamily="monospace", fontweight="bold")
        ax_srw.set_title("STORM-REL\nWIND", color=FG, fontsize=8,
                        fontfamily="monospace", fontweight="bold", pad=3)
    else:
        ax_srw.text(0.5, 0.5, "N/A", transform=ax_srw.transAxes,
                   color=FG_FAINT, ha="center", fontsize=12)
    
    ax_srw.tick_params(colors=FG_DIM, labelsize=7, width=1.2)
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
                              fontsize=7, color=color, fontweight="bold",
                              fontfamily="monospace")
        
        # Axis limits
        ax_sw.set_xlim(0, 105)
        ax_sw.set_ylim(0, min(sw_h.max(), 6))
        ax_sw.set_xlabel("Streamwiseness (%)", color=FG_DIM, fontsize=8,
                         fontfamily="monospace", fontweight="bold")
        ax_sw.set_ylabel("Height AGL (km)", color=FG_DIM, fontsize=8,
                         fontfamily="monospace", fontweight="bold")
        ax_sw.set_title("STREAM-\nWISENESS", color=FG, fontsize=8,
                        fontfamily="monospace", fontweight="bold", pad=3)
        
        # Legend
        leg = ax_sw.legend(loc="lower right", fontsize=6, facecolor=BG,
                           edgecolor=BORDER, labelcolor=FG, framealpha=0.9)
    else:
        ax_sw.text(0.5, 0.5, "N/A", transform=ax_sw.transAxes,
                   color=FG_FAINT, ha="center", fontsize=12)
        ax_sw.set_title("STREAM-\nWISENESS", color=FG, fontsize=8,
                        fontfamily="monospace", fontweight="bold", pad=3)
    
    ax_sw.tick_params(colors=FG_DIM, labelsize=7, width=1.2)
    for spine in ax_sw.spines.values():
        spine.set_color(BORDER)
    ax_sw.grid(True, alpha=0.25, color=GRID_CLR)
    
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
    lr03 = f"{params['lr_03']:.1f}" if params.get('lr_03') is not None else "---"
    lr36 = f"{params['lr_36']:.1f}" if params.get('lr_36') is not None else "---"
    
    pwat = fv(params.get("pwat"), "mm", 1) if params.get("pwat") else "---"
    frz = f"{params['frz_level']:.0f}m" if params.get("frz_level") else "---"
    wbo_v = f"{params['wbo']:.0f}m" if params.get('wbo') is not None else "---"
    rh_01 = f"{params.get('rh_0_1km', 0):.0f}%" if params.get('rh_0_1km') else "---"
    rh_13 = f"{params.get('rh_1_3km', 0):.0f}%" if params.get('rh_1_3km') else "---"
    rh_36 = f"{params.get('rh_3_6km', 0):.0f}%" if params.get('rh_3_6km') else "---"
    stp_v = fv(params.get("stp"), "", 1)
    
    # Rows: (text, color, y_position) — tightly spaced
    thermo_rows = [
        ("THERMODYNAMIC", ACCENT, 0.97),
        (header, FG_FAINT, 0.87),
        (f"   {'SB:':8s}  {sb_cape:>10s}  {sb_cin:>10s}  {sb_lcl:>8s}", "#ff8800", 0.77),
        (f"   {'MU:':8s}  {mu_cape:>10s}  {mu_cin:>10s}  {mu_lcl:>8s}", FG, 0.67),
        (f"   {'ML:':8s}  {ml_cape:>10s}  {ml_cin:>10s}  {ml_lcl:>8s}", "#dd44dd", 0.57),
        (f"   DCAPE: {dcape_v}", FG_DIM, 0.47),
        (f"   \u03930-3: {lr03} \u00b0C/km   \u03933-6: {lr36} \u00b0C/km", FG_DIM, 0.37),
        (f"   PWAT: {pwat}  |  FRZ: {frz}  |  WBO: {wbo_v}", FG_DIM, 0.27),
        (f"   RH  0-1km: {rh_01}  1-3km: {rh_13}  3-6km: {rh_36}", FG_DIM, 0.17),
        (f"   STP: {stp_v}", ACCENT, 0.07),
    ]
    
    for text, color, y_pos in thermo_rows:
        ax_params.text(0.01, y_pos, text,
                      transform=ax_params.transAxes,
                      fontsize=10, color=color, fontfamily="monospace",
                      fontweight="bold", va="top")
    
    # ── KINEMATIC TABLE ──
    ax_kin = fig.add_subplot(gs[1, 1:4])
    ax_kin.set_facecolor(BG)
    ax_kin.axis("off")
    
    kin_header = f"{'':3s}{'':8s} {'BWD':>8s}   {'SRH':>12s}"
    
    bwd1 = fv(params.get("bwd_1km"), "kt")
    bwd3 = fv(params.get("bwd_3km"), "kt")
    bwd6 = fv(params.get("bwd_6km"), "kt")
    srh1 = fv(params.get("srh_1km"), "m²/s²")
    srh3 = fv(params.get("srh_3km"), "m²/s²")
    
    kin_rows = [
        ("KINEMATIC", ACCENT, 0.97),
        (kin_header, FG_FAINT, 0.85),
        (f"   {'0-1km:':8s} {bwd1:>8s}   {srh1:>12s}", "#ff3333", 0.73),
        (f"   {'0-3km:':8s} {bwd3:>8s}   {srh3:>12s}", "#ff8800", 0.61),
        (f"   {'0-6km:':8s} {bwd6:>8s}", "#ffcc00", 0.49),
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
            FG_DIM, 0.23
        ))
    
    for text, color, y_pos in kin_rows:
        ax_kin.text(0.01, y_pos, text,
                   transform=ax_kin.transAxes,
                   fontsize=10, color=color, fontfamily="monospace",
                   fontweight="bold", va="top")
    
    # ════════════════════════════════════════════════════════════════
    # MINI MAP INSET (station location) — inside kinematic panel
    # ════════════════════════════════════════════════════════════════
    ax_map = ax_kin.inset_axes([0.55, 0.05, 0.44, 0.90])  # right side of kinematic panel
    ax_map.set_facecolor(BG_PANEL)
    ax_map.plot(_conus_lon, _conus_lat, color=FG_FAINT, linewidth=0.8, zorder=2)
    ax_map.plot(lon, lat, marker="o", color="#ff3333", markersize=6,
                markeredgecolor=FG, markeredgewidth=0.8, zorder=10)
    ax_map.set_xlim(-128, -65)
    ax_map.set_ylim(23, 51)
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
    elif "era5" in stn_name_lower:
        source_label = "ERA5 Reanalysis (ECMWF/CDS)"
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
             fontsize=8.5, color=FG_FAINT, fontfamily="monospace",
             fontweight="bold")
    fig.text(0.96, 0.012,
             f"Station: {station_name}  |  {lat:.4f}, {lon:.4f}",
             fontsize=8.5, color=ACCENT, ha="right", va="bottom",
             fontfamily="monospace", fontweight="bold")
    
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
  era5     ERA5 global reanalysis — any lat/lon, 1940-present (needs CDS key)
  acars    ACARS/AMDAR aircraft obs — airport-based (IEM)

Examples:
  python sounding.py --station OUN                                 # Observed
  python sounding.py --source rap --lat 35.2 --lon -97.4           # RAP at point
  python sounding.py --source bufkit --model hrrr --station OUN    # HRRR forecast
  python sounding.py --source era5 --lat 35.2 --lon -97.4 --date 2020050100
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
                        help="Latitude — required for rap/era5 if no --station")
    parser.add_argument("--lon", type=float, default=None,
                        help="Longitude — required for rap/era5 if no --station")
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
        print("  |  [4] ERA5 reanalysis       (any lat/lon, global)  |")
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
            "4": "era5", "5": "acars",
        }
        source = source_map.get(choice)
        if choice not in source_map:
            print(f"  Invalid choice '{choice}'. Defaulting to observed.")
            source = "obs"

        # Additional interactive prompts depending on source
        if source in ("rap", "era5"):
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
    if source in ("rap", "era5") and lat is not None:
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
        if source == "era5":
            print("  ERA5 tip: ensure CDS API key is configured "
                  "(~/.cdsapirc). See --list-sources.")
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
