"""
Data-fetching functions for every supported source,
plus the unified ``fetch_sounding`` dispatcher.
"""
import warnings
from datetime import datetime, timedelta, timezone
from io import StringIO

import numpy as np
import requests

import metpy.calc as mpcalc
from metpy.units import units

from .constants import (
    STATIONS, STATION_WMO, BUFKIT_MODELS, PSU_MODELS, DATA_SOURCES,
)

warnings.filterwarnings("ignore")


def fetch_iem_sounding(station_id, dt, quiet=False):
    """
    Fetch sounding data from the Iowa Environmental Mesonet (IEM).
    This is generally more reliable than the UWyo server.
    Returns parsed arrays with units.
    """
    # Try the 3-letter station ID first, then WMO number as fallback.
    # Some station IDs (e.g. MFL) are primarily WFO identifiers and
    # IEM may only recognise the WMO number for the radiosonde.
    ids_to_try = [station_id]
    wmo = STATION_WMO.get(station_id)
    if wmo and wmo != station_id:
        ids_to_try.append(wmo)

    profiles = None
    result = None
    last_url = None
    for sid in ids_to_try:
        url = (
            f"https://mesonet.agron.iastate.edu/json/raob.py?"
            f"ts={dt.strftime('%Y%m%d%H%M')}&station={sid}"
        )
        last_url = url
        if not quiet:
            print(f"  Fetching from IEM: {url}")
        resp = requests.get(url, timeout=12)
        resp.raise_for_status()
        result = resp.json()
        profiles = result.get("profiles", [])
        if profiles:
            break

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

    # Try the modern wsgi endpoint first (with required src=FM35),
    # then fall back to the legacy cgi-bin endpoint.
    html = None
    wsgi_url = "https://weather.uwyo.edu/wsgi/sounding"
    wsgi_params = {
        "type": "TEXT:LIST",
        "datetime": dt_str,
        "id": wmo,
        "src": "FM35",
    }
    print(f"  Fetching from UWyo (wsgi): {wsgi_url}  params={wsgi_params}")
    try:
        resp = requests.get(wsgi_url, params=wsgi_params, timeout=20)
        resp.raise_for_status()
        html = resp.text
    except requests.exceptions.RequestException as e:
        print(f"  UWyo wsgi failed ({e}), trying legacy cgi-bin...")
        # Legacy cgi-bin endpoint (widely used by Siphon and others)
        cgi_url = "https://weather.uwyo.edu/cgi-bin/sounding"
        cgi_params = {
            "region": "naconf",
            "TYPE": "TEXT:LIST",
            "YEAR": f"{dt.year}",
            "MONTH": f"{dt.month:02d}",
            "FROM": f"{dt.day:02d}{dt.hour:02d}",
            "TO": f"{dt.day:02d}{dt.hour:02d}",
            "STNM": wmo,
        }
        print(f"  Fetching from UWyo (cgi): {cgi_url}  params={cgi_params}")
        resp = requests.get(cgi_url, params=cgi_params, timeout=20)
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
    
    # UWyo TEXT:LIST has 11 fixed columns:
    # PRES HGHT TEMP DWPT RELH MIXR DRCT SKNT THTA THTE THTV
    # When DRCT+SKNT are missing (common at surface/tropo), the line has
    # 9 tokens instead of 11.  Using split() collapses blank columns, so
    # we must check the token count to avoid reading THTA/THTE as wind.
    pressure, height, temp, dewpoint, wind_dir, wind_spd = [], [], [], [], [], []
    
    for line in lines[data_start:]:
        line = line.strip()
        if not line or line.lower().startswith("</pre>"):
            break
        parts = line.split()
        if len(parts) < 7:       # need at least thermo + theta cols
            continue
        try:
            p  = float(parts[0])
            h  = float(parts[1])
            t  = float(parts[2])
            td = float(parts[3])
            if p < 10 or t < -999 or td < -999:
                continue
            # Wind is only at indices 6,7 when ALL 11 columns are present.
            if len(parts) == 11:
                wd = float(parts[6])
                ws = float(parts[7])
                # Sanity: reject obviously bogus wind
                if ws > 250 or wd > 360:
                    wd, ws = np.nan, np.nan
            else:
                wd, ws = np.nan, np.nan
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
    
    # Interpolate missing wind data (same approach as IEM parser).
    wd_arr = np.array(wind_dir)
    ws_arr = np.array(wind_spd)
    valid_wind = ~np.isnan(wd_arr) & ~np.isnan(ws_arr)
    has_wind = valid_wind.copy()
    if np.any(~valid_wind) and np.sum(valid_wind) >= 2:
        h_arr = np.array(height)
        u_valid, v_valid = mpcalc.wind_components(
            ws_arr[valid_wind] * units.knot,
            wd_arr[valid_wind] * units.degree,
        )
        u_all = np.interp(h_arr, h_arr[valid_wind],
                          u_valid.to("knot").magnitude)
        v_all = np.interp(h_arr, h_arr[valid_wind],
                          v_valid.to("knot").magnitude)
        wd_arr = mpcalc.wind_direction(
            u_all * units.knot, v_all * units.knot
        ).to("degree").magnitude
        ws_arr = mpcalc.wind_speed(
            u_all * units.knot, v_all * units.knot
        ).to("knot").magnitude

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
        "wind_direction": wd_arr * units.degree,
        "wind_speed": ws_arr * units.knot,
        "has_wind": has_wind,
        "station_info": station_info,
    }
    
    return data



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

    stn_lower = station_id.lower()

    # Iowa State BUFKIT archive uses {model}_{station}.buf naming
    url = (
        f"https://mtarchive.geol.iastate.edu/"
        f"{dt:%Y}/{dt:%m}/{dt:%d}/bufkit/{dt:%H}/{model}/{model}_{stn_lower}.buf"
    )

    print(f"  Fetching BUFKIT {model.upper()} f{fhour:03d} for {station_id} "
          f"init {dt:%Y-%m-%d %H}Z...")
    print(f"  URL: {url}")

    resp = requests.get(url, timeout=45)
    if resp.status_code == 404:
        # Some stations use ICAO prefix (k + 3-letter)
        stn_icao = f"k{stn_lower}" if len(stn_lower) == 3 else stn_lower
        url2 = (
            f"https://mtarchive.geol.iastate.edu/"
            f"{dt:%Y}/{dt:%m}/{dt:%d}/bufkit/{dt:%H}/{model}/{model}_{stn_icao}.buf"
        )
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


def fetch_psu_bufkit(station_id, model="rap", fhour=0):
    """
    Fetch latest-run BUFKIT sounding from Penn State's e-wall server.

    PSU provides real-time (latest model cycle) BUFKIT profiles for many
    NWS models.  Unlike the Iowa State archive, only the most recent
    model run is available — no date-based archive.

    URL pattern:
        https://www.meteo.psu.edu/bufkit/data/{MODEL_UPPER}/{CYCLE_HH}/{model_lower}_{station}.buf

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

    model_upper = model.upper()
    model_lower = model.lower()

    # Determine valid cycles for this model and try the most recent ones
    # RAP/HRRR run hourly, NAM/GFS/NAMNEST run every 6h, etc.
    _MODEL_CYCLES = {
        "rap":     list(range(24)),                 # 00-23Z hourly
        "hrrr":    list(range(24)),                 # 00-23Z hourly
        "nam":     [0, 6, 12, 18],                  # 4x daily
        "namnest": [0, 6, 12, 18],
        "gfs":     [0, 6, 12, 18],
        "nam4km":  [0, 6, 12, 18],
        "hiresw":  [0, 6, 12, 18],
        "sref":    [3, 9, 15, 21],
    }
    valid_cycles = _MODEL_CYCLES.get(model_lower, [0, 6, 12, 18])

    # Build candidate cycles: most recent first (based on current UTC hour)
    now_utc = datetime.now(timezone.utc)
    current_hour = now_utc.hour
    # Sort cycles so the most recently passed cycle comes first
    candidates = sorted(valid_cycles, key=lambda c: (current_hour - c) % 24)

    # Try up to 4 most recent cycles
    resp = None
    used_url = None
    for cycle in candidates[:4]:
        cycle_str = f"{cycle:02d}"
        url = (
            f"https://www.meteo.psu.edu/bufkit/data/"
            f"{model_upper}/{cycle_str}/{model_lower}_{stn_k}.buf"
        )
        print(f"  Fetching PSU BUFKIT {model_upper} {cycle_str}Z f{fhour:03d} for {station_id}...")
        print(f"  URL: {url}")
        r = requests.get(url, timeout=45)
        if r.status_code == 200:
            resp = r
            used_url = url
            break
        # Try alternate naming — 3-letter without K prefix
        url2 = url.replace(f"_{stn_k}.buf", f"_{stn_lower}.buf")
        r2 = requests.get(url2, timeout=45)
        if r2.status_code == 200:
            resp = r2
            used_url = url2
            break
        print(f"    → 404, trying next cycle...")

    if resp is None:
        raise ValueError(
            f"PSU BUFKIT data not found for {station_id} ({model_upper}). "
            f"The station may not be available in Penn State's feed."
        )
    resp.raise_for_status()

    print(f"  ✓ Got PSU data from: {used_url}")
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



def fetch_point_sounding(lat, lon, dt, model="gfs", fhour=0):
    """
    Fetch a model-profile point sounding at an arbitrary lat/lon
    from the Open-Meteo pressure-level API.

    Parameters
    ----------
    lat, lon : float
        WGS-84 coordinates.
    dt : datetime
        Model initialisation time (used to pick the correct forecast step).
    model : str
        Model name from our BUFKIT naming (hrrr, rap, gfs, nam, …).
    fhour : int
        Forecast hour offset from init time.
    """
    PRESSURE_LEVELS = [
        1000, 975, 950, 925, 900, 875, 850, 825, 800,
        775, 750, 725, 700, 675, 650, 625, 600, 575,
        550, 525, 500, 475, 450, 425, 400, 375, 350,
        325, 300, 275, 250, 225, 200, 150, 100, 70, 50,
    ]

    # Build variable list for each pressure level
    var_names = []
    for p in PRESSURE_LEVELS:
        var_names.extend([
            f"temperature_{p}hPa",
            f"dew_point_{p}hPa",
            f"wind_speed_{p}hPa",
            f"wind_direction_{p}hPa",
            f"geopotential_height_{p}hPa",
        ])

    # Use gfs_seamless (auto-selects HRRR for CONUS, GFS elsewhere)
    om_model = "gfs_seamless"

    # Compute target valid time
    valid_dt = dt + timedelta(hours=int(fhour))
    start_hour = valid_dt.strftime("%Y-%m-%dT%H:00")
    end_hour = start_hour

    url = "https://api.open-meteo.com/v1/gfs"
    params = {
        "latitude": round(float(lat), 4),
        "longitude": round(float(lon), 4),
        "hourly": ",".join(var_names),
        "models": om_model,
        "wind_speed_unit": "kn",
        "timezone": "GMT",
        "start_hour": start_hour,
        "end_hour": end_hour,
    }

    print(f"  Fetching point sounding at {lat:.2f}, {lon:.2f} "
          f"({model.upper()} F{fhour:03d}) via Open-Meteo...")
    try:
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()
    except requests.exceptions.Timeout:
        raise ValueError(
            "Open-Meteo API timed out. The service may be temporarily "
            "unavailable. Please try again in a few minutes."
        )
    except requests.exceptions.HTTPError as e:
        raise ValueError(
            f"Open-Meteo API error ({e.response.status_code}). "
            f"The service may be temporarily unavailable."
        )
    data = resp.json()

    hourly = data.get("hourly", {})
    times = hourly.get("time", [])
    if not times:
        raise ValueError(
            f"No data returned from Open-Meteo for {lat:.2f}, {lon:.2f} "
            f"at {valid_dt:%Y-%m-%d %HZ}"
        )

    # Use the first (and only) timestep
    idx = 0

    pressure, height, temp, dewpoint, wind_dir, wind_spd = [], [], [], [], [], []

    for p in PRESSURE_LEVELS:
        t_val = hourly.get(f"temperature_{p}hPa", [None])[idx]
        td_val = hourly.get(f"dew_point_{p}hPa", [None])[idx]
        ws_val = hourly.get(f"wind_speed_{p}hPa", [None])[idx]
        wd_val = hourly.get(f"wind_direction_{p}hPa", [None])[idx]
        gh_val = hourly.get(f"geopotential_height_{p}hPa", [None])[idx]

        # Skip levels with missing critical data
        if any(v is None for v in [t_val, td_val, gh_val]):
            continue

        pressure.append(float(p))
        height.append(float(gh_val))
        temp.append(float(t_val))
        dewpoint.append(float(td_val))
        wind_dir.append(float(wd_val) if wd_val is not None else np.nan)
        wind_spd.append(float(ws_val) if ws_val is not None else np.nan)

    if len(pressure) < 5:
        raise ValueError(
            f"Insufficient point sounding data ({len(pressure)} levels) "
            f"from Open-Meteo for {lat:.2f}, {lon:.2f}"
        )

    lat_str = f"{abs(float(lat)):.2f}°{'N' if float(lat) >= 0 else 'S'}"
    lon_str = f"{abs(float(lon)):.2f}°{'W' if float(lon) < 0 else 'E'}"

    return {
        "pressure": np.array(pressure) * units.hPa,
        "height": np.array(height) * units.meter,
        "temperature": np.array(temp) * units.degC,
        "dewpoint": np.array(dewpoint) * units.degC,
        "wind_direction": np.array(wind_dir) * units.degree,
        "wind_speed": np.array(wind_spd) * units.knot,
        "has_wind": np.ones(len(pressure), dtype=bool),
        "station_info": {
            "lat": float(lat),
            "lon": float(lon),
            "name": f"Point Sounding ({lat_str}, {lon_str})",
        },
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
        One of: obs, bufkit, psu.
    lat, lon : float or None
        For point soundings at arbitrary coordinates.
    model : str
        BUFKIT model name (rap, hrrr, nam, …).  Ignored for other sources.
    fhour : int
        BUFKIT forecast hour (0 = analysis).     Ignored for other sources.
    """
    source = source.lower()
    errors = []

    # ── Observed (default) ──────────────────────────────────────────
    if source == "obs":
        # Try the requested time first, then nearby standard cycles
        # (00Z / 12Z) so we can still return data if the exact time
        # is missing from the archive.
        times_to_try = [dt]
        for offset_h in [12, -12, 24, -24]:
            alt = dt + timedelta(hours=offset_h)
            # Only add standard sounding hours (00Z, 12Z)
            if alt.hour in (0, 12) and alt not in times_to_try:
                times_to_try.append(alt)

        for try_dt in times_to_try:
            if try_dt != dt:
                print(f"  Trying nearby cycle: {try_dt}")
            try:
                data = fetch_iem_sounding(station_id, try_dt)
                if try_dt != dt:
                    print(f"  NOTE: Using data from {try_dt} "
                          f"(requested {dt} was unavailable)")
                return data
            except Exception as e:
                errors.append(f"IEM@{try_dt:%Y-%m-%d %HZ}: {e}")
                if not any("UWyo" in err for err in errors):
                    # Only print IEM failure for first attempt
                    if try_dt == dt:
                        print(f"  IEM fetch failed: {e}")
            try:
                data = fetch_wyoming_sounding(station_id, try_dt)
                if try_dt != dt:
                    print(f"  NOTE: Using data from {try_dt} "
                          f"(requested {dt} was unavailable)")
                return data
            except Exception as e:
                errors.append(f"UWyo@{try_dt:%Y-%m-%d %HZ}: {e}")
                if try_dt == dt:
                    print(f"  UWyo fetch failed: {e}")

        # Build a clean, user-friendly error message
        station_name = STATIONS.get(station_id.upper(), (station_id,))[0] if station_id else "Unknown"
        raise ValueError(
            f"No observed sounding data available for {station_id} ({station_name}) "
            f"near {dt:%Y-%m-%d %HZ}. This station may not launch radiosondes regularly, "
            f"or the data hasn't been ingested yet. Try a forecast model (HRRR/RAP/NAM) instead."
        )

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

    raise ValueError(
        f"Unknown source '{source}'. Options: {list(DATA_SOURCES.keys())}"
    )

