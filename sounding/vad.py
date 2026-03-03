"""
NEXRAD VAD Wind Profile (VWP) fetching, parsing, and plotting.
"""
import numpy as np
import requests


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

