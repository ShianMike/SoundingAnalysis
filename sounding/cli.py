"""
Command-line interface for the sounding analysis tool.
"""
import argparse
import sys
from datetime import datetime, timedelta, timezone

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .constants import (
    STATIONS, STATION_WMO, DATA_SOURCES, BUFKIT_MODELS,
)
from .utils import get_latest_sounding_time, find_nearest_station
from .fetchers import fetch_sounding
from .tornado import find_highest_tornado_risk
from .parameters import compute_parameters
from .plotting import plot_sounding


def main():
    parser = argparse.ArgumentParser(
        description="Sounding Analysis Tool — fetch & plot real upper-air data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Data sources (--source):
  obs      Observed radiosonde (IEM / UWyo) — default
  bufkit   BUFKIT forecast models — station-based (Iowa State)
  psu      PSU BUFKIT feed — latest run (Penn State)

Examples:
  python sounding.py --station OUN                                 # Observed
  python sounding.py --source bufkit --model hrrr --station OUN    # HRRR forecast
  python sounding.py --source bufkit --model rap --station OUN     # RAP analysis
""",
    )
    parser.add_argument("--station", type=str, default=None,
                        help="3-letter station ID (e.g. OUN)")
    parser.add_argument("--date", type=str, default=None,
                        help="Date/time as YYYYMMDDHH (e.g. 2024061200). "
                             "Defaults to most recent sounding time.")
    parser.add_argument("--lat", type=float, default=None,
                        help="Latitude (for nearest station lookup)")
    parser.add_argument("--lon", type=float, default=None,
                        help="Longitude (for nearest station lookup)")
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
        print("  |  [2] BUFKIT forecast model (station-based)        |")
        print("  |                                                   |")
        print("  |  [0] Auto-select (tornado risk scan, observed)    |")
        print("  +====================================================+")
        try:
            choice = input("\n  Enter choice [0-2]: ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            return

        source_map = {
            "0": None, "1": "obs", "2": "bufkit",
        }
        source = source_map.get(choice)
        if choice not in source_map:
            print(f"  Invalid choice '{choice}'. Defaulting to observed.")
            source = "obs"

        # Additional interactive prompts depending on source
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
    if station:
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


