"""
Automated tornado/severe-weather risk scanning across station networks.
"""
import numpy as np

import metpy.calc as mpcalc
from metpy.units import units

from .constants import STATIONS, STATION_WMO, TORNADO_SCAN_STATIONS
from .fetchers import fetch_iem_sounding, fetch_wyoming_sounding
from .parameters import _mixing_ratio_from_dewpoint, compute_parameters


def _quick_tornado_score(station_id, dt):
    """
    Fetch sounding for a station and compute quick severe weather composite scores.
    Returns (stp_score, raw_score, cape, srh, bwd, scp, ship, dcp) or None on failure.
    """
    try:
        # Try IEM with both 3-letter ID and WMO number fallback
        data = None
        ids_to_try = [station_id]
        wmo = STATION_WMO.get(station_id)
        if wmo and wmo != station_id:
            ids_to_try.append(wmo)
        for sid in ids_to_try:
            try:
                data = fetch_iem_sounding(sid, dt, quiet=True)
                break
            except Exception:
                pass
        if data is None:
            # Last resort: try UWyo
            try:
                data = fetch_wyoming_sounding(station_id, dt)
            except Exception:
                return None
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

