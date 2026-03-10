"""
Analysis routes: ensemble plume, compare, composite, merge profiles.
"""
import base64
import io
import os
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta, timezone

import numpy as np
from flask import Blueprint, jsonify, request
from metpy.units import units as mpu

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sounding import (
    STATIONS, STATION_WMO,
    fetch_sounding, fetch_psu_bufkit, fetch_bufkit_sounding,
    compute_parameters, plot_sounding, plot_composite_sounding,
    merge_profiles, get_latest_sounding_time, find_nearest_station,
)
from .helpers import _serialize_params, _nan_safe

bp = Blueprint("analysis", __name__)


# ─── Ensemble Sounding Plume ────────────────────────────────────────
@bp.route("/api/ensemble-plume", methods=["POST", "OPTIONS"])
def ensemble_plume():
    """
    Generate an ensemble sounding plume by fetching multiple forecast hours
    from a BUFKIT model and overlaying them on a single Skew-T.
    """
    if request.method == "OPTIONS":
        return "", 204

    body = request.get_json(force=True) if request.data else {}
    station = body.get("station", "OUN")
    model = body.get("model", "rap")
    src = body.get("source", "psu")
    date_str = body.get("date")
    hours = body.get("hours", [0, 1, 2, 3, 6, 9, 12])
    theme = body.get("theme", "dark")
    cb = body.get("colorblind", False)

    if date_str:
        try:
            dt = datetime.strptime(str(date_str), "%Y%m%d%H").replace(tzinfo=timezone.utc)
        except ValueError:
            return jsonify({"error": f"Invalid date '{date_str}'."}), 400
    else:
        dt = get_latest_sounding_time()

    profiles = []
    errors = []
    used_source = src

    for fh in hours:
        try:
            if src == "psu":
                data = fetch_psu_bufkit(station, model=model, fhour=int(fh))
            else:
                data = fetch_bufkit_sounding(station, dt, model=model, fhour=int(fh))
            profiles.append({"fhour": fh, "data": data})
        except Exception as e:
            errors.append(f"f{fh:03d}: {e}")

    # Auto-fallback: try the other source
    if len(profiles) < 2 and src == "bufkit":
        print(f"  [ensemble-plume] Iowa State failed ({len(profiles)} profiles), trying PSU fallback...")
        profiles.clear()
        errors_fb = []
        for fh in hours:
            try:
                data = fetch_psu_bufkit(station, model=model, fhour=int(fh))
                profiles.append({"fhour": fh, "data": data})
            except Exception as e:
                errors_fb.append(f"f{fh:03d}: {e}")
        if len(profiles) >= 2:
            errors = errors_fb
            used_source = "psu"
        else:
            errors.extend(errors_fb)

    elif len(profiles) < 2 and src == "psu":
        print(f"  [ensemble-plume] PSU failed ({len(profiles)} profiles), trying Iowa State fallback...")
        profiles.clear()
        errors_fb = []
        for fh in hours:
            try:
                data = fetch_bufkit_sounding(station, dt, model=model, fhour=int(fh))
                profiles.append({"fhour": fh, "data": data})
            except Exception as e:
                errors_fb.append(f"f{fh:03d}: {e}")
        if len(profiles) >= 2:
            errors = errors_fb
            used_source = "bufkit"
        else:
            errors.extend(errors_fb)

    # Try previous init cycle (12h earlier) if still not enough
    if len(profiles) < 2 and dt:
        prev_dt = dt - timedelta(hours=12)
        print(f"  [ensemble-plume] Trying previous init cycle: {prev_dt:%Y-%m-%d %H}Z...")
        for try_src_label in ("psu", "bufkit"):
            profiles.clear()
            errors_prev = []
            for fh in hours:
                try:
                    if try_src_label == "psu":
                        data = fetch_psu_bufkit(station, model=model, fhour=int(fh))
                    else:
                        data = fetch_bufkit_sounding(station, prev_dt, model=model, fhour=int(fh))
                    profiles.append({"fhour": fh, "data": data})
                except Exception as e:
                    errors_prev.append(f"f{fh:03d}: {e}")
            if len(profiles) >= 2:
                errors = errors_prev
                dt = prev_dt
                used_source = try_src_label
                break

    if len(profiles) < 2:
        unique_reasons = []
        seen = set()
        for e in errors:
            reason = e.split(": ", 1)[1] if ": " in e else e
            short = reason.split(".")[0].strip()
            if short not in seen:
                seen.add(short)
                unique_reasons.append(short)

        tried_sources = []
        if src == "bufkit" or used_source == "bufkit":
            tried_sources.append("Iowa State")
        if src == "psu" or used_source == "psu":
            tried_sources.append("Penn State")
        tried_str = " and ".join(tried_sources) if tried_sources else src

        summary = (
            f"Could not load forecast data for {station.upper()} ({model.upper()}).\n"
            f"Tried {len(hours)} forecast hours from {tried_str} — none returned data."
        )
        if unique_reasons:
            summary += f"\n\nReason: {unique_reasons[0]}."

        suggestions = []
        if model == "sref":
            suggestions.append("SREF was discontinued in 2025 — switch to RAP, HRRR, or NAM.")
        if src == "bufkit":
            suggestions.append('Switch source to "Penn State (latest)" — it has the most recent model runs.')
        elif src == "psu":
            suggestions.append('Switch source to "Iowa State (archive)" for historical data.')
        suggestions.append("Try a different model (RAP and HRRR have the best station coverage).")
        suggestions.append("Try a different station — not all sites are available in all model feeds.")

        return jsonify({"error": summary, "suggestions": suggestions}), 400

    base_data = profiles[0]["data"]
    base_params = compute_parameters(base_data)

    try:
        from matplotlib.gridspec import GridSpec
        from metpy.plots import SkewT

        bg = "#f5f5f5" if theme == "light" else "#0d0d0d"
        fg = "#333333" if theme == "light" else "#e2e8f0"
        grid_c = "#cccccc" if theme == "light" else "#2a2a2a"

        fig = plt.figure(figsize=(12, 10), facecolor=bg)
        gs = GridSpec(1, 2, width_ratios=[3, 1], figure=fig)

        ax_skew = fig.add_subplot(gs[0, 0])
        skew = SkewT(fig, rotation=45, subplot=ax_skew)
        skew.ax.set_facecolor(bg)
        skew.ax.tick_params(colors=fg, labelsize=8)
        for spine in skew.ax.spines.values():
            spine.set_color(grid_c)
        skew.ax.set_ylim(1050, 100)
        skew.ax.set_xlim(-40, 50)

        n = len(profiles)
        if cb:
            colors = plt.cm.viridis(np.linspace(0.2, 0.9, n))
        else:
            colors = plt.cm.plasma(np.linspace(0.15, 0.85, n))

        # Interpolate each member to a common pressure grid for percentile bands.
        p_grid = np.arange(1000, 99, -25, dtype=float)  # hPa (surface to upper levels)
        t_members = []
        td_members = []

        def _interp_to_grid(d, field_name):
            p = np.asarray(d["pressure"].m, dtype=float)
            x = np.asarray(d[field_name].m, dtype=float)
            mask = np.isfinite(p) & np.isfinite(x)
            if mask.sum() < 3:
                return None
            p_valid = p[mask]
            x_valid = x[mask]
            order = np.argsort(p_valid)
            p_sorted = p_valid[order]
            x_sorted = x_valid[order]
            p_unique, idx = np.unique(p_sorted, return_index=True)
            x_unique = x_sorted[idx]
            if p_unique.size < 3:
                return None
            out = np.full_like(p_grid, np.nan, dtype=float)
            valid = (p_grid >= p_unique.min()) & (p_grid <= p_unique.max())
            if np.any(valid):
                out[valid] = np.interp(p_grid[valid], p_unique, x_unique)
            return out

        for i, pf in enumerate(profiles):
            d = pf["data"]
            p = d["pressure"].m
            t = d["temperature"].m
            td = d["dewpoint"].m
            c = colors[i]
            alpha = 0.3 if n > 5 else 0.5
            skew.plot(p, t, c, linewidth=1.0, alpha=alpha)
            skew.plot(p, td, c, linewidth=0.8, alpha=alpha * 0.7, linestyle="--")

            t_interp = _interp_to_grid(d, "temperature")
            td_interp = _interp_to_grid(d, "dewpoint")
            if t_interp is not None:
                t_members.append(t_interp)
            if td_interp is not None:
                td_members.append(td_interp)

        # Percentile plume (10/25/50/75/90) across all available members.
        if len(t_members) >= 3:
            t_arr = np.asarray(t_members)
            t10, t25, t50, t75, t90 = np.nanpercentile(t_arr, [10, 25, 50, 75, 90], axis=0)
            valid_wide = np.isfinite(t10) & np.isfinite(t90)
            valid_core = np.isfinite(t25) & np.isfinite(t75)
            valid_med = np.isfinite(t50)
            if np.any(valid_wide):
                skew.ax.fill_betweenx(
                    p_grid[valid_wide], t10[valid_wide], t90[valid_wide],
                    color="#ef4444", alpha=0.08, linewidth=0, label="T 10–90%",
                )
            if np.any(valid_core):
                skew.ax.fill_betweenx(
                    p_grid[valid_core], t25[valid_core], t75[valid_core],
                    color="#ef4444", alpha=0.16, linewidth=0, label="T 25–75%",
                )
            if np.any(valid_med):
                skew.plot(p_grid[valid_med], t50[valid_med], color="#dc2626", linewidth=2.0, alpha=0.9, label="T median")

        if len(td_members) >= 3:
            td_arr = np.asarray(td_members)
            td10, td25, td50, td75, td90 = np.nanpercentile(td_arr, [10, 25, 50, 75, 90], axis=0)
            valid_wide = np.isfinite(td10) & np.isfinite(td90)
            valid_core = np.isfinite(td25) & np.isfinite(td75)
            valid_med = np.isfinite(td50)
            if np.any(valid_wide):
                skew.ax.fill_betweenx(
                    p_grid[valid_wide], td10[valid_wide], td90[valid_wide],
                    color="#22c55e", alpha=0.06, linewidth=0, label="Td 10–90%",
                )
            if np.any(valid_core):
                skew.ax.fill_betweenx(
                    p_grid[valid_core], td25[valid_core], td75[valid_core],
                    color="#22c55e", alpha=0.12, linewidth=0, label="Td 25–75%",
                )
            if np.any(valid_med):
                skew.plot(
                    p_grid[valid_med], td50[valid_med],
                    color="#16a34a", linewidth=1.6, alpha=0.9, linestyle="--", label="Td median",
                )

        d0 = profiles[0]["data"]
        skew.plot(d0["pressure"].m, d0["temperature"].m, "r", linewidth=2, alpha=0.9, label="Analysis (f000)")
        skew.plot(d0["pressure"].m, d0["dewpoint"].m, "g", linewidth=1.5, alpha=0.8)

        handles, labels = skew.ax.get_legend_handles_labels()
        if handles:
            dedup = {}
            for h, lbl in zip(handles, labels):
                if lbl and lbl not in dedup:
                    dedup[lbl] = h
            skew.ax.legend(
                dedup.values(), dedup.keys(),
                loc="lower left", fontsize=7, framealpha=0.22,
                facecolor=bg, edgecolor=grid_c,
            )

        skew.ax.set_title(
            f"Ensemble Plume: {station.upper()} {model.upper()}\n"
            f"{n} members (f{min(hours):03d} – f{max(hours):03d})",
            fontsize=12, color=fg, pad=10,
        )

        ax_hodo = fig.add_subplot(gs[0, 1])
        ax_hodo.set_facecolor(bg)
        ax_hodo.set_aspect("equal")
        ax_hodo.tick_params(colors=fg, labelsize=7)
        ax_hodo.set_title("Hodograph Spread", fontsize=10, color=fg)

        for i, pf in enumerate(profiles):
            d = pf["data"]
            try:
                u_w = -d["wind_speed"].m * np.sin(np.radians(d["wind_direction"].m)) * 0.514444
                v_w = -d["wind_speed"].m * np.cos(np.radians(d["wind_direction"].m)) * 0.514444
                h = d["height"].m
                mask = h <= h[0] + 10000
                ax_hodo.plot(u_w[mask], v_w[mask], color=colors[i], alpha=0.3, linewidth=0.8)
            except Exception:
                pass

        ax_hodo.axhline(0, color=grid_c, linewidth=0.5)
        ax_hodo.axvline(0, color=grid_c, linewidth=0.5)
        ax_hodo.set_xlabel("u (m/s)", fontsize=8, color=fg)
        ax_hodo.set_ylabel("v (m/s)", fontsize=8, color=fg)

        fig.tight_layout()

        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=180, facecolor=bg)
        plt.close(fig)
        buf.seek(0)
        image_b64 = base64.b64encode(buf.read()).decode("utf-8")

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": f"Plume plot failed: {e}"}), 500

    serialized = _serialize_params(base_params, base_data, station, dt, "bufkit")

    return jsonify({
        "image": image_b64,
        "params": serialized,
        "members": len(profiles),
        "hours": [p["fhour"] for p in profiles],
        "errors": errors if errors else None,
        "meta": {
            "station": station.upper(),
            "model": model.upper(),
            "source": used_source,
            "date": dt.strftime("%Y-%m-%d %HZ") if dt else "latest",
        },
    })


# ─── Forecast Profiles (raw data for animation) ─────────────────────
@bp.route("/api/forecast-profiles", methods=["POST", "OPTIONS"])
def forecast_profiles():
    """
    Return raw profile + params data for multiple BUFKIT forecast hours.
    Used by the client-side sounding animation feature.
    """
    if request.method == "OPTIONS":
        return "", 204

    import math
    body = request.get_json(force=True) if request.data else {}
    station = body.get("station", "OUN")
    model = body.get("model", "rap")
    src = body.get("source", "psu")
    date_str = body.get("date")
    hours = body.get("hours", [0, 1, 2, 3, 6, 9, 12])

    if date_str:
        try:
            dt = datetime.strptime(str(date_str), "%Y%m%d%H").replace(tzinfo=timezone.utc)
        except ValueError:
            return jsonify({"error": f"Invalid date '{date_str}'."}), 400
    else:
        dt = get_latest_sounding_time()

    results = []
    errors = []

    def _safe_round(val, decimals=1):
        try:
            v = float(val.magnitude) if hasattr(val, "magnitude") else float(val)
            return None if (math.isnan(v) or math.isinf(v)) else round(v, decimals)
        except Exception:
            return None

    for fh in hours:
        try:
            if src == "psu":
                data = fetch_psu_bufkit(station, model=model, fhour=int(fh))
            else:
                data = fetch_bufkit_sounding(station, dt, model=model, fhour=int(fh))

            params = compute_parameters(data)
            n = len(data["pressure"])
            profile_rows = []
            for i in range(n):
                profile_rows.append({
                    "p": _safe_round(data["pressure"][i]),
                    "h": _safe_round(data["height"][i]),
                    "t": _safe_round(data["temperature"][i]),
                    "td": _safe_round(data["dewpoint"][i]),
                    "wd": _safe_round(data["wind_direction"][i]),
                    "ws": _safe_round(data["wind_speed"][i]),
                })

            # Parcel profiles
            def _parcel_arr(key):
                prof = params.get(key)
                if prof is None:
                    return None
                try:
                    arr = prof.to("degC").magnitude if hasattr(prof, "magnitude") else prof
                    return [round(float(v), 1) if not (math.isnan(v) or math.isinf(v)) else None for v in arr]
                except Exception:
                    return None

            ser = _serialize_params(params, data, station, dt, "bufkit")
            results.append({
                "fhour": fh,
                "profile": profile_rows,
                "sbParcel": _parcel_arr("sb_profile"),
                "mlParcel": _parcel_arr("ml_profile"),
                "params": ser,
            })
        except Exception as e:
            errors.append(f"f{fh:03d}: {e}")

    if len(results) == 0:
        return jsonify({"error": "No forecast profiles could be loaded.", "errors": errors}), 400

    return _nan_safe(jsonify({
        "frames": results,
        "errors": errors if errors else None,
        "meta": {
            "station": station.upper(),
            "model": model.upper(),
            "source": src,
            "date": dt.strftime("%Y-%m-%d %HZ") if dt else "latest",
        },
    }))


# ── Compare endpoint ────────────────────────────────────────────────
@bp.route("/api/compare", methods=["POST", "OPTIONS"])
def compare_soundings():
    """Fetch multiple soundings in parallel for side-by-side comparison."""
    if request.method == "OPTIONS":
        return "", 204

    body = request.get_json(force=True)
    items = body.get("soundings", [])

    if not items or not isinstance(items, list):
        return jsonify({"error": "Provide a 'soundings' array."}), 400
    if len(items) > 4:
        return jsonify({"error": "Maximum 4 soundings per comparison."}), 400

    compare_workers = int(os.environ.get("COMPARE_WORKERS", "2"))

    def _fetch_one_compare(item):
        src = item.get("source", "obs").lower()
        stn = item.get("station")
        date_str = item.get("date")
        lat_v = item.get("lat")
        lon_v = item.get("lon")
        mdl = item.get("model", "rap")
        fhr = item.get("fhour", 0)

        if date_str:
            dt = datetime.strptime(str(date_str), "%Y%m%d%H").replace(tzinfo=timezone.utc)
        else:
            dt = get_latest_sounding_time()

        if src == "obs" and stn is None and lat_v is not None and lon_v is not None:
            stn = find_nearest_station(lat_v, lon_v)
        if stn:
            stn = stn.upper()
        if src in ("rap", "era5") and lat_v is None and lon_v is None:
            if stn and stn in STATIONS:
                _, lat_v, lon_v = STATIONS[stn]

        data = fetch_sounding(
            station_id=stn, dt=dt, source=src,
            lat=lat_v, lon=lon_v, model=mdl, fhour=fhr,
        )
        params = compute_parameters(data)

        plot_id = stn or src.upper()
        fig = plot_sounding(data, params, plot_id, dt)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=180, facecolor="#0d0d0d")
        plt.close(fig)
        buf.seek(0)
        image_b64 = base64.b64encode(buf.read()).decode("utf-8")

        serialized = _serialize_params(params, data, stn, dt, src)

        return {
            "image": image_b64,
            "params": serialized,
            "meta": {
                "station": stn,
                "stationName": STATIONS.get(stn, (stn or src.upper(),))[0],
                "source": src,
                "date": dt.strftime("%Y-%m-%d %HZ"),
                "levels": len(data["pressure"]),
                "sfcPressure": round(float(data["pressure"][0].magnitude)),
                "topPressure": round(float(data["pressure"][-1].magnitude)),
            },
        }

    try:
        results = []
        with ThreadPoolExecutor(max_workers=compare_workers) as executor:
            futures = {executor.submit(_fetch_one_compare, item): i for i, item in enumerate(items)}
            for future in as_completed(futures):
                idx = futures[future]
                try:
                    results.append((idx, future.result()))
                except Exception as e:
                    results.append((idx, {"error": str(e), "meta": items[idx]}))

        results.sort(key=lambda x: x[0])
        return jsonify({"soundings": [r[1] for r in results]})

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


# ── Composite endpoint ──────────────────────────────────────────────
@bp.route("/api/composite", methods=["POST", "OPTIONS"])
def composite_sounding():
    """Generate a composite overlay plot with multiple soundings on one Skew-T."""
    if request.method == "OPTIONS":
        return "", 204

    body = request.get_json(force=True)
    items = body.get("soundings", [])

    if not items or not isinstance(items, list):
        return jsonify({"error": "Provide a 'soundings' array."}), 400
    if len(items) > 6:
        return jsonify({"error": "Maximum 6 soundings per composite."}), 400

    comp_profiles = []
    comp_errors = []
    for item in items:
        try:
            src = item.get("source", "obs").lower()
            stn = item.get("station")
            date_str = item.get("date")
            lat_v = item.get("lat")
            lon_v = item.get("lon")
            mdl = item.get("model", "rap")
            fhr = item.get("fhour", 0)

            if date_str:
                dt = datetime.strptime(str(date_str), "%Y%m%d%H").replace(tzinfo=timezone.utc)
            else:
                dt = get_latest_sounding_time()

            if src == "obs" and stn is None and lat_v is not None and lon_v is not None:
                stn = find_nearest_station(lat_v, lon_v)
            if stn:
                stn = stn.upper()
            if src in ("rap", "era5") and lat_v is None and lon_v is None:
                if stn and stn in STATIONS:
                    _, lat_v, lon_v = STATIONS[stn]

            data = fetch_sounding(
                station_id=stn, dt=dt, source=src,
                lat=lat_v, lon=lon_v, model=mdl, fhour=fhr,
            )
            params = compute_parameters(data)

            label = f"{stn or src.upper()} {dt.strftime('%d/%HZ')}"
            comp_profiles.append({"data": data, "params": params, "label": label})
        except Exception as ex:
            comp_errors.append({"item": item, "error": str(ex)})

    if not comp_profiles:
        return jsonify({"error": "All soundings failed to fetch.", "details": comp_errors}), 500

    try:
        fig = plot_composite_sounding(comp_profiles, title="Composite Sounding Overlay")
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=180, facecolor="#0d0d0d")
        plt.close(fig)
        buf.seek(0)
        image_b64 = base64.b64encode(buf.read()).decode("utf-8")

        return jsonify({
            "image": image_b64,
            "count": len(comp_profiles),
            "errors": comp_errors if comp_errors else None,
        })

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


# ── Merge Profiles endpoint ─────────────────────────────────────────
@bp.route("/api/merge-profiles", methods=["POST", "OPTIONS"])
def merge_profiles_endpoint():
    """Merge two soundings into a weighted-average blended profile."""
    if request.method == "OPTIONS":
        return "", 204

    body = request.get_json(force=True)
    items = body.get("soundings", [])
    weight = float(body.get("weight", 0.5))
    theme = body.get("theme", "dark")
    cb = body.get("colorblind", False)

    if len(items) != 2:
        return jsonify({"error": "Exactly 2 soundings required for merging."}), 400

    try:
        datasets = []
        labels = []
        for item in items:
            src = item.get("source", "obs").lower()
            stn = item.get("station")
            date_str = item.get("date")
            lat_v = item.get("lat")
            lon_v = item.get("lon")
            mdl = item.get("model", "rap")
            fhr = item.get("fhour", 0)

            if date_str:
                dt = datetime.strptime(str(date_str), "%Y%m%d%H").replace(tzinfo=timezone.utc)
            else:
                dt = get_latest_sounding_time()

            if src == "obs" and stn is None and lat_v is not None and lon_v is not None:
                stn = find_nearest_station(lat_v, lon_v)
            if stn:
                stn = stn.upper()
            if src in ("rap", "era5") and lat_v is None and lon_v is None:
                if stn and stn in STATIONS:
                    _, lat_v, lon_v = STATIONS[stn]

            data = fetch_sounding(
                station_id=stn, dt=dt, source=src,
                lat=lat_v, lon=lon_v, model=mdl, fhour=fhr,
            )
            datasets.append(data)
            labels.append(f"{stn or src.upper()} {dt.strftime('%d/%HZ')}")

        merged_data = merge_profiles(datasets[0], datasets[1], weight_a=weight)
        merged_params = compute_parameters(merged_data)

        merge_label = f"MERGE: {int(weight*100)}% {labels[0]} + {int((1-weight)*100)}% {labels[1]}"
        fig = plot_sounding(merged_data, merged_params, "MERGED", dt,
                            theme=theme, colorblind=cb)

        _facecolor = "#f5f5f5" if theme == "light" else "#0d0d0d"
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=180, facecolor=_facecolor)
        plt.close(fig)
        buf.seek(0)
        image_b64 = base64.b64encode(buf.read()).decode("utf-8")

        serialized = _serialize_params(merged_params, merged_data, "MERGED", dt, "merge")

        return jsonify({
            "image": image_b64,
            "params": serialized,
            "label": merge_label,
            "meta": {
                "station": "MERGED",
                "stationName": merge_label,
                "source": "merge",
                "date": dt.strftime("%Y-%m-%d %HZ"),
                "levels": len(merged_data["pressure"]),
                "weight": weight,
                "profileA": labels[0],
                "profileB": labels[1],
            },
        })

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500
