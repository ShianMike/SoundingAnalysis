"""
Core sounding routes: fetch, custom upload, file upload.
"""
import base64
import io
import math
import traceback
from datetime import datetime, timezone

from flask import Blueprint, jsonify, request

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sounding import (
    STATIONS, STATION_WMO,
    fetch_sounding, fetch_psu_bufkit, fetch_bufkit_sounding,
    fetch_vad_data,
    compute_parameters, plot_sounding,
    get_latest_sounding_time, find_nearest_station,
)
from .helpers import _fmt, _nan_safe, _serialize_params

bp = Blueprint("sounding_routes", __name__)


@bp.route("/api/sounding", methods=["POST", "OPTIONS"])
def get_sounding():
    if request.method == "OPTIONS":
        return "", 204

    body = request.get_json(force=True)

    source = body.get("source", "obs").lower()
    station = body.get("station")
    date_str = body.get("date")
    lat = body.get("lat")
    lon = body.get("lon")
    model = body.get("model", "rap")
    fhour = body.get("fhour", 0)
    storm_motion = body.get("stormMotion")
    surface_mod = body.get("surfaceMod")
    vad_radar = body.get("vad")
    smoothing = body.get("smoothing")
    sr_hodograph = body.get("srHodograph", False)
    theme = body.get("theme", "dark")
    colorblind = body.get("colorblind", False)
    boundary_orientation = body.get("boundaryOrientation")
    map_zoom = body.get("mapZoom", 1.0)

    # Parse date
    if date_str:
        try:
            dt = datetime.strptime(str(date_str), "%Y%m%d%H").replace(tzinfo=timezone.utc)
        except ValueError:
            return jsonify({"error": f"Invalid date format '{date_str}'. Use YYYYMMDDHH."}), 400
    else:
        dt = get_latest_sounding_time()

    if source == "obs" and station is None and lat is not None and lon is not None:
        station = find_nearest_station(lat, lon)
    if station:
        station = station.upper()

    if source == "obs" and (not station or station not in STATION_WMO):
        return jsonify({"error": f"Unknown station '{station}'. Provide a valid 3-letter ID."}), 400

    if source == "rap" and lat is None and lon is None:
        if station and station in STATIONS:
            _, lat, lon = STATIONS[station]
        else:
            return jsonify({"error": "RAP source requires lat/lon or a known station."}), 400

    try:
        data = fetch_sounding(
            station_id=station, dt=dt, source=source,
            lat=lat, lon=lon, model=model, fhour=fhour,
        )
    except (ValueError, Exception) as fetch_err:
        fallback_data = None
        if source == "bufkit":
            try:
                print(f"  [sounding] BUFKIT failed, trying PSU fallback for {station}...")
                fallback_data = fetch_psu_bufkit(station, model=model, fhour=fhour)
                source = "psu"
            except Exception:
                pass
        elif source == "psu":
            try:
                print(f"  [sounding] PSU failed, trying BUFKIT fallback for {station}...")
                fallback_data = fetch_bufkit_sounding(station, dt, model=model, fhour=fhour)
                source = "bufkit"
            except Exception:
                pass

        if fallback_data is None:
            if isinstance(fetch_err, ValueError):
                return jsonify({"error": str(fetch_err)}), 400
            traceback.print_exc()
            return jsonify({"error": str(fetch_err)}), 500
        data = fallback_data

    try:
        params = compute_parameters(data, storm_motion=storm_motion,
                                    surface_mod=surface_mod, smoothing=smoothing)

        vad_result = None
        if vad_radar:
            try:
                vad_result = fetch_vad_data(vad_radar)
            except Exception as ve:
                print(f"[VAD] Non-fatal error fetching VAD: {ve}")

        plot_id = station or source.upper()
        fig = plot_sounding(data, params, plot_id, dt, vad_data=vad_result,
                            sr_hodograph=sr_hodograph, theme=theme,
                            colorblind=colorblind,
                            boundary_orientation=boundary_orientation,
                            map_zoom=map_zoom)

        _facecolor = "#f5f5f5" if theme == "light" else "#0d0d0d"
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=180, facecolor=_facecolor)
        plt.close(fig)
        buf.seek(0)
        image_b64 = base64.b64encode(buf.read()).decode("utf-8")

        serialized = _serialize_params(params, data, station, dt, source)

        n = len(data["pressure"])
        profile_rows = []
        for i in range(n):
            def _safe_round(val, decimals=1):
                try:
                    v = float(val.magnitude) if hasattr(val, "magnitude") else float(val)
                    return None if (math.isnan(v) or math.isinf(v)) else round(v, decimals)
                except Exception:
                    return None
            row = {
                "p": _safe_round(data["pressure"][i]),
                "h": _safe_round(data["height"][i]),
                "t": _safe_round(data["temperature"][i]),
                "td": _safe_round(data["dewpoint"][i]),
                "wd": _safe_round(data["wind_direction"][i]),
                "ws": _safe_round(data["wind_speed"][i]),
            }
            profile_rows.append(row)

        # Parcel profiles for interactive Skew-T
        def _parcel_array(key):
            prof = params.get(key)
            if prof is None:
                return None
            try:
                arr = prof.to("degC").magnitude if hasattr(prof, "magnitude") else prof
                return [round(float(v), 1) if not (math.isnan(v) or math.isinf(v)) else None for v in arr]
            except Exception:
                return None

        sb_parcel = _parcel_array("sb_profile")
        ml_parcel = _parcel_array("ml_profile")

        stn_info = data.get("station_info", {})
        sfc_elev = stn_info.get("elev", data["height"][0].magnitude if n > 0 else 0)

        return jsonify(_nan_safe({
            "image": image_b64,
            "params": serialized,
            "profile": profile_rows,
            "sbParcel": sb_parcel,
            "mlParcel": ml_parcel,
            "meta": {
                "station": station,
                "stationName": STATIONS.get(station, (station or source.upper(),))[0],
                "source": source,
                "date": dt.strftime("%Y-%m-%d %HZ"),
                "levels": len(data["pressure"]),
                "sfcPressure": round(float(data["pressure"][0].magnitude)),
                "topPressure": round(float(data["pressure"][-1].magnitude)),
                "sfcElevation": round(float(sfc_elev)) if sfc_elev is not None else 0,
                "vadRadar": vad_result.get("radar") if vad_result and isinstance(vad_result, dict) else None,
                "vadTime": vad_result.get("meta", {}).get("time") if vad_result and isinstance(vad_result, dict) else None,
            },
        }))

    except ValueError as e:
        print(f"  [sounding] ValueError: {e}")
        return jsonify({"error": str(e)}), 400
    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


@bp.route("/api/custom-sounding", methods=["POST", "OPTIONS"])
def custom_sounding():
    """Accept raw sounding profile data (text) and return analysis."""
    if request.method == "OPTIONS":
        return "", 204

    import numpy as np
    from metpy.units import units as mpu

    body = request.get_json(force=True)
    raw_text = body.get("text", "").strip()
    fmt = body.get("format", "auto").lower()
    theme = body.get("theme", "dark")
    colorblind = body.get("colorblind", False)

    if not raw_text:
        return jsonify({"error": "No sounding data provided."}), 400

    try:
        lines = [l.strip() for l in raw_text.splitlines() if l.strip()]

        if fmt == "auto":
            first_data = [l for l in lines if not l.startswith("%") and not l.startswith("#")]
            if first_data and len(first_data[0].split()) == 3 and first_data[0].replace(".", "").replace("-", "").replace(" ", "").isdigit():
                fmt = "cm1"
            elif any("," in l for l in first_data[:5]):
                fmt = "csv"
            else:
                fmt = "sharppy"

        p_arr, h_arr, t_arr, td_arr, wd_arr, ws_arr = [], [], [], [], [], []

        if fmt == "csv":
            for line in lines:
                if line.startswith("#") or line.startswith("P") or line.startswith("p"):
                    continue
                parts = [x.strip() for x in line.split(",")]
                if len(parts) < 6:
                    continue
                try:
                    p_arr.append(float(parts[0]))
                    h_arr.append(float(parts[1]))
                    t_arr.append(float(parts[2]))
                    td_arr.append(float(parts[3]))
                    wd_arr.append(float(parts[4]))
                    ws_arr.append(float(parts[5]))
                except ValueError:
                    continue

        elif fmt == "cm1":
            header_parts = lines[0].split()
            sfc_p_hPa = float(header_parts[0])
            sfc_theta = float(header_parts[1])
            sfc_qv = float(header_parts[2]) / 1000.0

            sfc_T_K = sfc_theta * (sfc_p_hPa / 1000.0) ** 0.286
            sfc_T_C = sfc_T_K - 273.15
            e_sfc = sfc_qv * sfc_p_hPa / (0.622 + sfc_qv)
            if e_sfc > 0:
                sfc_Td_C = (243.5 * math.log(e_sfc / 6.112)) / (17.67 - math.log(e_sfc / 6.112))
            else:
                sfc_Td_C = sfc_T_C - 20

            p_arr.append(sfc_p_hPa)
            h_arr.append(0)
            t_arr.append(sfc_T_C)
            td_arr.append(sfc_Td_C)
            wd_arr.append(0)
            ws_arr.append(0)

            for line in lines[1:]:
                parts = line.split()
                if len(parts) < 5:
                    continue
                try:
                    z_m = float(parts[0])
                    theta_K = float(parts[1])
                    qv_kgkg = float(parts[2]) / 1000.0
                    u_ms = float(parts[3])
                    v_ms = float(parts[4])
                    p_est = sfc_p_hPa * math.exp(-z_m / 8500.0)
                    T_K = theta_K * (p_est / 1000.0) ** 0.286
                    T_C = T_K - 273.15
                    e = qv_kgkg * p_est / (0.622 + qv_kgkg)
                    if e > 0:
                        Td_C = (243.5 * math.log(e / 6.112)) / (17.67 - math.log(e / 6.112))
                    else:
                        Td_C = T_C - 20
                    spd = math.sqrt(u_ms**2 + v_ms**2) * 1.94384
                    dirn = (math.degrees(math.atan2(-u_ms, -v_ms)) + 360) % 360
                    p_arr.append(p_est)
                    h_arr.append(z_m)
                    t_arr.append(T_C)
                    td_arr.append(Td_C)
                    wd_arr.append(dirn)
                    ws_arr.append(spd)
                except ValueError:
                    continue
        else:
            for line in lines:
                if line.startswith("%") or line.startswith("#") or line.startswith("SNPARM"):
                    continue
                parts = line.split()
                if len(parts) < 6:
                    continue
                try:
                    vals = [float(x) for x in parts[:6]]
                    if any(v <= -9998 for v in vals):
                        continue
                    p_arr.append(vals[0])
                    h_arr.append(vals[1])
                    t_arr.append(vals[2])
                    td_arr.append(vals[3])
                    wd_arr.append(vals[4])
                    ws_arr.append(vals[5])
                except ValueError:
                    continue

        if len(p_arr) < 5:
            return jsonify({"error": f"Could not parse enough levels ({len(p_arr)} found). "
                           "Need at least 5.  Check format: P(hPa), H(m), T(°C), Td(°C), "
                           "WD(°), WS(kt)."}), 400

        data = {
            "pressure": np.array(p_arr) * mpu.hPa,
            "height": np.array(h_arr) * mpu.meter,
            "temperature": np.array(t_arr) * mpu.degC,
            "dewpoint": np.array(td_arr) * mpu.degC,
            "wind_direction": np.array(wd_arr) * mpu.degree,
            "wind_speed": np.array(ws_arr) * mpu.knot,
            "station_info": {"name": "Custom Upload", "lat": 35.0, "lon": -97.0, "elev": h_arr[0]},
        }

        params = compute_parameters(data)
        fig = plot_sounding(data, params, "CUSTOM", datetime.now(timezone.utc),
                            theme=theme, colorblind=colorblind)

        _facecolor = "#f5f5f5" if theme == "light" else "#0d0d0d"
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=180, facecolor=_facecolor)
        plt.close(fig)
        buf.seek(0)
        image_b64 = base64.b64encode(buf.read()).decode("utf-8")

        serialized = _serialize_params(params, data, "CUSTOM",
                                       datetime.now(timezone.utc), "custom")

        return jsonify({
            "image": image_b64,
            "params": serialized,
            "meta": {
                "station": "CUSTOM",
                "stationName": "Custom Upload",
                "source": "custom",
                "date": datetime.now(timezone.utc).strftime("%Y-%m-%d %HZ"),
                "levels": len(p_arr),
                "sfcPressure": round(p_arr[0]),
                "topPressure": round(p_arr[-1]),
            },
        })

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


@bp.route("/api/upload-file", methods=["POST", "OPTIONS"])
def upload_file():
    """Accept binary file uploads (WRF netCDF, BUFKIT text, etc.)."""
    if request.method == "OPTIONS":
        return "", 204

    import numpy as np
    from metpy.units import units as mpu

    uploaded = request.files.get("file")
    if not uploaded:
        return jsonify({"error": "No file uploaded."}), 400

    fmt = request.form.get("format", "auto").lower()
    theme = request.form.get("theme", "dark")
    colorblind = request.form.get("colorblind", "false").lower() == "true"
    req_lat = request.form.get("lat")
    req_lon = request.form.get("lon")

    filename = uploaded.filename or ""
    file_bytes = uploaded.read()

    if fmt == "auto":
        ext = filename.rsplit(".", 1)[-1].lower() if "." in filename else ""
        if ext in ("nc", "nc4", "ncf", "netcdf") or filename.startswith("wrfout"):
            fmt = "wrf"
        elif ext in ("buf", "bufkit"):
            fmt = "bufkit_file"
        else:
            if file_bytes[:4] in (b'\x89HDF', b'CDF\x01', b'CDF\x02'):
                fmt = "wrf"
            else:
                return jsonify({"error": "Could not auto-detect file format. "
                               "Please select WRF (.nc) or specify the format."}), 400

    try:
        if fmt == "wrf":
            try:
                import netCDF4  # type: ignore[reportMissingImports]
            except ImportError:
                return jsonify({"error": "netCDF4 library is not installed on the server."}), 500

            import tempfile
            with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as tmp:
                tmp.write(file_bytes)
                tmp_path = tmp.name

            try:
                ds = netCDF4.Dataset(tmp_path, "r")
                lats = ds.variables.get("XLAT") or ds.variables.get("XLAT_M")
                lons = ds.variables.get("XLONG") or ds.variables.get("XLONG_M")

                if lats is None or lons is None:
                    ds.close()
                    return jsonify({"error": "WRF file does not contain XLAT/XLONG variables."}), 400

                lat_2d = np.array(lats[0])
                lon_2d = np.array(lons[0])

                if req_lat is not None and req_lon is not None:
                    target_lat, target_lon = float(req_lat), float(req_lon)
                else:
                    target_lat = float(lat_2d[lat_2d.shape[0]//2, lat_2d.shape[1]//2])
                    target_lon = float(lon_2d[lon_2d.shape[0]//2, lon_2d.shape[1]//2])

                dist = (lat_2d - target_lat)**2 + (lon_2d - target_lon)**2
                j, i = np.unravel_index(dist.argmin(), dist.shape)
                t_idx = 0

                if "PB" in ds.variables and "P" in ds.variables:
                    p_pa = np.array(ds.variables["PB"][t_idx,:,j,i]) + np.array(ds.variables["P"][t_idx,:,j,i])
                elif "P" in ds.variables:
                    p_pa = np.array(ds.variables["P"][t_idx,:,j,i])
                else:
                    ds.close()
                    return jsonify({"error": "WRF file missing pressure variables (P, PB)."}), 400
                p_hpa = p_pa / 100.0

                if "PHB" in ds.variables and "PH" in ds.variables:
                    geopot = (np.array(ds.variables["PHB"][t_idx,:,j,i]) + np.array(ds.variables["PH"][t_idx,:,j,i])) / 9.81
                    h_m = 0.5 * (geopot[:-1] + geopot[1:])
                elif "Z" in ds.variables:
                    h_m = np.array(ds.variables["Z"][t_idx,:,j,i])
                else:
                    h_m = 44330.0 * (1.0 - (p_hpa / p_hpa[0])**0.19)

                if "T" in ds.variables:
                    theta_k = np.array(ds.variables["T"][t_idx,:,j,i]) + 300.0
                else:
                    ds.close()
                    return jsonify({"error": "WRF file missing temperature variable (T)."}), 400

                t_k = theta_k * (p_hpa / 1000.0)**0.286
                t_c = t_k - 273.15

                if "QVAPOR" in ds.variables:
                    qv = np.array(ds.variables["QVAPOR"][t_idx,:,j,i])
                    e = qv * p_hpa / (0.622 + qv)
                    e = np.maximum(e, 0.001)
                    td_c = (243.5 * np.log(e / 6.112)) / (17.67 - np.log(e / 6.112))
                else:
                    td_c = t_c - 15.0

                if "U" in ds.variables and "V" in ds.variables:
                    u_stag = np.array(ds.variables["U"][t_idx,:,j,i:i+2])
                    v_stag = np.array(ds.variables["V"][t_idx,:,j:j+2,i])
                    u_ms = 0.5*(u_stag[:,0]+u_stag[:,1]) if u_stag.shape[1]==2 else u_stag[:,0]
                    v_ms = 0.5*(v_stag[:,0]+v_stag[:,1]) if v_stag.shape[1]==2 else v_stag[:,0]
                    ws_kt = np.sqrt(u_ms**2 + v_ms**2) * 1.94384
                    wd_deg = (np.degrees(np.arctan2(-u_ms, -v_ms)) + 360) % 360
                else:
                    ws_kt = np.zeros_like(t_c)
                    wd_deg = np.zeros_like(t_c)

                nz = min(len(p_hpa), len(h_m), len(t_c), len(td_c), len(ws_kt))
                p_hpa, h_m, t_c, td_c, ws_kt, wd_deg = [a[:nz] for a in [p_hpa, h_m, t_c, td_c, ws_kt, wd_deg]]

                grid_lat = float(lat_2d[j, i])
                grid_lon = float(lon_2d[j, i])
                ds.close()

                p_arr, h_arr, t_arr, td_arr, wd_arr, ws_arr = (
                    p_hpa.tolist(), h_m.tolist(), t_c.tolist(), td_c.tolist(), wd_deg.tolist(), ws_kt.tolist()
                )
            finally:
                import os as _os
                try:
                    _os.unlink(tmp_path)
                except OSError:
                    pass
        else:
            return jsonify({"error": f"Unsupported file format: '{fmt}'."}), 400

        if len(p_arr) < 5:
            return jsonify({"error": f"Only {len(p_arr)} levels extracted. Need at least 5."}), 400

        data = {
            "pressure": np.array(p_arr) * mpu.hPa,
            "height": np.array(h_arr) * mpu.meter,
            "temperature": np.array(t_arr) * mpu.degC,
            "dewpoint": np.array(td_arr) * mpu.degC,
            "wind_direction": np.array(wd_arr) * mpu.degree,
            "wind_speed": np.array(ws_arr) * mpu.knot,
            "station_info": {
                "name": f"WRF ({grid_lat:.2f}, {grid_lon:.2f})" if fmt == "wrf" else "File Upload",
                "lat": grid_lat if fmt == "wrf" else 35.0,
                "lon": grid_lon if fmt == "wrf" else -97.0,
                "elev": h_arr[0],
            },
        }

        params = compute_parameters(data)
        label = f"WRF ({grid_lat:.1f}, {grid_lon:.1f})" if fmt == "wrf" else "FILE"
        fig = plot_sounding(data, params, label, datetime.now(timezone.utc),
                            theme=theme, colorblind=colorblind)

        _facecolor = "#f5f5f5" if theme == "light" else "#0d0d0d"
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=180, facecolor=_facecolor)
        plt.close(fig)
        buf.seek(0)
        image_b64 = base64.b64encode(buf.read()).decode("utf-8")

        serialized = _serialize_params(params, data, "WRF" if fmt == "wrf" else "FILE",
                                       datetime.now(timezone.utc), "custom")

        n = len(p_arr)
        profile_rows = [{"p": round(p_arr[k],1), "h": round(h_arr[k],1),
                         "t": round(t_arr[k],1), "td": round(td_arr[k],1),
                         "wd": round(wd_arr[k],1), "ws": round(ws_arr[k],1)} for k in range(n)]

        # Parcel profiles for interactive Skew-T
        try:
            _sp = params.get("sb_profile")
            sb_parcel2 = [round(float(v.magnitude if hasattr(v, 'magnitude') else v), 1) for v in _sp.to("degC")] if _sp is not None else None
        except Exception:
            sb_parcel2 = None
        try:
            _mp = params.get("ml_profile")
            ml_parcel2 = [round(float(v.magnitude if hasattr(v, 'magnitude') else v), 1) for v in _mp.to("degC")] if _mp is not None else None
        except Exception:
            ml_parcel2 = None

        return jsonify(_nan_safe({
            "image": image_b64,
            "params": serialized,
            "profile": profile_rows,
            "sbParcel": sb_parcel2,
            "mlParcel": ml_parcel2,
            "meta": {
                "station": label,
                "stationName": f"WRF Output ({grid_lat:.2f}°N, {grid_lon:.2f}°E)" if fmt == "wrf" else "File Upload",
                "source": "wrf" if fmt == "wrf" else "custom",
                "date": datetime.now(timezone.utc).strftime("%Y-%m-%d %HZ"),
                "levels": len(p_arr),
                "sfcPressure": round(p_arr[0]),
                "topPressure": round(p_arr[-1]),
            },
        }))

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500
