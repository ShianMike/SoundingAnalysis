"""
Sounding analysis plot generation (Skew-T, hodograph, parameter tables).
"""
import warnings
from datetime import datetime, timezone

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as path_effects
import numpy as np

import metpy.calc as mpcalc
from metpy.plots import SkewT, Hodograph
from metpy.units import units

from .constants import STATIONS

warnings.filterwarnings("ignore")


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
        CLR_DEW       = "#22dd22"
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
        height_ratios=[2.7, 1.3],
        hspace=0.10, wspace=0.12,
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
    skew.plot(p, T, color=CLR_TEMP, linewidth=4.0, zorder=6, label="TEMPERATURE")
    # Dewpoint (green, solid thick)
    skew.plot(p, Td, color=CLR_DEW, linewidth=4.0, zorder=6, label="DEWPOINT")
    
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
    skew.ax.set_xlim(-50, 50)
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
    
    # Rows: (text, color, y_position) — evenly spaced for readability
    thermo_rows = [
        ("THERMODYNAMIC", ACCENT, 0.97),
        (header, FG_FAINT, 0.90),
        (f"   {'SB:':8s}  {sb_cape:>10s}  {sb_cin:>10s}  {sb_lcl:>8s}", CLR_SB_PARCEL, 0.84),
        (f"   {'MU:':8s}  {mu_cape:>10s}  {mu_cin:>10s}  {mu_lcl:>8s}", CLR_MU_PARCEL, 0.77),
        (f"   {'ML:':8s}  {ml_cape:>10s}  {ml_cin:>10s}  {ml_lcl:>8s}", CLR_ML_PARCEL, 0.70),
        (f"   ML LFC: {ml_lfc_v}  |  ML EL: {ml_el_v}  |  WCD: {wcd_v}", FG_DIM, 0.63),
        (f"   DCAPE: {dcape_v}  |  DCIN: {dcin_v}  |  ECAPE: {ecape_v}", FG_DIM, 0.57),
        (f"   3CAPE: {cape3_v}  |  6CAPE: {cape6_v}  |  NCAPE: {ncape_v}", FG_DIM, 0.50),
        (f"   \u03930-3: {lr03} \u00b0C/km   \u03933-6: {lr36} \u00b0C/km", FG_DIM, 0.43),
        (f"   PWAT: {pwat}  |  FRZ: {frz}  |  WBO: {wbo_v}", FG_DIM, 0.37),
        (f"   RH  0-1km: {rh_01}  1-3km: {rh_13}  3-6km: {rh_36}", FG_DIM, 0.30),
        (f"   Precip Type: {precip_type_v}", "#22d3ee", 0.23),
        (f"   STP: {stp_v}  STPeff: {stp_eff_v}  |  SCP: {fv(params.get('scp'), '', 1)}  |  SHIP: {fv(params.get('ship'), '', 1)}", ACCENT, 0.17),
        (f"   DCP: {fv(params.get('dcp'), '', 1)}  |  FireWx: FWI {fosberg_v}  Haines {haines_v}", ACCENT, 0.10),
        (f"   {'[SURFACE MODIFIED]' if params.get('surface_modified') else ''}{'  [CUSTOM SM]' if params.get('custom_storm_motion') else ''}", "#ff5555" if params.get('surface_modified') or params.get('custom_storm_motion') else BG, 0.03),
    ]
    
    for text, color, y_pos in thermo_rows:
        ax_params.text(0.01, y_pos, text,
                      transform=ax_params.transAxes,
                      fontsize=11, color=color, fontfamily="monospace",
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
                   fontsize=11, color=color, fontfamily="monospace",
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

