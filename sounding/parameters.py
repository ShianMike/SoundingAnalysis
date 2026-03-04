"""
Comprehensive thermodynamic & kinematic parameter computation.
"""
import warnings

import numpy as np
from scipy.ndimage import gaussian_filter1d

import metpy.calc as mpcalc
from metpy.units import units
import metpy.constants as mpconst

warnings.filterwarnings("ignore")


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

    # ── Piecewise CAPE (50 hPa layers from LFC to EL) ──────────────
    try:
        _pw_layers = []
        _mu_cape_pw = float(params.get("mu_cape", 0 * units("J/kg")).magnitude)
        _mu_si_pw = params.get("mu_start_idx", 0) or 0
        if _mu_cape_pw > 0 and params.get("mu_profile") is not None:
            # Slice environment arrays to match the MU profile (starts at mu_start_idx)
            _T_env_pw = T[_mu_si_pw:].magnitude
            _T_parcel_pw = params["mu_profile"].magnitude
            _p_arr_pw = p[_mu_si_pw:].magnitude
            _h_arr_pw = h[_mu_si_pw:].magnitude
            # Define layers by pressure: each ~50 hPa thick from 900 to 200 hPa
            _layer_edges = list(range(900, 150, -50))
            for i in range(len(_layer_edges) - 1):
                _p_top = _layer_edges[i + 1]
                _p_bot = _layer_edges[i]
                _mask_pw = (_p_arr_pw <= _p_bot) & (_p_arr_pw >= _p_top)
                if np.sum(_mask_pw) < 2:
                    continue
                _buoy = _T_parcel_pw[_mask_pw] - _T_env_pw[_mask_pw]
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

    # ── Bulk Richardson Number (BRN) & Convective Mode ──────────────────
    try:
        _brn_cape = float(params.get("sb_cape", 0 * units("J/kg")).magnitude)
        _brn_shear = float(params.get("bwd_6km", 0 * units.knot).to("m/s").magnitude)
        if _brn_shear > 0 and _brn_cape > 0:
            _brn = _brn_cape / (0.5 * _brn_shear ** 2)
        else:
            _brn = None
        params["brn"] = round(_brn, 1) if _brn is not None else None

        # Convective mode estimate (Thompson et al. 2007 framework)
        _bwd6_kt = float(params.get("bwd_6km", 0 * units.knot).magnitude)
        _scp_val = float(params.get("scp", 0))
        _srh1_val = float(params.get("srh_1km", 0 * units("m**2/s**2")).magnitude)

        if _brn is not None and _brn < 10 and _bwd6_kt > 50:
            _conv_mode = "Discrete Supercell"
        elif _brn is not None and 10 <= _brn < 45 and _bwd6_kt >= 35:
            _conv_mode = "Discrete / Supercell"
        elif _bwd6_kt >= 30 and _scp_val >= 1:
            _conv_mode = "Supercell likely"
        elif _bwd6_kt >= 25 and _srh1_val >= 100:
            _conv_mode = "Rotating Storms"
        elif _bwd6_kt >= 20:
            _conv_mode = "Multicell / Clusters"
        elif _bwd6_kt >= 10:
            _conv_mode = "Weak Multicell"
        else:
            _conv_mode = "Single Cell / Pulse"

        params["convective_mode"] = _conv_mode
    except Exception as e:
        print(f"  Warning: BRN / convective mode calc failed: {e}")
        params["brn"] = None
        params["convective_mode"] = None

    return params


