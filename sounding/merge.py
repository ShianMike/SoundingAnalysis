"""
Profile merging — weighted blend of two sounding profiles.
"""
import numpy as np

import metpy.calc as mpcalc
from metpy.units import units


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


