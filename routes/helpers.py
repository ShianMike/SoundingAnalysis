"""
Shared helpers used by multiple route modules.
"""
import math


def _fmt(val, unit_str="", decimals=0):
    """Format a pint-style value, returning None for missing/NaN data."""
    if val is None:
        return None
    try:
        v = val.magnitude if hasattr(val, "magnitude") else val
        fv = float(v)
        if math.isnan(fv) or math.isinf(fv):
            return None
        if decimals == 0:
            return round(fv)
        return round(fv, decimals)
    except Exception:
        return None


def _nan_safe(obj):
    """Recursively replace NaN/Inf float values with None for JSON safety."""
    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return None
        return obj
    if isinstance(obj, dict):
        return {k: _nan_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_nan_safe(item) for item in obj]
    return obj


def _serialize_params(params, data, station, dt, source):
    """Extract key computed parameters into a JSON-friendly dict."""
    return {
        # Thermodynamic
        "sbCape": _fmt(params.get("sb_cape")),
        "sbCin": _fmt(params.get("sb_cin")),
        "sbLclM": round(params["sb_lcl_m"]) if params.get("sb_lcl_m") is not None else None,
        "sbLclP": _fmt(params.get("sb_lcl_p")),
        "sbLfcP": _fmt(params.get("sb_lfc_p")),
        "sbElP": _fmt(params.get("sb_el_p")),
        "muCape": _fmt(params.get("mu_cape")),
        "muCin": _fmt(params.get("mu_cin")),
        "muLclM": round(params["mu_lcl_m"]) if params.get("mu_lcl_m") is not None else None,
        "mlCape": _fmt(params.get("ml_cape")),
        "mlCin": _fmt(params.get("ml_cin")),
        "mlLclM": round(params["ml_lcl_m"]) if params.get("ml_lcl_m") is not None else None,
        "dcape": _fmt(params.get("dcape")),
        "lr03": round(params["lr_03"], 1) if params.get("lr_03") is not None else None,
        "lr36": round(params["lr_36"], 1) if params.get("lr_36") is not None else None,
        "pwat": _fmt(params.get("pwat"), decimals=1),
        "frzLevel": round(params["frz_level"]) if params.get("frz_level") is not None else None,
        "wbo": round(params["wbo"]) if params.get("wbo") is not None else None,
        "stp": round(params.get("stp", 0), 1),
        "stpEff": round(params.get("stp_eff", 0), 1),
        "scp": round(params.get("scp", 0), 1),
        "ship": round(params.get("ship", 0), 1),
        "dcp": round(params.get("dcp", 0), 1),
        "ecape": round(params.get("ecape", 0), 1),
        "piecewiseCape": params.get("piecewise_cape", []),
        "cape3km": params.get("cape_3km", 0),
        "cape6km": params.get("cape_6km", 0),
        "dcin": params.get("dcin", 0),
        "ncape": params.get("ncape", 0),
        "mlLfcM": round(params["ml_lfc_m"]) if params.get("ml_lfc_m") is not None else None,
        "mlLfcP": _fmt(params.get("ml_lfc_p")),
        "mlElM": round(params["ml_el_m"]) if params.get("ml_el_m") is not None else None,
        "mlElP": _fmt(params.get("ml_el_p")),
        "wcd": params.get("wcd"),
        "surfaceModified": params.get("surface_modified", False),
        "customStormMotion": params.get("custom_storm_motion", False),
        "smoothingApplied": params.get("smoothing_applied", False),
        "rh01": round(params["rh_0_1km"]) if params.get("rh_0_1km") is not None else None,
        "rh13": round(params["rh_1_3km"]) if params.get("rh_1_3km") is not None else None,
        "rh36": round(params["rh_3_6km"]) if params.get("rh_3_6km") is not None else None,
        # Winter weather / precip type
        "precipType": params.get("precip_type", "N/A"),
        "warmLayerEnergy": params.get("warm_layer_energy", 0),
        "coldLayerEnergy": params.get("cold_layer_energy", 0),
        # Microburst / downburst composites
        "wmsi": params.get("wmsi", 0),
        "mdpi": params.get("mdpi", 0),
        "maxGust": params.get("max_gust", 0),
        # Corfidi MCS motion vectors
        "corfidiUpSpd": params.get("corfidi_up_spd"),
        "corfidiDnSpd": params.get("corfidi_dn_spd"),
        "corfidiUpU": params.get("corfidi_up_u"),
        "corfidiUpV": params.get("corfidi_up_v"),
        "corfidiDnU": params.get("corfidi_dn_u"),
        "corfidiDnV": params.get("corfidi_dn_v"),
        # Fire weather indices
        "fosbergFwi": params.get("fosberg_fwi"),
        "haines": params.get("haines"),
        "hdw": params.get("hdw"),
        # Hazard classification
        "hazards": params.get("hazards", []),
        # Temperature advection
        "tempAdvection": params.get("temp_advection", []),
        # Kinematic
        "bwd500m": _fmt(params.get("bwd_500m")),
        "bwd1km": _fmt(params.get("bwd_1km")),
        "bwd3km": _fmt(params.get("bwd_3km")),
        "bwd6km": _fmt(params.get("bwd_6km")),
        "srh500m": _fmt(params.get("srh_500m")),
        "srh1km": _fmt(params.get("srh_1km")),
        "srh3km": _fmt(params.get("srh_3km")),
        "esrh": _fmt(params.get("esrh")),
        "ebwd": _fmt(params.get("ebwd")),
        "eilBot": round(params["eil_bot_h"]) if params.get("eil_bot_h") is not None else None,
        "eilTop": round(params["eil_top_h"]) if params.get("eil_top_h") is not None else None,
        # Convective mode
        "brn": params.get("brn"),
        "convectiveMode": params.get("convective_mode"),
        # Bunkers storm motion vectors (m/s)
        "rmU": _fmt(params.get("rm_u")),
        "rmV": _fmt(params.get("rm_v")),
        "lmU": _fmt(params.get("lm_u")),
        "lmV": _fmt(params.get("lm_v")),
        "mwU": _fmt(params.get("mw_u")),
        "mwV": _fmt(params.get("mw_v")),
        # Shear vector components (m/s) for hodograph rendering
        "bwdU500m": _fmt(params.get("bwd_u_500m"), decimals=1),
        "bwdV500m": _fmt(params.get("bwd_v_500m"), decimals=1),
        "bwdU1km": _fmt(params.get("bwd_u_1km"), decimals=1),
        "bwdV1km": _fmt(params.get("bwd_v_1km"), decimals=1),
        "bwdU3km": _fmt(params.get("bwd_u_3km"), decimals=1),
        "bwdV3km": _fmt(params.get("bwd_v_3km"), decimals=1),
        "bwdU6km": _fmt(params.get("bwd_u_6km"), decimals=1),
        "bwdV6km": _fmt(params.get("bwd_v_6km"), decimals=1),
        # Critical angle (degrees)
        "criticalAngle": params.get("critical_angle"),
        # Streamwiseness profile
        "streamwiseness": _nan_safe(params.get("streamwiseness", []).tolist()) if hasattr(params.get("streamwiseness", []), "tolist") else None,
        "streamwisenessSigned": _nan_safe(params.get("streamwiseness_signed", []).tolist()) if hasattr(params.get("streamwiseness_signed", []), "tolist") else None,
        "streamwisenessHeight": _nan_safe(params.get("streamwiseness_height", []).tolist()) if hasattr(params.get("streamwiseness_height", []), "tolist") else None,
        # Energy-Helicity Index & Vorticity Generation Parameter
        "ehi01": params.get("ehi_01", 0),
        "ehi03": params.get("ehi_03", 0),
        "vgp": params.get("vgp", 0),
    }
