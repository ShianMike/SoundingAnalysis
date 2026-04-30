/**
 * Pure helpers + constants for NWS-feed warning classification & styling.
 *
 * No DOM, no React, no Leaflet — safe for unit tests in Node/jsdom.
 *
 * Classification follows a Weatherwise-style sub-category model so the map
 * can surface the underlying NWS hazard tier:
 *   Tornado:                    Radar-Indicated → Observed → PDS → Emergency
 *   Severe Thunderstorm:        Base → Considerable → Destructive
 *   Flash Flood / Watches / SPS:single tier each.
 */

/**
 * Per-sub-category style values consumed by both Leaflet (for the polygon
 * stroke / fill) and the legend chip rendered in the toolbar.
 */
export const WARNING_SUBCAT_STYLES = {
  // Tornado sub-categories (highest → lowest)
  TOR_EMG: { color: "#ff00ff", weight: 5,   fill: 0.30, label: "TORNADO EMERGENCY",        cls: "smap-warn-emg" },
  TOR_PDS: { color: "#ff1493", weight: 4,   fill: 0.22, label: "PDS TORNADO",              cls: "smap-warn-pds" },
  TOR_OBS: { color: "#dc143c", weight: 3.5, fill: 0.20, label: "OBSERVED TORNADO",         cls: "smap-warn-obs" },
  TOR_RI:  { color: "#ff0000", weight: 3,   fill: 0.18, label: "RADAR-INDICATED TORNADO",  cls: "smap-warn-ri"  },
  // Severe-thunderstorm sub-categories
  SVR_DST: { color: "#ff4500", weight: 3.5, fill: 0.22, label: "DESTRUCTIVE SVR",          cls: "smap-warn-svr-dst" },
  SVR_CON: { color: "#ff8c00", weight: 3,   fill: 0.20, label: "CONSIDERABLE SVR",         cls: "smap-warn-svr-con" },
  SVR:     { color: "#ffa500", weight: 2.5, fill: 0.16, label: "SVR",                      cls: "smap-warn-svr" },
  // Other event types
  FFW: { color: "#00ff00", weight: 2,   fill: 0.18, label: "FFW",                          cls: "smap-warn-ffw" },
  TOA: { color: "#ffff00", weight: 2,   fill: 0.08, label: "TOA", dash: "8 4",             cls: "smap-warn-toa" },
  SVA: { color: "#ffa500", weight: 2,   fill: 0.08, label: "SVA", dash: "8 4",             cls: "smap-warn-sva" },
  SPS: { color: "#ffe4b5", weight: 1.5, fill: 0.10, label: "SPS",                          cls: "smap-warn-sps" },
  FLW: { color: "#228b22", weight: 1.5, fill: 0.10, label: "FLW",                          cls: "smap-warn-flw" },
};

/** Set of NWS `event` strings the UI accepts; everything else is dropped. */
export const VALID_WARNING_EVENTS = new Set([
  "Tornado Warning",
  "Particularly Dangerous Situation Tornado Warning",
  "Severe Thunderstorm Warning",
  "Flash Flood Warning",
  "Tornado Watch",
  "Severe Thunderstorm Watch",
  "Special Weather Statement",
  "Flood Warning",
]);

/** Display priority for stacking / counting (highest first). */
export const WARNING_PRIORITY = [
  "TOR_EMG", "TOR_PDS", "TOR_OBS", "TOR_RI",
  "SVR_DST", "SVR_CON", "SVR",
  "FFW", "TOA", "SVA", "SPS", "FLW",
];

/**
 * Safely pull the first element of an NWS `parameters` array field.
 * Returns `""` when the field is missing or empty.
 *
 * @param {Record<string, unknown>} params
 * @param {string} key
 * @returns {string}
 */
export function firstParam(params, key) {
  const arr = params?.[key];
  if (!Array.isArray(arr) || arr.length === 0) return "";
  return String(arr[0] || "").trim();
}

/**
 * Classify an NWS GeoJSON feature into one of the keys of
 * `WARNING_SUBCAT_STYLES`. Returns `null` when the event is not in
 * `VALID_WARNING_EVENTS` (so callers can drop the feature entirely).
 *
 * Tornado warnings are decomposed via the `tornadoDamageThreat` /
 * `tornadoDetection` parameters and the headline text. Severe storms
 * follow the same pattern via `thunderstormDamageThreat`.
 *
 * @param {{ properties?: Record<string, unknown> }} feature
 * @returns {string | null}
 */
export function classifyWarning(feature) {
  const p = feature?.properties || {};
  const ev = p.event || "";
  const params = p.parameters || {};
  const headline = (p.headline || "").toUpperCase();
  const description = (p.description || "").toUpperCase();
  const text = `${headline} ${description}`;

  if (ev === "Tornado Warning" || ev === "Particularly Dangerous Situation Tornado Warning") {
    const damage = firstParam(params, "tornadoDamageThreat").toUpperCase();
    const detect = firstParam(params, "tornadoDetection").toUpperCase();
    if (damage === "CATASTROPHIC" || text.includes("TORNADO EMERGENCY")) return "TOR_EMG";
    if (damage === "CONSIDERABLE" || ev.startsWith("Particularly Dangerous") ||
        text.includes("PARTICULARLY DANGEROUS SITUATION")) return "TOR_PDS";
    if (detect === "OBSERVED" || text.includes("CONFIRMED TORNADO") ||
        text.includes("OBSERVED TORNADO")) return "TOR_OBS";
    return "TOR_RI";
  }
  if (ev === "Severe Thunderstorm Warning") {
    const damage = firstParam(params, "thunderstormDamageThreat").toUpperCase();
    if (damage === "DESTRUCTIVE") return "SVR_DST";
    if (damage === "CONSIDERABLE") return "SVR_CON";
    return "SVR";
  }
  if (ev === "Flash Flood Warning")        return "FFW";
  if (ev === "Tornado Watch")              return "TOA";
  if (ev === "Severe Thunderstorm Watch")  return "SVA";
  if (ev === "Special Weather Statement")  return "SPS";
  if (ev === "Flood Warning")              return "FLW";
  return null;
}

/**
 * Parse the NWS `parameters.eventMotionDescription` string format:
 *
 *   `"2025-04-30T01:50:00-05:00...260DEG...45KT...3848 9421"`
 *   → `{ dirDeg: 260, speedKt: 45, mph: 52 }`
 *
 * Any field may be `NaN` if the corresponding token is missing. Returns
 * `null` when no recognizable token is present.
 *
 * @param {string} motionStr
 * @returns {{ dirDeg: number, speedKt: number, mph: number } | null}
 */
export function parseEventMotion(motionStr) {
  if (!motionStr || typeof motionStr !== "string") return null;
  const dirMatch   = motionStr.match(/(\d{1,3})\s*DEG/i);
  const speedMatch = motionStr.match(/(\d{1,3})\s*KT/i);
  if (!dirMatch && !speedMatch) return null;
  const dirDeg  = dirMatch ? parseInt(dirMatch[1], 10) : NaN;
  const speedKt = speedMatch ? parseInt(speedMatch[1], 10) : NaN;
  const mph     = Number.isFinite(speedKt) ? Math.round(speedKt * 1.15078) : NaN;
  return { dirDeg, speedKt, mph };
}

/**
 * HTML-escape a string for safe injection into popup template literals.
 * Treats `null` / `undefined` as empty.
 *
 * @param {unknown} s
 * @returns {string}
 */
export function escapeHtml(s) {
  return String(s ?? "")
    .replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;").replace(/'/g, "&#39;");
}

/**
 * Build the Leaflet style object for an NWS warning polygon, derived from
 * `classifyWarning(feature)`. Falls back to `SPS` styling for unclassifiable
 * features so they still render (faintly) instead of disappearing.
 */
export function warningStyle(feature) {
  const key = classifyWarning(feature) || "SPS";
  const s = WARNING_SUBCAT_STYLES[key] || { color: "#ccc", weight: 1, fill: 0.1 };
  return {
    color: s.color,
    weight: s.weight,
    fillColor: s.color,
    fillOpacity: s.fill,
    dashArray: s.dash || null,
    className: `smap-warn-poly ${s.cls || ""}`.trim(),
  };
}
