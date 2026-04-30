/**
 * Pure geometry helpers used by the map for nearest-radar lookup, view
 * bounds decimation, and storm-bearing labels.
 *
 * No DOM, no React, no Leaflet — safe to unit-test in jsdom or Node.
 */

const _toRad = (deg) => (deg * Math.PI) / 180;
const _toDeg = (rad) => (rad * 180) / Math.PI;

/**
 * Clamp a number into the inclusive `[min, max]` range.
 *
 * @param {number} v
 * @param {number} min
 * @param {number} max
 * @returns {number}
 */
export function clampNum(v, min, max) {
  return Math.max(min, Math.min(max, v));
}

/**
 * Initial-bearing (forward azimuth) from `(lat1,lon1)` to `(lat2,lon2)` in
 * degrees, normalized to `[0, 360)`. 0° is true north.
 *
 * Uses the standard great-circle formula; accurate for any pair of points
 * outside the polar singularities.
 */
export function bearingDeg(lat1, lon1, lat2, lon2) {
  const phi1 = _toRad(lat1);
  const phi2 = _toRad(lat2);
  const dLon = _toRad(lon2 - lon1);
  const y = Math.sin(dLon) * Math.cos(phi2);
  const x = Math.cos(phi1) * Math.sin(phi2) -
            Math.sin(phi1) * Math.cos(phi2) * Math.cos(dLon);
  return (_toDeg(Math.atan2(y, x)) + 360) % 360;
}

/**
 * Great-circle distance between two points in kilometers. Uses the
 * haversine formula with mean Earth radius 6371 km.
 */
export function haversineKm(lat1, lon1, lat2, lon2) {
  const R = 6371.0;
  const dLat = _toRad(lat2 - lat1);
  const dLon = _toRad(lon2 - lon1);
  const a = Math.sin(dLat / 2) ** 2 +
            Math.cos(_toRad(lat1)) * Math.cos(_toRad(lat2)) * Math.sin(dLon / 2) ** 2;
  return 2 * R * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
}

/** Maximum hour-offset the wind-field overlay scrubs across when animating. */
export const WIND_MAX_HOUR_OFFSET = 12;

/**
 * Map a frame index (0…frameCount-1) onto an integer hour-offset between
 * 0 and `maxOffset`. Used by the wind-field overlay so the scrub bar
 * advances uniformly across whatever number of frames is loaded.
 */
export function frameToHourOffset(frameIdx, frameCount, maxOffset = WIND_MAX_HOUR_OFFSET) {
  if (frameCount <= 1) return 0;
  const ratio = clampNum(frameIdx, 0, frameCount - 1) / (frameCount - 1);
  return Math.round(ratio * maxOffset);
}
