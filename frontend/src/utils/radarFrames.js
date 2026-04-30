/**
 * Radar-frame URL builders for the IEM RIDGE archive used by the map's
 * radar reflectivity / velocity overlays.
 *
 * Pure helpers — no DOM, no React. Each returns an `oldest → newest`
 * array of frame descriptors `{ time (unix s), ts ("YYYYMMDDHHmm"), url? }`.
 */

/**
 * Format a UNIX-seconds timestamp as a short UTC `HH:MMZ` label
 * (e.g. `"14:05Z"`). Returns `"--:--Z"` for non-finite input.
 *
 * @param {number} ts UNIX seconds
 * @returns {string}
 */
export function fmtRadarTime(ts) {
  if (!Number.isFinite(ts)) return "--:--Z";
  const d = new Date(ts * 1000);
  return d.toISOString().slice(11, 16) + "Z";
}

/**
 * Build an IEM mosaic frame list — last `count × 5 minutes`, snapped to the
 * nearest 5-minute clock. IEM RIDGE tile URLs accept timestamps as
 * `YYYYMMDDHHmm` in the path, so the returned `ts` is in that exact form.
 *
 * @param {number} [count=24]
 * @returns {{time:number,ts:string}[]} oldest → newest
 */
export function buildMosaicFrames(count = 24) {
  const now = new Date();
  // Round down to nearest 5 minutes
  now.setUTCSeconds(0, 0);
  now.setUTCMinutes(Math.floor(now.getUTCMinutes() / 5) * 5);
  const frames = [];
  for (let i = count - 1; i >= 0; i--) {
    const t = new Date(now.getTime() - i * 5 * 60_000);
    const ts = t.getUTCFullYear().toString()
      + String(t.getUTCMonth() + 1).padStart(2, "0")
      + String(t.getUTCDate()).padStart(2, "0")
      + String(t.getUTCHours()).padStart(2, "0")
      + String(t.getUTCMinutes()).padStart(2, "0");
    frames.push({ time: Math.floor(t.getTime() / 1000), ts });
  }
  return frames;
}

/**
 * Build synthetic single-site archive frames at estimated 5-minute scan
 * intervals. Used as a fallback when IEM's scan-list API returns no rows
 * (e.g. for very recent times that haven't been indexed yet).
 *
 * The returned frames carry a fully-qualified IEM archive URL so callers
 * can drop them straight into a `<RadarLayer/>` without a follow-up
 * lookup.
 */
export function buildSingleSiteFrames(radar, product, count = 24) {
  const now = new Date();
  now.setUTCSeconds(0, 0);
  now.setUTCMinutes(Math.floor(now.getUTCMinutes() / 5) * 5);
  const frames = [];
  for (let i = count - 1; i >= 0; i--) {
    const t = new Date(now.getTime() - i * 5 * 60_000);
    const yyyy = t.getUTCFullYear().toString();
    const mm   = String(t.getUTCMonth() + 1).padStart(2, "0");
    const dd   = String(t.getUTCDate()).padStart(2, "0");
    const hh   = String(t.getUTCHours()).padStart(2, "0");
    const mi   = String(t.getUTCMinutes()).padStart(2, "0");
    const ts   = `${yyyy}${mm}${dd}${hh}${mi}`;
    const url  = `https://mesonet.agron.iastate.edu/archive/data/${yyyy}/${mm}/${dd}/GIS/ridge/${radar}/${product}/${radar}_${product}_${ts}.png`;
    frames.push({ url, time: Math.floor(t.getTime() / 1000), ts });
  }
  return frames;
}
