const API_BASE = import.meta.env.VITE_API_URL || "";  // Empty in dev (Vite proxy), full URL in prod

/**
 * Wrapper around fetch that adds a timeout (default 10 s).
 * Throws on timeout so callers can retry or show an error.
 */
function fetchWithTimeout(url, options = {}, timeoutMs = 10000) {
  const controller = new AbortController();
  const timer = setTimeout(() => controller.abort(), timeoutMs);
  return fetch(url, { ...options, signal: controller.signal })
    .catch((err) => {
      if (err.name === "AbortError") throw new Error("Request timed out");
      throw err;
    })
    .finally(() => clearTimeout(timer));
}

/**
 * Retry a function up to `retries` times with `delayMs` between attempts.
 */
async function withRetry(fn, retries = 2, delayMs = 1500) {
  for (let i = 0; i <= retries; i++) {
    try {
      return await fn();
    } catch (err) {
      if (i === retries) throw err;
      await new Promise((r) => setTimeout(r, delayMs));
    }
  }
}

export async function fetchStations() {
  return withRetry(async () => {
    const res = await fetchWithTimeout(`${API_BASE}/api/stations`);
    if (!res.ok) throw new Error("Failed to fetch stations");
    return res.json();
  });
}

export async function fetchSources() {
  return withRetry(async () => {
    const res = await fetchWithTimeout(`${API_BASE}/api/sources`);
    if (!res.ok) throw new Error("Failed to fetch sources");
    return res.json();
  });
}

export async function fetchRiskScan(date) {
  const res = await fetchWithTimeout(`${API_BASE}/api/risk-scan`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(date ? { date } : {}),
  }, 180000);
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || "Risk scan failed");
  return data;
}

export async function fetchSounding(params) {
  const res = await fetchWithTimeout(`${API_BASE}/api/sounding`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(params),
  }, 120000);
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || "Request failed");
  return data;
}

export async function fetchTimeSeries(params) {
  const res = await fetchWithTimeout(`${API_BASE}/api/time-series`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(params),
  }, 300000);   // 5 min timeout — multiple fetches
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || "Time-series request failed");
  return data;
}

export async function fetchCompare(soundings) {
  const res = await fetchWithTimeout(`${API_BASE}/api/compare`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ soundings }),
  }, 300000);
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || "Comparison request failed");
  return data;
}

/**
 * Fetch a composite overlay plot (multiple profiles on one Skew-T).
 */
export async function fetchComposite(soundings) {
  const res = await fetchWithTimeout(`${API_BASE}/api/composite`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ soundings }),
  }, 300000);
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || "Composite request failed");
  return data;
}

/**
 * Fetch international (IGRAv2) station list, optionally filtered by region.
 */
export async function fetchIntlStations(region = "") {
  const url = region
    ? `${API_BASE}/api/stations/intl?region=${encodeURIComponent(region)}`
    : `${API_BASE}/api/stations/intl`;
  const res = await fetchWithTimeout(url);
  if (!res.ok) throw new Error("Failed to fetch international stations");
  return res.json();
}

/**
 * Fetch available international station regions with counts.
 */
export async function fetchIntlRegions() {
  const res = await fetchWithTimeout(`${API_BASE}/api/stations/intl/regions`);
  if (!res.ok) throw new Error("Failed to fetch regions");
  return res.json();
}

/**
 * Fetch SPC convective outlook GeoJSON for a given day (1, 2, or 3).
 */
export async function fetchSpcOutlook(day = 1) {
  const res = await fetchWithTimeout(`${API_BASE}/api/spc-outlook?day=${day}`, {}, 15000);
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || "Failed to fetch SPC outlook");
  return data;
}
