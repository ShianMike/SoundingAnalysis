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
 * Fetch SPC convective outlook GeoJSON for a given day (1, 2, or 3).
 */
export async function fetchSpcOutlook(day = 1) {
  const res = await fetchWithTimeout(`${API_BASE}/api/spc-outlook?day=${day}`, {}, 15000);
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || "Failed to fetch SPC outlook");
  return data;
}

/**
 * Fetch VWP time-height display image for a given NEXRAD radar.
 */
export async function fetchVwpDisplay(radar, hours = 12) {
  const res = await fetchWithTimeout(
    `${API_BASE}/api/vwp-display?radar=${encodeURIComponent(radar)}&hours=${hours}`,
    {},
    120000 // 2 min timeout — fetches many files
  );
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || "Failed to fetch VWP display");
  return data;
}

/**
 * Upload custom sounding text and receive analysis.
 */
export async function fetchCustomSounding({ text, format = "auto", theme = "dark", colorblind = false }) {
  const res = await fetchWithTimeout(`${API_BASE}/api/custom-sounding`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ text, format, theme, colorblind }),
  }, 120000);
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || "Custom sounding analysis failed");
  return data;
}

/**
 * Upload a binary file (WRF netCDF etc.) for analysis.
 */
export async function fetchUploadFile(formData) {
  const res = await fetchWithTimeout(`${API_BASE}/api/upload-file`, {
    method: "POST",
    body: formData,   // FormData — browser sets multipart boundary
  }, 180000);
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || "File upload analysis failed");
  return data;
}

/**
 * Merge two soundings into a weighted-average blended profile.
 */
export async function fetchMergeProfiles({ soundings, weight = 0.5, theme = "dark", colorblind = false }) {
  const res = await fetchWithTimeout(`${API_BASE}/api/merge-profiles`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ soundings, weight, theme, colorblind }),
  }, 300000);
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || "Profile merge failed");
  return data;
}

/**
 * Fetch ensemble sounding plume (multiple BUFKIT forecast hours overlaid).
 */
export async function fetchEnsemblePlume(params) {
  const res = await fetchWithTimeout(`${API_BASE}/api/ensemble-plume`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(params),
  }, 300000);
  const data = await res.json();
  if (!res.ok) {
    const err = new Error(data.error || "Ensemble plume request failed");
    err.suggestions = data.suggestions || [];
    throw err;
  }
  return data;
}
