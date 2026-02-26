const API_BASE = import.meta.env.VITE_API_URL || "";  // Empty in dev (Vite proxy), full URL in prod

export async function fetchStations() {
  const res = await fetch(`${API_BASE}/api/stations`);
  if (!res.ok) throw new Error("Failed to fetch stations");
  return res.json();
}

export async function fetchSources() {
  const res = await fetch(`${API_BASE}/api/sources`);
  if (!res.ok) throw new Error("Failed to fetch sources");
  return res.json();
}

export async function fetchRiskScan(date) {
  const res = await fetch(`${API_BASE}/api/risk-scan`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(date ? { date } : {}),
  });
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || "Risk scan failed");
  return data;
}

export async function fetchSounding(params) {
  const res = await fetch(`${API_BASE}/api/sounding`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(params),
  });
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || "Request failed");
  return data;
}
