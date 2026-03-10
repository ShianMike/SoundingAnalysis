import { useState, useEffect, useRef, useCallback, useMemo } from "react";
import { MapContainer, TileLayer, WMSTileLayer, ImageOverlay, CircleMarker, Marker, Circle, Popup, Tooltip, Pane, GeoJSON, Polyline, useMapEvents, useMap } from "react-leaflet";
import L from "leaflet";
import { X, Crosshair, CloudLightning, Wind, Maximize2, Minimize2, Layers, Play, Pause, Zap, RefreshCw, AlertTriangle, Tornado, Binoculars, FileWarning, Radio, ChevronDown, Star, Plane, Navigation, Search, FileText, Minus, Compass } from "lucide-react";
import { fetchSpcOutlook, fetchSpcOutlookStations, fetchSpcDiscussion, fetchWindField, fetchAcarsAirports } from "../api";
import { getFavorites } from "../favorites";
import WindCanvas from "./WindCanvas";
import "leaflet/dist/leaflet.css";
import "./StationMap.css";

/* ── Wrapper that keeps a SINGLE TileLayer alive and swaps URL via
   setUrl() instead of destroying / recreating the Leaflet layer.
   This avoids re-fetching every visible tile on each frame change. ── */
function RadarLayer({ url, opacity, ...rest }) {
  const ref = useRef(null);
  const prevUrl = useRef(url);
  useEffect(() => {
    if (ref.current) ref.current.setOpacity(opacity);
  }, [opacity]);
  useEffect(() => {
    if (ref.current && url !== prevUrl.current) {
      ref.current.setUrl(url);
      prevUrl.current = url;
    }
  }, [url]);
  return <TileLayer ref={ref} url={url} opacity={opacity} {...rest} />;
}

/* ── base-map tile options ───────────────────────────────────── */
const BASE_MAPS = [
  { id: "dark",      name: "Dark",      url: "https://{s}.basemaps.cartocdn.com/dark_nolabels/{z}/{x}/{y}{r}.png",  labels: "https://{s}.basemaps.cartocdn.com/dark_only_labels/{z}/{x}/{y}{r}.png" },
  { id: "light",     name: "Light",     url: "https://{s}.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}{r}.png", labels: "https://{s}.basemaps.cartocdn.com/light_only_labels/{z}/{x}/{y}{r}.png" },
  { id: "satellite", name: "Satellite", url: "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}", labels: "https://{s}.basemaps.cartocdn.com/dark_only_labels/{z}/{x}/{y}{r}.png" },
  { id: "terrain",   name: "Terrain",   url: "https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png", labels: "https://{s}.basemaps.cartocdn.com/dark_only_labels/{z}/{x}/{y}{r}.png" },
];

/* ── animated radar via RainViewer API (live composite, real timestamps) ──── */
const RADAR_INTERVAL_MS = 1200;         // ms between animation frames (slower to avoid 429s)
const RAINVIEWER_API    = "https://api.rainviewer.com/public/weather-maps.json";
const RAINVIEWER_TILE   = "https://tilecache.rainviewer.com";
const RAINVIEWER_OPTS   = "/256/{z}/{x}/{y}/1/1_1.png";   // size/color(1=transparent-bg)/smooth/snow

/** Fetch available radar frame paths from RainViewer (with retry). */
async function fetchRadarFrames(retries = 2) {
  for (let attempt = 0; attempt <= retries; attempt++) {
    try {
      const res = await fetch(RAINVIEWER_API);
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      const data = await res.json();
      const past    = (data.radar?.past    || []).filter(f => f.path).slice(-12);  // last 12 frames
      const nowcast = (data.radar?.nowcast || []).filter(f => f.path).slice(0, 3); // up to 3 forecast frames
      return { past, nowcast, generated: data.generated };
    } catch (e) {
      if (attempt < retries) {
        await new Promise((r) => setTimeout(r, 2000 * (attempt + 1))); // back off
        continue;
      }
      console.error("RainViewer API error:", e);
      return { past: [], nowcast: [], generated: 0 };
    }
  }
}

/** Format a UNIX timestamp to a short UTC label like "14:05Z". */
function fmtRadarTime(ts) {
  const d = new Date(ts * 1000);
  return d.toISOString().slice(11, 16) + "Z";
}

/** Build IEM mosaic frame list — last ~2 hours at 5-min intervals.
 *  IEM RIDGE tiles accept timestamps as YYYYMMDDHHmm in the URL path.
 *  Returns [{time (unix), ts (YYYYMMDDHHmm string)}] oldest→newest. */
function buildMosaicFrames(count = 24) {
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

const WIND_MAX_HOUR_OFFSET = 12;

function clampNum(v, min, max) {
  return Math.max(min, Math.min(max, v));
}

function frameToHourOffset(frameIdx, frameCount, maxOffset = WIND_MAX_HOUR_OFFSET) {
  if (frameCount <= 1) return 0;
  const ratio = clampNum(frameIdx, 0, frameCount - 1) / (frameCount - 1);
  return Math.round(ratio * maxOffset);
}

/* ── NWS Active Warnings ─────────────────────────────────────── */
const NWS_ALERTS_API = "https://api.weather.gov/alerts/active?status=actual";

/* ── SPC Mesoscale Discussions & Watches ──────────────────────── */
const SPC_WATCH_URL =
  "https://mapservices.weather.noaa.gov/eventdriven/rest/services/WWA/watch_warn_adv/MapServer/1/query" +
  "?where=(phenom%3D%27TO%27+OR+phenom%3D%27SV%27)+AND+sig%3D%27A%27" +
  "&f=geojson&outFields=prod_type,event,sig,phenom,wfo,expiration,issuance" +
  "&returnGeometry=true&outSR=4326";

const WATCH_STYLES = {
  "Tornado Watch":           { color: "#ff0000", fill: "#ff000022", label: "TOR" },
  "Severe Thunderstorm Watch": { color: "#ffa500", fill: "#ffa50022", label: "SVR" },
};

async function fetchSpcMds() {
  // IEM removed the spc_mcd.geojson endpoint; no alternative GeoJSON source exists
  return { type: "FeatureCollection", features: [] };
}

async function fetchSpcWatches() {
  try {
    const res = await fetch(SPC_WATCH_URL);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const data = await res.json();
    // Map NWS MapServer property names to what MdWatchLayer expects
    for (const f of data.features || []) {
      f.properties.type = f.properties.prod_type;
      f.properties.number = f.properties.event;
    }
    return data;
  } catch (e) {
    console.error("SPC Watch fetch error:", e);
    return { type: "FeatureCollection", features: [] };
  }
}

function mdStyle() {
  return {
    color: "#00e5ff",
    weight: 2,
    fillColor: "#00e5ff",
    fillOpacity: 0.08,
    dashArray: "6 4",
  };
}

function watchStyle(feature) {
  const ev = feature.properties?.type || feature.properties?.event || "";
  const s = Object.entries(WATCH_STYLES).find(([k]) => ev.includes(k));
  if (s) return { color: s[1].color, weight: 2.5, fillColor: s[1].color, fillOpacity: 0.1, dashArray: "10 5" };
  return { color: "#ffa500", weight: 2, fillColor: "#ffa500", fillOpacity: 0.08, dashArray: "10 5" };
}

function MdWatchLayer({ mdData, watchData }) {
  const hasMd = mdData?.features?.length > 0;
  const hasWatch = watchData?.features?.length > 0;
  if (!hasMd && !hasWatch) return null;
  return (
    <Pane name="spc-md-watch" style={{ zIndex: 415, pointerEvents: "auto" }}>
      {hasWatch && (
        <GeoJSON
          key={"w-" + watchData.features.length + (watchData.features[0]?.properties?.number || "")}
          data={watchData}
          style={watchStyle}
          onEachFeature={(feature, layer) => {
            const p = feature.properties;
            const ev = p.type || p.event || "Watch";
            const isTor = ev.includes("Tornado");
            const num = p.number || p.wn || "";
            const expires = p.expires || p.expiration || "";
            const expStr = expires ? new Date(expires).toLocaleString() : "";
            layer.bindPopup(
              `<div class="smap-warn-popup">
                <div class="smap-warn-badge" style="background:${isTor ? '#ff000022' : '#ffa50022'};color:${isTor ? '#ff0000' : '#ffa500'};border-color:${isTor ? '#ff000055' : '#ffa50055'}">
                  ${isTor ? 'TOR' : 'SVR'} WATCH ${num}
                </div>
                <div class="smap-warn-headline">${ev}${num ? ` #${num}` : ''}</div>
                ${expStr ? `<div class="smap-warn-time">Expires: ${expStr}</div>` : ""}
              </div>`,
              { className: "smap-popup smap-warn-popup-wrap", minWidth: 220, maxWidth: 320 }
            );
          }}
        />
      )}
      {hasMd && (
        <GeoJSON
          key={"md-" + mdData.features.length + (mdData.features[0]?.properties?.number || "")}
          data={mdData}
          style={mdStyle}
          onEachFeature={(feature, layer) => {
            const p = feature.properties;
            const num = p.number || p.md_num || "";
            const concerning = p.concerning || p.hazard || "";
            const areas = p.areas || p.states || "";
            layer.bindPopup(
              `<div class="smap-warn-popup">
                <div class="smap-warn-badge" style="background:#00e5ff22;color:#00e5ff;border-color:#00e5ff55">
                  MD ${num}
                </div>
                <div class="smap-warn-headline">Mesoscale Discussion ${num}</div>
                ${concerning ? `<div class="smap-warn-area">Concerning: ${concerning}</div>` : ""}
                ${areas ? `<div class="smap-warn-area">${areas}</div>` : ""}
              </div>`,
              { className: "smap-popup smap-warn-popup-wrap", minWidth: 220, maxWidth: 320 }
            );
          }}
        />
      )}
    </Pane>
  );
}

/* ── Lightning overlay (Blitzortung) ─────────────────────────── */
const BLITZ_WS_SERVERS = [
  "wss://ws1.blitzortung.org/",
];
const LIGHTNING_BOUNDS = { west: -130, east: -60, north: 55, south: 20 };

function decodeBlitzPayload(ciphertext) {
  // Mirrors the LZW-like decoder used by known Blitzortung clients.
  const chars = [...String(ciphertext || "")];
  if (!chars.length) return "";

  let c = chars[0];
  let prev = c;
  let out = c;
  const dict = [];

  for (let i = 1; i < chars.length; i += 1) {
    const ch = chars[i];
    const code = ch.codePointAt(0);
    const a = code < 256
      ? ch
      : (dict[code - 256] != null ? dict[code - 256] : `${prev}${c}`);

    out += a;
    c = [...a][0] || "";
    dict.push(`${prev}${c}`);
    prev = a;
  }

  return out;
}

function normalizeStrikeTime(rawTime) {
  if (rawTime == null) return Date.now();
  const t = Number(rawTime);
  if (!Number.isFinite(t)) return Date.now();
  if (t > 1e15) return t / 1e6; // unix nanos -> millis
  if (t > 1e12) return t / 1e3; // unix micros -> millis
  if (t > 1e10) return t;       // unix millis
  return t * 1000;              // unix seconds
}

function parseLightningMessage(raw) {
  if (typeof raw !== "string") return null;

  const parseJson = (txt) => {
    try { return JSON.parse(txt); } catch { return null; }
  };
  const parsed = parseJson(raw) || parseJson(decodeBlitzPayload(raw));
  if (!parsed) return null;

  const lat = Number(parsed.lat);
  const lon = Number(parsed.lon);
  if (!Number.isFinite(lat) || !Number.isFinite(lon)) return null;

  return {
    lat,
    lon,
    time: normalizeStrikeTime(parsed.time),
  };
}

function useLightningData(enabled) {
  const [strikes, setStrikes] = useState([]);
  const wsRef = useRef(null);
  const strikesRef = useRef([]);

  useEffect(() => {
    if (!enabled) {
      setStrikes([]);
      strikesRef.current = [];
      if (wsRef.current) { wsRef.current.close(); wsRef.current = null; }
      return;
    }

    let reconnectTimer = null;
    const MAX_STRIKES = 2000;
    const STRIKE_TTL = 300_000; // 5 min

    function connect() {
      if (wsRef.current) return;
      try {
        const url = BLITZ_WS_SERVERS[Math.floor(Math.random() * BLITZ_WS_SERVERS.length)];
        const ws = new WebSocket(url);
        wsRef.current = ws;

        ws.onopen = () => {
          // Required handshake for the live strike stream.
          ws.send(JSON.stringify({ a: 111 }));
          ws.send(JSON.stringify(LIGHTNING_BOUNDS));
        };

        ws.onmessage = (event) => {
          const strike = parseLightningMessage(event.data);
          if (!strike) return;
          if (
            strike.lat < LIGHTNING_BOUNDS.south ||
            strike.lat > LIGHTNING_BOUNDS.north ||
            strike.lon < LIGHTNING_BOUNDS.west ||
            strike.lon > LIGHTNING_BOUNDS.east
          ) return;

          const now = Date.now();
          strikesRef.current = [
            { ...strike, id: `${Math.round(strike.time)}_${strike.lat.toFixed(3)}_${strike.lon.toFixed(3)}` },
            ...strikesRef.current.filter((s) => now - s.time < STRIKE_TTL),
          ].slice(0, MAX_STRIKES);
        };

        ws.onclose = () => {
          wsRef.current = null;
          reconnectTimer = setTimeout(connect, 5000);
        };
        ws.onerror = () => ws.close();
      } catch {
        reconnectTimer = setTimeout(connect, 5000);
      }
    }

    connect();

    // Batch update the React state every 2s to avoid per-strike renders
    const renderInterval = setInterval(() => {
      setStrikes([...strikesRef.current]);
    }, 2000);

    // Prune old strikes every 30s
    const pruneInterval = setInterval(() => {
      const now = Date.now();
      strikesRef.current = strikesRef.current.filter((s) => now - s.time < STRIKE_TTL);
    }, 30_000);

    return () => {
      if (wsRef.current) { wsRef.current.close(); wsRef.current = null; }
      clearTimeout(reconnectTimer);
      clearInterval(renderInterval);
      clearInterval(pruneInterval);
    };
  }, [enabled]);

  return strikes;
}

// Event → colour + short label
const WARNING_STYLES = {
  "Tornado Warning":              { color: "#ff0000", weight: 3, label: "TOR" },
  "Particularly Dangerous Situation Tornado Warning": { color: "#ff0000", weight: 4, label: "PDS TOR" },
  "Severe Thunderstorm Warning":   { color: "#ffa500", weight: 2.5, label: "SVR" },
  "Flash Flood Warning":           { color: "#00ff00", weight: 2, label: "FFW" },
  "Tornado Watch":                 { color: "#ffff00", weight: 2, label: "TOA", dash: "8 4" },
  "Severe Thunderstorm Watch":     { color: "#ffa500", weight: 2, label: "SVA", dash: "8 4" },
  "Special Weather Statement":     { color: "#ffe4b5", weight: 1.5, label: "SPS" },
  "Flood Warning":                 { color: "#228b22", weight: 1.5, label: "FLW" },
};
// Ordered priority for display (highest first)
const WARNING_PRIORITY = [
  "Tornado Warning", "Particularly Dangerous Situation Tornado Warning",
  "Severe Thunderstorm Warning", "Flash Flood Warning",
  "Tornado Watch", "Severe Thunderstorm Watch",
  "Special Weather Statement", "Flood Warning",
];

async function fetchNwsWarnings() {
  try {
    const res = await fetch(NWS_ALERTS_API, {
      headers: { "Accept": "application/geo+json", "User-Agent": "SoundingAnalysis/1.0" },
    });
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const data = await res.json();
    // Keep only events we have styles for, and that have geometry
    const validEvents = new Set(Object.keys(WARNING_STYLES));
    const features = (data.features || []).filter(
      (f) => f.geometry && validEvents.has(f.properties?.event)
    );
    return { type: "FeatureCollection", features };
  } catch (e) {
    console.error("NWS warnings fetch error:", e);
    return { type: "FeatureCollection", features: [] };
  }
}

function warningStyle(feature) {
  const ev = feature.properties?.event || "";
  const s = WARNING_STYLES[ev] || { color: "#ccc", weight: 1 };
  return {
    color: s.color,
    weight: s.weight,
    fillColor: s.color,
    fillOpacity: ev.includes("Watch") ? 0.08 : 0.18,
    dashArray: s.dash || null,
  };
}

function WarningsLayer({ data }) {
  if (!data || !data.features || data.features.length === 0) return null;
  return (
    <Pane name="nws-warnings" style={{ zIndex: 420, pointerEvents: "auto" }}>
      <GeoJSON
        key={data.features.length + (data.features[0]?.properties?.id || "")}
        data={data}
        style={warningStyle}
        onEachFeature={(feature, layer) => {
          const p = feature.properties;
          const ev = p.event || "Alert";
          const s = WARNING_STYLES[ev];
          const headline = p.headline || ev;
          const areas = p.areaDesc || "";
          const onset = p.onset ? new Date(p.onset).toLocaleString() : "";
          const expires = p.expires ? new Date(p.expires).toLocaleString() : "";
          layer.bindPopup(
            `<div class="smap-warn-popup">
              <div class="smap-warn-badge" style="background:${s?.color || '#ccc'}22;color:${s?.color || '#ccc'};border-color:${s?.color || '#ccc'}55">
                ${s?.label || "WARN"}
              </div>
              <div class="smap-warn-headline">${headline}</div>
              <div class="smap-warn-area">${areas}</div>
              ${onset ? `<div class="smap-warn-time">Onset: ${onset}</div>` : ""}
              ${expires ? `<div class="smap-warn-time">Expires: ${expires}</div>` : ""}
            </div>`,
            { className: "smap-popup smap-warn-popup-wrap", minWidth: 240, maxWidth: 340 }
          );
        }}
      />
    </Pane>
  );
}

/* ── TVS / Mesocyclone custom map icons ─────────────────────── */
const TVS_ICON = L.divIcon({
  className: "smap-tvs-icon",
  iconSize: [28, 28],
  iconAnchor: [14, 24],
  html: `<svg width="28" height="28" viewBox="0 0 28 28" xmlns="http://www.w3.org/2000/svg">
    <polygon points="14,26 1,3 27,3" fill="#dc2626" stroke="#fff" stroke-width="2" stroke-linejoin="round"/>
    <line x1="14" y1="8" x2="14" y2="17" stroke="#fff" stroke-width="2.5" stroke-linecap="round"/>
    <circle cx="14" cy="21" r="1.5" fill="#fff"/>
  </svg>`,
});

function makeMesoIcon(color, rank) {
  const sz = Math.max(22, 18 + rank * 2);
  const half = sz / 2;
  const r = half - 2;
  const ir = r * 0.55;
  const dash = Math.max(2.5, r * 0.28);
  return L.divIcon({
    className: "smap-meso-icon",
    iconSize: [sz, sz],
    iconAnchor: [half, half],
    html: `<svg width="${sz}" height="${sz}" viewBox="0 0 ${sz} ${sz}" xmlns="http://www.w3.org/2000/svg">
      <circle cx="${half}" cy="${half}" r="${r}" fill="${color}" fill-opacity="0.55" stroke="#fff" stroke-width="2"/>
      <circle cx="${half}" cy="${half}" r="${ir}" fill="none" stroke="#fff" stroke-width="1.5" stroke-dasharray="${dash},${dash * 0.9}" opacity="0.85"/>
      <circle cx="${half}" cy="${half}" r="1.5" fill="#fff" opacity="0.9"/>
    </svg>`,
  });
}

/* ── NEXRAD storm attribute (TVS / Mesocyclone) feed ───────── */
const IEM_STORM_ATTR = "https://mesonet.agron.iastate.edu/geojson/nexrad_attr.geojson";

async function fetchStormAttr() {
  try {
    const res = await fetch(IEM_STORM_ATTR);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const data = await res.json();
    // Only keep features that have TVS or Meso detection
    const features = (data.features || []).filter(
      (f) => f.properties?.tvs !== "NONE" || f.properties?.meso !== "NONE"
    );
    return { type: "FeatureCollection", features };
  } catch (e) {
    console.error("Storm attr fetch error:", e);
    return { type: "FeatureCollection", features: [] };
  }
}

/* ── Spotter Network live chaser positions + reports ─────── */
const SN_POSITIONS = "https://www.spotternetwork.org/feeds/gr.txt";
const SN_REPORTS   = "https://www.spotternetwork.org/feeds/reports.txt";

// Report icon-index → type mapping (from SN sprite sheet)
const SN_REPORT_TYPES = { 1: "Tornado", 2: "Funnel Cloud", 3: "Rotating Wall Cloud", 4: "Hail", 5: "Wind Damage", 6: "Flooding", 7: "Other" };
const SN_REPORT_COLORS = {
  Tornado:                "#ff0000",
  "Funnel Cloud":         "#ff6600",
  "Rotating Wall Cloud":  "#ff9900",
  Hail:                   "#00ccff",
  "Wind Damage":          "#ffaa00",
  Flooding:               "#00ff66",
  Other:                  "#cccccc",
};

/** Parse Spotter Network GR-format position feed → [{lat,lon,name,time,heading}] */
function parseSpotterPositions(text) {
  const lines = text.split(/\r?\n/);
  const spotters = [];
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i].trim();
    if (!line.startsWith("Object:")) continue;
    const coords = line.slice(7).trim().split(",");
    const lat = parseFloat(coords[0]), lon = parseFloat(coords[1]);
    if (isNaN(lat) || isNaN(lon)) continue;
    let name = "", time = "", heading = "";
    // Look ahead for Icon line with tooltip
    for (let j = i + 1; j < Math.min(i + 5, lines.length); j++) {
      const l = lines[j].trim();
      if (l.startsWith("Icon:") && l.includes('"')) {
        const m = l.match(/"([^"]*)"/); // tooltip in quotes
        if (m) {
          const parts = m[1].split("\\n");
          name = parts[0] || "";
          time = parts[1] || "";
          heading = parts[2] || "";
        }
        break;
      }
    }
    spotters.push({ lat, lon, name, time, heading });
  }
  return spotters;
}

/** Parse Spotter Network report feed → [{lat,lon,reporter,type,time,notes,sheet}] */
function parseSpotterReports(text) {
  const lines = text.split(/\r?\n/);
  const reports = [];
  for (const line of lines) {
    const l = line.trim();
    if (!l.startsWith("Icon:")) continue;
    const m = l.match(/^Icon:\s*([\d.-]+),([\d.-]+),\d+,(\d+),(\d+),"(.+)"$/);
    if (!m) continue;
    const lat = parseFloat(m[1]), lon = parseFloat(m[2]);
    const sheet = parseInt(m[3], 10);    // 3=recent, 4=older, 5=oldest
    const idx = parseInt(m[4], 10);      // report type index
    const tooltip = m[5];
    const parts = tooltip.split("\\n");
    const reporter = (parts[0] || "").replace(/^Reported By:\s*/, "");
    const type = parts[1] || "Other";
    const timePart = (parts[2] || "").replace(/^Time:\s*/, "");
    // Remaining parts may have Size / Notes
    const rest = parts.slice(3).join(" · ").replace(/^(Size:|Notes:)\s*/g, "");
    reports.push({ lat, lon, reporter, type, time: timePart, notes: rest, sheet, idx });
  }
  return reports;
}

async function fetchSpotterNetwork() {
  try {
    const [posRes, repRes] = await Promise.all([
      fetch(SN_POSITIONS),
      fetch(SN_REPORTS),
    ]);
    const positions = posRes.ok ? parseSpotterPositions(await posRes.text()) : [];
    const reports = repRes.ok ? parseSpotterReports(await repRes.text()) : [];
    return { positions, reports };
  } catch (e) {
    console.error("Spotter Network fetch error:", e);
    return { positions: [], reports: [] };
  }
}

/* ── colours ────────────────────────────────────────────────── */
const RISK_COLORS = {
  high: "#ef4444",   // red
  med:  "#facc15",   // yellow
  low:  "#22c55e",   // green
  none: "#555",      // neutral
};

function riskColor(stp) {
  if (stp == null) return RISK_COLORS.none;
  if (stp >= 1)    return RISK_COLORS.high;
  if (stp >= 0.3)  return RISK_COLORS.med;
  return RISK_COLORS.low;
}

/* ── Velocity product definitions ─────────────────────────────── */
/* IEM RIDGE only archives N0S (Storm-Relative Mean Velocity) for NEXRAD sites */
const VELOCITY_PRODUCTS = [
  { id: "N0S", label: "SRM 0.5°", desc: "Storm-Relative Mean Velocity (0.5° tilt) — removes storm motion, best for mesocyclone rotation couplets and tornado signatures" },
];

const VEL_ANIM_INTERVAL_MS = 800;
const VEL_FRAME_HOURS = 1; // how many hours of scans to fetch

/** Parse IEM scan timestamps → archive image URLs + bounds from world file.
 *  World file format: pixelSizeX, 0, 0, -pixelSizeY, topLeftX, topLeftY */
function parseWorldFile(text) {
  const lines = text.trim().split(/\r?\n/).map(Number);
  return { dx: lines[0], dy: lines[3], x0: lines[4], y0: lines[5] };
}

function velImageBounds(wld, width = 1000, height = 1000) {
  // wld: {dx, dy, x0, y0} — dy is negative
  const south = wld.y0 + wld.dy * height;
  const north = wld.y0;
  const west = wld.x0;
  const east = wld.x0 + wld.dx * width;
  return [[south, west], [north, east]];
}

/** Fetch velocity scan list + world file for a radar/product, return frames. */
async function fetchVelFrames(radar, product, hours = VEL_FRAME_HOURS) {
  const end = new Date();
  const start = new Date(end.getTime() - hours * 3600_000);
  const iso = (d) => d.toISOString().replace(/\.\d+Z/, "Z");

  // 1. Get scan list first
  const listRes = await fetch(`https://mesonet.agron.iastate.edu/json/radar.py?operation=list&radar=${radar}&product=${product}&start=${iso(start)}&end=${iso(end)}`);
  if (!listRes.ok) throw new Error(`HTTP ${listRes.status}`);
  const listData = await listRes.json();
  const scans = listData.scans || [];
  if (!scans.length) return { frames: [], bounds: null };

  // 2. Use the first actual scan's timestamp to fetch the world file (always exists)
  const firstTs = scans[0].ts.replace(/[-:TZ]/g, "").slice(0, 12);
  const firstDate = scans[0].ts.slice(0, 10).split("-");
  const wldUrl = `https://mesonet.agron.iastate.edu/archive/data/${firstDate[0]}/${firstDate[1]}/${firstDate[2]}/GIS/ridge/${radar}/${product}/${radar}_${product}_${firstTs}.wld`;
  const wldRes = await fetch(wldUrl).catch(() => null);
  if (!wldRes || !wldRes.ok) return { frames: [], bounds: null };
  const wld = parseWorldFile(await wldRes.text());
  const bounds = velImageBounds(wld);

  // 3. Build frame objects with archive image URLs
  const frames = scans.map((s) => {
    const ts = s.ts.replace(/[-:TZ]/g, "").slice(0, 12); // YYYYMMDDHHmm
    const d = s.ts.slice(0, 10).split("-"); // [YYYY, MM, DD]
    const url = `https://mesonet.agron.iastate.edu/archive/data/${d[0]}/${d[1]}/${d[2]}/GIS/ridge/${radar}/${product}/${radar}_${product}_${ts}.png`;
    const time = Math.floor(new Date(s.ts.endsWith("Z") ? s.ts : s.ts + "Z").getTime() / 1000);
    return { url, time, ts };
  });
  return { frames, bounds };
}

/* ── NEXRAD WSR-88D sites (CONUS) for velocity overlay ──────── */
const NEXRAD_SITES = [
  ["ABR",45.46,-98.41],["ABX",35.15,-106.82],["AKQ",36.98,-77.01],
  ["AMA",35.23,-101.71],["AMX",25.61,-80.41],["APX",44.91,-84.72],
  ["ARX",43.82,-91.19],["ATX",48.19,-122.50],["BBX",39.50,-121.63],
  ["BGM",42.20,-75.98],["BMX",33.17,-86.77],["BOX",41.96,-71.14],
  ["BRO",25.92,-97.42],["BUF",42.95,-78.74],["BYX",24.60,-81.70],
  ["CAE",33.95,-81.12],["CBW",46.04,-67.81],["CBX",43.49,-116.24],
  ["CCX",40.92,-78.00],["CLE",41.41,-81.86],["CLX",32.66,-81.04],
  ["CRP",27.78,-97.51],["CXX",44.51,-73.17],["CYS",41.15,-104.81],
  ["DAX",38.50,-121.68],["DDC",37.76,-99.97],["DFX",29.27,-100.28],
  ["DGX",32.28,-89.98],["DIX",39.95,-74.41],["DLH",46.84,-92.21],
  ["DMX",41.73,-93.72],["DOX",38.83,-75.44],["DTX",42.70,-83.47],
  ["DVN",41.61,-90.58],["DYX",32.54,-99.25],["EAX",38.81,-94.26],
  ["EMX",31.89,-110.63],["ENX",42.59,-74.06],["EOX",31.46,-85.46],
  ["EPZ",31.87,-106.70],["ESX",35.70,-114.89],["EVX",30.56,-85.92],
  ["EWX",29.70,-98.03],["EYX",35.10,-117.56],["FCX",37.02,-80.27],
  ["FDR",34.36,-98.98],["FDX",34.64,-103.63],["FFC",33.36,-84.57],
  ["FSD",43.59,-96.73],["FSX",34.57,-111.20],["FTG",39.79,-104.55],
  ["FWS",32.57,-97.30],["GGW",48.21,-106.63],["GJX",39.06,-108.21],
  ["GLD",39.37,-101.70],["GRB",44.50,-88.11],["GRK",30.72,-97.38],
  ["GRR",42.89,-85.54],["GSP",34.88,-82.22],["GWX",33.90,-88.33],
  ["GYX",43.89,-70.26],["HDX",33.08,-106.12],["HGX",29.47,-95.08],
  ["HNX",36.31,-119.63],["HPX",36.74,-87.28],["HTX",34.93,-86.08],
  ["HWA",38.51,-82.97],["ICT",37.65,-97.44],["ICX",37.59,-112.86],
  ["ILN",39.42,-83.82],["ILX",40.15,-89.34],["IND",39.71,-86.28],
  ["INX",36.18,-95.56],["IWA",33.29,-111.67],["IWX",41.36,-85.70],
  ["JAX",30.48,-81.70],["JGX",32.68,-83.35],["JKL",37.59,-83.31],
  ["KEY",24.55,-81.78],["KSG",31.48,-82.31],["LBB",33.65,-101.81],
  ["LCH",30.13,-93.22],["LIX",30.34,-89.83],["LNX",41.96,-100.58],
  ["LOT",41.60,-88.08],["LRX",40.74,-116.80],["LSX",38.70,-90.68],
  ["LTX",33.99,-78.43],["LVX",37.98,-85.94],["LWX",38.98,-77.48],
  ["LZK",34.84,-92.26],["MAF",31.94,-102.19],["MAX",42.08,-122.72],
  ["MBX",48.39,-100.86],["MHX",34.78,-76.88],["MKX",42.97,-88.55],
  ["MLB",28.11,-80.65],["MOB",30.68,-88.24],["MPX",44.85,-93.57],
  ["MQT",46.53,-87.55],["MRX",36.17,-83.40],["MSX",47.04,-113.99],
  ["MTX",41.26,-112.45],["MUX",37.16,-121.90],["MVX",47.53,-97.33],
  ["MXX",32.54,-85.79],["NKX",32.92,-117.04],["NQA",35.34,-89.87],
  ["OAX",41.32,-96.37],["OHX",36.25,-86.56],["OKX",40.87,-72.86],
  ["OTX",47.68,-117.63],["PAH",37.07,-88.77],["PBZ",40.53,-80.22],
  ["PDT",45.69,-118.85],["POE",34.41,-116.16],["PUX",38.46,-104.18],
  ["RAX",35.67,-78.49],["RGX",39.75,-119.46],["RIW",43.07,-108.48],
  ["RLX",38.31,-81.72],["RTX",45.71,-122.97],["SFX",43.11,-112.69],
  ["SGF",37.24,-93.40],["SHV",32.45,-93.84],["SJT",31.37,-100.49],
  ["SOX",33.82,-117.64],["SRX",35.29,-94.36],["TBW",27.71,-82.40],
  ["TFX",47.46,-111.39],["TLH",30.40,-84.33],["TLX",35.33,-97.28],
  ["TWX",38.99,-96.23],["TYX",43.76,-75.68],["UDX",44.13,-102.83],
  ["UEX",40.32,-98.44],["VAX",30.89,-83.00],["VBX",34.84,-120.40],
  ["VNX",36.74,-98.13],["VTX",34.41,-119.18],["VWX",38.26,-87.72],
  ["YUX",32.50,-114.66],
];

function nearestNexrad(lat, lon) {
  let best = NEXRAD_SITES[0], bestD = Infinity;
  for (const s of NEXRAD_SITES) {
    const d = (s[1] - lat) ** 2 + (s[2] - lon) ** 2;
    if (d < bestD) { bestD = d; best = s; }
  }
  return { id: best[0], lat: best[1], lon: best[2] };
}

function _toRad(deg) { return (deg * Math.PI) / 180; }
function _toDeg(rad) { return (rad * 180) / Math.PI; }

function bearingDeg(lat1, lon1, lat2, lon2) {
  const phi1 = _toRad(lat1);
  const phi2 = _toRad(lat2);
  const dLon = _toRad(lon2 - lon1);
  const y = Math.sin(dLon) * Math.cos(phi2);
  const x = Math.cos(phi1) * Math.sin(phi2) -
            Math.sin(phi1) * Math.cos(phi2) * Math.cos(dLon);
  return (_toDeg(Math.atan2(y, x)) + 360) % 360;
}

function haversineKm(lat1, lon1, lat2, lon2) {
  const R = 6371.0;
  const dLat = _toRad(lat2 - lat1);
  const dLon = _toRad(lon2 - lon1);
  const a = Math.sin(dLat / 2) ** 2 +
            Math.cos(_toRad(lat1)) * Math.cos(_toRad(lat2)) * Math.sin(dLon / 2) ** 2;
  return 2 * R * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
}

/* ── Track map center + view bounds for nearest-radar and decimation ───────────── */
function MapViewTracker({ onCenterChange, onViewChange }) {
  const map = useMap();
  const report = useCallback(() => {
    const c = map.getCenter();
    const b = map.getBounds();
    onCenterChange?.(c.lat, c.lng);
    onViewChange?.({
      zoom: map.getZoom(),
      south: b.getSouth(),
      west: b.getWest(),
      north: b.getNorth(),
      east: b.getEast(),
    });
  }, [map, onCenterChange, onViewChange]);

  useMapEvents({
    moveend() {
      report();
    },
    zoomend() {
      report();
    },
  });

  // Fire once on mount
  useEffect(() => {
    report();
  }, [report]);

  return null;
}

/* ── SPC outlook category styling ───────────────────────────── */
const SPC_CATEGORIES = [
  { label: "TSTM", name: "General Thunder", color: "#55BB55", fill: "#C1E9C1",
    info: "Non-severe thunderstorms possible" },
  { label: "MRGL", name: "Marginal",        color: "#005500", fill: "#66A366",
    info: "5% severe · 2% tornado" },
  { label: "SLGT", name: "Slight",          color: "#DDB600", fill: "#FFE066",
    info: "15% severe · 5% tornado" },
  { label: "ENH",  name: "Enhanced",        color: "#FF6600", fill: "#FFA366",
    info: "30% severe · 10% tornado" },
  { label: "MDT",  name: "Moderate",        color: "#CC0000", fill: "#E06666",
    info: "45% severe · 15% tornado" },
  { label: "HIGH", name: "High",            color: "#FF00FF", fill: "#FF99FF",
    info: "60%+ severe · 30%+ tornado" },
];

/* ── SPC probabilistic / conditional intensity outlook styling ── */
const PROB_COLORS = {
  torn: [
    { label: "0.02", pct: "2%",  color: "#008B00", fill: "#66CC66" },
    { label: "0.05", pct: "5%",  color: "#8B4513", fill: "#CD853F" },
    { label: "0.10", pct: "10%", color: "#FFD700", fill: "#FFEE88" },
    { label: "0.15", pct: "15%", color: "#FF0000", fill: "#FF6666" },
    { label: "0.30", pct: "30%", color: "#FF00FF", fill: "#FF88FF" },
    { label: "0.45", pct: "45%", color: "#9400D3", fill: "#C488FF" },
    { label: "0.60", pct: "60%", color: "#000080", fill: "#6666CC" },
    { label: "SIGN", pct: "Sig", color: "#000000", fill: "#000000" },
    { label: "CIG1", pct: "CIG1", color: "#005500", fill: "#66A366" },
    { label: "CIG2", pct: "CIG2", color: "#FF6600", fill: "#FFA366" },
    { label: "CIG3", pct: "CIG3", color: "#CC0000", fill: "#E06666" },
  ],
  wind: [
    { label: "0.05", pct: "5%",  color: "#8B4513", fill: "#CD853F" },
    { label: "0.15", pct: "15%", color: "#FFD700", fill: "#FFEE88" },
    { label: "0.30", pct: "30%", color: "#FF0000", fill: "#FF6666" },
    { label: "0.45", pct: "45%", color: "#FF00FF", fill: "#FF88FF" },
    { label: "0.60", pct: "60%", color: "#9400D3", fill: "#C488FF" },
    { label: "SIGN", pct: "Sig", color: "#000000", fill: "#000000" },
    { label: "CIG1", pct: "CIG1", color: "#005500", fill: "#66A366" },
    { label: "CIG2", pct: "CIG2", color: "#FF6600", fill: "#FFA366" },
    { label: "CIG3", pct: "CIG3", color: "#CC0000", fill: "#E06666" },
  ],
  hail: [
    { label: "0.05", pct: "5%",  color: "#8B4513", fill: "#CD853F" },
    { label: "0.15", pct: "15%", color: "#FFD700", fill: "#FFEE88" },
    { label: "0.30", pct: "30%", color: "#FF0000", fill: "#FF6666" },
    { label: "0.45", pct: "45%", color: "#FF00FF", fill: "#FF88FF" },
    { label: "0.60", pct: "60%", color: "#9400D3", fill: "#C488FF" },
    { label: "SIGN", pct: "Sig", color: "#000000", fill: "#000000" },
    { label: "CIG1", pct: "CIG1", color: "#005500", fill: "#66A366" },
    { label: "CIG2", pct: "CIG2", color: "#FF6600", fill: "#FFA366" },
  ],
  /* Day 4-8 extended-range probability contours */
  prob: [
    { label: "0.15", pct: "15%", color: "#DDB600", fill: "#FFE066",
      info: "15% probability of severe weather within 25 mi" },
    { label: "0.30", pct: "30%", color: "#FF6600", fill: "#FFA366",
      info: "30% probability of severe weather within 25 mi" },
    { label: "CATEGORICAL", pct: "Cat", color: "#DDB600", fill: "#FFE066",
      info: "General severe weather area" },
  ],
};

const OUTLOOK_TYPES = [
  { id: "cat",  name: "Categorical" },
  { id: "torn", name: "Tornado", title: "Probabilistic Tornado Outlook" },
  { id: "wind", name: "Wind",    title: "Probabilistic Wind Outlook" },
  { id: "hail", name: "Hail",    title: "Probabilistic Hail Outlook" },
];

function spcStyle(feature, outlookType) {
  const label = feature.properties?.LABEL || feature.properties?.LABEL2 || "";
  if (outlookType === "cat" || !outlookType) {
    const cat = SPC_CATEGORIES.find((c) => c.label === label);
    return {
      color: cat?.color || feature.properties?.stroke || "#888",
      fillColor: cat?.fill || feature.properties?.fill || "#888",
      fillOpacity: 0.25,
      weight: 1.5,
    };
  }
  // Probabilistic layers
  const probList = PROB_COLORS[outlookType] || [];
  const prob = probList.find((p) => p.label === label);
  if (label === "SIGN") {
    return {
      color: "#111",
      fillColor: "transparent",
      fillOpacity: 0,
      weight: 2.5,
      dashArray: "8 4",
    };
  }
  return {
    color: prob?.color || feature.properties?.stroke || "#888",
    fillColor: prob?.fill || feature.properties?.fill || "#888",
    fillOpacity: 0.3,
    weight: 1.5,
  };
}

/* ── GeoJSON overlay in a non-interactive pane below markers ── */
function OutlookLayer({ data, outlookType }) {
  if (!data || !data.features || data.features.length === 0) return null;
  return (
    <Pane name="spc-outlook" style={{ zIndex: 250, pointerEvents: "none" }}>
      <GeoJSON
        key={JSON.stringify(data).slice(0, 100) + outlookType}
        data={data}
        style={(f) => spcStyle(f, outlookType)}
        interactive={false}
      />
    </Pane>
  );
}

/* ── click-for-coords layer ─────────────────────────────────── */
function MapClickHandler({ enabled, onLatLonSelect, onRightClick, onMapClick, onMouseMove, disableRightClick = false }) {
  useMapEvents({
    click(e) {
      if (onMapClick && onMapClick(e.latlng.lat, e.latlng.lng, e) === true) return;
      if (enabled && onLatLonSelect) {
        onLatLonSelect(e.latlng.lat.toFixed(2), e.latlng.lng.toFixed(2));
      }
    },
    mousemove(e) {
      onMouseMove?.(e.latlng.lat, e.latlng.lng);
    },
    contextmenu(e) {
      if (disableRightClick) {
        e.originalEvent.preventDefault();
        return;
      }
      // Right-click → proximity search
      if (onRightClick) {
        e.originalEvent.preventDefault();
        onRightClick(e.latlng.lat, e.latlng.lng);
      }
    },
  });
  return null;
}

/* ── fly-to when selected station changes ───────────────────── */
function FlyToStation({ station, stations }) {
  const map = useMap();
  useEffect(() => {
    if (!station) return;
    const s = stations.find((st) => st.id === station);
    if (s) map.flyTo([s.lat, s.lon], 7, { duration: 0.6 });
  }, [station]); // eslint-disable-line react-hooks/exhaustive-deps
  return null;
}

/* ── Invalidate map size when container changes ────────────── */
function MapResizeHandler({ trigger }) {
  const map = useMap();
  useEffect(() => {
    // Small delay lets the CSS transition finish before recalculating
    const id = setTimeout(() => map.invalidateSize(), 350);
    return () => clearTimeout(id);
  }, [trigger, map]);
  return null;
}

/* ── main component ─────────────────────────────────────────── */
export default function StationMap({
  stations = [],
  riskData,
  selectedStation,
  onStationSelect,
  onLatLonSelect,
  onStormMotionSelect,
  onBoundaryOrientationSelect,
  latLonMode,   // true when source is rap
  onClose,
  onFetchLatest,            // (stationId) => void  — triggers auto-fetch with best source
  onCompareStations,        // (stationIds[]) => void — compare nearest stations
}) {
  const [isFullscreen, setIsFullscreen] = useState(false);
  const [utcNow, setUtcNow] = useState(() => new Date());
  const [outlookDay, setOutlookDay] = useState(1);
  const [outlookType, setOutlookType] = useState("cat");
  const [outlookData, setOutlookData] = useState(null);
  const [outlookLoading, setOutlookLoading] = useState(false);
  const [showOutlook, setShowOutlook] = useState(true);
  const [showRadar, setShowRadar] = useState(false);
  const [radarSource, setRadarSource] = useState("mosaic"); // "composite" | "mosaic" | "singlesite"
  const [ssRadarBounds, setSsRadarBounds] = useState(null); // [[south,west],[north,east]] for single-site ImageOverlay
  const [showVelocity, setShowVelocity] = useState(false);
  const [velocityRadar, setVelocityRadar] = useState({ id: "TLX", lat: 35.33, lon: -97.28 });   // nearest NEXRAD
  const [velCacheBust, setVelCacheBust] = useState(() => Date.now());
  const [velScanTime, setVelScanTime] = useState(null);        // latest volume scan UTC string
  const [velProduct, setVelProduct] = useState("N0S");          // active velocity product
  const [velAvailProducts, setVelAvailProducts] = useState([]);  // IEM available products for current radar

  // Animated velocity state
  const [velFrames, setVelFrames] = useState([]);     // [{url, time, ts}]
  const [velBounds, setVelBounds] = useState(null);   // [[south,west],[north,east]]
  const [velFrame, setVelFrame] = useState(0);
  const [velPlaying, setVelPlaying] = useState(false);
  const velTimerRef = useRef(null);
  const velFrameCount = velFrames.length;
  const [baseMap, setBaseMap] = useState("dark");
  const [showBaseMapPicker, setShowBaseMapPicker] = useState(false);

  // NWS active warnings overlay
  const [showWarnings, setShowWarnings] = useState(false);
  const [warningsData, setWarningsData] = useState(null);
  const [warningsLoading, setWarningsLoading] = useState(false);
  const warningCount = warningsData?.features?.length || 0;

  // TVS / Mesocyclone storm attributes overlay
  const [showStorms, setShowStorms] = useState(false);
  const [stormData, setStormData] = useState(null);
  const [stormsLoading, setStormsLoading] = useState(false);
  const tvsCount = stormData?.features?.filter((f) => f.properties?.tvs === "TVS").length || 0;
  const mesoCount = stormData?.features?.length || 0;

  // Spotter Network live chasers
  const [showSpotters, setShowSpotters] = useState(false);
  const [spotterData, setSpotterData] = useState(null);    // {positions, reports}
  const [spottersLoading, setSpottersLoading] = useState(false);
  const spotterCount = spotterData?.positions?.length || 0;
  const spotterReportCount = spotterData?.reports?.length || 0;

  // Animated wind flow overlay
  const [showWind, setShowWind] = useState(false);
  const [windLevel, setWindLevel] = useState("500");  // "500" = steering, "surface" = 10m
  const [windHourOffset, setWindHourOffset] = useState(0);
  const [windData, setWindData] = useState(null);
  const [syncTimeline, setSyncTimeline] = useState(true);
  const windFetchDebounceRef = useRef(null);
  const windReqSeqRef = useRef(0);

  // SPC Mesoscale Discussions & Watches overlay
  const [showMdWatch, setShowMdWatch] = useState(false);
  const [mdData, setMdData] = useState(null);
  const [watchData, setWatchData] = useState(null);
  const [mdWatchLoading, setMdWatchLoading] = useState(false);
  const mdCount = mdData?.features?.length || 0;
  const watchCount = watchData?.features?.length || 0;

  // Lightning overlay (Blitzortung WebSocket)
  const [showLightning, setShowLightning] = useState(false);
  const lightningStrikes = useLightningData(showLightning);

  // SPC outlook refresh counter
  const [outlookRefresh, setOutlookRefresh] = useState(0);

  // Toolbar group dropdown (only one open at a time)
  const [openGroup, setOpenGroup] = useState(null); // "forecasts" | "radar" | "overlays" | "tools" | null
  const toggleGroup = useCallback((g) => setOpenGroup((prev) => (prev === g ? null : g)), []);

  // Draw tools: map-driven storm motion and boundary orientation
  const [drawMode, setDrawMode] = useState(null); // "stormMotion" | "boundary" | null
  const [drawPoints, setDrawPoints] = useState([]); // [[lat, lon]]
  const [stormMotionGuide, setStormMotionGuide] = useState(null); // {start, end, direction, speed, toward}
  const [boundaryGuide, setBoundaryGuide] = useState(null); // {start, end, orientation}
  const [drawPreview, setDrawPreview] = useState(null); // live preview while cursor moves after first click

  // Outlook stations (stations within SPC outlook polygons)
  const [outlookStations, setOutlookStations] = useState(null);
  const [outlookStationsLoading, setOutlookStationsLoading] = useState(false);
  const [showOutlookStations, setShowOutlookStations] = useState(false);
  const [outlookError, setOutlookError] = useState(null);
  const [popupStationId, setPopupStationId] = useState(null);
  const [showDiscussion, setShowDiscussion] = useState(false);
  const [discussionText, setDiscussionText] = useState("");
  const [discussionLoading, setDiscussionLoading] = useState(false);

  // ── Proximity search (right-click map) ────────────────────
  const [proximityResult, setProximityResult] = useState(null); // {lat, lon, stations: [{id, name, dist}]}

  // ── Wind barbs from risk-scan data on map ─────────────────
  const [showWindBarbs, setShowWindBarbs] = useState(false);

  // ── Favorite stations as distinct markers ─────────────────
  const [showFavorites, setShowFavorites] = useState(true);
  const [favorites, setFavorites] = useState(() => getFavorites());
  // Refresh favorites when window regains focus (in case changed in ControlPanel)
  useEffect(() => {
    const handler = () => setFavorites(getFavorites());
    window.addEventListener("focus", handler);
    return () => window.removeEventListener("focus", handler);
  }, []);

  // ── ACARS airports overlay ────────────────────────────────
  const [showAcars, setShowAcars] = useState(false);
  const [acarsAirports, setAcarsAirports] = useState(null);
  const [acarsLoading, setAcarsLoading] = useState(false);

  // Animated radar state — unified for both composite (RainViewer) and mosaic (IEM)
  const [radarPlaying, setRadarPlaying] = useState(false);
  const [radarFrame, setRadarFrame] = useState(0);
  const [radarFrames, setRadarFrames] = useState([]);   // composite: [{time, path, forecast}]  mosaic: [{time, ts}]  singlesite: [{url, time, ts}]
  const radarTimerRef = useRef(null);
  const radarFrameCount = radarFrames.length;

  // Fetch radar frames — RainViewer for composite, generated timestamps for mosaic, IEM archive for singlesite
  useEffect(() => {
    if (!showRadar) return;
    let cancelled = false;
    // Clear stale frames immediately so the wrong source type isn't rendered
    setRadarFrames([]);
    setRadarFrame(0);
    setRadarPlaying(false);
    setSsRadarBounds(null);
    if (radarSource === "composite") {
      const load = async () => {
        const { past, nowcast } = await fetchRadarFrames();
        if (cancelled) return;
        const frames = [
          ...past.map((f) => ({ ...f, forecast: false })),
          ...nowcast.map((f) => ({ ...f, forecast: true })),
        ];
        setRadarFrames(frames);
        setRadarFrame(Math.max(0, past.length - 1));
      };
      load();
      const id = setInterval(load, 300_000);
      return () => { cancelled = true; clearInterval(id); };
    } else if (radarSource === "singlesite") {
      // Single-site via NWS WMS — no frames needed, layer is always live
      return;
    } else {
      // Mosaic: build frames from timestamps (last 2 hrs, 5-min intervals)
      const frames = buildMosaicFrames(24);
      setRadarFrames(frames);
      setRadarFrame(frames.length - 1); // default to latest
      // Rebuild every 5 min to pick up new data
      const id = setInterval(() => {
        if (cancelled) return;
        const f = buildMosaicFrames(24);
        setRadarFrames(f);
        // Keep at latest if user hasn't scrubbed
        setRadarFrame((prev) => prev >= f.length - 2 ? f.length - 1 : prev);
      }, 300_000);
      return () => { cancelled = true; clearInterval(id); };
    }
  }, [showRadar, radarSource, velocityRadar.id]);

  // Advance radar frame
  useEffect(() => {
    if (!radarPlaying || !showRadar || radarFrameCount === 0) {
      clearInterval(radarTimerRef.current);
      return;
    }
    radarTimerRef.current = setInterval(() => {
      setRadarFrame((f) => (f + 1) % radarFrameCount);
    }, RADAR_INTERVAL_MS);
    return () => clearInterval(radarTimerRef.current);
  }, [radarPlaying, showRadar, radarFrameCount]);

  // Reset frame when radar is toggled off
  useEffect(() => {
    if (!showRadar) { setRadarPlaying(false); setRadarFrame(0); setRadarFrames([]); }
  }, [showRadar]);

  // Auto-refresh velocity tiles every 2 min (NEXRAD volume scan ~4-5 min)
  useEffect(() => {
    if (!showVelocity) return;
    const id = setInterval(() => setVelCacheBust(Date.now()), 120_000);
    return () => clearInterval(id);
  }, [showVelocity]);

  // Fetch velocity animation frames from IEM archive
  useEffect(() => {
    if (!showVelocity) {
      setVelScanTime(null); setVelFrames([]); setVelFrame(0); setVelPlaying(false);
      return;
    }
    let cancelled = false;
    setVelAvailProducts(VELOCITY_PRODUCTS);
    const load = async () => {
      try {
        const { frames, bounds } = await fetchVelFrames(velocityRadar.id, velProduct);
        if (cancelled) return;
        setVelFrames(frames);
        setVelBounds(bounds);
        setVelFrame(Math.max(0, frames.length - 1)); // default to latest
        if (frames.length) {
          const last = frames[frames.length - 1];
          setVelScanTime(new Date(last.time * 1000).toISOString().slice(0, 16).replace("T", " ") + "Z");
        } else {
          setVelScanTime(null);
        }
      } catch (e) {
        console.warn("IEM velocity frame error:", e);
        if (!cancelled) setVelScanTime(null);
      }
    };
    load();
    const id = setInterval(load, 120_000); // refresh every 2 min
    return () => { cancelled = true; clearInterval(id); };
  }, [showVelocity, velocityRadar.id, velProduct, velCacheBust]);

  // Advance velocity animation frame
  useEffect(() => {
    if (!velPlaying || !showVelocity || velFrameCount === 0) {
      clearInterval(velTimerRef.current);
      return;
    }
    velTimerRef.current = setInterval(() => {
      setVelFrame((f) => (f + 1) % velFrameCount);
    }, VEL_ANIM_INTERVAL_MS);
    return () => clearInterval(velTimerRef.current);
  }, [velPlaying, showVelocity, velFrameCount]);

  // Optional timeline sync: map radar/velocity frame position to wind forecast offset.
  useEffect(() => {
    if (!syncTimeline || !showWind) return;
    let nextOffset = null;
    if (showVelocity && velFrameCount > 0) {
      nextOffset = frameToHourOffset(velFrame, velFrameCount);
    } else if (showRadar && radarFrameCount > 0) {
      nextOffset = frameToHourOffset(radarFrame, radarFrameCount);
    }
    if (nextOffset == null) return;
    setWindHourOffset((prev) => (prev === nextOffset ? prev : nextOffset));
  }, [
    syncTimeline,
    showWind,
    showVelocity,
    velFrame,
    velFrameCount,
    showRadar,
    radarFrame,
    radarFrameCount,
  ]);

  const handleCenterChange = useCallback((lat, lng) => {
    const nr = nearestNexrad(lat, lng);
    setVelocityRadar((prev) => prev.id === nr.id ? prev : nr);
  }, []);

  const activeTimeline = useMemo(() => {
    if (showVelocity && velFrameCount > 0) {
      return {
        kind: "velocity",
        frame: velFrame,
        count: velFrameCount,
        playing: velPlaying,
        label: velFrames[velFrame] ? fmtRadarTime(velFrames[velFrame].time) : "Animate",
      };
    }
    if (showRadar && radarFrameCount > 0) {
      return {
        kind: "radar",
        frame: radarFrame,
        count: radarFrameCount,
        playing: radarPlaying,
        label: radarFrames[radarFrame]
          ? `${radarFrames[radarFrame].forecast ? "\u2601 " : ""}${fmtRadarTime(radarFrames[radarFrame].time)}`
          : "Animate",
      };
    }
    return null;
  }, [
    showVelocity, velFrameCount, velFrame, velPlaying, velFrames,
    showRadar, radarFrameCount, radarFrame, radarPlaying, radarFrames,
  ]);

  const toggleTimelinePlayback = useCallback(() => {
    if (!activeTimeline) return;
    if (activeTimeline.kind === "velocity") {
      setVelPlaying((v) => !v);
    } else {
      setRadarPlaying((v) => !v);
    }
  }, [activeTimeline]);

  const scrubTimelineFrame = useCallback((rawFrame) => {
    const frame = Number(rawFrame);
    if (showVelocity && velFrameCount > 0) {
      const idx = clampNum(frame, 0, velFrameCount - 1);
      setVelFrame(idx);
      setVelPlaying(false);
      if (syncTimeline && showWind) {
        setWindHourOffset(frameToHourOffset(idx, velFrameCount));
      }
      return;
    }
    if (showRadar && radarFrameCount > 0) {
      const idx = clampNum(frame, 0, radarFrameCount - 1);
      setRadarFrame(idx);
      setRadarPlaying(false);
      if (syncTimeline && showWind) {
        setWindHourOffset(frameToHourOffset(idx, radarFrameCount));
      }
    }
  }, [
    showVelocity, velFrameCount, syncTimeline, showWind,
    showRadar, radarFrameCount,
  ]);

  // ── Proximity search: right-click map → find nearest N stations ──
  const handleProximitySearch = useCallback((lat, lng) => {
    if (!stations.length) return;
    const nearest = stations
      .map((s) => {
        const dLat = s.lat - lat;
        const dLon = (s.lon - lng) * Math.cos((lat * Math.PI) / 180);
        const dist = Math.sqrt(dLat * dLat + dLon * dLon) * 111.32; // approx km
        return { id: s.id, name: s.name, lat: s.lat, lon: s.lon, dist };
      })
      .sort((a, b) => a.dist - b.dist)
      .slice(0, 5);
    setProximityResult({ lat, lon: lng, stations: nearest });
  }, [stations]);

  const toggleDrawMode = useCallback((mode) => {
    setDrawPoints([]);
    setDrawPreview(null);
    setDrawMode((prev) => (prev === mode ? null : mode));
  }, []);

  const clearDrawGuides = useCallback(() => {
    setDrawMode(null);
    setDrawPoints([]);
    setDrawPreview(null);
    setStormMotionGuide(null);
    setBoundaryGuide(null);
  }, []);

  const handleMapToolClick = useCallback((lat, lon) => {
    if (!drawMode) return false;

    const clickPoint = [Number(lat), Number(lon)];

    // First click picks the start point; second click completes the vector/line.
    if (drawPoints.length === 0) {
      setDrawPoints([clickPoint]);
      return true;
    }

    const [start] = drawPoints;
    const [sLat, sLon] = start;
    const [eLat, eLon] = clickPoint;
    const distKm = haversineKm(sLat, sLon, eLat, eLon);

    // Ignore accidental near-zero second clicks.
    if (distKm < 0.5) return true;

    if (drawMode === "stormMotion") {
      // User draws where the storm is moving toward.
      const toward = bearingDeg(sLat, sLon, eLat, eLon);
      // Backend expects meteorological "from" direction for wind_components.
      const direction = (toward + 180) % 360;
      // Interpret line length as 1-hour displacement (nm/hr = kt).
      const speed = Math.max(5, Math.min(120, Math.round(distKm * 0.539957)));
      const payload = {
        start,
        end: clickPoint,
        direction: Math.round(direction),
        speed,
        toward: Math.round(toward),
      };
      setStormMotionGuide(payload);
      onStormMotionSelect?.({ direction: payload.direction, speed: payload.speed });
    } else if (drawMode === "boundary") {
      const rawOrient = bearingDeg(sLat, sLon, eLat, eLon);
      // Boundary orientation is an axis, so 045° and 225° are equivalent.
      const orientation = Math.round(((rawOrient % 180) + 180) % 180);
      const payload = {
        start,
        end: clickPoint,
        orientation,
      };
      setBoundaryGuide(payload);
      onBoundaryOrientationSelect?.(orientation);
    }

    setDrawMode(null);
    setDrawPoints([]);
    setDrawPreview(null);
    return true;
  }, [drawMode, drawPoints, onStormMotionSelect, onBoundaryOrientationSelect]);

  const handleDrawMouseMove = useCallback((lat, lon) => {
    if (!drawMode || drawPoints.length !== 1) { setDrawPreview(null); return; }
    const [sLat, sLon] = drawPoints[0];
    const eLat = Number(lat), eLon = Number(lon);
    const distKm = haversineKm(sLat, sLon, eLat, eLon);
    if (distKm < 0.5) { setDrawPreview(null); return; }
    const cursorPos = [eLat, eLon];
    if (drawMode === "stormMotion") {
      const toward = bearingDeg(sLat, sLon, eLat, eLon);
      const direction = (toward + 180) % 360;
      const speed = Math.max(5, Math.min(120, Math.round(distKm * 0.539957)));
      setDrawPreview({ cursorPos, direction: Math.round(direction), speed, toward: Math.round(toward) });
    } else if (drawMode === "boundary") {
      const rawOrient = bearingDeg(sLat, sLon, eLat, eLon);
      const orientation = Math.round(((rawOrient % 180) + 180) % 180);
      setDrawPreview({ cursorPos, orientation });
    }
  }, [drawMode, drawPoints]);

  // ── Fetch ACARS airports when toggled on ─────────────────────
  useEffect(() => {
    if (!showAcars) { setAcarsAirports(null); return; }
    let cancelled = false;
    setAcarsLoading(true);
    fetchAcarsAirports()
      .then((data) => { if (!cancelled) setAcarsAirports(data); })
      .catch((e) => console.error("ACARS airports fetch error:", e))
      .finally(() => { if (!cancelled) setAcarsLoading(false); });
    return () => { cancelled = true; };
  }, [showAcars]);

  // Fetch NWS active warnings on mount + every 2 min when toggle is on
  useEffect(() => {
    if (!showWarnings) { setWarningsData(null); return; }
    let cancelled = false;
    const load = async () => {
      setWarningsLoading(true);
      const data = await fetchNwsWarnings();
      if (!cancelled) { setWarningsData(data); setWarningsLoading(false); }
    };
    load();
    const id = setInterval(load, 120_000); // refresh every 2 min
    return () => { cancelled = true; clearInterval(id); };
  }, [showWarnings]);

  // Fetch Spotter Network positions + reports every 60s when toggle is on
  useEffect(() => {
    if (!showSpotters) { setSpotterData(null); return; }
    let cancelled = false;
    const load = async () => {
      setSpottersLoading(true);
      const data = await fetchSpotterNetwork();
      if (!cancelled) { setSpotterData(data); setSpottersLoading(false); }
    };
    load();
    const id = setInterval(load, 60_000);
    return () => { cancelled = true; clearInterval(id); };
  }, [showSpotters]);

  // Fetch storm attributes (TVS/Meso) every 60s when toggle is on
  useEffect(() => {
    if (!showStorms) { setStormData(null); return; }
    let cancelled = false;
    const load = async () => {
      setStormsLoading(true);
      const data = await fetchStormAttr();
      if (!cancelled) { setStormData(data); setStormsLoading(false); }
    };
    load();
    const id = setInterval(load, 60_000); // refresh every 1 min
    return () => { cancelled = true; clearInterval(id); };
  }, [showStorms]);

  // Fetch SPC MDs & Watches every 2 min when toggle is on
  useEffect(() => {
    if (!showMdWatch) { setMdData(null); setWatchData(null); return; }
    let cancelled = false;
    const load = async () => {
      setMdWatchLoading(true);
      const [mds, watches] = await Promise.all([fetchSpcMds(), fetchSpcWatches()]);
      if (!cancelled) { setMdData(mds); setWatchData(watches); setMdWatchLoading(false); }
    };
    load();
    const id = setInterval(load, 120_000);
    return () => { cancelled = true; clearInterval(id); };
  }, [showMdWatch]);

  // Fetch wind field with scrub-friendly debounce to avoid UI jitter.
  useEffect(() => {
    clearTimeout(windFetchDebounceRef.current);

    if (!showWind) {
      return;
    }

    let cancelled = false;
    const reqSeq = ++windReqSeqRef.current;

    windFetchDebounceRef.current = setTimeout(() => {
      fetchWindField(windLevel, windHourOffset)
        .then((d) => {
          if (!cancelled && reqSeq === windReqSeqRef.current) setWindData(d);
        })
        .catch((e) => {
          if (!cancelled && reqSeq === windReqSeqRef.current) {
            console.error("Wind field error:", e);
          }
        });
    }, 120);

    return () => {
      cancelled = true;
      clearTimeout(windFetchDebounceRef.current);
    };
  }, [showWind, windLevel, windHourOffset]);

  // Fetch SPC outlook on mount and when day/type changes
  useEffect(() => {
    if (!showOutlook) return;
    // Day 3 only has categorical (no prob outlooks)
    // Day 4-8 only has extended-range prob
    let effectiveT = outlookType;
    if (outlookDay >= 4) {
      effectiveT = "prob";
    } else if (outlookDay === 3 && ["torn", "wind", "hail"].includes(outlookType)) {
      effectiveT = "cat";
    }
    let cancelled = false;
    setOutlookLoading(true);
    setOutlookError(null);
    fetchSpcOutlook(outlookDay, effectiveT)
      .then((data) => { if (!cancelled) { setOutlookData(data); setOutlookError(null); } })
      .catch(() => {
        if (!cancelled) {
          setOutlookData(null);
          setOutlookError(null);
        }
      })
      .finally(() => { if (!cancelled) setOutlookLoading(false); });
    return () => { cancelled = true; };
  }, [outlookDay, outlookType, showOutlook, outlookRefresh]);

  // Fetch stations within outlook when toggled
  useEffect(() => {
    if (!showOutlookStations || !showOutlook) {
      setOutlookStations(null);
      return;
    }
    let effectiveT = outlookType;
    if (outlookDay >= 4) {
      effectiveT = "prob";
    } else if (outlookDay === 3 && ["torn", "wind", "hail"].includes(outlookType)) {
      effectiveT = "cat";
    }
    let cancelled = false;
    setOutlookStationsLoading(true);
    fetchSpcOutlookStations(outlookDay, effectiveT)
      .then((data) => { if (!cancelled) setOutlookStations(data); })
      .catch(() => { if (!cancelled) setOutlookStations(null); })
      .finally(() => { if (!cancelled) setOutlookStationsLoading(false); });
    return () => { cancelled = true; };
  }, [showOutlookStations, showOutlook, outlookDay, outlookType, outlookRefresh]);

  // Determine which SPC categories are present in the current data
  const effectiveType = (() => {
    if (outlookDay >= 4) return "prob";
    if (outlookDay === 3 && ["torn", "wind", "hail"].includes(outlookType)) return "cat";
    return outlookType;
  })();
  const activeCategories = useMemo(() => {
    if (!outlookData?.features) return [];
    const labels = new Set(outlookData.features.map((f) => f.properties?.LABEL || f.properties?.LABEL2));
    if (effectiveType === "cat") {
      return SPC_CATEGORIES.filter((c) => labels.has(c.label));
    }
    const probList = PROB_COLORS[effectiveType] || [];
    return probList.filter((p) => labels.has(p.label));
  }, [outlookData, effectiveType]);

  // Get outlook metadata (valid/expire/forecaster)
  const outlookMeta = useMemo(() => {
    if (!outlookData?.features?.length) return null;
    const p = outlookData.features[0].properties;
    return {
      valid: p.VALID_ISO ? new Date(p.VALID_ISO).toUTCString().slice(0, 22) : "",
      expire: p.EXPIRE_ISO ? new Date(p.EXPIRE_ISO).toUTCString().slice(0, 22) : "",
      forecaster: p.FORECASTER || "",
    };
  }, [outlookData]);

  // Fetch SPC discussion when toggled
  useEffect(() => {
    if (!showDiscussion || !showOutlook) { setDiscussionText(""); return; }
    let cancelled = false;
    setDiscussionLoading(true);
    fetchSpcDiscussion(outlookDay)
      .then((data) => { if (!cancelled) setDiscussionText(data.text || ""); })
      .catch(() => { if (!cancelled) setDiscussionText("Failed to load discussion."); })
      .finally(() => { if (!cancelled) setDiscussionLoading(false); });
    return () => { cancelled = true; };
  }, [showDiscussion, showOutlook, outlookDay, outlookRefresh]);

  // Set of station IDs within the outlook area
  const outlookStationIds = useMemo(() => {
    if (!outlookStations?.stations) return new Set();
    return new Set(outlookStations.stations.map((s) => s.id));
  }, [outlookStations]);

  const riskMap = useMemo(() => {
    const m = {};
    if (riskData) {
      riskData.stations.forEach((s) => {
        m[s.id] = s;
      });
    }
    return m;
  }, [riskData]);

  const CONUS_CENTER = [39.0, -96.0];
  const CONUS_ZOOM = 4;

  const mapCenter = CONUS_CENTER;
  const mapZoom = CONUS_ZOOM;
  const [mapView, setMapView] = useState({
    zoom: CONUS_ZOOM,
    south: 22,
    west: -128,
    north: 52,
    east: -66,
  });

  const stationRender = useMemo(() => {
    if (!stations.length) {
      return { stations: [], clusters: [], hiddenCount: 0 };
    }

    const view = mapView || { zoom: 4, south: 20, west: -130, north: 55, east: -60 };
    const zoom = view.zoom ?? 4;
    const pad = zoom <= 4 ? 2.5 : zoom <= 6 ? 1.5 : 0.8;
    const south = view.south - pad;
    const north = view.north + pad;
    const west = view.west;
    const east = view.east;

    const lonInBounds = (lon) => (
      west <= east ? (lon >= west && lon <= east) : (lon >= west || lon <= east)
    );
    const inBounds = (s) => (
      s.lat >= south && s.lat <= north && lonInBounds(s.lon)
    );

    const priorityIds = new Set();
    if (selectedStation) priorityIds.add(selectedStation);
    if (popupStationId) priorityIds.add(popupStationId);
    favorites.forEach((fid) => priorityIds.add(fid));
    if (showOutlookStations) {
      outlookStationIds.forEach((sid) => priorityIds.add(sid));
    }

    const candidates = stations.filter((s) => inBounds(s) || priorityIds.has(s.id));
    if (zoom >= 7) {
      return { stations: candidates, clusters: [], hiddenCount: 0 };
    }

    const cellDeg = zoom <= 3 ? 2.4 : zoom <= 4 ? 1.8 : zoom <= 5 ? 1.2 : 0.9;
    const buckets = new Map();
    const kept = new Map();

    for (const s of candidates) {
      if (priorityIds.has(s.id)) {
        kept.set(s.id, s);
        continue;
      }

      const wrappedLon = ((s.lon + 540) % 360) - 180;
      const key = `${Math.floor((s.lat + 90) / cellDeg)}_${Math.floor((wrappedLon + 180) / cellDeg)}`;
      let bucket = buckets.get(key);
      if (!bucket) {
        bucket = { members: [], sumLat: 0, sumLon: 0, maxStp: null };
        buckets.set(key, bucket);
      }
      bucket.members.push(s);
      bucket.sumLat += s.lat;
      bucket.sumLon += s.lon;
      const stp = riskMap[s.id]?.stp;
      if (stp != null && (bucket.maxStp == null || stp > bucket.maxStp)) {
        bucket.maxStp = stp;
      }
    }

    const renderedStations = [...kept.values()];
    const clusters = [];
    let hiddenCount = 0;

    for (const [key, bucket] of buckets.entries()) {
      if (bucket.members.length === 1) {
        renderedStations.push(bucket.members[0]);
        continue;
      }
      hiddenCount += bucket.members.length;
      clusters.push({
        id: key,
        lat: bucket.sumLat / bucket.members.length,
        lon: bucket.sumLon / bucket.members.length,
        count: bucket.members.length,
        color: riskColor(bucket.maxStp),
      });
    }

    return { stations: renderedStations, clusters, hiddenCount };
  }, [
    stations,
    mapView,
    selectedStation,
    popupStationId,
    favorites,
    showOutlookStations,
    outlookStationIds,
    riskMap,
  ]);

  // Live UTC clock
  useEffect(() => {
    const id = setInterval(() => setUtcNow(new Date()), 1000);
    return () => clearInterval(id);
  }, []);

  // Escape key exits fullscreen
  useEffect(() => {
    if (!isFullscreen) return;
    const handler = (e) => { if (e.key === "Escape") setIsFullscreen(false); };
    window.addEventListener("keydown", handler);
    return () => window.removeEventListener("keydown", handler);
  }, [isFullscreen]);

  // Toggle body scroll lock when fullscreen
  useEffect(() => {
    if (isFullscreen) {
      document.body.classList.add("smap-fullscreen-active");
    } else {
      document.body.classList.remove("smap-fullscreen-active");
    }
    return () => document.body.classList.remove("smap-fullscreen-active");
  }, [isFullscreen]);

  return (
    <section className={`station-map-panel${isFullscreen ? " smap-fullscreen" : ""}`}>
      {/* Two-tier toolbar matching rv-meta-bar layout */}
      <div className="smap-toolbar">
        <div className="smap-toolbar-bottom">
          {/* ── Group header buttons + controls ── */}
          <div className="smap-groups">
            {latLonMode && <span className="smap-tag">Click map to set coordinates</span>}
            <button
              className={`smap-group-btn${openGroup === "forecasts" ? " open" : ""}${showOutlook ? " has-active" : ""}`}
              onClick={() => toggleGroup("forecasts")}
            >
              <CloudLightning size={11} />
              Forecasts
              {showOutlook && <span className="smap-group-dot" />}
              {showOutlookStations && outlookStations && <span className="smap-radar-badge">{outlookStations.count}</span>}
              <ChevronDown size={10} className="smap-group-chevron" />
            </button>
            <button
              className={`smap-group-btn${openGroup === "radar" ? " open" : ""}${showRadar || showVelocity ? " has-active" : ""}`}
              onClick={() => toggleGroup("radar")}
            >
              <Crosshair size={11} />
              Radar
              {showRadar && <span className="smap-group-dot" />}
              {showVelocity && <span className="smap-group-dot smap-dot-vel" />}
              <ChevronDown size={10} className="smap-group-chevron" />
            </button>
            <button
              className={`smap-group-btn${openGroup === "overlays" ? " open" : ""}${(showStorms || showWarnings || showMdWatch || showSpotters || showLightning || showWind) ? " has-active" : ""}`}
              onClick={() => toggleGroup("overlays")}
            >
              <Layers size={11} />
              Overlays
              {warningCount > 0 && <span className="smap-radar-badge smap-warn-count">{warningCount}</span>}
              {tvsCount > 0 && <span className="smap-radar-badge smap-warn-count">{tvsCount} TVS</span>}
              {lightningStrikes.length > 0 && <span className="smap-radar-badge">{lightningStrikes.length}</span>}
              <ChevronDown size={10} className="smap-group-chevron" />
            </button>
            <div className="smap-toolbar-end">
              <button
                className={`smap-group-btn${openGroup === "tools" ? " open" : ""}${(stormMotionGuide || boundaryGuide || drawMode) ? " has-active" : ""}`}
                onClick={() => toggleGroup("tools")}
              >
                <Compass size={11} />
                Tools
                {stormMotionGuide && <span className="smap-group-dot" />}
                {boundaryGuide && <span className="smap-group-dot smap-dot-vel" />}
                <ChevronDown size={10} className="smap-group-chevron" />
              </button>
              <span className="smap-utc-clock" title="Current UTC / Zulu time">
                {utcNow.toISOString().slice(11, 16)}Z
              </span>
              <button
                className={`smap-expand${isFullscreen ? " active" : ""}`}
                onClick={() => setIsFullscreen((v) => !v)}
                title={isFullscreen ? "Exit fullscreen (Esc)" : "Expand map fullscreen"}
              >
                {isFullscreen ? <Minimize2 size={14} /> : <Maximize2 size={14} />}
              </button>
              <button className="smap-expand" onClick={() => { setIsFullscreen(false); onClose(); }} title="Close map">
                <X size={14} />
              </button>
            </div>
          </div>

          {/* ── Tools panel ── */}
          {openGroup === "tools" && (
            <div className="smap-group-panel">
              <div className="smap-controls">
                <button
                  className={`smap-tbtn smap-draw-btn ${drawMode === "stormMotion" ? "active" : ""}`}
                  onClick={() => toggleDrawMode("stormMotion")}
                  title="Draw storm-motion vector: click start then end point (1-hour displacement)"
                >
                  <Navigation size={11} />
                  Draw Storm Motion
                </button>
                <button
                  className={`smap-tbtn smap-draw-btn smap-draw-boundary ${drawMode === "boundary" ? "active" : ""}`}
                  onClick={() => toggleDrawMode("boundary")}
                  title="Draw boundary orientation line: click two points along the boundary"
                >
                  <Minus size={11} />
                  Draw Boundary
                </button>
                {(drawMode || drawPoints.length > 0) && (
                  <span className="smap-tag smap-tag-draw">
                    {drawPoints.length === 0
                      ? "Click start point"
                      : drawPreview
                        ? drawMode === "stormMotion"
                          ? `${drawPreview.direction}°/${drawPreview.speed} kt — click to set`
                          : `${drawPreview.orientation}° — click to set`
                        : "Move cursor…"}
                  </span>
                )}
                {(stormMotionGuide || boundaryGuide) && (
                  <button
                    className="smap-tbtn smap-draw-clear"
                    onClick={clearDrawGuides}
                    title="Clear drawn storm-motion and boundary guides"
                  >
                    <X size={10} />
                    Clear Guides
                  </button>
                )}
                {stormMotionGuide && (
                  <span className="smap-tag smap-tag-sm">
                    SM {stormMotionGuide.direction}°/{stormMotionGuide.speed} kt
                  </span>
                )}
                {boundaryGuide && (
                  <span className="smap-tag smap-tag-bdry">
                    BDRY {boundaryGuide.orientation}°
                  </span>
                )}
              </div>
            </div>
          )}

          {/* ── Forecasts panel ── */}
          {openGroup === "forecasts" && (
            <>
            <div className="smap-group-panel">
              <div className="smap-controls">
                <button
                  className={`smap-tbtn ${showOutlook ? "active" : ""}`}
                  onClick={() => setShowOutlook((v) => !v)}
                  title="Toggle SPC outlook overlay"
                >
                  <CloudLightning size={11} />
                  SPC Outlook
                </button>
                {showOutlook && (
                  <button
                    className="smap-tbtn smap-refresh-btn"
                    onClick={() => setOutlookRefresh((r) => r + 1)}
                    title="Refresh SPC outlook data"
                  >
                    <RefreshCw size={10} />
                  </button>
                )}
                {showOutlook && (
                  <button
                    className={`smap-tbtn ${showOutlookStations ? "active" : ""}`}
                    onClick={() => setShowOutlookStations((v) => !v)}
                    title="Highlight stations within outlook area — click to fetch forecast soundings"
                  >
                    <Zap size={11} />
                    Outlook Sndgs
                    {outlookStationsLoading && <span className="smap-outlook-loading-dot">...</span>}
                    {outlookStations && <span className="smap-radar-badge">{outlookStations.count}</span>}
                  </button>
                )}
                {showOutlook && (
                  <>
                    <div className="smap-day-btns">
                      {[1, 2, 3].map((d) => (
                        <button
                          key={d}
                          className={`smap-day-btn ${outlookDay === d ? "active" : ""}`}
                          onClick={() => setOutlookDay(d)}
                        >
                          D{d}
                        </button>
                      ))}
                      <span className="smap-day-sep">|</span>
                      {[4, 5, 6, 7, 8].map((d) => (
                        <button
                          key={d}
                          className={`smap-day-btn smap-day-ext ${outlookDay === d ? "active" : ""}`}
                          onClick={() => setOutlookDay(d)}
                          title={`Day ${d} extended-range outlook (probability only)`}
                        >
                          D{d}
                        </button>
                      ))}
                    </div>
                    {outlookDay <= 3 && (
                    <div className="smap-day-btns smap-type-btns">
                      {OUTLOOK_TYPES.map((t) => {
                        const isDisabledD3 = outlookDay === 3 && ["torn", "wind", "hail"].includes(t.id);
                        return (
                          <button
                            key={t.id}
                            className={`smap-day-btn ${effectiveType === t.id ? "active" : ""}${isDisabledD3 ? " smap-btn-disabled" : ""}`}
                            onClick={() => { if (!isDisabledD3) setOutlookType(t.id); }}
                            title={isDisabledD3 ? "Day 3: only categorical available" : (t.title || t.name)}
                          >
                            {t.name}
                          </button>
                        );
                      })}
                    </div>
                    )}
                    {outlookDay >= 4 && (
                      <span className="smap-ext-label">Extended Range – Probability Only</span>
                    )}
                  </>
                )}
              </div>
              {showOutlook && (
                <button
                  className={`smap-tbtn ${showDiscussion ? "active" : ""}`}
                  onClick={() => setShowDiscussion((v) => !v)}
                  title="Show/hide forecaster discussion narrative"
                >
                  <FileText size={11} />
                  Discussion
                </button>
              )}
              {showOutlook && outlookLoading && (
                <span className="smap-outlook-loading">Loading…</span>
              )}
              {showOutlook && outlookMeta && !outlookLoading && (
                <span className="smap-outlook-meta">
                  {outlookMeta.forecaster} · {outlookMeta.valid}
                </span>
              )}
              {showOutlook && outlookError && !outlookLoading && (
                <span className="smap-outlook-error">{outlookError}</span>
              )}
            </div>
            {showOutlook && showDiscussion && (
              <div className="smap-discussion-panel">
                <div className="smap-discussion-header">
                  <FileText size={12} />
                  <span>Day {outlookDay} Convective Outlook Discussion</span>
                  <button className="smap-discussion-close" onClick={() => setShowDiscussion(false)}><X size={12} /></button>
                </div>
                <div className="smap-discussion-body">
                  {discussionLoading ? "Loading discussion…" : discussionText || "No discussion available."}
                </div>
              </div>
            )}
          </>
          )}

          {/* ── Radar panel ── */}
          {openGroup === "radar" && (
            <div className="smap-group-panel">
              <div className="smap-controls">
                <button
                  className={`smap-tbtn ${showRadar ? "active" : ""}`}
                  onClick={() => { setShowRadar((v) => !v); setShowVelocity(false); }}
                  title="Toggle radar reflectivity overlay"
                >
                  <Crosshair size={11} />
                  Radar
                </button>
                {showRadar && (
                  <div className="smap-day-btns smap-radar-source-btns">
                    <button
                      className={`smap-day-btn${radarSource === "composite" ? " active" : ""}`}
                      onClick={() => setRadarSource("composite")}
                      title="RainViewer global composite — animated multi-frame"
                    >
                      Composite
                    </button>
                    <button
                      className={`smap-day-btn${radarSource === "mosaic" ? " active" : ""}`}
                      onClick={() => { setRadarSource("mosaic"); }}
                      title="IEM US national composite — CONUS-wide animated mosaic"
                    >
                      US Mosaic
                    </button>
                    <button
                      className={`smap-day-btn${radarSource === "singlesite" ? " active" : ""}`}
                      onClick={() => setRadarSource("singlesite")}
                      title="Single-site NEXRAD super-res base reflectivity — gate-level detail for hook echoes and storm structure"
                    >
                      Single Site{radarSource === "singlesite" && <span className="smap-radar-badge">{velocityRadar.id}</span>}
                    </button>
                  </div>
                )}
                <span className="smap-group-sep" />
                <button
                  className={`smap-tbtn ${showVelocity ? "active" : ""}`}
                  onClick={() => { setShowVelocity((v) => { if (!v) setVelCacheBust(Date.now()); return !v; }); setShowRadar(false); }}
                  title="Toggle NEXRAD velocity overlay (single-site) — super-res base velocity for tornado rotation signatures"
                >
                  <Wind size={11} />
                  Velocity{showVelocity && <span className="smap-radar-badge">{velocityRadar.id}</span>}
                </button>
                {showVelocity && velAvailProducts.length > 0 && (
                  <div className="smap-day-btns smap-vel-product-btns">
                    {velAvailProducts.map((vp) => (
                      <button
                        key={vp.id}
                        className={`smap-day-btn${velProduct === vp.id ? " active" : ""}${vp.id.endsWith("S") ? " smap-vel-srm" : ""}`}
                        onClick={() => { setVelProduct(vp.id); setVelCacheBust(Date.now()); }}
                        title={vp.desc}
                      >
                        {vp.label}
                      </button>
                    ))}
                  </div>
                )}
                {(showRadar || showVelocity) && showWind && (
                  <button
                    className={`smap-tbtn ${syncTimeline ? "active" : ""}`}
                    onClick={() => setSyncTimeline((v) => !v)}
                    title="Sync wind-flow forecast hour to radar/velocity timeline position"
                  >
                    <Navigation size={11} />
                    Sync Wind
                  </button>
                )}
                {activeTimeline && (
                  <>
                    <button
                      className={`smap-tbtn smap-anim-btn ${activeTimeline.playing ? "active" : ""}`}
                      onClick={toggleTimelinePlayback}
                      title={activeTimeline.playing
                        ? `Pause ${activeTimeline.kind} animation`
                        : `Animate ${activeTimeline.kind} frames`}
                    >
                      {activeTimeline.playing ? <Pause size={10} /> : <Play size={10} />}
                      {activeTimeline.label}
                    </button>
                    <input
                      type="range"
                      className="smap-radar-slider"
                      min={0}
                      max={activeTimeline.count - 1}
                      step={1}
                      value={activeTimeline.frame}
                      onChange={(e) => scrubTimelineFrame(e.target.value)}
                      title={`Scrub through ${activeTimeline.kind} frames`}
                    />
                  </>
                )}
              </div>
            </div>
          )}

          {/* ── Overlays panel ── */}
          {openGroup === "overlays" && (
            <div className="smap-group-panel">
              <div className="smap-controls">
                <button
                  className={`smap-tbtn ${showStorms ? "active" : ""}`}
                  onClick={() => setShowStorms((v) => !v)}
                  title="Toggle TVS (Tornado Vortex Signature) and Mesocyclone detections from NEXRAD storm attributes"
                >
                  <Tornado size={11} />
                  TVS/Meso
                  {stormsLoading && <span className="smap-outlook-loading-dot">...</span>}
                  {tvsCount > 0 && <span className="smap-radar-badge smap-warn-count">{tvsCount} TVS</span>}
                  {mesoCount > 0 && <span className="smap-radar-badge">{mesoCount}</span>}
                </button>
                <button
                  className={`smap-tbtn ${showWarnings ? "active" : ""}`}
                  onClick={() => setShowWarnings((v) => !v)}
                  title="Toggle NWS active warnings — tornado, severe, flash flood, watches"
                >
                  <AlertTriangle size={11} />
                  Warnings
                  {warningsLoading && <span className="smap-outlook-loading-dot">...</span>}
                  {warningCount > 0 && <span className="smap-radar-badge smap-warn-count">{warningCount}</span>}
                </button>
                <button
                  className={`smap-tbtn ${showMdWatch ? "active" : ""}`}
                  onClick={() => setShowMdWatch((v) => !v)}
                  title="Toggle SPC Mesoscale Discussions and active Tornado/Severe Thunderstorm Watches"
                >
                  <FileWarning size={11} />
                  MDs/Watch
                  {mdWatchLoading && <span className="smap-outlook-loading-dot">...</span>}
                  {watchCount > 0 && <span className="smap-radar-badge smap-warn-count">{watchCount}</span>}
                  {mdCount > 0 && <span className="smap-radar-badge">{mdCount} MD</span>}
                </button>
                <button
                  className={`smap-tbtn ${showSpotters ? "active" : ""}`}
                  onClick={() => setShowSpotters((v) => !v)}
                  title="Toggle Spotter Network — live storm chaser positions and field reports (tornado, hail, funnel cloud)"
                >
                  <Binoculars size={11} />
                  Spotters
                  {spottersLoading && <span className="smap-outlook-loading-dot">...</span>}
                  {spotterReportCount > 0 && <span className="smap-radar-badge smap-warn-count">{spotterReportCount}</span>}
                  {spotterCount > 0 && <span className="smap-radar-badge">{spotterCount}</span>}
                </button>
                <button
                  className={`smap-tbtn ${showLightning ? "active" : ""}`}
                  onClick={() => setShowLightning((v) => !v)}
                  title="Toggle real-time lightning strikes (Blitzortung network)"
                >
                  <Zap size={11} />
                  Lightning
                  {lightningStrikes.length > 0 && <span className="smap-radar-badge">{lightningStrikes.length}</span>}
                </button>
                <button
                  className={`smap-tbtn ${showWind ? "active" : ""}`}
                  onClick={() => setShowWind((v) => !v)}
                  title="Toggle animated wind flow overlay"
                >
                  <Wind size={11} />
                  Wind Flow
                </button>
                {showWind && (
                  <>
                    <div className="smap-day-btns smap-wind-level-btns">
                      <button
                        className={`smap-day-btn${windLevel === "500" ? " active" : ""}`}
                        onClick={() => setWindLevel("500")}
                        title="500 hPa steering flow — matches radar echo motion"
                      >
                        Steering
                      </button>
                      <button
                        className={`smap-day-btn${windLevel === "surface" ? " active" : ""}`}
                        onClick={() => setWindLevel("surface")}
                        title="10 m surface wind"
                      >
                        Surface
                      </button>
                    </div>
                    <div className="smap-wind-time">
                      <input
                        type="range"
                        className="smap-radar-slider"
                        min={0}
                        max={WIND_MAX_HOUR_OFFSET}
                        step={1}
                        value={windHourOffset}
                        disabled={syncTimeline && Boolean(activeTimeline)}
                        onChange={(e) => setWindHourOffset(Number(e.target.value))}
                        title={syncTimeline && activeTimeline
                          ? "Wind timeline synced to radar/velocity scrubber"
                          : "Set wind forecast hour offset"}
                      />
                      <span className="smap-radar-badge">{`+${windHourOffset}h`}</span>
                      {syncTimeline && activeTimeline && (
                        <span className="smap-radar-badge">Synced</span>
                      )}
                    </div>
                  </>
                )}
                <span className="smap-group-sep" />
                <button
                  className={`smap-tbtn ${showWindBarbs ? "active" : ""}`}
                  onClick={() => setShowWindBarbs((v) => !v)}
                  title={riskData ? "Show 0-6 km bulk shear vectors from risk-scan data on station markers" : "Run a Risk Scan first to enable shear vectors"}
                  disabled={!riskData}
                >
                  <Navigation size={11} />
                  Shear Vectors
                </button>
                <button
                  className={`smap-tbtn ${showFavorites ? "active" : ""}`}
                  onClick={() => setShowFavorites((v) => !v)}
                  title="Highlight favorite stations with star markers"
                >
                  <Star size={11} />
                  Favorites
                  {favorites.length > 0 && <span className="smap-radar-badge">{favorites.length}</span>}
                </button>
                <button
                  className={`smap-tbtn ${showAcars ? "active" : ""}`}
                  onClick={() => setShowAcars((v) => !v)}
                  title="Show ACARS/AMDAR aircraft observation airports — click to fetch profile"
                >
                  <Plane size={11} />
                  ACARS
                  {acarsLoading && <span className="smap-outlook-loading-dot">...</span>}
                  {acarsAirports && <span className="smap-radar-badge">{acarsAirports.length}</span>}
                </button>
              </div>
            </div>
          )}
        </div>
      </div>

      <div className="smap-container">
        {/* Base map switcher — below zoom controls */}
        <div className="smap-basemap-wrapper">
          <button
            className="smap-basemap-btn"
            onClick={() => setShowBaseMapPicker((v) => !v)}
            title="Change base map"
          >
            <Layers size={14} />
          </button>
          {showBaseMapPicker && (
            <div className="smap-basemap-picker">
              {BASE_MAPS.map((bm) => (
                <button
                  key={bm.id}
                  className={`smap-basemap-opt ${baseMap === bm.id ? "active" : ""}`}
                  onClick={() => { setBaseMap(bm.id); setShowBaseMapPicker(false); }}
                >
                  {bm.name}
                </button>
              ))}
            </div>
          )}
        </div>

        <MapContainer
          center={mapCenter}
          zoom={mapZoom}
          className="smap-leaflet"
          zoomControl={true}
          attributionControl={false}
        >
          <TileLayer
            key={baseMap}
            url={BASE_MAPS.find((b) => b.id === baseMap)?.url}
            maxZoom={18}
          />

          {/* Draw-tool guides: storm-motion vector and boundary line */}
          <Pane name="draw-guides" style={{ zIndex: 465, pointerEvents: "none" }}>
            {drawPoints.length === 1 && (
              <>
                <CircleMarker
                  center={drawPoints[0]}
                  radius={6}
                  pathOptions={{
                    color: drawMode === "boundary" ? "#ff44ff" : "#22c55e",
                    fillColor: drawMode === "boundary" ? "#ff44ff" : "#22c55e",
                    fillOpacity: 0.9,
                    weight: 2,
                  }}
                />
                {drawPreview && (
                  <Polyline
                    positions={[drawPoints[0], drawPreview.cursorPos]}
                    pathOptions={{
                      color: drawMode === "boundary" ? "#ff44ff" : "#22c55e",
                      weight: 2,
                      opacity: 0.5,
                      dashArray: "6 4",
                    }}
                  />
                )}
              </>
            )}
            {stormMotionGuide && (
              <>
                <Polyline
                  positions={[stormMotionGuide.start, stormMotionGuide.end]}
                  pathOptions={{
                    color: "#22c55e",
                    weight: 3,
                    opacity: 0.95,
                  }}
                />
                <CircleMarker
                  center={stormMotionGuide.end}
                  radius={5}
                  pathOptions={{
                    color: "#22c55e",
                    fillColor: "#22c55e",
                    fillOpacity: 0.95,
                    weight: 2,
                  }}
                />
              </>
            )}
            {boundaryGuide && (
              <Polyline
                positions={[boundaryGuide.start, boundaryGuide.end]}
                pathOptions={{
                  color: "#ff44ff",
                  weight: 3,
                  opacity: 0.9,
                  dashArray: "8 6",
                }}
              />
            )}
          </Pane>

          {/* Animated wind particle field — lowest overlay (z-index 250, below tiles/SPC) */}
          {showWind && windData && (
            <WindCanvas
              key={`wind-${windData.level || "500"}-${windData.hourOffset ?? 0}-${Math.round(windData.validAt || 0)}`}
              data={windData}
            />
          )}

          {/* Radar reflectivity — Pane z-index 300 renders ABOVE SPC outlook (250) */}
          <Pane name="radar-tiles" style={{ zIndex: 300 }}>
            {/* Composite (RainViewer) — single persistent layer
                whose URL is swapped via setUrl() to avoid tile-request floods. */}
            {showRadar && radarSource === "composite" && radarFrames[radarFrame]?.path && (
              <RadarLayer
                pane="radar-tiles"
                url={`${RAINVIEWER_TILE}${radarFrames[radarFrame].path}${RAINVIEWER_OPTS}`}
                opacity={0.85}
                maxNativeZoom={7}
                maxZoom={18}
                attribution="Radar via RainViewer"
              />
            )}
            {showRadar && radarSource === "mosaic" && radarFrames[radarFrame] && (
              <RadarLayer
                pane="radar-tiles"
                url={`https://mesonet.agron.iastate.edu/cache/tile.py/1.0.0/ridge::USCOMP-N0Q-${radarFrames[radarFrame].ts}/{z}/{x}/{y}.png`}
                opacity={0.85}
                maxZoom={18}
                attribution="US National Composite via IEM"
              />
            )}
          </Pane>

          {/* Single-site NEXRAD super-res base reflectivity via NWS WMS */}
          {showRadar && radarSource === "singlesite" && (
            <WMSTileLayer
              key={`ss-bref-${velocityRadar.id}`}
              pane="radar-tiles"
              url={`https://opengeo.ncep.noaa.gov/geoserver/k${velocityRadar.id.toLowerCase()}/k${velocityRadar.id.toLowerCase()}_sr_bref/ows`}
              layers={`k${velocityRadar.id.toLowerCase()}_sr_bref`}
              format="image/png"
              transparent={true}
              opacity={0.85}
              maxZoom={18}
              attribution="NWS NEXRAD"
            />
          )}

          {/* Range ring + site marker for single-site radar */}
          {showRadar && radarSource === "singlesite" && (
            <Pane name="ss-range" style={{ zIndex: 310, pointerEvents: "none" }}>
              <Circle
                center={[velocityRadar.lat, velocityRadar.lon]}
                radius={230000}
                pathOptions={{ color: "#00e5ff", weight: 1.2, fillColor: "transparent", fillOpacity: 0, dashArray: "6 4", opacity: 0.5 }}
              />
              <CircleMarker
                center={[velocityRadar.lat, velocityRadar.lon]}
                radius={4}
                pathOptions={{ color: "#00e5ff", fillColor: "#00e5ff", fillOpacity: 1, weight: 1 }}
              >
                <Tooltip direction="top" offset={[0, -6]} permanent className="smap-vel-site-label">
                  {velocityRadar.id}
                </Tooltip>
              </CircleMarker>
            </Pane>
          )}

          {/* NEXRAD Velocity overlay — NWS WMS super-res base velocity */}
          {showVelocity && (
            <WMSTileLayer
              key={`vel-bvel-${velocityRadar.id}`}
              url={`https://opengeo.ncep.noaa.gov/geoserver/k${velocityRadar.id.toLowerCase()}/k${velocityRadar.id.toLowerCase()}_sr_bvel/ows`}
              layers={`k${velocityRadar.id.toLowerCase()}_sr_bvel`}
              format="image/png"
              transparent={true}
              opacity={0.85}
              zIndex={300}
              maxZoom={18}
              attribution="NWS NEXRAD"
            />
          )}

          {/* Range ring + site marker for velocity radar */}
          {showVelocity && (
            <Pane name="vel-range" style={{ zIndex: 310, pointerEvents: "none" }}>
              <Circle
                center={[velocityRadar.lat, velocityRadar.lon]}
                radius={230000}
                pathOptions={{ color: "#00e5ff", weight: 1.2, fillColor: "transparent", fillOpacity: 0, dashArray: "6 4", opacity: 0.5 }}
              />
              <CircleMarker
                center={[velocityRadar.lat, velocityRadar.lon]}
                radius={4}
                pathOptions={{ color: "#00e5ff", fillColor: "#00e5ff", fillOpacity: 1, weight: 1 }}
              >
                <Tooltip direction="top" offset={[0, -6]} permanent className="smap-vel-site-label">
                  {velocityRadar.id}
                </Tooltip>
              </CircleMarker>
            </Pane>
          )}

          {/* SPC outlook polygons (render first so stations sit on top) */}
          {showOutlook && outlookData && <OutlookLayer data={outlookData} outlookType={effectiveType} />}

          {/* NWS active warnings (tornado, severe, flash flood, watches) */}
          {showWarnings && warningsData && <WarningsLayer data={warningsData} />}

          {/* SPC Mesoscale Discussions & Watches */}
          {showMdWatch && <MdWatchLayer mdData={mdData} watchData={watchData} />}

          {/* Lightning strikes */}
          {showLightning && lightningStrikes.length > 0 && (
            <Pane name="lightning-layer" style={{ zIndex: 425 }}>
              {lightningStrikes.map((s) => {
                const age = Date.now() - s.time;
                const opacity = Math.max(0.2, 1 - age / 300_000);
                return (
                  <CircleMarker
                    key={s.id}
                    center={[s.lat, s.lon]}
                    radius={3}
                    pathOptions={{
                      color: "#fff",
                      fillColor: age < 10_000 ? "#ffffff" : age < 60_000 ? "#ffe066" : "#ffaa00",
                      fillOpacity: opacity,
                      weight: age < 10_000 ? 2 : 1,
                      className: age < 10_000 ? "smap-lightning-flash" : "",
                    }}
                  />
                );
              })}
            </Pane>
          )}

          {/* TVS / Mesocyclone storm attribute markers */}
          {showStorms && stormData && stormData.features.length > 0 && (
            <Pane name="storm-attr" style={{ zIndex: 430 }}>
              {stormData.features.map((f, i) => {
                const p = f.properties;
                const [lng, lat] = f.geometry.coordinates;
                const isTvs = p.tvs === "TVS";
                const mesoRank = parseInt(p.meso, 10) || 0;
                const color = isTvs ? "#dc2626" : mesoRank >= 5 ? "#ff4400" : mesoRank >= 3 ? "#ff8800" : "#ffcc00";
                const icon = isTvs ? TVS_ICON : makeMesoIcon(color, mesoRank);
                return (
                  <Marker
                    key={`storm-${p.nexrad}-${p.storm_id}-${i}`}
                    position={[lat, lng]}
                    icon={icon}
                  >
                    <Tooltip direction="top" offset={[0, -12]} className="smap-tooltip">
                      <div className="smap-tt-inner">
                        <div className="smap-tt-header">
                          <span className="smap-tt-id" style={{ color }}>{isTvs ? "\u26a0 TVS" : `Meso ${mesoRank}`}</span>
                          <span className="smap-tt-name">{p.nexrad} {p.storm_id}</span>
                        </div>
                        <span className="smap-tt-type">VIL {p.vil} · {p.max_dbz} dBZ · Top {p.top} kft</span>
                      </div>
                    </Tooltip>
                  </Marker>
                );
              })}
            </Pane>
          )}

          {/* Road / city labels on top of all weather overlays but below popups */}
          <Pane name="top-labels" style={{ zIndex: 399, pointerEvents: "none" }}>
            <TileLayer
              url={BASE_MAPS.find((b) => b.id === baseMap)?.labels}
              maxZoom={18}
              opacity={0.9}
              className="smap-label-tiles"
            />
          </Pane>

          {/* Spotter Network — live chaser positions + field reports */}
          {showSpotters && spotterData && (
            <Pane name="spotter-layer" style={{ zIndex: 435 }}>
              {/* Chaser positions — small green dots */}
              {spotterData.positions.map((s, i) => (
                <CircleMarker
                  key={`sp-${i}`}
                  center={[s.lat, s.lon]}
                  radius={4}
                  pathOptions={{ color: "#00cc44", fillColor: "#00cc44", fillOpacity: 0.7, weight: 1 }}
                >
                  <Tooltip direction="top" offset={[0, -6]} className="smap-tooltip">
                    <div className="smap-tt-inner">
                      <div className="smap-tt-header">
                        <span className="smap-tt-id" style={{ color: "#00cc44" }}>&#x1f441; Spotter</span>
                        <span className="smap-tt-name">{s.name}</span>
                      </div>
                      <span className="smap-tt-type">{s.time}{s.heading ? ` · ${s.heading}` : ""}</span>
                    </div>
                  </Tooltip>
                </CircleMarker>
              ))}
              {/* Field reports — larger color-coded markers */}
              {spotterData.reports.map((r, i) => {
                const rColor = SN_REPORT_COLORS[r.type] || SN_REPORT_COLORS.Other;
                const isTornado = r.type === "Tornado";
                const age = r.sheet === 3 ? "Recent" : r.sheet === 4 ? "Older" : "Old";
                return (
                  <CircleMarker
                    key={`sr-${i}`}
                    center={[r.lat, r.lon]}
                    radius={isTornado ? 11 : 8}
                    pathOptions={{
                      color: rColor,
                      fillColor: rColor,
                      fillOpacity: isTornado ? 0.9 : 0.7,
                      weight: isTornado ? 3 : 2,
                      className: isTornado ? "smap-tvs-pulse" : "",
                    }}
                  >
                    <Tooltip direction="top" offset={[0, -8]} className="smap-tooltip">
                      <div className="smap-tt-inner">
                        <div className="smap-tt-header">
                          <span className="smap-tt-id" style={{ color: rColor }}>{isTornado ? "\u26a0" : "\u25cf"} {r.type}</span>
                          <span className="smap-tt-name">{r.reporter}</span>
                        </div>
                        <span className="smap-tt-type">{r.time} · {age}</span>
                        {r.notes && <span className="smap-tt-type">{r.notes}</span>}
                      </div>
                    </Tooltip>
                  </CircleMarker>
                );
              })}
            </Pane>
          )}

          <MapClickHandler
            enabled={latLonMode}
            onLatLonSelect={onLatLonSelect}
            onRightClick={handleProximitySearch}
            onMapClick={handleMapToolClick}
            onMouseMove={handleDrawMouseMove}
            disableRightClick={Boolean(drawMode)}
          />
          <MapViewTracker onCenterChange={handleCenterChange} onViewChange={setMapView} />
          <FlyToStation station={selectedStation} stations={stations} />
          <MapResizeHandler trigger={isFullscreen} />

          {/* ACARS airports overlay */}
          {showAcars && acarsAirports && (
            <Pane name="acars-layer" style={{ zIndex: 440 }}>
              {acarsAirports.map((ap) => (
                <CircleMarker
                  key={ap.id}
                  center={[ap.lat, ap.lon]}
                  radius={6}
                  pathOptions={{
                    color: "#a855f7",
                    fillColor: "#a855f7",
                    fillOpacity: 0.7,
                    weight: 2,
                  }}
                  eventHandlers={{
                    click: () => {
                      // Select best nearby sounding station or use ACARS source
                      onStationSelect(ap.id.replace(/^K/, ""));
                    },
                  }}
                >
                  <Tooltip direction="top" offset={[0, -8]} className="smap-tooltip">
                    <div className="smap-tt-inner">
                      <div className="smap-tt-header">
                        <span className="smap-tt-id" style={{ color: "#a855f7" }}>✈ {ap.id}</span>
                        <span className="smap-tt-name">{ap.name}</span>
                      </div>
                      <span className="smap-tt-type">ACARS aircraft profiles</span>
                    </div>
                  </Tooltip>
                  <Popup className="smap-popup" maxWidth={260} minWidth={200}>
                    <div className="smap-popup-inner">
                      <div className="smap-popup-header">
                        <span className="smap-popup-id" style={{ color: "#a855f7" }}>✈ {ap.id}</span>
                        <span className="smap-popup-name">{ap.name}</span>
                      </div>
                      <div className="smap-popup-meta">
                        <span className="smap-popup-coords">{ap.lat.toFixed(2)}°, {ap.lon.toFixed(2)}°</span>
                      </div>
                      <div className="smap-popup-outlook">
                        <span className="smap-outlook-badge" style={{ background: "#a855f722", color: "#a855f7", borderColor: "#a855f755" }}>ACARS Airport</span>
                        <div className="smap-forecast-btns">
                          <button className="smap-forecast-btn smap-forecast-btn-auto" onClick={(e) => { e.stopPropagation(); onFetchLatest?.(ap.id.replace(/^K/, "")); }}>Fetch ACARS Profile</button>
                        </div>
                      </div>
                    </div>
                  </Popup>
                </CircleMarker>
              ))}
            </Pane>
          )}

          {/* Proximity search result markers + lines */}
          {proximityResult && (
            <Pane name="proximity-layer" style={{ zIndex: 445 }}>
              <CircleMarker
                center={[proximityResult.lat, proximityResult.lon]}
                radius={6}
                pathOptions={{ color: "#f43f5e", fillColor: "#f43f5e", fillOpacity: 0.9, weight: 2 }}
              />
            </Pane>
          )}

          {stationRender.clusters.length > 0 && (
            <Pane name="station-clusters" style={{ zIndex: 448 }}>
              {stationRender.clusters.map((c) => (
                <CircleMarker
                  key={`cluster-${c.id}`}
                  center={[c.lat, c.lon]}
                  radius={7}
                  pathOptions={{
                    color: c.color || "#60a5fa",
                    fillColor: c.color || "#60a5fa",
                    fillOpacity: 0.28,
                    weight: 1.5,
                    className: "smap-cluster-marker",
                  }}
                  eventHandlers={{
                    click: (e) => {
                      const map = e?.target?._map;
                      if (map) {
                        map.flyTo([c.lat, c.lon], Math.min((map.getZoom() || 4) + 2, 10), { duration: 0.45 });
                      }
                    },
                  }}
                >
                  <Tooltip direction="top" offset={[0, -8]} className="smap-tooltip">
                    <div className="smap-tt-inner">
                      <div className="smap-tt-header">
                        <span className="smap-tt-id">Cluster</span>
                        <span className="smap-tt-name">{c.count} station{c.count !== 1 ? "s" : ""}</span>
                      </div>
                      <span className="smap-tt-type">Zoom in to view individual sites</span>
                    </div>
                  </Tooltip>
                </CircleMarker>
              ))}
            </Pane>
          )}

          {stationRender.stations.map((s) => {
            const risk = riskMap[s.id];
            const stp = risk?.stp;
            const color = riskColor(stp);
            const isSelected = s.id === selectedStation;
            const inOutlook = showOutlookStations && outlookStationIds.has(s.id);
            const outlookInfo = inOutlook
              ? outlookStations?.stations?.find((os) => os.id === s.id)
              : null;
            const radius = isSelected ? 8 : inOutlook ? 7 : stp != null ? (stp >= 1 ? 7 : 5) : 4;

            // Orange glow for stations in the SPC outlook area
            const markerColor = isSelected ? "#3b82f6" : inOutlook ? "#f59e0b" : color;

            return (
              <CircleMarker
                key={s.id}
                center={[s.lat, s.lon]}
                radius={radius}
                pathOptions={{
                  color: markerColor,
                  fillColor: markerColor,
                  fillOpacity: isSelected ? 0.9 : inOutlook ? 0.85 : 0.7,
                  weight: isSelected ? 3 : inOutlook ? 2.5 : 1.5,
                  className: isSelected ? "smap-marker-pulse" : inOutlook ? "smap-marker-outlook" : "",
                }}
                eventHandlers={{
                  click: () => onStationSelect(s.id),
                  popupopen: () => setPopupStationId(s.id),
                  popupclose: () => setPopupStationId((prev) => prev === s.id ? null : prev),
                }}
              >
                {popupStationId !== s.id && (
                  <Tooltip
                    direction="top"
                    offset={[0, -8]}
                    className="smap-tooltip"
                  >
                    <div className="smap-tt-inner">
                      <div className="smap-tt-header">
                        <span className="smap-tt-id">{s.id}</span>
                        <span className="smap-tt-name">{s.name}</span>
                      </div>
                      <span className="smap-tt-type">{s.wmo ? "Radiosonde" : "Model-only"}</span>
                    </div>
                  </Tooltip>
                )}
                <Popup className="smap-popup" maxWidth={260} minWidth={200}>
                  <div className="smap-popup-inner">
                    <div className="smap-popup-header">
                      <span className="smap-popup-id">{s.id}</span>
                      <span className="smap-popup-name">{s.name}</span>
                    </div>
                    <div className="smap-popup-meta">
                      <span className="smap-popup-coords">{s.lat.toFixed(2)}°, {s.lon.toFixed(2)}°</span>
                      {s.wmo && <span className="smap-popup-wmo">WMO {s.wmo}</span>}
                    </div>
                    {risk && (
                      <div className="smap-popup-risk">
                        <div className="smap-popup-risk-item">
                          <span className="smap-popup-risk-label">STP</span>
                          <span className="smap-popup-risk-val">{stp.toFixed(2)}</span>
                        </div>
                        {risk.cape != null && (
                          <div className="smap-popup-risk-item">
                            <span className="smap-popup-risk-label">CAPE</span>
                            <span className="smap-popup-risk-val">{Math.round(risk.cape)}</span>
                          </div>
                        )}
                        {risk.srh != null && (
                          <div className="smap-popup-risk-item">
                            <span className="smap-popup-risk-label">SRH</span>
                            <span className="smap-popup-risk-val">{Math.round(risk.srh)}</span>
                          </div>
                        )}
                      </div>
                    )}
                    {inOutlook && (
                      <div className="smap-popup-outlook">
                        <span className="smap-outlook-badge">{outlookInfo?.riskLabel} Outlook</span>
                        <div className="smap-forecast-btns">
                          <button className="smap-forecast-btn smap-forecast-btn-auto" onClick={(e) => { e.stopPropagation(); onFetchLatest?.(s.id); }}>Fetch Latest</button>
                        </div>
                      </div>
                    )}
                  </div>
                </Popup>
              </CircleMarker>
            );
          })}

          {/* Wind barb arrows from risk-scan shear data */}
          {showWindBarbs && riskData && (
            <Pane name="wind-barbs-layer" style={{ zIndex: 455, pointerEvents: "none" }}>
              {riskData.stations.filter(s => s.bwdDir != null && s.bwd > 5).map((s) => {
                const rad = ((s.bwdDir + 180) % 360) * Math.PI / 180; // point in shear direction
                const len = Math.min(0.8, s.bwd / 60); // degrees, scaled by magnitude
                const endLat = s.lat + len * Math.cos(rad);
                const endLon = s.lon + len * Math.sin(rad) / Math.cos(s.lat * Math.PI / 180);
                const color = s.bwd >= 50 ? "#ef4444" : s.bwd >= 35 ? "#f59e0b" : s.bwd >= 20 ? "#22c55e" : "#888";
                // SVG arrow rendered as polyline overlay
                return (
                  <Pane key={`wb-${s.id}`} name={`wb-${s.id}`} style={{ zIndex: 455, pointerEvents: "none" }}>
                    <CircleMarker
                      center={[endLat, endLon]}
                      radius={3}
                      pathOptions={{ color, fillColor: color, fillOpacity: 0.9, weight: 1 }}
                    />
                    {/* Line from station to arrow tip */}
                    <GeoJSON
                      data={{ type: "Feature", geometry: { type: "LineString", coordinates: [[s.lon, s.lat], [endLon, endLat]] }, properties: {} }}
                      style={() => ({ color, weight: 2, opacity: 0.8 })}
                    />
                  </Pane>
                );
              })}
            </Pane>
          )}

          {/* Favorite station star markers */}
          {showFavorites && favorites.length > 0 && (
            <Pane name="favorites-layer" style={{ zIndex: 460 }}>
              {favorites.map((fid) => {
                const s = stations.find((st) => st.id === fid);
                if (!s) return null;
                return (
                  <Marker
                    key={`fav-${fid}`}
                    position={[s.lat, s.lon]}
                    icon={L.divIcon({
                      className: "smap-fav-icon",
                      iconSize: [20, 20],
                      iconAnchor: [10, 10],
                      html: `<svg width="20" height="20" viewBox="0 0 24 24" fill="#f59e0b" stroke="#fff" stroke-width="1.5"><polygon points="12,2 15.09,8.26 22,9.27 17,14.14 18.18,21.02 12,17.77 5.82,21.02 7,14.14 2,9.27 8.91,8.26"/></svg>`,
                    })}
                    eventHandlers={{
                      click: () => onStationSelect(fid),
                    }}
                  >
                    <Tooltip direction="top" offset={[0, -12]} className="smap-tooltip">
                      <div className="smap-tt-inner">
                        <div className="smap-tt-header">
                          <span className="smap-tt-id" style={{ color: "#f59e0b" }}>★ {s.id}</span>
                          <span className="smap-tt-name">{s.name}</span>
                        </div>
                        <span className="smap-tt-type">Favorite station</span>
                      </div>
                    </Tooltip>
                  </Marker>
                );
              })}
            </Pane>
          )}
        </MapContainer>

        {/* Floating legend overlays — each in its own container */}
        <div className="smap-legend-stack">
          {showSpotters && spotterData && (
            <div className="smap-legend-float">
              <div className="smap-legend smap-legend-spotters">
                <span className="smap-legend-item">
                  <span className="smap-dot" style={{ background: "#00cc44" }} />
                  {spotterCount} chaser{spotterCount !== 1 ? "s" : ""}
                </span>
                {spotterReportCount > 0 && (
                  <span className="smap-legend-item">
                    <span className="smap-dot" style={{ background: "#ff0000" }} />
                    {spotterReportCount} report{spotterReportCount !== 1 ? "s" : ""}
                  </span>
                )}
              </div>
            </div>
          )}
          {showRadar && radarSource === "singlesite" && (
            <div className="smap-legend-float smap-legend-vel">
              <div className="smap-legend smap-legend-vel-row">
                <span className="smap-legend-vel-title">
                  {velocityRadar.id} SR Bref
                </span>
              </div>
            </div>
          )}
          {showVelocity && (
            <div className="smap-legend-float smap-legend-vel">
              <div className="smap-legend smap-legend-vel-row">
                <span className="smap-legend-vel-title">
                  {velocityRadar.id} {(() => { const vp = VELOCITY_PRODUCTS.find((p) => p.id === velProduct); return vp ? vp.label : velProduct; })()}
                </span>
                {velFrames[velFrame] ? (
                  <span className="smap-legend-vel-time">{fmtRadarTime(velFrames[velFrame].time)}</span>
                ) : velScanTime ? (
                  <span className="smap-legend-vel-time">{velScanTime}</span>
                ) : null}
              </div>
              <div className="smap-legend smap-legend-vel-bar">
                <span className="smap-vel-neg">Toward</span>
                <span className="smap-vel-gradient" />
                <span className="smap-vel-pos">Away</span>
              </div>
            </div>
          )}
          {riskData && (
            <div className="smap-legend-float">
              <div className="smap-legend">
                <span className="smap-legend-item"><span className="smap-dot" style={{ background: RISK_COLORS.high }} />STP &ge; 1</span>
                <span className="smap-legend-item"><span className="smap-dot" style={{ background: RISK_COLORS.med }} />0.3–1</span>
                <span className="smap-legend-item"><span className="smap-dot" style={{ background: RISK_COLORS.low }} />&lt; 0.3</span>
              </div>
            </div>
          )}
          {stationRender.hiddenCount > 0 && (
            <div className="smap-legend-float">
              <div className="smap-legend smap-legend-decimation">
                <span className="smap-legend-item">
                  <span className="smap-dot" style={{ background: "#60a5fa" }} />
                  {stationRender.hiddenCount} station{stationRender.hiddenCount !== 1 ? "s" : ""} grouped at this zoom
                </span>
                <span className="smap-legend-hint">Click a cluster or zoom in for individual markers</span>
              </div>
            </div>
          )}
          {showOutlookStations && outlookStations && (
            <div className="smap-legend-float">
              <div className="smap-legend smap-legend-outlook-stations">
                <span className="smap-legend-item">
                  <span className="smap-dot" style={{ background: "#f59e0b" }} />
                  {outlookStations.count} station{outlookStations.count !== 1 ? "s" : ""} in outlook
                </span>
                <span className="smap-legend-hint">Click station for forecast sounding</span>
              </div>
            </div>
          )}
          {showOutlook && activeCategories.length > 0 && (
            <div className="smap-legend-float">
              <div className="smap-legend smap-legend-outlook">
                {["torn", "wind", "hail", "prob"].includes(effectiveType) && (
                  <span className="smap-legend-ci-label">
                    {effectiveType === "prob" ? `Day ${outlookDay} Extended` : "Probabilistic"}
                  </span>
                )}
                {effectiveType === "cat" ? (
                  activeCategories.map((c) => (
                    <span key={c.label} className="smap-legend-item smap-outlook-cat" title={c.info}>
                      <span className="smap-dot-rect" style={{ background: c.fill, borderColor: c.color }} />
                      <span className="smap-cat-label">{c.name}</span>
                    </span>
                  ))
                ) : (
                  activeCategories.map((p) => (
                    <span key={p.label} className="smap-legend-item smap-outlook-cat">
                      {p.label === "SIGN" ? (
                        <span className="smap-dot-rect smap-sig-hatch" style={{ borderColor: "#111" }} />
                      ) : (
                        <span className="smap-dot-rect" style={{ background: p.fill, borderColor: p.color }} />
                      )}
                      <span className="smap-cat-label">{p.label === "SIGN" ? "Sig" : p.pct}</span>
                    </span>
                  ))
                )}
              </div>
            </div>
          )}
          {showAcars && acarsAirports && (
            <div className="smap-legend-float">
              <div className="smap-legend">
                <span className="smap-legend-item"><span className="smap-dot" style={{ background: "#a855f7" }} />{acarsAirports.length} ACARS airports</span>
                <span className="smap-legend-hint">Click airport for aircraft profile</span>
              </div>
            </div>
          )}
          {showWindBarbs && riskData && (
            <div className="smap-legend-float">
              <div className="smap-legend">
                <span className="smap-legend-item"><span className="smap-dot" style={{ background: "#ef4444" }} />BWD ≥ 50 kt</span>
                <span className="smap-legend-item"><span className="smap-dot" style={{ background: "#f59e0b" }} />35–50 kt</span>
                <span className="smap-legend-item"><span className="smap-dot" style={{ background: "#22c55e" }} />20–35 kt</span>
                <span className="smap-legend-hint">0-6 km bulk shear vectors</span>
              </div>
            </div>
          )}
        </div>
        {/* Proximity panel — rendered OUTSIDE Leaflet so labels can't overlap */}
        {proximityResult && (
          <div className="smap-prox-overlay">
            <div className="smap-prox-panel">
              <div className="smap-prox-header">
                <Crosshair size={14} />
                <span>Nearest Stations</span>
                <button className="smap-prox-close" onClick={() => setProximityResult(null)}><X size={14} /></button>
              </div>
              <div className="smap-prox-coords">
                {proximityResult.lat.toFixed(2)}°, {proximityResult.lon.toFixed(2)}°
              </div>
              <div className="smap-prox-list">
                {proximityResult.stations.map((ps, i) => (
                  <button key={ps.id} className="smap-prox-item" onClick={() => onStationSelect(ps.id)}>
                    <span className="smap-prox-rank">{i + 1}</span>
                    <span className="smap-prox-sid">{ps.id}</span>
                    <span className="smap-prox-name">{ps.name}</span>
                    <span className="smap-prox-dist">{ps.dist.toFixed(0)} km</span>
                  </button>
                ))}
              </div>
              <div className="smap-prox-actions">
                <select className="smap-prox-select" id="prox-compare-count" defaultValue="4">
                  <option value="2">Top 2</option>
                  <option value="3">Top 3</option>
                  <option value="4">Top 4</option>
                </select>
                <button className="smap-prox-compare" onClick={() => {
                  const n = parseInt(document.getElementById('prox-compare-count').value, 10);
                  onCompareStations?.(proximityResult.stations.slice(0, n).map(s => s.id));
                }}>Compare</button>
                <button className="smap-prox-dismiss" onClick={() => setProximityResult(null)}>Dismiss</button>
              </div>
            </div>
          </div>
        )}
        <div className="smap-hint-bar">
          {drawMode
            ? "Draw mode active: click two points on the map"
            : "Right-click map for nearest stations"}
        </div>
      </div>
    </section>
  );
}
