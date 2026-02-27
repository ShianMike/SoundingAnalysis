import { useState, useEffect, useRef, useCallback, useMemo } from "react";
import { MapContainer, TileLayer, CircleMarker, Popup, Pane, GeoJSON, useMapEvents, useMap } from "react-leaflet";
import { X, Crosshair, CloudLightning, Wind } from "lucide-react";
import { fetchSpcOutlook } from "../api";
import "leaflet/dist/leaflet.css";
import "./StationMap.css";

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
  return best[0];
}

/* ── Track map center for nearest-radar selection ───────────── */
function MapCenterTracker({ onCenterChange }) {
  const map = useMap();
  useMapEvents({
    moveend() {
      const c = map.getCenter();
      onCenterChange(c.lat, c.lng);
    },
  });
  // Fire once on mount
  useEffect(() => {
    const c = map.getCenter();
    onCenterChange(c.lat, c.lng);
  }, []); // eslint-disable-line react-hooks/exhaustive-deps
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

function spcStyle(feature) {
  const label = feature.properties?.LABEL || "";
  const cat = SPC_CATEGORIES.find((c) => c.label === label);
  return {
    color: cat?.color || feature.properties?.stroke || "#888",
    fillColor: cat?.fill || feature.properties?.fill || "#888",
    fillOpacity: 0.25,
    weight: 1.5,
  };
}

/* ── GeoJSON overlay in a non-interactive pane below markers ── */
function OutlookLayer({ data }) {
  if (!data || !data.features || data.features.length === 0) return null;
  return (
    <Pane name="spc-outlook" style={{ zIndex: 250, pointerEvents: "none" }}>
      <GeoJSON
        key={JSON.stringify(data).slice(0, 100)}
        data={data}
        style={spcStyle}
        interactive={false}
      />
    </Pane>
  );
}

/* ── click-for-coords layer ─────────────────────────────────── */
function MapClickHandler({ enabled, onLatLonSelect }) {
  useMapEvents({
    click(e) {
      if (enabled && onLatLonSelect) {
        onLatLonSelect(e.latlng.lat.toFixed(2), e.latlng.lng.toFixed(2));
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

/* ── main component ─────────────────────────────────────────── */
export default function StationMap({
  stations = [],
  riskData,
  selectedStation,
  onStationSelect,
  onLatLonSelect,
  latLonMode,   // true when source is rap
  onClose,
}) {
  const [outlookDay, setOutlookDay] = useState(1);
  const [outlookData, setOutlookData] = useState(null);
  const [outlookLoading, setOutlookLoading] = useState(false);
  const [showOutlook, setShowOutlook] = useState(true);
  const [showRadar, setShowRadar] = useState(false);
  const [showVelocity, setShowVelocity] = useState(false);
  const [velocityRadar, setVelocityRadar] = useState("TLX");   // nearest NEXRAD ID

  const handleCenterChange = useCallback((lat, lng) => {
    setVelocityRadar(nearestNexrad(lat, lng));
  }, []);

  // Fetch SPC outlook on mount and when day changes
  useEffect(() => {
    if (!showOutlook) return;
    let cancelled = false;
    setOutlookLoading(true);
    fetchSpcOutlook(outlookDay)
      .then((data) => { if (!cancelled) setOutlookData(data); })
      .catch(() => { if (!cancelled) setOutlookData(null); })
      .finally(() => { if (!cancelled) setOutlookLoading(false); });
    return () => { cancelled = true; };
  }, [outlookDay, showOutlook]);

  // Determine which SPC categories are present in the current data
  const activeCategories = useMemo(() => {
    if (!outlookData?.features) return [];
    const labels = new Set(outlookData.features.map((f) => f.properties?.LABEL));
    return SPC_CATEGORIES.filter((c) => labels.has(c.label));
  }, [outlookData]);

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

  return (
    <section className="station-map-panel">
      <header className="smap-header">
        <div className="smap-title-group">
          <Crosshair size={14} />
          <span className="smap-title">Station Map</span>
          {latLonMode && (
            <span className="smap-tag">Click map to set coordinates</span>
          )}
        </div>
        <button className="smap-close" onClick={onClose} title="Close map">
          <X size={16} />
        </button>
      </header>

      {/* SPC Outlook & Radar controls */}
      <div className="smap-outlook-bar">
        <button
          className={`smap-outlook-toggle ${showOutlook ? "active" : ""}`}
          onClick={() => setShowOutlook((v) => !v)}
          title="Toggle SPC outlook overlay"
        >
          <CloudLightning size={12} />
          SPC Outlook
        </button>
        <button
          className={`smap-outlook-toggle ${showRadar ? "active" : ""}`}
          onClick={() => { setShowRadar((v) => !v); setShowVelocity(false); }}
          title="Toggle NEXRAD radar reflectivity overlay"
        >
          <Crosshair size={12} />
          Radar
        </button>
        <button
          className={`smap-outlook-toggle ${showVelocity ? "active" : ""}`}
          onClick={() => { setShowVelocity((v) => !v); setShowRadar(false); }}
          title="Toggle NEXRAD storm-relative velocity overlay"
        >
          <Wind size={12} />
          Velocity{showVelocity && <span className="smap-radar-badge">{velocityRadar}</span>}
        </button>
        {showOutlook && (
          <div className="smap-day-btns">
            {[1, 2, 3].map((d) => (
              <button
                key={d}
                className={`smap-day-btn ${outlookDay === d ? "active" : ""}`}
                onClick={() => setOutlookDay(d)}
              >
                Day {d}
              </button>
            ))}
          </div>
        )}
        {showOutlook && outlookLoading && (
          <span className="smap-outlook-loading">Loading…</span>
        )}
        {showOutlook && outlookMeta && !outlookLoading && (
          <span className="smap-outlook-meta">
            {outlookMeta.forecaster} · Valid {outlookMeta.valid}
          </span>
        )}
      </div>

      {/* Legends */}
      <div className="smap-legends">
        {riskData && (
          <div className="smap-legend">
            <span className="smap-legend-item">
              <span className="smap-dot" style={{ background: RISK_COLORS.high }} />
              STP &ge; 1
            </span>
            <span className="smap-legend-item">
              <span className="smap-dot" style={{ background: RISK_COLORS.med }} />
              0.3 &ndash; 1
            </span>
            <span className="smap-legend-item">
              <span className="smap-dot" style={{ background: RISK_COLORS.low }} />
              &lt; 0.3
            </span>
          </div>
        )}
        {showOutlook && activeCategories.length > 0 && (
          <div className="smap-legend smap-legend-outlook">
            {activeCategories.map((c) => (
              <span key={c.label} className="smap-legend-item smap-outlook-cat" title={c.info}>
                <span className="smap-dot-rect" style={{ background: c.fill, borderColor: c.color }} />
                <span className="smap-cat-label">{c.name}</span>
                <span className="smap-cat-info">{c.info}</span>
              </span>
            ))}
          </div>
        )}
      </div>

      <div className="smap-container">
        <MapContainer
          center={mapCenter}
          zoom={mapZoom}
          className="smap-leaflet"
          zoomControl={true}
          attributionControl={false}
        >
          <TileLayer
            url="https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png"
            maxZoom={18}
          />

          {/* NEXRAD Radar reflectivity overlay (IEM) */}
          {showRadar && (
            <TileLayer
              url="https://mesonet.agron.iastate.edu/cache/tile.py/1.0.0/nexrad-n0q-900913/{z}/{x}/{y}.png"
              opacity={0.55}
              maxZoom={12}
              attribution="NEXRAD via IEM"
            />
          )}

          {/* NEXRAD Storm-Relative Velocity overlay (RIDGE single-site) */}
          {showVelocity && (
            <TileLayer
              key={`vel-${velocityRadar}`}
              url={`https://mesonet.agron.iastate.edu/cache/tile.py/1.0.0/ridge::${velocityRadar}-N0U-0/{z}/{x}/{y}.png`}
              opacity={0.55}
              maxZoom={12}
              attribution={`SRV ${velocityRadar} via IEM`}
            />
          )}

          {/* SPC outlook polygons (render first so stations sit on top) */}
          {showOutlook && outlookData && <OutlookLayer data={outlookData} />}

          <MapClickHandler enabled={latLonMode} onLatLonSelect={onLatLonSelect} />
          <MapCenterTracker onCenterChange={handleCenterChange} />
          <FlyToStation station={selectedStation} stations={stations} />

          {stations.map((s) => {
            const risk = riskMap[s.id];
            const stp = risk?.stp;
            const color = riskColor(stp);
            const isSelected = s.id === selectedStation;
            const radius = isSelected ? 8 : stp != null ? (stp >= 1 ? 7 : 5) : 4;

            return (
              <CircleMarker
                key={s.id}
                center={[s.lat, s.lon]}
                radius={radius}
                pathOptions={{
                  color: isSelected ? "#3b82f6" : color,
                  fillColor: isSelected ? "#3b82f6" : color,
                  fillOpacity: isSelected ? 0.9 : 0.7,
                  weight: isSelected ? 3 : 1.5,
                }}
                eventHandlers={{
                  click: () => onStationSelect(s.id),
                }}
              >
                <Popup className="smap-popup">
                  <div className="smap-popup-inner">
                    <strong>{s.id}</strong> &mdash; {s.name}
                    <br />
                    <span className="smap-popup-coords">
                      {s.lat.toFixed(2)}, {s.lon.toFixed(2)}
                    </span>
                    {risk && (
                      <div className="smap-popup-risk">
                        <span>STP: <b>{stp.toFixed(2)}</b></span>
                        {risk.cape != null && <span>CAPE: <b>{Math.round(risk.cape)}</b></span>}
                        {risk.srh != null && <span>SRH: <b>{Math.round(risk.srh)}</b></span>}
                      </div>
                    )}
                  </div>
                </Popup>
              </CircleMarker>
            );
          })}
        </MapContainer>
      </div>
    </section>
  );
}
