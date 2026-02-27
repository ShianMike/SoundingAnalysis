import { useState, useEffect, useRef, useCallback, useMemo } from "react";
import { MapContainer, TileLayer, CircleMarker, Popup, Pane, GeoJSON, useMapEvents, useMap } from "react-leaflet";
import { X, Crosshair, CloudLightning } from "lucide-react";
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
  latLonMode,   // true when source is rap/era5
  onClose,
}) {
  const [outlookDay, setOutlookDay] = useState(1);
  const [outlookData, setOutlookData] = useState(null);
  const [outlookLoading, setOutlookLoading] = useState(false);
  const [showOutlook, setShowOutlook] = useState(true);

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

      {/* SPC Outlook controls */}
      <div className="smap-outlook-bar">
        <button
          className={`smap-outlook-toggle ${showOutlook ? "active" : ""}`}
          onClick={() => setShowOutlook((v) => !v)}
          title="Toggle SPC outlook overlay"
        >
          <CloudLightning size={12} />
          SPC Outlook
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
          center={CONUS_CENTER}
          zoom={CONUS_ZOOM}
          className="smap-leaflet"
          zoomControl={true}
          attributionControl={false}
        >
          <TileLayer
            url="https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png"
            maxZoom={18}
          />

          {/* SPC outlook polygons (render first so stations sit on top) */}
          {showOutlook && outlookData && <OutlookLayer data={outlookData} />}

          <MapClickHandler enabled={latLonMode} onLatLonSelect={onLatLonSelect} />
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
