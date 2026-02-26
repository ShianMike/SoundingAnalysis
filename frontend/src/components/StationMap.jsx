import { useEffect, useRef, useCallback, useMemo } from "react";
import { MapContainer, TileLayer, CircleMarker, Popup, useMapEvents, useMap } from "react-leaflet";
import { X, Crosshair } from "lucide-react";
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

      {/* Legend */}
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
