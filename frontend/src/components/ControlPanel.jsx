import { useState, useRef, useEffect } from "react";
import {
  Search,
  MapPin,
  Calendar,
  Database,
  Layers,
  Clock,
  Loader2,
  ChevronDown,
  Zap,
  ArrowUpDown,
  History,
  Map,
  Star,
  TrendingUp,
  GitCompareArrows,
  MessageSquarePlus,
  Github,
  Thermometer,
  Wind,
} from "lucide-react";
import { fetchRiskScan } from "../api";
import { getFavorites, toggleFavorite } from "../favorites";
import "./ControlPanel.css";

const SOURCE_META = {
  obs: {
    label: "Observed Radiosonde",
    desc: "Real observed upper-air data from the Iowa Environmental Mesonet and University of Wyoming archives.",
  },
  rap: {
    label: "RAP Model Analysis",
    desc: "Rapid Refresh model analysis at any lat/lon point over CONUS via NCEI THREDDS. ~13 km resolution.",
  },
  bufkit: {
    label: "BUFKIT Forecast",
    desc: "Station-based forecast soundings from HRRR, RAP, NAM, GFS, and other models via Iowa State archive.",
  },
  era5: {
    label: "ERA5 Reanalysis",
    desc: "ECMWF ERA5 global reanalysis at any lat/lon. Covers 1940 to near-present. Requires CDS API key.",
  },
  acars: {
    label: "ACARS Aircraft Obs",
    desc: "ACARS/AMDAR aircraft observation profiles at major airports from the IEM archive.",
  },
  igrav2: {
    label: "IGRAv2 Global",
    desc: "IGRA v2 global radiosonde archive via UWyo. Enter a WMO station ID for any station worldwide.",
  },
};

export default function ControlPanel({
  stations,
  sources,
  models,
  onSubmit,
  loading,
  initialLoading,
  onRetry,
  connectError,
  riskData,
  onRiskDataChange,
  showHistory,
  onToggleHistory,
  showMap,
  onToggleMap,
  showTimeSeries,
  onToggleTimeSeries,
  showCompare,
  onToggleCompare,
  selectedStation,
  onStationChange,
  onSourceChange,
  mapLatLon,
  onFeedbackClick,
  showFeedback: feedbackActive,
}) {
  const [source, setSourceLocal] = useState("obs");
  const [station, setStationLocal] = useState("OUN");
  const [date, setDate] = useState("");
  const [lat, setLat] = useState("");
  const [lon, setLon] = useState("");
  const [model, setModel] = useState("hrrr");
  const [fhour, setFhour] = useState("0");
  const [stationSearch, setStationSearch] = useState("");
  const [hoveredSource, setHoveredSource] = useState(null);
  const [scanning, setScanning] = useState(false);
  const [sortMode, setSortMode] = useState("az");
  const [favorites, setFavorites] = useState(() => getFavorites());
  const [soundingHour, setSoundingHour] = useState("12");
  const listRef = useRef(null);

  // Surface modification state
  const [sfcModEnabled, setSfcModEnabled] = useState(false);
  const [sfcModT, setSfcModT] = useState("");
  const [sfcModTd, setSfcModTd] = useState("");
  const [sfcModWspd, setSfcModWspd] = useState("");
  const [sfcModWdir, setSfcModWdir] = useState("");

  // Custom storm motion state
  const [smEnabled, setSmEnabled] = useState(false);
  const [smDirection, setSmDirection] = useState("");
  const [smSpeed, setSmSpeed] = useState("");

  // IGRAv2 WMO ID input
  const [wmoId, setWmoId] = useState("");

  // Sync source to parent
  const setSource = (src) => {
    setSourceLocal(src);
    if (onSourceChange) onSourceChange(src);
  };

  // Sync station to parent
  const setStation = (id) => {
    setStationLocal(id);
    if (onStationChange) onStationChange(id);
  };

  // Sync station from parent (map click)
  useEffect(() => {
    if (selectedStation && selectedStation !== station) {
      setStationLocal(selectedStation);
      const stn = stations.find((s) => s.id === selectedStation);
      if (stn) {
        setLat(String(stn.lat));
        setLon(String(stn.lon));
      }
    }
  }, [selectedStation]); // eslint-disable-line react-hooks/exhaustive-deps

  // Sync lat/lon from map click
  useEffect(() => {
    if (mapLatLon) {
      setLat(String(mapLatLon.lat));
      setLon(String(mapLatLon.lon));
    }
  }, [mapLatLon]);

  // Scroll selected station into view on mount
  useEffect(() => {
    if (listRef.current) {
      const active = listRef.current.querySelector(".cp-station-item.active");
      if (active) active.scrollIntoView({ block: "center" });
    }
  }, [stations]);

  const needsLatLon = source === "rap" || source === "era5";
  const needsStation = source === "obs" || source === "bufkit" || source === "acars" || source === "igrav2";
  const needsModel = source === "bufkit";
  const needsWmoId = source === "igrav2";

  // Build risk lookup from scan results
  const riskMap = {};
  if (riskData) {
    riskData.stations.forEach((s, i) => {
      riskMap[s.id] = { rank: i + 1, stp: s.stp, raw: s.raw, cape: s.cape, srh: s.srh, bwd: s.bwd };
    });
  }

  // Merge station data with risk scores
  const mergedStations = stations.map((s) => ({
    ...s,
    risk: riskMap[s.id] || null,
  }));

  // Sort based on selected mode
  const sortedStations = [...mergedStations].sort((a, b) => {
    switch (sortMode) {
      case "za":
        return b.id.localeCompare(a.id);
      case "favs": {
        const aFav = favorites.includes(a.id) ? 0 : 1;
        const bFav = favorites.includes(b.id) ? 0 : 1;
        if (aFav !== bFav) return aFav - bFav;
        return a.id.localeCompare(b.id);
      }
      case "risk-high":
        // Stations with risk first (highest STP first), then unscanned at end
        if (a.risk && !b.risk) return -1;
        if (!a.risk && b.risk) return 1;
        if (a.risk && b.risk) return b.risk.stp - a.risk.stp;
        return a.id.localeCompare(b.id);
      case "risk-low":
        if (a.risk && !b.risk) return -1;
        if (!a.risk && b.risk) return 1;
        if (a.risk && b.risk) return a.risk.stp - b.risk.stp;
        return a.id.localeCompare(b.id);
      case "az":
      default:
        return a.id.localeCompare(b.id);
    }
  });

  const filteredStations = sortedStations.filter(
    (s) =>
      s.id.toLowerCase().includes(stationSearch.toLowerCase()) ||
      s.name.toLowerCase().includes(stationSearch.toLowerCase())
  );

  const handleRiskScan = async () => {
    setScanning(true);
    try {
      const dateParam = date ? date.replace(/[-T:]/g, "").slice(0, 10) : undefined;
      const data = await fetchRiskScan(dateParam);
      onRiskDataChange(data);
      setSortMode("risk-high");
      // Auto-open the map
      if (!showMap && onToggleMap) onToggleMap();
      // Auto-select the highest risk station
      if (data.stations.length > 0) {
        handleStationSelect(data.stations[0].id);
      }
    } catch (e) {
      console.error("Risk scan failed:", e);
    } finally {
      setScanning(false);
    }
  };

  const handleSubmit = (e) => {
    e.preventDefault();
    const params = { source };

    if (needsStation) params.station = station;
    if (needsWmoId && wmoId) params.station = wmoId;  // WMO override
    if (needsLatLon) {
      params.lat = parseFloat(lat);
      params.lon = parseFloat(lon);
    }
    if (date) params.date = date.replace(/[-T:]/g, "").slice(0, 10);
    if (needsModel) {
      params.model = model;
      params.fhour = parseInt(fhour) || 0;
    }

    // Surface modification
    if (sfcModEnabled) {
      const mod = {};
      if (sfcModT !== "") mod.temperature = parseFloat(sfcModT);
      if (sfcModTd !== "") mod.dewpoint = parseFloat(sfcModTd);
      if (sfcModWspd !== "") mod.wind_speed = parseFloat(sfcModWspd);
      if (sfcModWdir !== "") mod.wind_direction = parseFloat(sfcModWdir);
      if (Object.keys(mod).length > 0) params.surfaceMod = mod;
    }

    // Custom storm motion
    if (smEnabled && smDirection !== "" && smSpeed !== "") {
      params.stormMotion = {
        direction: parseFloat(smDirection),
        speed: parseFloat(smSpeed),
      };
    }

    onSubmit(params);
  };

  const handleStationSelect = (id, autoFetch = false) => {
    setStation(id);
    setStationSearch("");
    const stn = stations.find((s) => s.id === id);
    if (stn) {
      setLat(String(stn.lat));
      setLon(String(stn.lon));
      // Auto-fill WMO ID for IGRAv2
      if (stn.wmo) setWmoId(stn.wmo);
    }
    if (autoFetch && !loading) {
      const params = { source, station: id };
      if (date) params.date = date.replace(/[-T:]/g, "").slice(0, 10);
      if (source === "bufkit") {
        params.model = model;
        params.fhour = parseInt(fhour) || 0;
      }
      if (source === "rap" || source === "era5") {
        const s = stations.find((st) => st.id === id);
        if (s) {
          params.lat = s.lat;
          params.lon = s.lon;
        }
      }
      onSubmit(params);
    }
  };

  if (initialLoading || connectError) {
    return (
      <aside className="control-panel">
        <div className="cp-brand">
          <svg width="28" height="28" viewBox="0 0 28 28" fill="none" xmlns="http://www.w3.org/2000/svg">
            <rect x="2" y="2" width="24" height="24" rx="4" fill="rgba(59,130,246,0.12)" stroke="#3b82f6" strokeWidth="1.5"/>
            <path d="M8 22 L10 18 L11 16 L12 13 L14 10 L16 8 L19 5" stroke="#ef4444" strokeWidth="1.8" strokeLinecap="round" strokeLinejoin="round" fill="none"/>
            <path d="M6 22 L8 19 L9 17 L9.5 15 L10 13 L10 11 L10.5 9 L11 6" stroke="#22c55e" strokeWidth="1.8" strokeLinecap="round" strokeLinejoin="round" fill="none"/>
            <line x1="21" y1="10" x2="21" y2="22" stroke="#60a5fa" strokeWidth="1.2" strokeLinecap="round"/>
            <line x1="21" y1="10" x2="24" y2="8" stroke="#60a5fa" strokeWidth="1.2" strokeLinecap="round"/>
            <line x1="21" y1="13" x2="24" y2="11" stroke="#60a5fa" strokeWidth="1.2" strokeLinecap="round"/>
          </svg>
          <div>
            <h1 className="cp-brand-title">Sounding Analysis</h1>
            <p className="cp-brand-sub">Atmospheric Profile Tool</p>
          </div>
        </div>
        <div className="cp-loading">
          {initialLoading ? (
            <>
              <Loader2 className="spin" size={20} />
              <span>Connecting to API...</span>
            </>
          ) : (
            <>
              <span style={{ color: "var(--danger, #e74c3c)" }}>{connectError}</span>
              <button
                type="button"
                onClick={onRetry}
                style={{
                  marginTop: 8,
                  padding: "6px 16px",
                  cursor: "pointer",
                  borderRadius: 6,
                  border: "1px solid var(--border, #555)",
                  background: "var(--surface, #2a2a2a)",
                  color: "inherit",
                }}
              >
                Retry
              </button>
            </>
          )}
        </div>
      </aside>
    );
  }

  return (
    <aside className="control-panel">
      {/* Brand */}
      <div className="cp-brand">
        <svg width="28" height="28" viewBox="0 0 28 28" fill="none" xmlns="http://www.w3.org/2000/svg">
          <rect x="2" y="2" width="24" height="24" rx="4" fill="rgba(59,130,246,0.12)" stroke="#3b82f6" strokeWidth="1.5"/>
          <path d="M8 22 L10 18 L11 16 L12 13 L14 10 L16 8 L19 5" stroke="#ef4444" strokeWidth="1.8" strokeLinecap="round" strokeLinejoin="round" fill="none"/>
          <path d="M6 22 L8 19 L9 17 L9.5 15 L10 13 L10 11 L10.5 9 L11 6" stroke="#22c55e" strokeWidth="1.8" strokeLinecap="round" strokeLinejoin="round" fill="none"/>
          <line x1="21" y1="10" x2="21" y2="22" stroke="#60a5fa" strokeWidth="1.2" strokeLinecap="round"/>
          <line x1="21" y1="10" x2="24" y2="8" stroke="#60a5fa" strokeWidth="1.2" strokeLinecap="round"/>
          <line x1="21" y1="13" x2="24" y2="11" stroke="#60a5fa" strokeWidth="1.2" strokeLinecap="round"/>
        </svg>
        <div>
          <h1 className="cp-brand-title">Sounding Analysis</h1>
          <p className="cp-brand-sub">Atmospheric Profile Tool</p>
        </div>
      </div>

      <form onSubmit={handleSubmit} className="cp-form">
        {/* Source */}
        <div className="cp-section">
          <label className="cp-label">
            <Database size={14} />
            Data Source
          </label>
          <div className="cp-source-grid">
            {sources.map((s) => (
              <div
                key={s.id}
                className="cp-source-btn-wrap"
                onMouseEnter={() => setHoveredSource(s.id)}
                onMouseLeave={() => setHoveredSource(null)}
              >
                <button
                  type="button"
                  className={`cp-source-btn ${source === s.id ? "active" : ""}`}
                  onClick={() => setSource(s.id)}
                >
                  <span className="cp-source-id">{s.id.toUpperCase()}</span>
                </button>
              </div>
            ))}
          </div>
          {SOURCE_META[hoveredSource || source] && (
            <div className="cp-source-info">
              <span className="cp-source-info-title">{SOURCE_META[hoveredSource || source].label}</span>
              <span className="cp-source-info-desc">{SOURCE_META[hoveredSource || source].desc}</span>
            </div>
          )}
        </div>

        {/* Station */}
        {needsStation && (
          <div className="cp-section">
            <label className="cp-label">
              <MapPin size={14} />
              Station
            </label>
            <div className="cp-station-picker">
              <button
                type="button"
                className="cp-risk-btn"
                onClick={handleRiskScan}
                disabled={scanning}
              >
                {scanning ? (
                  <>
                    <Loader2 size={14} className="spin" />
                    Scanning stations...
                  </>
                ) : (
                  <>
                    <Zap size={14} />
                    {riskData ? "Rescan Tornado Risk" : "Scan Tornado Risk"}
                  </>
                )}
              </button>
              {riskData && (
                <p className="cp-risk-hint">
                  Scanned {riskData.stations.length} stations at {riskData.date}
                </p>
              )}
              <div className="cp-station-toolbar">
                <div className="cp-input-wrap cp-search-flex">
                  <Search size={14} className="cp-input-icon" />
                  <input
                    type="text"
                    className="cp-input"
                    placeholder="Filter..."
                    value={stationSearch}
                    onChange={(e) => setStationSearch(e.target.value)}
                  />
                </div>
                <div className="cp-sort-wrap">
                  <ArrowUpDown size={12} className="cp-sort-icon" />
                  <select
                    className="cp-sort-select"
                    value={sortMode}
                    onChange={(e) => setSortMode(e.target.value)}
                  >
                    <option value="az">A → Z</option>
                    <option value="za">Z → A</option>
                    <option value="favs">★ Favs</option>
                    <option value="risk-high">Risk ↓</option>
                    <option value="risk-low">Risk ↑</option>
                  </select>
                </div>
              </div>
              <div className="cp-station-list" ref={listRef}>
                {filteredStations.length === 0 && (
                  <div className="cp-station-empty">No stations found</div>
                )}
                {filteredStations.map((s) => (
                  <button
                    key={s.id}
                    type="button"
                    className={`cp-station-item ${station === s.id ? "active" : ""}`}
                    onClick={() => handleStationSelect(s.id, !!riskData)}
                  >
                    <span
                      className={`cp-fav-star ${favorites.includes(s.id) ? "faved" : ""}`}
                      onClick={(e) => {
                        e.stopPropagation();
                        const { favorites: newFavs } = toggleFavorite(s.id);
                        setFavorites(newFavs);
                      }}
                      title={favorites.includes(s.id) ? "Remove from favorites" : "Add to favorites"}
                    >
                      <Star size={12} fill={favorites.includes(s.id) ? "currentColor" : "none"} />
                    </span>
                    <span className="cp-station-item-id">{s.id}</span>
                    <span className="cp-station-item-name">{s.name}</span>
                    {s.risk ? (
                      <span className={`cp-risk-score ${s.risk.stp >= 1 ? "high" : s.risk.stp >= 0.3 ? "med" : "low"}`}>
                        {s.risk.stp.toFixed(1)}
                      </span>
                    ) : (
                      <span className="cp-station-item-coords">
                        {s.lat.toFixed(1)}, {s.lon.toFixed(1)}
                      </span>
                    )}
                  </button>
                ))}
              </div>
            </div>
          </div>
        )}

        {/* Lat/Lon */}
        {needsLatLon && (
          <div className="cp-section">
            <label className="cp-label">
              <MapPin size={14} />
              Coordinates
            </label>
            <div className="cp-row">
              <div className="cp-field">
                <span className="cp-field-label">Lat</span>
                <input
                  type="number"
                  step="0.01"
                  className="cp-input cp-input-sm"
                  placeholder="35.22"
                  value={lat}
                  onChange={(e) => setLat(e.target.value)}
                  required
                />
              </div>
              <div className="cp-field">
                <span className="cp-field-label">Lon</span>
                <input
                  type="number"
                  step="0.01"
                  className="cp-input cp-input-sm"
                  placeholder="-97.46"
                  value={lon}
                  onChange={(e) => setLon(e.target.value)}
                  required
                />
              </div>
            </div>
          </div>
        )}

        {/* BUFKIT Model */}
        {needsModel && (
          <div className="cp-section">
            <label className="cp-label">
              <Layers size={14} />
              Model
            </label>
            <div className="cp-select-wrap">
              <select
                className="cp-select"
                value={model}
                onChange={(e) => setModel(e.target.value)}
              >
                {models.map((m) => (
                  <option key={m.id} value={m.id}>
                    {m.id.toUpperCase()} — {m.name}
                  </option>
                ))}
              </select>
              <ChevronDown size={14} className="cp-select-icon" />
            </div>
            <div className="cp-field" style={{ marginTop: 8 }}>
              <span className="cp-field-label">Forecast Hour</span>
              <input
                type="number"
                min="0"
                max="384"
                className="cp-input cp-input-sm"
                placeholder="0"
                value={fhour}
                onChange={(e) => setFhour(e.target.value)}
              />
            </div>
          </div>
        )}

        {/* Date */}
        <div className="cp-section">
          <label className="cp-label">
            <Calendar size={14} />
            Date / Time (UTC)
          </label>
          {source === "obs" || source === "acars" ? (
            <>
              <div className="cp-date-row">
                <input
                  type="date"
                  className="cp-input cp-calendar-input"
                  value={date ? date.slice(0, 10) : ""}
                  onChange={(e) => {
                    const d = e.target.value;
                    if (d) {
                      const hour = soundingHour === "12" ? "12:00" : "00:00";
                      setDate(`${d}T${hour}`);
                    } else {
                      setDate("");
                    }
                  }}
                />
                {date && (
                  <button
                    type="button"
                    className="cp-date-clear"
                    onClick={() => setDate("")}
                    title="Clear date"
                  >
                    ✕
                  </button>
                )}
              </div>
              <div className="cp-sounding-hour-row">
                <button
                  type="button"
                  className={`cp-hour-btn ${soundingHour === "00" ? "active" : ""}`}
                  onClick={() => {
                    setSoundingHour("00");
                    if (date) setDate(`${date.slice(0, 10)}T00:00`);
                  }}
                >
                  00Z
                </button>
                <button
                  type="button"
                  className={`cp-hour-btn ${soundingHour === "12" ? "active" : ""}`}
                  onClick={() => {
                    setSoundingHour("12");
                    if (date) setDate(`${date.slice(0, 10)}T12:00`);
                  }}
                >
                  12Z
                </button>
              </div>
              <p className="cp-hint">
                Radiosondes launch at 00Z and 12Z only
              </p>
            </>
          ) : (
            <>
              <div className="cp-date-row">
                <input
                  type="datetime-local"
                  className="cp-input cp-calendar-input"
                  value={date}
                  onChange={(e) => setDate(e.target.value)}
                />
                {date && (
                  <button
                    type="button"
                    className="cp-date-clear"
                    onClick={() => setDate("")}
                    title="Clear date"
                  >
                    ✕
                  </button>
                )}
              </div>
              <p className="cp-hint">
                Leave blank to use the most recent time
              </p>
            </>
          )}
        </div>

        {/* IGRAv2 WMO Station ID */}
        {needsWmoId && (
          <div className="cp-section">
            <label className="cp-label">
              <Database size={14} />
              WMO Station ID
            </label>
            <input
              type="text"
              className="cp-input"
              placeholder="e.g. 72451, 47646"
              value={wmoId}
              onChange={(e) => setWmoId(e.target.value)}
              required
            />
            <p className="cp-hint">
              Enter any WMO station number for global radiosonde data
            </p>
          </div>
        )}

        {/* Surface Modification */}
        <div className="cp-section">
          <button
            type="button"
            className={`cp-toggle-inline ${sfcModEnabled ? "active" : ""}`}
            onClick={() => setSfcModEnabled((v) => !v)}
          >
            <Thermometer size={12} />
            Surface Modification {sfcModEnabled ? "ON" : "OFF"}
          </button>
          {sfcModEnabled && (
            <div className="cp-mod-grid">
              <div className="cp-field">
                <span className="cp-field-label">T (°C)</span>
                <input type="number" step="0.1" className="cp-input cp-input-sm" placeholder="Sfc T" value={sfcModT} onChange={(e) => setSfcModT(e.target.value)} />
              </div>
              <div className="cp-field">
                <span className="cp-field-label">Td (°C)</span>
                <input type="number" step="0.1" className="cp-input cp-input-sm" placeholder="Sfc Td" value={sfcModTd} onChange={(e) => setSfcModTd(e.target.value)} />
              </div>
              <div className="cp-field">
                <span className="cp-field-label">Wind (kt)</span>
                <input type="number" step="1" className="cp-input cp-input-sm" placeholder="Speed" value={sfcModWspd} onChange={(e) => setSfcModWspd(e.target.value)} />
              </div>
              <div className="cp-field">
                <span className="cp-field-label">Dir (°)</span>
                <input type="number" step="1" min="0" max="360" className="cp-input cp-input-sm" placeholder="Dir" value={sfcModWdir} onChange={(e) => setSfcModWdir(e.target.value)} />
              </div>
            </div>
          )}
        </div>

        {/* Custom Storm Motion */}
        <div className="cp-section">
          <button
            type="button"
            className={`cp-toggle-inline ${smEnabled ? "active" : ""}`}
            onClick={() => setSmEnabled((v) => !v)}
          >
            <Wind size={12} />
            Custom Storm Motion {smEnabled ? "ON" : "OFF"}
          </button>
          {smEnabled && (
            <div className="cp-mod-grid">
              <div className="cp-field">
                <span className="cp-field-label">Dir (°)</span>
                <input type="number" step="1" min="0" max="360" className="cp-input cp-input-sm" placeholder="Direction" value={smDirection} onChange={(e) => setSmDirection(e.target.value)} />
              </div>
              <div className="cp-field">
                <span className="cp-field-label">Speed (kt)</span>
                <input type="number" step="1" className="cp-input cp-input-sm" placeholder="Speed" value={smSpeed} onChange={(e) => setSmSpeed(e.target.value)} />
              </div>
            </div>
          )}
        </div>

        {/* Submit */}
        <button
          type="submit"
          className="cp-submit"
          disabled={loading}
        >
          {loading ? (
            <>
              <Loader2 size={16} className="spin" />
              Fetching & Analyzing...
            </>
          ) : (
            <>
              <Search size={16} />
              Generate Sounding
            </>
          )}
        </button>

        {/* Map / Trends / History */}
        <div className="cp-toggle-row">
          <button
            type="button"
            className={`cp-toggle-btn ${showMap ? "active" : ""}`}
            onClick={onToggleMap}
          >
            <Map size={14} />
            {showMap ? "Hide Map" : "Map"}
          </button>
          <button
            type="button"
            className={`cp-toggle-btn ${showTimeSeries ? "active" : ""}`}
            onClick={onToggleTimeSeries}
          >
            <TrendingUp size={14} />
            {showTimeSeries ? "Hide" : "Trends"}
          </button>
          <button
            type="button"
            className={`cp-toggle-btn ${showCompare ? "active" : ""}`}
            onClick={onToggleCompare}
          >
            <GitCompareArrows size={14} />
            {showCompare ? "Hide" : "Compare"}
          </button>
          <button
            type="button"
            className={`cp-toggle-btn ${showHistory ? "active" : ""}`}
            onClick={onToggleHistory}
          >
            <History size={14} />
            {showHistory ? "Hide" : "History"}
          </button>
        </div>
      </form>

      {/* Footer actions */}
      <div className="cp-footer">
        <button
          type="button"
          className={`cp-footer-btn ${feedbackActive ? "active" : ""}`}
          onClick={onFeedbackClick}
          title="Send feedback"
        >
          <MessageSquarePlus size={14} />
          <span>Feedback</span>
        </button>
        <a
          href="https://github.com/ShianMike/SoundingAnalysis"
          target="_blank"
          rel="noopener noreferrer"
          className="cp-footer-btn"
          title="View on GitHub"
        >
          <Github size={14} />
          <span>GitHub</span>
        </a>
      </div>
    </aside>
  );
}
