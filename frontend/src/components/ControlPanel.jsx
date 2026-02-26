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
  selectedStation,
  onStationChange,
  onSourceChange,
  mapLatLon,
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
  const listRef = useRef(null);

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
  const needsStation = source === "obs" || source === "bufkit" || source === "acars";
  const needsModel = source === "bufkit";

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
    if (needsLatLon) {
      params.lat = parseFloat(lat);
      params.lon = parseFloat(lon);
    }
    if (date) params.date = date.replace(/[-T:]/g, "").slice(0, 10);
    if (needsModel) {
      params.model = model;
      params.fhour = parseInt(fhour) || 0;
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
                {hoveredSource === s.id && (
                  <div className="cp-source-card">
                    <span className="cp-source-card-title">{SOURCE_META[s.id]?.label}</span>
                    <span className="cp-source-card-desc">{SOURCE_META[s.id]?.desc}</span>
                  </div>
                )}
              </div>
            ))}
          </div>
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
            Leave blank to use the most recent sounding time
          </p>
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
            className={`cp-toggle-btn ${showHistory ? "active" : ""}`}
            onClick={onToggleHistory}
          >
            <History size={14} />
            {showHistory ? "Hide" : "History"}
          </button>
        </div>
      </form>
    </aside>
  );
}
