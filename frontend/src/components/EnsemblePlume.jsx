import { useState, useMemo, useRef, useEffect } from "react";
import { Layers, ArrowLeft, Loader2, AlertTriangle, Info, RefreshCw, Lightbulb, HelpCircle, Zap, MapPin, Search, X } from "lucide-react";
import { fetchEnsemblePlume } from "../api";
import "./EnsemblePlume.css";

const MODELS = [
  { id: "rap", label: "RAP", desc: "Hourly analysis, 13 km – best overall coverage" },
  { id: "hrrr", label: "HRRR", desc: "Hourly, 3 km – high-res CONUS" },
  { id: "nam", label: "NAM", desc: "Hourly, 12 km" },
  { id: "namnest", label: "NAM Nest", desc: "Hourly, 3 km" },
  { id: "gfs", label: "GFS", desc: "3-hourly, global" },
  { id: "sref", label: "SREF ⚠", desc: "Discontinued — data unavailable" },
];

const SOURCES = [
  { id: "psu", label: "Penn State (latest)", desc: "Most recent model run" },
  { id: "bufkit", label: "Iowa State (archive)", desc: "Historical archive" },
];

const HOUR_PRESETS = [
  { id: "short",    label: "Short (0–6 h)",    hours: [0, 1, 2, 3, 4, 5, 6],          desc: "Near-term evolution" },
  { id: "medium",   label: "Medium (0–12 h)",   hours: [0, 1, 2, 3, 6, 9, 12],         desc: "Half-day outlook" },
  { id: "long",     label: "Long (0–24 h)",     hours: [0, 3, 6, 9, 12, 15, 18, 21, 24], desc: "Full day outlook" },
  { id: "extended", label: "Extended (0–48 h)",  hours: [0, 6, 12, 18, 24, 30, 36, 42, 48], desc: "2-day outlook" },
];

/* Stations with reliable BUFKIT data on Penn State */
const QUICK_COMBOS = [
  { station: "OUN", model: "rap", source: "psu", label: "OUN", region: "Oklahoma" },
  { station: "TBW", model: "rap", source: "psu", label: "TBW", region: "Tampa Bay" },
  { station: "ILX", model: "rap", source: "psu", label: "ILX", region: "Lincoln, IL" },
  { station: "OAX", model: "rap", source: "psu", label: "OAX", region: "Omaha" },
  { station: "SGF", model: "rap", source: "psu", label: "SGF", region: "Springfield" },
  { station: "DVN", model: "rap", source: "psu", label: "DVN", region: "Davenport" },
];

export default function EnsemblePlume({ station, stations = [], onBack, theme, colorblind }) {
  const [localStation, setLocalStation] = useState(station || "OUN");
  const [model, setModel] = useState("rap");
  const [presetId, setPresetId] = useState("medium");
  const [source, setSource] = useState("psu");
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState(null);
  const [error, setError] = useState(null);
  const [suggestions, setSuggestions] = useState([]);
  const [showHelp, setShowHelp] = useState(false);
  const [stationSearch, setStationSearch] = useState("");
  const [stationDropdownOpen, setStationDropdownOpen] = useState(false);
  const stationRef = useRef(null);

  // Close dropdown on outside click
  useEffect(() => {
    const handler = (e) => {
      if (stationRef.current && !stationRef.current.contains(e.target)) {
        setStationDropdownOpen(false);
      }
    };
    document.addEventListener("mousedown", handler);
    return () => document.removeEventListener("mousedown", handler);
  }, []);

  // Sync if parent station changes
  useEffect(() => {
    if (station) setLocalStation(station);
  }, [station]);

  const filteredStations = useMemo(() => {
    if (!stations.length) return [];
    const q = stationSearch.trim().toUpperCase();
    if (!q) return stations.slice(0, 40);
    return stations.filter((s) => {
      const id = (s.id || "").toUpperCase();
      const name = (s.name || s.label || "").toUpperCase();
      return id.includes(q) || name.includes(q);
    }).slice(0, 40);
  }, [stations, stationSearch]);

  const preset = HOUR_PRESETS.find((p) => p.id === presetId) || HOUR_PRESETS[1];

  const handleFetch = async (overrides = {}) => {
    setLoading(true);
    setError(null);
    setSuggestions([]);
    setResult(null);
    const usedStation = overrides.station || localStation;
    const usedModel = overrides.model || model;
    const usedSource = overrides.source || source;
    const usedHours = preset.hours;
    try {
      const data = await fetchEnsemblePlume({
        station: usedStation,
        model: usedModel,
        source: usedSource,
        hours: usedHours,
        theme: theme || "dark",
        colorblind: colorblind || false,
      });
      setResult(data);
    } catch (e) {
      setError(e.message);
      setSuggestions(e.suggestions || []);
    } finally {
      setLoading(false);
    }
  };

  const handleQuickFix = (action) => {
    if (action === "switch-source") {
      const newSrc = source === "psu" ? "bufkit" : "psu";
      setSource(newSrc);
      handleFetch({ source: newSrc });
    } else if (action === "try-rap") {
      setModel("rap");
      handleFetch({ model: "rap" });
    } else if (action === "try-hrrr") {
      setModel("hrrr");
      handleFetch({ model: "hrrr" });
    }
  };

  const handleQuickCombo = (combo) => {
    setLocalStation(combo.station);
    setModel(combo.model);
    setSource(combo.source);
    handleFetch({ model: combo.model, source: combo.source, station: combo.station });
  };

  const modelInfo = MODELS.find((m) => m.id === model);
  const sourceInfo = SOURCES.find((s) => s.id === source);

  return (
    <div className="ens-page">
      <div className="ens-page-inner">
        {/* Header card */}
        <div className="ens-header-card">
          <div className="ens-header">
            <div className="ens-header-left">
              <button className="ens-back" onClick={onBack} title="Back to main"><ArrowLeft size={16} /></button>
              <div className="ens-header-icon-wrap"><Layers size={18} className="ens-header-icon" /></div>
              <div>
                <span className="ens-title">Forecast Sounding Plume</span>
                <span className="ens-subtitle">
                  Overlay multiple forecast hours on one Skew-T
                </span>
              </div>
            </div>
            <button className="ens-help-btn" onClick={() => setShowHelp((v) => !v)} title="How it works">
              <HelpCircle size={15} />
            </button>
          </div>

          {/* Help explainer */}
          {showHelp && (
            <div className="ens-help-box">
              <p><strong>What is this?</strong> A forecast plume overlays soundings from multiple forecast hours (e.g. f000, f003, f006…) onto one diagram, showing how the atmosphere is predicted to evolve over time.</p>
              <p><strong>How to use:</strong> Pick a station, model, data source, and forecast range, then click "Generate Plume". The first forecast hour (f000 = analysis) is drawn in bold; later hours are semi-transparent to show the spread.</p>
              <p><strong>Source tips:</strong> "Penn State (latest)" has the most recent model run and works best for current data. "Iowa State (archive)" is better for past dates.</p>
            </div>
          )}
        </div>

        {/* Controls card — grid layout */}
        <div className="ens-controls-card">
          <div className="ens-controls-grid">
            <div className="ens-ctrl-cell" ref={stationRef}>
              <span className="ens-ctrl-label">Station</span>
              <button
                className="ens-station-btn"
                onClick={() => setStationDropdownOpen((v) => !v)}
                type="button"
              >
                <MapPin size={12} className="ens-station-icon" />
                {localStation.toUpperCase()}
                <span className="ens-station-arrow">▾</span>
              </button>
              {stationDropdownOpen && (
                <div className="ens-station-dropdown">
                  <div className="ens-station-search-wrap">
                    <Search size={12} />
                    <input
                      className="ens-station-search"
                      placeholder="Search stations…"
                      value={stationSearch}
                      onChange={(e) => setStationSearch(e.target.value)}
                      autoFocus
                    />
                    {stationSearch && (
                      <button className="ens-station-search-clear" onClick={() => setStationSearch("")}><X size={10} /></button>
                    )}
                  </div>
                  <div className="ens-station-list">
                    {filteredStations.length === 0 && (
                      <div className="ens-station-empty">No stations found</div>
                    )}
                    {filteredStations.map((s) => (
                      <button
                        key={s.id}
                        className={`ens-station-option${s.id === localStation ? " active" : ""}`}
                        onClick={() => {
                          setLocalStation(s.id);
                          setStationDropdownOpen(false);
                          setStationSearch("");
                        }}
                      >
                        <span className="ens-station-option-id">{s.id}</span>
                        {(s.name || s.label) && <span className="ens-station-option-name">{s.name || s.label}</span>}
                      </button>
                    ))}
                  </div>
                </div>
              )}
              <span className="ens-ctrl-hint">{stations.length > 0 ? `${stations.length} stations available` : "Type to search"}</span>
            </div>
            <div className="ens-ctrl-cell">
              <span className="ens-ctrl-label">Model</span>
              <select value={model} onChange={(e) => setModel(e.target.value)}>
                {MODELS.map((m) => (
                  <option key={m.id} value={m.id}>{m.label}</option>
                ))}
              </select>
              <span className="ens-ctrl-hint">{modelInfo?.desc}</span>
            </div>
            <div className="ens-ctrl-cell">
              <span className="ens-ctrl-label">Source</span>
              <select value={source} onChange={(e) => setSource(e.target.value)}>
                {SOURCES.map((s) => (
                  <option key={s.id} value={s.id}>{s.label}</option>
                ))}
              </select>
              <span className="ens-ctrl-hint">{sourceInfo?.desc}</span>
            </div>
            <div className="ens-ctrl-cell">
              <span className="ens-ctrl-label">Forecast Range</span>
              <select value={presetId} onChange={(e) => setPresetId(e.target.value)}>
                {HOUR_PRESETS.map((p) => (
                  <option key={p.id} value={p.id}>{p.label}</option>
                ))}
              </select>
              <span className="ens-ctrl-hint">{preset.hours.length} forecast hours · {preset.desc}</span>
            </div>
          </div>
          <div className="ens-controls-action">
            <button className="ens-fetch-btn" onClick={() => handleFetch()} disabled={loading}>
              {loading ? <><Loader2 size={14} className="spin" /> Generating…</> : <><Zap size={14} /> Generate Plume</>}
            </button>
          </div>

          {/* SREF warning */}
          {model === "sref" && !result && !loading && (
            <div className="ens-hint ens-hint-warn">
              <AlertTriangle size={13} />
              <span>SREF was discontinued in 2025. Data is unavailable for recent dates. <button className="ens-link-btn" onClick={() => handleQuickFix("try-rap")}>Switch to RAP</button></span>
            </div>
          )}
        </div>

        {/* Error */}
        {error && (
          <div className="ens-error-card">
            <div className="ens-error-header">
              <AlertTriangle size={18} />
              <span>Data Not Available</span>
            </div>
            <p className="ens-error-msg">{error}</p>

            {/* Quick-fix actions */}
            <div className="ens-error-actions">
              <div className="ens-error-actions-row">
                <button className="ens-suggestion-btn" onClick={() => handleQuickFix("switch-source")}>
                  <RefreshCw size={12} /> Switch to {source === "psu" ? "Iowa State" : "Penn State"}
                </button>
                {model !== "rap" && (
                  <button className="ens-suggestion-btn" onClick={() => handleQuickFix("try-rap")}>
                    Try RAP model
                  </button>
                )}
                {model !== "hrrr" && (
                  <button className="ens-suggestion-btn" onClick={() => handleQuickFix("try-hrrr")}>
                    Try HRRR model
                  </button>
                )}
              </div>
            </div>

            {/* Quick combos — known working */}
            <div className="ens-quick-combos">
              <div className="ens-quick-combos-label"><Zap size={11} /> Try a station that's available now:</div>
              <div className="ens-quick-combos-grid">
                {QUICK_COMBOS.map((c) => (
                  <button key={c.station} className="ens-combo-chip" onClick={() => handleQuickCombo(c)}>
                    <span className="ens-combo-station">{c.label}</span>
                    <span className="ens-combo-region">{c.region}</span>
                  </button>
                ))}
              </div>
            </div>

            {suggestions.length > 0 && (
              <ul className="ens-suggestion-list">
                {suggestions.map((s, i) => <li key={i}>{s}</li>)}
              </ul>
            )}
          </div>
        )}

        {/* Loading state */}
        {loading && (
          <div className="ens-loading">
            <div className="ens-loading-spinner"><Loader2 size={22} className="spin" /></div>
            <p>Fetching {preset.hours.length} forecast hours for {localStation.toUpperCase()} ({model.toUpperCase()})…</p>
            <span className="ens-loading-sub">This typically takes 10–30 seconds</span>
            <div className="ens-loading-bar"><div className="ens-loading-bar-fill" /></div>
          </div>
        )}

        {/* Empty state — show when no result, no error, not loading */}
        {!result && !error && !loading && (
          <div className="ens-empty-state">
            <div className="ens-empty-icon"><Layers size={24} /></div>
            <div className="ens-empty-title">No plume generated yet</div>
            <div className="ens-empty-desc">Configure the settings above and click "Generate Plume" to overlay forecast soundings on one Skew-T diagram.</div>
          </div>
        )}

        {/* Results */}
        {result && (
          <div className="ens-result-card">
            <div className="ens-result">
              <div className="ens-plot-wrap">
                <img
                  src={`data:image/png;base64,${result.image}`}
                  alt="Forecast sounding plume"
                  className="ens-plot-img"
                />
              </div>
              <div className="ens-meta-bar">
                <span className="ens-meta-badge">{result.members} forecast hours</span>
                <span className="ens-meta-badge">f{Math.min(...result.hours).toString().padStart(3,"0")} – f{Math.max(...result.hours).toString().padStart(3,"0")}</span>
                {result.meta?.source && <span className="ens-meta-badge">{result.meta.source === "psu" ? "Penn State" : "Iowa State"}</span>}
                {result.meta?.date && <span className="ens-meta-badge">{result.meta.date}</span>}
                {result.errors?.length > 0 && <span className="ens-meta-badge ens-meta-warn">{result.errors.length} hours failed</span>}
              </div>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
