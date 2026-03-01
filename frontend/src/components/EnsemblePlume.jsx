import { useState } from "react";
import { Layers, ArrowLeft, Loader2, AlertTriangle, Info } from "lucide-react";
import { fetchEnsemblePlume } from "../api";
import "./EnsemblePlume.css";

const MODELS = [
  { id: "rap", label: "RAP", desc: "Hourly, 13 km" },
  { id: "hrrr", label: "HRRR", desc: "Hourly, 3 km" },
  { id: "nam", label: "NAM", desc: "Hourly, 12 km" },
  { id: "namnest", label: "NAM Nest", desc: "Hourly, 3 km" },
  { id: "gfs", label: "GFS", desc: "3-hourly, global" },
  { id: "sref", label: "SREF", desc: "Ensemble (legacy)" },
];

const HOUR_PRESETS = {
  "Short (0-6h)": [0, 1, 2, 3, 4, 5, 6],
  "Medium (0-12h)": [0, 1, 2, 3, 6, 9, 12],
  "Long (0-24h)": [0, 3, 6, 9, 12, 15, 18, 21, 24],
  "Extended (0-48h)": [0, 6, 12, 18, 24, 30, 36, 42, 48],
};

export default function EnsemblePlume({ station, onBack, theme, colorblind }) {
  const [model, setModel] = useState("rap");
  const [hourPreset, setHourPreset] = useState("Medium (0-12h)");
  const [source, setSource] = useState("psu");
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState(null);
  const [error, setError] = useState(null);

  const handleFetch = async () => {
    setLoading(true);
    setError(null);
    setResult(null);
    try {
      const hours = HOUR_PRESETS[hourPreset] || [0, 1, 2, 3, 6, 9, 12];
      const data = await fetchEnsemblePlume({
        station: station || "OUN",
        model,
        source,
        hours,
        theme: theme || "dark",
        colorblind: colorblind || false,
      });
      setResult(data);
    } catch (e) {
      setError(e.message);
    } finally {
      setLoading(false);
    }
  };

  const modelInfo = MODELS.find((m) => m.id === model);

  return (
    <div className="ens-page">
      <div className="ens-page-inner">
        {/* Header */}
        <div className="ens-header">
          <div className="ens-header-left">
            <button className="ens-back" onClick={onBack} title="Back to main"><ArrowLeft size={16} /></button>
            <Layers size={15} className="ens-header-icon" />
            <div>
              <span className="ens-title">Ensemble Sounding Plume</span>
              <span className="ens-subtitle">
                {station || "OUN"} &middot; {modelInfo?.label || model.toUpperCase()}
              </span>
            </div>
          </div>
        </div>

      {/* Controls */}
      <div className="ens-controls">
        <div className="ens-ctrl-group">
          <span className="ens-ctrl-label">Model</span>
          <select value={model} onChange={(e) => setModel(e.target.value)}>
            {MODELS.map((m) => (
              <option key={m.id} value={m.id}>{m.label}</option>
            ))}
          </select>
        </div>
        <div className="ens-ctrl-group">
          <span className="ens-ctrl-label">Source</span>
          <select value={source} onChange={(e) => setSource(e.target.value)}>
            <option value="psu">Penn State (latest)</option>
            <option value="bufkit">Iowa State (archive)</option>
          </select>
        </div>
        <div className="ens-ctrl-group">
          <span className="ens-ctrl-label">Fcst Range</span>
          <select value={hourPreset} onChange={(e) => setHourPreset(e.target.value)}>
            {Object.keys(HOUR_PRESETS).map((k) => (
              <option key={k} value={k}>{k}</option>
            ))}
          </select>
        </div>
        <button className="ens-fetch-btn" onClick={handleFetch} disabled={loading}>
          {loading ? <><Loader2 size={13} className="spin" /> Generating...</> : "Generate Plume"}
        </button>
      </div>

      {/* Hint for SREF */}
      {model === "sref" && !result && !loading && (
        <div className="ens-hint">
          <Info size={12} />
          SREF was discontinued. Data may be unavailable for recent dates. Consider using RAP or HRRR instead.
        </div>
      )}

      {/* Error */}
      {error && (
        <div className="ens-error">
          <AlertTriangle size={14} />
          <div className="ens-error-text">{error}</div>
        </div>
      )}

      {/* Results */}
      {result && (
        <div className="ens-result">
          <div className="ens-plot-wrap">
            <img
              src={`data:image/png;base64,${result.image}`}
              alt="Ensemble sounding plume"
              className="ens-plot-img"
            />
          </div>
          <div className="ens-meta-bar">
            <span className="ens-meta-badge">{result.members} members</span>
            <span className="ens-meta-badge">f{Math.min(...result.hours).toString().padStart(3,"0")} – f{Math.max(...result.hours).toString().padStart(3,"0")}</span>
            {result.meta?.source && <span className="ens-meta-badge">{result.meta.source === "psu" ? "Penn State" : "Iowa State"}</span>}
            {result.meta?.date && <span className="ens-meta-badge">{result.meta.date}</span>}
            {result.errors?.length > 0 && <span className="ens-meta-badge ens-meta-warn">{result.errors.length} failed</span>}
          </div>
        </div>
      )}
      </div>
    </div>
  );
}
