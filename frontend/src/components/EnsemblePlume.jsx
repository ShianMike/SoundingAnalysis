import { useState } from "react";
import { Layers, X, Loader2 } from "lucide-react";
import { fetchEnsemblePlume } from "../api";
import "./EnsemblePlume.css";

const MODELS = [
  { id: "sref", label: "SREF" },
  { id: "rap", label: "RAP" },
  { id: "hrrr", label: "HRRR" },
  { id: "nam", label: "NAM" },
  { id: "namnest", label: "NAM Nest" },
  { id: "gfs", label: "GFS" },
];

const HOUR_PRESETS = {
  "Short (0-6h)": [0, 1, 2, 3, 4, 5, 6],
  "Medium (0-12h)": [0, 1, 2, 3, 6, 9, 12],
  "Long (0-24h)": [0, 3, 6, 9, 12, 15, 18, 21, 24],
  "Extended (0-48h)": [0, 6, 12, 18, 24, 30, 36, 42, 48],
};

export default function EnsemblePlume({ station, onClose, theme, colorblind }) {
  const [model, setModel] = useState("sref");
  const [hourPreset, setHourPreset] = useState("Medium (0-12h)");
  const [source, setSource] = useState("bufkit");
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

  return (
    <div className="ens-panel">
      <div className="ens-header">
        <span className="ens-title"><Layers size={14} /> Ensemble Sounding Plume</span>
        <button className="ens-close" onClick={onClose}><X size={14} /></button>
      </div>

      <div className="ens-controls">
        <label>
          Station:
          <input type="text" value={station || "OUN"} readOnly style={{ width: 50 }} />
        </label>
        <label>
          Model:
          <select value={model} onChange={(e) => setModel(e.target.value)}>
            {MODELS.map((m) => (
              <option key={m.id} value={m.id}>{m.label}</option>
            ))}
          </select>
        </label>
        <label>
          Source:
          <select value={source} onChange={(e) => setSource(e.target.value)}>
            <option value="bufkit">Iowa State</option>
            <option value="psu">Penn State</option>
          </select>
        </label>
        <label>
          Hours:
          <select value={hourPreset} onChange={(e) => setHourPreset(e.target.value)}>
            {Object.keys(HOUR_PRESETS).map((k) => (
              <option key={k} value={k}>{k}</option>
            ))}
          </select>
        </label>
        <button className="ens-fetch-btn" onClick={handleFetch} disabled={loading}>
          {loading ? <><Loader2 size={12} className="spin" /> Loading...</> : "Generate Plume"}
        </button>
      </div>

      {error && <div className="ens-error">{error}</div>}

      {result && (
        <>
          <div className="ens-plot-wrap">
            <img
              src={`data:image/png;base64,${result.image}`}
              alt="Ensemble sounding plume"
              className="ens-plot-img"
            />
          </div>
          <div className="ens-info">
            {result.members} members loaded
            {result.hours && ` • f${result.hours.join(", f")}`}
            {result.errors?.length > 0 && ` • ${result.errors.length} failed`}
          </div>
        </>
      )}
    </div>
  );
}
