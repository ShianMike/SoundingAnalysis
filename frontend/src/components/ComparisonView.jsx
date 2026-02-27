import { useState, useEffect } from "react";
import {
  GitCompareArrows,
  Plus,
  X,
  Loader2,
  AlertTriangle,
  MapPin,
  Calendar,
  Database,
  ChevronDown,
  Search,
  Download,
} from "lucide-react";
import { fetchCompare, fetchComposite } from "../api";
import { saveCompareToHistory } from "../history";
import "./ComparisonView.css";

/* ── parameter labels for the comparison table ────────────────── */
const PARAM_GROUPS = [
  {
    title: "Severe Indices",
    params: [
      { key: "stp", label: "STP" },
      { key: "scp", label: "SCP" },
      { key: "ship", label: "SHIP" },
      { key: "dcp", label: "DCP" },
      { key: "ecape", label: "ECAPE", unit: "J/kg" },
    ],
  },
  {
    title: "CAPE / CIN",
    params: [
      { key: "sbCape", label: "SB CAPE", unit: "J/kg" },
      { key: "sbCin", label: "SB CIN", unit: "J/kg" },
      { key: "muCape", label: "MU CAPE", unit: "J/kg" },
      { key: "mlCape", label: "ML CAPE", unit: "J/kg" },
      { key: "mlCin", label: "ML CIN", unit: "J/kg" },
      { key: "dcape", label: "DCAPE", unit: "J/kg" },
    ],
  },
  {
    title: "Shear / Helicity",
    params: [
      { key: "bwd1km", label: "BWD 0-1 km", unit: "kt" },
      { key: "bwd3km", label: "BWD 0-3 km", unit: "kt" },
      { key: "bwd6km", label: "BWD 0-6 km", unit: "kt" },
      { key: "srh500m", label: "SRH 500m", unit: "m²/s²" },
      { key: "srh1km", label: "SRH 0-1 km", unit: "m²/s²" },
      { key: "srh3km", label: "SRH 0-3 km", unit: "m²/s²" },
    ],
  },
  {
    title: "Lapse Rates / Moisture",
    params: [
      { key: "lr03", label: "LR 0-3 km", unit: "°C/km" },
      { key: "lr36", label: "LR 3-6 km", unit: "°C/km" },
      { key: "pwat", label: "PWAT", unit: "mm" },
      { key: "frzLevel", label: "FRZ Level", unit: "m" },
      { key: "rh01", label: "RH 0-1 km", unit: "%" },
    ],
  },
];

/* ── color palette for each slot ──────────────────────────────── */
const SLOT_COLORS = ["#60a5fa", "#f59e0b", "#10b981", "#a78bfa"];

function SlotCard({ slot, index, stations, onUpdate, onRemove, canRemove }) {
  const [search, setSearch] = useState("");

  const filtered = stations.filter(
    (s) =>
      s.id.toLowerCase().includes(search.toLowerCase()) ||
      s.name.toLowerCase().includes(search.toLowerCase())
  );

  return (
    <div className="cv-slot" style={{ borderColor: SLOT_COLORS[index] }}>
      <div className="cv-slot-header">
        <span className="cv-slot-badge" style={{ background: SLOT_COLORS[index] }}>
          {index + 1}
        </span>
        <span className="cv-slot-title">
          {slot.station || "Select station"}
        </span>
        {canRemove && (
          <button className="cv-slot-remove" onClick={onRemove} title="Remove">
            <X size={14} />
          </button>
        )}
      </div>

      {/* Source selector */}
      <div className="cv-slot-field">
        <Database size={12} />
        <select
          className="cv-slot-select"
          value={slot.source}
          onChange={(e) => onUpdate({ ...slot, source: e.target.value })}
        >
          <option value="obs">OBS</option>
          <option value="rap">RAP</option>
          <option value="bufkit">BUFKIT</option>
          <option value="acars">ACARS</option>
        </select>
      </div>

      {/* Station picker */}
      <div className="cv-slot-field">
        <MapPin size={12} />
        <div className="cv-station-pick">
          <input
            type="text"
            className="cv-slot-input"
            placeholder="Filter station..."
            value={search}
            onChange={(e) => setSearch(e.target.value)}
          />
          <div className="cv-station-dropdown">
            {filtered.slice(0, 30).map((s) => (
              <button
                key={s.id}
                type="button"
                className={`cv-station-opt ${slot.station === s.id ? "active" : ""}`}
                onClick={() => {
                  onUpdate({ ...slot, station: s.id });
                  setSearch("");
                }}
              >
                <span className="cv-opt-id">{s.id}</span>
                <span className="cv-opt-name">{s.name}</span>
              </button>
            ))}
          </div>
        </div>
      </div>

      {/* Date picker */}
      <div className="cv-slot-field">
        <Calendar size={12} />
        {slot.source === "obs" || slot.source === "acars" ? (
          <div className="cv-date-obs">
            <input
              type="date"
              className="cv-slot-input"
              value={slot.date ? slot.date.slice(0, 10) : ""}
              onChange={(e) => {
                const d = e.target.value;
                if (d) {
                  onUpdate({ ...slot, date: `${d}T${slot.hour || "12"}:00` });
                } else {
                  onUpdate({ ...slot, date: "" });
                }
              }}
            />
            <div className="cv-hour-toggle">
              <button
                type="button"
                className={`cv-hour-btn ${(slot.hour || "12") === "00" ? "active" : ""}`}
                onClick={() => {
                  const d = slot.date ? slot.date.slice(0, 10) : "";
                  onUpdate({ ...slot, hour: "00", date: d ? `${d}T00:00` : "" });
                }}
              >
                00Z
              </button>
              <button
                type="button"
                className={`cv-hour-btn ${(slot.hour || "12") === "12" ? "active" : ""}`}
                onClick={() => {
                  const d = slot.date ? slot.date.slice(0, 10) : "";
                  onUpdate({ ...slot, hour: "12", date: d ? `${d}T12:00` : "" });
                }}
              >
                12Z
              </button>
            </div>
          </div>
        ) : (
          <input
            type="datetime-local"
            className="cv-slot-input"
            value={slot.date || ""}
            onChange={(e) => onUpdate({ ...slot, date: e.target.value })}
          />
        )}
      </div>
    </div>
  );
}

/* For highlighting best/worst values */
function getHighlightClass(key, value, allValues) {
  if (value == null || allValues.length < 2) return "";
  const nums = allValues.filter((v) => v != null);
  if (nums.length < 2) return "";
  const max = Math.max(...nums);
  const min = Math.min(...nums);
  if (max === min) return "";

  // Higher is "more significant" for these params
  const higherIsBolder = [
    "stp", "scp", "ship", "dcp", "ecape", "sbCape", "muCape", "mlCape", "dcape",
    "bwd1km", "bwd3km", "bwd6km", "srh500m", "srh1km", "srh3km",
    "lr03", "lr36", "pwat",
  ];

  if (higherIsBolder.includes(key)) {
    if (value === max) return "cv-cell-high";
    if (value === min) return "cv-cell-low";
  } else {
    if (value === min) return "cv-cell-high";
    if (value === max) return "cv-cell-low";
  }
  return "";
}

export default function ComparisonView({ stations, onClose, historyData, onHistoryConsumed }) {
  const makeSlot = (station = "", source = "obs", date = "", hour = "12") => ({
    source,
    station,
    date,
    hour,
  });

  const [slots, setSlots] = useState([makeSlot("OUN"), makeSlot("FWD")]);
  const [results, setResults] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [compositeImage, setCompositeImage] = useState(null);
  const [compositeLoading, setCompositeLoading] = useState(false);

  // Load from history if provided
  useEffect(() => {
    if (historyData && historyData.slots && historyData.results) {
      setSlots(historyData.slots.map((s) => makeSlot(s.station, s.source, s.date, s.hour)));
      setResults(historyData.results);
      setError(null);
      if (onHistoryConsumed) onHistoryConsumed();
    }
  }, [historyData]);

  const updateSlot = (idx, slot) => {
    setSlots((prev) => prev.map((s, i) => (i === idx ? slot : s)));
  };

  const removeSlot = (idx) => {
    setSlots((prev) => prev.filter((_, i) => i !== idx));
  };

  const addSlot = () => {
    if (slots.length < 4) setSlots((prev) => [...prev, makeSlot()]);
  };

  const handleCompare = async () => {
    const valid = slots.filter((s) => s.station);
    if (valid.length < 2) {
      setError("Select at least 2 stations to compare.");
      return;
    }

    setLoading(true);
    setError(null);
    setResults(null);

    try {
      const soundings = valid.map((s) => {
        const params = { source: s.source, station: s.station };
        if (s.date) {
          params.date = s.date.replace(/[-T:]/g, "").slice(0, 10);
        }
        return params;
      });
      const data = await fetchCompare(soundings);
      setResults(data.soundings);
      // Save to history
      saveCompareToHistory(valid, data.soundings);
    } catch (e) {
      setError(e.message);
    } finally {
      setLoading(false);
    }
  };

  const handleComposite = async () => {
    const valid = slots.filter((s) => s.station);
    if (valid.length < 2) {
      setError("Select at least 2 stations for the overlay.");
      return;
    }

    setCompositeLoading(true);
    setError(null);
    setCompositeImage(null);

    try {
      const soundings = valid.map((s) => {
        const params = { source: s.source, station: s.station };
        if (s.date) {
          params.date = s.date.replace(/[-T:]/g, "").slice(0, 10);
        }
        return params;
      });
      const data = await fetchComposite(soundings);
      setCompositeImage(data.image);
    } catch (e) {
      setError(e.message);
    } finally {
      setCompositeLoading(false);
    }
  };

  /* Download comparison as a composite image */
  const handleDownload = () => {
    if (!results) return;
    const validResults = results.filter((r) => !r.error && r.image);
    if (validResults.length === 0) return;

    // Load all images, then tile them side-by-side on a canvas
    const imgs = validResults.map((r) => {
      const img = new Image();
      img.src = `data:image/png;base64,${r.image}`;
      return img;
    });

    Promise.all(imgs.map((img) => new Promise((res) => {
      if (img.complete) return res(img);
      img.onload = () => res(img);
      img.onerror = () => res(img);
    }))).then((loaded) => {
      const maxH = Math.max(...loaded.map((i) => i.naturalHeight || 800));
      const totalW = loaded.reduce((s, i) => s + (i.naturalWidth || 600), 0);
      const canvas = document.createElement("canvas");
      canvas.width = totalW;
      canvas.height = maxH;
      const ctx = canvas.getContext("2d");
      ctx.fillStyle = "#0e1117";
      ctx.fillRect(0, 0, totalW, maxH);
      let x = 0;
      loaded.forEach((img) => {
        const w = img.naturalWidth || 600;
        const h = img.naturalHeight || 800;
        ctx.drawImage(img, x, 0, w, h);
        x += w;
      });
      const link = document.createElement("a");
      const stations = validResults.map((r) => r.meta?.station || "?").join("_vs_");
      link.download = `comparison_${stations}.png`;
      link.href = canvas.toDataURL("image/png");
      link.click();
    });
  };

  return (
    <div className="cv-wrap">
      {/* Header */}
      <div className="cv-header">
        <div className="cv-header-left">
          <GitCompareArrows size={16} />
          <h3>Multi-Sounding Comparison</h3>
        </div>
        <button className="cv-close-btn" onClick={onClose}>
          <X size={16} />
        </button>
      </div>

      {/* Slot cards */}
      <div className="cv-slots">
        {slots.map((slot, i) => (
          <SlotCard
            key={i}
            slot={slot}
            index={i}
            stations={stations}
            onUpdate={(s) => updateSlot(i, s)}
            onRemove={() => removeSlot(i)}
            canRemove={slots.length > 2}
          />
        ))}
        {slots.length < 4 && (
          <button className="cv-add-slot" onClick={addSlot}>
            <Plus size={16} />
            Add Sounding
          </button>
        )}
      </div>

      {/* Compare and Composite buttons */}
      <div className="cv-action-row" style={{ display: "flex", gap: "8px" }}>
        <button
          className="cv-compare-btn"
          onClick={handleCompare}
          disabled={loading || compositeLoading || slots.filter((s) => s.station).length < 2}
        >
          {loading ? (
            <>
              <Loader2 size={16} className="spin" />
              Fetching soundings...
            </>
          ) : (
            <>
              <GitCompareArrows size={16} />
              Compare ({slots.filter((s) => s.station).length})
            </>
          )}
        </button>
        <button
          className="cv-compare-btn cv-composite-btn"
          onClick={handleComposite}
          disabled={loading || compositeLoading || slots.filter((s) => s.station).length < 2}
          title="Overlay all profiles on a single Skew-T"
        >
          {compositeLoading ? (
            <>
              <Loader2 size={16} className="spin" />
              Generating overlay...
            </>
          ) : (
            <>
              <GitCompareArrows size={16} />
              Overlay
            </>
          )}
        </button>
      </div>

      {error && (
        <div className="cv-error">
          <AlertTriangle size={14} />
          <span>{error}</span>
        </div>
      )}

      {/* Results */}
      {results && (
        <div className="cv-results">
          {/* Download button */}
          <div className="cv-results-actions">
            <button className="cv-download-btn" onClick={handleDownload} title="Download comparison as PNG">
              <Download size={14} />
              Download Comparison
            </button>
          </div>

          {/* Side-by-side plots */}
          <div className="cv-plots" style={{ gridTemplateColumns: `repeat(${results.length}, 1fr)` }}>
            {results.map((r, i) => (
              <div key={i} className="cv-plot-card" style={{ borderTopColor: SLOT_COLORS[i] }}>
                {r.error ? (
                  <div className="cv-plot-error">
                    <AlertTriangle size={16} />
                    <span>{r.error}</span>
                  </div>
                ) : (
                  <>
                    <div className="cv-plot-meta">
                      <span className="cv-plot-station" style={{ color: SLOT_COLORS[i] }}>
                        {r.meta.station}
                      </span>
                      <span className="cv-plot-date">{r.meta.date}</span>
                      <span className="cv-plot-src">{r.meta.source.toUpperCase()}</span>
                    </div>
                    <img
                      src={`data:image/png;base64,${r.image}`}
                      alt={`Sounding ${r.meta.station}`}
                      className="cv-plot-img"
                    />
                  </>
                )}
              </div>
            ))}
          </div>

          {/* Parameter comparison table */}
          <div className="cv-table-wrap">
            <table className="cv-table">
              <thead>
                <tr>
                  <th className="cv-th-param">Parameter</th>
                  {results.map((r, i) => (
                    <th key={i} style={{ color: SLOT_COLORS[i] }}>
                      {r.meta?.station || r.meta?.source || `#${i + 1}`}
                    </th>
                  ))}
                  {results.length >= 2 && results.every((r) => !r.error) && (
                    <th className="cv-th-diff">Δ</th>
                  )}
                </tr>
              </thead>
              <tbody>
                {PARAM_GROUPS.map((group) => (
                  <>
                    <tr key={group.title} className="cv-group-row">
                      <td colSpan={results.length + 2} className="cv-group-label">
                        {group.title}
                      </td>
                    </tr>
                    {group.params.map((p) => {
                      const values = results.map((r) => (r.params ? r.params[p.key] : null));
                      const nums = values.filter((v) => v != null);
                      const diff =
                        nums.length >= 2 ? (Math.max(...nums) - Math.min(...nums)).toFixed(1) : null;

                      return (
                        <tr key={p.key}>
                          <td className="cv-td-label">
                            {p.label}
                            {p.unit && <span className="cv-td-unit">{p.unit}</span>}
                          </td>
                          {results.map((r, i) => {
                            const val = r.params ? r.params[p.key] : null;
                            return (
                              <td
                                key={i}
                                className={`cv-td-val ${getHighlightClass(p.key, val, values)}`}
                              >
                                {val != null ? val : "—"}
                              </td>
                            );
                          })}
                          {results.length >= 2 && results.every((r) => !r.error) && (
                            <td className="cv-td-diff">{diff ?? "—"}</td>
                          )}
                        </tr>
                      );
                    })}
                  </>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}

      {/* Composite overlay image */}
      {compositeImage && (
        <div className="cv-results">
          <div className="cv-results-actions">
            <span style={{ color: "#60a5fa", fontWeight: 600, fontSize: 13 }}>
              Composite Overlay
            </span>
            <button
              className="cv-download-btn"
              onClick={() => {
                const link = document.createElement("a");
                link.href = `data:image/png;base64,${compositeImage}`;
                link.download = "composite_sounding_overlay.png";
                link.click();
              }}
              title="Download composite overlay"
            >
              <Download size={14} />
              Download Overlay
            </button>
          </div>
          <div className="cv-composite-plot">
            <img
              src={`data:image/png;base64,${compositeImage}`}
              alt="Composite sounding overlay"
              style={{ width: "100%", borderRadius: 8 }}
            />
          </div>
        </div>
      )}
    </div>
  );
}
