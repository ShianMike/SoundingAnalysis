import { useState, useMemo } from "react";
import { X, ChevronDown, ChevronUp, AlertTriangle, Zap, Wind, CloudHail } from "lucide-react";
import "./OutlookPanel.css";

/* ── Colour & label maps ─────────────────────────────── */
const CAT_META = {
  NONE:  { color: "#666",    bg: "rgba(102,102,102,0.10)", label: "None" },
  TSTM:  { color: "#78b878", bg: "rgba(120,184,120,0.12)", label: "Thunder" },
  MRGL:  { color: "#66a366", bg: "rgba(102,163,102,0.14)", label: "Marginal" },
  SLGT:  { color: "#ffe066", bg: "rgba(255,224,102,0.14)", label: "Slight" },
  ENH:   { color: "#ffa500", bg: "rgba(255,165,0,0.14)",   label: "Enhanced" },
  MDT:   { color: "#ff4444", bg: "rgba(255,68,68,0.14)",   label: "Moderate" },
  HIGH:  { color: "#ff00ff", bg: "rgba(255,0,255,0.14)",   label: "High" },
};

const HAZ_META = {
  NONE:    { color: "#666",    label: "None" },
  LOW:     { color: "#66a366", label: "Low" },
  MOD:     { color: "#ffa500", label: "Moderate" },
  HIGH:    { color: "#ff4444", label: "High" },
  EXTREME: { color: "#ff00ff", label: "Extreme" },
};

function Badge({ level, meta }) {
  const m = meta[level] || meta.NONE;
  return (
    <span className="olp-badge" style={{ color: m.color, borderColor: m.color }}>
      {m.label}
    </span>
  );
}

export default function OutlookPanel({ data, onClose, onStationSelect }) {
  const [expanded, setExpanded] = useState(true);
  const [filter, setFilter] = useState("ALL");

  const summary = data?.summary;
  const stations = useMemo(() => data?.stations || [], [data?.stations]);

  const filtered = useMemo(() => {
    if (filter === "ALL") return stations.filter((s) => s.categorical !== "NONE");
    return stations.filter((s) => s.categorical === filter);
  }, [stations, filter]);

  if (!data) return null;

  const catM = CAT_META[summary?.categorical] || CAT_META.NONE;

  return (
    <div className="olp-panel">
      {/* Header */}
      <div className="olp-header" style={{ borderLeftColor: catM.color }}>
        <div className="olp-header-left">
          <AlertTriangle size={16} style={{ color: catM.color }} />
          <h3 className="olp-title">Custom Severe Weather Outlook</h3>
        </div>
        <div className="olp-header-right">
          <span className="olp-valid">
            {data.model && `${data.model} F${data.fhour} · `}
            Valid {data.validTime}
          </span>
          <button className="olp-close" onClick={onClose}><X size={14} /></button>
        </div>
      </div>

      {/* Summary cards */}
      <div className="olp-summary">
        <div className="olp-card olp-card-cat" style={{ background: catM.bg, borderColor: catM.color }}>
          <span className="olp-card-label">Categorical</span>
          <span className="olp-card-value" style={{ color: catM.color }}>
            {catM.label}
          </span>
        </div>
        <div className="olp-card">
          <Zap size={13} className="olp-card-icon" />
          <span className="olp-card-label">Tornado</span>
          <Badge level={summary?.tornado} meta={HAZ_META} />
        </div>
        <div className="olp-card">
          <Wind size={13} className="olp-card-icon" />
          <span className="olp-card-label">Wind</span>
          <Badge level={summary?.wind} meta={HAZ_META} />
        </div>
        <div className="olp-card">
          <CloudHail size={13} className="olp-card-icon" />
          <span className="olp-card-label">Hail</span>
          <Badge level={summary?.hail} meta={HAZ_META} />
        </div>
      </div>

      {/* Category counts bar */}
      <div className="olp-counts">
        <span className="olp-counts-label">{summary?.total} stations scanned</span>
        <div className="olp-count-chips">
          {["HIGH", "MDT", "ENH", "SLGT", "MRGL", "TSTM"].map((cat) => {
            const cnt = summary?.counts?.[cat] || 0;
            if (cnt === 0) return null;
            const m = CAT_META[cat];
            return (
              <button
                key={cat}
                className={`olp-chip${filter === cat ? " active" : ""}`}
                style={{ color: m.color, borderColor: m.color }}
                onClick={() => setFilter(filter === cat ? "ALL" : cat)}
              >
                {cat} ({cnt})
              </button>
            );
          })}
        </div>
      </div>

      {/* Station table */}
      <button
        className="olp-toggle"
        onClick={() => setExpanded((v) => !v)}
      >
        {expanded ? <ChevronUp size={12} /> : <ChevronDown size={12} />}
        {expanded ? "Hide stations" : `Show stations (${filtered.length})`}
      </button>

      {expanded && filtered.length > 0 && (
        <div className="olp-table-wrap">
          <table className="olp-table">
            <thead>
              <tr>
                <th>Station</th>
                <th>Cat</th>
                <th>Tor</th>
                <th>Wind</th>
                <th>Hail</th>
                <th>STP</th>
                <th>SCP</th>
                <th>SHIP</th>
                <th>CAPE</th>
                <th>SRH</th>
                <th>BWD</th>
              </tr>
            </thead>
            <tbody>
              {filtered.map((s) => {
                const cm = CAT_META[s.categorical] || CAT_META.NONE;
                return (
                  <tr
                    key={s.id}
                    className="olp-row"
                    onClick={() => onStationSelect?.(s.id)}
                    style={{ cursor: onStationSelect ? "pointer" : undefined }}
                  >
                    <td className="olp-cell-station">
                      <span className="olp-dot" style={{ background: cm.color }} />
                      <span className="olp-sid">{s.id}</span>
                      <span className="olp-sname">{s.name}</span>
                    </td>
                    <td><span className="olp-cat-tag" style={{ color: cm.color, borderColor: cm.color }}>{s.categorical}</span></td>
                    <td><Badge level={s.tornado} meta={HAZ_META} /></td>
                    <td><Badge level={s.wind} meta={HAZ_META} /></td>
                    <td><Badge level={s.hail} meta={HAZ_META} /></td>
                    <td>{s.stp}</td>
                    <td>{s.scp}</td>
                    <td>{s.ship}</td>
                    <td>{s.cape}</td>
                    <td>{s.srh}</td>
                    <td>{s.bwd}</td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
      )}

      {expanded && filtered.length === 0 && (
        <p className="olp-empty">No stations meet the selected filter criteria.</p>
      )}
    </div>
  );
}
