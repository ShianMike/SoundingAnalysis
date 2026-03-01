import { useState, useMemo } from "react";
import { Layers, X, ChevronUp, ChevronDown, AlertTriangle, TrendingUp, Search } from "lucide-react";
import "./MesoPanel.css";

const PARAMS = [
  { id: "stp",  label: "STP",    unit: "",       thresholds: [0.5, 1, 4],    desc: "Sig. Tornado" },
  { id: "scp",  label: "SCP",    unit: "",       thresholds: [1, 4, 10],     desc: "Supercell" },
  { id: "cape", label: "CAPE",   unit: "J/kg",   thresholds: [500, 1500, 3000], desc: "Instability" },
  { id: "srh",  label: "0-1 SRH", unit: "m²/s²", thresholds: [100, 200, 400], desc: "Helicity" },
  { id: "bwd",  label: "0-6 BWD", unit: "kt",    thresholds: [20, 35, 50],   desc: "Deep Shear" },
  { id: "ship", label: "SHIP",   unit: "",       thresholds: [0.5, 1, 2.5],  desc: "Sig. Hail" },
  { id: "dcp",  label: "DCP",    unit: "",       thresholds: [2, 4, 6],      desc: "Dmg. Wind" },
];

function chipColor(val, thresholds) {
  if (val == null) return "transparent";
  if (val >= thresholds[2]) return "#ef4444";
  if (val >= thresholds[1]) return "#f59e0b";
  if (val >= thresholds[0]) return "#22c55e";
  return "transparent";
}

function riskLevel(val, thresholds) {
  if (val == null) return "";
  if (val >= thresholds[2]) return "high";
  if (val >= thresholds[1]) return "mod";
  if (val >= thresholds[0]) return "low";
  return "";
}

/**
 * Mesoscale Analysis Dashboard — shows a sortable station table from risk scan data.
 */
export default function MesoPanel({ riskData, onStationSelect, onClose }) {
  const [activeParam, setActiveParam] = useState("stp");
  const [sortCol, setSortCol] = useState("stp");
  const [sortDir, setSortDir] = useState("desc");
  const [search, setSearch] = useState("");

  const paramDef = PARAMS.find((p) => p.id === activeParam) || PARAMS[0];

  const rows = useMemo(() => {
    if (!riskData?.stations) return [];
    return riskData.stations.map((s) => ({
      id: s.id,
      name: s.name || s.id,
      stp: s.stp ?? null,
      scp: s.scp ?? null,
      cape: s.cape ?? null,
      srh: s.srh ?? null,
      bwd: s.bwd ?? null,
      ship: s.ship ?? null,
      dcp: s.dcp ?? null,
    }));
  }, [riskData]);

  const sorted = useMemo(() => {
    let arr = [...rows];
    if (search) {
      const q = search.toLowerCase();
      arr = arr.filter((r) => r.id.toLowerCase().includes(q) || r.name.toLowerCase().includes(q));
    }
    arr.sort((a, b) => {
      if (sortCol === "id") {
        return sortDir === "desc" ? b.id.localeCompare(a.id) : a.id.localeCompare(b.id);
      }
      const va = a[sortCol] ?? -9999;
      const vb = b[sortCol] ?? -9999;
      return sortDir === "desc" ? vb - va : va - vb;
    });
    return arr;
  }, [rows, sortCol, sortDir, search]);

  const handleSort = (col) => {
    if (col === sortCol) {
      setSortDir((d) => (d === "desc" ? "asc" : "desc"));
    } else {
      setSortCol(col);
      setSortDir("desc");
    }
  };

  const SortIcon = ({ col }) => {
    if (col !== sortCol) return null;
    return sortDir === "desc" ? <ChevronDown size={10} /> : <ChevronUp size={10} />;
  };

  // Count stations with notable values for the active parameter
  const notableCount = rows.filter((r) => r[activeParam] != null && r[activeParam] >= paramDef.thresholds[0]).length;
  const highCount = rows.filter((r) => r[activeParam] != null && r[activeParam] >= paramDef.thresholds[2]).length;

  if (!riskData) {
    return (
      <div className="meso-panel">
        <div className="meso-header">
          <div className="meso-header-left">
            <Layers size={15} className="meso-header-icon" />
            <span className="meso-title">Mesoscale Dashboard</span>
          </div>
          <button className="meso-close" onClick={onClose}><X size={14} /></button>
        </div>
        <div className="meso-no-data">
          <AlertTriangle size={16} />
          <p>Run a <strong>Risk Scan</strong> first to populate the mesoscale dashboard.</p>
        </div>
      </div>
    );
  }

  const fmt = (v, digits = 1) => (v != null ? Number(v).toFixed(digits) : "—");

  return (
    <div className="meso-panel">
      {/* ── Header ── */}
      <div className="meso-header">
        <div className="meso-header-left">
          <Layers size={15} className="meso-header-icon" />
          <div>
            <span className="meso-title">Mesoscale Dashboard</span>
            <span className="meso-header-meta">
              {rows.length} stations &middot; {riskData.date || "latest"}
            </span>
          </div>
        </div>
        <button className="meso-close" onClick={onClose} title="Close"><X size={14} /></button>
      </div>

      {/* ── Parameter selector ── */}
      <div className="meso-param-row">
        {PARAMS.map((p) => (
          <button
            key={p.id}
            className={`meso-param-btn ${activeParam === p.id ? "active" : ""}`}
            onClick={() => { setActiveParam(p.id); setSortCol(p.id); setSortDir("desc"); }}
            title={p.desc}
          >
            {p.label}
          </button>
        ))}
      </div>

      {/* ── Summary badges ── */}
      <div className="meso-summary-bar">
        <span className="meso-summary-badge">{paramDef.label}: {paramDef.desc}</span>
        {highCount > 0 && (
          <span className="meso-summary-badge meso-badge-high">
            <TrendingUp size={10} /> {highCount} high-risk
          </span>
        )}
        {notableCount > 0 && (
          <span className="meso-summary-badge meso-badge-notable">
            {notableCount} notable
          </span>
        )}
        <div className="meso-search-wrap">
          <Search size={11} />
          <input
            type="text"
            className="meso-search"
            placeholder="Filter stations…"
            value={search}
            onChange={(e) => setSearch(e.target.value)}
          />
        </div>
      </div>

      {/* ── Station table ── */}
      <div className="meso-table-wrap">
        <table className="meso-table">
          <thead>
            <tr>
              <th className="meso-th-rank">#</th>
              <th onClick={() => handleSort("id")} className={`meso-th-stn ${sortCol === "id" ? "sorted" : ""}`}>
                Station <SortIcon col="id" />
              </th>
              <th className="meso-th-name">Name</th>
              {PARAMS.map((p) => (
                <th
                  key={p.id}
                  onClick={() => handleSort(p.id)}
                  className={`${sortCol === p.id ? "sorted" : ""} ${activeParam === p.id ? "meso-th-active" : ""}`}
                >
                  {p.label} <SortIcon col={p.id} />
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {sorted.map((r, idx) => {
              const highlight = r[activeParam] != null &&
                r[activeParam] >= paramDef.thresholds[2];
              return (
                <tr key={r.id} className={highlight ? "meso-row-high" : ""}>
                  <td className="meso-rank-cell">{idx + 1}</td>
                  <td
                    className="meso-stn-cell"
                    onClick={() => onStationSelect?.(r.id)}
                  >
                    {r.id}
                  </td>
                  <td className="meso-name-cell">{r.name}</td>
                  {PARAMS.map((p) => {
                    const level = riskLevel(r[p.id], p.thresholds);
                    return (
                      <td key={p.id} className={`meso-val-cell ${level ? `meso-val-${level}` : ""}`}>
                        {chipColor(r[p.id], p.thresholds) !== "transparent" && (
                          <span
                            className="meso-chip"
                            style={{ background: chipColor(r[p.id], p.thresholds) }}
                          />
                        )}
                        {fmt(r[p.id], p.id === "cape" ? 0 : (p.id === "bwd" || p.id === "srh" ? 0 : 2))}
                      </td>
                    );
                  })}
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>

      <div className="meso-info">
        Showing {sorted.length} of {rows.length} stations &middot; Click station ID to load sounding
      </div>
    </div>
  );
}
