import { useState, useMemo } from "react";
import { Layers, X, ChevronUp, ChevronDown } from "lucide-react";
import "./MesoPanel.css";

const PARAMS = [
  { id: "stp",  label: "STP",    unit: "",     thresholds: [0.5, 1, 4] },
  { id: "scp",  label: "SCP",    unit: "",     thresholds: [1, 4, 10] },
  { id: "cape", label: "CAPE",   unit: "J/kg", thresholds: [500, 1500, 3000] },
  { id: "srh",  label: "0-1 SRH", unit: "m²/s²", thresholds: [100, 200, 400] },
  { id: "shear",label: "0-6 BWD", unit: "kt",  thresholds: [20, 35, 50] },
  { id: "ship", label: "SHIP",   unit: "",     thresholds: [0.5, 1, 2.5] },
  { id: "dcape",label: "DCAPE",  unit: "J/kg", thresholds: [500, 1000, 1500] },
  { id: "pwat", label: "PWAT",   unit: "mm",   thresholds: [20, 35, 50] },
];

function chipColor(val, thresholds) {
  if (val == null) return "#555";
  if (val >= thresholds[2]) return "#ef4444"; // high
  if (val >= thresholds[1]) return "#f59e0b"; // moderate
  if (val >= thresholds[0]) return "#22c55e"; // notable
  return "#6b7280"; // low
}

/**
 * Mesoscale Analysis Dashboard — shows a sortable station table from risk scan data.
 */
export default function MesoPanel({ riskData, onStationSelect, onClose }) {
  const [activeParam, setActiveParam] = useState("stp");
  const [sortCol, setSortCol] = useState("stp");
  const [sortDir, setSortDir] = useState("desc");

  const paramDef = PARAMS.find((p) => p.id === activeParam) || PARAMS[0];

  const rows = useMemo(() => {
    if (!riskData?.stations) return [];
    return riskData.stations.map((s) => ({
      id: s.station,
      stp: s.stp ?? null,
      scp: s.scp ?? null,
      cape: s.cape ?? null,
      srh: s.srh ?? null,
      shear: s.shear ?? null,
      ship: s.ship ?? null,
      dcape: s.dcape ?? null,
      pwat: s.pwat ?? null,
    }));
  }, [riskData]);

  const sorted = useMemo(() => {
    const arr = [...rows];
    arr.sort((a, b) => {
      const va = a[sortCol] ?? -9999;
      const vb = b[sortCol] ?? -9999;
      return sortDir === "desc" ? vb - va : va - vb;
    });
    return arr;
  }, [rows, sortCol, sortDir]);

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

  if (!riskData) {
    return (
      <div className="meso-panel">
        <div className="meso-header">
          <span className="meso-title"><Layers size={14} /> Mesoscale Dashboard</span>
          <button className="meso-close" onClick={onClose}><X size={14} /></button>
        </div>
        <div className="meso-no-data">
          Run a Risk Scan first to populate the mesoscale dashboard.
        </div>
      </div>
    );
  }

  const fmt = (v, digits = 1) => (v != null ? Number(v).toFixed(digits) : "—");

  return (
    <div className="meso-panel">
      <div className="meso-header">
        <span className="meso-title"><Layers size={14} /> Mesoscale Dashboard</span>
        <button className="meso-close" onClick={onClose}><X size={14} /></button>
      </div>

      {/* Parameter highlight selector */}
      <div className="meso-param-row">
        {PARAMS.map((p) => (
          <button
            key={p.id}
            className={`meso-param-btn ${activeParam === p.id ? "active" : ""}`}
            onClick={() => { setActiveParam(p.id); setSortCol(p.id); setSortDir("desc"); }}
          >
            {p.label}
          </button>
        ))}
      </div>

      {/* Station table */}
      <div className="meso-table-wrap">
        <table className="meso-table">
          <thead>
            <tr>
              <th onClick={() => handleSort("id")} className={sortCol === "id" ? "sorted" : ""}>
                Stn <SortIcon col="id" />
              </th>
              {PARAMS.map((p) => (
                <th
                  key={p.id}
                  onClick={() => handleSort(p.id)}
                  className={sortCol === p.id ? "sorted" : ""}
                >
                  {p.label} <SortIcon col={p.id} />
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {sorted.map((r) => {
              const highlight = r[activeParam] != null &&
                r[activeParam] >= paramDef.thresholds[2];
              return (
                <tr key={r.id} className={highlight ? "meso-row-high" : ""}>
                  <td
                    className="meso-stn-cell"
                    onClick={() => onStationSelect?.(r.id)}
                  >
                    {r.id}
                  </td>
                  {PARAMS.map((p) => (
                    <td key={p.id} className="meso-val-cell">
                      <span
                        className="meso-chip"
                        style={{ background: chipColor(r[p.id], p.thresholds) }}
                      />
                      {fmt(r[p.id], p.id === "cape" || p.id === "dcape" ? 0 : 1)}
                    </td>
                  ))}
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>

      <div className="meso-info">
        {sorted.length} stations &bull; Scan: {riskData.scanDate || "latest"} &bull;
        Click station to select
      </div>
    </div>
  );
}
