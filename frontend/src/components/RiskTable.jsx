import { Zap, Download } from "lucide-react";
import { exportRiskTablePng } from "../utils/exportImage";

/**
 * Severe-Weather Risk Scan summary table.
 *
 * Renders a sortable list of stations ranked by severe-weather composite
 * indices (STP / SCP / SHIP / DCP / CAPE / SRH / BWD). Clicking a row
 * dispatches `onStationSelect(stationId, riskData)` so the caller can load
 * that station's full sounding.
 *
 * Extracted from `ResultsView.jsx` as part of the Phase-1 refactor; behavior
 * and styling are unchanged.
 */
export default function RiskTable({ riskData, onStationSelect }) {
  if (!riskData || !riskData.stations || riskData.stations.length === 0) return null;

  const handleExportPng = () => exportRiskTablePng(riskData);

  return (
    <div className="rv-risk-table-wrap">
      <div className="rv-risk-table-header">
        <Zap size={14} />
        <h3>Severe Weather Risk Scan</h3>
        {riskData.model && (
          <span className="rv-risk-table-model">{riskData.model} F{riskData.fhour}</span>
        )}
        <span className="rv-risk-table-date">{riskData.date}</span>
        <span className="rv-risk-table-count">{riskData.stations.length} stations</span>
        <button
          type="button"
          className="rv-risk-export-btn"
          onClick={handleExportPng}
          title="Export as PNG"
        >
          <Download size={13} />
        </button>
      </div>
      <div className="rv-risk-table-scroll">
        <table className="rv-risk-table">
          <thead>
            <tr>
              <th>#</th>
              <th>Station</th>
              <th>Name</th>
              <th className="rv-rt-num">STP</th>
              <th className="rv-rt-num">SCP</th>
              <th className="rv-rt-num">SHIP</th>
              <th className="rv-rt-num">DCP</th>
              <th className="rv-rt-num">CAPE</th>
              <th className="rv-rt-num">SRH</th>
              <th className="rv-rt-num">BWD</th>
            </tr>
          </thead>
          <tbody>
            {riskData.stations.map((s, i) => (
              <tr
                key={s.id}
                className={`${s.stp >= 1 ? "rv-rt-high" : s.stp >= 0.3 ? "rv-rt-med" : ""} rv-rt-clickable`}
                onClick={() => onStationSelect?.(s.id, riskData)}
                title={`Load ${s.id} sounding${riskData.model ? ` (${riskData.model} F${riskData.fhour})` : ""}`}
              >
                <td className="rv-rt-rank">{i + 1}</td>
                <td className="rv-rt-id">{s.id}</td>
                <td className="rv-rt-name">{s.name}</td>
                <td className="rv-rt-num">
                  <span className={`rv-rt-stp ${s.stp >= 1 ? "high" : s.stp >= 0.3 ? "med" : "low"}`}>
                    {s.stp.toFixed(2)}
                  </span>
                </td>
                <td className="rv-rt-num">
                  <span className={`rv-rt-stp ${s.scp >= 4 ? "high" : s.scp >= 1 ? "med" : "low"}`}>
                    {s.scp.toFixed(2)}
                  </span>
                </td>
                <td className="rv-rt-num">
                  <span className={`rv-rt-stp ${s.ship >= 1.5 ? "high" : s.ship >= 0.5 ? "med" : "low"}`}>
                    {s.ship.toFixed(2)}
                  </span>
                </td>
                <td className="rv-rt-num">
                  <span className={`rv-rt-stp ${s.dcp >= 4 ? "high" : s.dcp >= 2 ? "med" : "low"}`}>
                    {s.dcp.toFixed(2)}
                  </span>
                </td>
                <td className="rv-rt-num">{s.cape}</td>
                <td className="rv-rt-num">{s.srh}</td>
                <td className="rv-rt-num">{s.bwd}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}
