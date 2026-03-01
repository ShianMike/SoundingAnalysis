import { useState, useRef, useCallback, useEffect } from "react";
import {
  ImageIcon,
  BarChart3,
  Info,
  Loader2,
  AlertTriangle,
  Download,
  Maximize2,
  ZoomIn,
  ZoomOut,
  Wind,
  Thermometer,
  ArrowUpDown,
  Droplets,
  Gauge,
  Zap,
  FileSpreadsheet,
  Link2,
  Check,
  FileText,
  ChevronDown,
} from "lucide-react";
import StationMap from "./StationMap";
import TimeSeriesChart from "./TimeSeriesChart";
import ComparisonView from "./ComparisonView";
import VwpDisplay from "./VwpDisplay";
import "./ResultsView.css";

/* ── Parameter severity thresholds ─────────────────────────── */
const THRESHOLDS = {
  "SB CAPE":  [{ v: 3000, c: "extreme" }, { v: 1500, c: "high" }, { v: 500, c: "mod" }],
  "MU CAPE":  [{ v: 3000, c: "extreme" }, { v: 1500, c: "high" }, { v: 500, c: "mod" }],
  "ML CAPE":  [{ v: 3000, c: "extreme" }, { v: 1500, c: "high" }, { v: 500, c: "mod" }],
  "ECAPE":    [{ v: 2000, c: "extreme" }, { v: 1000, c: "high" }, { v: 400, c: "mod" }],
  "DCAPE":    [{ v: 1500, c: "extreme" }, { v: 800, c: "high" }],
  "3CAPE":    [{ v: 100, c: "extreme" }, { v: 50, c: "high" }],
  "STP":      [{ v: 4, c: "extreme" }, { v: 1, c: "high" }, { v: 0.5, c: "mod" }],
  "SCP":      [{ v: 8, c: "extreme" }, { v: 4, c: "high" }, { v: 1, c: "mod" }],
  "SHIP":     [{ v: 3, c: "extreme" }, { v: 1.5, c: "high" }, { v: 0.5, c: "mod" }],
  "DCP":      [{ v: 6, c: "extreme" }, { v: 4, c: "high" }, { v: 2, c: "mod" }],
  "LR 0-3 km":[{ v: 9, c: "extreme" }, { v: 8, c: "high" }, { v: 7, c: "mod" }],
  "LR 3-6 km":[{ v: 9, c: "extreme" }, { v: 8, c: "high" }, { v: 7, c: "mod" }],
  "BWD 0-6 km":[{ v: 60, c: "extreme" }, { v: 40, c: "high" }, { v: 25, c: "mod" }],
  "BWD 0-1 km":[{ v: 30, c: "extreme" }, { v: 20, c: "high" }, { v: 15, c: "mod" }],
  "SRH 0-1 km":[{ v: 300, c: "extreme" }, { v: 150, c: "high" }, { v: 100, c: "mod" }],
  "SRH 0-3 km":[{ v: 400, c: "extreme" }, { v: 200, c: "high" }, { v: 100, c: "mod" }],
  "Eff. SRH": [{ v: 300, c: "extreme" }, { v: 150, c: "high" }, { v: 100, c: "mod" }],
};

function getAlertClass(label, value) {
  if (value == null || typeof value !== "number") return "";
  const rules = THRESHOLDS[label];
  if (!rules) return "";
  for (const r of rules) {
    if (value >= r.v) return `param-alert-${r.c}`;
  }
  return "";
}

/* ── Sounding text summary generator ───────────────────────── */
function generateSoundingSummary(params, meta) {
  const lines = [];
  const n = (v) => (v != null && typeof v === "number" ? v : null);

  // ── Instability ──
  const sbCape = n(params.sbCape);
  const mlCape = n(params.mlCape);
  const muCape = n(params.muCape);
  const ecape  = n(params.ecape);
  const bestCape = Math.max(sbCape ?? 0, mlCape ?? 0, muCape ?? 0);

  if (bestCape < 100) {
    lines.push("The atmosphere is currently stable with negligible buoyancy. No significant thunderstorm development is expected from this profile alone.");
  } else if (bestCape < 500) {
    lines.push(`Marginal instability is present (MUCAPE ${muCape ?? "--"} J/kg). Weak updrafts are possible if a lifting mechanism is available, but storm intensity would be limited.`);
  } else if (bestCape < 1500) {
    lines.push(`Moderate instability exists with MUCAPE of ${muCape ?? "--"} J/kg${ecape != null ? ` (ECAPE ${ecape} J/kg)` : ""}. This supports organized thunderstorm updrafts if convective initiation occurs.`);
  } else if (bestCape < 3000) {
    lines.push(`Large instability is present — MUCAPE ${muCape ?? "--"} J/kg${ecape != null ? `, ECAPE ${ecape} J/kg` : ""}. This provides substantial energy for strong to severe updrafts.`);
  } else {
    lines.push(`Extreme instability is present with MUCAPE of ${muCape ?? "--"} J/kg${ecape != null ? ` (ECAPE ${ecape} J/kg)` : ""}. This is a high-end thermodynamic environment supporting violent updrafts and very large hail potential.`);
  }

  // ── CIN / Cap ──
  const mlCin = n(params.mlCin);
  if (mlCin != null) {
    if (mlCin > -25) {
      lines.push("The convective cap is weak or absent — storms could initiate easily with minimal forcing.");
    } else if (mlCin > -100) {
      lines.push(`A moderate cap is in place (MLCIN ${mlCin} J/kg). Mesoscale forcing (boundaries, terrain, or outflow) may be needed to initiate convection.`);
    } else {
      lines.push(`A strong cap is present (MLCIN ${mlCin} J/kg). Convective initiation is unlikely without strong synoptic-scale or mesoscale forcing, but storms that do break through may be explosive.`);
    }
  }

  // ── LCL height ──
  const mlLcl = n(params.mlLclM);
  if (mlLcl != null) {
    if (mlLcl < 1000) {
      lines.push(`Cloud bases are very low (MLLCL ${mlLcl} m AGL), which is favorable for tornado development by reducing the sub-cloud dry layer that would otherwise disrupt low-level rotation.`);
    } else if (mlLcl < 1500) {
      lines.push(`Cloud bases are relatively low (MLLCL ${mlLcl} m AGL), marginally supportive of tornado potential.`);
    } else {
      lines.push(`Cloud bases are elevated (MLLCL ${mlLcl} m AGL). This reduces tornado probability but may enhance downdraft wind potential and dry-slot entrainment.`);
    }
  }

  // ── Lapse rates ──
  const lr03 = n(params.lr03);
  const lr36 = n(params.lr36);
  if (lr03 != null || lr36 != null) {
    const parts = [];
    if (lr03 != null) {
      if (lr03 >= 8) parts.push(`very steep 0–3 km lapse rates (${lr03} °C/km) indicating near-adiabatic boundary layer conditions`);
      else if (lr03 >= 7) parts.push(`moderately steep 0–3 km lapse rates (${lr03} °C/km)`);
    }
    if (lr36 != null) {
      if (lr36 >= 8) parts.push(`extreme mid-level lapse rates (${lr36} °C/km) enhancing CAPE depth`);
      else if (lr36 >= 7) parts.push(`steep mid-level lapse rates (${lr36} °C/km)`);
    }
    if (parts.length > 0) {
      lines.push(`Lapse rates show ${parts.join(" and ")}.`);
    }
  }

  // ── Deep shear ──
  const bwd6 = n(params.bwd6km);
  const bwd1 = n(params.bwd1km);
  if (bwd6 != null) {
    if (bwd6 >= 60) {
      lines.push(`Deep-layer shear is extreme (0–6 km BWD ${bwd6} kt), strongly favoring discrete supercells with long-lived mesocyclones.`);
    } else if (bwd6 >= 40) {
      lines.push(`Strong deep-layer shear is present (0–6 km BWD ${bwd6} kt), largely supporting supercellular convective mode.`);
    } else if (bwd6 >= 25) {
      lines.push(`Moderate deep-layer shear (0–6 km BWD ${bwd6} kt) supports organized multicell to supercell storms.`);
    } else {
      lines.push(`Weak deep-layer shear (0–6 km BWD ${bwd6} kt) favors pulse or weakly organized multicellular storms.`);
    }
  }

  // ── Low-level shear & SRH ──
  const srh1 = n(params.srh1km);
  const srh3 = n(params.srh3km);
  if (srh1 != null) {
    if (srh1 >= 300) {
      lines.push(`Low-level storm-relative helicity is extreme (0–1 km SRH ${srh1} m²/s²), highly favorable for violent tornadoes.`);
    } else if (srh1 >= 150) {
      lines.push(`Significant low-level helicity (0–1 km SRH ${srh1} m²/s²) supports strong mesocyclones and significant tornado potential.`);
    } else if (srh1 >= 100) {
      lines.push(`Moderate low-level helicity (0–1 km SRH ${srh1} m²/s²) supports tornado development in supercells.`);
    }
  }

  // ── Composite indices ──
  const stp = n(params.stp);
  const scp = n(params.scp);
  const ship = n(params.ship);
  const dcp = n(params.dcp);
  const composites = [];

  if (stp != null && stp >= 1) {
    composites.push(stp >= 4
      ? `STP of ${stp.toFixed(1)} signals a **high-end significant tornado environment**`
      : `STP of ${stp.toFixed(1)} indicates a favorable environment for significant (EF2+) tornadoes`);
  }
  if (scp != null && scp >= 1) {
    composites.push(scp >= 8
      ? `SCP of ${scp.toFixed(1)} strongly favors long-lived discrete supercells`
      : `SCP of ${scp.toFixed(1)} supports supercell development`);
  }
  if (ship != null && ship >= 0.5) {
    composites.push(ship >= 1.5
      ? `SHIP of ${ship.toFixed(1)} suggests a significant hail (≥2 in.) environment`
      : `SHIP of ${ship.toFixed(1)} indicates some potential for large hail`);
  }
  if (dcp != null && dcp >= 2) {
    composites.push(dcp >= 6
      ? `DCP of ${dcp.toFixed(1)} signals an extreme derecho/long-lived wind event risk`
      : `DCP of ${dcp.toFixed(1)} suggests potential for organized damaging wind events`);
  }

  if (composites.length > 0) {
    lines.push("Composite parameters: " + composites.join(". ") + ".");
  } else if (bestCape >= 500) {
    lines.push("Composite severe parameters (STP, SCP, SHIP) are below significant thresholds, suggesting isolated or sub-severe convection is most likely.");
  }

  // ── DCAPE / downburst ──
  const dcape = n(params.dcape);
  if (dcape != null && dcape >= 800) {
    lines.push(`DCAPE of ${dcape} J/kg${dcape >= 1500 ? " is extreme and" : ""} supports strong to damaging outflow winds at the surface.`);
  }

  // ── PWAT ──
  const pwat = n(params.pwat);
  if (pwat != null && pwat >= 40) {
    lines.push(`Precipitable water of ${pwat} mm is anomalously high, elevating the risk of flash flooding from training or slow-moving convection.`);
  }

  // ── Overall threat assessment ──
  const threats = [];
  if (stp != null && stp >= 1) threats.push("tornadoes");
  if (ship != null && ship >= 1) threats.push("significant hail");
  else if (bestCape >= 2000 && bwd6 >= 30) threats.push("large hail");
  if (dcp != null && dcp >= 2) threats.push("damaging straight-line winds");
  else if (dcape != null && dcape >= 800 && bwd6 >= 25) threats.push("damaging outflow winds");
  if (pwat != null && pwat >= 40) threats.push("flash flooding");

  if (threats.length > 0) {
    lines.push(`Primary threats from this profile: ${threats.join(", ")}.`);
  } else if (bestCape >= 500) {
    lines.push("Overall, the environment supports convection but severe weather parameters remain below high-confidence thresholds for significant threats.");
  }

  return lines.join("\n\n");
}

function ParamCard({ label, value, unit, color, desc }) {
  const alertCls = getAlertClass(label, value);
  return (
    <div className={`param-card ${alertCls}`} title="">
      <span className="param-label">{label}</span>
      <span className="param-value" style={color ? { color } : {}}>
        {value ?? "---"}
      </span>
      {unit && value != null && <span className="param-unit">{unit}</span>}
      {desc && <span className="param-tooltip">{desc}</span>}
    </div>
  );
}

function ParamSection({ title, icon, children }) {
  return (
    <div className="param-section">
      <div className="param-section-header">
        {icon}
        <h3>{title}</h3>
      </div>
      <div className="param-grid">{children}</div>
    </div>
  );
}

function RiskTable({ riskData }) {
  if (!riskData || !riskData.stations || riskData.stations.length === 0) return null;

  return (
    <div className="rv-risk-table-wrap">
      <div className="rv-risk-table-header">
        <Zap size={14} />
        <h3>Severe Weather Risk Scan</h3>
        <span className="rv-risk-table-date">{riskData.date}</span>
        <span className="rv-risk-table-count">{riskData.stations.length} stations</span>
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
              <tr key={s.id} className={s.stp >= 1 ? "rv-rt-high" : s.stp >= 0.3 ? "rv-rt-med" : ""}>
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

export default function ResultsView({ result, loading, error, riskData, showRisk, showMap, mapProps, showTimeSeries, onCloseTimeSeries, showCompare, onCloseCompare, showVwp, onCloseVwp, compareHistoryData, onCompareHistoryConsumed, stations, selectedStation, source, lastParams }) {
  if (error) {
    return (
      <div className="results-view">
        {showMap && mapProps && <StationMap {...mapProps} />}
        {showRisk && <RiskTable riskData={riskData} />}
        {showTimeSeries && (
          <TimeSeriesChart station={selectedStation} source={source} onClose={onCloseTimeSeries} />
        )}
        {showCompare && (
          <ComparisonView stations={stations || []} onClose={onCloseCompare} historyData={compareHistoryData} onHistoryConsumed={onCompareHistoryConsumed} />
        )}
        {showVwp && (
          <VwpDisplay stations={stations || []} selectedStation={selectedStation} onClose={onCloseVwp} />
        )}
        <div className="rv-state rv-error">
          <AlertTriangle size={24} />
          <div>
            <h3>Analysis Failed</h3>
            <p>{error}</p>
          </div>
        </div>
      </div>
    );
  }

  if (loading) {
    return (
      <div className="results-view">
        {showMap && mapProps && <StationMap {...mapProps} />}
        {showRisk && <RiskTable riskData={riskData} />}
        {showTimeSeries && (
          <TimeSeriesChart station={selectedStation} source={source} onClose={onCloseTimeSeries} />
        )}
        {showCompare && (
          <ComparisonView stations={stations || []} onClose={onCloseCompare} historyData={compareHistoryData} onHistoryConsumed={onCompareHistoryConsumed} />
        )}
        {showVwp && (
          <VwpDisplay stations={stations || []} selectedStation={selectedStation} onClose={onCloseVwp} />
        )}
        <div className="rv-state rv-loading">
          <Loader2 size={24} className="spin" />
          <div>
            <h3>Fetching & Analyzing</h3>
            <p>
              Retrieving sounding data, computing thermodynamic and kinematic
              parameters, and generating the analysis plot...
            </p>
          </div>
        </div>
      </div>
    );
  }

  if (!result) {
    return (
      <div className="results-view">
        {showMap && mapProps && <StationMap {...mapProps} />}
        {showRisk && <RiskTable riskData={riskData} />}
        {showTimeSeries && (
          <TimeSeriesChart station={selectedStation} source={source} onClose={onCloseTimeSeries} />
        )}
        {showCompare && (
          <ComparisonView stations={stations || []} onClose={onCloseCompare} historyData={compareHistoryData} onHistoryConsumed={onCompareHistoryConsumed} />
        )}
        {showVwp && (
          <VwpDisplay stations={stations || []} selectedStation={selectedStation} onClose={onCloseVwp} />
        )}
        {!riskData && (
          <div className="rv-state rv-empty">
            <Wind size={32} />
            <div>
              <h3>No Sounding Loaded</h3>
              <p>
                Select a data source and station, then click Generate Sounding
                to fetch and analyze upper-air data.
              </p>
            </div>
          </div>
        )}
      </div>
    );
  }

  const { image, params, meta } = result;
  const [zoomed, setZoomed] = useState(false);
  const [linkCopied, setLinkCopied] = useState(false);
  const [exportOpen, setExportOpen] = useState(false);
  const [showSummary, setShowSummary] = useState(false);
  const exportRef = useRef(null);
  const plotRef = useRef(null);
  const dragRef = useRef({ dragging: false, startX: 0, startY: 0, scrollLeft: 0, scrollTop: 0 });

  // Close export menu on outside click
  useEffect(() => {
    if (!exportOpen) return;
    const handler = (e) => {
      if (exportRef.current && !exportRef.current.contains(e.target)) setExportOpen(false);
    };
    document.addEventListener("mousedown", handler);
    return () => document.removeEventListener("mousedown", handler);
  }, [exportOpen]);

  const handleMouseDown = useCallback((e) => {
    if (!zoomed) return;
    const el = plotRef.current;
    if (!el) return;
    dragRef.current = {
      dragging: true,
      startX: e.clientX,
      startY: e.clientY,
      scrollLeft: el.scrollLeft,
      scrollTop: el.scrollTop,
    };
    el.style.cursor = "grabbing";
    e.preventDefault();
  }, [zoomed]);

  const handleMouseMove = useCallback((e) => {
    const d = dragRef.current;
    if (!d.dragging) return;
    const el = plotRef.current;
    if (!el) return;
    el.scrollLeft = d.scrollLeft - (e.clientX - d.startX);
    el.scrollTop = d.scrollTop - (e.clientY - d.startY);
  }, []);

  const handleMouseUp = useCallback(() => {
    dragRef.current.dragging = false;
    const el = plotRef.current;
    if (el && zoomed) el.style.cursor = "grab";
  }, [zoomed]);

  const handleDownload = () => {
    const link = document.createElement("a");
    link.href = `data:image/png;base64,${image}`;
    link.download = `sounding_${meta.station || "analysis"}_${meta.date.replace(/\s/g, "_")}.png`;
    link.click();
  };

  const handleFullscreen = () => {
    const w = window.open();
    w.document.write(
      `<html><head><title>Sounding - ${meta.station}</title>
       <style>body{margin:0;background:#0a0a0a;display:flex;align-items:center;justify-content:center;min-height:100vh}
       img{max-width:100%;height:auto}</style></head>
       <body><img src="data:image/png;base64,${image}" /></body></html>`
    );
  };

  const handleCsvExport = () => {
    const rows = [["Parameter", "Value", "Unit"]];
    const entries = [
      ["SB CAPE", params.sbCape, "J/kg"],
      ["SB CIN", params.sbCin, "J/kg"],
      ["SB LCL", params.sbLclM, "m AGL"],
      ["MU CAPE", params.muCape, "J/kg"],
      ["MU CIN", params.muCin, "J/kg"],
      ["MU LCL", params.muLclM, "m AGL"],
      ["ML CAPE", params.mlCape, "J/kg"],
      ["ML CIN", params.mlCin, "J/kg"],
      ["ML LCL", params.mlLclM, "m AGL"],
      ["DCAPE", params.dcape, "J/kg"],
      ["ECAPE", params.ecape, "J/kg"],
      ["STP", params.stp, ""],
      ["SCP", params.scp, ""],
      ["SHIP", params.ship, ""],
      ["DCP", params.dcp, ""],
      ["LR 0-3 km", params.lr03, "C/km"],
      ["LR 3-6 km", params.lr36, "C/km"],
      ["PWAT", params.pwat, "mm"],
      ["FRZ Level", params.frzLevel, "m AGL"],
      ["WB Zero", params.wbo, "m AGL"],
      ["RH 0-1 km", params.rh01, "%"],
      ["RH 1-3 km", params.rh13, "%"],
      ["RH 3-6 km", params.rh36, "%"],
      ["BWD 0-500m", params.bwd500m, "kt"],
      ["BWD 0-1 km", params.bwd1km, "kt"],
      ["BWD 0-3 km", params.bwd3km, "kt"],
      ["BWD 0-6 km", params.bwd6km, "kt"],
      ["SRH 500m", params.srh500m, "m2/s2"],
      ["SRH 0-1 km", params.srh1km, "m2/s2"],
      ["SRH 0-3 km", params.srh3km, "m2/s2"],
    ];
    entries.forEach(([name, val, unit]) => {
      rows.push([name, val != null ? String(val) : "", unit]);
    });

    const csv = rows.map((r) => r.map((c) => `"${c}"`).join(",")).join("\n");
    const blob = new Blob([csv], { type: "text/csv" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = `sounding_${meta.station || "analysis"}_${meta.date.replace(/\s/g, "_")}.csv`;
    link.click();
    URL.revokeObjectURL(url);
  };

  // ── SHARPpy format export ──────────────────────────────────
  const handleSharppy = () => {
    const profile = result.profile;
    if (!profile || profile.length === 0) return;
    const lines = [];
    lines.push(`%TITLE%`);
    lines.push(` ${meta.station || "XXXX"}   ${meta.date.replace(/\s/g, "").replace("-", "").replace("Z", "")}`);
    lines.push(``);
    lines.push(`   LEVEL       HGHT       TEMP       DWPT       WDIR       WSPD`);
    lines.push(`-------------------------------------------------------------------`);
    lines.push(`%RAW%`);
    for (const lv of profile) {
      const p = lv.p != null ? lv.p.toFixed(2) : "9999.00";
      const h = lv.h != null ? lv.h.toFixed(2) : "9999.00";
      const t = lv.t != null ? lv.t.toFixed(2) : "9999.00";
      const td = lv.td != null ? lv.td.toFixed(2) : "9999.00";
      const wd = lv.wd != null ? lv.wd.toFixed(2) : "9999.00";
      const ws = lv.ws != null ? lv.ws.toFixed(2) : "9999.00";
      lines.push(`${p},${h},${t},${td},${wd},${ws}`);
    }
    lines.push(`%END%`);
    const text = lines.join("\n");
    const blob = new Blob([text], { type: "text/plain" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = `${meta.station || "sounding"}_${meta.date.replace(/\s/g, "_")}.sharppy`;
    link.click();
    URL.revokeObjectURL(url);
  };

  // ── CM1 input_sounding export ──────────────────────────────
  const handleCm1 = () => {
    const profile = result.profile;
    if (!profile || profile.length === 0) return;
    const lines = [];
    // First line: sfc pressure (mb), sfc theta (K), sfc mixing ratio (g/kg)
    const sfcP = profile[0].p;
    const sfcT = profile[0].t + 273.15; // K
    const sfcTd = profile[0].td;
    // Compute sfc potential temperature: theta = T * (1000/p)^0.286
    const sfcTheta = sfcT * Math.pow(1000.0 / sfcP, 0.286);
    // Compute sfc mixing ratio from dewpoint and pressure (Bolton 1980)
    const es = 6.112 * Math.exp((17.67 * sfcTd) / (sfcTd + 243.5));
    const sfcQv = (621.97 * es) / (sfcP - es); // g/kg
    lines.push(`${sfcP.toFixed(2)}  ${sfcTheta.toFixed(2)}  ${sfcQv.toFixed(2)}`);
    // Subsequent lines: height AGL (m), theta (K), qv (g/kg), u (m/s), v (m/s)
    const sfcH = profile[0].h;
    for (const lv of profile) {
      const hAgl = lv.h - sfcH;
      const tk = lv.t + 273.15;
      const theta = tk * Math.pow(1000.0 / lv.p, 0.286);
      const e = 6.112 * Math.exp((17.67 * lv.td) / (lv.td + 243.5));
      const qv = (621.97 * e) / (lv.p - e);
      let u = 0, v = 0;
      if (lv.wd != null && lv.ws != null) {
        const wsMs = lv.ws * 0.51444; // kt → m/s
        const wdRad = (lv.wd * Math.PI) / 180;
        u = -wsMs * Math.sin(wdRad);
        v = -wsMs * Math.cos(wdRad);
      }
      lines.push(`${hAgl.toFixed(1)}  ${theta.toFixed(2)}  ${qv.toFixed(2)}  ${u.toFixed(2)}  ${v.toFixed(2)}`);
    }
    const text = lines.join("\n");
    const blob = new Blob([text], { type: "text/plain" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = `input_sounding_${meta.station || "sounding"}_${meta.date.replace(/\s/g, "_")}`;
    link.click();
    URL.revokeObjectURL(url);
  };

  const handleCopyLink = async () => {
    try {
      await navigator.clipboard.writeText(window.location.href);
      setLinkCopied(true);
      setTimeout(() => setLinkCopied(false), 2000);
    } catch {
      // Fallback for non-HTTPS
      const ta = document.createElement("textarea");
      ta.value = window.location.href;
      document.body.appendChild(ta);
      ta.select();
      document.execCommand("copy");
      document.body.removeChild(ta);
      setLinkCopied(true);
      setTimeout(() => setLinkCopied(false), 2000);
    }
  };

  return (
    <div className="results-view">
      {/* Meta bar */}
      <div className="rv-meta-bar">
        <div className="rv-meta-items">
          <span className="rv-meta-tag rv-meta-station">{meta.station || meta.source.toUpperCase()}</span>
          <span className="rv-meta-text">{meta.stationName}</span>
          <span className="rv-meta-sep" />
          <span className="rv-meta-text">{meta.date}</span>
          <span className="rv-meta-sep" />
          <span className="rv-meta-text">{meta.levels} levels</span>
          <span className="rv-meta-text rv-meta-dim">
            {meta.sfcPressure}&ndash;{meta.topPressure} hPa
          </span>
          {meta.vadRadar && (
            <span className="rv-meta-tag rv-meta-vad" title={meta.vadTime ? `VAD valid: ${meta.vadTime}` : ""}>
              VAD: {meta.vadRadar}
            </span>
          )}
        </div>
        <div className="rv-meta-actions">
          <button
            className={`rv-btn ${zoomed ? "rv-btn-active" : ""}`}
            onClick={() => setZoomed((z) => !z)}
            title={zoomed ? "Fit to width" : "Zoom to full resolution"}
          >
            {zoomed ? <ZoomOut size={14} /> : <ZoomIn size={14} />}
          </button>
          <button className="rv-btn" onClick={handleDownload} title="Download PNG">
            <Download size={14} />
          </button>
          <div className="rv-export-wrap" ref={exportRef}>
            <button className={`rv-btn ${exportOpen ? "rv-btn-active" : ""}`} onClick={() => setExportOpen((v) => !v)} title="Export sounding data">
              <Download size={14} />
              <span className="rv-btn-label">Export</span>
              <ChevronDown size={10} />
            </button>
            {exportOpen && (
              <div className="rv-export-menu">
                <button onClick={() => { handleCsvExport(); setExportOpen(false); }}>
                  <FileSpreadsheet size={13} />
                  <div>
                    <span className="rv-export-title">CSV Parameters</span>
                    <span className="rv-export-desc">All computed parameters in spreadsheet format</span>
                  </div>
                </button>
                <button onClick={() => { handleSharppy(); setExportOpen(false); }}>
                  <FileText size={13} />
                  <div>
                    <span className="rv-export-title">SHARPpy Format</span>
                    <span className="rv-export-desc">Raw profile for SHARPpy analysis software</span>
                  </div>
                </button>
                <button onClick={() => { handleCm1(); setExportOpen(false); }}>
                  <FileText size={13} />
                  <div>
                    <span className="rv-export-title">CM1 input_sounding</span>
                    <span className="rv-export-desc">Profile for Cloud Model 1 simulations</span>
                  </div>
                </button>
              </div>
            )}
          </div>
          <button className={`rv-btn ${linkCopied ? "rv-btn-active" : ""}`} onClick={handleCopyLink} title={linkCopied ? "Link copied!" : "Copy shareable link"}>
            {linkCopied ? <Check size={14} /> : <Link2 size={14} />}
          </button>
          <button className="rv-btn" onClick={handleFullscreen} title="Fullscreen">
            <Maximize2 size={14} />
          </button>
        </div>
      </div>

      {/* Map + Risk scan table */}
      {showMap && mapProps && <StationMap {...mapProps} />}
      {showRisk && <RiskTable riskData={riskData} />}

      {/* Time-Series Chart */}
      {showTimeSeries && (
        <TimeSeriesChart station={selectedStation} source={source} onClose={onCloseTimeSeries} />
      )}

      {/* Comparison View */}
      {showCompare && (
        <ComparisonView stations={stations || []} onClose={onCloseCompare} historyData={compareHistoryData} onHistoryConsumed={onCompareHistoryConsumed} />
      )}

      {/* VWP Display */}
      {showVwp && (
        <VwpDisplay stations={stations || []} selectedStation={selectedStation} onClose={onCloseVwp} />
      )}

      {/* Plot */}
      <div
        ref={plotRef}
        className={`rv-plot-wrap ${zoomed ? "rv-plot-zoomed" : ""}`}
        onMouseDown={handleMouseDown}
        onMouseMove={handleMouseMove}
        onMouseUp={handleMouseUp}
        onMouseLeave={handleMouseUp}
      >
        <img
          src={`data:image/png;base64,${image}`}
          alt="Sounding analysis plot"
          className="rv-plot-img"
          draggable={false}
        />
      </div>

      {/* Parameters */}
      <div className="rv-params">
        <ParamSection
          title="Thermodynamic"
          icon={<Thermometer size={14} />}
        >
          <ParamCard label="SB CAPE" value={params.sbCape} unit="J/kg" color="#f97316" desc="Surface-Based CAPE — total buoyant energy for a parcel lifted from the surface. Higher values indicate stronger updraft potential. >1000 notable, >3000 extreme." />
          <ParamCard label="SB CIN" value={params.sbCin} unit="J/kg" desc="Surface-Based CIN — energy needed to lift a surface parcel to its LFC. More negative = stronger cap. Values < -50 often inhibit convection initiation." />
          <ParamCard label="SB LCL" value={params.sbLclM} unit="m AGL" desc="Surface-Based Lifted Condensation Level — height where a surface parcel reaches saturation. Lower LCL (<1000m) favors tornadoes." />
          <ParamCard label="MU CAPE" value={params.muCape} unit="J/kg" desc="Most-Unstable CAPE — CAPE computed for the parcel with highest θe in the lowest 300 hPa. Represents the maximum buoyancy available." />
          <ParamCard label="MU CIN" value={params.muCin} unit="J/kg" desc="Most-Unstable CIN — inhibition for the most-unstable parcel. Useful when elevated convection is possible." />
          <ParamCard label="MU LCL" value={params.muLclM} unit="m AGL" desc="Most-Unstable LCL — condensation height for the MU parcel. May differ from SB LCL when the most unstable parcel is elevated." />
          <ParamCard label="ML CAPE" value={params.mlCape} unit="J/kg" color="#d946ef" desc="Mixed-Layer CAPE — CAPE for a parcel representing the mean of the lowest 100 hPa. Best estimate for surface-based storms in a well-mixed boundary layer." />
          <ParamCard label="ML CIN" value={params.mlCin} unit="J/kg" desc="Mixed-Layer CIN — inhibition for the ML parcel. More representative than SB CIN in the afternoon when the boundary layer is mixed." />
          <ParamCard label="ML LCL" value={params.mlLclM} unit="m AGL" desc="Mixed-Layer LCL — cloud base height for the ML parcel. Lower values increase tornado probability; <1000m is favorable." />
          <ParamCard label="DCAPE" value={params.dcape} unit="J/kg" desc="Downdraft CAPE — energy available for downdrafts. Higher values (>800) indicate strong outflow winds and potential for damaging gusts." />
          <ParamCard label="ECAPE" value={params.ecape} unit="J/kg" color="#06b6d4" desc="Entraining CAPE — CAPE adjusted for entrainment of dry environmental air (Peters et al. 2023). More physically realistic than standard CAPE." />
          <ParamCard label="3CAPE" value={params.cape3km} unit="J/kg" color="#fb923c" desc="0–3 km CAPE — buoyant energy in the lowest 3 km (MU parcel). Higher values indicate stronger low-level accelerations; key for tornado intensity." />
          <ParamCard label="6CAPE" value={params.cape6km} unit="J/kg" color="#facc15" desc="0–6 km CAPE — buoyant energy in the lowest 6 km (MU parcel). Indicates how quickly an updraft develops in the mid-levels." />
          <ParamCard label="DCIN" value={params.dcin} unit="J/kg" color="#818cf8" desc="Downdraft CIN — inhibition of downdrafts reaching the surface. More negative = stronger capping of downdrafts. Near 0 means downdrafts easily penetrate to the surface." />
          <ParamCard label="NCAPE" value={params.ncape} unit="J/kg/m" color="#38bdf8" desc="Normalized CAPE — MUCAPE divided by the LFC-to-EL depth. Measures buoyancy intensity per unit depth. >0.3 is very buoyant; indicates narrow, intense updrafts." />
          <ParamCard label="STP" value={params.stp} unit="" color="#60a5fa" desc="Significant Tornado Parameter — composite index combining CAPE, SRH, shear, and LCL. Values ≥1 suggest significant (EF2+) tornado environment. Higher = more favorable." />
          <ParamCard label="SCP" value={params.scp} unit="" color="#f59e0b" desc="Supercell Composite Parameter — combines CAPE, deep shear, and SRH. Values ≥1 support supercells; >4 strongly favors discrete supercells." />
          <ParamCard label="SHIP" value={params.ship} unit="" color="#10b981" desc="Significant Hail Parameter — composite for significant hail (≥2 in.). Values ≥1 indicate potential; >1.5 strongly favors significant hail." />
          <ParamCard label="DCP" value={params.dcp} unit="" color="#a78bfa" desc="Derecho Composite Parameter — combines DCAPE, MUCAPE, shear, and mean wind. Values ≥2 suggest potential for long-lived wind events (derechos)." />
        </ParamSection>

        <ParamSection
          title="Lapse Rates & Moisture"
          icon={<Droplets size={14} />}
        >
          <ParamCard label="LR 0-3 km" value={params.lr03} unit="C/km" desc="0–3 km Lapse Rate — temperature decrease per km in the lowest 3 km. Values >7°C/km are steep; >8°C/km nearly absolute unstable. Key for tornado environments." />
          <ParamCard label="LR 3-6 km" value={params.lr36} unit="C/km" desc="3–6 km Lapse Rate — mid-level lapse rate. Steeper values (>7°C/km) enhance CAPE and updraft strength. >8°C/km is extreme instability." />
          <ParamCard label="PWAT" value={params.pwat} unit="mm" desc="Precipitable Water — total column water vapor. Higher values mean more moisture available for heavy rainfall. >50mm is extremely high for CONUS." />
          <ParamCard label="FRZ Level" value={params.frzLevel} unit="m AGL" desc="Freezing Level — height of the 0°C isotherm AGL. Affects hail size (higher FRZ = more melting) and snow levels." />
          <ParamCard label="WB Zero" value={params.wbo} unit="m AGL" desc="Wet-Bulb Zero Height — where the wet-bulb temperature crosses 0°C. Better predictor of hail reaching the surface than the freezing level. <2500m favors large hail." />
          <ParamCard label="RH 0-1 km" value={params.rh01} unit="%" desc="Relative Humidity 0–1 km — low-level moisture. Higher values (>80%) favor tornado development by reducing evaporative cooling of rain." />
          <ParamCard label="RH 1-3 km" value={params.rh13} unit="%" desc="Relative Humidity 1–3 km — mid-low moisture. Dry layers here (<50%) enhance DCAPE and outflow potential." />
          <ParamCard label="RH 3-6 km" value={params.rh36} unit="%" desc="Relative Humidity 3–6 km — mid-level moisture. Dry air here entrains into storms, increasing evaporative cooling and downdraft strength." />
        </ParamSection>

        <ParamSection
          title="Kinematic"
          icon={<ArrowUpDown size={14} />}
        >
          <ParamCard label="BWD 0-1 km" value={params.bwd1km} unit="kt" color="#ef4444" desc="0–1 km Bulk Wind Difference — magnitude of the wind shear vector over the lowest 1 km. >15 kt supports organized storms; >20 kt favors tornadoes." />
          <ParamCard label="BWD 0-3 km" value={params.bwd3km} unit="kt" color="#f97316" desc="0–3 km Bulk Wind Difference — shear in the low-to-mid levels. Important for mesocyclone development. >30 kt favors supercells." />
          <ParamCard label="BWD 0-6 km" value={params.bwd6km} unit="kt" color="#eab308" desc="0–6 km Bulk Wind Difference — deep-layer shear. The primary discriminator between organized and disorganized convection. >40 kt strongly favors supercells." />
          <ParamCard label="SRH 500m" value={params.srh500m} unit="m²/s²" desc="0–500m Storm-Relative Helicity — streamwise vorticity in the lowest 500m relative to storm motion. Key for tornado potential. >150 is significant." />
          <ParamCard label="SRH 0-1 km" value={params.srh1km} unit="m²/s²" color="#ef4444" desc="0–1 km Storm-Relative Helicity — measures the rotational potential of a storm's updraft in the lowest 1 km. >100 favors mesocyclones; >300 strongly favors tornadoes." />
          <ParamCard label="SRH 0-3 km" value={params.srh3km} unit="m²/s²" color="#f97316" desc="0–3 km Storm-Relative Helicity — total low-level rotational potential. >200 favors strong mesocyclones; >400 is extreme. Used in STP and SCP composites." />
          <ParamCard label="Eff. SRH" value={params.esrh} unit="m²/s²" color="#2dd4bf" desc="Effective SRH — storm-relative helicity computed within the effective inflow layer (where CAPE ≥ 100 and CIN > -250). More physically meaningful than fixed-layer SRH." />
          <ParamCard label="Eff. BWD" value={params.ebwd} unit="kt" color="#34d399" desc="Effective Bulk Wind Difference — shear from the effective inflow base to half the MU EL height. Better discriminator for supercells than fixed 0-6 km shear." />
          <ParamCard label="EIL Base" value={params.eilBot} unit="m AGL" color="#a7f3d0" desc="Effective Inflow Layer base — lowest level where CAPE ≥ 100 J/kg and CIN > -250 J/kg. Identifies the true inflow layer for storms." />
          <ParamCard label="EIL Top" value={params.eilTop} unit="m AGL" color="#a7f3d0" desc="Effective Inflow Layer top — highest contiguous level meeting the CAPE/CIN thresholds. Deeper EIL = deeper inflow available for storms." />
        </ParamSection>
      </div>

      {/* Text Summary */}
      <div className="rv-summary-section">
        <button className="rv-summary-toggle" onClick={() => setShowSummary((s) => !s)}>
          <FileText size={14} />
          <span>Sounding Text Summary</span>
          <ChevronDown size={12} className={showSummary ? "rv-chev-open" : ""} />
        </button>
        {showSummary && (
          <div className="rv-summary-body">
            {generateSoundingSummary(params, meta).split("\n\n").map((p, i) => (
              <p key={i}>{p}</p>
            ))}
          </div>
        )}
      </div>
    </div>
  );
}
