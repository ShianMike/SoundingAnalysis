import { useState, useRef, useCallback } from "react";
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
} from "lucide-react";
import "./ResultsView.css";

function ParamCard({ label, value, unit, color }) {
  return (
    <div className="param-card">
      <span className="param-label">{label}</span>
      <span className="param-value" style={color ? { color } : {}}>
        {value ?? "---"}
      </span>
      {unit && value != null && <span className="param-unit">{unit}</span>}
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
        <h3>Tornado Risk Scan</h3>
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
              <th className="rv-rt-num">Raw</th>
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
                <td className="rv-rt-num">{s.raw.toFixed(2)}</td>
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

export default function ResultsView({ result, loading, error, riskData }) {
  if (error) {
    return (
      <div className="results-view">
        <RiskTable riskData={riskData} />
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
        <RiskTable riskData={riskData} />
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
        <RiskTable riskData={riskData} />
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
  const plotRef = useRef(null);
  const dragRef = useRef({ dragging: false, startX: 0, startY: 0, scrollLeft: 0, scrollTop: 0 });

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
          <button className="rv-btn" onClick={handleFullscreen} title="Fullscreen">
            <Maximize2 size={14} />
          </button>
        </div>
      </div>

      {/* Risk scan table */}
      <RiskTable riskData={riskData} />

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
          <ParamCard label="SB CAPE" value={params.sbCape} unit="J/kg" color="#f97316" />
          <ParamCard label="SB CIN" value={params.sbCin} unit="J/kg" />
          <ParamCard label="SB LCL" value={params.sbLclM} unit="m AGL" />
          <ParamCard label="MU CAPE" value={params.muCape} unit="J/kg" />
          <ParamCard label="MU CIN" value={params.muCin} unit="J/kg" />
          <ParamCard label="MU LCL" value={params.muLclM} unit="m AGL" />
          <ParamCard label="ML CAPE" value={params.mlCape} unit="J/kg" color="#d946ef" />
          <ParamCard label="ML CIN" value={params.mlCin} unit="J/kg" />
          <ParamCard label="ML LCL" value={params.mlLclM} unit="m AGL" />
          <ParamCard label="DCAPE" value={params.dcape} unit="J/kg" />
          <ParamCard label="STP" value={params.stp} unit="" color="#60a5fa" />
        </ParamSection>

        <ParamSection
          title="Lapse Rates & Moisture"
          icon={<Droplets size={14} />}
        >
          <ParamCard label="LR 0-3 km" value={params.lr03} unit="C/km" />
          <ParamCard label="LR 3-6 km" value={params.lr36} unit="C/km" />
          <ParamCard label="PWAT" value={params.pwat} unit="mm" />
          <ParamCard label="FRZ Level" value={params.frzLevel} unit="m AGL" />
          <ParamCard label="WB Zero" value={params.wbo} unit="m AGL" />
          <ParamCard label="RH 0-1 km" value={params.rh01} unit="%" />
          <ParamCard label="RH 1-3 km" value={params.rh13} unit="%" />
          <ParamCard label="RH 3-6 km" value={params.rh36} unit="%" />
        </ParamSection>

        <ParamSection
          title="Kinematic"
          icon={<ArrowUpDown size={14} />}
        >
          <ParamCard label="BWD 0-1 km" value={params.bwd1km} unit="kt" color="#ef4444" />
          <ParamCard label="BWD 0-3 km" value={params.bwd3km} unit="kt" color="#f97316" />
          <ParamCard label="BWD 0-6 km" value={params.bwd6km} unit="kt" color="#eab308" />
          <ParamCard label="SRH 500m" value={params.srh500m} unit="m²/s²" />
          <ParamCard label="SRH 0-1 km" value={params.srh1km} unit="m²/s²" color="#ef4444" />
          <ParamCard label="SRH 0-3 km" value={params.srh3km} unit="m²/s²" color="#f97316" />
        </ParamSection>
      </div>
    </div>
  );
}
