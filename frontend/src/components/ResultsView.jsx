import { useState, useRef, useCallback, useEffect } from "react";
import {
  BarChart3,
  Loader2,
  AlertTriangle,
  Download,
  Maximize2,
  ZoomIn,
  ZoomOut,
  Wind,
  Thermometer,
  Droplets,
  FileSpreadsheet,
  Link2,
  Check,
  FileText,
  ChevronDown,
  Printer,
  ArrowUp,
  ArrowDown,
  ShieldOff,
  Cloud,
  TrendingUp,
  RotateCcw,
  Target,
  Snowflake,
  Play,
  RefreshCw,
} from "lucide-react";

/* ── Icon map for summary section cards ─────────────── */
const SECTION_ICONS = {
  thermo: Thermometer,
  llcape: ArrowUp,
  cin: ShieldOff,
  lcl: Cloud,
  lapse: TrendingUp,
  moisture: Droplets,
  dlshear: Wind,
  llshear: RotateCcw,
  composite: BarChart3,
  dcape: ArrowDown,
  hail: Snowflake,
  mode: Target,
};
import { lazy, Suspense } from "react";
const StationMap       = lazy(() => import("./StationMap"));
const TimeSeriesChart  = lazy(() => import("./TimeSeriesChart"));
const ComparisonView   = lazy(() => import("./ComparisonView"));
const VwpDisplay       = lazy(() => import("./VwpDisplay"));
const InteractiveSkewT = lazy(() => import("./InteractiveSkewT"));
const InteractiveHodograph = lazy(() => import("./InteractiveHodograph"));
const SoundingAnimator = lazy(() => import("./SoundingAnimator"));
const SoundingTimeline = lazy(() => import("./SoundingTimeline"));
import "./ResultsView.css";
import { useAppStore } from "../store/useAppStore";
import { useShallow } from "zustand/shallow";
import RiskTable from "./RiskTable";
import { generateSoundingSummary } from "../utils/soundingSummary";
import { exportSoundingReportPng } from "../utils/exportImage";

/* Original `generateSoundingSummary` and inline `RiskTable` blocks were
   moved to `utils/soundingSummary.js` and `components/RiskTable.jsx` in
   the Phase 1 refactor. */

export default function ResultsView({
  autoRefresh, onToggleAutoRefresh, refreshInterval, onRefreshIntervalChange,
  onTimelineSelect, onRiskStationSelect, mapProps
}) {
  // Subscribe with shallow equality so the (very large) ResultsView only
  // re-renders when one of these specific fields actually changes — not on
  // every store mutation (e.g. toggling unrelated UI panels).
  const {
    stations, error, loading, result, riskData, showRisk, showMap,
    showTimeSeries, setShowTimeSeries, showCompare, setShowCompare,
    showVwp, setShowVwp, compareHistoryData, setCompareHistoryData,
    selectedStation, source, theme, lastParams,
  } = useAppStore(useShallow((s) => ({
    stations: s.stations,
    error: s.error,
    loading: s.loading,
    result: s.result,
    riskData: s.riskData,
    showRisk: s.showRisk,
    showMap: s.showMap,
    showTimeSeries: s.showTimeSeries,
    setShowTimeSeries: s.setShowTimeSeries,
    showCompare: s.showCompare,
    setShowCompare: s.setShowCompare,
    showVwp: s.showVwp,
    setShowVwp: s.setShowVwp,
    compareHistoryData: s.compareHistoryData,
    setCompareHistoryData: s.setCompareHistoryData,
    selectedStation: s.selectedStation,
    source: s.source,
    theme: s.theme,
    lastParams: s.lastParams,
  })));

  const onCloseTimeSeries = () => setShowTimeSeries(false);
  const onCloseCompare = () => setShowCompare(false);
  const onCloseVwp = () => setShowVwp(false);
  const onCompareHistoryConsumed = () => setCompareHistoryData(null);

  /* ── Hooks — must be called unconditionally before any return ── */
  const [zoomed, setZoomed] = useState(false);
  const [linkCopied, setLinkCopied] = useState(false);
  const [exportOpen, setExportOpen] = useState(false);
  const [showSummary, setShowSummary] = useState(false);
  const [interactiveMode, setInteractiveMode] = useState(false);
  const [showAnimator, setShowAnimator] = useState(false);
  const [reportBusy, setReportBusy] = useState(false);
  const exportRef = useRef(null);
  const plotRef = useRef(null);
  const dragRef = useRef({ dragging: false, startX: 0, startY: 0, scrollLeft: 0, scrollTop: 0 });

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

  if (error) {
    return (
      <div className="results-view">
        <Suspense fallback={null}>
          {showMap && mapProps && <div id="section-map"><StationMap {...mapProps} /></div>}
        </Suspense>
        {showRisk && <div id="section-risk"><RiskTable riskData={riskData} onStationSelect={onRiskStationSelect} /></div>}
        <Suspense fallback={null}>
          {showTimeSeries && (
            <div id="section-timeseries"><TimeSeriesChart station={selectedStation} source={source} onClose={onCloseTimeSeries} /></div>
          )}
          {showCompare && (
            <div id="section-compare"><ComparisonView stations={stations || []} onClose={onCloseCompare} historyData={compareHistoryData} onHistoryConsumed={onCompareHistoryConsumed} /></div>
          )}
          {showVwp && (
            <div id="section-vwp"><VwpDisplay stations={stations || []} selectedStation={selectedStation} onClose={onCloseVwp} /></div>
          )}
        </Suspense>
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
        <Suspense fallback={null}>
          {showMap && mapProps && <div id="section-map"><StationMap {...mapProps} /></div>}
        </Suspense>
        {showRisk && <div id="section-risk"><RiskTable riskData={riskData} onStationSelect={onRiskStationSelect} /></div>}
        <Suspense fallback={null}>
          {showTimeSeries && (
            <div id="section-timeseries"><TimeSeriesChart station={selectedStation} source={source} onClose={onCloseTimeSeries} /></div>
          )}
          {showCompare && (
            <div id="section-compare"><ComparisonView stations={stations || []} onClose={onCloseCompare} historyData={compareHistoryData} onHistoryConsumed={onCompareHistoryConsumed} /></div>
          )}
          {showVwp && (
            <div id="section-vwp"><VwpDisplay stations={stations || []} selectedStation={selectedStation} onClose={onCloseVwp} /></div>
          )}
        </Suspense>
        <div id="section-sounding" className="rv-state rv-loading">
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
        <Suspense fallback={null}>
          {showMap && mapProps && <div id="section-map"><StationMap {...mapProps} /></div>}
        </Suspense>
        {showRisk && <div id="section-risk"><RiskTable riskData={riskData} onStationSelect={onRiskStationSelect} /></div>}
        <Suspense fallback={null}>
          {showTimeSeries && (
            <div id="section-timeseries"><TimeSeriesChart station={selectedStation} source={source} onClose={onCloseTimeSeries} /></div>
          )}
          {showCompare && (
            <div id="section-compare"><ComparisonView stations={stations || []} onClose={onCloseCompare} historyData={compareHistoryData} onHistoryConsumed={onCompareHistoryConsumed} /></div>
          )}
          {showVwp && (
            <div id="section-vwp"><VwpDisplay stations={stations || []} selectedStation={selectedStation} onClose={onCloseVwp} /></div>
          )}
        </Suspense>
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

  // ── JSON export ─────────────────────────────────────────────
  const handleJsonExport = () => {
    const payload = {
      meta: { station: meta.station, date: meta.date, source: meta.source, stationName: meta.stationName },
      params,
      profile: result.profile,
    };
    const json = JSON.stringify(payload, null, 2);
    const blob = new Blob([json], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = `sounding_${meta.station || "analysis"}_${meta.date.replace(/\s/g, "_")}.json`;
    link.click();
    URL.revokeObjectURL(url);
  };

  const handleFullReport = async () => {
    if (reportBusy) return;
    setReportBusy(true);
    try {
      await exportSoundingReportPng({ image, params, meta });
    } catch (err) {
      console.error("Report export failed:", err);
    }
    setReportBusy(false);
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
      {/* Map + Risk scan table */}
      <Suspense fallback={null}>
        {showMap && mapProps && <div id="section-map"><StationMap {...mapProps} /></div>}
      </Suspense>
      {showRisk && <div id="section-risk"><RiskTable riskData={riskData} onStationSelect={onRiskStationSelect} /></div>}

      {/* Time-Series Chart / Comparison / VWP */}
      <Suspense fallback={null}>
        {showTimeSeries && (
          <div id="section-timeseries"><TimeSeriesChart station={selectedStation} source={source} onClose={onCloseTimeSeries} /></div>
        )}
        {showCompare && (
          <div id="section-compare"><ComparisonView stations={stations || []} onClose={onCloseCompare} historyData={compareHistoryData} onHistoryConsumed={onCompareHistoryConsumed} /></div>
        )}
        {showVwp && (
          <div id="section-vwp"><VwpDisplay stations={stations || []} selectedStation={selectedStation} onClose={onCloseVwp} /></div>
        )}
      </Suspense>

      {/* Meta bar — above the sounding plot */}
      <div id="section-sounding" className="rv-meta-bar">
        <div className="rv-meta-row-top">
          <span className="rv-meta-station">{meta.station || meta.source.toUpperCase()}</span>
          <span className="rv-meta-name">{meta.stationName}</span>
          <span className="rv-meta-date">{meta.date}</span>
        </div>
        <div className="rv-meta-row-bottom">
          <div className="rv-meta-details">
            <span className="rv-meta-chip">{meta.levels} <span className="rv-meta-chip-value">levels</span></span>
            <span className="rv-meta-chip">{meta.sfcPressure}&ndash;{meta.topPressure} <span className="rv-meta-chip-value">hPa</span></span>
            {meta.vadRadar && (
              <span className="rv-meta-chip rv-meta-vad" title={meta.vadTime ? `VAD valid: ${meta.vadTime}` : ""}>
                VAD <span className="rv-meta-chip-value">{meta.vadRadar}</span>
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
                <button onClick={() => { handleJsonExport(); setExportOpen(false); }}>
                  <FileText size={13} />
                  <div>
                    <span className="rv-export-title">JSON</span>
                    <span className="rv-export-desc">Full params + profile data as JSON</span>
                  </div>
                </button>
                <button onClick={() => { handleFullReport(); setExportOpen(false); }} disabled={reportBusy}>
                  <Printer size={13} />
                  <div>
                    <span className="rv-export-title">{reportBusy ? "Generating…" : "Full Report PNG"}</span>
                    <span className="rv-export-desc">Screenshot of sounding + all parameters</span>
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
          <button className="rv-btn" onClick={() => window.print()} title="Print multi-panel layout">
            <Printer size={14} />
          </button>
          <button
            className={`rv-btn ${interactiveMode ? "rv-btn-active" : ""}`}
            onClick={() => setInteractiveMode((v) => !v)}
            title={interactiveMode ? "Switch to static plot" : "Switch to interactive Skew-T"}
          >
            <BarChart3 size={14} />
            <span className="rv-btn-label">{interactiveMode ? "Static" : "Interactive"}</span>
          </button>
          {(source === "bufkit" || source === "psu") && (
            <button
              className={`rv-btn ${showAnimator ? "rv-btn-active" : ""}`}
              onClick={() => setShowAnimator((v) => !v)}
              title="Animate through forecast hours"
            >
              <Play size={14} />
              <span className="rv-btn-label">Animate</span>
            </button>
          )}
          {onToggleAutoRefresh && (
            <div className="rv-autorefresh-wrap">
              <button
                className={`rv-btn ${autoRefresh ? "rv-btn-active" : ""}`}
                onClick={onToggleAutoRefresh}
                title={autoRefresh ? `Auto-refreshing every ${refreshInterval / 60000} min` : "Enable auto-refresh"}
              >
                <RefreshCw size={14} className={autoRefresh ? "spin-slow" : ""} />
              </button>
              {autoRefresh && (
                <select
                  className="rv-autorefresh-sel"
                  value={refreshInterval}
                  onChange={(e) => onRefreshIntervalChange(Number(e.target.value))}
                  title="Refresh interval"
                >
                  <option value={60000}>1 min</option>
                  <option value={120000}>2 min</option>
                  <option value={300000}>5 min</option>
                  <option value={600000}>10 min</option>
                  <option value={900000}>15 min</option>
                </select>
              )}
            </div>
          )}
        </div>
        </div>
      </div>
      {/* Sounding Timeline */}
      {onTimelineSelect && (
        <Suspense fallback={null}>
          <SoundingTimeline
            station={selectedStation}
            currentDate={meta?.date}
            onSelectTime={onTimelineSelect}
            source={source}
          />
        </Suspense>
      )}
      {/* Sounding Animator */}
      {showAnimator && (
        <Suspense fallback={<div className="rv-state rv-loading"><Loader2 size={20} className="spin" /><span>Loading animator…</span></div>}>
          <SoundingAnimator
            station={lastParams?.station || selectedStation}
            model={lastParams?.model}
            source={lastParams?.source || source}
            date={lastParams?.date}
            theme={theme || "dark"}
            onClose={() => setShowAnimator(false)}
          />
        </Suspense>
      )}
      {/* Plot — static PNG or interactive Skew-T */}
      {interactiveMode ? (
        <Suspense fallback={<div className="rv-state rv-loading"><Loader2 size={20} className="spin" /><span>Loading interactive Skew-T…</span></div>}>
          <InteractiveSkewT
            profile={result.profile}
            sbParcel={result.sbParcel}
            mlParcel={result.mlParcel}
            params={params}
            theme={theme || "dark"}
          />
          <InteractiveHodograph
            profile={result.profile}
            params={params}
            theme={theme || "dark"}
          />
          {/* Piecewise CAPE strip */}
          {params.piecewiseCape && params.piecewiseCape.length > 0 && (
            <div className="rv-piecewise-cape">
              <div className="rv-pw-title">Piecewise CAPE (50 hPa layers)</div>
              <div className="rv-pw-strip">
                {params.piecewiseCape.map((layer, i) => {
                  const maxCape = Math.max(...params.piecewiseCape.map(l => l.cape), 1);
                  const barW = Math.max(2, (layer.cape / maxCape) * 100);
                  const alpha = Math.min(1, layer.cape / maxCape * 0.8 + 0.2);
                  return (
                    <div key={i} className="rv-pw-row" title={`${layer.p_bot}\u2013${layer.p_top} hPa: CAPE ${layer.cape} J/kg`}>
                      <span className="rv-pw-label">{layer.p_bot}</span>
                      <div className="rv-pw-bar-bg">
                        <div className="rv-pw-bar" style={{ width: `${barW}%`, opacity: alpha }} />
                      </div>
                      <span className="rv-pw-val">{Math.round(layer.cape)}</span>
                    </div>
                  );
                })}
              </div>
            </div>
          )}
        </Suspense>
      ) : (
      <div
        ref={plotRef}
        className={`rv-plot-wrap ${zoomed ? "rv-plot-zoomed" : ""}`}
        onMouseDown={handleMouseDown}
        onMouseMove={handleMouseMove}
        onMouseUp={handleMouseUp}
        onMouseLeave={handleMouseUp}
        onTouchStart={(e) => {
          if (e.touches.length === 1 && zoomed) {
            const t = e.touches[0];
            dragRef.current = {
              dragging: true,
              startX: t.clientX,
              startY: t.clientY,
              scrollLeft: plotRef.current?.scrollLeft || 0,
              scrollTop: plotRef.current?.scrollTop || 0,
            };
          }
        }}
        onTouchMove={(e) => {
          const d = dragRef.current;
          if (d.dragging && e.touches.length === 1) {
            const t = e.touches[0];
            const el = plotRef.current;
            if (el) {
              el.scrollLeft = d.scrollLeft - (t.clientX - d.startX);
              el.scrollTop = d.scrollTop - (t.clientY - d.startY);
            }
          }
        }}
        onTouchEnd={() => { dragRef.current.dragging = false; }}
        style={{ touchAction: zoomed ? "none" : "auto" }}
      >
        <img
          src={`data:image/png;base64,${image}`}
          alt="Sounding analysis plot"
          className="rv-plot-img"
          draggable={false}
        />
      </div>
      )}

      {/* Text Summary */}
      <div className="rv-summary-section">
        <button className="rv-summary-toggle" onClick={() => setShowSummary((s) => !s)}>
          <FileText size={14} />
          <span>Sounding Text Summary</span>
          <ChevronDown size={12} className={showSummary ? "rv-chev-open" : ""} />
        </button>
        {showSummary && (() => {
          const summary = generateSoundingSummary(params, meta);
          return (
            <div className="rv-summary-body">
              {/* ── Header banner ── */}
              <div className="rv-sum-banner">
                <div className="rv-sum-banner-left">
                  <span className="rv-sum-station">{summary.station}</span>
                  {summary.date && <span className="rv-sum-date">{summary.date}</span>}
                </div>
                {summary.threatLevel !== "none" && (
                  <span className="rv-sum-threat-badge" style={{ background: summary.threatColor + "22", color: summary.threatColor, borderColor: summary.threatColor + "44" }}>
                    {summary.threatLevel} RISK
                  </span>
                )}
              </div>

              {/* ── Threat callout ── */}
              {summary.threats.length > 0 && (
                <div className="rv-sum-callout" style={{ borderLeftColor: summary.threatColor }}>
                  <span className="rv-sum-callout-label">Primary Hazards</span>
                  <div className="rv-sum-callout-tags">
                    {summary.threats.map((t, i) => (
                      <span key={i} className="rv-sum-hazard-tag" style={{ background: summary.threatColor + "18", color: summary.threatColor }}>{t}</span>
                    ))}
                  </div>
                </div>
              )}

              {/* ── Analysis sections ── */}
              <div className="rv-sum-sections">
                {summary.sections.map((sec) => {
                  const IconComp = SECTION_ICONS[sec.id];
                  return (
                  <div key={sec.id} className="rv-sum-card">
                    <div className="rv-sum-card-hdr">
                      {IconComp && <IconComp size={14} className="rv-sum-card-icon" />}
                      <span className="rv-sum-card-title">{sec.title}</span>
                    </div>
                    {sec.vals && sec.vals.length > 0 && (
                      <div className="rv-sum-card-vals">
                        {sec.vals.map((v) => (
                          <span key={v.k} className="rv-sum-kv">
                            <span className="rv-sum-kv-k">{v.k}</span>
                            <span className="rv-sum-kv-v">{v.v}</span>
                          </span>
                        ))}
                      </div>
                    )}
                    <p className="rv-sum-card-text">{sec.text}</p>
                  </div>
                  );
                })}
              </div>
            </div>
          );
        })()}
      </div>


    </div>
  );
}
