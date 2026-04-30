import { useState, useRef, useEffect } from "react";
import {
  Database,
  Loader2,
  ChevronDown,
  MessageSquarePlus,
  Github,
  Keyboard,
  Sun,
  Moon,
  GripVertical,
  ChevronUp,
  Sliders,
  Wrench,
  AlertTriangle,
} from "lucide-react";
import { fetchRiskScan, fetchForecastRiskScan } from "../api";
import { getFavorites, toggleFavorite } from "../favorites";
import "./ControlPanel.css";
import { useAppStore } from "../store/useAppStore";
import { useShallow } from "zustand/shallow";
import { useDraggable } from "../hooks/useDraggable";
import ToolsTab from "./ControlPanel/ToolsTab";
import ModifyTab from "./ControlPanel/ModifyTab";
import DataTab from "./ControlPanel/DataTab";
import { nearestNexrad } from "../config/constants";

/** Local thin wrapper preserving the existing `"KXXX"` return shape used here. */
function nearestNexradForVad(lat, lon) {
  const r = nearestNexrad(lat, lon);
  return r ? r.idK : null;
}

export default function ControlPanel({
  onSubmit,
  onRetry,
  connectError,
  onNavigateEnsemble,
  mapLatLon,
  mapStormMotion,
  mapBoundaryOrientation,
  urlParams,
  onNavigateUpload,
  onShowShortcuts,
}) {
  // Subscribe to only this component's slice of the store with shallow
  // equality so unrelated store mutations (result, error, lastParams, page,
  // compareHistoryData, …) don't trigger re-renders here.
  const {
    stations, sources, models, psuModels, loading, initialLoading,
    riskData, setRiskData, showRisk, toggleRisk,
    showHistory, toggleHistory,
    showMap, toggleMap,
    showTimeSeries, toggleTimeSeries,
    showCompare, toggleCompare,
    showVwp, toggleVwp,
    station, setStation,
    source, setSource,
    feedbackActive, toggleFeedback,
    theme, toggleTheme, colorblind, toggleColorblind,
  } = useAppStore(useShallow((s) => ({
    stations: s.stations,
    sources: s.sources,
    models: s.models,
    psuModels: s.psuModels,
    loading: s.loading,
    initialLoading: s.initialLoading,
    riskData: s.riskData,
    setRiskData: s.setRiskData,
    showRisk: s.showRisk,
    toggleRisk: s.toggleRisk,
    showHistory: s.showHistory,
    toggleHistory: s.toggleHistory,
    showMap: s.showMap,
    toggleMap: s.toggleMap,
    showTimeSeries: s.showTimeSeries,
    toggleTimeSeries: s.toggleTimeSeries,
    showCompare: s.showCompare,
    toggleCompare: s.toggleCompare,
    showVwp: s.showVwp,
    toggleVwp: s.toggleVwp,
    station: s.selectedStation,
    setStation: s.setSelectedStation,
    source: s.source,
    setSource: s.setSource,
    feedbackActive: s.showFeedback,
    toggleFeedback: s.toggleFeedback,
    theme: s.theme,
    toggleTheme: s.toggleTheme,
    colorblind: s.colorblind,
    toggleColorblind: s.toggleColorblind,
  })));

  const onFeedbackClick = toggleFeedback;
  const onToggleColorblind = toggleColorblind;
  const onToggleTheme = toggleTheme;
  const onToggleRisk = toggleRisk;
  const onToggleHistory = toggleHistory;
  const onToggleMap = toggleMap;
  const onToggleTimeSeries = toggleTimeSeries;
  const onToggleCompare = toggleCompare;
  const onToggleVwp = toggleVwp;
  const onRiskDataChange = setRiskData;
  /* ── Drag & collapse state for fullscreen-float mode ───────
     Logic lives in `hooks/useDraggable.js`. */
  const {
    dragPos,
    floatCollapsed,
    toggleFloatCollapsed,
    dragRef,
    dragHandlers,
  } = useDraggable();

  const [date, setDate] = useState(() => {
    // Convert YYYYMMDDHH back to datetime-local value for the input
    if (urlParams?.date && /^\d{10}$/.test(urlParams.date)) {
      const d = urlParams.date;
      return `${d.slice(0,4)}-${d.slice(4,6)}-${d.slice(6,8)}T${d.slice(8,10)}:00`;
    }
    return "";
  });
  const [lat, setLat] = useState(urlParams?.lat != null ? String(urlParams.lat) : "");
  const [lon, setLon] = useState(urlParams?.lon != null ? String(urlParams.lon) : "");
  const [model, setModel] = useState(urlParams?.model || "hrrr");
  const [fhour, setFhour] = useState(urlParams?.fhour != null ? String(urlParams.fhour) : "0");
  const [stationSearch, setStationSearch] = useState("");
  const [scanning, setScanning] = useState(false);
  const [scanMode, setScanMode] = useState("obs"); // "obs" | "forecast"
  const [fcstModel, setFcstModel] = useState("hrrr");
  const [fcstFhour, setFcstFhour] = useState("12");
  const [sortMode, setSortMode] = useState("az");
  const [favorites, setFavorites] = useState(() => getFavorites());
  const [soundingHour, setSoundingHour] = useState("latest");
  const listRef = useRef(null);

  // Surface modification state
  const [sfcModEnabled, setSfcModEnabled] = useState(false);
  const [sfcModT, setSfcModT] = useState("");
  const [sfcModTd, setSfcModTd] = useState("");
  const [sfcModWspd, setSfcModWspd] = useState("");
  const [sfcModWdir, setSfcModWdir] = useState("");

  // Custom storm motion state
  const [smEnabled, setSmEnabled] = useState(false);
  const [smDirection, setSmDirection] = useState("");
  const [smSpeed, setSmSpeed] = useState("");

  // VAD Wind Profile overlay state
  const [vadEnabled, setVadEnabled] = useState(false);

  // Storm-relative hodograph state
  const [srHodoEnabled, setSrHodoEnabled] = useState(false);

  // Profile smoothing state
  const [smoothEnabled, setSmoothEnabled] = useState(false);
  const [smoothSigma, setSmoothSigma] = useState("3");

  // Boundary line state
  const [boundaryEnabled, setBoundaryEnabled] = useState(false);
  const [boundaryOrientation, setBoundaryOrientation] = useState("");

  // Navigation tab state
  const [navTab, setNavTab] = useState("data");




  // Point sounding mode — set when user clicks map for arbitrary lat/lon
  const [pointMode, setPointMode] = useState(false);

  // Sync lat/lon and scroll into view when station changes (e.g. map click)
  useEffect(() => {
    if (!station) return;
    setPointMode(false);
    const stn = stations.find((s) => s.id === station);
    if (stn) {
      setLat(String(stn.lat));
      setLon(String(stn.lon));
    }
    // Scroll the selected station into view in the list
    if (listRef.current) {
      requestAnimationFrame(() => {
        const el = listRef.current?.querySelector(`.cp-station-item[data-id="${station}"]`);
        if (el) el.scrollIntoView({ block: "center", behavior: "smooth" });
      });
    }
  }, [station]); // eslint-disable-line react-hooks/exhaustive-deps

  // Sync lat/lon from map click — enter point sounding mode
  useEffect(() => {
    if (mapLatLon) {
      setLat(String(mapLatLon.lat));
      setLon(String(mapLatLon.lon));
      setPointMode(true);
    }
  }, [mapLatLon]);

  // Sync custom storm motion from map draw tool
  useEffect(() => {
    if (
      mapStormMotion &&
      mapStormMotion.direction != null &&
      mapStormMotion.speed != null
    ) {
      setSmEnabled(true);
      setSmDirection(String(Math.round(mapStormMotion.direction)));
      setSmSpeed(String(Math.round(mapStormMotion.speed)));
    }
  }, [mapStormMotion]);

  // Sync boundary orientation from map draw tool
  useEffect(() => {
    if (mapBoundaryOrientation != null) {
      setBoundaryEnabled(true);
      setBoundaryOrientation(String(Math.round(mapBoundaryOrientation)));
    }
  }, [mapBoundaryOrientation]);

  // Scroll selected station into view on mount
  useEffect(() => {
    if (listRef.current) {
      const active = listRef.current.querySelector(".cp-station-item.active");
      if (active) active.scrollIntoView({ block: "center" });
    }
  }, [stations]);

  const needsStation = source === "obs" || source === "bufkit" || source === "psu";
  const needsModel = source === "bufkit" || source === "psu" || pointMode;

  // Build risk lookup from scan results
  const riskMap = {};
  if (riskData) {
    riskData.stations.forEach((s, i) => {
      riskMap[s.id] = { rank: i + 1, stp: s.stp, raw: s.raw, cape: s.cape, srh: s.srh, bwd: s.bwd };
    });
  }

  // Merge station data with risk scores
  const mergedStations = stations.map((s) => ({
    ...s,
    risk: riskMap[s.id] || null,
  }));

  // Sort based on selected mode
  const sortedStations = [...mergedStations].sort((a, b) => {
    switch (sortMode) {
      case "za":
        return b.id.localeCompare(a.id);
      case "favs": {
        const aFav = favorites.includes(a.id) ? 0 : 1;
        const bFav = favorites.includes(b.id) ? 0 : 1;
        if (aFav !== bFav) return aFav - bFav;
        return a.id.localeCompare(b.id);
      }
      case "risk-high":
        // Stations with risk first (highest STP first), then unscanned at end
        if (a.risk && !b.risk) return -1;
        if (!a.risk && b.risk) return 1;
        if (a.risk && b.risk) return b.risk.stp - a.risk.stp;
        return a.id.localeCompare(b.id);
      case "risk-low":
        if (a.risk && !b.risk) return -1;
        if (!a.risk && b.risk) return 1;
        if (a.risk && b.risk) return a.risk.stp - b.risk.stp;
        return a.id.localeCompare(b.id);
      case "az":
      default:
        return a.id.localeCompare(b.id);
    }
  });

  const filteredStations = sortedStations.filter(
    (s) => {
      const matchesSearch = s.id.toLowerCase().includes(stationSearch.toLowerCase()) ||
        s.name.toLowerCase().includes(stationSearch.toLowerCase());
      return matchesSearch;
    }
  );

  const handleRiskScan = async () => {
    setScanning(true);
    try {
      let data;
      if (scanMode === "forecast") {
        data = await fetchForecastRiskScan({
          model: fcstModel,
          fhour: parseInt(fcstFhour) || 0,
        });
      } else {
        const dateParam = date ? date.replace(/[-T:]/g, "").slice(0, 10) : undefined;
        data = await fetchRiskScan(dateParam);
      }
      onRiskDataChange(data);
      setSortMode("risk-high");
      // Auto-open the map
      if (!showMap && onToggleMap) onToggleMap();
      // Auto-select the highest risk station
      if (data.stations.length > 0) {
        handleStationSelect(data.stations[0].id);
      }
    } catch (e) {
      console.error("Risk scan failed:", e);
    } finally {
      setScanning(false);
    }
  };

  const handleSubmit = (e) => {
    e.preventDefault();
    const params = { source };

    // Point sounding: send lat/lon only, no station
    if (pointMode && lat && lon) {
      params.lat = parseFloat(lat);
      params.lon = parseFloat(lon);
      // Point soundings need a model source
      if (source === "obs") params.source = "bufkit";
    } else if (needsStation) {
      params.station = station;
    }
    // Include lat/lon when available (for obs nearest-station lookup)
    if (!pointMode && lat && lon) {
      params.lat = parseFloat(lat);
      params.lon = parseFloat(lon);
    }
    if (date) params.date = date.replace(/[-T:]/g, "").slice(0, 10);
    if (needsModel) {
      params.model = model;
      params.fhour = parseInt(fhour) || 0;
    }

    // Surface modification
    if (sfcModEnabled) {
      const mod = {};
      if (sfcModT !== "") mod.temperature = parseFloat(sfcModT);
      if (sfcModTd !== "") mod.dewpoint = parseFloat(sfcModTd);
      if (sfcModWspd !== "") mod.wind_speed = parseFloat(sfcModWspd);
      if (sfcModWdir !== "") mod.wind_direction = parseFloat(sfcModWdir);
      if (Object.keys(mod).length > 0) params.surfaceMod = mod;
    }

    // Custom storm motion
    if (smEnabled && smDirection !== "" && smSpeed !== "") {
      params.stormMotion = {
        direction: parseFloat(smDirection),
        speed: parseFloat(smSpeed),
      };
    }

    // VAD Wind Profile overlay
    if (vadEnabled) {
      // Determine nearest NEXRAD radar for the selected station or lat/lon
      let vadLat = parseFloat(lat), vadLon = parseFloat(lon);
      if ((!vadLat || !vadLon) && station) {
        const stn = stations.find((s) => s.id === station);
        if (stn) { vadLat = stn.lat; vadLon = stn.lon; }
      }
      if (vadLat && vadLon) {
        params.vad = nearestNexradForVad(vadLat, vadLon);
      }
    }

    // Storm-relative hodograph
    if (srHodoEnabled) {
      params.srHodograph = true;
    }

    // Profile smoothing
    if (smoothEnabled) {
      const sigma = parseFloat(smoothSigma);
      if (sigma > 0) params.smoothing = sigma;
    }

    // Boundary line
    if (boundaryEnabled && boundaryOrientation !== "") {
      const deg = parseFloat(boundaryOrientation);
      if (!isNaN(deg)) params.boundaryOrientation = deg;
    }

    onSubmit(params);
  };

  const handleStationSelect = (id, autoFetch = false) => {
    setStation(id);
    setStationSearch("");
    setPointMode(false);  // Exit point mode when selecting a station
    const stn = stations.find((s) => s.id === id);
    if (stn) {
      setLat(String(stn.lat));
      setLon(String(stn.lon));
    }
    if (autoFetch && !loading) {
      let params;
      if (scanMode === "forecast" && riskData?.model) {
        // Forecast risk scan — fetch forecast sounding with the scanned model/fhour
        params = {
          source: "psu",
          station: id,
          model: riskData.model.toLowerCase(),
          fhour: riskData.fhour || 0,
        };
      } else {
        params = { source, station: id };
        if (date) params.date = date.replace(/[-T:]/g, "").slice(0, 10);
        if (source === "bufkit" || source === "psu") {
          params.model = model;
          params.fhour = parseInt(fhour) || 0;
        }
      }
      onSubmit(params);
    }
  };

  if (initialLoading || connectError) {
    return (
      <aside className="control-panel">
        <div className="cp-brand">
          <svg width="28" height="28" viewBox="0 0 28 28" fill="none" xmlns="http://www.w3.org/2000/svg">
            <rect x="2" y="2" width="24" height="24" rx="4" fill="rgba(59,130,246,0.12)" stroke="#3b82f6" strokeWidth="1.5"/>
            <path d="M8 22 L10 18 L11 16 L12 13 L14 10 L16 8 L19 5" stroke="#ef4444" strokeWidth="1.8" strokeLinecap="round" strokeLinejoin="round" fill="none"/>
            <path d="M6 22 L8 19 L9 17 L9.5 15 L10 13 L10 11 L10.5 9 L11 6" stroke="#22c55e" strokeWidth="1.8" strokeLinecap="round" strokeLinejoin="round" fill="none"/>
            <line x1="21" y1="10" x2="21" y2="22" stroke="#60a5fa" strokeWidth="1.2" strokeLinecap="round"/>
            <line x1="21" y1="10" x2="24" y2="8" stroke="#60a5fa" strokeWidth="1.2" strokeLinecap="round"/>
            <line x1="21" y1="13" x2="24" y2="11" stroke="#60a5fa" strokeWidth="1.2" strokeLinecap="round"/>
          </svg>
          <div>
            <h1 className="cp-brand-title">Sounding Analysis</h1>
            <p className="cp-brand-sub">Atmospheric Profile Tool</p>
          </div>
        </div>
        <div className="cp-loading-state">
          {initialLoading ? (
            <>
              <Loader2 className="spin cp-loading-icon" size={32} />
              <span className="cp-loading-title">Connecting to API&hellip;</span>
              <span className="cp-loading-hint">May take a moment, be patient</span>
              <div className="cp-loading-skeleton">
                <div className="cp-skel-bar" style={{ width: "80%" }} />
                <div className="cp-skel-bar" style={{ width: "60%" }} />
                <div className="cp-skel-bar" style={{ width: "90%" }} />
                <div className="cp-skel-bar" style={{ width: "45%" }} />
              </div>
            </>
          ) : (
            <>
              <AlertTriangle size={32} style={{ color: "var(--danger, #e74c3c)" }} />
              <span className="cp-loading-title" style={{ color: "var(--danger, #e74c3c)" }}>Connection Failed</span>
              <span className="cp-loading-hint">{connectError}</span>
              <button type="button" className="cp-retry-btn" onClick={onRetry}>Retry</button>
            </>
          )}
        </div>
      </aside>
    );
  }

  const isFloat = typeof window !== "undefined" && document.body.classList.contains("smap-fullscreen-active");

  return (
    <aside
      className={`control-panel${floatCollapsed ? " cp-float-collapsed" : ""}`}
      style={isFloat ? { left: dragPos.x, top: dragPos.y, bottom: 'auto', margin: 0 } : undefined}
      ref={dragRef}
    >
      {/* Drag handle – only visible in float mode */}
      <div
        className="cp-drag-handle"
        {...dragHandlers}
      >
        <GripVertical size={14} />
        <span className="cp-drag-label">Sounding Analysis</span>
        <div className="cp-drag-actions">
          <button
            type="button"
            className="cp-drag-btn"
            onClick={toggleFloatCollapsed}
            title={floatCollapsed ? "Expand panel" : "Collapse panel"}
          >
            {floatCollapsed ? <ChevronDown size={14} /> : <ChevronUp size={14} />}
          </button>
        </div>
      </div>

      {/* Brand + Tab Navigation — compact single header */}
      <div className="cp-brand">
        <span className="cp-brand-title">Sounding Analysis</span>
        <div className="cp-brand-actions">
          <button
            type="button"
            className="cp-brand-icon-btn"
            onClick={onToggleTheme}
            title={theme === "dark" ? "Switch to light theme" : "Switch to dark theme"}
          >
            {theme === "dark" ? <Sun size={14} /> : <Moon size={14} />}
          </button>
        </div>
      </div>

      {/* Tab Navigation */}
      <nav className="cp-nav">
        <button
          type="button"
          className={`cp-nav-tab${navTab === "data" ? " active" : ""}`}
          onClick={() => setNavTab("data")}
        >
          <Database size={14} />
          <span>Data</span>
        </button>
        <button
          type="button"
          className={`cp-nav-tab${navTab === "modify" ? " active" : ""}`}
          onClick={() => setNavTab("modify")}
        >
          <Sliders size={14} />
          <span>Modify</span>
        </button>
        <button
          type="button"
          className={`cp-nav-tab${navTab === "tools" ? " active" : ""}`}
          onClick={() => setNavTab("tools")}
        >
          <Wrench size={14} />
          <span>Tools</span>
        </button>
      </nav>

      <div className="cp-tab-content">

      {/* ═══ DATA TAB ═══ */}
      {navTab === "data" && (
        <DataTab
          sources={sources}
          models={models}
          psuModels={psuModels}
          filteredStations={filteredStations}
          riskData={riskData}
          stations={stations}
          source={source} setSource={setSource}
          station={station}
          lat={lat} setLat={setLat}
          lon={lon} setLon={setLon}
          date={date} setDate={setDate}
          model={model} setModel={setModel}
          fhour={fhour} setFhour={setFhour}
          scanMode={scanMode} setScanMode={setScanMode}
          fcstModel={fcstModel} setFcstModel={setFcstModel}
          fcstFhour={fcstFhour} setFcstFhour={setFcstFhour}
          soundingHour={soundingHour} setSoundingHour={setSoundingHour}
          stationSearch={stationSearch} setStationSearch={setStationSearch}
          sortMode={sortMode} setSortMode={setSortMode}
          favorites={favorites}
          onToggleFavorite={(id) => {
            const { favorites: newFavs } = toggleFavorite(id);
            setFavorites(newFavs);
          }}
          pointMode={pointMode} setPointMode={setPointMode}
          needsStation={needsStation}
          needsModel={needsModel}
          scanning={scanning}
          loading={loading}
          listRef={listRef}
          onSubmit={handleSubmit}
          onStationSelect={handleStationSelect}
          onRiskScan={handleRiskScan}
        />
      )}

      {/* ═══ MODIFY TAB ═══ */}
      {navTab === "modify" && (() => {
        // Build the VAD button hint: append the nearest WSR-88D radar id
        // when we have a usable lat/lon (either from the form, or by
        // looking up the selected station).
        let vLat = parseFloat(lat), vLon = parseFloat(lon);
        if ((!vLat || !vLon) && station) {
          const stn = stations.find((s) => s.id === station);
          if (stn) { vLat = stn.lat; vLon = stn.lon; }
        }
        const vadSuffix = (vLat && vLon) ? ` (${nearestNexradForVad(vLat, vLon)})` : "";
        const vadTitleHint = `Overlay NEXRAD VAD winds on the hodograph from the nearest WSR-88D radar${vadSuffix}`;
        return (
          <ModifyTab
            sfcModEnabled={sfcModEnabled} setSfcModEnabled={setSfcModEnabled}
            sfcModT={sfcModT} setSfcModT={setSfcModT}
            sfcModTd={sfcModTd} setSfcModTd={setSfcModTd}
            sfcModWspd={sfcModWspd} setSfcModWspd={setSfcModWspd}
            sfcModWdir={sfcModWdir} setSfcModWdir={setSfcModWdir}
            smEnabled={smEnabled} setSmEnabled={setSmEnabled}
            smDirection={smDirection} setSmDirection={setSmDirection}
            smSpeed={smSpeed} setSmSpeed={setSmSpeed}
            vadEnabled={vadEnabled} setVadEnabled={setVadEnabled}
            vadTitleHint={vadTitleHint}
            srHodoEnabled={srHodoEnabled} setSrHodoEnabled={setSrHodoEnabled}
            boundaryEnabled={boundaryEnabled} setBoundaryEnabled={setBoundaryEnabled}
            boundaryOrientation={boundaryOrientation} setBoundaryOrientation={setBoundaryOrientation}
            smoothEnabled={smoothEnabled} setSmoothEnabled={setSmoothEnabled}
            smoothSigma={smoothSigma} setSmoothSigma={setSmoothSigma}
          />
        );
      })()}

      {/* ═══ TOOLS TAB ═══ */}
      {navTab === "tools" && (
        <ToolsTab
          showMap={showMap}
          showRisk={showRisk}
          showTimeSeries={showTimeSeries}
          showCompare={showCompare}
          showVwp={showVwp}
          showHistory={showHistory}
          riskData={riskData}
          colorblind={colorblind}
          onToggleMap={onToggleMap}
          onToggleRisk={onToggleRisk}
          onToggleTimeSeries={onToggleTimeSeries}
          onToggleCompare={onToggleCompare}
          onToggleVwp={onToggleVwp}
          onToggleHistory={onToggleHistory}
          onNavigateEnsemble={onNavigateEnsemble}
          onToggleColorblind={onToggleColorblind}
          onNavigateUpload={onNavigateUpload}
        />
      )}

      </div>{/* end cp-tab-content */}

      {/* Links — always visible at bottom */}
      <div className="cp-section-group cp-bottom-links">
        <span className="cp-group-label">Links</span>
        <div className="cp-footer">
          <button
            type="button"
            className={`cp-footer-btn ${feedbackActive ? "active" : ""}`}
            onClick={onFeedbackClick}
            title="Send feedback"
          >
            <MessageSquarePlus size={14} />
            <span>Feedback</span>
          </button>
          <a
            href="https://github.com/ShianMike/SoundingAnalysis"
            target="_blank"
            rel="noopener noreferrer"
            className="cp-footer-btn"
            title="View on GitHub"
          >
            <Github size={14} />
            <span>GitHub</span>
          </a>
          {onShowShortcuts && (
            <button
              type="button"
              className="cp-footer-btn"
              onClick={onShowShortcuts}
              title="Keyboard shortcuts (?)"
            >
              <Keyboard size={14} />
              <span>Shortcuts</span>
            </button>
          )}
        </div>
      </div>
    </aside>
  );
}
