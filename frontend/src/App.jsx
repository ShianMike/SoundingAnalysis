import { useState, useEffect, useCallback, useRef, lazy, Suspense } from "react";
import Header from "./components/Header";
import ControlPanel from "./components/ControlPanel";
import ResultsView from "./components/ResultsView";

const HistoryPanel   = lazy(() => import("./components/HistoryPanel"));
const CustomUpload   = lazy(() => import("./components/CustomUpload"));
const EnsemblePlume  = lazy(() => import("./components/EnsemblePlume"));
import { fetchStations, fetchSources, fetchSounding } from "./api";
import { saveToHistory } from "./history";
import "./App.css";

/* ── URL ↔ params helpers ────────────────────────────────── */
function parseUrlParams() {
  const sp = new URLSearchParams(window.location.search);
  if (sp.size === 0) return null;
  const p = {};
  if (sp.has("source")) p.source = sp.get("source");
  if (sp.has("station")) p.station = sp.get("station").toUpperCase();
  if (sp.has("date")) p.date = sp.get("date");
  if (sp.has("lat")) p.lat = parseFloat(sp.get("lat"));
  if (sp.has("lon")) p.lon = parseFloat(sp.get("lon"));
  if (sp.has("model")) p.model = sp.get("model");
  if (sp.has("fhour")) p.fhour = parseInt(sp.get("fhour"), 10);
  return Object.keys(p).length ? p : null;
}

function updateUrl(params) {
  const sp = new URLSearchParams();
  if (params.source) sp.set("source", params.source);
  if (params.station) sp.set("station", params.station);
  if (params.date) sp.set("date", params.date);
  if (params.lat != null) sp.set("lat", String(params.lat));
  if (params.lon != null) sp.set("lon", String(params.lon));
  if (params.model) sp.set("model", params.model);
  if (params.fhour != null) sp.set("fhour", String(params.fhour));
  const qs = sp.toString();
  const newUrl = qs ? `${window.location.pathname}?${qs}` : window.location.pathname;
  window.history.replaceState(null, "", newUrl);
}

/* ── Persisted preferences (theme / colorblind) ──────────── */
function loadPref(key, fallback) {
  try { return localStorage.getItem(key) || fallback; } catch { return fallback; }
}
function savePref(key, val) {
  try { localStorage.setItem(key, val); } catch { /* noop */ }
}

export default function App() {
  const [stations, setStations] = useState([]);
  const [sources, setSources] = useState([]);
  const [models, setModels] = useState([]);
  const [psuModels, setPsuModels] = useState([]);
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState(null);
  const [error, setError] = useState(null);
  const [initialLoading, setInitialLoading] = useState(true);
  const [riskData, setRiskData] = useState(null);
  const [showHistory, setShowHistory] = useState(false);
  const [showMap, setShowMap] = useState(true);
  const [showRisk, setShowRisk] = useState(false);
  const [showTimeSeries, setShowTimeSeries] = useState(false);
  const [showCompare, setShowCompare] = useState(false);
  const [showVwp, setShowVwp] = useState(false);

  const [compareHistoryData, setCompareHistoryData] = useState(null);
  const [lastParams, setLastParams] = useState(null);
  const [selectedStation, setSelectedStation] = useState("OUN");
  const [source, setSource] = useState("obs");
  const [showFeedback, setShowFeedback] = useState(false);

  // ── Page routing (main / upload) ──────────────────────────
  const [page, setPage] = useState("main");

  // ── Theme & accessibility prefs ───────────────────────────
  const [theme, setThemeState] = useState(() => loadPref("sa_theme", "dark"));
  const [colorblind, setColorblindState] = useState(() => loadPref("sa_cb", "false") === "true");

  // Apply to <html> on mount & change
  useEffect(() => {
    document.documentElement.setAttribute("data-theme", theme);
    savePref("sa_theme", theme);
  }, [theme]);
  useEffect(() => {
    document.documentElement.setAttribute("data-cb", String(colorblind));
    savePref("sa_cb", String(colorblind));
  }, [colorblind]);

  const toggleTheme = () => setThemeState((t) => (t === "dark" ? "light" : "dark"));
  const toggleColorblind = () => setColorblindState((v) => !v);

  /** Scroll a section into view after React re-renders */
  const scrollToSection = (id) => {
    requestAnimationFrame(() => {
      setTimeout(() => {
        const el = document.getElementById(id);
        if (el) {
          const y = el.getBoundingClientRect().top + window.scrollY - 20;
          window.scrollTo({ top: Math.max(0, y), behavior: "smooth" });
        }
      }, 120);
    });
  };

  // URL-based initial params (parsed once on mount)
  const urlParamsRef = useRef(parseUrlParams());

  const loadInitialData = useCallback(() => {
    setInitialLoading(true);
    setError(null);
    Promise.all([fetchStations(), fetchSources()])
      .then(([stationsData, sourcesData]) => {
        setStations(stationsData);
        setSources(sourcesData.sources);
        setModels(sourcesData.models);
        setPsuModels(sourcesData.psuModels || []);
      })
      .catch(() => setError("Failed to connect to API. Is the backend running?"))
      .finally(() => setInitialLoading(false));
  }, []);

  useEffect(() => {
    loadInitialData();
  }, [loadInitialData]);

  // Auto-fetch from URL params (or latest for default station) once stations are loaded
  const autoFetchedRef = useRef(false);
  useEffect(() => {
    if (autoFetchedRef.current || initialLoading || stations.length === 0) return;
    autoFetchedRef.current = true;
    const up = urlParamsRef.current;
    if (up) {
      // Sync parent state from URL params
      if (up.source) setSource(up.source);
      if (up.station) setSelectedStation(up.station);
      handleSubmit(up);
    } else {
      // No URL params — auto-fetch latest sounding for default station
      handleSubmit({ source, station: selectedStation, date: "", soundingHour: "latest" });
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [initialLoading, stations]);

  const handleSubmit = async (params) => {
    setLoading(true);
    setError(null);
    setLastParams(params);
    updateUrl(params);
    scrollToSection("section-sounding");
    try {
      const data = await fetchSounding({ ...params, theme, colorblind });
      setResult(data);
      saveToHistory(params, data);
      scrollToSection("section-sounding");
    } catch (e) {
      setError(e.message);
      setResult(null);
    } finally {
      setLoading(false);
    }
  };

  const handleLoadHistory = (historyResult) => {
    setResult(historyResult);
    setError(null);
    setShowHistory(false);
  };

  const handleLoadCompareHistory = (data) => {
    // data = { slots, results } — open Compare view and pre-load results
    setShowCompare(true);
    setShowHistory(false);
    // Store in a ref so ComparisonView can pick it up
    setCompareHistoryData(data);
  };

  const handleCompareStations = (stationIds) => {
    // Open Compare view with pre-filled station slots from proximity search
    const slots = stationIds.map((id) => ({ station: id, source: "obs", date: "", hour: "12" }));
    setShowCompare(true);
    setCompareHistoryData({ slots, results: [] });
    scrollToSection("section-compare");
  };

  const handleMapStationSelect = (stationId) => {
    setSelectedStation(stationId);
  };

  const handleTimelineSelect = (dateKey) => {
    // Re-fetch sounding for the same station with a different date
    handleSubmit({
      source: source || "obs",
      station: selectedStation,
      date: dateKey,
    });
  };

  const handleRiskStationSelect = (stationId, riskMeta) => {
    // Load sounding for a station from risk scan results
    if (riskMeta?.model) {
      // Forecast scan — load PSU sounding (backend falls back to BUFKIT)
      handleSubmit({
        source: "psu",
        station: stationId,
        model: riskMeta.model.toLowerCase(),
        fhour: riskMeta.fhour || 0,
      });
    } else {
      // Observed scan — load observed sounding
      handleSubmit({
        source: "obs",
        station: stationId,
        date: riskMeta?.date?.replace(/[^0-9]/g, "").slice(0, 10) || "",
      });
    }
    setSelectedStation(stationId);
  };

  const handleMapLatLonSelect = (lat, lon) => {
    // Store lat/lon for ControlPanel to pick up — point sounding mode
    setLastParams((prev) => ({ ...prev, _mapLat: lat, _mapLon: lon }));
  };

  const handleMapStormMotionSelect = ({ direction, speed }) => {
    // stored for ControlPanel to pick up via props
    setLastParams((prev) => ({
      ...prev,
      _mapStormDirection: direction,
      _mapStormSpeed: speed,
    }));
  };

  const handleMapBoundaryOrientationSelect = (orientation) => {
    // stored for ControlPanel to pick up via props
    setLastParams((prev) => ({ ...prev, _mapBoundaryOrientation: orientation }));
  };

  const handleFetchLatest = (stationId) => {
    // Auto-select the best source: try PSU (latest model run) first,
    // the backend auto-falls back to BUFKIT / OBS if PSU fails
    setSelectedStation(stationId);
    handleSubmit({
      source: "psu",
      station: stationId,
      model: "hrrr",
      fhour: 0,
    });
  };

  const handleSourceChange = (src) => {
    setSource(src);
  };

  const handleStationChange = (id) => {
    setSelectedStation(id);
  };

  // ── Auto-refresh polling ──────────────────────────────────
  const [autoRefresh, setAutoRefresh] = useState(false);
  const [refreshInterval, setRefreshInterval] = useState(300000); // 5 min default
  const autoRefreshRef = useRef(null);

  useEffect(() => {
    clearInterval(autoRefreshRef.current);
    if (autoRefresh && lastParams) {
      autoRefreshRef.current = setInterval(() => {
        handleSubmit(lastParams);
      }, refreshInterval);
    }
    return () => clearInterval(autoRefreshRef.current);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [autoRefresh, refreshInterval, lastParams]);

  // ── Keyboard shortcuts ────────────────────────────────────
  useEffect(() => {
    const handler = (e) => {
      // Ignore when typing in inputs
      const tag = e.target.tagName;
      if (tag === "INPUT" || tag === "TEXTAREA" || tag === "SELECT") return;
      if (e.ctrlKey || e.metaKey || e.altKey) return;

      switch (e.key.toLowerCase()) {
        case "h": e.preventDefault(); setShowHistory((v) => !v); break;
        case "c": e.preventDefault(); setShowCompare((v) => !v); break;
        case "m": e.preventDefault(); setShowMap((v) => !v); break;
        case "t": e.preventDefault(); setShowTimeSeries((v) => !v); break;
        case "v": e.preventDefault(); setShowVwp((v) => !v); break;
        case "?": e.preventDefault(); setShowShortcuts((v) => !v); break;
        case "escape": setShowShortcuts(false); break;
        default: break;
      }
    };
    window.addEventListener("keydown", handler);
    return () => window.removeEventListener("keydown", handler);
  }, []);

  const [showShortcuts, setShowShortcuts] = useState(false);

  return (
    <div className="app">
      <Header showFeedback={showFeedback} onCloseFeedback={() => setShowFeedback(false)} />

      {/* Keyboard shortcuts modal */}
      {showShortcuts && (
        <div className="shortcuts-overlay" onClick={() => setShowShortcuts(false)}>
          <div className="shortcuts-modal" onClick={(e) => e.stopPropagation()}>
            <h3>Keyboard Shortcuts</h3>
            <div className="shortcuts-grid">
              <kbd>H</kbd><span>Toggle History panel</span>
              <kbd>C</kbd><span>Toggle Compare view</span>
              <kbd>M</kbd><span>Toggle Station Map</span>
              <kbd>T</kbd><span>Toggle Time-Series</span>
              <kbd>V</kbd><span>Toggle VWP Display</span>
              <kbd>?</kbd><span>Show this help</span>
              <kbd>Esc</kbd><span>Close dialogs</span>
            </div>
            <button className="shortcuts-close" onClick={() => setShowShortcuts(false)}>Close</button>
          </div>
        </div>
      )}
      {page === "upload" ? (
        <Suspense fallback={<div className="loading-placeholder">Loading…</div>}>
          <CustomUpload
            onBack={() => setPage("main")}
            theme={theme}
            colorblind={colorblind}
          />
        </Suspense>
      ) : page === "ensemble" ? (
        <Suspense fallback={<div className="loading-placeholder">Loading…</div>}>
          <EnsemblePlume
            station={selectedStation}
            stations={stations}
            onBack={() => setPage("main")}
            theme={theme}
            colorblind={colorblind}
          />
        </Suspense>
      ) : (
        <main className="app-main">
          <ControlPanel
            stations={stations}
            sources={sources}
            models={models}
            psuModels={psuModels}
            onSubmit={handleSubmit}
            loading={loading}
            initialLoading={initialLoading}
            onRetry={loadInitialData}
            connectError={initialLoading ? null : (stations.length === 0 ? error : null)}
            riskData={riskData}
            onRiskDataChange={(data) => { setRiskData(data); if (data) { setShowRisk(true); scrollToSection("section-risk"); } }}
            showRisk={showRisk}
            onToggleRisk={() => setShowRisk((v) => { if (!v) scrollToSection("section-risk"); return !v; })}
            showHistory={showHistory}
            onToggleHistory={() => setShowHistory((v) => !v)}
            showMap={showMap}
            onToggleMap={() => setShowMap((v) => { if (!v) scrollToSection("section-map"); return !v; })}
            showTimeSeries={showTimeSeries}
            onToggleTimeSeries={() => setShowTimeSeries((v) => { if (!v) scrollToSection("section-timeseries"); return !v; })}
            showCompare={showCompare}
            onToggleCompare={() => setShowCompare((v) => { if (!v) scrollToSection("section-compare"); return !v; })}
            showVwp={showVwp}
            onToggleVwp={() => setShowVwp((v) => { if (!v) scrollToSection("section-vwp"); return !v; })}
            onNavigateEnsemble={() => setPage("ensemble")}
            selectedStation={selectedStation}
            onStationChange={handleStationChange}
            onSourceChange={handleSourceChange}
            mapLatLon={lastParams?._mapLat != null ? { lat: lastParams._mapLat, lon: lastParams._mapLon } : null}
            mapStormMotion={lastParams?._mapStormDirection != null ? {
              direction: lastParams._mapStormDirection,
              speed: lastParams._mapStormSpeed,
            } : null}
            mapBoundaryOrientation={lastParams?._mapBoundaryOrientation != null ? lastParams._mapBoundaryOrientation : null}
            onFeedbackClick={() => setShowFeedback((v) => !v)}
            showFeedback={showFeedback}
            urlParams={urlParamsRef.current}
            theme={theme}
            onToggleTheme={toggleTheme}
            colorblind={colorblind}
            onToggleColorblind={toggleColorblind}
            onNavigateUpload={() => setPage("upload")}
            onShowShortcuts={() => setShowShortcuts(true)}
          />
          {showHistory && (
            <Suspense fallback={null}>
              <HistoryPanel
                onLoad={handleLoadHistory}
                onLoadCompare={handleLoadCompareHistory}
                onClose={() => setShowHistory(false)}
              />
            </Suspense>
          )}

          <ResultsView
            result={result}
            loading={loading}
            error={error}
            riskData={riskData}
            showRisk={showRisk}
            showMap={showMap}
            showTimeSeries={showTimeSeries}
            onCloseTimeSeries={() => setShowTimeSeries(false)}
            showCompare={showCompare}
            onCloseCompare={() => setShowCompare(false)}
            showVwp={showVwp}
            onCloseVwp={() => setShowVwp(false)}

            compareHistoryData={compareHistoryData}
            onCompareHistoryConsumed={() => setCompareHistoryData(null)}
            stations={stations}
            selectedStation={selectedStation}
            source={source}
            lastParams={lastParams}
            autoRefresh={autoRefresh}
            onToggleAutoRefresh={() => setAutoRefresh((v) => !v)}
            refreshInterval={refreshInterval}
            onRefreshIntervalChange={setRefreshInterval}
            theme={theme}
            onTimelineSelect={handleTimelineSelect}
            onRiskStationSelect={handleRiskStationSelect}
            mapProps={{
              stations,
              riskData,
              selectedStation,
              onStationSelect: handleMapStationSelect,
              onLatLonSelect: handleMapLatLonSelect,
              onStormMotionSelect: handleMapStormMotionSelect,
              onBoundaryOrientationSelect: handleMapBoundaryOrientationSelect,
              onClose: () => setShowMap(false),
              onFetchLatest: handleFetchLatest,
              onCompareStations: handleCompareStations,

            }}
          />
        </main>
      )}
    </div>
  );
}
