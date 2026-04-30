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
import { useAppStore } from "./store/useAppStore";
import { useShallow } from "zustand/shallow";

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

export default function App() {
  // Only subscribe to the slice of store this component actually uses.
  // `useShallow` keeps the destructure shape but skips re-renders when
  // unrelated fields change.
  const {
    stations, setGlobalData,
    setLoading,
    setResult,
    error, setError,
    initialLoading, setInitialLoading,
    setRiskData,
    showHistory, toggleHistory, setShowHistory,
    toggleMap,
    toggleTimeSeries,
    toggleCompare, setShowCompare,
    toggleVwp,
    showFeedback, setShowFeedback,
    lastParams, setLastParams,
    selectedStation, setSelectedStation,
    source, setSource,
    page, setPage,
    theme,
    colorblind,
    setCompareHistoryData,
  } = useAppStore(useShallow((s) => ({
    stations: s.stations,
    setGlobalData: s.setGlobalData,
    setLoading: s.setLoading,
    setResult: s.setResult,
    error: s.error,
    setError: s.setError,
    initialLoading: s.initialLoading,
    setInitialLoading: s.setInitialLoading,
    setRiskData: s.setRiskData,
    showHistory: s.showHistory,
    toggleHistory: s.toggleHistory,
    setShowHistory: s.setShowHistory,
    toggleMap: s.toggleMap,
    toggleTimeSeries: s.toggleTimeSeries,
    toggleCompare: s.toggleCompare,
    setShowCompare: s.setShowCompare,
    toggleVwp: s.toggleVwp,
    showFeedback: s.showFeedback,
    setShowFeedback: s.setShowFeedback,
    lastParams: s.lastParams,
    setLastParams: s.setLastParams,
    selectedStation: s.selectedStation,
    setSelectedStation: s.setSelectedStation,
    source: s.source,
    setSource: s.setSource,
    page: s.page,
    setPage: s.setPage,
    theme: s.theme,
    colorblind: s.colorblind,
    setCompareHistoryData: s.setCompareHistoryData,
  })));

  // Apply to <html> on mount & change
  useEffect(() => {
    document.documentElement.setAttribute("data-theme", theme);
  }, [theme]);
  useEffect(() => {
    document.documentElement.setAttribute("data-cb", String(colorblind));
  }, [colorblind]);

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
        setGlobalData({
          stations: stationsData,
          sources: sourcesData.sources,
          models: sourcesData.models,
          psuModels: sourcesData.psuModels || []
        });
      })
      .catch(() => setError("Failed to connect to API. Is the backend running?"))
      .finally(() => setInitialLoading(false));
  }, [setError, setGlobalData, setInitialLoading]);

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
    setRiskData(riskMeta);
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
        case "h": e.preventDefault(); toggleHistory(); break;
        case "c": e.preventDefault(); toggleCompare(); break;
        case "m": e.preventDefault(); toggleMap(); break;
        case "t": e.preventDefault(); toggleTimeSeries(); break;
        case "v": e.preventDefault(); toggleVwp(); break;
        case "?": e.preventDefault(); setShowShortcuts((v) => !v); break;
        case "escape": setShowShortcuts(false); break;
        default: break;
      }
    };
    window.addEventListener("keydown", handler);
    return () => window.removeEventListener("keydown", handler);
  }, [toggleCompare, toggleHistory, toggleMap, toggleTimeSeries, toggleVwp]);

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
            onSubmit={handleSubmit}
            onRetry={loadInitialData}
            connectError={initialLoading ? null : (stations.length === 0 ? error : null)}
            onNavigateEnsemble={() => setPage("ensemble")}
            mapLatLon={lastParams?._mapLat != null ? { lat: lastParams._mapLat, lon: lastParams._mapLon } : null}
            mapStormMotion={lastParams?._mapStormDirection != null ? {
              direction: lastParams._mapStormDirection,
              speed: lastParams._mapStormSpeed,
            } : null}
            mapBoundaryOrientation={lastParams?._mapBoundaryOrientation != null ? lastParams._mapBoundaryOrientation : null}
            urlParams={urlParamsRef.current}
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
            autoRefresh={autoRefresh}
            onToggleAutoRefresh={() => setAutoRefresh((v) => !v)}
            refreshInterval={refreshInterval}
            onRefreshIntervalChange={setRefreshInterval}
            onTimelineSelect={handleTimelineSelect}
            onRiskStationSelect={handleRiskStationSelect}
            mapProps={{
              onLatLonSelect: handleMapLatLonSelect,
              onStormMotionSelect: handleMapStormMotionSelect,
              onBoundaryOrientationSelect: handleMapBoundaryOrientationSelect,
              onFetchLatest: handleFetchLatest,
              onCompareStations: handleCompareStations,
            }}
          />
        </main>
      )}
    </div>
  );
}
