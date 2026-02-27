import { useState, useEffect, useCallback } from "react";
import Header from "./components/Header";
import ControlPanel from "./components/ControlPanel";
import ResultsView from "./components/ResultsView";
import HistoryPanel from "./components/HistoryPanel";
import { fetchStations, fetchSources, fetchSounding } from "./api";
import { saveToHistory } from "./history";
import "./App.css";

export default function App() {
  const [stations, setStations] = useState([]);
  const [sources, setSources] = useState([]);
  const [models, setModels] = useState([]);
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState(null);
  const [error, setError] = useState(null);
  const [initialLoading, setInitialLoading] = useState(true);
  const [riskData, setRiskData] = useState(null);
  const [showHistory, setShowHistory] = useState(false);
  const [showMap, setShowMap] = useState(true);
  const [showTimeSeries, setShowTimeSeries] = useState(false);
  const [showCompare, setShowCompare] = useState(false);
  const [compareHistoryData, setCompareHistoryData] = useState(null);
  const [lastParams, setLastParams] = useState(null);
  const [selectedStation, setSelectedStation] = useState("OUN");
  const [source, setSource] = useState("obs");
  const [showFeedback, setShowFeedback] = useState(false);

  const loadInitialData = useCallback(() => {
    setInitialLoading(true);
    setError(null);
    Promise.all([fetchStations(), fetchSources()])
      .then(([stationsData, sourcesData]) => {
        setStations(stationsData);
        setSources(sourcesData.sources);
        setModels(sourcesData.models);
      })
      .catch(() => setError("Failed to connect to API. Is the backend running?"))
      .finally(() => setInitialLoading(false));
  }, []);

  useEffect(() => {
    loadInitialData();
  }, [loadInitialData]);

  const handleSubmit = async (params) => {
    setLoading(true);
    setError(null);
    setLastParams(params);
    try {
      const data = await fetchSounding(params);
      setResult(data);
      saveToHistory(params, data);
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

  const handleMapStationSelect = (stationId) => {
    setSelectedStation(stationId);
  };

  const handleMapLatLonSelect = (lat, lon) => {
    // stored for ControlPanel to pick up via props
    setLastParams((prev) => ({ ...prev, _mapLat: lat, _mapLon: lon }));
  };

  const handleSourceChange = (src) => {
    setSource(src);
  };

  const handleStationChange = (id) => {
    setSelectedStation(id);
  };

  return (
    <div className="app">
      <Header showFeedback={showFeedback} onCloseFeedback={() => setShowFeedback(false)} />
      <main className="app-main">
        <ControlPanel
          stations={stations}
          sources={sources}
          models={models}
          onSubmit={handleSubmit}
          loading={loading}
          initialLoading={initialLoading}
          onRetry={loadInitialData}
          connectError={initialLoading ? null : (stations.length === 0 ? error : null)}
          riskData={riskData}
          onRiskDataChange={setRiskData}
          showHistory={showHistory}
          onToggleHistory={() => setShowHistory((v) => !v)}
          showMap={showMap}
          onToggleMap={() => setShowMap((v) => !v)}
          showTimeSeries={showTimeSeries}
          onToggleTimeSeries={() => setShowTimeSeries((v) => !v)}
          showCompare={showCompare}
          onToggleCompare={() => setShowCompare((v) => !v)}
          selectedStation={selectedStation}
          onStationChange={handleStationChange}
          onSourceChange={handleSourceChange}
          mapLatLon={lastParams?._mapLat ? { lat: lastParams._mapLat, lon: lastParams._mapLon } : null}
          onFeedbackClick={() => setShowFeedback((v) => !v)}
          showFeedback={showFeedback}
        />
        {showHistory && (
          <HistoryPanel
            onLoad={handleLoadHistory}
            onLoadCompare={handleLoadCompareHistory}
            onClose={() => setShowHistory(false)}
          />
        )}
        <ResultsView
          result={result}
          loading={loading}
          error={error}
          riskData={riskData}
          showMap={showMap}
          showTimeSeries={showTimeSeries}
          onCloseTimeSeries={() => setShowTimeSeries(false)}
          showCompare={showCompare}
          onCloseCompare={() => setShowCompare(false)}
          compareHistoryData={compareHistoryData}
          onCompareHistoryConsumed={() => setCompareHistoryData(null)}
          stations={stations}
          selectedStation={selectedStation}
          source={source}
          mapProps={{
            stations,
            riskData,
            selectedStation,
            onStationSelect: handleMapStationSelect,
            onLatLonSelect: handleMapLatLonSelect,
            latLonMode: source === "rap",
            onClose: () => setShowMap(false),
          }}
        />
      </main>
    </div>
  );
}
