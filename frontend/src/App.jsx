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
  const [lastParams, setLastParams] = useState(null);

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

  return (
    <div className="app">
      <Header />
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
        />
        {showHistory && (
          <HistoryPanel
            onLoad={handleLoadHistory}
            onClose={() => setShowHistory(false)}
          />
        )}
        <ResultsView result={result} loading={loading} error={error} riskData={riskData} />
      </main>
    </div>
  );
}
