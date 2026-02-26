import { useState, useEffect, useCallback } from "react";
import Header from "./components/Header";
import ControlPanel from "./components/ControlPanel";
import ResultsView from "./components/ResultsView";
import { fetchStations, fetchSources, fetchSounding } from "./api";
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
    try {
      const data = await fetchSounding(params);
      setResult(data);
    } catch (e) {
      setError(e.message);
      setResult(null);
    } finally {
      setLoading(false);
    }
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
        />
        <ResultsView result={result} loading={loading} error={error} riskData={riskData} />
      </main>
    </div>
  );
}
