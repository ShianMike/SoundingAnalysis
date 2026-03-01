import { useState } from "react";
import { Clock, Trash2, X, RotateCcw, GitCompareArrows } from "lucide-react";
import {
  getHistory, loadFromHistory, deleteFromHistory, clearHistory,
  getCompareHistory, loadCompareFromHistory, deleteCompareFromHistory, clearCompareHistory,
} from "../history";
import "./HistoryPanel.css";

export default function HistoryPanel({ onLoad, onLoadCompare, onClose }) {
  const [tab, setTab] = useState("soundings"); // "soundings" | "comparisons"
  const [entries, setEntries] = useState(() => getHistory());
  const [compareEntries, setCompareEntries] = useState(() => getCompareHistory());

  const handleLoad = async (id) => {
    const result = await loadFromHistory(id);
    if (result) {
      onLoad(result);
    }
  };

  const handleDelete = (e, id) => {
    e.stopPropagation();
    deleteFromHistory(id);
    setEntries(getHistory());
  };

  const handleClear = () => {
    clearHistory();
    setEntries([]);
  };

  const handleLoadCompare = async (id) => {
    const data = await loadCompareFromHistory(id);
    if (data && onLoadCompare) {
      onLoadCompare(data);
    }
  };

  const handleDeleteCompare = (e, id) => {
    e.stopPropagation();
    deleteCompareFromHistory(id);
    setCompareEntries(getCompareHistory());
  };

  const handleClearCompare = () => {
    clearCompareHistory();
    setCompareEntries([]);
  };

  const formatTime = (ts) => {
    const d = new Date(ts);
    const now = new Date();
    const diff = now - d;
    if (diff < 60_000) return "Just now";
    if (diff < 3_600_000) return `${Math.floor(diff / 60_000)}m ago`;
    if (diff < 86_400_000) return `${Math.floor(diff / 3_600_000)}h ago`;
    if (diff < 604_800_000) return `${Math.floor(diff / 86_400_000)}d ago`;
    return d.toLocaleDateString();
  };

  const activeEntries = tab === "soundings" ? entries : compareEntries;
  const activeClear = tab === "soundings" ? handleClear : handleClearCompare;

  return (
    <div className="history-panel">
      <div className="hp-header">
        <div className="hp-header-left">
          <Clock size={14} />
          <h3>History</h3>
          <span className="hp-count">{activeEntries.length}</span>
        </div>
        <div className="hp-header-right">
          {activeEntries.length > 0 && (
            <button className="hp-clear-btn" onClick={activeClear} title="Clear all">
              <Trash2 size={12} />
              Clear
            </button>
          )}
          <button className="hp-close-btn" onClick={onClose} title="Close">
            <X size={14} />
          </button>
        </div>
      </div>

      {/* Tabs */}
      <div className="hp-tabs">
        <button
          className={`hp-tab ${tab === "soundings" ? "active" : ""}`}
          onClick={() => setTab("soundings")}
        >
          Soundings
          {entries.length > 0 && <span className="hp-tab-count">{entries.length}</span>}
        </button>
        <button
          className={`hp-tab ${tab === "comparisons" ? "active" : ""}`}
          onClick={() => setTab("comparisons")}
        >
          <GitCompareArrows size={12} />
          Comparisons
          {compareEntries.length > 0 && <span className="hp-tab-count">{compareEntries.length}</span>}
        </button>
      </div>

      {tab === "soundings" && (
        <>
          {entries.length === 0 ? (
            <div className="hp-empty">
              <RotateCcw size={20} />
              <p>No soundings in history yet.</p>
              <p className="hp-empty-hint">Generated soundings will appear here automatically.</p>
            </div>
          ) : (
            <div className="hp-list">
              {entries.map((entry) => (
                <button
                  key={entry.id}
                  className="hp-item"
                  onClick={() => handleLoad(entry.id)}
                >
                  <div className="hp-item-main">
                    <span className="hp-item-station">
                      {entry.meta.station || entry.meta.source?.toUpperCase()}
                    </span>
                    <span className="hp-item-name">{entry.meta.stationName}</span>
                  </div>
                  <div className="hp-item-details">
                    <span className="hp-item-date">{entry.meta.date}</span>
                    <span className="hp-item-source">{entry.requestParams?.source?.toUpperCase()}</span>
                  </div>
                  <div className="hp-item-right">
                    <span className="hp-item-time">{formatTime(entry.timestamp)}</span>
                    <button
                      className="hp-item-delete"
                      onClick={(e) => handleDelete(e, entry.id)}
                      title="Remove"
                    >
                      <X size={12} />
                    </button>
                  </div>
                </button>
              ))}
            </div>
          )}
        </>
      )}

      {tab === "comparisons" && (
        <>
          {compareEntries.length === 0 ? (
            <div className="hp-empty">
              <GitCompareArrows size={20} />
              <p>No comparisons in history yet.</p>
              <p className="hp-empty-hint">Comparisons will appear here automatically.</p>
            </div>
          ) : (
            <div className="hp-list">
              {compareEntries.map((entry) => (
                <button
                  key={entry.id}
                  className="hp-item hp-item-compare"
                  onClick={() => handleLoadCompare(entry.id)}
                >
                  <div className="hp-item-main">
                    <GitCompareArrows size={12} className="hp-compare-icon" />
                    <span className="hp-item-station">
                      {entry.summary.map((s) => s.meta?.station || "?").join(" vs ")}
                    </span>
                  </div>
                  <div className="hp-item-details">
                    <span className="hp-item-date">
                      {entry.summary[0]?.meta?.date || ""}
                    </span>
                    <span className="hp-item-source">
                      {entry.slots.length} soundings
                    </span>
                  </div>
                  <div className="hp-item-right">
                    <span className="hp-item-time">{formatTime(entry.timestamp)}</span>
                    <button
                      className="hp-item-delete"
                      onClick={(e) => handleDeleteCompare(e, entry.id)}
                      title="Remove"
                    >
                      <X size={12} />
                    </button>
                  </div>
                </button>
              ))}
            </div>
          )}
        </>
      )}
    </div>
  );
}
