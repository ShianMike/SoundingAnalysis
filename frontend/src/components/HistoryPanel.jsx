import { useState } from "react";
import { Clock, Trash2, X, RotateCcw } from "lucide-react";
import { getHistory, loadFromHistory, deleteFromHistory, clearHistory } from "../history";
import "./HistoryPanel.css";

export default function HistoryPanel({ onLoad, onClose }) {
  const [entries, setEntries] = useState(() => getHistory());

  const handleLoad = (id) => {
    const result = loadFromHistory(id);
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

  return (
    <div className="history-panel">
      <div className="hp-header">
        <div className="hp-header-left">
          <Clock size={14} />
          <h3>History</h3>
          <span className="hp-count">{entries.length}</span>
        </div>
        <div className="hp-header-right">
          {entries.length > 0 && (
            <button className="hp-clear-btn" onClick={handleClear} title="Clear all">
              <Trash2 size={12} />
              Clear
            </button>
          )}
          <button className="hp-close-btn" onClick={onClose} title="Close">
            <X size={14} />
          </button>
        </div>
      </div>

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
    </div>
  );
}
