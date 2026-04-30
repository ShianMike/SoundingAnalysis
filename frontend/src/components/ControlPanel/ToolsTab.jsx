import {
  Map,
  Zap,
  TrendingUp,
  GitCompareArrows,
  Radio,
  History,
  Layers,
  Eye,
  Upload,
} from "lucide-react";

/**
 * `ToolsTab` — sidebar tab housing view toggles + cross-app navigation.
 *
 * Pure presentational. All state and handlers come from props; the parent
 * (ControlPanel) wires these to Zustand selectors and route navigation.
 */
export default function ToolsTab({
  // View toggles
  showMap,
  showRisk,
  showTimeSeries,
  showCompare,
  showVwp,
  showHistory,
  riskData,
  // Setting state
  colorblind,
  // Handlers
  onToggleMap,
  onToggleRisk,
  onToggleTimeSeries,
  onToggleCompare,
  onToggleVwp,
  onToggleHistory,
  onNavigateEnsemble,
  onToggleColorblind,
  onNavigateUpload,
}) {
  return (
    <div className="cp-form">
      <div className="cp-section-group">
        <span className="cp-group-label">Views</span>
        <div className="cp-tools-grid">
          <button
            type="button"
            className={`cp-tool-btn ${showMap ? "active" : ""}`}
            onClick={onToggleMap}
          >
            <Map size={13} />
            {showMap ? "Hide" : "Map"}
          </button>
          {riskData && (
            <button
              type="button"
              className={`cp-tool-btn ${showRisk ? "active" : ""}`}
              onClick={onToggleRisk}
            >
              <Zap size={13} />
              {showRisk ? "Hide" : "Risk"}
            </button>
          )}
          <button
            type="button"
            className={`cp-tool-btn ${showTimeSeries ? "active" : ""}`}
            onClick={onToggleTimeSeries}
          >
            <TrendingUp size={13} />
            Trends
          </button>
          <button
            type="button"
            className={`cp-tool-btn ${showCompare ? "active" : ""}`}
            onClick={onToggleCompare}
          >
            <GitCompareArrows size={13} />
            Compare
          </button>
          <button
            type="button"
            className={`cp-tool-btn ${showVwp ? "active" : ""}`}
            onClick={onToggleVwp}
          >
            <Radio size={13} />
            VWP
          </button>
          <button
            type="button"
            className={`cp-tool-btn ${showHistory ? "active" : ""}`}
            onClick={onToggleHistory}
          >
            <History size={13} />
            History
          </button>
          <button
            type="button"
            className="cp-tool-btn"
            onClick={onNavigateEnsemble}
          >
            <Layers size={13} />
            Plume
          </button>
        </div>
      </div>

      <div className="cp-section-group">
        <span className="cp-group-label">Settings</span>
        <div className="cp-settings-grid">
          <button
            type="button"
            className={`cp-settings-btn ${colorblind ? "active" : ""}`}
            onClick={onToggleColorblind}
            title="Toggle color-blind safe palette"
          >
            <Eye size={14} />
            <span>CB Mode</span>
          </button>
          <button
            type="button"
            className="cp-settings-btn"
            onClick={onNavigateUpload}
            title="Upload custom sounding data"
          >
            <Upload size={14} />
            <span>Upload</span>
          </button>
          <a
            href="https://modelforecastpy.app"
            target="_blank"
            rel="noopener noreferrer"
            className="cp-settings-btn"
            title="Open gridded model forecast maps"
          >
            <Map size={14} />
            <span>Forecast</span>
          </a>
        </div>
      </div>
    </div>
  );
}
