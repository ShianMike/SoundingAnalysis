import {
  Thermometer,
  ChevronRight,
  Wind,
  RotateCcw,
  Layers,
  Crosshair,
  Minus,
  Waves,
} from "lucide-react";

/**
 * `ModifyTab` — sidebar tab containing the modification accordions:
 * surface mod, custom storm motion, VAD overlay, SR hodograph,
 * boundary line, profile smoothing.
 *
 * Pure presentational. All state + setters come from props (so the parent
 * can decide where to keep them — currently in `ControlPanel.jsx`).
 */
export default function ModifyTab({
  // Surface mod
  sfcModEnabled, setSfcModEnabled,
  sfcModT, setSfcModT,
  sfcModTd, setSfcModTd,
  sfcModWspd, setSfcModWspd,
  sfcModWdir, setSfcModWdir,
  // Custom storm motion
  smEnabled, setSmEnabled,
  smDirection, setSmDirection,
  smSpeed, setSmSpeed,
  // VAD overlay
  vadEnabled, setVadEnabled,
  vadTitleHint, // pre-built title hint built by parent (knows nearestNexrad)
  // SR hodograph
  srHodoEnabled, setSrHodoEnabled,
  // Boundary
  boundaryEnabled, setBoundaryEnabled,
  boundaryOrientation, setBoundaryOrientation,
  // Smoothing
  smoothEnabled, setSmoothEnabled,
  smoothSigma, setSmoothSigma,
}) {
  return (
    <div className="cp-form">
      <div className="cp-section-group">
        <span className="cp-group-label">Modifications</span>

        {/* Surface Modification */}
        <div className={`cp-accordion ${sfcModEnabled ? "cp-accordion--active" : ""}`}>
          <button
            type="button"
            className="cp-accordion-header"
            onClick={() => setSfcModEnabled((v) => !v)}
          >
            <div className="cp-accordion-left">
              <Thermometer size={14} className="cp-accordion-icon" />
              <span className="cp-accordion-title">Surface Modification</span>
            </div>
            <div className="cp-accordion-right">
              <span className={`cp-toggle-chip ${sfcModEnabled ? "on" : ""}`}>
                {sfcModEnabled ? "ON" : "OFF"}
              </span>
              <ChevronRight size={14} className={`cp-accordion-chevron ${sfcModEnabled ? "open" : ""}`} />
            </div>
          </button>
          <div className={`cp-accordion-body ${sfcModEnabled ? "expanded" : ""}`}>
            <div className="cp-accordion-content">
              <div className="cp-input-row">
                <div className="cp-input-group">
                  <label className="cp-input-group-label">Temperature</label>
                  <div className="cp-input-with-unit">
                    <input type="number" step="0.1" className="cp-input cp-input-sm" placeholder="—" value={sfcModT} onChange={(e) => setSfcModT(e.target.value)} />
                    <span className="cp-unit-badge">°C</span>
                  </div>
                </div>
                <div className="cp-input-group">
                  <label className="cp-input-group-label">Dewpoint</label>
                  <div className="cp-input-with-unit">
                    <input type="number" step="0.1" className="cp-input cp-input-sm" placeholder="—" value={sfcModTd} onChange={(e) => setSfcModTd(e.target.value)} />
                    <span className="cp-unit-badge">°C</span>
                  </div>
                </div>
              </div>
              <div className="cp-input-row">
                <div className="cp-input-group">
                  <label className="cp-input-group-label">Wind Speed</label>
                  <div className="cp-input-with-unit">
                    <input type="number" step="1" className="cp-input cp-input-sm" placeholder="—" value={sfcModWspd} onChange={(e) => setSfcModWspd(e.target.value)} />
                    <span className="cp-unit-badge">kt</span>
                  </div>
                </div>
                <div className="cp-input-group">
                  <label className="cp-input-group-label">Wind Dir</label>
                  <div className="cp-input-with-unit">
                    <input type="number" step="1" min="0" max="360" className="cp-input cp-input-sm" placeholder="—" value={sfcModWdir} onChange={(e) => setSfcModWdir(e.target.value)} />
                    <span className="cp-unit-badge">°</span>
                  </div>
                </div>
              </div>
              <button
                type="button"
                className="cp-reset-btn"
                onClick={() => { setSfcModT(""); setSfcModTd(""); setSfcModWspd(""); setSfcModWdir(""); }}
              >
                <RotateCcw size={11} /> Reset values
              </button>
            </div>
          </div>
        </div>

        {/* Custom Storm Motion */}
        <div className={`cp-accordion ${smEnabled ? "cp-accordion--active" : ""}`}>
          <button
            type="button"
            className="cp-accordion-header"
            onClick={() => setSmEnabled((v) => !v)}
          >
            <div className="cp-accordion-left">
              <Wind size={14} className="cp-accordion-icon" />
              <span className="cp-accordion-title">Custom Storm Motion</span>
            </div>
            <div className="cp-accordion-right">
              <span className={`cp-toggle-chip ${smEnabled ? "on" : ""}`}>
                {smEnabled ? "ON" : "OFF"}
              </span>
              <ChevronRight size={14} className={`cp-accordion-chevron ${smEnabled ? "open" : ""}`} />
            </div>
          </button>
          <div className={`cp-accordion-body ${smEnabled ? "expanded" : ""}`}>
            <div className="cp-accordion-content">
              <div className="cp-input-row">
                <div className="cp-input-group">
                  <label className="cp-input-group-label">Direction</label>
                  <div className="cp-input-with-unit">
                    <input type="number" step="1" min="0" max="360" className="cp-input cp-input-sm" placeholder="—" value={smDirection} onChange={(e) => setSmDirection(e.target.value)} />
                    <span className="cp-unit-badge">°</span>
                  </div>
                </div>
                <div className="cp-input-group">
                  <label className="cp-input-group-label">Speed</label>
                  <div className="cp-input-with-unit">
                    <input type="number" step="1" className="cp-input cp-input-sm" placeholder="—" value={smSpeed} onChange={(e) => setSmSpeed(e.target.value)} />
                    <span className="cp-unit-badge">kt</span>
                  </div>
                </div>
              </div>
              <button
                type="button"
                className="cp-reset-btn"
                onClick={() => { setSmDirection(""); setSmSpeed(""); }}
              >
                <RotateCcw size={11} /> Reset values
              </button>
            </div>
          </div>
        </div>

        {/* VAD Wind Profile Overlay */}
        <button
          type="button"
          className={`cp-toggle-btn ${vadEnabled ? "cp-toggle-btn--active" : ""}`}
          onClick={() => setVadEnabled((v) => !v)}
          title={vadTitleHint}
        >
          <Layers size={14} />
          <span>VAD Wind Profile</span>
          <span className={`cp-toggle-chip ${vadEnabled ? "on" : ""}`}>
            {vadEnabled ? "ON" : "OFF"}
          </span>
        </button>

        {/* Storm-Relative Hodograph */}
        <button
          type="button"
          className={`cp-toggle-btn ${srHodoEnabled ? "cp-toggle-btn--active" : ""}`}
          onClick={() => setSrHodoEnabled((v) => !v)}
          title="Plot hodograph in storm-relative frame — subtracts Bunkers RM (or custom SM) from all winds so storm motion is at the origin"
        >
          <Crosshair size={14} />
          <span>SR Hodograph</span>
          <span className={`cp-toggle-chip ${srHodoEnabled ? "on" : ""}`}>
            {srHodoEnabled ? "ON" : "OFF"}
          </span>
        </button>

        {/* Boundary Line */}
        <div className="cp-accordion-section">
          <button
            type="button"
            className={`cp-toggle-btn ${boundaryEnabled ? "cp-toggle-btn--active" : ""}`}
            onClick={() => setBoundaryEnabled((v) => !v)}
            title="Draw a boundary orientation line on the hodograph — represents an outflow boundary, front, or dryline orientation"
          >
            <Minus size={14} />
            <span>Boundary Line</span>
            <span className={`cp-toggle-chip ${boundaryEnabled ? "on" : ""}`}>
              {boundaryEnabled ? "ON" : "OFF"}
            </span>
          </button>
          {boundaryEnabled && (
            <div style={{ padding: "6px 10px" }}>
              <div className="cp-input-row">
                <div className="cp-input-group">
                  <label className="cp-input-group-label">Orientation</label>
                  <div className="cp-input-with-unit">
                    <input
                      type="number"
                      step="5"
                      min="0"
                      max="360"
                      className="cp-input cp-input-sm"
                      placeholder="e.g. 210"
                      value={boundaryOrientation}
                      onChange={(e) => setBoundaryOrientation(e.target.value)}
                    />
                    <span className="cp-unit-badge">°</span>
                  </div>
                </div>
              </div>
              <p style={{ margin: "4px 0 0", fontSize: 10, color: "var(--fg-faint, #707070)" }}>
                Direction the boundary runs (0–360°). e.g. 210° = SW–NE oriented.
              </p>
            </div>
          )}
        </div>

        {/* Profile Smoothing */}
        <div className="cp-accordion-section">
          <button
            type="button"
            className={`cp-toggle-btn ${smoothEnabled ? "cp-toggle-btn--active" : ""}`}
            onClick={() => setSmoothEnabled((v) => !v)}
            title="Apply Gaussian smoothing to T, Td, and wind profiles — reduces noise in model data"
          >
            <Waves size={14} />
            <span>Profile Smoothing</span>
            <span className={`cp-toggle-chip ${smoothEnabled ? "on" : ""}`}>
              {smoothEnabled ? "ON" : "OFF"}
            </span>
          </button>
          {smoothEnabled && (
            <div style={{ padding: "6px 10px" }}>
              <label className="cp-field-label" style={{ marginBottom: 0, display: "flex", alignItems: "center", gap: 8 }}>
                <span style={{ minWidth: 50 }}>σ = {smoothSigma}</span>
                <input
                  type="range"
                  min="1"
                  max="10"
                  step="0.5"
                  value={smoothSigma}
                  onChange={(e) => setSmoothSigma(e.target.value)}
                  style={{ flex: 1 }}
                  title="Gaussian sigma in data levels (higher = smoother). Typical: 2-5"
                />
              </label>
            </div>
          )}
        </div>

      </div>{/* end Modifications */}

      <p className="cp-hint" style={{ marginTop: 4 }}>
        These options modify the next sounding fetch. Go to Data tab and click Generate to apply.
      </p>
    </div>
  );
}
