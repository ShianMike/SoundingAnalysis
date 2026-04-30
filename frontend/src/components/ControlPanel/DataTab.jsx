import {
  Database,
  MapPin,
  Calendar,
  Layers,
  Clock,
  Loader2,
  Search,
  Eye,
  TrendingUp,
  Crosshair,
  X,
} from "lucide-react";
import StationPicker from "./StationPicker";
import { SOURCE_META, MODEL_META, FHOUR_PRESETS } from "../../config/constants";

/**
 * `DataTab` — sidebar tab housing the main data-fetch form: source picker,
 * station picker (with risk-scan controls), date/time picker, model picker
 * for forecast sources, and the Generate-Sounding submit button.
 *
 * Pure presentational. Owns no state — everything is plumbed from the
 * parent `ControlPanel`. Extracted from `ControlPanel.jsx` in the Phase 2
 * refactor to keep the parent's render under ~500 lines.
 */
export default function DataTab({
  // Data the form needs
  sources,
  models,
  psuModels,
  filteredStations,
  riskData,
  stations,
  // Form state
  source, setSource,
  station,
  lat, setLat,
  lon, setLon,
  date, setDate,
  model, setModel,
  fhour, setFhour,
  scanMode, setScanMode,
  fcstModel, setFcstModel,
  fcstFhour, setFcstFhour,
  soundingHour, setSoundingHour,
  stationSearch, setStationSearch,
  sortMode, setSortMode,
  favorites, onToggleFavorite,
  pointMode, setPointMode,
  // Flags
  needsStation,
  needsModel,
  scanning,
  loading,
  // Refs / handlers
  listRef,
  onSubmit,
  onStationSelect,
  onRiskScan,
}) {
  return (
    <form onSubmit={onSubmit} className="cp-form">
      {/* Source */}
      <div className="cp-section">
        <label className="cp-label">
          <Database size={14} />
          Data Source
        </label>
        <div className="cp-source-grid">
          {sources.map((s) => {
            const meta = SOURCE_META[s.id];
            return (
              <div key={s.id} className="cp-source-btn-wrap">
                <button
                  type="button"
                  className={`cp-source-btn ${source === s.id ? "active" : ""}`}
                  onClick={() => setSource(s.id)}
                  title={meta ? `${meta.label}\n${meta.desc}` : s.id.toUpperCase()}
                >
                  <span className="cp-source-id">{s.id.toUpperCase()}</span>
                </button>
              </div>
            );
          })}
        </div>
      </div>

      {/* Station */}
      {needsStation && (
        <div className="cp-section">
          <label className="cp-label">
            <MapPin size={14} />
            Station
          </label>
          <div className="cp-station-picker">
            {/* Scan mode toggle */}
            <div className="cp-scan-tabs">
              <button
                type="button"
                className={`cp-scan-tab${scanMode === "obs" ? " active" : ""}`}
                onClick={() => setScanMode("obs")}
              >
                <Eye size={12} /> Observed
              </button>
              <button
                type="button"
                className={`cp-scan-tab${scanMode === "forecast" ? " active" : ""}`}
                onClick={() => setScanMode("forecast")}
              >
                <TrendingUp size={12} /> Forecast
              </button>
            </div>

            {/* Forecast model + fhour pickers */}
            {scanMode === "forecast" && (() => {
              const meta = MODEL_META[fcstModel] || { maxF: 48, step: 1, lag: 3, interval: 1 };
              const fh = parseInt(fcstFhour) || 0;
              /* Mirror backend cycle logic: candidate = now - lag, snap to interval */
              const lag = meta.lag || 3;
              const interval = meta.interval || 1;
              const candidate = new Date(Date.now() - lag * 3600000);
              let initHour = candidate.getUTCHours();
              if (interval > 1) initHour = Math.floor(initHour / interval) * interval;
              const initDate = new Date(Date.UTC(candidate.getUTCFullYear(), candidate.getUTCMonth(), candidate.getUTCDate(), initHour));
              const validDate = new Date(initDate.getTime() + fh * 3600000);
              const months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
              const zLabel = `${validDate.getUTCDate()} ${months[validDate.getUTCMonth()]} ${String(validDate.getUTCHours()).padStart(2,"0")}Z`;
              return (
                <div className="cp-scan-forecast-opts">
                  <div className="cp-scan-row">
                    <select
                      className="cp-input cp-scan-select"
                      value={fcstModel}
                      onChange={(e) => {
                        setFcstModel(e.target.value);
                        const m = MODEL_META[e.target.value] || { maxF: 384, step: 1 };
                        if (parseInt(fcstFhour) > m.maxF) setFcstFhour(String(m.maxF));
                      }}
                    >
                      {["hrrr", "rap", "nam", "namnest", "gfs"].map((m) => (
                        <option key={m} value={m}>
                          {MODEL_META[m]?.short || m} (0–{MODEL_META[m]?.maxF}h)
                        </option>
                      ))}
                    </select>
                    <span className="cp-scan-valid-tag">F{fcstFhour} · {zLabel}</span>
                  </div>
                  <div className="cp-scan-slider-row">
                    <span className="cp-scan-edge-label">0h</span>
                    <input
                      type="range"
                      className="cp-scan-slider"
                      min={0}
                      max={meta.maxF}
                      step={meta.step}
                      value={fcstFhour}
                      onChange={(e) => setFcstFhour(e.target.value)}
                    />
                    <span className="cp-scan-edge-label">{meta.maxF}h</span>
                  </div>
                </div>
              );
            })()}

            <button
              type="button"
              className="cp-risk-btn"
              onClick={onRiskScan}
              disabled={scanning}
            >
              {scanning ? (
                <>
                  <Loader2 size={14} className="spin" />
                  Analyzing...
                </>
              ) : (
                <>
                  <Search size={14} />
                  {riskData ? "Rescan Severe Risk" : "Scan Severe Risk"}
                </>
              )}
            </button>
            {riskData && (
              <p className="cp-risk-hint">
                {riskData.model
                  ? `${(MODEL_META[riskData.model]?.short || riskData.model).toUpperCase()} F${riskData.fhour} · ${riskData.stations.length} stations`
                  : `Scanned ${riskData.stations.length} stations at ${riskData.date}`}
              </p>
            )}
            <StationPicker
              filteredStations={filteredStations}
              selectedStation={station}
              onStationSelect={(id) => onStationSelect(id, !!riskData)}
              stationSearch={stationSearch}
              onStationSearchChange={setStationSearch}
              sortMode={sortMode}
              onSortModeChange={setSortMode}
              favorites={favorites}
              onToggleFavorite={onToggleFavorite}
              riskData={riskData}
              listRef={listRef}
            />
          </div>
          {lat && lon && !pointMode && (
            <div className="cp-coords-preview">
              <MapPin size={12} />
              <span>{parseFloat(lat).toFixed(2)}°N, {Math.abs(parseFloat(lon)).toFixed(2)}°{parseFloat(lon) < 0 ? "W" : "E"}</span>
              <span className="cp-coords-hint">· station selected</span>
            </div>
          )}
          {pointMode && lat && lon && (
            <div className="cp-point-banner">
              <div className="cp-point-row">
                <Crosshair size={13} />
                <span className="cp-point-tag">POINT SOUNDING</span>
                <button
                  type="button"
                  className="cp-point-clear"
                  onClick={() => {
                    setPointMode(false);
                    const stn = stations.find((s) => s.id === station);
                    if (stn) {
                      setLat(String(stn.lat));
                      setLon(String(stn.lon));
                    }
                  }}
                  title="Switch back to station mode"
                >
                  <X size={12} />
                </button>
              </div>
              <span className="cp-point-coords">
                {parseFloat(lat).toFixed(2)}°N, {Math.abs(parseFloat(lon)).toFixed(2)}°{parseFloat(lon) < 0 ? "W" : "E"}
              </span>
            </div>
          )}
          <p className="cp-map-click-hint">
            <Crosshair size={11} />
            Click map for point sounding
          </p>
        </div>
      )}

      {/* BUFKIT Model */}
      {needsModel && (() => {
        const modelList = source === "psu" ? (psuModels || []) : models;
        const meta = MODEL_META[model] || { maxF: 384, step: 1 };
        const maxF = meta.maxF;
        const step = meta.step;
        const currentFhour = Math.min(parseInt(fhour) || 0, maxF);
        return (
          <div className="cp-section">
            <label className="cp-label">
              <Layers size={14} />
              Model
            </label>
            <div className="cp-model-grid">
              {modelList.map((m) => {
                const mm = MODEL_META[m.id];
                return (
                  <button
                    key={m.id}
                    type="button"
                    className={`cp-model-btn ${model === m.id ? "active" : ""}`}
                    onClick={() => {
                      setModel(m.id);
                      // Clamp forecast hour to new model's max
                      const newMeta = MODEL_META[m.id] || { maxF: 384, step: 1 };
                      const cur = parseInt(fhour) || 0;
                      if (cur > newMeta.maxF) setFhour(String(newMeta.maxF));
                    }}
                    title={m.name}
                  >
                    <span className="cp-model-id">{mm?.short || m.id.toUpperCase()}</span>
                    {mm && <span className="cp-model-res">{mm.res}</span>}
                  </button>
                );
              })}
            </div>

            {/* Forecast hour */}
            <div className="cp-fhour-section">
              <div className="cp-fhour-header">
                <Clock size={13} />
                <span>Forecast Hour</span>
                <span className="cp-fhour-value">F{String(currentFhour).padStart(2, "0")}</span>
              </div>
              <input
                type="range"
                min="0"
                max={maxF}
                step={step}
                className="cp-fhour-slider"
                value={currentFhour}
                onChange={(e) => setFhour(e.target.value)}
              />
              <div className="cp-fhour-range">
                <span>F00</span>
                <span>F{String(maxF).padStart(2, "0")}</span>
              </div>
              <div className="cp-fhour-presets">
                {FHOUR_PRESETS.filter((h) => h <= maxF).map((h) => (
                  <button
                    key={h}
                    type="button"
                    className={`cp-fhour-preset ${currentFhour === h ? "active" : ""}`}
                    onClick={() => setFhour(String(h))}
                  >
                    {h === 0 ? "Anl" : `F${String(h).padStart(2, "0")}`}
                  </button>
                ))}
              </div>
            </div>
          </div>
        );
      })()}

      {/* Date */}
      <div className="cp-section">
        <label className="cp-label">
          <Calendar size={14} />
          Date / Time (UTC)
        </label>
        {source === "obs" ? (
          <div className="cp-dt-picker">
            <div className="cp-dt-hours">
              {[
                { id: "latest", label: "Latest" },
                { id: "00", label: "00Z" },
                { id: "12", label: "12Z" },
              ].map((h) => (
                <button
                  key={h.id}
                  type="button"
                  className={`cp-dt-hour${soundingHour === h.id ? " active" : ""}`}
                  onClick={() => {
                    setSoundingHour(h.id);
                    if (h.id === "latest") {
                      setDate("");
                    } else if (date) {
                      setDate(`${date.slice(0, 10)}T${h.id}:00`);
                    }
                  }}
                >
                  {soundingHour === h.id && <Clock size={11} />}
                  {h.label}
                </button>
              ))}
            </div>
            <div className="cp-dt-date-wrap">
              <Calendar size={13} className="cp-dt-date-icon" />
              <input
                type="date"
                className="cp-dt-date-input"
                value={date ? date.slice(0, 10) : ""}
                onChange={(e) => {
                  const d = e.target.value;
                  if (d) {
                    const hour = soundingHour === "12" ? "12:00" : "00:00";
                    if (soundingHour === "latest") setSoundingHour("00");
                    setDate(`${d}T${hour}`);
                  } else {
                    setDate("");
                  }
                }}
              />
              {date && (
                <button type="button" className="cp-dt-clear" onClick={() => { setDate(""); setSoundingHour("latest"); }} title="Reset to latest">
                  <X size={12} />
                </button>
              )}
            </div>
            {soundingHour === "latest" && (
              <p className="cp-hint" style={{ margin: 0 }}>Auto-selects the most recent available sounding</p>
            )}
          </div>
        ) : (
          <div className="cp-dt-picker">
            <div className="cp-dt-date-wrap">
              <Clock size={13} className="cp-dt-date-icon" />
              <input
                type="datetime-local"
                className="cp-dt-date-input"
                value={date}
                onChange={(e) => setDate(e.target.value)}
              />
              {date && (
                <button type="button" className="cp-dt-clear" onClick={() => setDate("")} title="Clear date">
                  <X size={12} />
                </button>
              )}
            </div>
            <p className="cp-hint" style={{ margin: 0 }}>Leave blank for the most recent time</p>
          </div>
        )}
      </div>

      {/* Submit */}
      <button
        type="submit"
        className="cp-submit"
        disabled={loading}
      >
        {loading ? (
          <>
            <Loader2 size={16} className="spin" />
            Analyzing...
          </>
        ) : (
          <>
            <Search size={16} />
            Generate Sounding
          </>
        )}
      </button>
    </form>
  );
}
