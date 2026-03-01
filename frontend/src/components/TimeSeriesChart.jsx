import { useState, useEffect, useRef, useCallback } from "react";
import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer,
  ReferenceLine,
} from "recharts";
import { TrendingUp, Loader2, AlertTriangle, X, ChevronDown, Info, Calendar, Play, Pause, SkipBack, SkipForward, RotateCcw } from "lucide-react";
import { fetchTimeSeries } from "../api";
import "./TimeSeriesChart.css";

/** Available parameter groups for charting */
const PARAM_GROUPS = {
  "Severe Weather Indices": {
    stp: { label: "STP", color: "#60a5fa", unit: "" },
    scp: { label: "SCP", color: "#f59e0b", unit: "" },
    ship: { label: "SHIP", color: "#10b981", unit: "" },
    dcp: { label: "DCP", color: "#a78bfa", unit: "" },
  },
  "CAPE / CIN": {
    sbCape: { label: "SB CAPE", color: "#f97316", unit: "J/kg" },
    muCape: { label: "MU CAPE", color: "#ef4444", unit: "J/kg" },
    mlCape: { label: "ML CAPE", color: "#d946ef", unit: "J/kg" },
    dcape: { label: "DCAPE", color: "#facc15", unit: "J/kg" },
  },
  "Shear / Helicity": {
    bwd1km: { label: "BWD 0-1km", color: "#ef4444", unit: "kt" },
    bwd3km: { label: "BWD 0-3km", color: "#f97316", unit: "kt" },
    bwd6km: { label: "BWD 0-6km", color: "#eab308", unit: "kt" },
    srh1km: { label: "SRH 0-1km", color: "#22d3ee", unit: "m²/s²" },
    srh3km: { label: "SRH 0-3km", color: "#818cf8", unit: "m²/s²" },
  },
  "Lapse Rates & Moisture": {
    lr03: { label: "LR 0-3km", color: "#fb923c", unit: "°C/km" },
    lr36: { label: "LR 3-6km", color: "#a3e635", unit: "°C/km" },
    pwat: { label: "PWAT", color: "#38bdf8", unit: "mm" },
  },
};

// Flatten for lookup
const ALL_PARAMS = {};
Object.values(PARAM_GROUPS).forEach((group) => {
  Object.entries(group).forEach(([key, cfg]) => {
    ALL_PARAMS[key] = cfg;
  });
});

const DEFAULT_SELECTED = ["stp", "scp", "ship"];

/** Custom tooltip */
function CustomTooltip({ active, payload, label }) {
  if (!active || !payload?.length) return null;
  return (
    <div className="ts-tooltip">
      <div className="ts-tooltip-date">{label}</div>
      {payload.map((p) => (
        <div key={p.dataKey} className="ts-tooltip-row">
          <span className="ts-tooltip-dot" style={{ background: p.color }} />
          <span className="ts-tooltip-label">{ALL_PARAMS[p.dataKey]?.label || p.dataKey}</span>
          <span className="ts-tooltip-value">
            {p.value != null ? p.value : "---"}
            {ALL_PARAMS[p.dataKey]?.unit ? ` ${ALL_PARAMS[p.dataKey].unit}` : ""}
          </span>
        </div>
      ))}
    </div>
  );
}

export default function TimeSeriesChart({ station, source, onClose }) {
  const [data, setData] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [selected, setSelected] = useState(DEFAULT_SELECTED);
  const [expandedGroup, setExpandedGroup] = useState("Severe Weather Indices");

  // Playback state
  const [playIdx, setPlayIdx] = useState(-1); // -1 = off
  const [playing, setPlaying] = useState(false);
  const [playSpeed, setPlaySpeed] = useState(1000); // ms per frame
  const playTimer = useRef(null);

  // Date range: default to last 7 days
  const todayStr = new Date().toISOString().slice(0, 10);
  const weekAgoStr = new Date(Date.now() - 7 * 86400000).toISOString().slice(0, 10);
  const [startDate, setStartDate] = useState(weekAgoStr);
  const [endDate, setEndDate] = useState(todayStr);

  // Enforce max 14-day span
  const maxStartDate = (() => {
    const d = new Date(endDate);
    d.setDate(d.getDate() - 14);
    return d.toISOString().slice(0, 10);
  })();

  const spanDays = Math.round(
    (new Date(endDate) - new Date(startDate)) / 86400000
  );

  const canFetch = !!station;

  const handleFetch = async () => {
    if (!canFetch) return;
    setLoading(true);
    setError(null);
    try {
      const fmt = (d) => d.replace(/-/g, "") + "00"; // YYYYMMDD00
      const result = await fetchTimeSeries({
        station,
        source: source || "obs",
        startDate: fmt(startDate),
        endDate: fmt(endDate),
      });
      setData(result);
      if (result.points.length === 0) {
        setError("No sounding data found for this date range.");
      }
    } catch (e) {
      if (e.message === "Request timed out" || e.message.includes("ERR_FAILED")) {
        setError("Request timed out. Try a shorter date range.");
      } else {
        setError(e.message);
      }
    } finally {
      setLoading(false);
    }
  };

  const toggleParam = (key) => {
    setSelected((prev) =>
      prev.includes(key) ? prev.filter((k) => k !== key) : [...prev, key]
    );
  };

  /** Clicking a group header selects that group exclusively */
  const selectGroup = (groupName) => {
    if (expandedGroup === groupName) {
      setExpandedGroup(null);
      return;
    }
    setExpandedGroup(groupName);
    const groupKeys = Object.keys(PARAM_GROUPS[groupName]);
    setSelected(groupKeys);
  };

  // Playback controls
  const totalFrames = data?.points?.length || 0;

  const stopPlayback = useCallback(() => {
    setPlaying(false);
    if (playTimer.current) { clearInterval(playTimer.current); playTimer.current = null; }
  }, []);

  const startPlayback = useCallback(() => {
    if (totalFrames < 2) return;
    setPlaying(true);
    setPlayIdx((prev) => (prev < 0 || prev >= totalFrames - 1) ? 0 : prev);
  }, [totalFrames]);

  const togglePlayback = useCallback(() => {
    if (playing) stopPlayback(); else startPlayback();
  }, [playing, stopPlayback, startPlayback]);

  const stepForward = useCallback(() => {
    setPlayIdx((prev) => Math.min((prev < 0 ? 0 : prev) + 1, totalFrames - 1));
  }, [totalFrames]);

  const stepBack = useCallback(() => {
    setPlayIdx((prev) => Math.max((prev < 0 ? 0 : prev) - 1, 0));
  }, []);

  const resetPlayback = useCallback(() => {
    stopPlayback();
    setPlayIdx(-1);
  }, [stopPlayback]);

  // Auto-advance when playing
  useEffect(() => {
    if (!playing) return;
    playTimer.current = setInterval(() => {
      setPlayIdx((prev) => {
        const next = prev + 1;
        if (next >= totalFrames) {
          setPlaying(false);
          clearInterval(playTimer.current);
          playTimer.current = null;
          return totalFrames - 1;
        }
        return next;
      });
    }, playSpeed);
    return () => { if (playTimer.current) clearInterval(playTimer.current); };
  }, [playing, playSpeed, totalFrames]);

  // Transform data for Recharts
  const chartData =
    data?.points?.map((p) => ({
      date: p.date,
      ...p.params,
    })) || [];

  return (
    <div className="ts-panel">
      {/* Header */}
      <div className="ts-header">
        <div className="ts-header-left">
          <TrendingUp size={16} />
          <h3>Parameter Time-Series</h3>
          {data && (
            <span className="ts-header-meta">
              {data.stationName} ({data.station}) · {data.points.length} soundings
            </span>
          )}
        </div>
        <button className="ts-close" onClick={onClose} title="Close">
          <X size={16} />
        </button>
      </div>

      {/* Controls */}
      <div className="ts-controls">
        <div className="ts-date-range">
          <div className="ts-date-field">
            <label className="ts-label"><Calendar size={11} /> From</label>
            <input
              type="date"
              className="ts-date-input"
              value={startDate}
              min={maxStartDate}
              max={endDate}
              onChange={(e) => setStartDate(e.target.value)}
            />
          </div>
          <span className="ts-date-arrow">→</span>
          <div className="ts-date-field">
            <label className="ts-label"><Calendar size={11} /> To</label>
            <input
              type="date"
              className="ts-date-input"
              value={endDate}
              max={todayStr}
              min={startDate}
              onChange={(e) => setEndDate(e.target.value)}
            />
          </div>
        </div>

        <button
          className="ts-fetch-btn"
          onClick={handleFetch}
          disabled={loading || !canFetch}
        >
          {loading ? (
            <>
              <Loader2 size={14} className="spin" />
              Fetching...
            </>
          ) : (
            <>
              <TrendingUp size={14} />
              {data ? "Refresh" : "Load Time-Series"}
            </>
          )}
        </button>
      </div>

      {/* Info banner */}
      <div className="ts-info">
        <Info size={13} />
        <span>
          {spanDays <= 7
            ? `${spanDays}-day range — sampling at 00Z & 12Z (${spanDays * 2} soundings).`
            : `${spanDays}-day range — sampling at 12Z only (${spanDays} soundings). Ranges ≤ 7 days include both 00Z & 12Z.`}
          {" "}Max range: 14 days.
        </span>
      </div>

      {/* Parameter selector */}
      <div className="ts-param-selector">
        {Object.entries(PARAM_GROUPS).map(([groupName, params]) => {
          const groupKeys = Object.keys(params);
          const activeCount = groupKeys.filter((k) => selected.includes(k)).length;
          return (
            <div key={groupName} className="ts-param-group">
              <button
                className={`ts-param-group-header ${expandedGroup === groupName ? "expanded" : ""} ${activeCount > 0 ? "has-active" : ""}`}
                onClick={() => selectGroup(groupName)}
              >
                <ChevronDown size={12} />
                <span>{groupName}</span>
                {activeCount > 0 && (
                  <span className="ts-group-badge">{activeCount}</span>
                )}
              </button>
              {expandedGroup === groupName && (
                <div className="ts-param-chips">
                  {Object.entries(params).map(([key, cfg]) => (
                    <button
                      key={key}
                      className={`ts-chip ${selected.includes(key) ? "active" : ""}`}
                      onClick={() => toggleParam(key)}
                      style={
                        selected.includes(key)
                          ? { borderColor: cfg.color, color: cfg.color }
                          : {}
                      }
                    >
                      <span
                        className="ts-chip-dot"
                        style={{ background: selected.includes(key) ? cfg.color : "#666" }}
                      />
                      {cfg.label}
                    </button>
                  ))}
                </div>
              )}
            </div>
          );
        })}
      </div>

      {/* Error */}
      {error && (
        <div className="ts-error">
          <AlertTriangle size={14} />
          <span>{error}</span>
        </div>
      )}

      {/* Playback controls */}
      {data && chartData.length > 1 && (
        <div className="ts-playback">
          <button className="ts-play-btn" onClick={resetPlayback} title="Reset">
            <RotateCcw size={13} />
          </button>
          <button className="ts-play-btn" onClick={stepBack} title="Step back" disabled={playIdx <= 0}>
            <SkipBack size={13} />
          </button>
          <button className={`ts-play-btn ts-play-main ${playing ? "active" : ""}`} onClick={togglePlayback} title={playing ? "Pause" : "Play"}>
            {playing ? <Pause size={14} /> : <Play size={14} />}
          </button>
          <button className="ts-play-btn" onClick={stepForward} title="Step forward" disabled={playIdx >= totalFrames - 1}>
            <SkipForward size={13} />
          </button>
          <select className="ts-speed-select" value={playSpeed} onChange={(e) => setPlaySpeed(Number(e.target.value))}>
            <option value={2000}>0.5x</option>
            <option value={1000}>1x</option>
            <option value={500}>2x</option>
            <option value={250}>4x</option>
          </select>
          {playIdx >= 0 && (
            <span className="ts-play-label">
              {chartData[playIdx]?.date} ({playIdx + 1}/{totalFrames})
            </span>
          )}
        </div>
      )}

      {/* Chart */}
      {data && chartData.length > 0 && (
        <div className="ts-chart-wrap">
          <ResponsiveContainer width="100%" height={360}>
            <LineChart data={chartData} margin={{ top: 8, right: 24, left: 4, bottom: 8 }}>
              <CartesianGrid strokeDasharray="3 3" stroke="#333" />
              <XAxis
                dataKey="date"
                tick={{ fill: "#999", fontSize: 11 }}
                tickLine={{ stroke: "#555" }}
                axisLine={{ stroke: "#555" }}
                angle={-30}
                textAnchor="end"
                height={60}
              />
              <YAxis
                tick={{ fill: "#999", fontSize: 11 }}
                tickLine={{ stroke: "#555" }}
                axisLine={{ stroke: "#555" }}
              />
              <Tooltip content={<CustomTooltip />} />
              <Legend
                wrapperStyle={{ paddingTop: 8, fontSize: 11, color: "#ccc" }}
                iconSize={10}
                formatter={(value) => <span style={{ color: "#aaa", fontSize: 11 }}>{value}</span>}
              />
              <ReferenceLine y={0} stroke="#555" strokeDasharray="2 2" />
              {playIdx >= 0 && chartData[playIdx] && (
                <ReferenceLine x={chartData[playIdx].date} stroke="#60a5fa" strokeWidth={2} strokeDasharray="4 2" />
              )}
              {selected.map(
                (key) =>
                  ALL_PARAMS[key] && (
                    <Line
                      key={key}
                      type="monotone"
                      dataKey={key}
                      name={ALL_PARAMS[key].label}
                      stroke={ALL_PARAMS[key].color}
                      strokeWidth={2}
                      dot={{ r: 3, fill: ALL_PARAMS[key].color }}
                      activeDot={{ r: 5 }}
                      connectNulls
                    />
                  )
              )}
            </LineChart>
          </ResponsiveContainer>
        </div>
      )}

      {/* Empty state */}
      {!data && !loading && !error && (
        <div className="ts-empty">
          <TrendingUp size={24} />
          <p>
            Choose a date range (up to 14 days) and click <strong>Load Time-Series</strong> to
            see how sounding parameters evolve over time.
          </p>
        </div>
      )}

      {/* No data */}
      {data && chartData.length === 0 && (
        <div className="ts-empty">
          <AlertTriangle size={24} />
          <p>No sounding data could be retrieved for this station and time range.</p>
        </div>
      )}
    </div>
  );
}
