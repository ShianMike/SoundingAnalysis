import { useState, useEffect, useRef, useCallback, lazy, Suspense } from "react";
import { Play, Pause, SkipForward, SkipBack, Loader2 } from "lucide-react";
import { fetchForecastProfiles } from "../api";
import "./SoundingAnimator.css";

const InteractiveSkewT = lazy(() => import("./InteractiveSkewT"));
const InteractiveHodograph = lazy(() => import("./InteractiveHodograph"));

const DEFAULT_HOURS = [0, 1, 2, 3, 6, 9, 12, 15, 18, 21, 24];

export default function SoundingAnimator({ station, model, source, date, theme, onClose }) {
  const [frames, setFrames] = useState([]);
  const [idx, setIdx] = useState(0);
  const [playing, setPlaying] = useState(false);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [speed, setSpeed] = useState(1000); // ms per frame
  const intervalRef = useRef(null);

  /* Fetch all forecast profiles on mount / params change */
  useEffect(() => {
    let cancelled = false;
    (async () => {
      setLoading(true);
      setError(null);
      setFrames([]);
      setIdx(0);
      setPlaying(false);
      try {
        const res = await fetchForecastProfiles({
          station,
          model: model || "rap",
          source: source || "psu",
          date,
          hours: DEFAULT_HOURS,
        });
        if (!cancelled) setFrames(res.frames || []);
      } catch (e) {
        if (!cancelled) setError(e.message);
      } finally {
        if (!cancelled) setLoading(false);
      }
    })();
    return () => { cancelled = true; };
  }, [station, model, source, date]);

  /* Play / pause timer */
  useEffect(() => {
    if (playing && frames.length > 1) {
      intervalRef.current = setInterval(() => {
        setIdx((prev) => (prev + 1) % frames.length);
      }, speed);
    }
    return () => clearInterval(intervalRef.current);
  }, [playing, frames.length, speed]);

  const step = useCallback(
    (dir) => {
      setPlaying(false);
      setIdx((prev) => {
        const next = prev + dir;
        if (next < 0) return frames.length - 1;
        if (next >= frames.length) return 0;
        return next;
      });
    },
    [frames.length]
  );

  const frame = frames[idx];

  return (
    <div className="animator-panel">
      <div className="animator-header">
        <span className="animator-title">Sounding Animation</span>
        <button className="animator-close" onClick={onClose} title="Close animator">✕</button>
      </div>

      {loading && (
        <div className="animator-loading">
          <Loader2 size={20} className="spin" />
          <span>Loading forecast profiles…</span>
        </div>
      )}
      {error && <div className="animator-error">{error}</div>}

      {!loading && frames.length > 0 && frame && (
        <>
          {/* Controls */}
          <div className="animator-controls">
            <button onClick={() => step(-1)} title="Previous frame"><SkipBack size={16} /></button>
            <button onClick={() => setPlaying((v) => !v)} title={playing ? "Pause" : "Play"}>
              {playing ? <Pause size={16} /> : <Play size={16} />}
            </button>
            <button onClick={() => step(1)} title="Next frame"><SkipForward size={16} /></button>
            <input
              type="range"
              min={0}
              max={frames.length - 1}
              value={idx}
              onChange={(e) => { setPlaying(false); setIdx(Number(e.target.value)); }}
              className="animator-slider"
            />
            <span className="animator-label">f{String(frame.fhour).padStart(3, "0")}</span>
            <select
              className="animator-speed"
              value={speed}
              onChange={(e) => setSpeed(Number(e.target.value))}
              title="Animation speed"
            >
              <option value={2000}>Slow</option>
              <option value={1000}>Normal</option>
              <option value={500}>Fast</option>
              <option value={250}>Very Fast</option>
            </select>
          </div>

          {/* Skew-T + Hodograph for current frame */}
          <Suspense fallback={null}>
            <InteractiveSkewT
              profile={frame.profile}
              sbParcel={frame.sbParcel}
              mlParcel={frame.mlParcel}
              params={frame.params}
              theme={theme || "dark"}
            />
            <InteractiveHodograph
              profile={frame.profile}
              params={frame.params}
              theme={theme || "dark"}
            />
          </Suspense>
        </>
      )}

      {!loading && frames.length === 0 && !error && (
        <div className="animator-empty">No forecast frames available.</div>
      )}
    </div>
  );
}
