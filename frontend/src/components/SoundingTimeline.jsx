import { useState, useEffect, useRef, useMemo } from "react";
import { Clock, ChevronLeft, ChevronRight } from "lucide-react";
import { getHistory } from "../history";
import "./SoundingTimeline.css";

/**
 * Build an array of standard sounding launch times:
 * 00Z and 12Z for the last `days` days, plus the upcoming time slot.
 */
function buildTimeline(days = 4) {
  const now = new Date();
  const slots = [];
  // Start from `days` ago at 00Z
  const start = new Date(Date.UTC(
    now.getUTCFullYear(), now.getUTCMonth(), now.getUTCDate() - days, 0, 0, 0
  ));
  const end = new Date(Date.UTC(
    now.getUTCFullYear(), now.getUTCMonth(), now.getUTCDate(), 
    now.getUTCHours() < 12 ? 12 : 24, 0, 0
  ));
  for (let t = new Date(start); t <= end; t.setUTCHours(t.getUTCHours() + 12)) {
    slots.push(new Date(t));
  }
  return slots;
}

/** Format date to YYYYMMDDHH for comparison with params. */
function toDateKey(d) {
  return d.getUTCFullYear().toString()
    + String(d.getUTCMonth() + 1).padStart(2, "0")
    + String(d.getUTCDate()).padStart(2, "0")
    + String(d.getUTCHours()).padStart(2, "0");
}

/** Short label like "Mar 7 12Z" */
function shortLabel(d) {
  const months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
  return `${months[d.getUTCMonth()]} ${d.getUTCDate()} ${String(d.getUTCHours()).padStart(2,"0")}Z`;
}

/** Day label like "Mar 7" */
function dayLabel(d) {
  const months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
  return `${months[d.getUTCMonth()]} ${d.getUTCDate()}`;
}

export default function SoundingTimeline({ station, currentDate, onSelectTime, source }) {
  const scrollRef = useRef(null);
  const activeRef = useRef(null);
  const [history, setHistory] = useState([]);

  // Load history for this station
  useEffect(() => {
    const all = getHistory();
    const filtered = all.filter(
      (h) => h.requestParams?.station === station || h.meta?.station === station
    );
    setHistory(filtered);
  }, [station]);

  const slots = useMemo(() => buildTimeline(4), []);

  // Set of YYYYMMDDHH that exist in history for this station
  const historyDateSet = useMemo(() => {
    const set = new Set();
    for (const h of history) {
      if (h.requestParams?.date) set.add(h.requestParams.date);
      // Also parse from meta.date (e.g., "2025-03-07 12:00Z")
      if (h.meta?.date) {
        const m = h.meta.date.match(/(\d{4})-(\d{2})-(\d{2})\s+(\d{2})/);
        if (m) set.add(m[1] + m[2] + m[3] + m[4]);
      }
    }
    return set;
  }, [history]);

  // Current active key from loaded sounding
  const activeKey = useMemo(() => {
    if (!currentDate) return null;
    // currentDate could be "2025-03-07 12:00Z" or "2025030712"
    if (/^\d{10}$/.test(currentDate)) return currentDate;
    const m = currentDate.match(/(\d{4})-(\d{2})-(\d{2})\s+(\d{2})/);
    if (m) return m[1] + m[2] + m[3] + m[4];
    return null;
  }, [currentDate]);

  // Scroll to active node on mount
  useEffect(() => {
    if (activeRef.current) {
      activeRef.current.scrollIntoView({ behavior: "smooth", inline: "center", block: "nearest" });
    }
  }, [activeKey]);

  const handleScroll = (dir) => {
    if (scrollRef.current) {
      scrollRef.current.scrollBy({ left: dir * 200, behavior: "smooth" });
    }
  };

  const handleClick = (dateKey) => {
    if (onSelectTime) onSelectTime(dateKey);
  };

  // Determine if we're past a slot time (in the past)
  const now = Date.now();

  return (
    <div className="stl-bar">
      <div className="stl-label">
        <Clock size={12} />
        <span>Timeline</span>
      </div>
      <button className="stl-arrow" onClick={() => handleScroll(-1)} title="Scroll left">
        <ChevronLeft size={14} />
      </button>
      <div className="stl-track" ref={scrollRef}>
        {slots.map((slot, i) => {
          const key = toDateKey(slot);
          const isActive = key === activeKey;
          const inHistory = historyDateSet.has(key);
          const isPast = slot.getTime() < now;
          const isFuture = !isPast;
          const showDay = i === 0 || slot.getUTCHours() === 0;

          return (
            <div key={key} className="stl-slot-group">
              {showDay && <div className="stl-day-label">{dayLabel(slot)}</div>}
              <button
                ref={isActive ? activeRef : null}
                className={`stl-node ${isActive ? "stl-active" : ""} ${inHistory ? "stl-has-data" : ""} ${isFuture ? "stl-future" : ""}`}
                onClick={() => !isFuture && handleClick(key)}
                disabled={isFuture}
                title={`${shortLabel(slot)}${inHistory ? " (in history)" : ""}${isActive ? " (current)" : ""}`}
              >
                <span className="stl-hour">{String(slot.getUTCHours()).padStart(2, "0")}Z</span>
                {inHistory && <span className="stl-dot" />}
              </button>
            </div>
          );
        })}
      </div>
      <button className="stl-arrow" onClick={() => handleScroll(1)} title="Scroll right">
        <ChevronRight size={14} />
      </button>
    </div>
  );
}
