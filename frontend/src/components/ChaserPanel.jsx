import { useState, useMemo } from "react";
import { X, Crosshair, AlertTriangle, Eye } from "lucide-react";

const SN_REPORT_COLORS = {
  Tornado: "#ff0000",
  "Funnel Cloud": "#ff6600",
  "Rotating Wall Cloud": "#ff9900",
  Hail: "#00ccff",
  "Wind Damage": "#ffaa00",
  Flooding: "#00ff66",
  Other: "#cccccc",
};

export default function ChaserPanel({ spotterData, onFlyTo, onClose }) {
  const [tab, setTab] = useState("all"); // "all" | "reports"
  const [search, setSearch] = useState("");

  const positions = spotterData?.positions || [];
  const reports = spotterData?.reports || [];

  const filtered = useMemo(() => {
    const q = search.toLowerCase().trim();

    const reportItems = reports.map((r) => ({
      kind: "report",
      name: r.reporter,
      detail: r.type,
      time: r.time,
      notes: r.notes,
      lat: r.lat,
      lon: r.lon,
      color: SN_REPORT_COLORS[r.type] || SN_REPORT_COLORS.Other,
      isTornado: r.type === "Tornado",
      age: r.sheet === 3 ? "Recent" : r.sheet === 4 ? "Older" : "Old",
    }));

    const posItems = positions.map((s) => ({
      kind: "position",
      name: s.name,
      detail: s.heading || "Active",
      time: s.time,
      lat: s.lat,
      lon: s.lon,
      color: "#00cc44",
    }));

    let items = tab === "reports" ? reportItems : [...reportItems, ...posItems];

    if (q) {
      items = items.filter(
        (it) =>
          it.name.toLowerCase().includes(q) ||
          it.detail.toLowerCase().includes(q) ||
          (it.notes || "").toLowerCase().includes(q)
      );
    }

    return items;
  }, [positions, reports, tab, search]);

  return (
    <div className="chaser-panel">
      <div className="chaser-panel-header">
        <span className="chaser-panel-title">
          <Eye size={14} />
          Live Chasers
        </span>
        <button className="chaser-panel-close" onClick={onClose}>
          <X size={14} />
        </button>
      </div>

      <div className="chaser-panel-tabs">
        <button
          className={`chaser-tab ${tab === "all" ? "active" : ""}`}
          onClick={() => setTab("all")}
        >
          All ({positions.length + reports.length})
        </button>
        <button
          className={`chaser-tab ${tab === "reports" ? "active" : ""}`}
          onClick={() => setTab("reports")}
        >
          Reports ({reports.length})
        </button>
      </div>

      <input
        className="chaser-search"
        type="text"
        placeholder="Search chasers..."
        value={search}
        onChange={(e) => setSearch(e.target.value)}
      />

      <div className="chaser-list">
        {filtered.length === 0 && (
          <div className="chaser-empty">No chasers found</div>
        )}
        {filtered.map((it, i) => (
          <div
            key={`${it.kind}-${it.lat}-${it.lon}-${i}`}
            className={`chaser-item ${it.kind === "report" ? "chaser-item-report" : ""} ${it.isTornado ? "chaser-item-tornado" : ""}`}
          >
            <div className="chaser-item-left" style={{ borderColor: it.color }}>
              {it.kind === "report" ? (
                <span className="chaser-badge chaser-badge-report" style={{ background: it.color }}>
                  {it.isTornado ? "TOR" : it.detail.slice(0, 4).toUpperCase()}
                </span>
              ) : (
                <span className="chaser-badge chaser-badge-active">ACTIVE</span>
              )}
            </div>
            <div className="chaser-item-info">
              <span className="chaser-name">{it.name || "Unknown"}</span>
              <span className="chaser-detail">
                {it.kind === "report" ? (
                  <>
                    <span style={{ color: it.color }}>{it.detail}</span>
                    {it.age && <span className="chaser-age"> · {it.age}</span>}
                  </>
                ) : (
                  <span>{it.detail}</span>
                )}
              </span>
              {it.time && <span className="chaser-time">{it.time}</span>}
              {it.notes && <span className="chaser-notes">{it.notes}</span>}
            </div>
            <div className="chaser-item-actions">
              <button
                className="chaser-fly-btn"
                onClick={() => onFlyTo(it.lat, it.lon)}
                title="Fly to location"
              >
                <Crosshair size={13} />
              </button>
            </div>
          </div>
        ))}
      </div>

      <div className="chaser-panel-footer">
        <span>Data: Spotter Network</span>
        <span>{filtered.length} chaser{filtered.length !== 1 ? "s" : ""}</span>
      </div>
    </div>
  );
}
