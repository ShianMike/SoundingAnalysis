import { useState, useEffect, useMemo, useCallback } from "react";
import { X, Crosshair, Eye, Radio, Wifi, MapPin, Users } from "lucide-react";

const LSC_API = "https://api.livestormchasing.com/api/v1/chasers";
const LSC_THUMB = "https://edge.livestormchasing.com/thumbnails/";
const LSC_STREAM = "https://livestormchasing.com/chasers/";
const LSC_HLS = "https://edge.livestormchasing.com/hls/";

function timeAgo(ts) {
  if (!ts) return "";
  const d = new Date(ts + " UTC");
  const s = Math.floor((Date.now() - d) / 1000);
  if (s < 60) return "just now";
  if (s < 3600) return `${Math.floor(s / 60)}m ago`;
  if (s < 86400) return `${Math.floor(s / 3600)}h ago`;
  return `${Math.floor(s / 86400)}d ago`;
}

export default function ChaserPanel({ onFlyTo, onClose }) {
  const [chasers, setChasers] = useState([]);
  const [loading, setLoading] = useState(true);
  const [search, setSearch] = useState("");
  const [activeStream, setActiveStream] = useState(null); // chaser id being watched

  const fetchChasers = useCallback(async () => {
    try {
      const res = await fetch(LSC_API);
      if (!res.ok) return;
      const data = await res.json();
      setChasers(data);
    } catch (e) {
      console.error("Live chasers fetch error:", e);
    } finally {
      setLoading(false);
    }
  }, []);

  useEffect(() => {
    fetchChasers();
    const id = setInterval(fetchChasers, 30_000);
    return () => clearInterval(id);
  }, [fetchChasers]);

  const { live, position } = useMemo(() => {
    const q = search.toLowerCase().trim();
    const all = chasers.map((c) => {
      const p = c.properties;
      const [lon, lat] = c.geometry?.coordinates || [0, 0];
      return {
        id: p.id,
        name: p.name,
        location: p.location,
        isLive: p.stream_status,
        viewers: p.viewers || 0,
        streamId: p.stream_id,
        username: p.username,
        gpsUpdate: p.gps_update,
        liveTimestamp: p.live_timestamp,
        lat,
        lon,
      };
    });
    const filtered = q
      ? all.filter((c) => c.name.toLowerCase().includes(q) || (c.location || "").toLowerCase().includes(q))
      : all;
    return {
      live: filtered.filter((c) => c.isLive).sort((a, b) => b.viewers - a.viewers),
      position: filtered.filter((c) => !c.isLive),
    };
  }, [chasers, search]);

  const handleWatch = (chaser) => {
    if (activeStream === chaser.id) {
      setActiveStream(null);
    } else {
      setActiveStream(chaser.id);
      onFlyTo(chaser.lat, chaser.lon);
    }
  };

  const openExternal = (chaser) => {
    const slug = chaser.username || chaser.id;
    window.open(LSC_STREAM + slug, "_blank", "noopener");
  };

  return (
    <div className="chaser-panel">
      <div className="chaser-header">
        <div className="chaser-header-left">
          <Radio size={14} className="chaser-header-icon" />
          <span className="chaser-header-title">Live Storm Chasers</span>
          {live.length > 0 && <span className="chaser-live-badge">{live.length} LIVE</span>}
        </div>
        <button className="chaser-close" onClick={onClose}><X size={14} /></button>
      </div>

      <input
        className="chaser-search"
        type="text"
        placeholder="Search chasers..."
        value={search}
        onChange={(e) => setSearch(e.target.value)}
      />

      <div className="chaser-list">
        {loading && <div className="chaser-empty"><span className="chaser-spinner" /> Loading chasers...</div>}
        {!loading && live.length === 0 && position.length === 0 && (
          <div className="chaser-empty">No chasers found</div>
        )}

        {live.length > 0 && (
          <div className="chaser-section">
            <div className="chaser-section-label"><Wifi size={10} /> Live Streams</div>
            {live.map((c) => (
              <div key={c.id} className={`chaser-card chaser-card-live ${activeStream === c.id ? "chaser-card-watching" : ""}`}>
                <div className="chaser-card-row" onClick={() => handleWatch(c)}>
                  <div className="chaser-avatar">
                    <div className="chaser-avatar-pulse" />
                    <Wifi size={12} />
                  </div>
                  <div className="chaser-card-info">
                    <div className="chaser-card-name">
                      {c.name}
                      <span className="chaser-card-live-tag">
                        <span className="chaser-live-dot" /> LIVE
                      </span>
                    </div>
                    <div className="chaser-card-meta">
                      <MapPin size={9} /> {c.location || "Unknown"}
                      <span className="chaser-card-sep">·</span>
                      <Users size={9} /> {c.viewers}
                    </div>
                    <div className="chaser-card-time">{timeAgo(c.liveTimestamp)}</div>
                  </div>
                  <div className="chaser-card-actions">
                    <button className="chaser-action-btn chaser-fly" onClick={(e) => { e.stopPropagation(); onFlyTo(c.lat, c.lon); }} title="Fly to"><Crosshair size={12} /></button>
                    <button className="chaser-action-btn chaser-watch" onClick={(e) => { e.stopPropagation(); openExternal(c); }} title="Watch on livestormchasing.com"><Eye size={12} /></button>
                  </div>
                </div>
                {activeStream === c.id && (
                  <div className="chaser-stream-embed">
                    <iframe
                      src={`${LSC_STREAM}${c.username || c.id}`}
                      title={`${c.name} live stream`}
                      allow="autoplay; fullscreen"
                      allowFullScreen
                      className="chaser-iframe"
                    />
                  </div>
                )}
              </div>
            ))}
          </div>
        )}

        {position.length > 0 && (
          <div className="chaser-section">
            <div className="chaser-section-label"><MapPin size={10} /> Position ({position.length})</div>
            {position.map((c) => (
              <div key={c.id} className="chaser-card chaser-card-pos">
                <div className="chaser-card-row">
                  <div className="chaser-avatar chaser-avatar-pos">
                    <MapPin size={12} />
                  </div>
                  <div className="chaser-card-info">
                    <div className="chaser-card-name">{c.name}</div>
                    <div className="chaser-card-meta">
                      <MapPin size={9} /> {c.location || "Unknown"}
                    </div>
                    <div className="chaser-card-time">{timeAgo(c.gpsUpdate)}</div>
                  </div>
                  <div className="chaser-card-actions">
                    <button className="chaser-action-btn chaser-fly" onClick={() => onFlyTo(c.lat, c.lon)} title="Fly to"><Crosshair size={12} /></button>
                  </div>
                </div>
              </div>
            ))}
          </div>
        )}
      </div>

      <div className="chaser-footer">
        <span>Data: livestormchasing.com</span>
        <span>{chasers.length} chaser{chasers.length !== 1 ? "s" : ""}</span>
      </div>
    </div>
  );
}
