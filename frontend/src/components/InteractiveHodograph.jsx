import { useRef, useEffect, useState, useCallback, useMemo } from "react";
import { Download, Crosshair } from "lucide-react";
import "./InteractiveHodograph.css";

/* ── Constants ────────────────────────────────────────── */
const MARGIN = { top: 20, right: 20, bottom: 30, left: 30 };
const SPEED_RINGS = [10, 20, 30, 40, 50, 60, 70, 80]; // kt
const HEIGHT_COLORS = [
  { maxKm: 1, color: "#ef4444", label: "0–1 km" },
  { maxKm: 3, color: "#f97316", label: "1–3 km" },
  { maxKm: 6, color: "#22c55e", label: "3–6 km" },
  { maxKm: 9, color: "#3b82f6", label: "6–9 km" },
  { maxKm: 12, color: "#a855f7", label: "9–12 km" },
  { maxKm: Infinity, color: "#888", label: ">12 km" },
];

function windComponents(wd, ws) {
  if (wd == null || ws == null) return null;
  const rad = ((270 - wd) * Math.PI) / 180;
  return { u: ws * Math.cos(rad), v: ws * Math.sin(rad) };
}

function heightColor(hAgl) {
  const km = hAgl / 1000;
  for (const hc of HEIGHT_COLORS) {
    if (km <= hc.maxKm) return hc.color;
  }
  return "#888";
}

export default function InteractiveHodograph({ profile, params, theme = "dark" }) {
  const canvasRef = useRef(null);
  const overlayRef = useRef(null);
  const containerRef = useRef(null);
  const [dims, setDims] = useState({ w: 500, h: 500 });
  const [cursor, setCursor] = useState(null);
  const [showCrosshair, setShowCrosshair] = useState(true);

  const colors = useMemo(
    () =>
      theme === "light"
        ? { bg: "#ffffff", grid: "#e0e0e0", ring: "#ccc", text: "#333", crosshair: "rgba(100,100,100,0.5)", crossText: "#333", zero: "#aaa" }
        : { bg: "#0d0d0d", grid: "#222", ring: "#333", text: "#ccc", crosshair: "rgba(200,200,200,0.4)", crossText: "#eee", zero: "#555" },
    [theme]
  );

  /* Responsive */
  useEffect(() => {
    const el = containerRef.current;
    if (!el) return;
    const ro = new ResizeObserver((entries) => {
      const { width } = entries[0].contentRect;
      const s = Math.max(300, Math.floor(width));
      setDims({ w: s, h: s });
    });
    ro.observe(el);
    return () => ro.disconnect();
  }, []);

  /* Build wind vectors with height info */
  const windData = useMemo(() => {
    if (!profile || profile.length === 0) return [];
    const sfcH = profile[0]?.h ?? 0;
    return profile
      .filter((lev) => lev.p != null && lev.wd != null && lev.ws != null && lev.h != null)
      .map((lev) => {
        const comp = windComponents(lev.wd, lev.ws);
        return { ...comp, hAgl: lev.h - sfcH, p: lev.p, wd: lev.wd, ws: lev.ws };
      });
  }, [profile]);

  /* Determine max speed for scale */
  const maxSpd = useMemo(() => {
    let m = 40; // minimum ring
    for (const w of windData) {
      const spd = Math.sqrt(w.u * w.u + w.v * w.v);
      if (spd > m) m = spd;
    }
    // Round up to next ring
    for (const r of SPEED_RINGS) {
      if (r >= m) return r;
    }
    return Math.ceil(m / 10) * 10;
  }, [windData]);

  const plotSize = dims.w - MARGIN.left - MARGIN.right;
  const cx = MARGIN.left + plotSize / 2;
  const cy = MARGIN.top + plotSize / 2;
  const scale = useCallback((spd) => (spd / maxSpd) * (plotSize / 2), [maxSpd, plotSize]);

  /* ── Main draw ── */
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas || windData.length === 0) return;
    const ctx = canvas.getContext("2d");
    const dpr = window.devicePixelRatio || 1;
    canvas.width = dims.w * dpr;
    canvas.height = dims.h * dpr;
    ctx.scale(dpr, dpr);

    ctx.fillStyle = colors.bg;
    ctx.fillRect(0, 0, dims.w, dims.h);

    // Speed rings
    ctx.strokeStyle = colors.ring;
    ctx.lineWidth = 0.5;
    ctx.setLineDash([3, 3]);
    ctx.font = "10px Inter, system-ui, sans-serif";
    ctx.fillStyle = colors.text;
    ctx.textAlign = "center";
    ctx.textBaseline = "bottom";
    for (const r of SPEED_RINGS) {
      if (r > maxSpd) break;
      const rad = scale(r);
      ctx.beginPath();
      ctx.arc(cx, cy, rad, 0, Math.PI * 2);
      ctx.stroke();
      ctx.fillText(`${r}`, cx, cy - rad - 2);
    }
    ctx.setLineDash([]);

    // Cross axes
    ctx.strokeStyle = colors.zero;
    ctx.lineWidth = 0.8;
    ctx.beginPath();
    ctx.moveTo(cx - plotSize / 2, cy);
    ctx.lineTo(cx + plotSize / 2, cy);
    ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(cx, cy - plotSize / 2);
    ctx.lineTo(cx, cy + plotSize / 2);
    ctx.stroke();

    // Axis labels
    ctx.fillStyle = colors.text;
    ctx.font = "11px Inter, system-ui, sans-serif";
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    ctx.fillText("u (kt)", cx, dims.h - 10);
    ctx.save();
    ctx.translate(10, cy);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "center";
    ctx.textBaseline = "bottom";
    ctx.fillText("v (kt)", 0, 0);
    ctx.restore();

    // Hodograph trace — color-coded by height
    ctx.lineWidth = 2.2;
    for (let i = 0; i < windData.length - 1; i++) {
      const w0 = windData[i], w1 = windData[i + 1];
      ctx.strokeStyle = heightColor(w0.hAgl);
      ctx.beginPath();
      ctx.moveTo(cx + scale(w0.u), cy - scale(w0.v));
      ctx.lineTo(cx + scale(w1.u), cy - scale(w1.v));
      ctx.stroke();
    }

    // Dots at each data point
    for (const w of windData) {
      ctx.fillStyle = heightColor(w.hAgl);
      ctx.beginPath();
      ctx.arc(cx + scale(w.u), cy - scale(w.v), 2, 0, Math.PI * 2);
      ctx.fill();
    }

    // Bunkers vectors
    const drawVector = (u, v, color, label) => {
      if (u == null || v == null) return;
      // Convert m/s to kt
      const uKt = u * 1.94384;
      const vKt = v * 1.94384;
      const px = cx + scale(uKt);
      const py = cy - scale(vKt);
      // Marker
      ctx.fillStyle = color;
      ctx.strokeStyle = color;
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.arc(px, py, 5, 0, Math.PI * 2);
      ctx.fill();
      // Label
      ctx.font = "bold 10px Inter, system-ui, sans-serif";
      ctx.textAlign = "left";
      ctx.textBaseline = "middle";
      ctx.fillText(label, px + 8, py);
    };

    if (params) {
      drawVector(params.rmU, params.rmV, "#ef4444", "RM");
      drawVector(params.lmU, params.lmV, "#3b82f6", "LM");
      drawVector(params.mwU, params.mwV, "#a3a3a3", "MW");
    }

    // Legend — height bands
    const legX = MARGIN.left + 6;
    let legY = MARGIN.top + 6;
    ctx.font = "10px Inter, system-ui, sans-serif";
    for (const hc of HEIGHT_COLORS) {
      if (hc.maxKm === Infinity && !windData.some((w) => w.hAgl > 12000)) continue;
      ctx.fillStyle = hc.color;
      ctx.fillRect(legX, legY, 12, 10);
      ctx.fillStyle = colors.text;
      ctx.textAlign = "left";
      ctx.textBaseline = "middle";
      ctx.fillText(hc.label, legX + 16, legY + 5);
      legY += 14;
    }
    // Storm motion legend
    [
      { color: "#ef4444", label: "RM" },
      { color: "#3b82f6", label: "LM" },
      { color: "#a3a3a3", label: "MW" },
    ].forEach(({ color, label }) => {
      ctx.fillStyle = color;
      ctx.beginPath();
      ctx.arc(legX + 6, legY + 5, 4, 0, Math.PI * 2);
      ctx.fill();
      ctx.fillStyle = colors.text;
      ctx.textAlign = "left";
      ctx.textBaseline = "middle";
      ctx.fillText(label, legX + 16, legY + 5);
      legY += 14;
    });
  }, [windData, dims, colors, scale, cx, cy, maxSpd, plotSize, params]);

  /* ── Crosshair overlay ── */
  useEffect(() => {
    const overlay = overlayRef.current;
    if (!overlay) return;
    const ctx = overlay.getContext("2d");
    const dpr = window.devicePixelRatio || 1;
    overlay.width = dims.w * dpr;
    overlay.height = dims.h * dpr;
    ctx.scale(dpr, dpr);
    ctx.clearRect(0, 0, dims.w, dims.h);
    if (!cursor || !showCrosshair) return;

    const { x, y } = cursor;
    // Find nearest data point
    let nearest = null;
    let minDist = Infinity;
    for (const w of windData) {
      const px = cx + scale(w.u);
      const py = cy - scale(w.v);
      const d = Math.hypot(px - x, py - y);
      if (d < minDist) { minDist = d; nearest = { ...w, px, py }; }
    }

    if (nearest && minDist < 40) {
      // Highlight point
      ctx.strokeStyle = colors.crosshair;
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.arc(nearest.px, nearest.py, 6, 0, Math.PI * 2);
      ctx.stroke();

      // Readout
      const hKm = (nearest.hAgl / 1000).toFixed(1);
      const lines = [
        `${nearest.wd}° / ${nearest.ws} kt`,
        `${hKm} km AGL`,
        `u=${nearest.u.toFixed(1)} v=${nearest.v.toFixed(1)}`,
      ];
      const boxW = 130, boxH = lines.length * 16 + 8;
      let bx = nearest.px + 12;
      let by = nearest.py - boxH / 2;
      if (bx + boxW > dims.w) bx = nearest.px - boxW - 12;
      if (by < 0) by = 4;
      if (by + boxH > dims.h) by = dims.h - boxH - 4;

      ctx.fillStyle = theme === "light" ? "rgba(255,255,255,0.92)" : "rgba(20,20,20,0.92)";
      ctx.strokeStyle = colors.crosshair;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.roundRect(bx, by, boxW, boxH, 4);
      ctx.fill();
      ctx.stroke();

      ctx.fillStyle = colors.crossText;
      ctx.font = "11px Inter, system-ui, sans-serif";
      ctx.textAlign = "left";
      ctx.textBaseline = "top";
      lines.forEach((line, i) => ctx.fillText(line, bx + 6, by + 4 + i * 16));
    }
  }, [cursor, dims, windData, showCrosshair, colors, scale, cx, cy, theme]);

  /* Mouse handlers */
  const handleMouseMove = useCallback((e) => {
    const rect = e.currentTarget.getBoundingClientRect();
    setCursor({ x: e.clientX - rect.left, y: e.clientY - rect.top });
  }, []);
  const handleMouseLeave = useCallback(() => setCursor(null), []);

  /* Download PNG */
  const handleDownload = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const link = document.createElement("a");
    link.download = "hodograph.png";
    link.href = canvas.toDataURL("image/png");
    link.click();
  }, []);

  if (!profile || profile.length === 0) return null;

  return (
    <div className="hodo-container" ref={containerRef}>
      <div className="hodo-toolbar">
        <button
          className={`hodo-btn${showCrosshair ? " active" : ""}`}
          onClick={() => setShowCrosshair((v) => !v)}
          title="Toggle crosshair"
        >
          <Crosshair size={16} />
        </button>
        <button className="hodo-btn" onClick={handleDownload} title="Download PNG">
          <Download size={16} />
        </button>
      </div>
      <div className="hodo-canvas-wrap">
        <canvas ref={canvasRef} style={{ width: dims.w, height: dims.h }} />
        <canvas
          ref={overlayRef}
          className="hodo-overlay"
          style={{ width: dims.w, height: dims.h }}
          onMouseMove={handleMouseMove}
          onMouseLeave={handleMouseLeave}
        />
      </div>
    </div>
  );
}
