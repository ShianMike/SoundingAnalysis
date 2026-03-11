import { useRef, useEffect, useState, useCallback, useMemo } from "react";
import { Download, Crosshair, Layers } from "lucide-react";
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

function streamwisenessColor(sw) {
  const hue = (1 - Math.max(0, Math.min(1, sw))) * 240;
  return `hsl(${hue}, 80%, 55%)`;
}

function interpWindAtHeight(data, targetH) {
  if (!data || data.length === 0) return null;
  for (let i = 0; i < data.length - 1; i++) {
    if (data[i].hAgl <= targetH && data[i + 1].hAgl >= targetH) {
      const frac = (targetH - data[i].hAgl) / (data[i + 1].hAgl - data[i].hAgl || 1);
      return {
        u: data[i].u + frac * (data[i + 1].u - data[i].u),
        v: data[i].v + frac * (data[i + 1].v - data[i].v),
      };
    }
  }
  if (targetH <= data[0].hAgl) return { u: data[0].u, v: data[0].v };
  return null;
}

function drawArrowhead(ctx, x, y, angle, len = 7) {
  ctx.beginPath();
  ctx.moveTo(x, y);
  ctx.lineTo(x - len * Math.cos(angle - 0.4), y - len * Math.sin(angle - 0.4));
  ctx.lineTo(x - len * Math.cos(angle + 0.4), y - len * Math.sin(angle + 0.4));
  ctx.closePath();
  ctx.fill();
}

export default function InteractiveHodograph({ profile, params, theme = "dark" }) {
  const canvasRef = useRef(null);
  const overlayRef = useRef(null);
  const containerRef = useRef(null);
  const [dims, setDims] = useState({ w: 500, h: 500 });
  const [cursor, setCursor] = useState(null);
  const [showCrosshair, setShowCrosshair] = useState(true);
  const [showStreamwise, setShowStreamwise] = useState(false);

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

    // ── Streamwiseness lookup ──
    const swArr = showStreamwise && params?.streamwiseness;
    const swH = showStreamwise && params?.streamwisenessHeight;
    const getStreamwise = (hKm) => {
      if (!swArr || !swH || swH.length === 0) return 0.5;
      if (hKm <= swH[0]) return swArr[0];
      if (hKm >= swH[swH.length - 1]) return swArr[swArr.length - 1];
      for (let j = 0; j < swH.length - 1; j++) {
        if (swH[j] <= hKm && swH[j + 1] >= hKm) {
          const f = (hKm - swH[j]) / (swH[j + 1] - swH[j] || 1);
          return swArr[j] + f * (swArr[j + 1] - swArr[j]);
        }
      }
      return 0.5;
    };

    // Hodograph trace — color-coded by height or streamwiseness
    ctx.lineWidth = 2.2;
    for (let i = 0; i < windData.length - 1; i++) {
      const w0 = windData[i], w1 = windData[i + 1];
      ctx.strokeStyle = swArr
        ? streamwisenessColor(getStreamwise((w0.hAgl + w1.hAgl) / 2000))
        : heightColor(w0.hAgl);
      ctx.beginPath();
      ctx.moveTo(cx + scale(w0.u), cy - scale(w0.v));
      ctx.lineTo(cx + scale(w1.u), cy - scale(w1.v));
      ctx.stroke();
    }

    // Dots at each data point
    for (const w of windData) {
      ctx.fillStyle = swArr
        ? streamwisenessColor(getStreamwise(w.hAgl / 1000))
        : heightColor(w.hAgl);
      ctx.beginPath();
      ctx.arc(cx + scale(w.u), cy - scale(w.v), 2, 0, Math.PI * 2);
      ctx.fill();
    }

    // ── Shear vectors (dashed arrows from surface to 1/3/6 km) ──
    const sfc = windData[0];
    const sfcPx = cx + scale(sfc.u);
    const sfcPy = cy - scale(sfc.v);
    const shearLayers = [
      { h: 1000, color: "#ef4444", label: "0\u20131" },
      { h: 3000, color: "#f97316", label: "0\u20133" },
      { h: 6000, color: "#22c55e", label: "0\u20136" },
    ];
    for (const sl of shearLayers) {
      const top = interpWindAtHeight(windData, sl.h);
      if (!top) continue;
      const topPx = cx + scale(top.u);
      const topPy = cy - scale(top.v);
      const mag = Math.round(Math.sqrt((top.u - sfc.u) ** 2 + (top.v - sfc.v) ** 2));
      ctx.save();
      ctx.globalAlpha = 0.6;
      ctx.setLineDash([5, 3]);
      ctx.strokeStyle = sl.color;
      ctx.lineWidth = 1.8;
      ctx.beginPath();
      ctx.moveTo(sfcPx, sfcPy);
      ctx.lineTo(topPx, topPy);
      ctx.stroke();
      ctx.setLineDash([]);
      const angle = Math.atan2(topPy - sfcPy, topPx - sfcPx);
      ctx.fillStyle = sl.color;
      drawArrowhead(ctx, topPx, topPy, angle);
      ctx.globalAlpha = 1;
      ctx.font = "bold 9px Inter, system-ui, sans-serif";
      ctx.fillStyle = sl.color;
      ctx.textAlign = "center";
      ctx.textBaseline = "bottom";
      ctx.fillText(`${sl.label}: ${mag}kt`, (sfcPx + topPx) / 2, (sfcPy + topPy) / 2 - 3);
      ctx.restore();
    }

    // ── Critical angle arc ──
    if (params?.criticalAngle != null && params.rmU != null && windData.length > 1) {
      const w500 = interpWindAtHeight(windData, 500);
      if (w500) {
        const shearAng = Math.atan2(-(w500.v - sfc.v), w500.u - sfc.u);
        const rmUKt = params.rmU * 1.94384;
        const rmVKt = params.rmV * 1.94384;
        const inflowAng = Math.atan2(-(sfc.v - rmVKt), sfc.u - rmUKt);
        let diff = inflowAng - shearAng;
        while (diff > Math.PI) diff -= 2 * Math.PI;
        while (diff < -Math.PI) diff += 2 * Math.PI;
        const arcR = 22;
        ctx.strokeStyle = "#fbbf24";
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        if (diff >= 0) ctx.arc(sfcPx, sfcPy, arcR, shearAng, shearAng + diff);
        else ctx.arc(sfcPx, sfcPy, arcR, shearAng + diff, shearAng);
        ctx.stroke();
        const midAng = shearAng + diff / 2;
        const labR = arcR + 11;
        ctx.font = "bold 9px Inter, system-ui, sans-serif";
        ctx.fillStyle = "#fbbf24";
        ctx.textAlign = "center";
        ctx.textBaseline = "middle";
        ctx.fillText(`${Math.round(params.criticalAngle)}\u00b0`, sfcPx + labR * Math.cos(midAng), sfcPy + labR * Math.sin(midAng));
      }
    }

    // ── Storm motion markers with direction/speed ──
    const drawStormMarker = (u, v, color, label, shape) => {
      if (u == null || v == null) return;
      const uKt = u * 1.94384;
      const vKt = v * 1.94384;
      const px = cx + scale(uKt);
      const py = cy - scale(vKt);
      const spd = Math.round(Math.sqrt(uKt * uKt + vKt * vKt));
      const heading = Math.round((90 - Math.atan2(vKt, uKt) * 180 / Math.PI + 360) % 360);
      ctx.strokeStyle = color;
      ctx.fillStyle = color;
      ctx.lineWidth = 2;
      if (shape === "cross") {
        const s = 5;
        ctx.beginPath();
        ctx.moveTo(px - s, py); ctx.lineTo(px + s, py);
        ctx.moveTo(px, py - s); ctx.lineTo(px, py + s);
        ctx.stroke();
      } else if (shape === "x") {
        const s = 4;
        ctx.beginPath();
        ctx.moveTo(px - s, py - s); ctx.lineTo(px + s, py + s);
        ctx.moveTo(px - s, py + s); ctx.lineTo(px + s, py - s);
        ctx.stroke();
      } else {
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.arc(px, py, 4, 0, Math.PI * 2);
        ctx.stroke();
      }
      ctx.font = "bold 10px Inter, system-ui, sans-serif";
      ctx.textAlign = "left";
      ctx.textBaseline = "middle";
      ctx.fillText(`${label} ${heading}\u00b0/${spd}kt`, px + 8, py);
    };

    if (params) {
      drawStormMarker(params.rmU, params.rmV, "#ef4444", "RM", "cross");
      drawStormMarker(params.lmU, params.lmV, "#3b82f6", "LM", "x");
      drawStormMarker(params.mwU, params.mwV, "#a3a3a3", "MW", "circle");
    }

    // ── Corfidi vectors ──
    if (params?.corfidiUpU != null && params?.corfidiUpV != null) {
      const px = cx + scale(params.corfidiUpU);
      const py = cy - scale(params.corfidiUpV);
      ctx.fillStyle = "#06b6d4";
      ctx.beginPath();
      ctx.moveTo(px, py - 5); ctx.lineTo(px - 4, py + 3); ctx.lineTo(px + 4, py + 3);
      ctx.closePath();
      ctx.fill();
      ctx.font = "bold 9px Inter, system-ui, sans-serif";
      ctx.textAlign = "left";
      ctx.textBaseline = "middle";
      ctx.fillText("CU", px + 7, py);
    }
    if (params?.corfidiDnU != null && params?.corfidiDnV != null) {
      const px = cx + scale(params.corfidiDnU);
      const py = cy - scale(params.corfidiDnV);
      ctx.fillStyle = "#8b5cf6";
      ctx.beginPath();
      ctx.moveTo(px, py + 5); ctx.lineTo(px - 4, py - 3); ctx.lineTo(px + 4, py - 3);
      ctx.closePath();
      ctx.fill();
      ctx.font = "bold 9px Inter, system-ui, sans-serif";
      ctx.textAlign = "left";
      ctx.textBaseline = "middle";
      ctx.fillText("CD", px + 7, py);
    }

    // ── Legend ──
    const legX = MARGIN.left + 6;
    let legY = MARGIN.top + 6;
    ctx.font = "10px Inter, system-ui, sans-serif";

    if (swArr) {
      // Streamwiseness gradient legend
      ctx.fillStyle = colors.text;
      ctx.textAlign = "left";
      ctx.textBaseline = "middle";
      ctx.fillText("Streamwiseness", legX, legY + 5);
      legY += 14;
      const gradH = 60;
      for (let i = 0; i < gradH; i++) {
        ctx.fillStyle = streamwisenessColor(1 - i / gradH);
        ctx.fillRect(legX, legY + i, 12, 1);
      }
      ctx.fillStyle = colors.text;
      ctx.font = "9px Inter, system-ui, sans-serif";
      ctx.textBaseline = "top";
      ctx.fillText("1.0 Streamwise", legX + 16, legY);
      ctx.textBaseline = "bottom";
      ctx.fillText("0.0 Crosswise", legX + 16, legY + gradH);
      legY += gradH + 6;
    } else {
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
    }

    // Storm motion legend
    legY += 2;
    ctx.font = "10px Inter, system-ui, sans-serif";
    [
      { color: "#ef4444", label: "RM", shape: "cross" },
      { color: "#3b82f6", label: "LM", shape: "x" },
      { color: "#a3a3a3", label: "MW", shape: "circle" },
    ].forEach(({ color, label, shape }) => {
      ctx.strokeStyle = color;
      ctx.fillStyle = color;
      ctx.lineWidth = 1.5;
      const mx = legX + 6, my = legY + 5;
      if (shape === "cross") {
        ctx.beginPath(); ctx.moveTo(mx - 4, my); ctx.lineTo(mx + 4, my);
        ctx.moveTo(mx, my - 4); ctx.lineTo(mx, my + 4); ctx.stroke();
      } else if (shape === "x") {
        ctx.beginPath(); ctx.moveTo(mx - 3, my - 3); ctx.lineTo(mx + 3, my + 3);
        ctx.moveTo(mx - 3, my + 3); ctx.lineTo(mx + 3, my - 3); ctx.stroke();
      } else {
        ctx.beginPath(); ctx.arc(mx, my, 3, 0, Math.PI * 2); ctx.stroke();
      }
      ctx.fillStyle = colors.text;
      ctx.textAlign = "left";
      ctx.textBaseline = "middle";
      ctx.fillText(label, legX + 16, legY + 5);
      legY += 14;
    });

    // Corfidi legend
    if (params?.corfidiUpU != null) {
      ctx.fillStyle = "#06b6d4";
      ctx.beginPath();
      ctx.moveTo(legX + 6, legY + 2); ctx.lineTo(legX + 2, legY + 8); ctx.lineTo(legX + 10, legY + 8);
      ctx.closePath(); ctx.fill();
      ctx.fillStyle = colors.text;
      ctx.fillText("CU", legX + 16, legY + 5);
      legY += 14;
    }
    if (params?.corfidiDnU != null) {
      ctx.fillStyle = "#8b5cf6";
      ctx.beginPath();
      ctx.moveTo(legX + 6, legY + 8); ctx.lineTo(legX + 2, legY + 2); ctx.lineTo(legX + 10, legY + 2);
      ctx.closePath(); ctx.fill();
      ctx.fillStyle = colors.text;
      ctx.fillText("CD", legX + 16, legY + 5);
      legY += 14;
    }

    // Critical angle legend
    if (params?.criticalAngle != null) {
      ctx.strokeStyle = "#fbbf24";
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.arc(legX + 6, legY + 5, 5, -0.5, 0.5);
      ctx.stroke();
      ctx.fillStyle = colors.text;
      ctx.fillText(`CA ${Math.round(params.criticalAngle)}\u00b0`, legX + 16, legY + 5);
      legY += 14;
    }
  }, [windData, dims, colors, scale, cx, cy, maxSpd, plotSize, params, showStreamwise]);

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
        <button
          className={`hodo-btn${showStreamwise ? " active" : ""}`}
          onClick={() => setShowStreamwise((v) => !v)}
          title="Toggle streamwiseness coloring"
        >
          <Layers size={16} />
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
