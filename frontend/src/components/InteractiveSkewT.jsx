import { useRef, useEffect, useState, useCallback, useMemo } from "react";
import * as d3 from "d3";
import {
  ZoomIn,
  ZoomOut,
  RotateCcw,
  Crosshair,
  Download,
  Layers,
} from "lucide-react";
import "./InteractiveSkewT.css";

/* ── Constants ────────────────────────────────────────── */
const MARGIN = { top: 30, right: 50, bottom: 50, left: 55 };
const P_MIN = 100;   // hPa top
const P_MAX = 1050;  // hPa bottom
const T_MIN = -50;   // °C left
const T_MAX = 50;    // °C right
const SKEW_ANGLE = 45; // degrees

/* Pressure levels for isotherms / isobars */
const ISOBARS = [100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 850, 900, 925, 950, 1000];
const ISOTHERMS = d3.range(-80, 60, 10);
const DRY_ADIABATS = d3.range(-40, 60, 10);
const MOIST_ADIABATS = [4, 8, 12, 16, 20, 24, 28, 32];
const MIXING_RATIOS = [1, 2, 4, 7, 10, 16, 24];

/* ── Thermo helpers ───────────────────────────────────── */
function skewX(t, p, yScale, width) {
  const y = yScale(p);
  const yBottom = yScale(P_MAX);
  const yFromBottom = yBottom - y;
  const xBase = ((t - T_MIN) / (T_MAX - T_MIN)) * width;
  return xBase + yFromBottom; // 1:1 pixel ratio → true 45° visual angle
}

/** Saturation vapor pressure (Bolton 1980) in hPa */
function es(tC) {
  return 6.112 * Math.exp((17.67 * tC) / (tC + 243.5));
}

/** Saturation mixing ratio in g/kg */
function satMixingRatio(tC, pHpa) {
  const e = es(tC);
  return (621.97 * e) / (pHpa - e);
}

/** Potential temperature (K) */
function theta(tC, pHpa) {
  return (tC + 273.15) * Math.pow(1000 / pHpa, 0.286);
}

/** Temperature from theta and pressure */
function tFromTheta(th, pHpa) {
  return th / Math.pow(1000 / pHpa, 0.286) - 273.15;
}

/** Moist (pseudo) adiabat via stepwise integration from the surface upward */
function moistAdiabat(thetaW) {
  const result = [];
  let t = thetaW; // temperature at ~1000 hPa in °C
  // Integrate from surface upward (decreasing pressure)
  for (let pp = 1050; pp >= 100; pp -= 5) {
    result.push({ p: pp, t });
    const rs = satMixingRatio(t, pp) / 1000; // g/kg → kg/kg
    const Lv = 2501000 - 2370 * t;
    const cp = 1005.7;
    const Tk = t + 273.15;
    const dTdp = (287 * Tk / (pp * 100)) *
      (1 + (Lv * rs) / (287 * Tk)) /
      (1 + (Lv * Lv * rs) / (cp * 461.5 * Tk * Tk));
    t -= dTdp * 5 * 100; // going up: dp = −5 hPa
  }
  return result;
}

/** Mixing ratio line: temperature at which ws(t,p) = w */
function mixingRatioT(w, pHpa) {
  // Invert satMixingRatio: e = w*p / (621.97+w), then Bolton inversion
  const e = (w * pHpa) / (621.97 + w);
  if (e <= 0) return -999;
  return (243.5 * Math.log(e / 6.112)) / (17.67 - Math.log(e / 6.112));
}

/* ── Wind barb drawing ────────────────────────────────── */
function drawWindBarb(ctx, cx, cy, dir, spd, scale = 1) {
  if (dir == null || spd == null || spd < 0) return;
  ctx.save();
  ctx.translate(cx, cy);
  ctx.rotate(((dir + 180) * Math.PI) / 180);

  const len = 25 * scale;
  ctx.beginPath();
  ctx.moveTo(0, 0);
  ctx.lineTo(0, -len);
  ctx.stroke();

  let remaining = spd;
  let y = -len;
  const barbSpacing = 3.5 * scale;
  const barbLen = 8 * scale;

  // Pennants (50 kt)
  while (remaining >= 50) {
    ctx.beginPath();
    ctx.moveTo(0, y);
    ctx.lineTo(barbLen, y + barbSpacing);
    ctx.lineTo(0, y + barbSpacing * 2);
    ctx.closePath();
    ctx.fill();
    y += barbSpacing * 2 + 1;
    remaining -= 50;
  }
  // Full barbs (10 kt)
  while (remaining >= 10) {
    ctx.beginPath();
    ctx.moveTo(0, y);
    ctx.lineTo(barbLen, y + barbSpacing);
    ctx.stroke();
    y += barbSpacing;
    remaining -= 10;
  }
  // Half barb (5 kt)
  if (remaining >= 5) {
    ctx.beginPath();
    ctx.moveTo(0, y);
    ctx.lineTo(barbLen * 0.6, y + barbSpacing * 0.6);
    ctx.stroke();
  }
  ctx.restore();
}

/* ═══════════════════════════════════════════════════════ */
/*  InteractiveSkewT Component                           */
/* ═══════════════════════════════════════════════════════ */
export default function InteractiveSkewT({
  profile,
  sbParcel,
  mlParcel,
  params,
  theme = "dark",
}) {
  const canvasRef = useRef(null);
  const overlayRef = useRef(null);
  const containerRef = useRef(null);
  const [dims, setDims] = useState({ w: 800, h: 500 });
  const [cursor, setCursor] = useState(null); // { x, y } in canvas coords
  const [showCrosshair, setShowCrosshair] = useState(true);
  const [showParcel, setShowParcel] = useState(true);
  const [showBarbs, setShowBarbs] = useState(true);
  const [zoom, setZoom] = useState(1);
  const [panOffset, setPanOffset] = useState({ x: 0, y: 0 });
  const isPanning = useRef(false);
  const panStart = useRef({ x: 0, y: 0 });
  const lastPanOffset = useRef({ x: 0, y: 0 });

  /* Colours based on theme */
  const colors = useMemo(
    () =>
      theme === "light"
        ? {
            bg: "#ffffff",
            grid: "#e0e0e0",
            isobar: "#cccccc",
            isotherm: "#dddddd",
            isotherm0: "#60a5fa",
            dryAdiabat: "rgba(220,100,100,0.2)",
            moistAdiabat: "rgba(80,180,80,0.2)",
            mixRatio: "rgba(100,200,100,0.25)",
            temp: "#ef4444",
            dewpt: "#22c55e",
            sbParcel: "#f97316",
            mlParcel: "#d946ef",
            barb: "#555",
            text: "#333",
            crosshair: "rgba(100,100,100,0.5)",
            crossText: "#333",
            cape: "rgba(239,68,68,0.12)",
            cin: "rgba(59,130,246,0.12)",
          }
        : {
            bg: "#0d0d0d",
            grid: "#222",
            isobar: "#333",
            isotherm: "#2a2a2a",
            isotherm0: "#3b82f6",
            dryAdiabat: "rgba(200,80,80,0.15)",
            moistAdiabat: "rgba(80,200,80,0.15)",
            mixRatio: "rgba(80,200,100,0.12)",
            temp: "#ef4444",
            dewpt: "#22c55e",
            sbParcel: "#f97316",
            mlParcel: "#d946ef",
            barb: "#aaa",
            text: "#ccc",
            crosshair: "rgba(200,200,200,0.4)",
            crossText: "#eee",
            cape: "rgba(239,68,68,0.15)",
            cin: "rgba(59,130,246,0.15)",
          },
    [theme]
  );

  /* Responsive sizing */
  useEffect(() => {
    const el = containerRef.current;
    if (!el) return;
    const ro = new ResizeObserver((entries) => {
      const { width } = entries[0].contentRect;
      const w = Math.max(400, Math.floor(width));
      const h = Math.max(500, Math.floor(w * 0.9));
      setDims({ w, h });
    });
    ro.observe(el);
    return () => ro.disconnect();
  }, []);

  /* Create scales */
  const plotW = dims.w - MARGIN.left - MARGIN.right;
  const plotH = dims.h - MARGIN.top - MARGIN.bottom;

  const yScale = useCallback(
    (p) => {
      // Log-pressure scale
      const logP = Math.log(p);
      const logPmin = Math.log(P_MIN);
      const logPmax = Math.log(P_MAX);
      return MARGIN.top + ((logP - logPmin) / (logPmax - logPmin)) * plotH * zoom - panOffset.y;
    },
    [plotH, zoom, panOffset.y]
  );

  const yScaleInv = useCallback(
    (y) => {
      const logPmin = Math.log(P_MIN);
      const logPmax = Math.log(P_MAX);
      const logP = ((y + panOffset.y - MARGIN.top) / (plotH * zoom)) * (logPmax - logPmin) + logPmin;
      return Math.exp(logP);
    },
    [plotH, zoom, panOffset.y]
  );

  const toX = useCallback(
    (t, p) => MARGIN.left + skewX(t, p, yScale, plotW * zoom) - panOffset.x,
    [yScale, plotW, zoom, panOffset.x]
  );

  /* ── Main draw ──────────────────────────────────────── */
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas || !profile || profile.length === 0) return;
    const ctx = canvas.getContext("2d");
    const dpr = window.devicePixelRatio || 1;
    canvas.width = dims.w * dpr;
    canvas.height = dims.h * dpr;
    ctx.scale(dpr, dpr);

    // Clear
    ctx.fillStyle = colors.bg;
    ctx.fillRect(0, 0, dims.w, dims.h);

    // Clip to plot area
    ctx.save();
    ctx.beginPath();
    ctx.rect(MARGIN.left, MARGIN.top, plotW, plotH);
    ctx.clip();

    // ── Dry adiabats ──
    ctx.strokeStyle = colors.dryAdiabat;
    ctx.lineWidth = 0.8;
    for (const th of DRY_ADIABATS) {
      const thK = th + 273.15;
      ctx.beginPath();
      let first = true;
      for (let pp = P_MIN; pp <= P_MAX; pp += 10) {
        const t = tFromTheta(thK, pp);
        const x = toX(t, pp);
        const y = yScale(pp);
        if (first) { ctx.moveTo(x, y); first = false; }
        else ctx.lineTo(x, y);
      }
      ctx.stroke();
    }

    // ── Moist adiabats ──
    ctx.strokeStyle = colors.moistAdiabat;
    ctx.lineWidth = 0.8;
    for (const tw of MOIST_ADIABATS) {
      const pts = moistAdiabat(tw);
      ctx.beginPath();
      let first = true;
      for (const pt of pts) {
        if (pt.p < P_MIN || pt.p > P_MAX) continue;
        const x = toX(pt.t, pt.p);
        const y = yScale(pt.p);
        if (first) { ctx.moveTo(x, y); first = false; }
        else ctx.lineTo(x, y);
      }
      ctx.stroke();
    }

    // ── Mixing ratio lines ──
    ctx.strokeStyle = colors.mixRatio;
    ctx.lineWidth = 0.6;
    ctx.setLineDash([4, 6]);
    for (const w of MIXING_RATIOS) {
      ctx.beginPath();
      let first = true;
      for (let pp = 400; pp <= P_MAX; pp += 10) {
        const t = mixingRatioT(w, pp);
        if (t < -100) continue;
        const x = toX(t, pp);
        const y = yScale(pp);
        if (first) { ctx.moveTo(x, y); first = false; }
        else ctx.lineTo(x, y);
      }
      ctx.stroke();
    }
    ctx.setLineDash([]);

    // ── Isotherms ──
    ctx.lineWidth = 0.6;
    for (const t of ISOTHERMS) {
      ctx.strokeStyle = t === 0 ? colors.isotherm0 : colors.isotherm;
      ctx.lineWidth = t === 0 ? 1.5 : 0.6;
      ctx.beginPath();
      ctx.moveTo(toX(t, P_MIN), yScale(P_MIN));
      ctx.lineTo(toX(t, P_MAX), yScale(P_MAX));
      ctx.stroke();
    }

    // ── Isobars ──
    ctx.strokeStyle = colors.isobar;
    ctx.lineWidth = 0.5;
    for (const pp of ISOBARS) {
      const y = yScale(pp);
      ctx.beginPath();
      ctx.moveTo(MARGIN.left, y);
      ctx.lineTo(MARGIN.left + plotW, y);
      ctx.stroke();
    }

    // ── CAPE/CIN shading between T and SB parcel ──
    if (showParcel && sbParcel && profile.length === sbParcel.length) {
      for (let i = 0; i < profile.length - 1; i++) {
        const p0 = profile[i].p, p1 = profile[i + 1].p;
        const t0 = profile[i].t, t1 = profile[i + 1].t;
        const sp0 = sbParcel[i], sp1 = sbParcel[i + 1];
        if (p0 == null || p1 == null || t0 == null || t1 == null || sp0 == null || sp1 == null) continue;

        const isCAPE = sp0 > t0 || sp1 > t1;
        ctx.fillStyle = isCAPE ? colors.cape : colors.cin;
        ctx.beginPath();
        ctx.moveTo(toX(t0, p0), yScale(p0));
        ctx.lineTo(toX(t1, p1), yScale(p1));
        ctx.lineTo(toX(sp1, p1), yScale(p1));
        ctx.lineTo(toX(sp0, p0), yScale(p0));
        ctx.closePath();
        ctx.fill();
      }
    }

    // ── Temperature line ──
    ctx.strokeStyle = colors.temp;
    ctx.lineWidth = 2.2;
    ctx.beginPath();
    let started = false;
    for (const lev of profile) {
      if (lev.p == null || lev.t == null) continue;
      const x = toX(lev.t, lev.p);
      const y = yScale(lev.p);
      if (!started) { ctx.moveTo(x, y); started = true; }
      else ctx.lineTo(x, y);
    }
    ctx.stroke();

    // ── Dewpoint line ──
    ctx.strokeStyle = colors.dewpt;
    ctx.lineWidth = 2.2;
    ctx.beginPath();
    started = false;
    for (const lev of profile) {
      if (lev.p == null || lev.td == null) continue;
      const x = toX(lev.td, lev.p);
      const y = yScale(lev.p);
      if (!started) { ctx.moveTo(x, y); started = true; }
      else ctx.lineTo(x, y);
    }
    ctx.stroke();

    // ── SB Parcel trace ──
    if (showParcel && sbParcel) {
      ctx.strokeStyle = colors.sbParcel;
      ctx.lineWidth = 1.8;
      ctx.setLineDash([6, 3]);
      ctx.beginPath();
      started = false;
      for (let i = 0; i < profile.length && i < sbParcel.length; i++) {
        if (profile[i].p == null || sbParcel[i] == null) continue;
        const x = toX(sbParcel[i], profile[i].p);
        const y = yScale(profile[i].p);
        if (!started) { ctx.moveTo(x, y); started = true; }
        else ctx.lineTo(x, y);
      }
      ctx.stroke();
      ctx.setLineDash([]);
    }

    // ── ML Parcel trace ──
    if (showParcel && mlParcel) {
      ctx.strokeStyle = colors.mlParcel;
      ctx.lineWidth = 1.5;
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      started = false;
      for (let i = 0; i < profile.length && i < mlParcel.length; i++) {
        if (profile[i].p == null || mlParcel[i] == null) continue;
        const x = toX(mlParcel[i], profile[i].p);
        const y = yScale(profile[i].p);
        if (!started) { ctx.moveTo(x, y); started = true; }
        else ctx.lineTo(x, y);
      }
      ctx.stroke();
      ctx.setLineDash([]);
    }

    // ── Wind barbs ──
    if (showBarbs) {
      ctx.strokeStyle = colors.barb;
      ctx.fillStyle = colors.barb;
      ctx.lineWidth = 1.2;
      const barbX = MARGIN.left + plotW + 20;
      // Thin out barbs to avoid overlap
      let lastBarbY = -Infinity;
      for (const lev of profile) {
        if (lev.p == null || lev.ws == null || lev.wd == null) continue;
        const y = yScale(lev.p);
        if (y - lastBarbY < 18) continue;
        drawWindBarb(ctx, barbX, y, lev.wd, lev.ws, 0.7);
        lastBarbY = y;
      }
    }

    // ── Significant level markers (LCL, LFC, EL, FRZ, WB0) ──
    {
      // Helper: interpolate pressure from height AGL using profile data
      const heightToPressure = (hAgl) => {
        if (hAgl == null || !profile || profile.length < 2) return null;
        const sfcH = profile[0]?.h ?? 0;
        const targetH = sfcH + hAgl;
        for (let i = 0; i < profile.length - 1; i++) {
          const h0 = profile[i].h, h1 = profile[i + 1].h;
          const p0 = profile[i].p, p1 = profile[i + 1].p;
          if (h0 == null || h1 == null || p0 == null || p1 == null) continue;
          if (targetH >= h0 && targetH <= h1) {
            const frac = (targetH - h0) / (h1 - h0);
            return p0 + frac * (p1 - p0);
          }
        }
        return null;
      };

      const levels = [];
      // LCL (use SB LCL pressure — directly available)
      if (params?.sbLclP != null) levels.push({ p: params.sbLclP, label: "LCL", color: "#facc15" });
      // LFC (use SB LFC pressure)
      if (params?.sbLfcP != null) levels.push({ p: params.sbLfcP, label: "LFC", color: "#fb923c" });
      // EL  (use SB EL pressure)
      if (params?.sbElP != null) levels.push({ p: params.sbElP, label: "EL", color: "#c084fc" });
      // Freezing level (height → pressure interpolation)
      const frzP = heightToPressure(params?.frzLevel);
      if (frzP != null) levels.push({ p: frzP, label: "FRZ", color: "#60a5fa" });
      // Wet-bulb zero (height → pressure interpolation)
      const wboP = heightToPressure(params?.wbo);
      if (wboP != null) levels.push({ p: wboP, label: "WB0", color: "#38bdf8" });

      ctx.font = "bold 10px Inter, system-ui, sans-serif";
      for (const lev of levels) {
        const y = yScale(lev.p);
        if (y < MARGIN.top || y > MARGIN.top + plotH) continue;
        // Dashed horizontal line
        ctx.strokeStyle = lev.color;
        ctx.lineWidth = 1.2;
        ctx.setLineDash([6, 4]);
        ctx.beginPath();
        ctx.moveTo(MARGIN.left, y);
        ctx.lineTo(MARGIN.left + plotW, y);
        ctx.stroke();
        ctx.setLineDash([]);
        // Label on right edge
        ctx.fillStyle = lev.color;
        ctx.textAlign = "left";
        ctx.textBaseline = "middle";
        ctx.fillText(lev.label, MARGIN.left + plotW + 2, y);
      }
    }

    ctx.restore(); // un-clip

    // ── Axis labels ──
    ctx.fillStyle = colors.text;
    ctx.font = "11px Inter, system-ui, sans-serif";

    // Pressure axis (left)
    ctx.textAlign = "right";
    ctx.textBaseline = "middle";
    for (const pp of ISOBARS) {
      const y = yScale(pp);
      if (y < MARGIN.top || y > MARGIN.top + plotH) continue;
      ctx.fillText(`${pp}`, MARGIN.left - 6, y);
    }

    // Temperature axis (bottom)
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    for (const t of ISOTHERMS) {
      const x = toX(t, P_MAX);
      if (x < MARGIN.left || x > MARGIN.left + plotW) continue;
      ctx.fillText(`${t}°`, x, MARGIN.top + plotH + 6);
    }

    // Axis titles
    ctx.save();
    ctx.font = "12px Inter, system-ui, sans-serif";
    ctx.fillStyle = colors.text;
    ctx.textAlign = "center";
    ctx.fillText("Temperature (°C)", MARGIN.left + plotW / 2, dims.h - 8);
    ctx.translate(12, MARGIN.top + plotH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText("Pressure (hPa)", 0, 0);
    ctx.restore();

    // ── Legend ──
    const legendX = MARGIN.left + 8;
    const legendY = MARGIN.top + 8;
    ctx.font = "10px Inter, system-ui, sans-serif";
    const legendItems = [
      { color: colors.temp, label: "Temperature", dash: false },
      { color: colors.dewpt, label: "Dewpoint", dash: false },
    ];
    if (showParcel) {
      legendItems.push({ color: colors.sbParcel, label: "SB Parcel", dash: true });
      legendItems.push({ color: colors.mlParcel, label: "ML Parcel", dash: true });
    }
    legendItems.forEach((item, i) => {
      const ly = legendY + i * 16;
      ctx.strokeStyle = item.color;
      ctx.lineWidth = 2;
      if (item.dash) ctx.setLineDash([5, 3]);
      else ctx.setLineDash([]);
      ctx.beginPath();
      ctx.moveTo(legendX, ly + 5);
      ctx.lineTo(legendX + 20, ly + 5);
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = colors.text;
      ctx.textAlign = "left";
      ctx.textBaseline = "middle";
      ctx.fillText(item.label, legendX + 24, ly + 5);
    });
  }, [profile, sbParcel, mlParcel, params, dims, colors, toX, yScale, showParcel, showBarbs, zoom, panOffset]);

  /* ── Crosshair overlay ──────────────────────────────── */
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
    if (x < MARGIN.left || x > MARGIN.left + plotW || y < MARGIN.top || y > MARGIN.top + plotH) return;

    // Crosshair lines
    ctx.strokeStyle = colors.crosshair;
    ctx.lineWidth = 1;
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    ctx.moveTo(MARGIN.left, y);
    ctx.lineTo(MARGIN.left + plotW, y);
    ctx.moveTo(x, MARGIN.top);
    ctx.lineTo(x, MARGIN.top + plotH);
    ctx.stroke();
    ctx.setLineDash([]);

    // Read-out values
    const pAtCursor = yScaleInv(y);
    // Find nearest profile level
    let nearest = null;
    let minDist = Infinity;
    for (const lev of profile) {
      if (lev.p == null) continue;
      const dist = Math.abs(lev.p - pAtCursor);
      if (dist < minDist) { minDist = dist; nearest = lev; }
    }

    // Draw readout box
    const boxW = 180;
    const boxH = nearest ? 78 : 30;
    let bx = x + 12;
    let by = y - 10;
    if (bx + boxW > dims.w - 10) bx = x - boxW - 12;
    if (by + boxH > dims.h - 10) by = y - boxH;
    if (by < 5) by = 5;

    ctx.fillStyle = theme === "light" ? "rgba(255,255,255,0.92)" : "rgba(20,20,20,0.92)";
    ctx.strokeStyle = theme === "light" ? "#ccc" : "#444";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.roundRect(bx, by, boxW, boxH, 6);
    ctx.fill();
    ctx.stroke();

    ctx.font = "11px Inter, system-ui, sans-serif";
    ctx.fillStyle = colors.crossText;
    ctx.textAlign = "left";
    ctx.textBaseline = "top";
    ctx.fillText(`P: ${pAtCursor.toFixed(0)} hPa`, bx + 8, by + 6);

    if (nearest) {
      ctx.fillText(`T: ${nearest.t != null ? nearest.t + "°C" : "—"}  Td: ${nearest.td != null ? nearest.td + "°C" : "—"}`, bx + 8, by + 22);
      ctx.fillText(`H: ${nearest.h != null ? Math.round(nearest.h) + " m" : "—"}`, bx + 8, by + 38);
      const windStr = nearest.wd != null && nearest.ws != null
        ? `${Math.round(nearest.wd)}° / ${Math.round(nearest.ws)} kt`
        : "—";
      ctx.fillText(`Wind: ${windStr}`, bx + 8, by + 54);
    }
  }, [cursor, dims, profile, colors, showCrosshair, yScaleInv, plotW, plotH, theme]);

  /* ── Mouse handlers ─────────────────────────────────── */
  const handleMouseMove = useCallback(
    (e) => {
      const rect = overlayRef.current?.getBoundingClientRect();
      if (!rect) return;
      if (isPanning.current) {
        const dx = e.clientX - panStart.current.x;
        const dy = e.clientY - panStart.current.y;
        setPanOffset({
          x: lastPanOffset.current.x - dx,
          y: lastPanOffset.current.y - dy,
        });
        return;
      }
      setCursor({ x: e.clientX - rect.left, y: e.clientY - rect.top });
    },
    []
  );

  const handleMouseDown = useCallback((e) => {
    if (e.button === 0 && zoom > 1) {
      isPanning.current = true;
      panStart.current = { x: e.clientX, y: e.clientY };
      lastPanOffset.current = { ...panOffset };
    }
  }, [zoom, panOffset]);

  const handleMouseUp = useCallback(() => {
    isPanning.current = false;
  }, []);

  const handleMouseLeave = useCallback(() => {
    setCursor(null);
    isPanning.current = false;
  }, []);

  const handleWheel = useCallback(
    (e) => {
      e.preventDefault();
      const newZoom = Math.min(4, Math.max(1, zoom - e.deltaY * 0.002));
      setZoom(newZoom);
      if (newZoom <= 1) setPanOffset({ x: 0, y: 0 });
    },
    [zoom]
  );

  // Attach wheel listener imperatively with { passive: false } so
  // preventDefault() works (React's onWheel is passive by default).
  useEffect(() => {
    const el = overlayRef.current;
    if (!el) return;
    el.addEventListener("wheel", handleWheel, { passive: false });
    return () => el.removeEventListener("wheel", handleWheel);
  }, [handleWheel]);

  /* ── Download ──────────────────────────────────────── */
  const handleDownload = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const link = document.createElement("a");
    link.download = "skew-t_interactive.png";
    link.href = canvas.toDataURL("image/png");
    link.click();
  }, []);

  /* ── Reset zoom ─────────────────────────────────────── */
  const resetZoom = useCallback(() => {
    setZoom(1);
    setPanOffset({ x: 0, y: 0 });
  }, []);

  if (!profile || profile.length === 0) {
    return <div className="skewt-empty">No profile data available</div>;
  }

  return (
    <div className="skewt-container" ref={containerRef}>
      <div className="skewt-toolbar">
        <button
          className={`skewt-btn ${showCrosshair ? "active" : ""}`}
          onClick={() => setShowCrosshair((v) => !v)}
          title="Toggle crosshair"
        >
          <Crosshair size={14} />
        </button>
        <button
          className={`skewt-btn ${showParcel ? "active" : ""}`}
          onClick={() => setShowParcel((v) => !v)}
          title="Toggle parcel traces"
        >
          <Layers size={14} />
        </button>
        <button
          className="skewt-btn"
          onClick={() => setZoom((z) => Math.min(4, z + 0.3))}
          title="Zoom in"
        >
          <ZoomIn size={14} />
        </button>
        <button
          className="skewt-btn"
          onClick={() => {
            const nz = Math.max(1, zoom - 0.3);
            setZoom(nz);
            if (nz <= 1) setPanOffset({ x: 0, y: 0 });
          }}
          title="Zoom out"
        >
          <ZoomOut size={14} />
        </button>
        <button className="skewt-btn" onClick={resetZoom} title="Reset zoom">
          <RotateCcw size={14} />
        </button>
        <button className="skewt-btn" onClick={handleDownload} title="Download PNG">
          <Download size={14} />
        </button>
        {zoom > 1 && (
          <span className="skewt-zoom-label">{zoom.toFixed(1)}×</span>
        )}
      </div>
      <div className="skewt-canvas-wrap">
        <canvas
          ref={canvasRef}
          style={{ width: dims.w, height: dims.h }}
          className="skewt-canvas"
        />
        <canvas
          ref={overlayRef}
          style={{ width: dims.w, height: dims.h }}
          className="skewt-overlay"
          onMouseMove={handleMouseMove}
          onMouseDown={handleMouseDown}
          onMouseUp={handleMouseUp}
          onMouseLeave={handleMouseLeave}
        />
      </div>
    </div>
  );
}
