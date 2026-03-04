/**
 * WindCanvas – animated wind-particle overlay for a Leaflet map.
 *
 * Renders ~3 500 lightweight particles that flow along a gridded 10 m
 * wind field (U/V components).  The canvas sits on top of the map
 * container (pointer-events: none) and re-samples the wind every frame
 * via bilinear interpolation, so it adjusts automatically as the user
 * pans or zooms.
 */
import { useEffect, useRef } from "react";
import { useMap } from "react-leaflet";

/* ── tunables ────────────────────────────────────────────────── */
const N_PARTICLES  = 3500;
const MAX_AGE      = 100;     // frames before particle resets
const FADE         = 0.97;    // trail fade per frame (0-1)
const LINE_W       = 0.75;
const BASE_SCALE   = 0.15;    // pixels per (m/s) at zoom 3
const COLOR_SLOW   = [140, 180, 255];  // light blue
const COLOR_FAST   = [255, 220, 120];  // warm yellow

export default function WindCanvas({ data }) {
  const map = useMap();
  const canvasRef   = useRef(null);
  const rafRef      = useRef(null);
  const particlesRef = useRef([]);

  /* ── main effect: create canvas, start animation ──────────── */
  useEffect(() => {
    if (!data) {
      cancelAnimationFrame(rafRef.current);
      if (canvasRef.current) canvasRef.current.style.display = "none";
      return;
    }

    const container = map.getContainer();
    let canvas = canvasRef.current;
    if (!canvas) {
      canvas = document.createElement("canvas");
      canvas.className = "wind-overlay-canvas";
      Object.assign(canvas.style, {
        position: "absolute",
        top: "0",
        left: "0",
        pointerEvents: "none",
        zIndex: "450",
      });
      container.appendChild(canvas);
      canvasRef.current = canvas;
    }

    /* sizing */
    const resize = () => {
      canvas.width  = container.clientWidth;
      canvas.height = container.clientHeight;
    };
    resize();
    canvas.style.display = "block";

    /* grid metadata */
    const { nx, ny, la1, lo1, dx, dy, uData, vData } = data;

    /* particle factory */
    const mkP = () => ({
      x:   Math.random() * canvas.width,
      y:   Math.random() * canvas.height,
      age: Math.floor(Math.random() * MAX_AGE),
      max: MAX_AGE + Math.floor(Math.random() * 40),
    });
    particlesRef.current = Array.from({ length: N_PARTICLES }, mkP);

    const ctx = canvas.getContext("2d");

    /* ── bilinear wind lookup ──────────────────────────────── */
    function windAt(px, py) {
      const bounds = map.getBounds();
      const west  = bounds.getWest(),  east  = bounds.getEast();
      const north = bounds.getNorth(), south = bounds.getSouth();
      const W = canvas.width, H = canvas.height;

      const lon = west  + (px / W) * (east  - west);
      const lat = north - (py / H) * (north - south);

      const gi = (lon - lo1) / dx;
      const gj = (lat - la1) / dy;
      const i0 = Math.floor(gi), j0 = Math.floor(gj);
      if (i0 < 0 || i0 + 1 >= nx || j0 < 0 || j0 + 1 >= ny) return null;

      const fi = gi - i0, fj = gj - j0;
      const idx = (j, i) => j * nx + i;

      const u =
        (1 - fi) * (1 - fj) * uData[idx(j0, i0)]     +
        fi       * (1 - fj) * uData[idx(j0, i0 + 1)] +
        (1 - fi) * fj       * uData[idx(j0 + 1, i0)] +
        fi       * fj       * uData[idx(j0 + 1, i0 + 1)];

      const v =
        (1 - fi) * (1 - fj) * vData[idx(j0, i0)]     +
        fi       * (1 - fj) * vData[idx(j0, i0 + 1)] +
        (1 - fi) * fj       * vData[idx(j0 + 1, i0)] +
        fi       * fj       * vData[idx(j0 + 1, i0 + 1)];

      return [u, v];
    }

    /* ── colour by speed ───────────────────────────────────── */
    function colorForSpeed(spd) {
      const t = Math.min(spd / 18, 1);           // 0-18 m/s range
      const r = Math.round(COLOR_SLOW[0] + t * (COLOR_FAST[0] - COLOR_SLOW[0]));
      const g = Math.round(COLOR_SLOW[1] + t * (COLOR_FAST[1] - COLOR_SLOW[1]));
      const b = Math.round(COLOR_SLOW[2] + t * (COLOR_FAST[2] - COLOR_SLOW[2]));
      return `rgba(${r},${g},${b},0.55)`;
    }

    /* ── animation frame ───────────────────────────────────── */
    function frame() {
      const zoom  = map.getZoom();
      const scale = BASE_SCALE * (1 + (zoom - 3) * 0.3);

      /* fade existing trails */
      ctx.globalCompositeOperation = "destination-in";
      ctx.fillStyle = `rgba(0,0,0,${FADE})`;
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      ctx.globalCompositeOperation = "source-over";

      /* batch particles by speed bucket for fewer style switches */
      const buckets = { slow: [], med: [], fast: [] };

      for (const p of particlesRef.current) {
        if (p.age >= p.max) { Object.assign(p, mkP()); continue; }

        const w = windAt(p.x, p.y);
        if (!w) { p.age = p.max; continue; }

        const [u, v] = w;
        const spd = Math.sqrt(u * u + v * v);
        if (spd < 0.2) { p.age++; continue; }

        const ddx = u * scale;
        const ddy = -v * scale;            // flip Y for screen coords
        const nx  = p.x + ddx;
        const ny  = p.y + ddy;

        const bucket = spd < 5 ? "slow" : spd < 12 ? "med" : "fast";
        buckets[bucket].push(p.x, p.y, nx, ny);

        p.x = nx;
        p.y = ny;
        p.age++;

        if (nx < 0 || nx > canvas.width || ny < 0 || ny > canvas.height) {
          Object.assign(p, mkP());
        }
      }

      /* draw each bucket */
      ctx.lineWidth = LINE_W;
      for (const [key, arr] of Object.entries(buckets)) {
        if (arr.length === 0) continue;
        ctx.strokeStyle =
          key === "slow" ? "rgba(140,180,255,0.4)" :
          key === "med"  ? "rgba(200,220,255,0.55)" :
                           "rgba(255,220,120,0.65)";
        ctx.beginPath();
        for (let k = 0; k < arr.length; k += 4) {
          ctx.moveTo(arr[k], arr[k + 1]);
          ctx.lineTo(arr[k + 2], arr[k + 3]);
        }
        ctx.stroke();
      }

      rafRef.current = requestAnimationFrame(frame);
    }

    frame();

    /* ── map event handlers ────────────────────────────────── */
    const clearTrails = () => ctx.clearRect(0, 0, canvas.width, canvas.height);
    const onResize    = () => { resize(); clearTrails(); };

    map.on("zoomstart", clearTrails);
    map.on("resize", onResize);

    return () => {
      cancelAnimationFrame(rafRef.current);
      map.off("zoomstart", clearTrails);
      map.off("resize", onResize);
    };
  }, [data, map]);

  /* ── cleanup on unmount ──────────────────────────────────── */
  useEffect(() => {
    return () => {
      cancelAnimationFrame(rafRef.current);
      if (canvasRef.current) {
        canvasRef.current.remove();
        canvasRef.current = null;
      }
    };
  }, []);

  return null;
}
