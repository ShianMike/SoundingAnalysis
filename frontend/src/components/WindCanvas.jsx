/**
 * WindCanvas – animated wind-particle overlay for a Leaflet map.
 *
 * Renders lightweight particles that flow along a gridded 10 m
 * wind field (U/V components).  The canvas sits on top of the map
 * container (pointer-events: none) and re-samples the wind every frame
 * via bilinear interpolation, so it adjusts automatically as the user
 * pans or zooms.
 *
 * Optimisations:
 *  • map.getBounds() called ONCE per frame, not per-particle
 *  • Pre-allocated Float64Array draw buffers — zero per-frame GC
 *  • Inline particle reset avoids Object.assign + ephemeral objects
 *  • Bucket draw uses 3 batched paths instead of per-particle styles
 *  • `desynchronized` canvas context skips compositor sync
 *  • Clears trails on movestart + zoomstart for correct panning
 */
import { useEffect, useRef } from "react";
import { useMap } from "react-leaflet";

/* ── tunables ────────────────────────────────────────────────── */
const N_PARTICLES  = 3500;
const MAX_AGE      = 100;     // frames before particle resets
const FADE         = 0.97;    // trail fade per frame (0-1)
const LINE_W       = 0.75;
const BASE_SCALE   = 0.15;    // pixels per (m/s) at zoom 3
const SPD_THRESH_LO = 0.2;   // ignore wind below this
const SPD_SQ_SLOW   = 25;    // 5² m/s — bucket boundary
const SPD_SQ_MED    = 144;   // 12² m/s — bucket boundary

/* bucket styling — defined once, avoids per-frame string creation */
const BUCKET_STYLES = [
  "rgba(140,180,255,0.4)",    // slow  (< 5 m/s)
  "rgba(200,220,255,0.55)",   // med   (5-12 m/s)
  "rgba(255,220,120,0.65)",   // fast  (> 12 m/s)
];

export default function WindCanvas({ data }) {
  const map = useMap();
  const canvasRef    = useRef(null);
  const rafRef       = useRef(null);
  const particlesRef = useRef(null);
  /* pre-allocated draw buffers: 4 floats per line (x1,y1,x2,y2) */
  const bufRef       = useRef(null);

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
    let W, H;
    const resize = () => {
      W = container.clientWidth;
      H = container.clientHeight;
      canvas.width  = W;
      canvas.height = H;
    };
    resize();
    canvas.style.display = "block";

    /* grid metadata — cache locally for hot-path access */
    const { nx, ny, la1, lo1, dx, dy } = data;
    const uGrid = data.uData;
    const vGrid = data.vData;
    const nxM1  = nx - 1;
    const nyM1  = ny - 1;

    /* particle array — flat typed struct-of-arrays for cache locality */
    const px  = new Float64Array(N_PARTICLES);
    const py  = new Float64Array(N_PARTICLES);
    const age = new Int32Array(N_PARTICLES);
    const maxAge = new Int32Array(N_PARTICLES);

    const resetParticle = (i) => {
      px[i]     = Math.random() * W;
      py[i]     = Math.random() * H;
      age[i]    = (Math.random() * MAX_AGE) | 0;
      maxAge[i] = MAX_AGE + ((Math.random() * 40) | 0);
    };
    for (let i = 0; i < N_PARTICLES; i++) resetParticle(i);
    particlesRef.current = { px, py, age, maxAge };

    /* pre-allocate draw buffers — 3 buckets × 4 floats × N_PARTICLES max */
    const maxLines = N_PARTICLES * 4;
    if (!bufRef.current || bufRef.current[0].length < maxLines) {
      bufRef.current = [
        new Float64Array(maxLines),
        new Float64Array(maxLines),
        new Float64Array(maxLines),
      ];
    }
    const bufs = bufRef.current;
    const bufLen = [0, 0, 0];   // current length per bucket

    const ctx = canvas.getContext("2d", { alpha: true, desynchronized: true });

    /* ── bilinear wind lookup (bounds cached per frame) ──── */
    let bWest, bEast, bNorth, bSouth, invW, invH;
    const cacheBounds = () => {
      const b = map.getBounds();
      bWest  = b.getWest();
      bEast  = b.getEast();
      bNorth = b.getNorth();
      bSouth = b.getSouth();
      invW   = (bEast  - bWest)  / W;
      invH   = (bNorth - bSouth) / H;
    };

    /** Returns [u, v] via bilinear interpolation, or null if off-grid. */
    const windAt = (ppx, ppy) => {
      const lon = bWest + ppx * invW;
      const lat = bNorth - ppy * invH;

      const gi = (lon - lo1) / dx;
      const gj = (lat - la1) / dy;
      const i0 = gi | 0, j0 = gj | 0;          // fast floor for positive values
      if (i0 < 0 || i0 >= nxM1 || j0 < 0 || j0 >= nyM1) return null;

      const fi = gi - i0, fj = gj - j0;
      const w00 = (1 - fi) * (1 - fj);
      const w10 = fi       * (1 - fj);
      const w01 = (1 - fi) * fj;
      const w11 = fi       * fj;

      const base = j0 * nx + i0;
      const row2 = base + nx;

      return [
        w00 * uGrid[base] + w10 * uGrid[base + 1] + w01 * uGrid[row2] + w11 * uGrid[row2 + 1],
        w00 * vGrid[base] + w10 * vGrid[base + 1] + w01 * vGrid[row2] + w11 * vGrid[row2 + 1],
      ];
    };

    /* ── animation frame ───────────────────────────────────── */
    const frame = () => {
      const zoom  = map.getZoom();
      const scale = BASE_SCALE * (1 + (zoom - 3) * 0.3);

      /* cache bounds once per frame */
      cacheBounds();

      /* fade existing trails */
      ctx.globalCompositeOperation = "destination-in";
      ctx.fillStyle = `rgba(0,0,0,${FADE})`;
      ctx.fillRect(0, 0, W, H);
      ctx.globalCompositeOperation = "source-over";

      /* reset bucket counters */
      bufLen[0] = bufLen[1] = bufLen[2] = 0;

      for (let i = 0; i < N_PARTICLES; i++) {
        if (age[i] >= maxAge[i]) { resetParticle(i); continue; }

        const w = windAt(px[i], py[i]);
        if (!w) { age[i] = maxAge[i]; continue; }

        const u = w[0], v = w[1];
        const spdSq = u * u + v * v;
        if (spdSq < SPD_THRESH_LO * SPD_THRESH_LO) { age[i]++; continue; }

        const ddx = u * scale;
        const ddy = -v * scale;
        const nx2 = px[i] + ddx;
        const ny2 = py[i] + ddy;

        /* classify into speed bucket (0=slow, 1=med, 2=fast) */
        const bk = spdSq < SPD_SQ_SLOW ? 0 : spdSq < SPD_SQ_MED ? 1 : 2;
        const off = bufLen[bk];
        bufs[bk][off]     = px[i];
        bufs[bk][off + 1] = py[i];
        bufs[bk][off + 2] = nx2;
        bufs[bk][off + 3] = ny2;
        bufLen[bk] = off + 4;

        px[i]  = nx2;
        py[i]  = ny2;
        age[i]++;

        if (nx2 < 0 || nx2 > W || ny2 < 0 || ny2 > H) resetParticle(i);
      }

      /* draw each bucket in a single batched path */
      ctx.lineWidth = LINE_W;
      for (let bk = 0; bk < 3; bk++) {
        const len = bufLen[bk];
        if (len === 0) continue;
        const buf = bufs[bk];
        ctx.strokeStyle = BUCKET_STYLES[bk];
        ctx.beginPath();
        for (let k = 0; k < len; k += 4) {
          ctx.moveTo(buf[k], buf[k + 1]);
          ctx.lineTo(buf[k + 2], buf[k + 3]);
        }
        ctx.stroke();
      }

      rafRef.current = requestAnimationFrame(frame);
    };

    frame();

    /* ── map event handlers ────────────────────────────────── */
    const clearTrails = () => ctx.clearRect(0, 0, W, H);
    const onResize    = () => { resize(); clearTrails(); };

    map.on("movestart", clearTrails);
    map.on("zoomstart", clearTrails);
    map.on("resize", onResize);

    return () => {
      cancelAnimationFrame(rafRef.current);
      map.off("movestart", clearTrails);
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
