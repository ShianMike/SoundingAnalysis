/**
 * WindCanvas – animated wind-particle overlay for a Leaflet map.
 *
 * Particles live in geographic (lat/lon) space and project to screen
 * pixels each frame via a cached Mercator mapping.
 *
 * The canvas is inserted into `.leaflet-map-pane` with z-index 250
 * (between tile-pane at 200 and overlay-pane at 400), so wind renders
 * behind SPC outlooks, radar, and markers.  We counter-translate for
 * the map-pane's CSS transform each frame.
 */
import { useEffect, useRef } from "react";
import { useMap } from "react-leaflet";

/* ── tunables ────────────────────────────────────────────────── */
const N_PARTICLES   = 800;
const MAX_AGE       = 100;
const FADE          = 0.96;
const LINE_W        = 0.7;
const BASE_SPEED    = 0.012;
const SPD_THRESH_SQ = 0.04;
const SPD_SQ_SLOW   = 25;
const SPD_SQ_MED    = 144;
const DEG2RAD       = Math.PI / 180;
const mercY = (lat) => Math.log(Math.tan(Math.PI * 0.25 + lat * DEG2RAD * 0.5));

const BUCKET_STYLES = [
  "rgba(140,180,255,0.25)",
  "rgba(200,220,255,0.35)",
  "rgba(255,220,120,0.45)",
];

export default function WindCanvas({ data }) {
  const map = useMap();
  const canvasRef = useRef(null);
  const rafRef    = useRef(null);
  const bufRef    = useRef(null);

  useEffect(() => {
    if (!data) {
      cancelAnimationFrame(rafRef.current);
      if (canvasRef.current) canvasRef.current.style.display = "none";
      return;
    }

    const isUpper  = data.level && data.level !== "surface";
    const spdSqSlow = isUpper ? 100 : SPD_SQ_SLOW;
    const spdSqMed  = isUpper ? 625 : SPD_SQ_MED;

    const container = map.getContainer();
    const mapPane   = map.getPane("mapPane");          // .leaflet-map-pane

    let canvas = canvasRef.current;
    if (!canvas) {
      canvas = document.createElement("canvas");
      canvas.className = "wind-overlay-canvas";
      Object.assign(canvas.style, {
        position: "absolute", top: "0", left: "0",
        pointerEvents: "none",
        /* Between tile-pane (200) and overlay-pane (400) */
        zIndex: "250",
      });
      mapPane.appendChild(canvas);
      canvasRef.current = canvas;
    }

    let W, H;
    const resize = () => {
      W = container.clientWidth;
      H = container.clientHeight;
      canvas.width  = W;
      canvas.height = H;
    };
    resize();
    canvas.style.display = "block";

    const { nx, ny, la1, lo1, dx, dy } = data;
    const uGrid = data.uData;
    const vGrid = data.vData;
    const nxM1 = nx - 1, nyM1 = ny - 1;
    const la2 = data.la2, lo2 = data.lo2;

    /* ── Particles in GEOGRAPHIC coords (lat, lon) ─────────── */
    const pLat  = new Float64Array(N_PARTICLES);
    const pLon  = new Float64Array(N_PARTICLES);
    const age   = new Int32Array(N_PARTICLES);
    const mxAge = new Int32Array(N_PARTICLES);

    const randomLat = () => la1 + Math.random() * (la2 - la1);
    const randomLon = () => lo1 + Math.random() * (lo2 - lo1);

    const resetParticle = (i) => {
      // Scatter within current map bounds intersected with the data domain
      const b = map.getBounds();
      const latLo = Math.max(la1, b.getSouth());
      const latHi = Math.min(la2, b.getNorth());
      const lonLo = Math.max(lo1, b.getWest());
      const lonHi = Math.min(lo2, b.getEast());
      if (latLo >= latHi || lonLo >= lonHi) {
        // Viewport is outside data domain — place anywhere in domain
        pLat[i] = randomLat();
        pLon[i] = randomLon();
      } else {
        pLat[i] = latLo + Math.random() * (latHi - latLo);
        pLon[i] = lonLo + Math.random() * (lonHi - lonLo);
      }
      age[i]   = (Math.random() * MAX_AGE) | 0;
      mxAge[i] = MAX_AGE + ((Math.random() * 40) | 0);
    };
    const resetAll = () => { for (let i = 0; i < N_PARTICLES; i++) resetParticle(i); };
    resetAll();

    /* draw buffers */
    const maxBuf = N_PARTICLES * 4;
    if (!bufRef.current || bufRef.current[0].length < maxBuf) {
      bufRef.current = [new Float64Array(maxBuf), new Float64Array(maxBuf), new Float64Array(maxBuf)];
    }
    const bufs = bufRef.current;
    const bufLen = [0, 0, 0];

    const ctx = canvas.getContext("2d", { alpha: true, desynchronized: true });

    /* ── bilinear wind lookup in lat/lon space ─────────────── */
    const windAt = (lat, lon) => {
      const gi = (lon - lo1) / dx;
      const gj = (lat - la1) / dy;
      const i0 = gi | 0, j0 = gj | 0;
      if (i0 < 0 || i0 >= nxM1 || j0 < 0 || j0 >= nyM1) return null;

      const fi = gi - i0, fj = gj - j0;
      const w00 = (1 - fi) * (1 - fj);
      const w10 = fi * (1 - fj);
      const w01 = (1 - fi) * fj;
      const w11 = fi * fj;
      const base = j0 * nx + i0;
      const row2 = base + nx;

      return [
        w00 * uGrid[base] + w10 * uGrid[base + 1] + w01 * uGrid[row2] + w11 * uGrid[row2 + 1],
        w00 * vGrid[base] + w10 * vGrid[base + 1] + w01 * vGrid[row2] + w11 * vGrid[row2 + 1],
      ];
    };

    /* ── Mercator projection (cached per frame) ──────────────── */
    let projOriginX, projScaleX, mercOY, mercSY;
    const cacheProjection = () => {
      const b = map.getBounds();
      const nw = map.latLngToContainerPoint(b.getNorthWest());
      const se = map.latLngToContainerPoint(b.getSouthEast());
      const west = b.getWest(), east = b.getEast();
      // X – linear in longitude
      projScaleX = (se.x - nw.x) / (east - west);
      projOriginX = nw.x - west * projScaleX;
      // Y – linear in Mercator-Y (correct Web Mercator projection)
      const mN = mercY(b.getNorth()), mS = mercY(b.getSouth());
      mercSY  = (se.y - nw.y) / (mS - mN);
      mercOY  = nw.y - mN * mercSY;
    };
    const toScreenX = (lon) => projOriginX + lon * projScaleX;
    const toScreenY = (lat) => mercOY + mercY(lat) * mercSY;

    /* ── animation frame ───────────────────────────────────── */
    const frame = () => {
      const zoom = map.getZoom();
      const speed = BASE_SPEED / Math.pow(2, (zoom - 4) * 0.5);

      /* Counteract the mapPane's CSS transform so the canvas stays
         aligned with the container (lat/lon→containerPoint coords). */
      const t = mapPane.style.transform;
      let px = 0, py = 0;
      let m = t.match(/translate3d\(([^,]+),\s*([^,]+)/);
      if (!m) m = t.match(/translate\(([^,]+),\s*([^,]+)/);
      if (m) { px = parseFloat(m[1]); py = parseFloat(m[2]); }
      canvas.style.transform = `translate(${-px}px, ${-py}px)`;

      cacheProjection();

      /* fade */
      ctx.globalCompositeOperation = "destination-in";
      ctx.fillStyle = `rgba(0,0,0,${FADE})`;
      ctx.fillRect(0, 0, W, H);
      ctx.globalCompositeOperation = "source-over";

      bufLen[0] = bufLen[1] = bufLen[2] = 0;

      for (let i = 0; i < N_PARTICLES; i++) {
        if (age[i] >= mxAge[i]) { resetParticle(i); continue; }

        const w = windAt(pLat[i], pLon[i]);
        if (!w) { age[i] = mxAge[i]; continue; }

        const u = w[0], v = w[1];
        const spdSq = u * u + v * v;
        if (spdSq < SPD_THRESH_SQ) { age[i]++; continue; }

        // Screen positions BEFORE move
        const sx1 = toScreenX(pLon[i]);
        const sy1 = toScreenY(pLat[i]);

        // Advance in geographic space (cos correction for longitude)
        const cosLat = Math.cos(pLat[i] * DEG2RAD);
        pLon[i] += u * speed / cosLat;
        pLat[i] += v * speed;
        age[i]++;

        // Screen positions AFTER move
        const sx2 = toScreenX(pLon[i]);
        const sy2 = toScreenY(pLat[i]);

        // Skip if both endpoints are off-screen
        if ((sx1 < -20 && sx2 < -20) || (sx1 > W + 20 && sx2 > W + 20) ||
            (sy1 < -20 && sy2 < -20) || (sy1 > H + 20 && sy2 > H + 20)) {
          resetParticle(i);
          continue;
        }

        // Bucket by speed (thresholds adapt to level)
        const bk = spdSq < spdSqSlow ? 0 : spdSq < spdSqMed ? 1 : 2;
        const off = bufLen[bk];
        bufs[bk][off]     = sx1;
        bufs[bk][off + 1] = sy1;
        bufs[bk][off + 2] = sx2;
        bufs[bk][off + 3] = sy2;
        bufLen[bk] = off + 4;

        // If particle left data domain, reset
        if (pLat[i] < la1 || pLat[i] > la2 || pLon[i] < lo1 || pLon[i] > lo2) {
          resetParticle(i);
        }
      }

      /* draw batched paths */
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

    /* ── map events ────────────────────────────────────────── */
    const clearTrails = () => ctx.clearRect(0, 0, W, H);
    const onResize = () => { resize(); clearTrails(); };
    // On zoom end, just clear old pixel trails — particles are in lat/lon, so
    // they'll naturally project to the correct new screen positions next frame.
    const onZoomEnd = () => { clearTrails(); };
    const onMoveEnd = () => { clearTrails(); };

    map.on("zoomstart", clearTrails);
    map.on("movestart", clearTrails);
    map.on("zoomend", onZoomEnd);
    map.on("moveend", onMoveEnd);
    map.on("resize", onResize);

    return () => {
      cancelAnimationFrame(rafRef.current);
      map.off("zoomstart", clearTrails);
      map.off("movestart", clearTrails);
      map.off("zoomend", onZoomEnd);
      map.off("moveend", onMoveEnd);
      map.off("resize", onResize);
    };
  }, [data, map]);

  /* cleanup on unmount */
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
