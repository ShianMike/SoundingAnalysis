import { useEffect } from "react";
import { useMap } from "react-leaflet";

/**
 * `FlyToCoords` — pans the map to a `{lat, lon}` pair whenever `coords`
 * changes (used by the chaser panel to focus on a selected report).
 * Holds at the current zoom unless it's < 9, in which case it zooms in.
 */
export function FlyToCoords({ coords, onDone }) {
  const map = useMap();
  useEffect(() => {
    if (!coords) return;
    map.flyTo([coords.lat, coords.lon], Math.max(map.getZoom(), 9), { duration: 0.5 });
    onDone?.();
  }, [coords]); // eslint-disable-line react-hooks/exhaustive-deps
  return null;
}

/**
 * `MapResizeHandler` — calls Leaflet's `invalidateSize()` shortly after
 * `trigger` changes. Use this whenever the surrounding container changes
 * dimensions (e.g. a sidebar collapses or fullscreen toggles) so Leaflet
 * recomputes its tile mask.
 */
export function MapResizeHandler({ trigger }) {
  const map = useMap();
  useEffect(() => {
    // Small delay lets the CSS transition finish before recalculating
    const id = setTimeout(() => map.invalidateSize(), 350);
    return () => clearTimeout(id);
  }, [trigger, map]);
  return null;
}
