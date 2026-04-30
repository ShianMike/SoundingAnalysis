import { useCallback, useEffect } from "react";
import { useMap, useMapEvents } from "react-leaflet";

/**
 * `MapViewTracker` — invisible component that surfaces map center +
 * viewport bounds back to the parent via callback props.
 *
 * Used by the parent for:
 *   - nearest-NEXRAD lookup (which radar's velocity to show)
 *   - decimation (skip points outside the visible bbox)
 *
 * Fires `onCenterChange(lat, lng)` and `onViewChange({zoom, south, west, north, east})`
 * on mount and on every `moveend` / `zoomend`.
 */
export default function MapViewTracker({ onCenterChange, onViewChange }) {
  const map = useMap();

  const report = useCallback(() => {
    const c = map.getCenter();
    const b = map.getBounds();
    onCenterChange?.(c.lat, c.lng);
    onViewChange?.({
      zoom: map.getZoom(),
      south: b.getSouth(),
      west: b.getWest(),
      north: b.getNorth(),
      east: b.getEast(),
    });
  }, [map, onCenterChange, onViewChange]);

  useMapEvents({
    moveend() { report(); },
    zoomend() { report(); },
  });

  // Fire once on mount so the parent has a starting point.
  useEffect(() => {
    report();
  }, [report]);

  return null;
}
