import { useMapEvents } from "react-leaflet";

/**
 * `MapClickHandler` — invisible interaction helper that turns Leaflet
 * pointer events into the parent's domain callbacks.
 *
 * Props:
 *   - `enabled` — when true, regular left-clicks trigger `onLatLonSelect`.
 *   - `onLatLonSelect(lat, lon)` — receives 2-decimal-place strings.
 *   - `onMapClick(lat, lng, event)` — generic click hook. If it returns
 *     `true`, the standard `onLatLonSelect` path is suppressed.
 *   - `onRightClick(lat, lng)` — context-menu (right-click) coords.
 *   - `onMouseMove(lat, lng)` — hover coords (cheap, fires often).
 *   - `disableRightClick` — when true, swallow context menu without
 *     dispatching `onRightClick`.
 *
 * Renders nothing.
 */
export default function MapClickHandler({
  enabled,
  onLatLonSelect,
  onRightClick,
  onMapClick,
  onMouseMove,
  disableRightClick = false,
}) {
  useMapEvents({
    click(e) {
      if (onMapClick && onMapClick(e.latlng.lat, e.latlng.lng, e) === true) return;
      if (enabled && onLatLonSelect) {
        onLatLonSelect(e.latlng.lat.toFixed(2), e.latlng.lng.toFixed(2));
      }
    },
    mousemove(e) {
      onMouseMove?.(e.latlng.lat, e.latlng.lng);
    },
    contextmenu(e) {
      if (disableRightClick) {
        e.originalEvent.preventDefault();
        return;
      }
      // Right-click → proximity search
      if (onRightClick) {
        e.originalEvent.preventDefault();
        onRightClick(e.latlng.lat, e.latlng.lng);
      }
    },
  });
  return null;
}
