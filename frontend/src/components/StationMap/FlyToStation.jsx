import { useEffect } from "react";
import { useMap } from "react-leaflet";

/**
 * `FlyToStation` — invisible helper that pans + zooms the Leaflet map to
 * the given station whenever the station id changes.
 *
 * Used by the parent to keep the map view synchronized with the user's
 * sidebar station selection.
 */
export default function FlyToStation({ station, stations }) {
  const map = useMap();
  useEffect(() => {
    if (!station) return;
    const s = stations.find((st) => st.id === station);
    if (s) map.flyTo([s.lat, s.lon], 7, { duration: 0.6 });
  }, [station]); // eslint-disable-line react-hooks/exhaustive-deps
  return null;
}
