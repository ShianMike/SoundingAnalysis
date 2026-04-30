import { useEffect, useRef } from "react";
import { TileLayer } from "react-leaflet";

/**
 * `RadarLayer` — wrapper that keeps a SINGLE Leaflet `<TileLayer/>` alive
 * across URL changes by mutating it via `setUrl()` instead of unmounting +
 * recreating it. This avoids re-fetching every visible tile on every
 * animation frame, which previously caused jank and burned through tile
 * provider rate limits.
 *
 * Drop-in replacement for `<TileLayer/>` — accepts the same props.
 */
export default function RadarLayer({ url, opacity, ...rest }) {
  const ref = useRef(null);
  const prevUrl = useRef(url);

  useEffect(() => {
    if (ref.current) ref.current.setOpacity(opacity);
  }, [opacity]);

  useEffect(() => {
    if (ref.current && url !== prevUrl.current) {
      ref.current.setUrl(url, false);
      ref.current.redraw();
      prevUrl.current = url;
    }
  }, [url]);

  return <TileLayer ref={ref} url={url} opacity={opacity} {...rest} />;
}
