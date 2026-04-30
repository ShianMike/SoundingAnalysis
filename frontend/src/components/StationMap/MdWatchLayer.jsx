import { Pane, GeoJSON } from "react-leaflet";

const WATCH_STYLES = {
  "Tornado Watch":           { color: "#ff0000", fill: "#ff000022", label: "TOR" },
  "Severe Thunderstorm Watch": { color: "#ffa500", fill: "#ffa50022", label: "SVR" },
};

function mdStyle() {
  return {
    color: "#00e5ff",
    weight: 2,
    fillColor: "#00e5ff",
    fillOpacity: 0.08,
    dashArray: "6 4",
  };
}

function watchStyle(feature) {
  const ev = feature.properties?.type || feature.properties?.event || "";
  const s = Object.entries(WATCH_STYLES).find(([k]) => ev.includes(k));
  if (s) return { color: s[1].color, weight: 2.5, fillColor: s[1].color, fillOpacity: 0.1, dashArray: "10 5" };
  return { color: "#ffa500", weight: 2, fillColor: "#ffa500", fillOpacity: 0.08, dashArray: "10 5" };
}

/**
 * `MdWatchLayer` — overlay for SPC Mesoscale Discussions and Watch boxes.
 *
 *   - Watch polygons (Tornado / Severe Thunderstorm) get a colored,
 *     dashed boundary with a popup listing the watch number + expiry.
 *   - MD polygons get a cyan dashed boundary with a popup listing the
 *     MD number, the "concerning" hazard string, and the affected areas.
 *
 * Both data feeds come from the backend proxy and pass through
 * `props.mdData` / `props.watchData` as plain GeoJSON objects.
 */
export default function MdWatchLayer({ mdData, watchData }) {
  const hasMd = mdData?.features?.length > 0;
  const hasWatch = watchData?.features?.length > 0;
  if (!hasMd && !hasWatch) return null;
  return (
    <Pane name="spc-md-watch" style={{ zIndex: 415, pointerEvents: "auto" }}>
      {hasWatch && (
        <GeoJSON
          key={"w-" + watchData.features.length + (watchData.features[0]?.properties?.number || "")}
          data={watchData}
          style={watchStyle}
          onEachFeature={(feature, layer) => {
            const p = feature.properties;
            const ev = p.type || p.event || "Watch";
            const isTor = ev.includes("Tornado");
            const num = p.number || p.wn || "";
            const expires = p.expires || p.expiration || "";
            const expStr = expires ? new Date(expires).toLocaleString() : "";
            layer.bindPopup(
              `<div class="smap-warn-popup">
                <div class="smap-warn-badge" style="background:${isTor ? '#ff000022' : '#ffa50022'};color:${isTor ? '#ff0000' : '#ffa500'};border-color:${isTor ? '#ff000055' : '#ffa50055'}">
                  ${isTor ? 'TOR' : 'SVR'} WATCH ${num}
                </div>
                <div class="smap-warn-headline">${ev}${num ? ` #${num}` : ''}</div>
                ${expStr ? `<div class="smap-warn-time">Expires: ${expStr}</div>` : ""}
              </div>`,
              { className: "smap-popup smap-warn-popup-wrap", minWidth: 220, maxWidth: 320 }
            );
          }}
        />
      )}
      {hasMd && (
        <GeoJSON
          key={"md-" + mdData.features.length + (mdData.features[0]?.properties?.number || "")}
          data={mdData}
          style={mdStyle}
          onEachFeature={(feature, layer) => {
            const p = feature.properties;
            const num = p.number || p.md_num || "";
            const concerning = p.concerning || p.hazard || "";
            const areas = p.areas || p.states || "";
            layer.bindPopup(
              `<div class="smap-warn-popup">
                <div class="smap-warn-badge" style="background:#00e5ff22;color:#00e5ff;border-color:#00e5ff55">
                  MD ${num}
                </div>
                <div class="smap-warn-headline">Mesoscale Discussion ${num}</div>
                ${concerning ? `<div class="smap-warn-area">Concerning: ${concerning}</div>` : ""}
                ${areas ? `<div class="smap-warn-area">${areas}</div>` : ""}
              </div>`,
              { className: "smap-popup smap-warn-popup-wrap", minWidth: 220, maxWidth: 320 }
            );
          }}
        />
      )}
    </Pane>
  );
}
