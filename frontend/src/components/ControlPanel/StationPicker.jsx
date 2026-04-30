import { Search, ArrowUpDown, Star } from "lucide-react";

/**
 * `StationPicker` — search bar, sort selector, and clickable station list.
 *
 * Pure presentational component extracted from `ControlPanel.jsx`. State
 * (search query, sort mode, favorites set, current selection, risk overlay)
 * is owned by the parent and passed in.
 *
 * Props:
 *   - filteredStations: pre-sorted/filtered list of stations to display.
 *       Each item may carry an optional `risk` field with shape
 *       `{ rank, stp, raw, cape, srh, bwd }`.
 *   - selectedStation: currently selected station id (string).
 *   - onStationSelect(id): called when a station row is clicked.
 *   - stationSearch / onStationSearchChange: filter input value + setter.
 *   - sortMode / onSortModeChange: current sort key + setter.
 *   - favorites: array of favorited station ids.
 *   - onToggleFavorite(id): called when the star icon is clicked.
 *   - riskData: pass-through (only used to color rows). May be null.
 *   - listRef: parent-owned ref attached to the list container so the parent
 *       can scroll the active row into view.
 */
export default function StationPicker({
  filteredStations,
  selectedStation,
  onStationSelect,
  stationSearch,
  onStationSearchChange,
  sortMode,
  onSortModeChange,
  favorites,
  onToggleFavorite,
  riskData,
  listRef,
}) {
  return (
    <>
      <div className="cp-station-toolbar">
        <div className="cp-input-wrap cp-search-flex">
          <Search size={14} className="cp-input-icon" />
          <input
            type="text"
            className="cp-input"
            placeholder="Filter..."
            value={stationSearch}
            onChange={(e) => onStationSearchChange(e.target.value)}
          />
        </div>
        <div className="cp-sort-wrap">
          <ArrowUpDown size={12} className="cp-sort-icon" />
          <select
            className="cp-sort-select"
            value={sortMode}
            onChange={(e) => onSortModeChange(e.target.value)}
          >
            <option value="az">A → Z</option>
            <option value="za">Z → A</option>
            <option value="favs">★ Favs</option>
            <option value="risk-high">Risk ↓</option>
            <option value="risk-low">Risk ↑</option>
          </select>
        </div>
      </div>

      <div className="cp-station-list" ref={listRef}>
        {filteredStations.length === 0 && (
          <div className="cp-station-empty">No stations found</div>
        )}
        {filteredStations.map((s) => (
          <button
            key={s.id}
            type="button"
            data-id={s.id}
            className={`cp-station-item ${selectedStation === s.id ? "active" : ""}`}
            onClick={() => onStationSelect(s.id, !!riskData)}
          >
            <span
              className={`cp-fav-star ${favorites.includes(s.id) ? "faved" : ""}`}
              onClick={(e) => {
                e.stopPropagation();
                onToggleFavorite(s.id);
              }}
              title={favorites.includes(s.id) ? "Remove from favorites" : "Add to favorites"}
            >
              <Star size={12} fill={favorites.includes(s.id) ? "currentColor" : "none"} />
            </span>
            <span className="cp-station-item-id">{s.id}</span>
            <span className="cp-station-item-name">{s.name}</span>
            {s.risk ? (
              <span className={`cp-risk-score ${s.risk.stp >= 1 ? "high" : s.risk.stp >= 0.3 ? "med" : "low"}`}>
                {s.risk.stp.toFixed(1)}
              </span>
            ) : (
              <span className="cp-station-item-coords">
                {s.lat.toFixed(1)}, {s.lon.toFixed(1)}
              </span>
            )}
          </button>
        ))}
      </div>
    </>
  );
}
