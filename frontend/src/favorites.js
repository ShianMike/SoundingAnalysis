const FAVORITES_KEY = "sounding_favorites";

export function getFavorites() {
  try {
    return JSON.parse(localStorage.getItem(FAVORITES_KEY)) || [];
  } catch {
    return [];
  }
}

export function addFavorite(stationId) {
  const favs = getFavorites();
  if (!favs.includes(stationId)) {
    favs.unshift(stationId);
    localStorage.setItem(FAVORITES_KEY, JSON.stringify(favs));
  }
  return favs;
}

export function removeFavorite(stationId) {
  const favs = getFavorites().filter((id) => id !== stationId);
  localStorage.setItem(FAVORITES_KEY, JSON.stringify(favs));
  return favs;
}

export function isFavorite(stationId) {
  return getFavorites().includes(stationId);
}

export function toggleFavorite(stationId) {
  if (isFavorite(stationId)) {
    return { favorites: removeFavorite(stationId), added: false };
  }
  return { favorites: addFavorite(stationId), added: true };
}

/* ── Station Groups ─────────────────────────────────────── */
const GROUPS_KEY = "sounding_station_groups";

export function getStationGroups() {
  try {
    return JSON.parse(localStorage.getItem(GROUPS_KEY)) || [];
  } catch {
    return [];
  }
}

export function saveStationGroup(name, stationIds) {
  const groups = getStationGroups();
  const idx = groups.findIndex((g) => g.name === name);
  if (idx >= 0) {
    groups[idx].stations = stationIds;
  } else {
    groups.push({ name, stations: stationIds });
  }
  localStorage.setItem(GROUPS_KEY, JSON.stringify(groups));
  return groups;
}

export function deleteStationGroup(name) {
  const groups = getStationGroups().filter((g) => g.name !== name);
  localStorage.setItem(GROUPS_KEY, JSON.stringify(groups));
  return groups;
}
