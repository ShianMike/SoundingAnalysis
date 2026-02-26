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
