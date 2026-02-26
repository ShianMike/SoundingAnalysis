/**
 * Sounding history stored in localStorage.
 * Each entry stores the request params, result metadata, and parameter values.
 * The base64 plot image is stored separately to allow listing without loading images.
 */

const HISTORY_KEY = "sounding_history";
const IMAGE_PREFIX = "sounding_img_";
const MAX_ENTRIES = 20;

/**
 * Get all history entries (without images).
 * Returns newest-first.
 */
export function getHistory() {
  try {
    const raw = localStorage.getItem(HISTORY_KEY);
    return raw ? JSON.parse(raw) : [];
  } catch {
    return [];
  }
}

/**
 * Save a sounding result to history.
 * @param {Object} requestParams - The params sent to the API (source, station, date, etc.)
 * @param {Object} result - The full API response { image, params, meta }
 */
export function saveToHistory(requestParams, result) {
  try {
    const entries = getHistory();
    const id = Date.now().toString(36) + Math.random().toString(36).slice(2, 6);

    const entry = {
      id,
      timestamp: Date.now(),
      requestParams,
      meta: result.meta,
      params: result.params,
    };

    // Store image separately (can be large)
    try {
      localStorage.setItem(IMAGE_PREFIX + id, result.image);
    } catch {
      // If storage is full, remove oldest images first
      _pruneOldImages(entries, 5);
      try {
        localStorage.setItem(IMAGE_PREFIX + id, result.image);
      } catch {
        // Still can't store — skip image
      }
    }

    // Add to front
    entries.unshift(entry);

    // Trim to max
    while (entries.length > MAX_ENTRIES) {
      const removed = entries.pop();
      localStorage.removeItem(IMAGE_PREFIX + removed.id);
    }

    localStorage.setItem(HISTORY_KEY, JSON.stringify(entries));
    return entry;
  } catch {
    return null;
  }
}

/**
 * Load a full history entry including the image.
 */
export function loadFromHistory(id) {
  const entries = getHistory();
  const entry = entries.find((e) => e.id === id);
  if (!entry) return null;

  const image = localStorage.getItem(IMAGE_PREFIX + id);
  return {
    image: image || null,
    params: entry.params,
    meta: entry.meta,
  };
}

/**
 * Delete a single history entry.
 */
export function deleteFromHistory(id) {
  try {
    const entries = getHistory().filter((e) => e.id !== id);
    localStorage.setItem(HISTORY_KEY, JSON.stringify(entries));
    localStorage.removeItem(IMAGE_PREFIX + id);
  } catch {
    // ignore
  }
}

/**
 * Clear all history.
 */
export function clearHistory() {
  try {
    const entries = getHistory();
    entries.forEach((e) => localStorage.removeItem(IMAGE_PREFIX + e.id));
    localStorage.removeItem(HISTORY_KEY);
  } catch {
    // ignore
  }
}

function _pruneOldImages(entries, count) {
  const toRemove = entries.slice(-count);
  toRemove.forEach((e) => localStorage.removeItem(IMAGE_PREFIX + e.id));
}
