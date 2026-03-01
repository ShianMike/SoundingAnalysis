/**
 * Sounding history — metadata in localStorage, images in IndexedDB.
 * IndexedDB has ~50 MB+ quota vs localStorage's ~5-10 MB,
 * so large base64 sounding PNGs survive across reloads reliably.
 */

const HISTORY_KEY = "sounding_history";
const IMAGE_PREFIX = "sounding_img_";
const MAX_ENTRIES = 20;

/* ── IndexedDB helpers ────────────────────────────────────────── */

const DB_NAME = "SoundingHistoryDB";
const DB_VERSION = 1;
const IMG_STORE = "images";

function _openDB() {
  return new Promise((resolve, reject) => {
    const req = indexedDB.open(DB_NAME, DB_VERSION);
    req.onupgradeneeded = (e) => {
      const db = e.target.result;
      if (!db.objectStoreNames.contains(IMG_STORE)) {
        db.createObjectStore(IMG_STORE);
      }
    };
    req.onsuccess = () => resolve(req.result);
    req.onerror = () => reject(req.error);
  });
}

async function _putImage(key, data) {
  try {
    const db = await _openDB();
    return new Promise((resolve, reject) => {
      const tx = db.transaction(IMG_STORE, "readwrite");
      tx.objectStore(IMG_STORE).put(data, key);
      tx.oncomplete = () => { db.close(); resolve(); };
      tx.onerror = () => { db.close(); reject(tx.error); };
    });
  } catch {
    // IndexedDB unavailable — silent fail
  }
}

async function _getImage(key) {
  try {
    const db = await _openDB();
    return new Promise((resolve, reject) => {
      const tx = db.transaction(IMG_STORE, "readonly");
      const req = tx.objectStore(IMG_STORE).get(key);
      req.onsuccess = () => { db.close(); resolve(req.result || null); };
      req.onerror = () => { db.close(); reject(req.error); };
    });
  } catch {
    return null;
  }
}

async function _deleteImage(key) {
  try {
    const db = await _openDB();
    return new Promise((resolve, reject) => {
      const tx = db.transaction(IMG_STORE, "readwrite");
      tx.objectStore(IMG_STORE).delete(key);
      tx.oncomplete = () => { db.close(); resolve(); };
      tx.onerror = () => { db.close(); reject(tx.error); };
    });
  } catch {
    // ignore
  }
}

/* ── Public API ───────────────────────────────────────────────── */

/**
 * Get all history entries (without images).  Returns newest-first.
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
 * Metadata → localStorage, image → IndexedDB.
 */
export async function saveToHistory(requestParams, result) {
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

    // Store image in IndexedDB (much larger quota than localStorage)
    await _putImage(IMAGE_PREFIX + id, result.image);

    // Add to front
    entries.unshift(entry);

    // Trim to max
    while (entries.length > MAX_ENTRIES) {
      const removed = entries.pop();
      _deleteImage(IMAGE_PREFIX + removed.id);          // async cleanup
      localStorage.removeItem(IMAGE_PREFIX + removed.id); // legacy cleanup
    }

    localStorage.setItem(HISTORY_KEY, JSON.stringify(entries));
    return entry;
  } catch {
    return null;
  }
}

/**
 * Load a full history entry including the image (from IndexedDB).
 * Falls back to localStorage for entries saved before the IndexedDB migration.
 */
export async function loadFromHistory(id) {
  const entries = getHistory();
  const entry = entries.find((e) => e.id === id);
  if (!entry) return null;

  // Try IndexedDB first
  let image = await _getImage(IMAGE_PREFIX + id);

  // Fall back to localStorage (legacy migration)
  if (!image) {
    image = localStorage.getItem(IMAGE_PREFIX + id);
    if (image) {
      // Migrate to IndexedDB, then remove from localStorage
      _putImage(IMAGE_PREFIX + id, image).then(() => {
        localStorage.removeItem(IMAGE_PREFIX + id);
      });
    }
  }

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
    _deleteImage(IMAGE_PREFIX + id);          // async fire-and-forget
    localStorage.removeItem(IMAGE_PREFIX + id); // legacy cleanup
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
    entries.forEach((e) => {
      _deleteImage(IMAGE_PREFIX + e.id);          // async fire-and-forget
      localStorage.removeItem(IMAGE_PREFIX + e.id); // legacy cleanup
    });
    localStorage.removeItem(HISTORY_KEY);
  } catch {
    // ignore
  }
}

/* ── Comparison History ─────────────────────────────────────────── */

const COMPARE_KEY = "comparison_history";
const COMPARE_IMG_PREFIX = "compare_img_";
const MAX_COMPARE = 15;

/**
 * Get all comparison history entries (without images). Newest-first.
 */
export function getCompareHistory() {
  try {
    const raw = localStorage.getItem(COMPARE_KEY);
    return raw ? JSON.parse(raw) : [];
  } catch {
    return [];
  }
}

/**
 * Save a comparison to history.
 * Metadata → localStorage, images → IndexedDB.
 */
export async function saveCompareToHistory(slots, results) {
  try {
    const entries = getCompareHistory();
    const id = Date.now().toString(36) + Math.random().toString(36).slice(2, 6);

    const summary = results.map((r) => ({
      meta: r.meta,
      params: r.params,
      error: r.error || null,
    }));

    const entry = {
      id,
      timestamp: Date.now(),
      slots: slots.filter((s) => s.station),
      summary,
    };

    // Store images in IndexedDB as JSON array
    const images = results.map((r) => r.image || null);
    await _putImage(COMPARE_IMG_PREFIX + id, JSON.stringify(images));

    entries.unshift(entry);
    while (entries.length > MAX_COMPARE) {
      const removed = entries.pop();
      _deleteImage(COMPARE_IMG_PREFIX + removed.id);
      localStorage.removeItem(COMPARE_IMG_PREFIX + removed.id);
    }

    localStorage.setItem(COMPARE_KEY, JSON.stringify(entries));
    return entry;
  } catch {
    return null;
  }
}

/**
 * Load a comparison entry including images (from IndexedDB).
 * Falls back to localStorage for legacy entries.
 */
export async function loadCompareFromHistory(id) {
  const entries = getCompareHistory();
  const entry = entries.find((e) => e.id === id);
  if (!entry) return null;

  // Try IndexedDB first
  let imagesRaw = await _getImage(COMPARE_IMG_PREFIX + id);

  // Fall back to localStorage (legacy)
  if (!imagesRaw) {
    imagesRaw = localStorage.getItem(COMPARE_IMG_PREFIX + id);
    if (imagesRaw) {
      _putImage(COMPARE_IMG_PREFIX + id, imagesRaw).then(() => {
        localStorage.removeItem(COMPARE_IMG_PREFIX + id);
      });
    }
  }

  let images = [];
  try {
    images = JSON.parse(imagesRaw || "[]");
  } catch {
    images = [];
  }

  // Reconstruct results array
  const results = entry.summary.map((s, i) => ({
    ...s,
    image: images[i] || null,
  }));

  return { slots: entry.slots, results };
}

/**
 * Delete a single comparison history entry.
 */
export function deleteCompareFromHistory(id) {
  try {
    const entries = getCompareHistory().filter((e) => e.id !== id);
    localStorage.setItem(COMPARE_KEY, JSON.stringify(entries));
    _deleteImage(COMPARE_IMG_PREFIX + id);
    localStorage.removeItem(COMPARE_IMG_PREFIX + id);
  } catch {
    // ignore
  }
}

/**
 * Clear all comparison history.
 */
export function clearCompareHistory() {
  try {
    const entries = getCompareHistory();
    entries.forEach((e) => {
      _deleteImage(COMPARE_IMG_PREFIX + e.id);
      localStorage.removeItem(COMPARE_IMG_PREFIX + e.id);
    });
    localStorage.removeItem(COMPARE_KEY);
  } catch {
    // ignore
  }
}
