// In dev mode (Vite), self-destruct so the SW doesn't intercept module requests
if (location.hostname === "localhost" || location.hostname === "127.0.0.1") {
  self.addEventListener("install", () => self.skipWaiting());
  self.addEventListener("activate", (event) => {
    event.waitUntil(
      caches.keys().then((keys) => Promise.all(keys.map((k) => caches.delete(k))))
        .then(() => self.registration.unregister())
        .then(() => self.clients.matchAll())
        .then((clients) => clients.forEach((c) => c.navigate(c.url)))
    );
  });
} else {

const CACHE_NAME = "sounding-v5";
const STATIC_ASSETS = [
  "./",
  "./index.html",
  "./favicon.svg",
  "./manifest.json",
];

// Install — cache static shell
self.addEventListener("install", (event) => {
  event.waitUntil(
    caches.open(CACHE_NAME).then((cache) => cache.addAll(STATIC_ASSETS))
  );
  self.skipWaiting();
});

// Activate — clean old caches
self.addEventListener("activate", (event) => {
  event.waitUntil(
    caches.keys().then((keys) =>
      Promise.all(keys.filter((k) => k !== CACHE_NAME).map((k) => caches.delete(k)))
    )
  );
  self.clients.claim();
});

// Fetch — network-first for API & navigation, cache-first for hashed assets
self.addEventListener("fetch", (event) => {
  const url = new URL(event.request.url);

  // Cross-origin requests (map tiles, fonts, CDN resources): let the browser
  // handle them directly — avoids CSP connect-src issues and opaque-response problems
  if (url.origin !== self.location.origin) return;

  // API calls: network-first, cache fallback for GET only
  if (url.pathname.startsWith("/api/")) {
    if (event.request.method !== "GET") return;
    event.respondWith(
      fetch(event.request)
        .then((res) => {
          if (res.ok) {
            const clone = res.clone();
            caches.open(CACHE_NAME).then((cache) => cache.put(event.request, clone));
          }
          return res;
        })
        .catch(() => caches.match(event.request))
    );
    return;
  }

  // Navigation requests (HTML pages): always network-first so deploys are seen immediately
  if (event.request.mode === "navigate") {
    event.respondWith(
      fetch(event.request)
        .then((res) => {
          const clone = res.clone();
          caches.open(CACHE_NAME).then((cache) => cache.put(event.request, clone));
          return res;
        })
        .catch(() => caches.match(event.request))
    );
    return;
  }

  // Hashed static assets (JS/CSS with content hash in filename): cache-first
  event.respondWith(
    caches.match(event.request).then((cached) => {
      if (cached) return cached;
      return fetch(event.request).then((res) => {
        if (res.ok) {
          const clone = res.clone();
          caches.open(CACHE_NAME).then((cache) => cache.put(event.request, clone));
        }
        return res;
      });
    })
  );
});

} // end else (production)
