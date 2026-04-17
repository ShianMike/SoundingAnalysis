"""
Flask API for the Sounding Analysis Tool.
Slim entry point — all route logic lives in the routes/ package.
"""

import os
import re

from flask import Flask, send_from_directory, request, abort
from flask_cors import CORS
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from flask_talisman import Talisman

import routes  # noqa: registers blueprints

FRONTEND_DIR = os.path.join(os.path.dirname(__file__), "frontend", "dist")

# ─── Allowed origins (lock down CORS) ──────────────────────────────
ALLOWED_ORIGINS = [
    "https://soundingscopepy.app",
    "https://www.soundingscopepy.app",
    "https://shianmike.github.io",
    "https://soundinganalysis-752306366750.asia-southeast1.run.app",
    "https://soundinganalysis-uvktu4ziyq-as.a.run.app",
]
# In local dev, also allow localhost origins
if os.environ.get("FLASK_DEBUG") or os.environ.get("FLASK_ENV") == "development":
    ALLOWED_ORIGINS += ["http://localhost:3000", "http://localhost:5000", "http://127.0.0.1:3000", "http://127.0.0.1:5000"]

app = Flask(__name__, static_folder=FRONTEND_DIR, static_url_path="")

# ─── CORS — restricted to known origins ────────────────────────────
CORS(app, resources={r"/api/*": {"origins": ALLOWED_ORIGINS}}, supports_credentials=False)

# ─── Rate limiting ─────────────────────────────────────────────────
limiter = Limiter(
    get_remote_address,
    app=app,
    default_limits=["200 per minute", "30 per second"],
    storage_uri="memory://",
)

# ─── Security headers via Talisman ─────────────────────────────────
# CSP must allow external tiles, fonts, and API calls the app uses.
_is_production = bool(os.environ.get("K_SERVICE"))  # Cloud Run sets K_SERVICE

csp = {
    "default-src": "'self'",
    "script-src":  "'self' 'unsafe-inline'",
    "style-src":   "'self' 'unsafe-inline' https://fonts.googleapis.com https://unpkg.com",
    "font-src":    "'self' https://fonts.gstatic.com data:",
    "img-src":     "'self' data: blob: https://*.basemaps.cartocdn.com https://server.arcgisonline.com "
                   "https://*.tile.opentopomap.org https://tilecache.rainviewer.com "
                   "https://mesonet.agron.iastate.edu https://*.openstreetmap.org",
    "connect-src": "'self' https://api.rainviewer.com https://mesonet.agron.iastate.edu "
                   "https://api.open-meteo.com "
                   "https://www.spc.noaa.gov https://spc.noaa.gov "
                   "https://api.weather.gov "
                   "https://api.livestormchasing.com https://edge.livestormchasing.com "
                   "https://soundingscopepy.app https://www.soundingscopepy.app "
                   "https://soundinganalysis-752306366750.asia-southeast1.run.app "
                   "https://soundinganalysis-uvktu4ziyq-as.a.run.app "
                   "https://*.run.app https://*.a.run.app",
    "media-src":   "'self' blob:",
    "frame-ancestors": "'none'",
    "base-uri":    "'self'",
    "form-action": "'self'",
    "object-src":  "'none'",
}

Talisman(
    app,
    force_https=_is_production,           # Only force HTTPS in prod (Cloud Run)
    force_https_permanent=False,
    strict_transport_security=True,
    strict_transport_security_max_age=63072000,  # 2 years (HSTS preload ready)
    strict_transport_security_include_subdomains=True,
    strict_transport_security_preload=True,
    content_security_policy=csp,
    content_security_policy_nonce_in=["script-src"],
    referrer_policy="strict-origin-when-cross-origin",
    frame_options="DENY",
    feature_policy={
        "geolocation":  "'none'",
        "camera":       "'none'",
        "microphone":   "'none'",
        "payment":      "'none'",
    },
    permissions_policy={
        "geolocation":       "()",
        "camera":            "()",
        "microphone":        "()",
        "payment":           "()",
        "usb":               "()",
        "magnetometer":      "()",
        "gyroscope":         "()",
        "accelerometer":     "()",
        "autoplay":          "()",
        "display-capture":   "()",
        "document-domain":   "()",
        "fullscreen":        "(self)",
        "interest-cohort":   "()",
    },
    session_cookie_secure=_is_production,
    session_cookie_http_only=True,
    session_cookie_samesite="Lax",
)


@app.after_request
def add_extra_security_headers(response):
    """Extra headers beyond what Talisman sets."""
    response.headers["X-Content-Type-Options"] = "nosniff"
    response.headers["X-Frame-Options"] = "DENY"
    response.headers["X-Permitted-Cross-Domain-Policies"] = "none"
    response.headers["Cross-Origin-Opener-Policy"] = "same-origin"
    response.headers["Cross-Origin-Resource-Policy"] = "same-origin"
    # Cache control for API responses
    if request.path.startswith("/api/"):
        response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate, max-age=0"
        response.headers["Pragma"] = "no-cache"
    return response


# ─── Request-size limit (16 MB — protects upload endpoints) ────────
app.config["MAX_CONTENT_LENGTH"] = 16 * 1024 * 1024


# ─── Input validation: reject suspicious path traversal ────────────
@app.before_request
def block_path_traversal():
    if ".." in request.path or re.search(r"[<>\"';\x00]", request.path):
        abort(400)


# Register all API blueprints
routes.register_all(app)


# ── Apply stricter rate limits to sensitive endpoints ──────────────
limiter.limit("10 per minute")(app.view_functions.get("feedback.submit_feedback", lambda: None))
limiter.limit("30 per minute")(app.view_functions.get("sounding_routes.get_sounding", lambda: None))


# ─── SPA catch-all ──────────────────────────────────────────────────
@app.route("/SoundingAnalysis/", defaults={"path": ""})
@app.route("/SoundingAnalysis/<path:path>")
def serve_spa_prefixed(path):
    """Serve the React SPA under the /SoundingAnalysis/ base path."""
    file_path = os.path.join(FRONTEND_DIR, path)
    if path and os.path.isfile(file_path):
        return send_from_directory(FRONTEND_DIR, path)
    return send_from_directory(FRONTEND_DIR, "index.html")


@app.route("/", defaults={"path": ""})
@app.route("/<path:path>")
def serve_spa(path):
    """Serve the React SPA. API routes are matched first by Flask."""
    file_path = os.path.join(FRONTEND_DIR, path)
    if path and os.path.isfile(file_path):
        return send_from_directory(FRONTEND_DIR, path)
    return send_from_directory(FRONTEND_DIR, "index.html")


if __name__ == "__main__":
    print("Starting Sounding API on http://localhost:5000")
    app.run(debug=True, port=5000)
