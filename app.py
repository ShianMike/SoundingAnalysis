"""
Flask API for the Sounding Analysis Tool.
Slim entry point — all route logic lives in the routes/ package.
"""

import os

from flask import Flask, send_from_directory
from flask_cors import CORS

import routes  # noqa: registers blueprints

FRONTEND_DIR = os.path.join(os.path.dirname(__file__), "frontend", "dist")

app = Flask(__name__, static_folder=FRONTEND_DIR, static_url_path="")
CORS(app, resources={r"/api/*": {"origins": "*"}}, supports_credentials=False)


@app.after_request
def add_cors_headers(response):
    response.headers["Access-Control-Allow-Origin"] = "*"
    response.headers["Access-Control-Allow-Methods"] = "GET, POST, OPTIONS"
    response.headers["Access-Control-Allow-Headers"] = "Content-Type"
    return response


# Register all API blueprints
routes.register_all(app)


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
