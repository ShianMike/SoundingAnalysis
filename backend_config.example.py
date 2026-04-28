# backend_config.example.py — template for backend_config.py (safe to commit)
# Copy this file to backend_config.py and fill in your own values.

# ── Allowed CORS origins ───────────────────────────────────────────────────
ALLOWED_ORIGINS = [
    "https://your-app.example.com",
]

# ── GCP / Cloud Run ────────────────────────────────────────────────────────
GCP_PROJECT  = "your-gcp-project-id"
GCP_SERVICE  = "your-cloud-run-service"
GCP_REGION   = "us-central1"
