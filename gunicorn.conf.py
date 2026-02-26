# Gunicorn configuration for deployment
# gunicorn auto-discovers this file by name.

import os

# Cloud Run sets PORT env var; Koyeb uses 8000
port = os.environ.get("PORT", "8000")
bind = f"0.0.0.0:{port}"

timeout = 300        # seconds – time-series can fetch many soundings
workers = int(os.environ.get("WEB_CONCURRENCY", "1"))
threads = int(os.environ.get("GUNICORN_THREADS", "4"))
graceful_timeout = 60
