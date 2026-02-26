# Gunicorn configuration for deployment
# gunicorn auto-discovers this file by name.

bind = "0.0.0.0:8000"
timeout = 300        # seconds – time-series can fetch many soundings
workers = 1          # keep memory usage low on free tier
threads = 4          # concurrency within the single worker
graceful_timeout = 60
