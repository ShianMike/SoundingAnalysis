# Gunicorn configuration for Render deployment
# gunicorn auto-discovers this file by name.

timeout = 180        # seconds – risk scan can take a while
workers = 1          # keep memory usage low on free tier
threads = 4          # concurrency within the single worker
graceful_timeout = 60
