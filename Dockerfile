FROM python:3.11-slim

WORKDIR /app

# Install system deps for matplotlib
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc g++ pkg-config libfreetype6-dev libpng-dev \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# ── Security: run as non-root user ──────────────────────────────
RUN addgroup --system appgroup && adduser --system --ingroup appgroup appuser \
    && chown -R appuser:appgroup /app
USER appuser

EXPOSE 8000

# gunicorn.conf.py is auto-discovered and handles PORT env var
CMD ["gunicorn", "app:app"]
