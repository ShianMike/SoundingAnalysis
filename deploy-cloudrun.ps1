#!/usr/bin/env pwsh
# deploy-cloudrun.ps1 - Deploy backend to Google Cloud Run
# Prerequisites:
#   1. Install Google Cloud SDK: https://cloud.google.com/sdk/docs/install
#   2. Run: gcloud auth login
#   3. Run: gcloud projects create sounding-analysis --name="Sounding Analysis"
#      (or use an existing project)
#   4. Run: gcloud config set project sounding-analysis
#   5. Enable billing on the project (required even for free tier)
#
# Usage: .\deploy-cloudrun.ps1

$ErrorActionPreference = "Stop"

# ── Configuration ──────────────────────────────────────────────
$PROJECT    = "sounding-analysis"      # your GCP project ID
$SERVICE    = "soundinganalysis"       # Cloud Run service name
$REGION     = "asia-southeast1"        # Singapore (closest to you)
$MEMORY     = "512Mi"                  # 512 MB RAM
$CPU        = "1"                      # 1 vCPU
$TIMEOUT    = "300"                    # 5 min request timeout
$MAX_INSTANCES = "2"                   # cap to stay in free tier
$CONCURRENCY   = "10"                  # concurrent requests per instance

# ── Ensure APIs are enabled ────────────────────────────────────
Write-Host "Enabling Cloud Run & Artifact Registry APIs..." -ForegroundColor Cyan
gcloud services enable run.googleapis.com artifactregistry.googleapis.com --project $PROJECT

# ── Deploy directly from source (Cloud Build) ─────────────────
Write-Host "`nDeploying to Cloud Run ($REGION)..." -ForegroundColor Cyan
gcloud run deploy $SERVICE `
    --source . `
    --project $PROJECT `
    --region $REGION `
    --platform managed `
    --allow-unauthenticated `
    --memory $MEMORY `
    --cpu $CPU `
    --timeout $TIMEOUT `
    --max-instances $MAX_INSTANCES `
    --concurrency $CONCURRENCY `
    --set-env-vars "WEB_CONCURRENCY=2,GUNICORN_THREADS=4,TS_WORKERS=4,TS_WALL_LIMIT=280,SCAN_WORKERS=6" `
    --port 8080

# ── Print the service URL ─────────────────────────────────────
$url = gcloud run services describe $SERVICE --region $REGION --project $PROJECT --format "value(status.url)"
Write-Host "`nDeployed! Backend URL:" -ForegroundColor Green
Write-Host $url -ForegroundColor Yellow
Write-Host "`nUpdate deploy.ps1 line 8 with this URL, then run .\deploy.ps1 to redeploy frontend." -ForegroundColor Cyan
