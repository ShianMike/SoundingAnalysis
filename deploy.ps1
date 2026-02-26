#!/usr/bin/env pwsh
# deploy.ps1 - Build frontend and deploy to GitHub Pages
# Usage: .\deploy.ps1

$ErrorActionPreference = "Stop"

$API_URL = "https://soundinganalysis-uvktu4ziyq-as.a.run.app"
$REPO    = "https://github.com/ShianMike/SoundingAnalysis.git"
$BASE    = "/SoundingAnalysis/"

Write-Host "Building frontend..." -ForegroundColor Cyan
Push-Location frontend
$env:VITE_API_URL = $API_URL
npx vite build --base=$BASE
Pop-Location

$tmpDir = "$env:TEMP\gh-pages-deploy"
Remove-Item -Recurse -Force $tmpDir -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Path $tmpDir -Force | Out-Null
Copy-Item "frontend\dist\*" $tmpDir -Recurse

Push-Location $tmpDir
git init
git checkout -b gh-pages
git add -A
git commit -m "Deploy $(Get-Date -Format 'yyyy-MM-dd HH:mm')"
git remote add origin $REPO
$null = git push origin gh-pages --force 2>&1
Pop-Location

Remove-Item -Recurse -Force $tmpDir -ErrorAction SilentlyContinue
Write-Host "`nDeployed to https://shianmike.github.io/SoundingAnalysis/" -ForegroundColor Green
