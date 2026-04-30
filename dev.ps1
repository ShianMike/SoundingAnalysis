<#
.SYNOPSIS
    Start backend (Flask) and frontend (Vite) for local development.
.DESCRIPTION
    - Activates the Python venv and runs Flask on port 5001  (Vite proxies /api → 5001)
    - Installs npm deps (if needed) and runs Vite dev server on port 3000
    - Watches Python/route/template files → auto-restarts Flask on changes
    - Vite handles frontend HMR natively
    - Press Ctrl+C to stop everything.
#>

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path

# ── Colors ──────────────────────────────────────────────────────────
function Write-Header($msg) { Write-Host "`n>> $msg" -ForegroundColor Cyan }

# ── 1. Python venv ─────────────────────────────────────────────────
$VenvDir = Join-Path $Root ".venv"
if (-not (Test-Path $VenvDir)) {
    Write-Header "Creating Python virtual environment..."
    python -m venv $VenvDir
}

$VenvPython = Join-Path $VenvDir "Scripts\python.exe"

# ── 2. Install Python dependencies ─────────────────────────────────
Write-Header "Installing Python dependencies..."
& $VenvPython -m pip install -q -r (Join-Path $Root "requirements.txt")

# ── 3. Install Node dependencies ───────────────────────────────────
$FrontendDir = Join-Path $Root "frontend"
$NodeModules = Join-Path $FrontendDir "node_modules"
if (-not (Test-Path $NodeModules)) {
    Write-Header "Installing Node dependencies..."
    Push-Location $FrontendDir
    npm install
    Pop-Location
} else {
    Write-Header "Node modules already installed. Skipping npm install."
}

# ── Helper: start / restart Flask backend job ───────────────────────
$script:backendJob = $null

function Start-Backend {
    if ($script:backendJob) {
        Stop-Job  $script:backendJob -ErrorAction SilentlyContinue
        Remove-Job $script:backendJob -Force -ErrorAction SilentlyContinue
    }
    Write-Host "[watcher] Starting Flask backend..." -ForegroundColor Yellow
    $script:backendJob = Start-Job -ScriptBlock {
        param($root, $venvPython)
        Set-Location $root
        $env:FLASK_DEBUG = "1"
        & $venvPython -m flask --app app run --host 127.0.0.1 --port 5001 --reload
    } -ArgumentList $Root, $VenvPython
}

# ── 4. Start Flask backend ─────────────────────────────────────────
Write-Header "Starting Flask backend on http://localhost:5001 ..."
Start-Backend

# ── 5. Start Vite frontend (background job) ────────────────────────
Write-Header "Starting Vite dev server on http://localhost:3000 ..."
$frontendJob = Start-Job -ScriptBlock {
    param($dir)
    Set-Location $dir
    npx vite --host 127.0.0.1 --port 3000
} -ArgumentList $FrontendDir

# ── 6. File watcher — restart backend on Python file changes ───────
#    Flask --reload handles most .py changes, but this also catches
#    new files, requirements.txt, and template changes that the
#    reloader can miss.
$watchers = @()

# Watch Python files in project root (app.py, etc.)
$pyWatcher = [System.IO.FileSystemWatcher]::new($Root, "*.py")
$pyWatcher.IncludeSubdirectories = $false
$pyWatcher.NotifyFilter = [System.IO.NotifyFilters]::LastWrite -bor [System.IO.NotifyFilters]::FileName
$watchers += $pyWatcher

# Watch routes/ folder
$routesDir = Join-Path $Root "routes"
if (Test-Path $routesDir) {
    $routeWatcher = [System.IO.FileSystemWatcher]::new($routesDir, "*.py")
    $routeWatcher.IncludeSubdirectories = $true
    $routeWatcher.NotifyFilter = [System.IO.NotifyFilters]::LastWrite -bor [System.IO.NotifyFilters]::FileName
    $watchers += $routeWatcher
}

# Watch sounding/ package
$soundingDir = Join-Path $Root "sounding"
if (Test-Path $soundingDir) {
    $soundingWatcher = [System.IO.FileSystemWatcher]::new($soundingDir, "*.py")
    $soundingWatcher.IncludeSubdirectories = $true
    $soundingWatcher.NotifyFilter = [System.IO.NotifyFilters]::LastWrite -bor [System.IO.NotifyFilters]::FileName
    $watchers += $soundingWatcher
}

# Watch forecast/ package
$forecastDir = Join-Path $Root "forecast"
if (Test-Path $forecastDir) {
    $fcWatcher = [System.IO.FileSystemWatcher]::new($forecastDir, "*.py")
    $fcWatcher.IncludeSubdirectories = $true
    $fcWatcher.NotifyFilter = [System.IO.NotifyFilters]::LastWrite -bor [System.IO.NotifyFilters]::FileName
    $watchers += $fcWatcher
}

# Watch requirements.txt
$reqWatcher = [System.IO.FileSystemWatcher]::new($Root, "requirements.txt")
$reqWatcher.NotifyFilter = [System.IO.NotifyFilters]::LastWrite
$watchers += $reqWatcher

# Enable all watchers
foreach ($w in $watchers) { $w.EnableRaisingEvents = $true }

# Debounce tracker — avoid rapid-fire restarts
$script:lastRestart = [datetime]::MinValue

# ── 7. Main loop — stream output + poll file watchers ──────────────
Write-Host ""
Write-Host "============================================" -ForegroundColor Green
Write-Host "  Backend:  http://localhost:5001"            -ForegroundColor Green
Write-Host "  Frontend: http://localhost:3000"            -ForegroundColor Green
Write-Host "  Auto-restart: watching .py + requirements"  -ForegroundColor Green
Write-Host "  Press Ctrl+C to stop everything."           -ForegroundColor Yellow
Write-Host "============================================" -ForegroundColor Green
Write-Host ""

try {
    while ($true) {
        # Stream output from both jobs
        Receive-Job $script:backendJob -ErrorAction SilentlyContinue |
            ForEach-Object { Write-Host "[backend]  $_" -ForegroundColor Magenta }
        Receive-Job $frontendJob -ErrorAction SilentlyContinue |
            ForEach-Object { Write-Host "[frontend] $_" -ForegroundColor Blue }

        # Poll each watcher for changes (non-blocking, 0ms timeout)
        $changed = $false
        $changedFile = ""
        foreach ($w in $watchers) {
            $evt = $w.WaitForChanged([System.IO.WatcherChangeTypes]::All, 0)
            if (-not $evt.TimedOut) {
                $changed = $true
                $changedFile = $evt.Name
                break
            }
        }

        if ($changed) {
            $now = [datetime]::Now
            # Debounce: ignore if less than 2 seconds since last restart
            if (($now - $script:lastRestart).TotalSeconds -ge 2) {
                $script:lastRestart = $now
                Write-Host ""
                Write-Host "[watcher] Change detected: $changedFile" -ForegroundColor Yellow

                # If requirements.txt changed, reinstall deps first
                if ($changedFile -eq "requirements.txt") {
                    Write-Host "[watcher] Reinstalling Python dependencies..." -ForegroundColor Yellow
                    & $VenvPython -m pip install -q -r (Join-Path $Root "requirements.txt")
                }

                Write-Host "[watcher] Restarting backend..." -ForegroundColor Yellow
                Start-Backend
            }
        }

        # Check if backend crashed (not just reloader restart)
        if ($script:backendJob.State -eq "Failed") {
            Write-Host "[backend]  Job failed — restarting in 3s..." -ForegroundColor Red
            Receive-Job $script:backendJob -ErrorAction SilentlyContinue |
                ForEach-Object { Write-Host "[backend]  $_" -ForegroundColor Red }
            Start-Sleep -Seconds 3
            Start-Backend
        }

        if ($frontendJob.State -eq "Failed") {
            Write-Host "[frontend] Job failed!" -ForegroundColor Red
            Receive-Job $frontendJob -ErrorAction SilentlyContinue |
                ForEach-Object { Write-Host "[frontend] $_" -ForegroundColor Red }
            break
        }

        Start-Sleep -Milliseconds 500
    }
} finally {
    Write-Header "Shutting down..."
    # Dispose file watchers
    foreach ($w in $watchers) { $w.EnableRaisingEvents = $false; $w.Dispose() }
    # Stop jobs
    Stop-Job  $script:backendJob -ErrorAction SilentlyContinue
    Stop-Job  $frontendJob       -ErrorAction SilentlyContinue
    Remove-Job $script:backendJob -Force -ErrorAction SilentlyContinue
    Remove-Job $frontendJob       -Force -ErrorAction SilentlyContinue
    Write-Host "Done." -ForegroundColor Green
}
