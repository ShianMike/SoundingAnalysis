# Sounding Analysis Tool

[![GitHub last commit](https://img.shields.io/github/last-commit/ShianMike/SoundingAnalysis?style=flat-square&color=blue)](https://github.com/ShianMike/SoundingAnalysis/commits/main)
[![GitHub stars](https://img.shields.io/github/stars/ShianMike/SoundingAnalysis?style=flat-square)](https://github.com/ShianMike/SoundingAnalysis/stargazers)
[![GitHub forks](https://img.shields.io/github/forks/ShianMike/SoundingAnalysis?style=flat-square)](https://github.com/ShianMike/SoundingAnalysis/forks)
[![Made with Python](https://img.shields.io/badge/Python-3.10+-3776AB?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
[![Made with React](https://img.shields.io/badge/React-18-61DAFB?style=flat-square&logo=react&logoColor=black)](https://react.dev/)
[![Deployed on Cloud Run](https://img.shields.io/badge/Cloud%20Run-deployed-4285F4?style=flat-square&logo=googlecloud&logoColor=white)](https://soundinganalysis-752306366750.asia-southeast1.run.app)
[![GitHub Pages](https://img.shields.io/badge/GitHub%20Pages-live-222?style=flat-square&logo=github&logoColor=white)](https://shianmike.github.io/SoundingAnalysis/)
[![License](https://img.shields.io/badge/license-Educational%20%2F%20Research-green?style=flat-square)](#license)
[![MetPy](https://img.shields.io/badge/MetPy-powered-orange?style=flat-square&logo=python&logoColor=white)](https://unidata.github.io/MetPy/)

A full-stack atmospheric sounding analysis application that fetches real upper-air data from multiple sources and produces comprehensive Skew-T / Hodograph analysis plots.

**Live site:** <https://shianmike.github.io/SoundingAnalysis/>

---

## Screenshots

### Dashboard
The main interface — sidebar control panel on the left with data source selection, station picker, date/time controls, and modification options. Results display on the right.

![Dashboard](docs/screenshots/dashboard.png)

### Sounding Plot (Dark Theme)
Full Skew-T Log-P diagram with hodograph, storm-relative wind/vorticity profiles, theta/theta-e panel, and computed parameters:

![Sounding Result](docs/screenshots/sounding-result.png)

### Parameters & Climatology
Computed thermodynamic, kinematic, and composite parameters displayed as color-coded cards with severity thresholds. Climatology percentile bars show how each parameter ranks against SPC severe-weather proximity sounding climatology:

![Parameters & Climatology](docs/screenshots/parameters-and-climo.png)

### Station Map with Radar & Warnings
Interactive dark-themed Leaflet map with station markers, animated radar overlays (composite + mosaic), NWS active weather warnings, SPC convective outlook polygons, and click-to-select functionality:

![Station Map with Radar](docs/screenshots/station-map-radar.png)

### VWP Time-Height Display
NEXRAD VAD wind barbs across time and height — searchable radar dropdown sorted by proximity, configurable time window (3–24h):

![VWP Display](docs/screenshots/vwp-display.png)

### Multi-Sounding Comparison
Compare up to 4 soundings side-by-side with full Skew-T plots and a parameter comparison table highlighting differences:

![Comparison View](docs/screenshots/comparison-view.png)

### Parameter Time-Series Charts
Plot CAPE, SRH, STP, shear, lapse rates, and more over a configurable date range (up to 14 days). Grouped parameter selector with dark-themed Recharts line charts:

![Time Series](docs/screenshots/time-series.png)

### Sounding History
Auto-saved last 20 soundings with one-click reload. Tabbed view for soundings and comparisons with relative timestamps:

![History Panel](docs/screenshots/history-panel.png)

### Light Theme
Toggle between dark and light themes — fully consistent styling across all panels and plots:

![Light Theme](docs/screenshots/example-sounding-light.png)

### Color-Blind Mode
Okabe-Ito / Wong 2011 color-safe palette for all plot traces:

![Color-Blind Mode](docs/screenshots/example-sounding-cb.png)

---

## Features

### Skew-T Log-P Diagram
- **Source-aware plot titles** — "OBSERVED UPPER-AIR SOUNDING", "HRRR FORECAST (INIT + FHOUR → VALID)", "RAP ANALYSIS", or "ACARS AIRCRAFT OBS" depending on data source
- **Temperature** & **Dewpoint** profiles
- **Wet-bulb temperature** trace
- **Virtual temperature** trace
- **SB / MU / ML parcel traces** with color-coded dashed lines
- **Downdraft (DCAPE) parcel** trace
- **CAPE/CIN shading** — red fill (CAPE) and blue fill (CIN) between SB parcel and environment
- Dry adiabats, moist adiabats, mixing ratio lines
- **0°C and -20°C highlighted isotherms**
- **Dendritic Growth Zone (DGZ)** shading (-12°C to -17°C)
- **Hail Growth Zone (HGZ)** shading (-10°C to -30°C)
- **PBL top** marker (mixed-layer depth)
- **Piecewise CAPE bars** — color-coded layer-by-layer CAPE visualization
- **OBS wind barbs** column (right edge)
- **VAD wind barbs** column (green, alongside OBS) — NEXRAD radar-derived winds
- Key level annotations (LCL, LFC, EL, Freezing level, WBZ)
- **Surface T/Td in °F** labels
- Height labels in km AGL

### Hodograph
- Color-coded by height (0–1, 1–3, 3–6, 6–9, 9–12 km layers)
- Bunkers right-mover, left-mover, and mean wind vectors
- Deviant Tornado Motion (DTM) marker
- Upshear / Downshear MCS motion markers
- **Critical angle** between 0–500m shear and storm-relative inflow
- VAD Wind Profiler overlay (green)
- Effective inflow layer SRH fill
- Dynamic bounds based on wind profile
- **Storm-relative mode** — shift all winds into the storm-relative frame (SM → origin crosshair)
- **Profile smoothing** — Gaussian filter (adjustable σ) to tame noisy ACARS profiles
- **Boundary line** — user-defined boundary orientation plotted on the hodograph

### Computed Parameters
| Category | Parameters |
|---|---|
| **Thermodynamic** | SB/ML/MU CAPE & CIN, DCAPE, DCIN, ECAPE, 3CAPE, 6CAPE, MU NCAPE, LCL height, LFC, EL, Lapse rates (0–3 km, 3–6 km), Precipitable water, Freezing level, Wet-bulb zero height, RH by layer |
| **Kinematic** | Bulk wind difference (500 m, 1 km, 3 km, 6 km), Effective BWD, Storm-relative helicity (500 m, 1 km, 3 km), Effective SRH, Storm-relative wind by layer, Bunkers storm motion, Critical angle, Corfidi Upwind/Downwind vectors |
| **Composite** | STP (fixed & effective), SCP, SHIP, DCP, BRN, ECAPE, NCAPE |
| **Hazard Assessment** | Tornado, Hail, Wind, Flood threat levels (HIGH / MOD / LOW / NONE) |
| **Convective Mode** | Predicted storm mode (Discrete Supercell → Multicell → Single Cell) via BRN + Thompson et al. 2007 framework |
| **Downburst** | WMSI, MDPI, Max Gust estimate |
| **Fire Weather** | Fosberg FWI, Haines Index, Hot-Dry-Windy Index |
| **Winter Wx** | Precip Type (Bourgouin 2000), Warm/Cold layer energy |
| **Temperature Advection** | WAA/CAA classification in 1 km layers (0–6 km) |

### Interactive Skew-T (Canvas)
- Full interactive Canvas-based Skew-T Log-P diagram alongside the static plot
- Pan, zoom, and hover with real-time readout of T, Td, pressure, and height
- Significant level markers: LCL, LFC, EL, Freezing Level, Wet-Bulb Zero
- Syncs with dark/light theme toggle

### Interactive Hodograph (Canvas)
- Canvas-based hodograph with height-colored wind trace (0–1, 1–3, 3–6, 6–9, 9–12 km)
- Bunkers RM/LM/MW storm motion vectors
- Hover readout showing wind speed, direction, and height at cursor
- Storm-relative mode support

### Sounding Animation
- Animate through BUFKIT forecast hours with play/pause/step controls
- Adjustable playback speed and timeline slider
- Watches Skew-T and parameters update frame-by-frame

### Additional Panels
- Storm-relative wind & streamwise vorticity profiles
- **Theta (θ) / Theta-e (θe) profile** — potential temperature and equivalent potential temperature vs height, with moisture gap fill and key height markers
- Comprehensive parameter text readout (thermodynamic + kinematic indices)
- **Climatology percentile comparison** — horizontal bar chart ranking each parameter against SPC severe-weather proximity sounding climatology (color-coded: grey <50th, green 50th–75th, orange 75th–90th, red >95th)
- **Predicted Convective Mode** — visual spectrum bar (Pulse → Multicell → Supercell) with active-mode highlight
- **Hazard Assessment** — always displays all 4 threat types (Tornado, Hail, Wind, Flood) with level or NONE

### VAD Wind Profile (VWP) Time-Height Display
- Standalone VWP panel showing wind barbs across time and height
- **Searchable radar dropdown** with all ~120 NEXRAD sites, sorted by proximity to selected station, with distance in km
- Fetches NEXRAD Level-III VAD data from Iowa Environmental Mesonet
- Configurable time range (3h, 6h, 12h, 24h)
- Dark-themed time-height section with color-coded wind speed

### Risk Scanner & Mesoscale Dashboard
- Scans all CONUS upper-air stations in parallel
- Scores sites by STP, SCP, SHIP, DCP / tornado potential
- Returns ranked results via API
- Auto-opens interactive station map with color-coded risk markers
- **Mesoscale panel:** sortable multi-station parameter table with color-coded threshold chips

### Interactive Station Map
- Dark-themed Leaflet map centered on CONUS
- Station markers color-coded by risk scan data (STP thresholds)
- Click station markers to select them
- Click anywhere on map to set lat/lon for RAP source
- Fly-to animation on station selection
- **Animated radar overlay:** composite reflectivity (RainViewer) and IEM US mosaic (N0Q) with play/pause/step controls and frame scrubbing (24 frames, 2-hour window)
- **Velocity overlay:** storm-relative velocity from nearest WSR-88D with dynamic product detection (prefers N0U, falls back to N0S), color legend (green toward / red away), radar site marker, 230 km range ring, and IEM scan timestamp display
- **NWS active warnings:** real-time Tornado/Severe Thunderstorm/Flash Flood/Special Weather warnings from the NWS API, rendered as color-coded polygons with click-for-details popups
- **SPC outlook overlays:** Day 1 through Day 8 convective outlook GeoJSON with color-coded risk categories and legend
- **Wind flow animation:** surface streamlines and 500 hPa steering flow overlays powered by Open-Meteo forecast data, rendered on a Canvas layer with animated particle trails
- **SPC watch boxes:** Tornado and Severe Thunderstorm watch outlines from NWS MapServer
- **Spotter Network:** real-time active storm spotter positions
- **Lightning overlay:** real-time lightning strikes via Blitzortung WebSocket feed
- **Proximity search:** click "Near me" to find the closest upper-air stations to your browser location
- **Wind barbs on map:** toggle surface/850 hPa/500 hPa Open-Meteo wind barbs over the map
- **Shear vectors:** 0-6 km bulk shear direction arrows on station markers (requires risk scan)
- **Favorite station stars:** gold star markers highlight your pinned stations on the map
- **ACARS airport markers:** blue aircraft icons on major ACARS-capable airports with click-to-fetch
- **Tools dropdown:** dedicated toolbar group (right side) with draw-on-map tools for storm motion vectors and boundary orientation lines
- **Draw Storm Motion:** two-click tool to define a custom storm-motion vector — first click sets start, cursor shows live speed (kt) and direction as you move, second click commits the vector and auto-fills the ControlPanel storm-motion fields
- **Draw Boundary:** two-click tool to define a boundary orientation line — live orientation angle preview while moving the cursor, auto-fills ControlPanel boundary orientation on commit
- **Live draw preview:** dashed guide line and real-time value readout (speed/direction or orientation) displayed on the map and toolbar between the first and second click

### Multi-Sounding Comparison
- Compare up to 4 soundings side-by-side
- Slot-based UI: pick station, source, and date for each sounding
- Side-by-side Skew-T plots in a responsive grid
- Full parameter comparison table with Δ (difference) column
- Highlights highest/lowest values across soundings
- **Download comparison** — tiles all plots into a single composite PNG
- **Comparison history** — auto-saved to localStorage; reload from History panel

### Parameter Time-Series Charts
- Plot CAPE, SRH, STP, shear, lapse rates over a date range
- Date range picker (up to 14 days)
- 00Z + 12Z resolution for ≤7 days, 12Z only for >7 days
- Grouped parameter selector with exclusive group selection
- Custom dark-themed Recharts line charts with tooltips

### Ensemble Sounding Plume
- Spaghetti plume by fetching multiple BUFKIT forecast hours
- T/Td traces at adjustable alpha for spread visualization
- Hodograph spread panel showing low-level wind variability
- Configurable hour presets: short (0–6h), medium (0–12h), long (0–24h), extended (0–48h)
- Supports both Iowa State archive and Penn State real-time feeds

### Sounding Modifications
- **Surface modification** — override surface T, Td, wind speed/direction and re-compute all parameters
- **Custom storm motion** — input direction + speed for SRH/SRW recalculation; also settable by drawing on the station map
- **Profile smoothing** — Gaussian filter with adjustable σ (great for noisy ACARS profiles)
- **Boundary orientation** — plot boundary line on hodograph at custom angle; also settable by drawing on the station map

### Sounding History & Favorites
- Auto-saves last 20 soundings to localStorage with one-click reload
- Relative timestamps ("3m ago", "2h ago")
- Tabbed view: Soundings and Comparisons tabs
- **Favorite stations** — pin frequently used stations; persisted in localStorage

### Sounding Timeline
- Horizontal scrollable bar showing available sounding times for the past 4 days
- 00Z and 12Z cycle buttons with one-click fetch
- Displayed above the results when a sounding is loaded

### Shareable Sounding Links
- Sounding parameters encoded in the URL query string
- Opening a link auto-fetches the sounding (station, source, date, model, etc.)
- "Copy Link" button in the results toolbar

### Export Formats
- **CSV** — all computed parameters in spreadsheet-ready format
- **JSON** — full sounding data and parameters in JSON format
- **SHARPpy** — raw profile data in SHARPpy-compatible format
- **CM1** — `input_sounding` format for Cloud Model 1 numerical simulations
- **Comparison CSV** — side-by-side parameter comparison export
- All exports available as one-click downloads from the results toolbar

### WRF / CM1 Data Ingestion
- Upload WRF netCDF output files (wrfout_d0x) directly from the Custom Upload page
- Automatic grid point extraction — specify target lat/lon or use domain center
- Supports WRF staggered grids: P/PB, PH/PHB, T (perturbation theta), QVAPOR, U/V
- Also accepts text-based formats: CSV, SHARPpy, CM1 input_sounding

### PSU BUFKIT Feed (Penn State)
- Real-time BUFKIT profiles from Penn State's e-wall server
- Latest model run (RAP, HRRR, NAM, NAM Nest, GFS, HiResW, SREF)
- Complementary to age-indexed Iowa State archive
- Separate "PSU" data source option in the control panel

### Auto-Refresh & Parameter Alerts
- **Auto-refresh polling** — configurable interval to re-fetch sounding data automatically
- **Parameter watch alerts** — highlight cards that changed since the last fetch with a yellow pulse animation

### PWA / Offline Support
- Installable as a Progressive Web App (manifest + service worker)
- Caches static assets for offline access to previously loaded data

### Theme & Accessibility
- **Dark theme** (default) and **Light theme** toggle with localStorage persistence
- **Color-blind mode** — Okabe-Ito/Wong 2011 color-safe palette
- **Keyboard shortcuts** — H=history, C=compare, M=map, T=trends, V=VWP, ?=help
- **Print layout** — optimized `@media print` stylesheet with 4-column compact parameter grid
- **Mobile-responsive** — breakpoints at 1024px, 768px, 480px; touch-friendly drag-to-pan on plots
- **Skew-T theme sync** — interactive Canvas plots automatically match dark/light mode

### Feedback System
- Built-in modal for Suggestion / Bug Report / Feature Request
- Submitted to backend API; logged to Cloud Run stdout

---

## Data Sources

| Source | Flag | Description |
|---|---|---|
| **Observed (IEM/UWyo)** | `--source obs` | Real radiosonde observations from CONUS upper-air sites |
| **RAP Analysis** | `--source rap` | Rapid Refresh model analysis at any lat/lon (requires `siphon`) |
| **BUFKIT Forecasts** | `--source bufkit` | HRRR, RAP, NAM, NAM-Nest, GFS, SREF forecasts from Iowa State |
| **ACARS/AMDAR** | `--source acars` | Aircraft observations at airports (IEM) |
| **PSU BUFKIT** | `--source psu` | Latest model run from Penn State's real-time feed (RAP, HRRR, NAM, GFS, etc.) |
| **Custom Upload** | (UI only) | Paste CSV/SHARPpy/CM1 text or upload WRF netCDF files |

---

## Project Structure

```
├── app.py               # Flask entry point (registers blueprints, SPA catch-all)
├── routes/              # Flask API route blueprints
│   ├── __init__.py        # Blueprint registration helper
│   ├── analysis.py        # /api/compare, /api/time-series, /api/ensemble-plume
│   ├── feedback.py        # /api/feedback (GET/POST)
│   ├── helpers.py         # Shared helpers (safe_round, parse_date, etc.)
│   ├── meta.py            # /api/stations, /api/sources, /api/acars-airports
│   ├── risk.py            # /api/risk-scan
│   ├── sounding_routes.py # /api/sounding (main endpoint)
│   ├── spc.py             # /api/spc-outlook
│   └── wind.py            # /api/vwp-display
├── sounding/            # Core sounding analysis package
│   ├── __init__.py        # Public API exports
│   ├── __main__.py        # CLI entry point (python -m sounding)
│   ├── cli.py             # Interactive CLI menu
│   ├── constants.py       # Station database, model metadata
│   ├── fetchers.py        # Data fetching (OBS, RAP, BUFKIT, ACARS, PSU)
│   ├── merge.py           # Profile merging utilities
│   ├── parameters.py      # Thermodynamic & kinematic parameter computation
│   ├── plotting.py        # Skew-T, hodograph, parameter panel plotting
│   ├── tornado.py         # Tornado risk scanning logic
│   ├── utils.py           # Time/station utilities
│   └── vad.py             # VAD Wind Profile fetching & display
├── Dockerfile           # Docker config for Cloud Run deployment
├── gunicorn.conf.py     # Gunicorn config (reads PORT from env)
├── requirements.txt     # Python dependencies
├── deploy.ps1           # Build frontend + deploy to GitHub Pages
├── deploy-cloudrun.ps1  # Deploy backend to Google Cloud Run
├── .gcloudignore        # Exclude dev files from Cloud Build uploads
└── frontend/            # React + Vite frontend
    ├── src/
    │   ├── App.jsx              # Main app shell
    │   ├── api.js               # API client with timeouts & retries
    │   ├── history.js           # localStorage sounding history
    │   ├── favorites.js         # localStorage station favorites
    │   └── components/
    │       ├── ControlPanel.jsx   # Dashboard sidebar
    │       ├── Header.jsx         # Feedback modal
    │       ├── ResultsView.jsx    # Plot image + parameter display + map
    │       ├── StationMap.jsx     # Interactive Leaflet station map
    │       ├── HistoryPanel.jsx   # Sounding history sidebar
    │       ├── TimeSeriesChart.jsx # Parameter time-series charts
    │       ├── ComparisonView.jsx  # Multi-sounding comparison
    │       ├── VwpDisplay.jsx     # VWP time-height display
    │       ├── MesoPanel.jsx      # Mesoscale analysis table
    │       ├── EnsemblePlume.jsx   # Ensemble sounding plume
    │       ├── CustomUpload.jsx   # WRF/CSV/SHARPpy upload page
    │       ├── InteractiveSkewT.jsx # Canvas interactive Skew-T
    │       ├── InteractiveHodograph.jsx # Canvas interactive hodograph
    │       ├── SoundingTimeline.jsx # Historical sounding time bar
    │       ├── SoundingAnimator.jsx # BUFKIT forecast animation
    │       └── *.css              # Component styles (dark/light theme)
    ├── public/
    │   ├── favicon.svg
    │   ├── manifest.json        # PWA manifest
    │   ├── sw.js                # Service worker for offline
    │   └── og-thumbnail.png     # Open Graph preview image
    ├── package.json
    └── vite.config.js
```

---

## Quick Start

### Backend (Python)

```bash
# Install dependencies
pip install -r requirements.txt

# Optional for RAP source:
pip install siphon
```

**CLI usage (standalone):**

```bash
python -m sounding                                            # Interactive menu
python -m sounding --station OUN                              # Latest observed
python -m sounding --station OUN --date 2024061200            # Specific date/time
python -m sounding --source rap --lat 36.4 --lon -99.4        # RAP at any point
python -m sounding --source bufkit --model hrrr --station OUN # HRRR forecast
python -m sounding --source acars --station KDFW              # Aircraft obs
python -m sounding --list-sources                             # Show all sources
```

**API server:**

```bash
python app.py
# → http://localhost:5000
```

### Frontend (React)

```bash
cd frontend
npm install
npm run dev
# → http://localhost:5173
```

Set `VITE_API_URL` to point to your backend (defaults to localhost:5000).

---

## Deployment

| Component | Platform | URL |
|---|---|---|
| **Backend API** | Google Cloud Run | `https://soundinganalysis-752306366750.asia-southeast1.run.app` |
| **Frontend** | GitHub Pages | `https://shianmike.github.io/SoundingAnalysis/` |

To redeploy frontend:

```powershell
.\deploy.ps1
```

To redeploy backend:

```powershell
.\deploy-cloudrun.ps1
```

The frontend build uses Vite and force-pushes to the `gh-pages` branch.
The backend deploys via Cloud Build from source to Cloud Run (Singapore region).

---

## API Endpoints

| Method | Endpoint | Description |
|---|---|---|
| `GET` | `/api/stations` | List all known sounding stations |
| `GET` | `/api/sources` | List available data sources and BUFKIT models |
| `GET` | `/api/acars-airports` | List ACARS-capable airport locations |
| `POST` | `/api/sounding` | Fetch sounding, compute params, return base64 plot + data |
| `POST` | `/api/risk-scan` | Scan stations and return tornado risk scores |
| `POST` | `/api/time-series` | Fetch parameter trends over a date range |
| `POST` | `/api/compare` | Fetch multiple soundings for side-by-side comparison |
| `POST` | `/api/feedback` | Submit user feedback/suggestion |
| `GET`  | `/api/feedback` | Retrieve all submitted feedback |
| `GET`  | `/api/vwp-display` | Fetch VAD Wind Profile time-height display image |

### `POST /api/sounding` — Request Body

```json
{
  "source": "obs",
  "station": "OUN",
  "date": "2024061200",
  "lat": 35.22,
  "lon": -97.46,
  "model": "hrrr",
  "fhour": 0
}
```

---

## Security

- **HTTPS:** forced in production via Flask-Talisman (HSTS 2-year preload)
- **Content Security Policy:** restrictive CSP with nonce-based script-src
- **CORS:** locked to production origins only (localhost allowed when `FLASK_DEBUG` is set)
- **Rate limiting:** Flask-Limiter — 200 req/min global, 30 req/sec burst, 10/min on feedback, 30/min on sounding
- **Security headers:** X-Content-Type-Options, X-Frame-Options DENY, COOP, CORP, Permissions-Policy, Referrer-Policy
- **Input validation:** path-traversal blocking, 16 MB request-size limit
- **Cloud Armor:** optional WAF setup script with OWASP CRS rules and IP rate limiting

---

## Dependencies

- **Python:** Flask, Flask-Talisman, Flask-Limiter, Flask-CORS, MetPy, Matplotlib, NumPy, Requests
- **Frontend:** React 18, Vite, Leaflet, React-Leaflet, Recharts, D3, Lucide React
- **Optional:** siphon (RAP source)

---

## License

This project is for educational and research purposes.
