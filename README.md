# Sounding Analysis Tool

A full-stack atmospheric sounding analysis application that fetches real upper-air data from multiple sources and produces comprehensive Skew-T / Hodograph analysis plots — similar to [SounderPy](https://github.com/kylejgillett/sounderpy).

**Live site:** <https://shianmike.github.io/SoundingAnalysis/>

---

## Features

### Skew-T Log-P Diagram
- **Temperature** & **Dewpoint** profiles
- **Wet-bulb temperature** trace
- **Virtual temperature** trace
- **Downdraft (DCAPE) parcel** trace
- Dry adiabats, moist adiabats, mixing ratio lines
- Wind barbs
- Key level annotations (LCL, LFC, EL, Freezing level)
- Height labels in km AGL

### Hodograph
- Color-coded by height (0–1, 1–3, 3–6, 6–9 km layers)
- Bunkers right-mover, left-mover, and mean wind vectors
- Dynamic bounds based on wind profile

### Computed Parameters
| Category | Parameters |
|---|---|
| **Thermodynamic** | SB/ML/MU CAPE & CIN, DCAPE, LCL height, LFC, EL, Lapse rates (0–3 km, 3–6 km), Precipitable water, Freezing level, Wet-bulb zero height, RH by layer |
| **Kinematic** | Bulk wind difference (500 m, 1 km, 3 km, 6 km), Storm-relative helicity (500 m, 1 km, 3 km), Bunkers storm motion |
| **Composite** | Significant Tornado Parameter (STP), Supercell Composite Parameter (SCP), Significant Hail Parameter (SHIP), Derecho Composite Parameter (DCP) |

### Additional Panels
- Storm-relative wind & streamwise vorticity profiles
- Comprehensive parameter text readout (thermodynamic + kinematic indices)

### Risk Scanner
- Scans all CONUS upper-air stations in parallel
- Scores sites by STP, SCP, SHIP, DCP / tornado potential
- Returns ranked results via API
- Auto-opens interactive station map with color-coded risk markers

### Interactive Station Map
- Dark-themed Leaflet map centered on CONUS
- Station markers color-coded by risk scan data (STP thresholds)
- Click station markers to select them
- Click anywhere on map to set lat/lon for RAP/ERA5 sources
- Fly-to animation on station selection

### SPC Convective Outlook Overlay
- Day 1, Day 2, and Day 3 outlook overlays on the station map
- Fetches live GeoJSON from SPC
- Color-coded risk categories (General Thunder, Marginal, Slight, Enhanced, Moderate, High)
- Legend with risk percentages (tornado, severe, wind)
- Non-blocking overlays — station markers remain clickable through outlook polygons

### Sounding History
- Automatically saves last 20 soundings to localStorage
- Quick-load previous soundings with one click
- Relative timestamps ("3m ago", "2h ago")
- Tabbed view: **Soundings** and **Comparisons** tabs

### Favorite Stations
- Star icon to pin frequently used stations
- Persisted in localStorage
- Sort stations by favorites

### Dashboard Sidebar Layout
- Fixed left sidebar with branding, controls, and footer
- Sidebar includes data source selection, station picker, date/time controls
- Footer with Feedback and GitHub links
- History panel slides in beside the sidebar without overlapping content
- Responsive: collapses to stacked layout on mobile (≤1024px)

### Feedback & Suggestions
- Built-in feedback modal (Suggestion / Bug Report / Feature Request)
- Submitted to backend API for developer review
- Feedback logged to Cloud Run stdout for persistence across ephemeral containers

### Parameter Time-Series Charts
- Plot CAPE, SRH, STP, shear, lapse rates over a date range
- Date range picker (up to 14 days)
- 00Z + 12Z resolution for ≤7 days, 12Z only for >7 days
- Grouped parameter selector with exclusive group selection
- Custom dark-themed Recharts line charts with tooltips

### Multi-Sounding Comparison
- Compare up to 4 soundings side-by-side
- Slot-based UI: pick station, source, and date for each sounding
- Side-by-side Skew-T plots in a responsive grid
- Full parameter comparison table with Δ (difference) column
- Highlights highest/lowest values across soundings
- 00Z/12Z toggle for observed/ACARS sources
- **Download comparison** — Tiles all sounding plots into a single composite PNG
- **Comparison history** — Automatically saved to localStorage; reload previous comparisons from the History panel's "Comparisons" tab

---

## Data Sources

| Source | Flag | Description |
|---|---|---|
| **Observed (IEM/UWyo)** | `--source obs` | Real radiosonde observations from CONUS upper-air sites |
| **RAP Analysis** | `--source rap` | Rapid Refresh model analysis at any lat/lon (requires `siphon`) |
| **BUFKIT Forecasts** | `--source bufkit` | HRRR, RAP, NAM, NAM-Nest, GFS, SREF forecasts from Iowa State |
| **ERA5 Reanalysis** | `--source era5` | Global reanalysis at any lat/lon (requires `cdsapi` + CDS key) |
| **ACARS/AMDAR** | `--source acars` | Aircraft observations at airports (IEM) |

---

## Project Structure

```
├── sounding.py          # Core: data fetching, parameter computation, plotting
├── app.py               # Flask API serving the React frontend
├── Dockerfile           # Docker config for Cloud Run deployment
├── gunicorn.conf.py     # Gunicorn config (reads PORT from env)
├── requirements.txt     # Python dependencies
├── deploy.ps1           # Build frontend + deploy to GitHub Pages
├── deploy-cloudrun.ps1  # Deploy backend to Google Cloud Run
├── .gcloudignore        # Exclude frontend from Cloud Build uploads
├── feedback.json        # Server-side feedback storage
└── frontend/            # React + Vite frontend
    ├── src/
    │   ├── App.jsx              # Main app shell
    │   ├── api.js               # API client with timeouts & retries
    │   ├── history.js           # localStorage sounding history
    │   ├── favorites.js         # localStorage station favorites
    │   └── components/
    │       ├── ControlPanel.jsx   # Dashboard sidebar: branding, controls, footer
    │       ├── Header.jsx         # Feedback modal (portal)
    │       ├── ResultsView.jsx    # Plot image + parameter display + map
    │       ├── StationMap.jsx     # Interactive Leaflet station map
    │       ├── HistoryPanel.jsx   # Sounding history sidebar
    │       ├── TimeSeriesChart.jsx # Recharts time-series parameter charts
    │       ├── ComparisonView.jsx  # Multi-sounding comparison view
    │       └── *.css              # Component styles (dark theme)
    ├── public/
    │   └── favicon.svg          # Custom skew-T favicon
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

# Optional for ERA5 source:
pip install cdsapi netCDF4
```

**CLI usage (standalone):**

```bash
python sounding.py                                            # Interactive menu
python sounding.py --station OUN                              # Latest observed
python sounding.py --station OUN --date 2024061200            # Specific date/time
python sounding.py --source rap --lat 36.4 --lon -99.4        # RAP at any point
python sounding.py --source bufkit --model hrrr --station OUN # HRRR forecast
python sounding.py --source era5 --lat 36.4 --lon -99.4 --date 2020050100
python sounding.py --source acars --station KDFW              # Aircraft obs
python sounding.py --list-sources                             # Show all sources
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
| `POST` | `/api/sounding` | Fetch sounding, compute params, return base64 plot + data |
| `POST` | `/api/risk-scan` | Scan stations and return tornado risk scores |
| `POST` | `/api/time-series` | Fetch parameter trends over a date range |
| `POST` | `/api/compare` | Fetch multiple soundings for side-by-side comparison |
| `POST` | `/api/feedback` | Submit user feedback/suggestion |
| `GET`  | `/api/feedback` | Retrieve all submitted feedback |

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

## Dependencies

- **Python:** Flask, MetPy, Matplotlib, NumPy, Requests
- **Frontend:** React 18, Vite, Lucide React, Leaflet, React-Leaflet, Recharts
- **Optional:** siphon (RAP), cdsapi + netCDF4 (ERA5)

---

## Feature Comparison: Ours vs SounderPy (Gillett)

A side-by-side look at what our Sounding Analysis Tool offers compared to [SounderPy](https://github.com/kylejgillett/sounderpy) (v3.1.0) by Kyle J. Gillett.

### Data Sources

| Capability | Ours | SounderPy |
|---|:---:|:---:|
| Observed radiosondes (UWyo / IEM) | ✅ | ✅ |
| RAP model analysis (THREDDS) | ✅ | ✅ |
| BUFKIT forecast models (ISU + PSU) | ✅ (ISU) | ✅ (ISU + PSU) |
| ERA5 global reanalysis (CDS) | ✅ | ✅ |
| ACARS aircraft observations | ✅ (IEM) | ✅ (OU archive) |
| RUC reanalysis (2005–2020) | ❌ | ✅ |
| NCEP-FNL reanalysis | ❌ | ✅ |
| IGRAv2 global obs archive (1905–present) | ❌ | ✅ |
| WRF output ingestion | ❌ | ✅ |
| CM1 input_sounding ingestion | ❌ | ✅ |
| Custom user-defined data dict | ❌ | ✅ |
| Omega / vertical velocity (model data) | ❌ | ✅ |

### Plotting & Visualization

| Capability | Ours | SounderPy |
|---|:---:|:---:|
| Skew-T Log-P diagram | ✅ | ✅ |
| Hodograph (color-coded by height) | ✅ | ✅ |
| Dark mode | ✅ (default) | ✅ (optional) |
| Light mode | ❌ | ✅ (default) |
| Color-blind mode | ❌ | ✅ |
| Wet-bulb temperature trace | ✅ | ✅ |
| Virtual temperature trace | ✅ | ✅ |
| Downdraft (DCAPE) parcel trace | ✅ | ❌ |
| Wind barbs on Skew-T | ✅ | ✅ |
| Key level annotations (LCL, LFC, EL) | ✅ | ✅ |
| Height labels (km AGL) | ✅ | ✅ |
| Radar reflectivity overlay (mosaic / single-site) | ❌ | ✅ |
| Map inset with station location | ✅ (CONUS outline) | ✅ (with radar) |
| Map inset zoom control | ❌ | ✅ |
| Storm-relative wind profile panel | ✅ | ✅ |
| Streamwise vorticity / streamwiseness panel | ✅ | ✅ |
| Piecewise CAPE / stepwise CAPE-CIN plot | ❌ | ✅ |
| Theta / Theta-e profile toggle | ❌ | ✅ |
| Hodograph boundary lines | ❌ | ✅ |
| Composite sounding (overlay multiple profiles) | ❌ | ✅ |
| VAD (radar wind profiler) hodograph | ❌ | ✅ |
| Storm-relative hodograph mode | ❌ | ✅ |

### Parcel & Parameter Calculations

| Capability | Ours | SounderPy |
|---|:---:|:---:|
| SB / MU / ML CAPE & CIN | ✅ | ✅ |
| DCAPE | ✅ | ✅ |
| ECAPE (entraining CAPE) | ❌ | ✅ |
| Irreversible adiabatic ascent parcels | ❌ | ✅ |
| Pseudoadiabatic ascent parcels | ✅ (implicit) | ✅ (explicit) |
| STP, SCP, SHIP, DCP composites | ✅ | ✅ (via SHARPpy) |
| SRH (500 m, 1 km, 3 km) | ✅ | ✅ |
| Bulk wind difference (shear layers) | ✅ | ✅ |
| Bunkers storm motion (RM, LM, MW) | ✅ | ✅ |
| Custom storm motion input | ❌ | ✅ |
| Surface modification (override sfc T/Td/wind) | ❌ | ✅ |
| Lapse rates (0–3, 3–6 km) | ✅ | ✅ |
| Precipitable water, freezing level, WBZ | ✅ | ✅ |
| RH by layer | ✅ | ✅ |

### Helper & Utility Tools

| Capability | Ours | SounderPy |
|---|:---:|:---:|
| Export to CSV | ❌ | ✅ |
| Export to CM1 format | ❌ | ✅ |
| Export to SHARPpy format | ❌ | ✅ |
| Profile merging (weighted average) | ❌ | ✅ |
| Profile smoothing (Gaussian) | ❌ | ✅ |
| Vertical interpolation utility | ✅ (internal, 100 m) | ✅ (exposed API) |
| Print variables to console | ✅ (CLI) | ✅ |
| Find nearest station | ✅ | ✅ |
| Lat/lon lookup by station type | ❌ | ✅ |

### Web Application Features (Ours Only)

These are features that exist in our web app but have **no equivalent in SounderPy** (which is a Python library, not a web app):

| Feature | Description |
|---|---|
| **Full-stack web app** | React frontend + Flask API, deployed on GitHub Pages + Cloud Run |
| **Interactive station map** | Dark Leaflet map with click-to-select stations & lat/lon picking |
| **SPC outlook overlays** | Day 1/2/3 convective outlook GeoJSON on the map |
| **Risk scanner** | Parallel CONUS-wide scan ranking stations by STP/SCP/SHIP/DCP |
| **Multi-sounding comparison** | Up to 4 side-by-side Skew-T plots with parameter diff table |
| **Parameter time-series charts** | Plot CAPE, SRH, STP, shear over a date range (up to 14 days) |
| **Sounding history** | Auto-saved last 20 soundings in localStorage, one-click reload |
| **Favorite stations** | Pinned stations persisted in localStorage |
| **Comparison history** | Auto-saved comparisons with reload from History panel |
| **Feedback system** | In-app bug report / feature request modal |
| **Download comparison PNG** | Tiles all compared sounding plots into a single image |

---

## Upcoming Features

### High Priority — Parity with SounderPy

- [x] **ECAPE (Entraining CAPE)** — Entraining CAPE calculations using the Peters et al. (2023) formulation; displayed on Skew-T and in parameter tables ✅
- [x] **IGRAv2 Global Obs Archive** — IGRAv2 global radiosonde data via UWyo with auto-region detection; 15+ international stations included ✅
- [x] **Radar Reflectivity Overlay** — NEXRAD mosaic reflectivity (IEM) toggle on the station map ✅
- [x] **Composite Sounding Overlay** — Overlay multiple T/Td profiles + hodographs on a single Skew-T via the Comparison panel ✅
- [x] **Surface Modification** — Override surface T, Td, wind speed/direction and re-compute all parameters ✅
- [x] **Custom Storm Motion** — Input custom storm motion (direction + speed) for SRH/SRW recalculation ✅
- [x] **Piecewise CAPE / Stepwise CAPE-CIN** — Layer-by-layer CAPE/CIN data computed and returned via API ✅

### Medium Priority — New Capabilities

- [x] **Export Parameters to CSV** — One-click CSV download of all computed thermodynamic, kinematic, and moisture parameters ✅
- [ ] **Export to SHARPpy / CM1 Format** — Save sounding data in formats compatible with SHARPpy and CM1
- [ ] **Profile Merging** — Weighted-average two soundings together to create a blended analysis profile
- [ ] **Profile Smoothing** — Apply Gaussian smoothing to noisy profiles (especially ACARS)
- [ ] **RUC & NCEP-FNL Reanalysis** — Add RUC (2005–2020) and NCEP-FNL (2005–2020) as additional reanalysis sources
- [ ] **VAD Wind Profiler Hodograph** — Fetch and plot NEXRAD VAD (Velocity Azimuth Display) wind data on the hodograph
- [ ] **Storm-Relative Hodograph Mode** — Toggle the hodograph between ground-relative and storm-relative frames
- [ ] **Hodograph Boundary Lines** — Draw user-defined boundary lines on the hodograph at custom angles
- [ ] **Color-Blind Mode** — Swap dewpoint trace from green/blue and adjust color palettes for accessibility
- [ ] **Light Theme Toggle** — Add a light-mode theme option alongside the current dark theme
- [ ] **Theta / Theta-e Profile Panel** — Add an optional panel showing potential temperature and equivalent potential temperature profiles

### Lower Priority — Quality of Life

- [ ] **Custom Station Groups** — Create and save named groups of stations for quick batch analysis
- [ ] **WRF / CM1 Data Ingestion** — Allow uploading WRF output or CM1 input_sounding files for analysis
- [ ] **Custom Data Upload** — Let users paste or upload raw sounding data (CSV/text) for plotting
- [ ] **PSU BUFKIT Feed** — Add Penn State's real-time BUFKIT feed as a secondary forecast model source
- [ ] **Map Inset Zoom Control** — Adjustable zoom level for the CONUS mini-map on sounding plots
- [ ] **Animated Time-Series Playback** — Step through time-series soundings as an animation
- [ ] **Mobile-Optimized Sounding View** — Pinch-zoom and swipe-friendly sounding plot rendering on mobile
- [ ] **Shareable Sounding Links** — Encode station/source/date in URL for direct sharing

---

## License

This project is for educational and research purposes.
