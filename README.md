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

## Upcoming Features

- [ ] **Export Parameters to CSV** — One-click download of computed parameters and risk scan tables as CSV
- [ ] **Custom Station Groups** — Create and save named groups of stations for quick batch analysis
- [ ] **Dark/Light Theme Toggle** — Switchable color theme

---

## License

This project is for educational and research purposes.
