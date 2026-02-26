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

### Sounding History
- Automatically saves last 20 soundings to localStorage
- Quick-load previous soundings with one click
- Relative timestamps ("3m ago", "2h ago")

### Favorite Stations
- Star icon to pin frequently used stations
- Persisted in localStorage
- Sort stations by favorites

### Feedback & Suggestions
- Built-in feedback modal (Suggestion / Bug Report / Feature Request)
- Submitted to backend API for developer review

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
├── Dockerfile           # Docker config for Koyeb deployment
├── gunicorn.conf.py     # Gunicorn config for production
├── requirements.txt     # Python dependencies
├── deploy.ps1           # Build frontend + deploy to GitHub Pages
├── feedback.json        # Server-side feedback storage
└── frontend/            # React + Vite frontend
    ├── src/
    │   ├── App.jsx              # Main app shell
    │   ├── api.js               # API client with timeouts & retries
    │   ├── history.js           # localStorage sounding history
    │   ├── favorites.js         # localStorage station favorites
    │   └── components/
    │       ├── ControlPanel.jsx  # Station/source/date controls + favorites
    │       ├── Header.jsx        # App header with feedback & GitHub links
    │       ├── ResultsView.jsx   # Plot image + parameter display + map
    │       ├── StationMap.jsx    # Interactive Leaflet station map
    │       ├── HistoryPanel.jsx  # Sounding history sidebar
    │       └── *.css             # Component styles (dark theme)
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
| **Backend API** | Koyeb | `https://constitutional-lissi-mypersonalprojs-de5b9491.koyeb.app` |
| **Frontend** | GitHub Pages | `https://shianmike.github.io/SoundingAnalysis/` |

To redeploy:

```powershell
.\deploy.ps1
```

This builds the frontend with Vite, then force-pushes to the `gh-pages` branch.

---

## API Endpoints

| Method | Endpoint | Description |
|---|---|---|
| `GET` | `/api/stations` | List all known sounding stations |
| `GET` | `/api/sources` | List available data sources and BUFKIT models |
| `POST` | `/api/sounding` | Fetch sounding, compute params, return base64 plot + data |
| `POST` | `/api/risk-scan` | Scan stations and return tornado risk scores |
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
- **Frontend:** React 18, Vite, Lucide React, Leaflet, React-Leaflet
- **Optional:** siphon (RAP), cdsapi + netCDF4 (ERA5)

---

## Upcoming Features

- [ ] **Parameter Time-Series Charts** — Plot CAPE, SRH, STP trends over multiple sounding times for a station (e.g. last 5 days of 00Z/12Z)
- [ ] **Multi-Sounding Comparison** — Load two soundings side-by-side (00Z vs 12Z, or two stations) with split-view parameter cards
- [ ] **Export Parameters to CSV** — One-click download of computed parameters and risk scan tables as CSV
- [ ] **Severe Weather Outlook Overlay** — Display SPC Day 1/2/3 convective outlooks on the station map
- [ ] **Custom Station Groups** — Create and save named groups of stations for quick batch analysis
- [ ] **Dark/Light Theme Toggle** — Switchable color theme

---

## License

This project is for educational and research purposes.
