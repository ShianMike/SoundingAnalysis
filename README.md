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
| **Composite** | Significant Tornado Parameter (STP) |

### Additional Panels
- Storm-relative wind & streamwise vorticity profiles
- Comprehensive parameter text readout (thermodynamic + kinematic indices)

### Risk Scanner
- Scans all CONUS upper-air stations in parallel
- Scores sites by STP / tornado potential
- Returns ranked results via API

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
├── gunicorn.conf.py     # Gunicorn config for Render deployment
├── requirements.txt     # Python dependencies
├── deploy.ps1           # Build frontend + deploy to GitHub Pages
└── frontend/            # React + Vite frontend
    ├── src/
    │   ├── App.jsx              # Main app shell
    │   ├── api.js               # API client
    │   └── components/
    │       ├── ControlPanel.jsx  # Station/source/date controls
    │       ├── Header.jsx        # App header
    │       └── ResultsView.jsx   # Plot image + parameter display
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
| **Backend API** | Render | `https://soundinganalysis-1.onrender.com` |
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
- **Frontend:** React 18, Vite, Lucide React
- **Optional:** siphon (RAP), cdsapi + netCDF4 (ERA5)

---

## License

This project is for educational and research purposes.
