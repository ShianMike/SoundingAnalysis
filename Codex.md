# Codex Project Context

## Overview

This repository is a full-stack atmospheric sounding analysis app:

- Backend: Flask API in `app.py` plus blueprints in `routes/`
- Frontend: React 18 + Vite app in `frontend/`
- Domain logic: meteorological fetching, merging, plotting, and parameter calculations in `sounding/`
- Production/static bundle: emitted into `static/` and served by Flask as a SPA

The app analyzes upper-air soundings, forecast soundings, point soundings, risk scans, SPC overlays, VWP displays, comparisons, and time-series views.

## Key Entry Points

- Backend app bootstrap: `app.py`
- Blueprint registration: `routes/__init__.py`
- Frontend root component: `frontend/src/App.jsx`
- Frontend API wrapper: `frontend/src/api.js`
- Local dev launcher: `dev.ps1`

## Architecture Notes

### Backend

- `app.py` is intentionally slim:
  - configures CORS, rate limiting, CSP/Talisman, upload size limits, and SPA file serving
  - registers all blueprints from `routes/`
- Route logic is split by concern:
  - `routes/sounding_routes.py`: sounding fetch/analysis flows
  - `routes/analysis.py`: comparison/time-series style analysis endpoints
  - `routes/risk.py`: risk scan endpoints
  - `routes/spc.py`, `routes/outlook.py`, `routes/overlays.py`: outlook and map overlay data
  - `routes/wind.py`: wind-field and VWP-related data
  - `routes/feedback.py`, `routes/meta.py`: supporting endpoints
- Core meteorology and plotting logic lives in `sounding/`, not in the route handlers.

### Frontend

- `frontend/src/App.jsx` is a large stateful orchestrator for:
  - initial data loading
  - URL param sync
  - section visibility
  - submit/fetch flows
  - handoff into major feature panels
- Components under `frontend/src/components/` contain most UI features.
- `frontend/src/api.js` is the network boundary for the SPA.
- Vite dev server proxies `/api` to `http://localhost:5001`.

### Build / Static Serving

- Flask serves the built SPA from `static/`.
- Do not manually edit hashed files under `static/assets/`; they are generated build artifacts.
- `static/index.html`, `static/sw.js`, and related files are generated from the frontend build and may change whenever the frontend is rebuilt.

## Local Development

Preferred dev workflow:

```powershell
.\dev.ps1
```

This script:

- creates/uses `.venv`
- installs Python requirements
- installs frontend npm dependencies if missing
- runs Flask on `http://localhost:5001`
- runs Vite on `http://localhost:3000`
- proxies frontend `/api` requests to the Flask backend

Useful direct commands:

```powershell
# Backend
.venv\Scripts\python -m flask --app app run --host 127.0.0.1 --port 5001 --reload

# Frontend
cd frontend
npm install
npm run dev
npm run build
npm run lint
```

## Dependencies

### Python

Important backend/runtime libraries from `requirements.txt`:

- Flask
- flask-cors
- flask-limiter
- flask-talisman
- MetPy
- matplotlib
- numpy
- scipy
- requests
- siphon
- netCDF4
- shapely

### Frontend

Important frontend dependencies from `frontend/package.json`:

- React 18
- Vite
- D3
- Leaflet / React Leaflet
- Recharts
- hls.js
- lucide-react

## Editing Guidance

- Prefer editing source files in `frontend/src/`, `routes/`, and `sounding/`.
- Treat `static/assets/*` as build output unless the task is explicitly about emitted artifacts.
- If a UI or API change is made, verify both sides of the contract:
  - request payload shape in `frontend/src/api.js` or the calling component
  - response handling in the matching Flask route
- Keep `app.py` focused on app wiring and security concerns; avoid moving business logic into it.
- Put meteorological calculations and transformations into `sounding/` when possible.

## Current Repo State

This repository may already contain user changes. At the time this context file was generated, `git status --short` showed modified files including:

- `app.py`
- `frontend/index.html`
- `frontend/public/sw.js`
- `static/index.html`
- `static/sw.js`

Do not overwrite unrelated user changes when making new edits.

## Verification

There is no obvious automated backend test suite in the repository root. For changes, prefer lightweight verification such as:

```powershell
cd frontend
npm run lint
```

and manual end-to-end validation through `.\dev.ps1`.

## Practical Rules For Future Agents

- Read `README.md` first for feature intent and user-facing behavior.
- Check `frontend/src/App.jsx` before changing cross-panel flows or shared state.
- Check `frontend/src/api.js` before changing endpoint names or payloads.
- Check `routes/__init__.py` and the relevant blueprint before adding a new endpoint.
- Avoid editing generated files in `static/` unless you are intentionally updating build output.
