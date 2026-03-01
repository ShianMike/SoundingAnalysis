import { useState, useRef, useEffect } from "react";
import {
  Search,
  MapPin,
  Calendar,
  Database,
  Layers,
  Clock,
  Loader2,
  ChevronDown,
  ChevronRight,
  Zap,
  ArrowUpDown,
  History,
  Map,
  Star,
  TrendingUp,
  GitCompareArrows,
  MessageSquarePlus,
  Github,
  Thermometer,
  Wind,
  RotateCcw,
  Crosshair,
  Waves,
  Radio,
  Sun,
  Moon,
  Eye,
  Upload,
  Minus,
  X,
} from "lucide-react";
import { fetchRiskScan } from "../api";
import { getFavorites, toggleFavorite, getStationGroups, saveStationGroup, deleteStationGroup } from "../favorites";
import "./ControlPanel.css";

/* ── NEXRAD radar sites for VAD nearest-radar lookup ──────── */
const NEXRAD_SITES = [
  ["KABR",45.46,-98.41],["KABX",35.15,-106.82],["KAKQ",36.98,-77.01],
  ["KAMA",35.23,-101.71],["KAMX",25.61,-80.41],["KAPX",44.91,-84.72],
  ["KARX",43.82,-91.19],["KATX",48.19,-122.50],["KBBX",39.50,-121.63],
  ["KBGM",42.20,-75.98],["KBMX",33.17,-86.77],["KBOX",41.96,-71.14],
  ["KBRO",25.92,-97.42],["KBUF",42.95,-78.74],["KBYX",24.60,-81.70],
  ["KCAE",33.95,-81.12],["KCBW",46.04,-67.81],["KCBX",43.49,-116.24],
  ["KCCX",40.92,-78.00],["KCLE",41.41,-81.86],["KCLX",32.66,-81.04],
  ["KCRP",27.78,-97.51],["KCXX",44.51,-73.17],["KCYS",41.15,-104.81],
  ["KDAX",38.50,-121.68],["KDDC",37.76,-99.97],["KDFX",29.27,-100.28],
  ["KDGX",32.28,-89.98],["KDIX",39.95,-74.41],["KDLH",46.84,-92.21],
  ["KDMX",41.73,-93.72],["KDOX",38.83,-75.44],["KDTX",42.70,-83.47],
  ["KDVN",41.61,-90.58],["KDYX",32.54,-99.25],["KEAX",38.81,-94.26],
  ["KEMX",31.89,-110.63],["KENX",42.59,-74.06],["KEOX",31.46,-85.46],
  ["KEPZ",31.87,-106.70],["KESX",35.70,-114.89],["KEVX",30.56,-85.92],
  ["KEWX",29.70,-98.03],["KEYX",35.10,-117.56],["KFCX",37.02,-80.27],
  ["KFDR",34.36,-98.98],["KFDX",34.64,-103.63],["KFFC",33.36,-84.57],
  ["KFSD",43.59,-96.73],["KFSX",34.57,-111.20],["KFTG",39.79,-104.55],
  ["KFWS",32.57,-97.30],["KGGW",48.21,-106.63],["KGJX",39.06,-108.21],
  ["KGLD",39.37,-101.70],["KGRB",44.50,-88.11],["KGRK",30.72,-97.38],
  ["KGRR",42.89,-85.54],["KGSP",34.88,-82.22],["KGWX",33.90,-88.33],
  ["KGYX",43.89,-70.26],["KHDX",33.08,-106.12],["KHGX",29.47,-95.08],
  ["KHNX",36.31,-119.63],["KHPX",36.74,-87.28],["KHTX",34.93,-86.08],
  ["KHWA",38.51,-82.97],["KICT",37.65,-97.44],["KICX",37.59,-112.86],
  ["KILN",39.42,-83.82],["KILX",40.15,-89.34],["KIND",39.71,-86.28],
  ["KINX",36.18,-95.56],["KIWA",33.29,-111.67],["KIWX",41.36,-85.70],
  ["KJAX",30.48,-81.70],["KJGX",32.68,-83.35],["KJKL",37.59,-83.31],
  ["KKEY",24.55,-81.78],["KLBB",33.65,-101.81],["KLCH",30.13,-93.22],
  ["KLIX",30.34,-89.83],["KLNX",41.96,-100.58],["KLOT",41.60,-88.08],
  ["KLRX",40.74,-116.80],["KLSX",38.70,-90.68],["KLTX",33.99,-78.43],
  ["KLVX",37.98,-85.94],["KLWX",38.98,-77.48],["KLZK",34.84,-92.26],
  ["KMAF",31.94,-102.19],["KMAX",42.08,-122.72],["KMBX",48.39,-100.86],
  ["KMHX",34.78,-76.88],["KMKX",42.97,-88.55],["KMLB",28.11,-80.65],
  ["KMOB",30.68,-88.24],["KMPX",44.85,-93.57],["KMQT",46.53,-87.55],
  ["KMRX",36.17,-83.40],["KMSX",47.04,-113.99],["KMTX",41.26,-112.45],
  ["KMUX",37.16,-121.90],["KMVX",47.53,-97.33],["KMXX",32.54,-85.79],
  ["KNKX",32.92,-117.04],["KNQA",35.34,-89.87],["KOAX",41.32,-96.37],
  ["KOHX",36.25,-86.56],["KOKX",40.87,-72.86],["KOTX",47.68,-117.63],
  ["KPAH",37.07,-88.77],["KPBZ",40.53,-80.22],["KPDT",45.69,-118.85],
  ["KPOE",34.41,-116.16],["KPUX",38.46,-104.18],["KRAX",35.67,-78.49],
  ["KRGX",39.75,-119.46],["KRIW",43.07,-108.48],["KRLX",38.31,-81.72],
  ["KRTX",45.71,-122.97],["KSFX",43.11,-112.69],["KSGF",37.24,-93.40],
  ["KSHV",32.45,-93.84],["KSJT",31.37,-100.49],["KSOX",33.82,-117.64],
  ["KSRX",35.29,-94.36],["KTBW",27.71,-82.40],["KTFX",47.46,-111.39],
  ["KTLH",30.40,-84.33],["KTLX",35.33,-97.28],["KTWX",38.99,-96.23],
  ["KTYX",43.76,-75.68],["KUDX",44.13,-102.83],["KUEX",40.32,-98.44],
  ["KVAX",30.89,-83.00],["KVBX",34.84,-120.40],["KVNX",36.74,-98.13],
  ["KVTX",34.41,-119.18],["KVWX",38.26,-87.72],["KYUX",32.50,-114.66],
];

function nearestNexradForVad(lat, lon) {
  let best = NEXRAD_SITES[0], bestD = Infinity;
  for (const s of NEXRAD_SITES) {
    const d = (s[1] - lat) ** 2 + (s[2] - lon) ** 2;
    if (d < bestD) { bestD = d; best = s; }
  }
  return best[0];  // Returns e.g. "KTLX"
}

const SOURCE_META = {
  obs: {
    label: "Observed Radiosonde",
    desc: "Real observed upper-air data from the Iowa Environmental Mesonet and University of Wyoming archives.",
  },
  rap: {
    label: "RAP Model Analysis",
    desc: "Rapid Refresh model analysis at any lat/lon point over CONUS via NCEI THREDDS. ~13 km resolution.",
  },
  bufkit: {
    label: "BUFKIT Forecast",
    desc: "Station-based forecast soundings from HRRR, RAP, NAM, GFS, and other models via Iowa State archive.",
  },
  psu: {
    label: "PSU BUFKIT (Latest)",
    desc: "Latest model run from Penn State's real-time BUFKIT feed. Supports RAP, HRRR, NAM, GFS, and more.",
  },
  acars: {
    label: "ACARS Aircraft Obs",
    desc: "ACARS/AMDAR aircraft observation profiles at major airports from the IEM archive.",
  },
};

export default function ControlPanel({
  stations,
  sources,
  models,
  psuModels,
  onSubmit,
  loading,
  initialLoading,
  onRetry,
  connectError,
  riskData,
  onRiskDataChange,
  showRisk,
  onToggleRisk,
  showHistory,
  onToggleHistory,
  showMap,
  onToggleMap,
  showTimeSeries,
  onToggleTimeSeries,
  showCompare,
  onToggleCompare,
  showVwp,
  onToggleVwp,
  showMeso,
  onToggleMeso,
  showEnsemble,
  onToggleEnsemble,
  selectedStation,
  onStationChange,
  onSourceChange,
  mapLatLon,
  onFeedbackClick,
  showFeedback: feedbackActive,
  urlParams,
  theme,
  onToggleTheme,
  colorblind,
  onToggleColorblind,
  onNavigateUpload,
}) {
  const [source, setSourceLocal] = useState(urlParams?.source || "obs");
  const [station, setStationLocal] = useState(urlParams?.station || "OUN");
  const [date, setDate] = useState(() => {
    // Convert YYYYMMDDHH back to datetime-local value for the input
    if (urlParams?.date && /^\d{10}$/.test(urlParams.date)) {
      const d = urlParams.date;
      return `${d.slice(0,4)}-${d.slice(4,6)}-${d.slice(6,8)}T${d.slice(8,10)}:00`;
    }
    return "";
  });
  const [lat, setLat] = useState(urlParams?.lat != null ? String(urlParams.lat) : "");
  const [lon, setLon] = useState(urlParams?.lon != null ? String(urlParams.lon) : "");
  const [model, setModel] = useState(urlParams?.model || "hrrr");
  const [fhour, setFhour] = useState(urlParams?.fhour != null ? String(urlParams.fhour) : "0");
  const [stationSearch, setStationSearch] = useState("");
  const [scanning, setScanning] = useState(false);
  const [sortMode, setSortMode] = useState("az");
  const [favorites, setFavorites] = useState(() => getFavorites());
  const [soundingHour, setSoundingHour] = useState("latest");
  const listRef = useRef(null);

  // Surface modification state
  const [sfcModEnabled, setSfcModEnabled] = useState(false);
  const [sfcModT, setSfcModT] = useState("");
  const [sfcModTd, setSfcModTd] = useState("");
  const [sfcModWspd, setSfcModWspd] = useState("");
  const [sfcModWdir, setSfcModWdir] = useState("");

  // Custom storm motion state
  const [smEnabled, setSmEnabled] = useState(false);
  const [smDirection, setSmDirection] = useState("");
  const [smSpeed, setSmSpeed] = useState("");

  // VAD Wind Profile overlay state
  const [vadEnabled, setVadEnabled] = useState(false);

  // Storm-relative hodograph state
  const [srHodoEnabled, setSrHodoEnabled] = useState(false);

  // Profile smoothing state
  const [smoothEnabled, setSmoothEnabled] = useState(false);
  const [smoothSigma, setSmoothSigma] = useState("3");

  // Boundary line state
  const [boundaryEnabled, setBoundaryEnabled] = useState(false);
  const [boundaryOrientation, setBoundaryOrientation] = useState("");

  // Station groups state
  const [stationGroups, setStationGroups] = useState(() => getStationGroups());
  const [activeGroup, setActiveGroup] = useState("");
  const [showGroupSave, setShowGroupSave] = useState(false);
  const [newGroupName, setNewGroupName] = useState("");

  // Map zoom state
  const [mapZoom, setMapZoom] = useState("1");

  // Sync source to parent
  const setSource = (src) => {
    setSourceLocal(src);
    if (onSourceChange) onSourceChange(src);
  };

  // Sync station to parent
  const setStation = (id) => {
    setStationLocal(id);
    if (onStationChange) onStationChange(id);
  };

  // Sync station from parent (map click)
  useEffect(() => {
    if (selectedStation && selectedStation !== station) {
      setStationLocal(selectedStation);
      const stn = stations.find((s) => s.id === selectedStation);
      if (stn) {
        setLat(String(stn.lat));
        setLon(String(stn.lon));
      }
    }
  }, [selectedStation]); // eslint-disable-line react-hooks/exhaustive-deps

  // Sync lat/lon from map click
  useEffect(() => {
    if (mapLatLon) {
      setLat(String(mapLatLon.lat));
      setLon(String(mapLatLon.lon));
    }
  }, [mapLatLon]);

  // Scroll selected station into view on mount
  useEffect(() => {
    if (listRef.current) {
      const active = listRef.current.querySelector(".cp-station-item.active");
      if (active) active.scrollIntoView({ block: "center" });
    }
  }, [stations]);

  const needsLatLon = source === "rap";
  const needsStation = source === "obs" || source === "bufkit" || source === "acars" || source === "psu";
  const needsModel = source === "bufkit" || source === "psu";

  // Station group filter
  const activeGroupStations = activeGroup
    ? (stationGroups.find((g) => g.name === activeGroup)?.stations || [])
    : null;

  const handleSaveGroup = () => {
    const name = newGroupName.trim();
    if (!name) return;
    // Save currently visible (filtered) station IDs
    const ids = filteredStations.map((s) => s.id);
    if (ids.length === 0) return;
    const updated = saveStationGroup(name, ids);
    setStationGroups(updated);
    setActiveGroup(name);
    setShowGroupSave(false);
    setNewGroupName("");
  };

  const handleDeleteGroup = (name) => {
    const updated = deleteStationGroup(name);
    setStationGroups(updated);
    if (activeGroup === name) setActiveGroup("");
  };

  // Build risk lookup from scan results
  const riskMap = {};
  if (riskData) {
    riskData.stations.forEach((s, i) => {
      riskMap[s.id] = { rank: i + 1, stp: s.stp, raw: s.raw, cape: s.cape, srh: s.srh, bwd: s.bwd };
    });
  }

  // Merge station data with risk scores
  const mergedStations = stations.map((s) => ({
    ...s,
    risk: riskMap[s.id] || null,
  }));

  // Sort based on selected mode
  const sortedStations = [...mergedStations].sort((a, b) => {
    switch (sortMode) {
      case "za":
        return b.id.localeCompare(a.id);
      case "favs": {
        const aFav = favorites.includes(a.id) ? 0 : 1;
        const bFav = favorites.includes(b.id) ? 0 : 1;
        if (aFav !== bFav) return aFav - bFav;
        return a.id.localeCompare(b.id);
      }
      case "risk-high":
        // Stations with risk first (highest STP first), then unscanned at end
        if (a.risk && !b.risk) return -1;
        if (!a.risk && b.risk) return 1;
        if (a.risk && b.risk) return b.risk.stp - a.risk.stp;
        return a.id.localeCompare(b.id);
      case "risk-low":
        if (a.risk && !b.risk) return -1;
        if (!a.risk && b.risk) return 1;
        if (a.risk && b.risk) return a.risk.stp - b.risk.stp;
        return a.id.localeCompare(b.id);
      case "az":
      default:
        return a.id.localeCompare(b.id);
    }
  });

  const filteredStations = sortedStations.filter(
    (s) => {
      const matchesSearch = s.id.toLowerCase().includes(stationSearch.toLowerCase()) ||
        s.name.toLowerCase().includes(stationSearch.toLowerCase());
      const matchesGroup = activeGroupStations ? activeGroupStations.includes(s.id) : true;
      return matchesSearch && matchesGroup;
    }
  );

  const handleRiskScan = async () => {
    setScanning(true);
    try {
      const dateParam = date ? date.replace(/[-T:]/g, "").slice(0, 10) : undefined;
      const data = await fetchRiskScan(dateParam);
      onRiskDataChange(data);
      setSortMode("risk-high");
      // Auto-open the map
      if (!showMap && onToggleMap) onToggleMap();
      // Auto-select the highest risk station
      if (data.stations.length > 0) {
        handleStationSelect(data.stations[0].id);
      }
    } catch (e) {
      console.error("Risk scan failed:", e);
    } finally {
      setScanning(false);
    }
  };

  const handleSubmit = (e) => {
    e.preventDefault();
    const params = { source };

    if (needsStation) params.station = station;
    if (needsLatLon) {
      params.lat = parseFloat(lat);
      params.lon = parseFloat(lon);
    }
    if (date) params.date = date.replace(/[-T:]/g, "").slice(0, 10);
    if (needsModel) {
      params.model = model;
      params.fhour = parseInt(fhour) || 0;
    }

    // Surface modification
    if (sfcModEnabled) {
      const mod = {};
      if (sfcModT !== "") mod.temperature = parseFloat(sfcModT);
      if (sfcModTd !== "") mod.dewpoint = parseFloat(sfcModTd);
      if (sfcModWspd !== "") mod.wind_speed = parseFloat(sfcModWspd);
      if (sfcModWdir !== "") mod.wind_direction = parseFloat(sfcModWdir);
      if (Object.keys(mod).length > 0) params.surfaceMod = mod;
    }

    // Custom storm motion
    if (smEnabled && smDirection !== "" && smSpeed !== "") {
      params.stormMotion = {
        direction: parseFloat(smDirection),
        speed: parseFloat(smSpeed),
      };
    }

    // VAD Wind Profile overlay
    if (vadEnabled) {
      // Determine nearest NEXRAD radar for the selected station or lat/lon
      let vadLat = parseFloat(lat), vadLon = parseFloat(lon);
      if ((!vadLat || !vadLon) && station) {
        const stn = stations.find((s) => s.id === station);
        if (stn) { vadLat = stn.lat; vadLon = stn.lon; }
      }
      if (vadLat && vadLon) {
        params.vad = nearestNexradForVad(vadLat, vadLon);
      }
    }

    // Storm-relative hodograph
    if (srHodoEnabled) {
      params.srHodograph = true;
    }

    // Profile smoothing
    if (smoothEnabled) {
      const sigma = parseFloat(smoothSigma);
      if (sigma > 0) params.smoothing = sigma;
    }

    // Boundary line
    if (boundaryEnabled && boundaryOrientation !== "") {
      const deg = parseFloat(boundaryOrientation);
      if (!isNaN(deg)) params.boundaryOrientation = deg;
    }

    // Map zoom
    const mz = parseFloat(mapZoom);
    if (mz > 1) params.mapZoom = mz;

    onSubmit(params);
  };

  const handleStationSelect = (id, autoFetch = false) => {
    setStation(id);
    setStationSearch("");
    const stn = stations.find((s) => s.id === id);
    if (stn) {
      setLat(String(stn.lat));
      setLon(String(stn.lon));
    }
    if (autoFetch && !loading) {
      const params = { source, station: id };
      if (date) params.date = date.replace(/[-T:]/g, "").slice(0, 10);
      if (source === "bufkit") {
        params.model = model;
        params.fhour = parseInt(fhour) || 0;
      }
      if (source === "rap") {
        const s = stations.find((st) => st.id === id);
        if (s) {
          params.lat = s.lat;
          params.lon = s.lon;
        }
      }
      onSubmit(params);
    }
  };

  if (initialLoading || connectError) {
    return (
      <aside className="control-panel">
        <div className="cp-brand">
          <svg width="28" height="28" viewBox="0 0 28 28" fill="none" xmlns="http://www.w3.org/2000/svg">
            <rect x="2" y="2" width="24" height="24" rx="4" fill="rgba(59,130,246,0.12)" stroke="#3b82f6" strokeWidth="1.5"/>
            <path d="M8 22 L10 18 L11 16 L12 13 L14 10 L16 8 L19 5" stroke="#ef4444" strokeWidth="1.8" strokeLinecap="round" strokeLinejoin="round" fill="none"/>
            <path d="M6 22 L8 19 L9 17 L9.5 15 L10 13 L10 11 L10.5 9 L11 6" stroke="#22c55e" strokeWidth="1.8" strokeLinecap="round" strokeLinejoin="round" fill="none"/>
            <line x1="21" y1="10" x2="21" y2="22" stroke="#60a5fa" strokeWidth="1.2" strokeLinecap="round"/>
            <line x1="21" y1="10" x2="24" y2="8" stroke="#60a5fa" strokeWidth="1.2" strokeLinecap="round"/>
            <line x1="21" y1="13" x2="24" y2="11" stroke="#60a5fa" strokeWidth="1.2" strokeLinecap="round"/>
          </svg>
          <div>
            <h1 className="cp-brand-title">Sounding Analysis</h1>
            <p className="cp-brand-sub">Atmospheric Profile Tool</p>
          </div>
        </div>
        <div className="cp-loading">
          {initialLoading ? (
            <>
              <Loader2 className="spin" size={20} />
              <span>Connecting to API...</span>
              <span style={{ fontSize: "0.85em", opacity: 0.7 }}>May take a moment, be patient</span>
            </>
          ) : (
            <>
              <span style={{ color: "var(--danger, #e74c3c)" }}>{connectError}</span>
              <button
                type="button"
                onClick={onRetry}
                style={{
                  marginTop: 8,
                  padding: "6px 16px",
                  cursor: "pointer",
                  borderRadius: 6,
                  border: "1px solid var(--border, #555)",
                  background: "var(--surface, #2a2a2a)",
                  color: "inherit",
                }}
              >
                Retry
              </button>
            </>
          )}
        </div>
      </aside>
    );
  }

  return (
    <aside className="control-panel">
      {/* Brand */}
      <div className="cp-brand">
        <svg width="28" height="28" viewBox="0 0 28 28" fill="none" xmlns="http://www.w3.org/2000/svg">
          <rect x="2" y="2" width="24" height="24" rx="4" fill="rgba(59,130,246,0.12)" stroke="#3b82f6" strokeWidth="1.5"/>
          <path d="M8 22 L10 18 L11 16 L12 13 L14 10 L16 8 L19 5" stroke="#ef4444" strokeWidth="1.8" strokeLinecap="round" strokeLinejoin="round" fill="none"/>
          <path d="M6 22 L8 19 L9 17 L9.5 15 L10 13 L10 11 L10.5 9 L11 6" stroke="#22c55e" strokeWidth="1.8" strokeLinecap="round" strokeLinejoin="round" fill="none"/>
          <line x1="21" y1="10" x2="21" y2="22" stroke="#60a5fa" strokeWidth="1.2" strokeLinecap="round"/>
          <line x1="21" y1="10" x2="24" y2="8" stroke="#60a5fa" strokeWidth="1.2" strokeLinecap="round"/>
          <line x1="21" y1="13" x2="24" y2="11" stroke="#60a5fa" strokeWidth="1.2" strokeLinecap="round"/>
        </svg>
        <div>
          <h1 className="cp-brand-title">Sounding Analysis</h1>
          <p className="cp-brand-sub">Atmospheric Profile Tool</p>
        </div>
      </div>

      <form onSubmit={handleSubmit} className="cp-form">
        {/* Source */}
        <div className="cp-section">
          <label className="cp-label">
            <Database size={14} />
            Data Source
          </label>
          <div className="cp-source-grid">
            {sources.map((s) => {
              const meta = SOURCE_META[s.id];
              return (
                <div key={s.id} className="cp-source-btn-wrap">
                  <button
                    type="button"
                    className={`cp-source-btn ${source === s.id ? "active" : ""}`}
                    onClick={() => setSource(s.id)}
                    title={meta ? `${meta.label}\n${meta.desc}` : s.id.toUpperCase()}
                  >
                    <span className="cp-source-id">{s.id.toUpperCase()}</span>
                  </button>
                </div>
              );
            })}
          </div>
        </div>

        {/* Station */}
        {needsStation && (
          <div className="cp-section">
            <label className="cp-label">
              <MapPin size={14} />
              Station
            </label>
            <div className="cp-station-picker">
              <button
                type="button"
                className="cp-risk-btn"
                onClick={handleRiskScan}
                disabled={scanning}
              >
                {scanning ? (
                  <>
                    <Loader2 size={14} className="spin" />
                    Scanning stations...
                  </>
                ) : (
                  <>
                    <Zap size={14} />
                    {riskData ? "Rescan Tornado Risk" : "Scan Tornado Risk"}
                  </>
                )}
              </button>
              {riskData && (
                <p className="cp-risk-hint">
                  Scanned {riskData.stations.length} stations at {riskData.date}
                </p>
              )}
              <div className="cp-station-toolbar">
                <div className="cp-input-wrap cp-search-flex">
                  <Search size={14} className="cp-input-icon" />
                  <input
                    type="text"
                    className="cp-input"
                    placeholder="Filter..."
                    value={stationSearch}
                    onChange={(e) => setStationSearch(e.target.value)}
                  />
                </div>
                <div className="cp-sort-wrap">
                  <ArrowUpDown size={12} className="cp-sort-icon" />
                  <select
                    className="cp-sort-select"
                    value={sortMode}
                    onChange={(e) => setSortMode(e.target.value)}
                  >
                    <option value="az">A → Z</option>
                    <option value="za">Z → A</option>
                    <option value="favs">★ Favs</option>
                    <option value="risk-high">Risk ↓</option>
                    <option value="risk-low">Risk ↑</option>
                  </select>
                </div>
              </div>

              {/* Station Groups */}
              <div className="cp-groups-row">
                <select
                  className="cp-group-select"
                  value={activeGroup}
                  onChange={(e) => setActiveGroup(e.target.value)}
                >
                  <option value="">All Stations</option>
                  {stationGroups.map((g) => (
                    <option key={g.name} value={g.name}>{g.name} ({g.stations.length})</option>
                  ))}
                </select>
                {activeGroup && (
                  <button type="button" className="cp-group-del" onClick={() => handleDeleteGroup(activeGroup)} title="Delete group">
                    <X size={11} />
                  </button>
                )}
                <button type="button" className="cp-group-save-btn" onClick={() => setShowGroupSave((v) => !v)} title="Save current filter as group">
                  +
                </button>
              </div>
              {showGroupSave && (
                <div className="cp-group-save-row">
                  <input
                    type="text"
                    className="cp-input cp-input-sm"
                    placeholder="Group name..."
                    value={newGroupName}
                    onChange={(e) => setNewGroupName(e.target.value)}
                    onKeyDown={(e) => { if (e.key === "Enter") { e.preventDefault(); handleSaveGroup(); } }}
                    autoFocus
                  />
                  <button type="button" className="cp-group-confirm" onClick={handleSaveGroup}>Save</button>
                </div>
              )}

              <div className="cp-station-list" ref={listRef}>
                {filteredStations.length === 0 && (
                  <div className="cp-station-empty">No stations found</div>
                )}
                {filteredStations.map((s) => (
                  <button
                    key={s.id}
                    type="button"
                    className={`cp-station-item ${station === s.id ? "active" : ""}`}
                    onClick={() => handleStationSelect(s.id, !!riskData)}
                  >
                    <span
                      className={`cp-fav-star ${favorites.includes(s.id) ? "faved" : ""}`}
                      onClick={(e) => {
                        e.stopPropagation();
                        const { favorites: newFavs } = toggleFavorite(s.id);
                        setFavorites(newFavs);
                      }}
                      title={favorites.includes(s.id) ? "Remove from favorites" : "Add to favorites"}
                    >
                      <Star size={12} fill={favorites.includes(s.id) ? "currentColor" : "none"} />
                    </span>
                    <span className="cp-station-item-id">{s.id}</span>
                    <span className="cp-station-item-name">{s.name}</span>
                    {s.risk ? (
                      <span className={`cp-risk-score ${s.risk.stp >= 1 ? "high" : s.risk.stp >= 0.3 ? "med" : "low"}`}>
                        {s.risk.stp.toFixed(1)}
                      </span>
                    ) : (
                      <span className="cp-station-item-coords">
                        {s.lat.toFixed(1)}, {s.lon.toFixed(1)}
                      </span>
                    )}
                  </button>
                ))}
              </div>
            </div>
          </div>
        )}

        {/* Lat/Lon */}
        {needsLatLon && (
          <div className="cp-section">
            <label className="cp-label">
              <MapPin size={14} />
              Coordinates
            </label>
            <div className="cp-row">
              <div className="cp-field">
                <span className="cp-field-label">Lat</span>
                <input
                  type="number"
                  step="0.01"
                  className="cp-input cp-input-sm"
                  placeholder="35.22"
                  value={lat}
                  onChange={(e) => setLat(e.target.value)}
                  required
                />
              </div>
              <div className="cp-field">
                <span className="cp-field-label">Lon</span>
                <input
                  type="number"
                  step="0.01"
                  className="cp-input cp-input-sm"
                  placeholder="-97.46"
                  value={lon}
                  onChange={(e) => setLon(e.target.value)}
                  required
                />
              </div>
            </div>
          </div>
        )}

        {/* BUFKIT Model */}
        {needsModel && (
          <div className="cp-section">
            <label className="cp-label">
              <Layers size={14} />
              Model
            </label>
            <div className="cp-select-wrap">
              <select
                className="cp-select"
                value={model}
                onChange={(e) => setModel(e.target.value)}
              >
                {(source === "psu" ? (psuModels || []) : models).map((m) => (
                  <option key={m.id} value={m.id}>
                    {m.id.toUpperCase()} — {m.name}
                  </option>
                ))}
              </select>
              <ChevronDown size={14} className="cp-select-icon" />
            </div>
            <div className="cp-field" style={{ marginTop: 8 }}>
              <span className="cp-field-label">Forecast Hour</span>
              <input
                type="number"
                min="0"
                max="384"
                className="cp-input cp-input-sm"
                placeholder="0"
                value={fhour}
                onChange={(e) => setFhour(e.target.value)}
              />
            </div>
          </div>
        )}

        {/* Date */}
        <div className="cp-section">
          <label className="cp-label">
            <Calendar size={14} />
            Date / Time (UTC)
          </label>
          {source === "obs" || source === "acars" ? (
            <>
              <div className="cp-date-row">
                <input
                  type="date"
                  className="cp-input cp-calendar-input"
                  value={date ? date.slice(0, 10) : ""}
                  onChange={(e) => {
                    const d = e.target.value;
                    if (d) {
                      const hour = soundingHour === "12" ? "12:00" : "00:00";
                      setDate(`${d}T${hour}`);
                    } else {
                      setDate("");
                    }
                  }}
                />
                {date && (
                  <button
                    type="button"
                    className="cp-date-clear"
                    onClick={() => setDate("")}
                    title="Clear date"
                  >
                    ✕
                  </button>
                )}
              </div>
              <div className="cp-sounding-hour-row">
                <button
                  type="button"
                  className={`cp-hour-btn ${soundingHour === "latest" ? "active" : ""}`}
                  onClick={() => {
                    setSoundingHour("latest");
                    setDate("");
                  }}
                >
                  Latest
                </button>
                <button
                  type="button"
                  className={`cp-hour-btn ${soundingHour === "00" ? "active" : ""}`}
                  onClick={() => {
                    setSoundingHour("00");
                    if (date) setDate(`${date.slice(0, 10)}T00:00`);
                  }}
                >
                  00Z
                </button>
                <button
                  type="button"
                  className={`cp-hour-btn ${soundingHour === "12" ? "active" : ""}`}
                  onClick={() => {
                    setSoundingHour("12");
                    if (date) setDate(`${date.slice(0, 10)}T12:00`);
                  }}
                >
                  12Z
                </button>
              </div>
              <p className="cp-hint">
                Latest fetches the most recent available sounding
              </p>
            </>
          ) : (
            <>
              <div className="cp-date-row">
                <input
                  type="datetime-local"
                  className="cp-input cp-calendar-input"
                  value={date}
                  onChange={(e) => setDate(e.target.value)}
                />
                {date && (
                  <button
                    type="button"
                    className="cp-date-clear"
                    onClick={() => setDate("")}
                    title="Clear date"
                  >
                    ✕
                  </button>
                )}
              </div>
              <p className="cp-hint">
                Leave blank to use the most recent time
              </p>
            </>
          )}
        </div>

        {/* ── Modifications ── */}
        <div className="cp-section-group">
        <span className="cp-group-label">Modifications</span>

        {/* Surface Modification */}
        <div className={`cp-accordion ${sfcModEnabled ? "cp-accordion--active" : ""}`}>
          <button
            type="button"
            className="cp-accordion-header"
            onClick={() => setSfcModEnabled((v) => !v)}
          >
            <div className="cp-accordion-left">
              <Thermometer size={14} className="cp-accordion-icon" />
              <span className="cp-accordion-title">Surface Modification</span>
            </div>
            <div className="cp-accordion-right">
              <span className={`cp-toggle-chip ${sfcModEnabled ? "on" : ""}`}>
                {sfcModEnabled ? "ON" : "OFF"}
              </span>
              <ChevronRight size={14} className={`cp-accordion-chevron ${sfcModEnabled ? "open" : ""}`} />
            </div>
          </button>
          <div className={`cp-accordion-body ${sfcModEnabled ? "expanded" : ""}`}>
            <div className="cp-accordion-content">
              <div className="cp-input-row">
                <div className="cp-input-group">
                  <label className="cp-input-group-label">Temperature</label>
                  <div className="cp-input-with-unit">
                    <input type="number" step="0.1" className="cp-input cp-input-sm" placeholder="—" value={sfcModT} onChange={(e) => setSfcModT(e.target.value)} />
                    <span className="cp-unit-badge">°C</span>
                  </div>
                </div>
                <div className="cp-input-group">
                  <label className="cp-input-group-label">Dewpoint</label>
                  <div className="cp-input-with-unit">
                    <input type="number" step="0.1" className="cp-input cp-input-sm" placeholder="—" value={sfcModTd} onChange={(e) => setSfcModTd(e.target.value)} />
                    <span className="cp-unit-badge">°C</span>
                  </div>
                </div>
              </div>
              <div className="cp-input-row">
                <div className="cp-input-group">
                  <label className="cp-input-group-label">Wind Speed</label>
                  <div className="cp-input-with-unit">
                    <input type="number" step="1" className="cp-input cp-input-sm" placeholder="—" value={sfcModWspd} onChange={(e) => setSfcModWspd(e.target.value)} />
                    <span className="cp-unit-badge">kt</span>
                  </div>
                </div>
                <div className="cp-input-group">
                  <label className="cp-input-group-label">Wind Dir</label>
                  <div className="cp-input-with-unit">
                    <input type="number" step="1" min="0" max="360" className="cp-input cp-input-sm" placeholder="—" value={sfcModWdir} onChange={(e) => setSfcModWdir(e.target.value)} />
                    <span className="cp-unit-badge">°</span>
                  </div>
                </div>
              </div>
              <button
                type="button"
                className="cp-reset-btn"
                onClick={() => { setSfcModT(""); setSfcModTd(""); setSfcModWspd(""); setSfcModWdir(""); }}
              >
                <RotateCcw size={11} /> Reset values
              </button>
            </div>
          </div>
        </div>

        {/* Custom Storm Motion */}
        <div className={`cp-accordion ${smEnabled ? "cp-accordion--active" : ""}`}>
          <button
            type="button"
            className="cp-accordion-header"
            onClick={() => setSmEnabled((v) => !v)}
          >
            <div className="cp-accordion-left">
              <Wind size={14} className="cp-accordion-icon" />
              <span className="cp-accordion-title">Custom Storm Motion</span>
            </div>
            <div className="cp-accordion-right">
              <span className={`cp-toggle-chip ${smEnabled ? "on" : ""}`}>
                {smEnabled ? "ON" : "OFF"}
              </span>
              <ChevronRight size={14} className={`cp-accordion-chevron ${smEnabled ? "open" : ""}`} />
            </div>
          </button>
          <div className={`cp-accordion-body ${smEnabled ? "expanded" : ""}`}>
            <div className="cp-accordion-content">
              <div className="cp-input-row">
                <div className="cp-input-group">
                  <label className="cp-input-group-label">Direction</label>
                  <div className="cp-input-with-unit">
                    <input type="number" step="1" min="0" max="360" className="cp-input cp-input-sm" placeholder="—" value={smDirection} onChange={(e) => setSmDirection(e.target.value)} />
                    <span className="cp-unit-badge">°</span>
                  </div>
                </div>
                <div className="cp-input-group">
                  <label className="cp-input-group-label">Speed</label>
                  <div className="cp-input-with-unit">
                    <input type="number" step="1" className="cp-input cp-input-sm" placeholder="—" value={smSpeed} onChange={(e) => setSmSpeed(e.target.value)} />
                    <span className="cp-unit-badge">kt</span>
                  </div>
                </div>
              </div>
              <button
                type="button"
                className="cp-reset-btn"
                onClick={() => { setSmDirection(""); setSmSpeed(""); }}
              >
                <RotateCcw size={11} /> Reset values
              </button>
            </div>
          </div>
        </div>

        {/* VAD Wind Profile Overlay */}
        <button
          type="button"
          className={`cp-toggle-btn ${vadEnabled ? "cp-toggle-btn--active" : ""}`}
          onClick={() => setVadEnabled((v) => !v)}
          title={`Overlay NEXRAD VAD winds on the hodograph from the nearest WSR-88D radar${(() => {
            let vLat = parseFloat(lat), vLon = parseFloat(lon);
            if ((!vLat || !vLon) && station) {
              const stn = stations.find((s) => s.id === station);
              if (stn) { vLat = stn.lat; vLon = stn.lon; }
            }
            if (vLat && vLon) return ` (${nearestNexradForVad(vLat, vLon)})`;
            return "";
          })()}`}
        >
          <Layers size={14} />
          <span>VAD Wind Profile</span>
          <span className={`cp-toggle-chip ${vadEnabled ? "on" : ""}`}>
            {vadEnabled ? "ON" : "OFF"}
          </span>
        </button>

        {/* Storm-Relative Hodograph */}
        <button
          type="button"
          className={`cp-toggle-btn ${srHodoEnabled ? "cp-toggle-btn--active" : ""}`}
          onClick={() => setSrHodoEnabled((v) => !v)}
          title="Plot hodograph in storm-relative frame — subtracts Bunkers RM (or custom SM) from all winds so storm motion is at the origin"
        >
          <Crosshair size={14} />
          <span>SR Hodograph</span>
          <span className={`cp-toggle-chip ${srHodoEnabled ? "on" : ""}`}>
            {srHodoEnabled ? "ON" : "OFF"}
          </span>
        </button>

        {/* Boundary Line */}
        <div className="cp-accordion-section">
          <button
            type="button"
            className={`cp-toggle-btn ${boundaryEnabled ? "cp-toggle-btn--active" : ""}`}
            onClick={() => setBoundaryEnabled((v) => !v)}
            title="Draw a boundary orientation line on the hodograph — represents an outflow boundary, front, or dryline orientation"
          >
            <Minus size={14} />
            <span>Boundary Line</span>
            <span className={`cp-toggle-chip ${boundaryEnabled ? "on" : ""}`}>
              {boundaryEnabled ? "ON" : "OFF"}
            </span>
          </button>
          {boundaryEnabled && (
            <div style={{ padding: "6px 10px" }}>
              <div className="cp-input-row">
                <div className="cp-input-group">
                  <label className="cp-input-group-label">Orientation</label>
                  <div className="cp-input-with-unit">
                    <input
                      type="number"
                      step="5"
                      min="0"
                      max="360"
                      className="cp-input cp-input-sm"
                      placeholder="e.g. 210"
                      value={boundaryOrientation}
                      onChange={(e) => setBoundaryOrientation(e.target.value)}
                    />
                    <span className="cp-unit-badge">°</span>
                  </div>
                </div>
              </div>
              <p style={{ margin: "4px 0 0", fontSize: 10, color: "var(--fg-faint, #707070)" }}>
                Direction the boundary runs (0–360°). e.g. 210° = SW–NE oriented.
              </p>
            </div>
          )}
        </div>

        {/* Profile Smoothing */}
        <div className="cp-accordion-section">
          <button
            type="button"
            className={`cp-toggle-btn ${smoothEnabled ? "cp-toggle-btn--active" : ""}`}
            onClick={() => setSmoothEnabled((v) => !v)}
            title="Apply Gaussian smoothing to T, Td, and wind profiles — reduces noise in ACARS and model data"
          >
            <Waves size={14} />
            <span>Profile Smoothing</span>
            <span className={`cp-toggle-chip ${smoothEnabled ? "on" : ""}`}>
              {smoothEnabled ? "ON" : "OFF"}
            </span>
          </button>
          {smoothEnabled && (
            <div style={{ padding: "6px 10px" }}>
              <label className="cp-field-label" style={{ marginBottom: 0, display: "flex", alignItems: "center", gap: 8 }}>
                <span style={{ minWidth: 50 }}>σ = {smoothSigma}</span>
                <input
                  type="range"
                  min="1"
                  max="10"
                  step="0.5"
                  value={smoothSigma}
                  onChange={(e) => setSmoothSigma(e.target.value)}
                  style={{ flex: 1 }}
                  title="Gaussian sigma in data levels (higher = smoother). Typical: 2-5"
                />
              </label>
            </div>
          )}
        </div>

        {/* Map Inset Zoom */}
        <div className="cp-accordion-section">
          <div style={{ padding: "6px 10px" }}>
            <label className="cp-field-label" style={{ marginBottom: 0, display: "flex", alignItems: "center", gap: 8 }}>
              <Map size={13} />
              <span style={{ minWidth: 70 }}>Map {mapZoom}x</span>
              <input
                type="range"
                min="1"
                max="8"
                step="0.5"
                value={mapZoom}
                onChange={(e) => setMapZoom(e.target.value)}
                style={{ flex: 1 }}
                title="Zoom level for the CONUS mini-map inset on the sounding plot (1x = full CONUS)"
              />
            </label>
          </div>
        </div>
        </div>{/* end Modifications */}

        {/* Submit */}
        <button
          type="submit"
          className="cp-submit"
          disabled={loading}
        >
          {loading ? (
            <>
              <Loader2 size={16} className="spin" />
              Fetching & Analyzing...
            </>
          ) : (
            <>
              <Search size={16} />
              Generate Sounding
            </>
          )}
        </button>
      </form>

      {/* ── Bottom group: toggles + settings + footer ── */}
      <div className="cp-bottom">
        <div className="cp-section-group">
          <span className="cp-group-label">Tools</span>
          <div className="cp-tools-grid">
            <button
              type="button"
              className={`cp-tool-btn ${showMap ? "active" : ""}`}
              onClick={onToggleMap}
            >
              <Map size={13} />
              {showMap ? "Hide" : "Map"}
            </button>
            {riskData && (
              <button
                type="button"
                className={`cp-tool-btn ${showRisk ? "active" : ""}`}
                onClick={onToggleRisk}
              >
                <Zap size={13} />
                {showRisk ? "Hide" : "Risk"}
              </button>
            )}
            <button
              type="button"
              className={`cp-tool-btn ${showTimeSeries ? "active" : ""}`}
              onClick={onToggleTimeSeries}
            >
              <TrendingUp size={13} />
              Trends
            </button>
            <button
              type="button"
              className={`cp-tool-btn ${showCompare ? "active" : ""}`}
              onClick={onToggleCompare}
            >
              <GitCompareArrows size={13} />
              Compare
            </button>
            <button
              type="button"
              className={`cp-tool-btn ${showVwp ? "active" : ""}`}
              onClick={onToggleVwp}
            >
              <Radio size={13} />
              VWP
            </button>
            <button
              type="button"
              className={`cp-tool-btn ${showHistory ? "active" : ""}`}
              onClick={onToggleHistory}
            >
              <History size={13} />
              History
            </button>
            <button
              type="button"
              className={`cp-tool-btn ${showMeso ? "active" : ""}`}
              onClick={onToggleMeso}
            >
              <Layers size={13} />
              Meso
            </button>
            <button
              type="button"
              className={`cp-tool-btn ${showEnsemble ? "active" : ""}`}
              onClick={onToggleEnsemble}
            >
              <Layers size={13} />
              Ensemble
            </button>
          </div>
        </div>

        <div className="cp-settings-row">
          <button
            type="button"
            className="cp-settings-btn"
            onClick={onToggleTheme}
            title={theme === "dark" ? "Switch to light theme" : "Switch to dark theme"}
          >
            {theme === "dark" ? <Sun size={14} /> : <Moon size={14} />}
            <span>{theme === "dark" ? "Light" : "Dark"}</span>
          </button>
          <button
            type="button"
            className={`cp-settings-btn ${colorblind ? "active" : ""}`}
            onClick={onToggleColorblind}
            title="Toggle color-blind safe palette"
          >
            <Eye size={14} />
            <span>CB Mode</span>
          </button>
          <button
            type="button"
            className="cp-settings-btn"
            onClick={onNavigateUpload}
            title="Upload custom sounding data"
          >
            <Upload size={14} />
            <span>Upload</span>
          </button>
        </div>

        <div className="cp-footer">
          <button
            type="button"
            className={`cp-footer-btn ${feedbackActive ? "active" : ""}`}
            onClick={onFeedbackClick}
            title="Send feedback"
          >
            <MessageSquarePlus size={14} />
            <span>Feedback</span>
          </button>
          <a
            href="https://github.com/ShianMike/SoundingAnalysis"
            target="_blank"
            rel="noopener noreferrer"
            className="cp-footer-btn"
            title="View on GitHub"
          >
            <Github size={14} />
            <span>GitHub</span>
          </a>
        </div>
      </div>
    </aside>
  );
}
