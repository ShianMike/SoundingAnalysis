import { useState, useRef, useEffect, useCallback } from "react";
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
  Keyboard,
  Waves,
  Radio,
  Sun,
  Moon,
  Eye,
  Upload,
  Minus,
  X,
  GripVertical,
  ChevronUp,
  PanelLeftClose,
  Sliders,
  Wrench,
  AlertTriangle,
} from "lucide-react";
import { fetchRiskScan, fetchForecastRiskScan } from "../api";
import { getFavorites, toggleFavorite } from "../favorites";
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

/* Model metadata: short label, resolution hint, max forecast hour, step */
const MODEL_META = {
  hrrr:    { short: "HRRR",     res: "3 km · hourly",     maxF: 48,  step: 1 },
  rap:     { short: "RAP",      res: "13 km · hourly",    maxF: 21,  step: 1 },
  nam:     { short: "NAM",      res: "12 km · hourly",    maxF: 84,  step: 1 },
  namnest: { short: "NAM Nest", res: "3 km · hourly",     maxF: 60,  step: 1 },
  nam4km:  { short: "NAM 4km",  res: "4 km · hourly",     maxF: 60,  step: 1 },
  gfs:     { short: "GFS",      res: "global · 3-hourly", maxF: 384, step: 3 },
  sref:    { short: "SREF",     res: "ens mean · 3-hr",   maxF: 87,  step: 3 },
  hiresw:  { short: "HiResW",   res: "NMMB / ARW",        maxF: 48,  step: 1 },
};

/* Common quick-pick forecast hours */
const FHOUR_PRESETS = [0, 3, 6, 12, 18, 24, 36, 48];

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
  onNavigateEnsemble,
  selectedStation,
  onStationChange,
  onSourceChange,
  mapLatLon,
  mapStormMotion,
  mapBoundaryOrientation,
  onFeedbackClick,
  showFeedback: feedbackActive,
  urlParams,
  theme,
  onToggleTheme,
  colorblind,
  onToggleColorblind,
  onNavigateUpload,
  onShowShortcuts,
}) {
  /* ── Drag & collapse state for fullscreen-float mode ─────── */
  const [floatCollapsed, setFloatCollapsed] = useState(false);
  const [dragPos, setDragPos] = useState({ x: 16, y: 16 });
  const dragRef = useRef(null);
  const isDragging = useRef(false);
  const dragStart = useRef({ x: 0, y: 0, posX: 0, posY: 0 });

  const onPointerDown = useCallback((e) => {
    if (!document.body.classList.contains("smap-fullscreen-active")) return;
    // Don't start drag if clicking a button (let the click event fire)
    if (e.target.closest("button")) return;
    isDragging.current = true;
    dragStart.current = { x: e.clientX, y: e.clientY, posX: dragPos.x, posY: dragPos.y };
    e.currentTarget.setPointerCapture(e.pointerId);
  }, [dragPos]);

  /** Minimum Y = bottom of the map toolbar + gap so sidebar never overlaps it */
  const getMinY = useCallback(() => {
    const tb = document.querySelector('.smap-fullscreen .smap-toolbar');
    if (tb) {
      const rect = tb.getBoundingClientRect();
      return rect.bottom + 8;        // 8px gap below toolbar
    }
    return 90;                         // safe fallback
  }, []);

  const onPointerMove = useCallback((e) => {
    if (!isDragging.current) return;
    const dx = e.clientX - dragStart.current.x;
    const dy = e.clientY - dragStart.current.y;
    const minY = getMinY();
    const pad = 12;
    const panelW = dragRef.current ? dragRef.current.offsetWidth : 270;
    setDragPos({
      x: Math.max(pad, Math.min(window.innerWidth - panelW - pad, dragStart.current.posX + dx)),
      y: Math.max(minY, Math.min(window.innerHeight - 40, dragStart.current.posY + dy)),
    });
  }, [getMinY]);

  const onPointerUp = useCallback(() => {
    isDragging.current = false;
  }, []);

  // Centre sidebar when entering fullscreen; reset when exiting
  // Only react when the fullscreen state actually *changes* (not every class mutation)
  const wasFullscreen = useRef(false);
  useEffect(() => {
    const obs = new MutationObserver(() => {
      const active = document.body.classList.contains("smap-fullscreen-active");
      if (active === wasFullscreen.current) return; // no change — ignore
      wasFullscreen.current = active;
      if (active) {
        // Position at absolute top-left corner
        setDragPos({ x: 0, y: 0 });
        setFloatCollapsed(true);
      } else {
        setDragPos({ x: 16, y: 16 });
        setFloatCollapsed(false);
      }
    });
    obs.observe(document.body, { attributes: true, attributeFilter: ["class"] });
    return () => obs.disconnect();
  }, []);

  // When a toolbar dropdown opens/closes, push sidebar below the toolbar
  useEffect(() => {
    const tb = document.querySelector('.smap-fullscreen .smap-toolbar');
    if (!tb) return;
    const ro = new ResizeObserver(() => {
      if (!document.body.classList.contains("smap-fullscreen-active")) return;
      const bottom = tb.getBoundingClientRect().bottom;
      setDragPos((prev) => ({ x: prev.x, y: (bottom > 0 ? bottom : 0) + 8 }));
    });
    ro.observe(tb);
    return () => ro.disconnect();
  });

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
  const [rapStationSearch, setRapStationSearch] = useState("");
  const [rapDropOpen, setRapDropOpen] = useState(false);
  const [scanning, setScanning] = useState(false);
  const [scanMode, setScanMode] = useState("obs"); // "obs" | "forecast"
  const [fcstModel, setFcstModel] = useState("hrrr");
  const [fcstFhour, setFcstFhour] = useState("12");
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

  // Navigation tab state
  const [navTab, setNavTab] = useState("data");




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

  // Sync station from parent (map click) and scroll into view
  useEffect(() => {
    if (selectedStation && selectedStation !== station) {
      setStationLocal(selectedStation);
      const stn = stations.find((s) => s.id === selectedStation);
      if (stn) {
        setLat(String(stn.lat));
        setLon(String(stn.lon));
      }
    }
    // Scroll the selected station into view in the list
    if (selectedStation && listRef.current) {
      requestAnimationFrame(() => {
        const el = listRef.current?.querySelector(`.cp-station-item[data-id="${selectedStation}"]`);
        if (el) el.scrollIntoView({ block: "center", behavior: "smooth" });
      });
    }
  }, [selectedStation]); // eslint-disable-line react-hooks/exhaustive-deps

  // Sync lat/lon from map click
  useEffect(() => {
    if (mapLatLon) {
      setLat(String(mapLatLon.lat));
      setLon(String(mapLatLon.lon));
    }
  }, [mapLatLon]);

  // Sync custom storm motion from map draw tool
  useEffect(() => {
    if (
      mapStormMotion &&
      mapStormMotion.direction != null &&
      mapStormMotion.speed != null
    ) {
      setSmEnabled(true);
      setSmDirection(String(Math.round(mapStormMotion.direction)));
      setSmSpeed(String(Math.round(mapStormMotion.speed)));
    }
  }, [mapStormMotion]);

  // Sync boundary orientation from map draw tool
  useEffect(() => {
    if (mapBoundaryOrientation != null) {
      setBoundaryEnabled(true);
      setBoundaryOrientation(String(Math.round(mapBoundaryOrientation)));
    }
  }, [mapBoundaryOrientation]);

  // Close RAP station dropdown on outside click
  useEffect(() => {
    if (!rapDropOpen) return;
    const handler = (e) => {
      if (!e.target.closest(".cp-rap-station-pick")) setRapDropOpen(false);
    };
    document.addEventListener("mousedown", handler);
    return () => document.removeEventListener("mousedown", handler);
  }, [rapDropOpen]);

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
      return matchesSearch;
    }
  );

  const handleRiskScan = async () => {
    setScanning(true);
    try {
      let data;
      if (scanMode === "forecast") {
        data = await fetchForecastRiskScan({
          model: fcstModel,
          fhour: parseInt(fcstFhour) || 0,
        });
      } else {
        const dateParam = date ? date.replace(/[-T:]/g, "").slice(0, 10) : undefined;
        data = await fetchRiskScan(dateParam);
      }
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
        <div className="cp-loading-state">
          {initialLoading ? (
            <>
              <Loader2 className="spin cp-loading-icon" size={32} />
              <span className="cp-loading-title">Connecting to API&hellip;</span>
              <span className="cp-loading-hint">May take a moment, be patient</span>
              <div className="cp-loading-skeleton">
                <div className="cp-skel-bar" style={{ width: "80%" }} />
                <div className="cp-skel-bar" style={{ width: "60%" }} />
                <div className="cp-skel-bar" style={{ width: "90%" }} />
                <div className="cp-skel-bar" style={{ width: "45%" }} />
              </div>
            </>
          ) : (
            <>
              <AlertTriangle size={32} style={{ color: "var(--danger, #e74c3c)" }} />
              <span className="cp-loading-title" style={{ color: "var(--danger, #e74c3c)" }}>Connection Failed</span>
              <span className="cp-loading-hint">{connectError}</span>
              <button type="button" className="cp-retry-btn" onClick={onRetry}>Retry</button>
            </>
          )}
        </div>
      </aside>
    );
  }

  const isFloat = typeof window !== "undefined" && document.body.classList.contains("smap-fullscreen-active");

  return (
    <aside
      className={`control-panel${floatCollapsed ? " cp-float-collapsed" : ""}`}
      style={isFloat ? { left: dragPos.x, top: dragPos.y, bottom: 'auto', margin: 0 } : undefined}
      ref={dragRef}
    >
      {/* Drag handle – only visible in float mode */}
      <div
        className="cp-drag-handle"
        onPointerDown={onPointerDown}
        onPointerMove={onPointerMove}
        onPointerUp={onPointerUp}
      >
        <GripVertical size={14} />
        <span className="cp-drag-label">Sounding Analysis</span>
        <div className="cp-drag-actions">
          <button
            type="button"
            className="cp-drag-btn"
            onClick={() => setFloatCollapsed((v) => !v)}
            title={floatCollapsed ? "Expand panel" : "Collapse panel"}
          >
            {floatCollapsed ? <ChevronDown size={14} /> : <ChevronUp size={14} />}
          </button>
        </div>
      </div>

      {/* Brand + Tab Navigation — compact single header */}
      <div className="cp-brand">
        <span className="cp-brand-title">Sounding Analysis</span>
        <div className="cp-brand-actions">
          <button
            type="button"
            className="cp-brand-icon-btn"
            onClick={onToggleTheme}
            title={theme === "dark" ? "Switch to light theme" : "Switch to dark theme"}
          >
            {theme === "dark" ? <Sun size={14} /> : <Moon size={14} />}
          </button>
        </div>
      </div>

      {/* Tab Navigation */}
      <nav className="cp-nav">
        <button
          type="button"
          className={`cp-nav-tab${navTab === "data" ? " active" : ""}`}
          onClick={() => setNavTab("data")}
        >
          <Database size={14} />
          <span>Data</span>
        </button>
        <button
          type="button"
          className={`cp-nav-tab${navTab === "modify" ? " active" : ""}`}
          onClick={() => setNavTab("modify")}
        >
          <Sliders size={14} />
          <span>Modify</span>
        </button>
        <button
          type="button"
          className={`cp-nav-tab${navTab === "tools" ? " active" : ""}`}
          onClick={() => setNavTab("tools")}
        >
          <Wrench size={14} />
          <span>Tools</span>
        </button>
      </nav>

      <div className="cp-tab-content">

      {/* ═══ DATA TAB ═══ */}
      {navTab === "data" && (
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
              {/* Scan mode toggle */}
              <div className="cp-scan-tabs">
                <button
                  type="button"
                  className={`cp-scan-tab${scanMode === "obs" ? " active" : ""}`}
                  onClick={() => setScanMode("obs")}
                >
                  <Eye size={12} /> Observed
                </button>
                <button
                  type="button"
                  className={`cp-scan-tab${scanMode === "forecast" ? " active" : ""}`}
                  onClick={() => setScanMode("forecast")}
                >
                  <TrendingUp size={12} /> Forecast
                </button>
              </div>

              {/* Forecast model + fhour pickers */}
              {scanMode === "forecast" && (
                <div className="cp-scan-forecast-opts">
                  <select
                    className="cp-input cp-scan-select"
                    value={fcstModel}
                    onChange={(e) => {
                      setFcstModel(e.target.value);
                      const meta = MODEL_META[e.target.value] || { maxF: 384, step: 1 };
                      if (parseInt(fcstFhour) > meta.maxF) setFcstFhour(String(meta.maxF));
                    }}
                  >
                    {["hrrr", "rap", "nam", "namnest", "gfs"].map((m) => (
                      <option key={m} value={m}>{MODEL_META[m]?.short || m}</option>
                    ))}
                  </select>
                  <div className="cp-scan-fhour">
                    <label className="cp-scan-fhour-label">F{fcstFhour}</label>
                    <input
                      type="range"
                      className="cp-scan-slider"
                      min={0}
                      max={(MODEL_META[fcstModel] || { maxF: 48 }).maxF}
                      step={(MODEL_META[fcstModel] || { step: 1 }).step}
                      value={fcstFhour}
                      onChange={(e) => setFcstFhour(e.target.value)}
                    />
                  </div>
                </div>
              )}

              <button
                type="button"
                className="cp-risk-btn"
                onClick={handleRiskScan}
                disabled={scanning}
              >
                {scanning ? (
                  <>
                    <Loader2 size={14} className="spin" />
                    Analyzing...
                  </>
                ) : (
                  <>
                    <Search size={14} />
                    {riskData ? "Rescan Severe Risk" : "Scan Severe Risk"}
                  </>
                )}
              </button>
              {riskData && (
                <p className="cp-risk-hint">
                  {riskData.model
                    ? `${(MODEL_META[riskData.model]?.short || riskData.model).toUpperCase()} F${riskData.fhour} · ${riskData.stations.length} stations`
                    : `Scanned ${riskData.stations.length} stations at ${riskData.date}`}
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

              <div className="cp-station-list" ref={listRef}>
                {filteredStations.length === 0 && (
                  <div className="cp-station-empty">No stations found</div>
                )}
                {filteredStations.map((s) => (
                  <button
                    key={s.id}
                    type="button"
                    data-id={s.id}
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

            {/* Quick-pick from station */}
            <div className="cp-rap-station-pick">
              <div className="cp-input-wrap cp-search-flex">
                <Search size={14} className="cp-input-icon" />
                <input
                  type="text"
                  className="cp-input"
                  placeholder="Jump to station..."
                  value={rapStationSearch}
                  onChange={(e) => setRapStationSearch(e.target.value)}
                  onFocus={() => setRapDropOpen(true)}
                />
              </div>
              {rapDropOpen && rapStationSearch.trim().length > 0 && (
                <div className="cp-rap-dropdown">
                  {stations
                    .filter((s) =>
                      s.id.toLowerCase().includes(rapStationSearch.toLowerCase()) ||
                      s.name.toLowerCase().includes(rapStationSearch.toLowerCase())
                    )
                    .slice(0, 8)
                    .map((s) => (
                      <button
                        key={s.id}
                        type="button"
                        className="cp-rap-drop-item"
                        onClick={() => {
                          setLat(String(s.lat));
                          setLon(String(s.lon));
                          setRapStationSearch("");
                          setRapDropOpen(false);
                        }}
                      >
                        <span className="cp-rap-drop-id">{s.id}</span>
                        <span className="cp-rap-drop-name">{s.name}</span>
                        <span className="cp-rap-drop-coords">{s.lat.toFixed(2)}, {s.lon.toFixed(2)}</span>
                      </button>
                    ))}
                  {stations.filter((s) =>
                    s.id.toLowerCase().includes(rapStationSearch.toLowerCase()) ||
                    s.name.toLowerCase().includes(rapStationSearch.toLowerCase())
                  ).length === 0 && (
                    <div className="cp-rap-drop-empty">No stations found</div>
                  )}
                </div>
              )}
            </div>

            <p className="cp-rap-hint">
              <Crosshair size={12} />
              Type a station or click anywhere on the map
            </p>

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
            {lat && lon && (
              <div className="cp-rap-coords-preview">
                <MapPin size={12} />
                <span>{parseFloat(lat).toFixed(2)}°N, {Math.abs(parseFloat(lon)).toFixed(2)}°{parseFloat(lon) < 0 ? "W" : "E"}</span>
              </div>
            )}
          </div>
        )}

        {/* BUFKIT Model */}
        {needsModel && (() => {
          const modelList = source === "psu" ? (psuModels || []) : models;
          const meta = MODEL_META[model] || { maxF: 384, step: 1 };
          const maxF = meta.maxF;
          const step = meta.step;
          const currentFhour = Math.min(parseInt(fhour) || 0, maxF);
          return (
          <div className="cp-section">
            <label className="cp-label">
              <Layers size={14} />
              Model
            </label>
            <div className="cp-model-grid">
              {modelList.map((m) => {
                const mm = MODEL_META[m.id];
                return (
                  <button
                    key={m.id}
                    type="button"
                    className={`cp-model-btn ${model === m.id ? "active" : ""}`}
                    onClick={() => {
                      setModel(m.id);
                      // Clamp forecast hour to new model's max
                      const newMeta = MODEL_META[m.id] || { maxF: 384, step: 1 };
                      const cur = parseInt(fhour) || 0;
                      if (cur > newMeta.maxF) setFhour(String(newMeta.maxF));
                    }}
                    title={m.name}
                  >
                    <span className="cp-model-id">{mm?.short || m.id.toUpperCase()}</span>
                    {mm && <span className="cp-model-res">{mm.res}</span>}
                  </button>
                );
              })}
            </div>

            {/* Forecast hour */}
            <div className="cp-fhour-section">
              <div className="cp-fhour-header">
                <Clock size={13} />
                <span>Forecast Hour</span>
                <span className="cp-fhour-value">F{String(currentFhour).padStart(2, "0")}</span>
              </div>
              <input
                type="range"
                min="0"
                max={maxF}
                step={step}
                className="cp-fhour-slider"
                value={currentFhour}
                onChange={(e) => setFhour(e.target.value)}
              />
              <div className="cp-fhour-range">
                <span>F00</span>
                <span>F{String(maxF).padStart(2, "0")}</span>
              </div>
              <div className="cp-fhour-presets">
                {FHOUR_PRESETS.filter((h) => h <= maxF).map((h) => (
                  <button
                    key={h}
                    type="button"
                    className={`cp-fhour-preset ${currentFhour === h ? "active" : ""}`}
                    onClick={() => setFhour(String(h))}
                  >
                    {h === 0 ? "Anl" : `F${String(h).padStart(2, "0")}`}
                  </button>
                ))}
              </div>
            </div>
          </div>
          );
        })()}

        {/* Date */}
        <div className="cp-section">
          <label className="cp-label">
            <Calendar size={14} />
            Date / Time (UTC)
          </label>
          {source === "obs" || source === "acars" ? (
            <div className="cp-dt-picker">
              <div className="cp-dt-hours">
                {[
                  { id: "latest", label: "Latest" },
                  { id: "00", label: "00Z" },
                  { id: "12", label: "12Z" },
                ].map((h) => (
                  <button
                    key={h.id}
                    type="button"
                    className={`cp-dt-hour${soundingHour === h.id ? " active" : ""}`}
                    onClick={() => {
                      setSoundingHour(h.id);
                      if (h.id === "latest") {
                        setDate("");
                      } else if (date) {
                        setDate(`${date.slice(0, 10)}T${h.id}:00`);
                      }
                    }}
                  >
                    {soundingHour === h.id && <Clock size={11} />}
                    {h.label}
                  </button>
                ))}
              </div>
              <div className="cp-dt-date-wrap">
                <Calendar size={13} className="cp-dt-date-icon" />
                <input
                  type="date"
                  className="cp-dt-date-input"
                  value={date ? date.slice(0, 10) : ""}
                  onChange={(e) => {
                    const d = e.target.value;
                    if (d) {
                      const hour = soundingHour === "12" ? "12:00" : "00:00";
                      if (soundingHour === "latest") setSoundingHour("00");
                      setDate(`${d}T${hour}`);
                    } else {
                      setDate("");
                    }
                  }}
                />
                {date && (
                  <button type="button" className="cp-dt-clear" onClick={() => { setDate(""); setSoundingHour("latest"); }} title="Reset to latest">
                    <X size={12} />
                  </button>
                )}
              </div>
              {soundingHour === "latest" && (
                <p className="cp-hint" style={{ margin: 0 }}>Auto-selects the most recent available sounding</p>
              )}
            </div>
          ) : (
            <div className="cp-dt-picker">
              <div className="cp-dt-date-wrap">
                <Clock size={13} className="cp-dt-date-icon" />
                <input
                  type="datetime-local"
                  className="cp-dt-date-input"
                  value={date}
                  onChange={(e) => setDate(e.target.value)}
                />
                {date && (
                  <button type="button" className="cp-dt-clear" onClick={() => setDate("")} title="Clear date">
                    <X size={12} />
                  </button>
                )}
              </div>
              <p className="cp-hint" style={{ margin: 0 }}>Leave blank for the most recent time</p>
            </div>
          )}
        </div>

        {/* Submit */}
        <button
          type="submit"
          className="cp-submit"
          disabled={loading}
        >
          {loading ? (
            <>
              <Loader2 size={16} className="spin" />
              Analyzing...
            </>
          ) : (
            <>
              <Search size={16} />
              Generate Sounding
            </>
          )}
        </button>
      </form>
      )}

      {/* ═══ MODIFY TAB ═══ */}
      {navTab === "modify" && (
      <div className="cp-form">
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

        </div>{/* end Modifications */}

        <p className="cp-hint" style={{ marginTop: 4 }}>
          These options modify the next sounding fetch. Go to Data tab and click Generate to apply.
        </p>
      </div>
      )}

      {/* ═══ TOOLS TAB ═══ */}
      {navTab === "tools" && (
      <div className="cp-form">
        <div className="cp-section-group">
          <span className="cp-group-label">Views</span>
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
              className="cp-tool-btn"
              onClick={onNavigateEnsemble}
            >
              <Layers size={13} />
              Plume
            </button>
          </div>
        </div>

        <div className="cp-section-group">
          <span className="cp-group-label">Settings</span>
          <div className="cp-settings-grid">
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
        </div>

      </div>
      )}

      </div>{/* end cp-tab-content */}

      {/* Links — always visible at bottom */}
      <div className="cp-section-group cp-bottom-links">
        <span className="cp-group-label">Links</span>
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
          {onShowShortcuts && (
            <button
              type="button"
              className="cp-footer-btn"
              onClick={onShowShortcuts}
              title="Keyboard shortcuts (?)"
            >
              <Keyboard size={14} />
              <span>Shortcuts</span>
            </button>
          )}
        </div>
      </div>
    </aside>
  );
}
