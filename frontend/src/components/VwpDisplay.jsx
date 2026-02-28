import { useState, useEffect, useCallback } from "react";
import { X, Loader2, Radio, RefreshCw } from "lucide-react";
import { fetchVwpDisplay } from "../api";
import "./VwpDisplay.css";

// NEXRAD sites for nearest radar lookup
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
  ["KICT",37.65,-97.44],["KICX",37.59,-112.86],
  ["KILN",39.42,-83.82],["KILX",40.15,-89.34],["KIND",39.71,-86.28],
  ["KINX",36.18,-95.56],["KIWA",33.29,-111.67],["KIWX",41.36,-85.70],
  ["KJAX",30.48,-81.70],["KJGX",32.68,-83.35],["KJKL",37.59,-83.31],
  ["KLBB",33.65,-101.81],["KLCH",30.13,-93.22],
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
  ["KPUX",38.46,-104.18],["KRAX",35.67,-78.49],
  ["KRGX",39.75,-119.46],["KRIW",43.07,-108.48],["KRLX",38.31,-81.72],
  ["KRTX",45.71,-122.97],["KSFX",43.11,-112.69],["KSGF",37.24,-93.40],
  ["KSHV",32.45,-93.84],["KSJT",31.37,-100.49],["KSOX",33.82,-117.64],
  ["KSRX",35.29,-94.36],["KTBW",27.71,-82.40],["KTFX",47.46,-111.39],
  ["KTLH",30.40,-84.33],["KTLX",35.33,-97.28],["KTWX",38.99,-96.23],
  ["KTYX",43.76,-75.68],["KUDX",44.13,-102.83],["KUEX",40.32,-98.44],
  ["KVAX",30.89,-83.00],["KVBX",34.84,-120.40],["KVNX",36.74,-98.13],
  ["KVTX",34.41,-119.18],["KVWX",38.26,-87.72],["KYUX",32.50,-114.66],
];

function nearestNexrad(lat, lon) {
  let best = NEXRAD_SITES[0], bestD = Infinity;
  for (const s of NEXRAD_SITES) {
    const d = (s[1] - lat) ** 2 + (s[2] - lon) ** 2;
    if (d < bestD) { bestD = d; best = s; }
  }
  return best[0];
}

export default function VwpDisplay({ stations, selectedStation, onClose }) {
  const [radar, setRadar] = useState("");
  const [hours, setHours] = useState(12);
  const [loading, setLoading] = useState(false);
  const [imageData, setImageData] = useState(null);
  const [error, setError] = useState(null);
  const [meta, setMeta] = useState(null);

  // Auto-set radar based on selected station
  useEffect(() => {
    if (selectedStation && stations?.length) {
      const st = stations.find((s) => s.id === selectedStation);
      if (st) {
        setRadar(nearestNexrad(st.lat, st.lon));
      }
    }
  }, [selectedStation, stations]);

  const handleFetch = useCallback(async () => {
    if (!radar) return;
    setLoading(true);
    setError(null);
    try {
      const data = await fetchVwpDisplay(radar, hours);
      setImageData(data.image);
      setMeta({ snapshots: data.snapshots, timeRange: data.timeRange, radar: data.radar });
    } catch (e) {
      setError(e.message);
      setImageData(null);
    } finally {
      setLoading(false);
    }
  }, [radar, hours]);

  // Auto-fetch on mount if radar is set
  useEffect(() => {
    if (radar) handleFetch();
  }, []); // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <div className="vwp-panel">
      <div className="vwp-header">
        <h3>
          <Radio size={14} />
          VWP Time-Height Display
        </h3>
        <div style={{ display: "flex", alignItems: "center", gap: 10 }}>
          {meta && (
            <span className="vwp-header-meta">
              {meta.snapshots} scans | {meta.timeRange?.start} → {meta.timeRange?.end}
            </span>
          )}
          <button className="vwp-close-btn" onClick={onClose} title="Close VWP">
            <X size={16} />
          </button>
        </div>
      </div>

      <div className="vwp-controls">
        <span className="vwp-radar-label">Radar:</span>
        <input
          type="text"
          className="vwp-radar-input"
          value={radar}
          onChange={(e) => setRadar(e.target.value.toUpperCase())}
          placeholder="KTLX"
          maxLength={4}
        />
        <span className="vwp-hours-label">Hours:</span>
        <select
          className="vwp-hours-select"
          value={hours}
          onChange={(e) => setHours(Number(e.target.value))}
        >
          <option value={3}>3h</option>
          <option value={6}>6h</option>
          <option value={12}>12h</option>
          <option value={24}>24h</option>
        </select>
        <button
          className="vwp-fetch-btn"
          onClick={handleFetch}
          disabled={loading || !radar}
        >
          {loading ? <Loader2 size={12} className="spin" /> : <RefreshCw size={12} />}
          {loading ? "Loading..." : "Fetch"}
        </button>
      </div>

      <div className="vwp-body">
        {loading && (
          <div className="vwp-loading">
            <Loader2 size={20} className="spin" />
            Fetching VWP data from {radar}... This may take a moment.
          </div>
        )}
        {error && !loading && (
          <div className="vwp-error">{error}</div>
        )}
        {imageData && !loading && (
          <div className="vwp-image-wrap">
            <img
              src={`data:image/png;base64,${imageData}`}
              alt={`VWP for ${meta?.radar || radar}`}
            />
          </div>
        )}
        {!imageData && !loading && !error && (
          <div className="vwp-empty">
            Enter a NEXRAD radar ID and click Fetch to generate VWP display.
          </div>
        )}
      </div>
    </div>
  );
}
