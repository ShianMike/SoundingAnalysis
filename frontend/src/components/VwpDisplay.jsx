import { useState, useEffect, useCallback, useMemo, useRef } from "react";
import { X, Loader2, Radio, RefreshCw, Search, ChevronDown } from "lucide-react";
import { fetchVwpDisplay } from "../api";
import "./VwpDisplay.css";

// NEXRAD sites: [id, lat, lon, location]
const NEXRAD_SITES = [
  ["KABR",45.46,-98.41,"Aberdeen, SD"],["KABX",35.15,-106.82,"Albuquerque, NM"],
  ["KAKQ",36.98,-77.01,"Wakefield, VA"],["KAMA",35.23,-101.71,"Amarillo, TX"],
  ["KAMX",25.61,-80.41,"Miami, FL"],["KAPX",44.91,-84.72,"Gaylord, MI"],
  ["KARX",43.82,-91.19,"La Crosse, WI"],["KATX",48.19,-122.50,"Seattle, WA"],
  ["KBBX",39.50,-121.63,"Beale AFB, CA"],["KBGM",42.20,-75.98,"Binghamton, NY"],
  ["KBMX",33.17,-86.77,"Birmingham, AL"],["KBOX",41.96,-71.14,"Boston, MA"],
  ["KBRO",25.92,-97.42,"Brownsville, TX"],["KBUF",42.95,-78.74,"Buffalo, NY"],
  ["KBYX",24.60,-81.70,"Key West, FL"],["KCAE",33.95,-81.12,"Columbia, SC"],
  ["KCBW",46.04,-67.81,"Caribou, ME"],["KCBX",43.49,-116.24,"Boise, ID"],
  ["KCCX",40.92,-78.00,"State College, PA"],["KCLE",41.41,-81.86,"Cleveland, OH"],
  ["KCLX",32.66,-81.04,"Charleston, SC"],["KCRP",27.78,-97.51,"Corpus Christi, TX"],
  ["KCXX",44.51,-73.17,"Burlington, VT"],["KCYS",41.15,-104.81,"Cheyenne, WY"],
  ["KDAX",38.50,-121.68,"Sacramento, CA"],["KDDC",37.76,-99.97,"Dodge City, KS"],
  ["KDFX",29.27,-100.28,"Laughlin AFB, TX"],["KDGX",32.28,-89.98,"Brandon, MS"],
  ["KDIX",39.95,-74.41,"Mt Holly, NJ"],["KDLH",46.84,-92.21,"Duluth, MN"],
  ["KDMX",41.73,-93.72,"Des Moines, IA"],["KDOX",38.83,-75.44,"Dover, DE"],
  ["KDTX",42.70,-83.47,"Detroit, MI"],["KDVN",41.61,-90.58,"Davenport, IA"],
  ["KDYX",32.54,-99.25,"Dyess AFB, TX"],["KEAX",38.81,-94.26,"Kansas City, MO"],
  ["KEMX",31.89,-110.63,"Tucson, AZ"],["KENX",42.59,-74.06,"Albany, NY"],
  ["KEOX",31.46,-85.46,"Fort Rucker, AL"],["KEPZ",31.87,-106.70,"El Paso, TX"],
  ["KESX",35.70,-114.89,"Las Vegas, NV"],["KEVX",30.56,-85.92,"Eglin AFB, FL"],
  ["KEWX",29.70,-98.03,"Austin/San Antonio, TX"],["KEYX",35.10,-117.56,"Edwards AFB, CA"],
  ["KFCX",37.02,-80.27,"Roanoke, VA"],["KFDR",34.36,-98.98,"Frederick, OK"],
  ["KFDX",34.64,-103.63,"Cannon AFB, NM"],["KFFC",33.36,-84.57,"Atlanta, GA"],
  ["KFSD",43.59,-96.73,"Sioux Falls, SD"],["KFSX",34.57,-111.20,"Flagstaff, AZ"],
  ["KFTG",39.79,-104.55,"Denver, CO"],["KFWS",32.57,-97.30,"Dallas/Fort Worth, TX"],
  ["KGGW",48.21,-106.63,"Glasgow, MT"],["KGJX",39.06,-108.21,"Grand Junction, CO"],
  ["KGLD",39.37,-101.70,"Goodland, KS"],["KGRB",44.50,-88.11,"Green Bay, WI"],
  ["KGRK",30.72,-97.38,"Central Texas, TX"],["KGRR",42.89,-85.54,"Grand Rapids, MI"],
  ["KGSP",34.88,-82.22,"Greenville, SC"],["KGWX",33.90,-88.33,"Columbus AFB, MS"],
  ["KGYX",43.89,-70.26,"Portland, ME"],["KHDX",33.08,-106.12,"Holloman AFB, NM"],
  ["KHGX",29.47,-95.08,"Houston, TX"],["KHNX",36.31,-119.63,"Hanford, CA"],
  ["KHPX",36.74,-87.28,"Fort Campbell, KY"],["KHTX",34.93,-86.08,"Huntsville, AL"],
  ["KICT",37.65,-97.44,"Wichita, KS"],["KICX",37.59,-112.86,"Cedar City, UT"],
  ["KILN",39.42,-83.82,"Wilmington, OH"],["KILX",40.15,-89.34,"Lincoln, IL"],
  ["KIND",39.71,-86.28,"Indianapolis, IN"],["KINX",36.18,-95.56,"Tulsa, OK"],
  ["KIWA",33.29,-111.67,"Phoenix, AZ"],["KIWX",41.36,-85.70,"Fort Wayne, IN"],
  ["KJAX",30.48,-81.70,"Jacksonville, FL"],["KJGX",32.68,-83.35,"Robins AFB, GA"],
  ["KJKL",37.59,-83.31,"Jackson, KY"],["KLBB",33.65,-101.81,"Lubbock, TX"],
  ["KLCH",30.13,-93.22,"Lake Charles, LA"],["KLIX",30.34,-89.83,"New Orleans, LA"],
  ["KLNX",41.96,-100.58,"North Platte, NE"],["KLOT",41.60,-88.08,"Chicago, IL"],
  ["KLRX",40.74,-116.80,"Elko, NV"],["KLSX",38.70,-90.68,"St. Louis, MO"],
  ["KLTX",33.99,-78.43,"Wilmington, NC"],["KLVX",37.98,-85.94,"Louisville, KY"],
  ["KLWX",38.98,-77.48,"Washington, DC"],["KLZK",34.84,-92.26,"Little Rock, AR"],
  ["KMAF",31.94,-102.19,"Midland/Odessa, TX"],["KMAX",42.08,-122.72,"Medford, OR"],
  ["KMBX",48.39,-100.86,"Minot, ND"],["KMHX",34.78,-76.88,"Morehead City, NC"],
  ["KMKX",42.97,-88.55,"Milwaukee, WI"],["KMLB",28.11,-80.65,"Melbourne, FL"],
  ["KMOB",30.68,-88.24,"Mobile, AL"],["KMPX",44.85,-93.57,"Minneapolis, MN"],
  ["KMQT",46.53,-87.55,"Marquette, MI"],["KMRX",36.17,-83.40,"Knoxville, TN"],
  ["KMSX",47.04,-113.99,"Missoula, MT"],["KMTX",41.26,-112.45,"Salt Lake City, UT"],
  ["KMUX",37.16,-121.90,"San Francisco, CA"],["KMVX",47.53,-97.33,"Fargo, ND"],
  ["KMXX",32.54,-85.79,"Maxwell AFB, AL"],["KNKX",32.92,-117.04,"San Diego, CA"],
  ["KNQA",35.34,-89.87,"Memphis, TN"],["KOAX",41.32,-96.37,"Omaha, NE"],
  ["KOHX",36.25,-86.56,"Nashville, TN"],["KOKX",40.87,-72.86,"New York, NY"],
  ["KOTX",47.68,-117.63,"Spokane, WA"],["KPAH",37.07,-88.77,"Paducah, KY"],
  ["KPBZ",40.53,-80.22,"Pittsburgh, PA"],["KPDT",45.69,-118.85,"Pendleton, OR"],
  ["KPUX",38.46,-104.18,"Pueblo, CO"],["KRAX",35.67,-78.49,"Raleigh, NC"],
  ["KRGX",39.75,-119.46,"Reno, NV"],["KRIW",43.07,-108.48,"Riverton, WY"],
  ["KRLX",38.31,-81.72,"Charleston, WV"],["KRTX",45.71,-122.97,"Portland, OR"],
  ["KSFX",43.11,-112.69,"Pocatello, ID"],["KSGF",37.24,-93.40,"Springfield, MO"],
  ["KSHV",32.45,-93.84,"Shreveport, LA"],["KSJT",31.37,-100.49,"San Angelo, TX"],
  ["KSOX",33.82,-117.64,"Santa Ana Mtns, CA"],["KSRX",35.29,-94.36,"Fort Smith, AR"],
  ["KTBW",27.71,-82.40,"Tampa Bay, FL"],["KTFX",47.46,-111.39,"Great Falls, MT"],
  ["KTLH",30.40,-84.33,"Tallahassee, FL"],["KTLX",35.33,-97.28,"Oklahoma City, OK"],
  ["KTWX",38.99,-96.23,"Topeka, KS"],["KTYX",43.76,-75.68,"Montague, NY"],
  ["KUDX",44.13,-102.83,"Rapid City, SD"],["KUEX",40.32,-98.44,"Hastings, NE"],
  ["KVAX",30.89,-83.00,"Valdosta, GA"],["KVBX",34.84,-120.40,"Vandenberg, CA"],
  ["KVNX",36.74,-98.13,"Vance AFB, OK"],["KVTX",34.41,-119.18,"Los Angeles, CA"],
  ["KVWX",38.26,-87.72,"Evansville, IN"],["KYUX",32.50,-114.66,"Yuma, AZ"],
];

function distSq(lat1, lon1, lat2, lon2) {
  return (lat1 - lat2) ** 2 + (lon1 - lon2) ** 2;
}

function nearestNexrad(lat, lon) {
  let best = NEXRAD_SITES[0], bestD = Infinity;
  for (const s of NEXRAD_SITES) {
    const d = distSq(s[1], s[2], lat, lon);
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
  const [radarSearch, setRadarSearch] = useState("");
  const [dropdownOpen, setDropdownOpen] = useState(false);
  const dropdownRef = useRef(null);

  // Close dropdown on outside click
  useEffect(() => {
    const handleClick = (e) => {
      if (dropdownRef.current && !dropdownRef.current.contains(e.target)) {
        setDropdownOpen(false);
      }
    };
    document.addEventListener("mousedown", handleClick);
    return () => document.removeEventListener("mousedown", handleClick);
  }, []);

  // Resolve selected station lat/lon for sorting
  const stationLatLon = useMemo(() => {
    if (selectedStation && stations?.length) {
      const st = stations.find((s) => s.id === selectedStation);
      if (st) return { lat: st.lat, lon: st.lon };
    }
    return null;
  }, [selectedStation, stations]);

  // Sorted & filtered radar list
  const radarOptions = useMemo(() => {
    let list = NEXRAD_SITES.map((s) => ({
      id: s[0], lat: s[1], lon: s[2], location: s[3],
      dist: stationLatLon ? distSq(s[1], s[2], stationLatLon.lat, stationLatLon.lon) : 0,
    }));
    // Sort by distance from selected station if available
    if (stationLatLon) {
      list.sort((a, b) => a.dist - b.dist);
    }
    // Filter by search
    if (radarSearch.trim()) {
      const q = radarSearch.trim().toUpperCase();
      list = list.filter((r) => r.id.includes(q) || r.location.toUpperCase().includes(q));
    }
    return list;
  }, [stationLatLon, radarSearch]);

  // Get display label for selected radar
  const selectedRadarInfo = useMemo(() => {
    const site = NEXRAD_SITES.find((s) => s[0] === radar);
    return site ? `${site[0]} – ${site[3]}` : radar || "Select radar...";
  }, [radar]);

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
        <div className="vwp-radar-dropdown" ref={dropdownRef}>
          <button
            className="vwp-radar-trigger"
            onClick={() => { setDropdownOpen((v) => !v); setRadarSearch(""); }}
            type="button"
          >
            <span className="vwp-radar-trigger-text">{selectedRadarInfo}</span>
            <ChevronDown size={12} className={dropdownOpen ? "vwp-chevron-open" : ""} />
          </button>
          {dropdownOpen && (
            <div className="vwp-radar-menu">
              <div className="vwp-radar-search-wrap">
                <Search size={12} />
                <input
                  type="text"
                  className="vwp-radar-search"
                  value={radarSearch}
                  onChange={(e) => setRadarSearch(e.target.value)}
                  placeholder="Search radar or city..."
                  autoFocus
                />
              </div>
              <div className="vwp-radar-list">
                {radarOptions.length === 0 && (
                  <div className="vwp-radar-no-match">No matching radars</div>
                )}
                {radarOptions.map((r) => (
                  <button
                    key={r.id}
                    className={`vwp-radar-option ${r.id === radar ? "active" : ""}`}
                    onClick={() => { setRadar(r.id); setDropdownOpen(false); setRadarSearch(""); }}
                    type="button"
                  >
                    <span className="vwp-radar-option-id">{r.id}</span>
                    <span className="vwp-radar-option-loc">{r.location}</span>
                    {stationLatLon && (
                      <span className="vwp-radar-option-dist">
                        {Math.round(Math.sqrt(r.dist) * 111)} km
                      </span>
                    )}
                  </button>
                ))}
              </div>
            </div>
          )}
        </div>
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
            Select a NEXRAD radar and click Fetch to generate VWP display.
          </div>
        )}
      </div>
    </div>
  );
}
