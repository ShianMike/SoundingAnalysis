/**
 * Centralised configuration constants for the frontend.
 *
 * Extracted from per-component literals so the data lives in one place and can
 * later be ported to TypeScript with shared types.
 */

/* ── NEXRAD WSR-88D sites (CONUS) ────────────────────────────────
 * Stored as 3-letter site codes (no "K" prefix). Helpers below add
 * the prefix when callers need the full ICAO identifier (e.g. "KTLX"). */
export const NEXRAD_SITES = [
  ["ABR",45.46,-98.41],["ABX",35.15,-106.82],["AKQ",36.98,-77.01],
  ["AMA",35.23,-101.71],["AMX",25.61,-80.41],["APX",44.91,-84.72],
  ["ARX",43.82,-91.19],["ATX",48.19,-122.50],["BBX",39.50,-121.63],
  ["BGM",42.20,-75.98],["BMX",33.17,-86.77],["BOX",41.96,-71.14],
  ["BRO",25.92,-97.42],["BUF",42.95,-78.74],["BYX",24.60,-81.70],
  ["CAE",33.95,-81.12],["CBW",46.04,-67.81],["CBX",43.49,-116.24],
  ["CCX",40.92,-78.00],["CLE",41.41,-81.86],["CLX",32.66,-81.04],
  ["CRP",27.78,-97.51],["CXX",44.51,-73.17],["CYS",41.15,-104.81],
  ["DAX",38.50,-121.68],["DDC",37.76,-99.97],["DFX",29.27,-100.28],
  ["DGX",32.28,-89.98],["DIX",39.95,-74.41],["DLH",46.84,-92.21],
  ["DMX",41.73,-93.72],["DOX",38.83,-75.44],["DTX",42.70,-83.47],
  ["DVN",41.61,-90.58],["DYX",32.54,-99.25],["EAX",38.81,-94.26],
  ["EMX",31.89,-110.63],["ENX",42.59,-74.06],["EOX",31.46,-85.46],
  ["EPZ",31.87,-106.70],["ESX",35.70,-114.89],["EVX",30.56,-85.92],
  ["EWX",29.70,-98.03],["EYX",35.10,-117.56],["FCX",37.02,-80.27],
  ["FDR",34.36,-98.98],["FDX",34.64,-103.63],["FFC",33.36,-84.57],
  ["FSD",43.59,-96.73],["FSX",34.57,-111.20],["FTG",39.79,-104.55],
  ["FWS",32.57,-97.30],["GGW",48.21,-106.63],["GJX",39.06,-108.21],
  ["GLD",39.37,-101.70],["GRB",44.50,-88.11],["GRK",30.72,-97.38],
  ["GRR",42.89,-85.54],["GSP",34.88,-82.22],["GWX",33.90,-88.33],
  ["GYX",43.89,-70.26],["HDX",33.08,-106.12],["HGX",29.47,-95.08],
  ["HNX",36.31,-119.63],["HPX",36.74,-87.28],["HTX",34.93,-86.08],
  ["HWA",38.51,-82.97],["ICT",37.65,-97.44],["ICX",37.59,-112.86],
  ["ILN",39.42,-83.82],["ILX",40.15,-89.34],["IND",39.71,-86.28],
  ["INX",36.18,-95.56],["IWA",33.29,-111.67],["IWX",41.36,-85.70],
  ["JAX",30.48,-81.70],["JGX",32.68,-83.35],["JKL",37.59,-83.31],
  ["KEY",24.55,-81.78],["KSG",31.48,-82.31],["LBB",33.65,-101.81],
  ["LCH",30.13,-93.22],["LIX",30.34,-89.83],["LNX",41.96,-100.58],
  ["LOT",41.60,-88.08],["LRX",40.74,-116.80],["LSX",38.70,-90.68],
  ["LTX",33.99,-78.43],["LVX",37.98,-85.94],["LWX",38.98,-77.48],
  ["LZK",34.84,-92.26],["MAF",31.94,-102.19],["MAX",42.08,-122.72],
  ["MBX",48.39,-100.86],["MHX",34.78,-76.88],["MKX",42.97,-88.55],
  ["MLB",28.11,-80.65],["MOB",30.68,-88.24],["MPX",44.85,-93.57],
  ["MQT",46.53,-87.55],["MRX",36.17,-83.40],["MSX",47.04,-113.99],
  ["MTX",41.26,-112.45],["MUX",37.16,-121.90],["MVX",47.53,-97.33],
  ["MXX",32.54,-85.79],["NKX",32.92,-117.04],["NQA",35.34,-89.87],
  ["OAX",41.32,-96.37],["OHX",36.25,-86.56],["OKX",40.87,-72.86],
  ["OTX",47.68,-117.63],["PAH",37.07,-88.77],["PBZ",40.53,-80.22],
  ["PDT",45.69,-118.85],["POE",34.41,-116.16],["PUX",38.46,-104.18],
  ["RAX",35.67,-78.49],["RGX",39.75,-119.46],["RIW",43.07,-108.48],
  ["RLX",38.31,-81.72],["RTX",45.71,-122.97],["SFX",43.11,-112.69],
  ["SGF",37.24,-93.40],["SHV",32.45,-93.84],["SJT",31.37,-100.49],
  ["SOX",33.82,-117.64],["SRX",35.29,-94.36],["TBW",27.71,-82.40],
  ["TFX",47.46,-111.39],["TLH",30.40,-84.33],["TLX",35.33,-97.28],
  ["TWX",38.99,-96.23],["TYX",43.76,-75.68],["UDX",44.13,-102.83],
  ["UEX",40.32,-98.44],["VAX",30.89,-83.00],["VBX",34.84,-120.40],
  ["VNX",36.74,-98.13],["VTX",34.41,-119.18],["VWX",38.26,-87.72],
  ["YUX",32.50,-114.66],
];

/**
 * Find the nearest WSR-88D site to the given lat/lon.
 * Returns `{ id, idK, lat, lon }` where `id` is the 3-letter code (e.g. "TLX")
 * and `idK` is the 4-letter ICAO form (e.g. "KTLX"). Returns `null` if the
 * input is not finite.
 */
export function nearestNexrad(lat, lon) {
  if (!Number.isFinite(lat) || !Number.isFinite(lon)) return null;
  let best = NEXRAD_SITES[0];
  let bestD = Infinity;
  for (const s of NEXRAD_SITES) {
    const d = (s[1] - lat) ** 2 + (s[2] - lon) ** 2;
    if (d < bestD) { bestD = d; best = s; }
  }
  return { id: best[0], idK: `K${best[0]}`, lat: best[1], lon: best[2] };
}

/* ── Sounding source metadata ────────────────────────────────── */
export const SOURCE_META = {
  obs: {
    label: "Observed Radiosonde",
    desc: "Real observed upper-air data from the Iowa Environmental Mesonet and University of Wyoming archives.",
  },
  bufkit: {
    label: "BUFKIT Forecast",
    desc: "Station-based forecast soundings from HRRR, RAP, NAM, GFS, and other models via Iowa State archive.",
  },
  psu: {
    label: "PSU BUFKIT (Latest)",
    desc: "Latest model run from Penn State's real-time BUFKIT feed. Supports RAP, HRRR, NAM, GFS, and more.",
  },
};

/* ── Forecast model metadata: short label, resolution hint, max forecast hour, step ── */
export const MODEL_META = {
  hrrr:    { short: "HRRR",     res: "3 km · hourly",     maxF: 48,  step: 1, lag: 3, interval: 1 },
  rap:     { short: "RAP",      res: "13 km · hourly",    maxF: 21,  step: 1, lag: 3, interval: 1 },
  nam:     { short: "NAM",      res: "12 km · hourly",    maxF: 84,  step: 1, lag: 5, interval: 6 },
  namnest: { short: "NAM Nest", res: "3 km · hourly",     maxF: 60,  step: 1, lag: 5, interval: 6 },
  nam4km:  { short: "NAM 4km",  res: "4 km · hourly",     maxF: 60,  step: 1, lag: 5, interval: 6 },
  gfs:     { short: "GFS",      res: "global · 3-hourly", maxF: 384, step: 3, lag: 6, interval: 6 },
  sref:    { short: "SREF",     res: "ens mean · 3-hr",   maxF: 87,  step: 3, lag: 6, interval: 6 },
  hiresw:  { short: "HiResW",   res: "NMMB / ARW",        maxF: 48,  step: 1, lag: 5, interval: 6 },
};

/* ── Common quick-pick forecast hours ────────────────────────── */
export const FHOUR_PRESETS = [0, 3, 6, 12, 18, 24, 36, 48];
