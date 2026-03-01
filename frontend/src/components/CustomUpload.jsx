import { useState, useRef } from "react";
import { ArrowLeft, Upload, FileText, Loader2, Zap, MapPin } from "lucide-react";
import { fetchCustomSounding, fetchUploadFile } from "../api";
import "./CustomUpload.css";

const FORMATS = [
  { id: "auto", label: "Auto-detect" },
  { id: "csv", label: "CSV" },
  { id: "sharppy", label: "SHARPpy" },
  { id: "cm1", label: "CM1" },
  { id: "wrf", label: "WRF netCDF" },
];

const PLACEHOLDER = `# Paste sounding data here (SHARPpy / CSV / CM1 format)
#
# CSV example:  P(hPa), H(m), T(°C), Td(°C), WD(°), WS(kt)
# 1013.0, 0.0, 30.0, 22.0, 180, 10
# 1000.0, 112.0, 29.0, 21.0, 185, 12
# 925.0, 766.0, 24.0, 18.0, 200, 18
# ...
#
# SHARPpy:  P  H  T  Td  WD  WS  (space-separated)
# CM1:  sfc_p theta qv  (header) then  z theta qv u v`;

const IS_BINARY_FORMAT = (f) => f === "wrf";

export default function CustomUpload({ onBack, theme, colorblind }) {
  const [text, setText] = useState("");
  const [format, setFormat] = useState("auto");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [result, setResult] = useState(null);
  const [dragover, setDragover] = useState(false);
  const [binaryFile, setBinaryFile] = useState(null);
  const [wrfLat, setWrfLat] = useState("");
  const [wrfLon, setWrfLon] = useState("");
  const fileRef = useRef(null);
  const binaryRef = useRef(null);

  const isBinary = IS_BINARY_FORMAT(format);

  const handleFile = (file) => {
    if (!file) return;
    // Auto-detect WRF from extension
    const ext = file.name.split(".").pop().toLowerCase();
    if (["nc", "nc4", "ncf", "netcdf"].includes(ext) || file.name.startsWith("wrfout")) {
      setFormat("wrf");
      setBinaryFile(file);
      return;
    }
    const reader = new FileReader();
    reader.onload = (e) => setText(e.target.result);
    reader.readAsText(file);
  };

  const handleBinaryFile = (file) => {
    if (!file) return;
    setBinaryFile(file);
  };

  const handleDrop = (e) => {
    e.preventDefault();
    setDragover(false);
    const file = e.dataTransfer?.files?.[0];
    if (isBinary) {
      handleBinaryFile(file);
    } else {
      handleFile(file);
    }
  };

  const handleAnalyze = async () => {
    if (isBinary) {
      if (!binaryFile) return;
      setLoading(true);
      setError(null);
      setResult(null);
      try {
        const fd = new FormData();
        fd.append("file", binaryFile);
        fd.append("format", format);
        fd.append("theme", theme || "dark");
        fd.append("colorblind", colorblind ? "true" : "false");
        if (wrfLat.trim()) fd.append("lat", wrfLat.trim());
        if (wrfLon.trim()) fd.append("lon", wrfLon.trim());
        const data = await fetchUploadFile(fd);
        setResult(data);
      } catch (e) {
        setError(e.message);
      } finally {
        setLoading(false);
      }
    } else {
      if (!text.trim()) return;
      setLoading(true);
      setError(null);
      setResult(null);
      try {
        const data = await fetchCustomSounding({ text, format, theme, colorblind });
        setResult(data);
      } catch (e) {
        setError(e.message);
      } finally {
        setLoading(false);
      }
    }
  };

  return (
    <div className="cu-page">
      <button className="cu-back" onClick={onBack}>
        <ArrowLeft size={14} /> Back to Main
      </button>

      <div className="cu-container">
        <h2 className="cu-title">Custom Sounding Upload</h2>
        <p className="cu-subtitle">
          Paste or upload raw sounding profile data for analysis
        </p>

        {/* Format selection */}
        <div className="cu-format-row">
          {FORMATS.map((f) => (
            <button
              key={f.id}
              className={`cu-format-btn ${format === f.id ? "active" : ""}`}
              onClick={() => setFormat(f.id)}
            >
              {f.label}
            </button>
          ))}
        </div>

        {/* WRF-specific options */}
        {isBinary && (
          <div className="cu-wrf-options">
            <div className="cu-wrf-info">
              Upload a WRF netCDF output file (wrfout_d0x_*). Optionally specify
              a lat/lon point to extract; otherwise the domain center is used.
            </div>
            <div className="cu-wrf-coords">
              <label>
                <MapPin size={12} /> Lat
                <input
                  type="number"
                  step="0.01"
                  placeholder="e.g. 35.22"
                  value={wrfLat}
                  onChange={(e) => setWrfLat(e.target.value)}
                />
              </label>
              <label>
                <MapPin size={12} /> Lon
                <input
                  type="number"
                  step="0.01"
                  placeholder="e.g. -97.44"
                  value={wrfLon}
                  onChange={(e) => setWrfLon(e.target.value)}
                />
              </label>
            </div>
            <div
              className={`cu-dropzone ${dragover ? "dragover" : ""}`}
              onClick={() => binaryRef.current?.click()}
              onDragOver={(e) => { e.preventDefault(); setDragover(true); }}
              onDragLeave={() => setDragover(false)}
              onDrop={handleDrop}
            >
              <Upload size={16} />
              {binaryFile ? binaryFile.name : "Drop a WRF netCDF file here or click to browse"}
              <input
                ref={binaryRef}
                type="file"
                accept=".nc,.nc4,.ncf,.netcdf"
                onChange={(e) => handleBinaryFile(e.target.files?.[0])}
              />
            </div>
          </div>
        )}

        {/* Text area (non-binary formats) */}
        {!isBinary && (
          <>
            <textarea
              className="cu-textarea"
              placeholder={PLACEHOLDER}
              value={text}
              onChange={(e) => setText(e.target.value)}
              spellCheck="false"
            />

            {/* File drop zone */}
            <div
              className={`cu-dropzone ${dragover ? "dragover" : ""}`}
              onClick={() => fileRef.current?.click()}
              onDragOver={(e) => { e.preventDefault(); setDragover(true); }}
              onDragLeave={() => setDragover(false)}
              onDrop={handleDrop}
            >
              <FileText size={16} />
              Drop a file here or click to browse (.txt, .csv)
              <input
                ref={fileRef}
                type="file"
                accept=".txt,.csv,.dat"
                onChange={(e) => handleFile(e.target.files?.[0])}
              />
            </div>
          </>
        )}

        {/* Analyze button */}
        <button
          className="cu-analyze"
          disabled={(isBinary ? !binaryFile : !text.trim()) || loading}
          onClick={handleAnalyze}
        >
          {loading ? (
            <><Loader2 size={16} className="spin" /> Analyzing...</>
          ) : (
            <><Zap size={16} /> Analyze Sounding</>
          )}
        </button>

        {error && <div className="cu-error">{error}</div>}

        {/* Result */}
        {result && (
          <div className="cu-result">
            <h3 className="cu-result-title">Analysis Result</h3>
            {result.image && (
              <img
                src={`data:image/png;base64,${result.image}`}
                alt="Custom sounding analysis"
              />
            )}

            {/* Key params */}
            {result.params && (
              <div className="cu-params">
                {[
                  ["SB CAPE", result.params.sbCape, "J/kg"],
                  ["SB CIN", result.params.sbCin, "J/kg"],
                  ["MU CAPE", result.params.muCape, "J/kg"],
                  ["ML CAPE", result.params.mlCape, "J/kg"],
                  ["0-6km BWD", result.params.bwd6km, "kt"],
                  ["0-1km SRH", result.params.srh1km, "m²/s²"],
                  ["0-3km SRH", result.params.srh3km, "m²/s²"],
                  ["STP", result.params.stp, ""],
                  ["PWAT", result.params.pwat, "mm"],
                  ["LR 0-3", result.params.lr03, "°C/km"],
                ].map(([label, val, unit]) => (
                  <div key={label} className="cu-param-card">
                    <div className="cu-param-label">{label}</div>
                    <div className="cu-param-value">
                      {val != null ? val : "---"}{val != null && unit ? ` ${unit}` : ""}
                    </div>
                  </div>
                ))}
              </div>
            )}
          </div>
        )}

        {/* Help */}
        <div className="cu-help">
          <h4>Supported Formats</h4>
          <pre>{`CSV:      P(hPa), H(m), T(°C), Td(°C), WD(°), WS(kt)
SHARPpy:  P  H  T  Td  WD  WS  (whitespace-separated)
CM1:      Line 1: sfc_p(hPa) theta(K) qv(g/kg)
          Subsequent: z(m) theta(K) qv(g/kg) u(m/s) v(m/s)
WRF:      Binary netCDF wrfout file — select WRF netCDF format`}</pre>
          <p>
            Lines starting with # or % are treated as comments and skipped.
            At least 5 valid levels are required. Wind direction in degrees,
            wind speed in knots (CSV/SHARPpy) or m/s (CM1).
          </p>
          <p>
            <strong>WRF netCDF:</strong> Upload a wrfout_d0x file. The nearest
            grid point to your specified lat/lon is extracted (first time step).
            Variables used: P, PB, T, QVAPOR, PH, PHB, U, V, XLAT, XLONG.
          </p>
        </div>
      </div>
    </div>
  );
}
