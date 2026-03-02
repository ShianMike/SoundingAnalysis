import { useState, useRef, useCallback, useEffect } from "react";
import {
  ImageIcon,
  BarChart3,
  Info,
  Loader2,
  AlertTriangle,
  Download,
  Maximize2,
  ZoomIn,
  ZoomOut,
  Wind,
  Thermometer,
  ArrowUpDown,
  Droplets,
  Gauge,
  Zap,
  FileSpreadsheet,
  Link2,
  Check,
  FileText,
  ChevronDown,
  Printer,
  ArrowUp,
  ArrowDown,
  ShieldOff,
  Cloud,
  TrendingUp,
  RotateCcw,
  Target,
  Snowflake,
  Flame,
  CloudRain,
  ShieldAlert,
} from "lucide-react";

/* ── Icon map for summary section cards ─────────────── */
const SECTION_ICONS = {
  thermo: Thermometer,
  llcape: ArrowUp,
  cin: ShieldOff,
  lcl: Cloud,
  lapse: TrendingUp,
  moisture: Droplets,
  dlshear: Wind,
  llshear: RotateCcw,
  composite: BarChart3,
  dcape: ArrowDown,
  hail: Snowflake,
  mode: Target,
};
import StationMap from "./StationMap";
import TimeSeriesChart from "./TimeSeriesChart";
import ComparisonView from "./ComparisonView";
import VwpDisplay from "./VwpDisplay";
import "./ResultsView.css";

/* ── Parameter severity thresholds ─────────────────────────── */
const THRESHOLDS = {
  "SB CAPE":  [{ v: 3000, c: "extreme" }, { v: 1500, c: "high" }, { v: 500, c: "mod" }],
  "MU CAPE":  [{ v: 3000, c: "extreme" }, { v: 1500, c: "high" }, { v: 500, c: "mod" }],
  "ML CAPE":  [{ v: 3000, c: "extreme" }, { v: 1500, c: "high" }, { v: 500, c: "mod" }],
  "ECAPE":    [{ v: 2000, c: "extreme" }, { v: 1000, c: "high" }, { v: 400, c: "mod" }],
  "DCAPE":    [{ v: 1500, c: "extreme" }, { v: 800, c: "high" }],
  "3CAPE":    [{ v: 100, c: "extreme" }, { v: 50, c: "high" }],
  "STP":      [{ v: 4, c: "extreme" }, { v: 1, c: "high" }, { v: 0.5, c: "mod" }],
  "SCP":      [{ v: 8, c: "extreme" }, { v: 4, c: "high" }, { v: 1, c: "mod" }],
  "SHIP":     [{ v: 3, c: "extreme" }, { v: 1.5, c: "high" }, { v: 0.5, c: "mod" }],
  "DCP":      [{ v: 6, c: "extreme" }, { v: 4, c: "high" }, { v: 2, c: "mod" }],
  "LR 0-3 km":[{ v: 9, c: "extreme" }, { v: 8, c: "high" }, { v: 7, c: "mod" }],
  "LR 3-6 km":[{ v: 9, c: "extreme" }, { v: 8, c: "high" }, { v: 7, c: "mod" }],
  "BWD 0-6 km":[{ v: 60, c: "extreme" }, { v: 40, c: "high" }, { v: 25, c: "mod" }],
  "BWD 0-1 km":[{ v: 30, c: "extreme" }, { v: 20, c: "high" }, { v: 15, c: "mod" }],
  "SRH 0-1 km":[{ v: 300, c: "extreme" }, { v: 150, c: "high" }, { v: 100, c: "mod" }],
  "SRH 0-3 km":[{ v: 400, c: "extreme" }, { v: 200, c: "high" }, { v: 100, c: "mod" }],
  "Eff. SRH": [{ v: 300, c: "extreme" }, { v: 150, c: "high" }, { v: 100, c: "mod" }],
  "WMSI":     [{ v: 5, c: "extreme" }, { v: 3, c: "high" }, { v: 1, c: "mod" }],
  "MDPI":     [{ v: 3, c: "extreme" }, { v: 2, c: "high" }, { v: 1, c: "mod" }],
  "Max Gust": [{ v: 80, c: "extreme" }, { v: 58, c: "high" }, { v: 40, c: "mod" }],
  "Fosberg FWI": [{ v: 75, c: "extreme" }, { v: 50, c: "high" }, { v: 35, c: "mod" }],
  "Haines":   [{ v: 6, c: "extreme" }, { v: 5, c: "high" }, { v: 4, c: "mod" }],
  "HDW":      [{ v: 400, c: "extreme" }, { v: 200, c: "high" }, { v: 100, c: "mod" }],
};

function getAlertClass(label, value) {
  if (value == null || typeof value !== "number") return "";
  const rules = THRESHOLDS[label];
  if (!rules) return "";
  for (const r of rules) {
    if (value >= r.v) return `param-alert-${r.c}`;
  }
  return "";
}

/* return short severity label for badge */
function getSeverityLabel(alertClass) {
  if (!alertClass) return null;
  if (alertClass.includes("extreme")) return "EXT";
  if (alertClass.includes("high")) return "HIGH";
  if (alertClass.includes("mod")) return "MOD";
  return null;
}

/* badge CSS class from alert class */
function getBadgeClass(alertClass) {
  if (!alertClass) return "";
  if (alertClass.includes("extreme")) return "param-severity-badge--extreme";
  if (alertClass.includes("high")) return "param-severity-badge--high";
  if (alertClass.includes("mod")) return "param-severity-badge--mod";
  return "";
}

/* ── Sounding text summary generator ───────────────────────── */
function generateSoundingSummary(params, meta) {
  const sections = [];
  const n = (v) => (v != null && typeof v === "number" ? v : null);
  const fmt = (v, d = 0) => v != null ? Number(v).toFixed(d) : "--";

  // ── Gather all values up front ──
  const sbCape = n(params.sbCape);
  const mlCape = n(params.mlCape);
  const muCape = n(params.muCape);
  const ecape  = n(params.ecape);
  const cape3  = n(params.cape3km);
  const mlCin  = n(params.mlCin);
  const mlLcl  = n(params.mlLclM);
  const dcape  = n(params.dcape);
  const lr03   = n(params.lr03);
  const lr36   = n(params.lr36);
  const bwd5   = n(params.bwd500m);
  const bwd1   = n(params.bwd1km);
  const bwd6   = n(params.bwd6km);
  const srh5   = n(params.srh500m);
  const srh1   = n(params.srh1km);
  const srh3   = n(params.srh3km);
  const esrh   = n(params.esrh);
  const ebwd   = n(params.ebwd);
  const stp    = n(params.stp);
  const scp    = n(params.scp);
  const ship   = n(params.ship);
  const dcp    = n(params.dcp);
  const pwat   = n(params.pwat);
  const rh01   = n(params.rh01);
  const rh36   = n(params.rh36);
  const frzLvl = n(params.frzLevel);
  const wbo    = n(params.wbo);
  const bestCape = Math.max(sbCape ?? 0, mlCape ?? 0, muCape ?? 0);

  // ── Header ──
  const stName = meta?.stationName || meta?.station || "this location";
  const dateStr = meta?.date || "";

  // — Overall threat level for banner —
  let threatLevel = "none";
  let threatColor = "#6b7280";
  if (stp >= 4 || (bestCape >= 4000 && bwd6 >= 50)) { threatLevel = "HIGH"; threatColor = "#ef4444"; }
  else if (stp >= 1 || scp >= 4 || ship >= 2 || dcp >= 4) { threatLevel = "MODERATE"; threatColor = "#f59e0b"; }
  else if (bestCape >= 500 && bwd6 >= 25) { threatLevel = "LOW"; threatColor = "#22c55e"; }
  else if (bestCape >= 250) { threatLevel = "MARGINAL"; threatColor = "#60a5fa"; }

  // ── 1) Thermodynamic profile overview ──
  {
    let text;
    const vals = [];
    if (muCape != null) vals.push({ k: "MUCAPE", v: `${fmt(muCape)} J/kg` });
    if (mlCape != null) vals.push({ k: "MLCAPE", v: `${fmt(mlCape)} J/kg` });
    if (ecape != null) vals.push({ k: "ECAPE", v: `${fmt(ecape)} J/kg` });

    if (bestCape < 50) {
      text = "The thermodynamic profile is stable with virtually no buoyancy. Thunderstorm development is not supported. Any convection would require significant mesoscale or synoptic forcing and would likely remain shallow and non-severe.";
    } else if (bestCape < 250) {
      text = "Marginal instability is present. Buoyancy is insufficient for robust updrafts, though isolated weak convection could occur along well-defined convergence zones or terrain-enhanced lifting. Severe weather is unlikely.";
    } else if (bestCape < 1000) {
      text = "Moderate instability exists. This level of buoyancy can support organized thunderstorms with modest updraft strength, particularly where mesoscale forcing provides reliable initiation.";
    } else if (bestCape < 2500) {
      text = "Substantial instability is present. This provides ample energy for vigorous updrafts, and severe weather is likely if storms develop within a supportive kinematic environment.";
    } else if (bestCape < 4000) {
      text = "Large instability is present. This is a high-end thermodynamic environment capable of producing strong to violent updrafts, very large hail, and extreme rainfall rates.";
    } else {
      text = "Extreme instability is evident. This is a rare, upper-tier thermodynamic environment where updrafts may exceed 50 m/s, supporting giant hail, violent tornadoes if shear is sufficient, and flash-flood-producing rainfall.";
    }
    sections.push({ id: "thermo", title: "Thermodynamic Overview", text, vals });
  }

  // ── 2) Low-level CAPE / updraft acceleration ──
  if (cape3 != null && bestCape >= 250 && cape3 >= 50) {
    const vals = [{ k: "0-3km CAPE", v: `${fmt(cape3)} J/kg` }];
    let text;
    if (cape3 >= 100) {
      text = "Strong low-level updraft acceleration enhances stretching of low-level vertical vorticity — a key ingredient for tornado intensity, dynamically stretching near-surface rotation into tighter, faster-spinning vortices.";
    } else {
      text = "Moderate low-level updraft acceleration provides meaningful support for efficient vortex stretching and tornado development.";
    }
    sections.push({ id: "llcape", title: "Low-Level Buoyancy", text, vals });
  }

  // ── 3) Cap / CIN analysis ──
  if (mlCin != null) {
    const vals = [{ k: "MLCIN", v: `${fmt(mlCin)} J/kg` }];
    let text;
    if (mlCin > -10) {
      text = "Convective inhibition is negligible. Storms could fire readily with minimal forcing — this is an uncapped environment where widespread initiation is possible.";
    } else if (mlCin > -50) {
      text = "A weak cap is present. Convective initiation should occur fairly easily along boundaries, outflow, or modest terrain features.";
    } else if (mlCin > -150) {
      text = "A moderate cap is in place. Initiation requires focused mesoscale forcing. Storms that breach the cap may be explosive, and discrete supercells become more likely than cluster modes.";
    } else {
      text = "A strong cap exists, significantly inhibiting deep convection. Only intense forcing is likely to break through this inversion. If a storm initiates, the sudden release of stored energy could produce a violently explosive updraft.";
    }
    sections.push({ id: "cin", title: "Convective Inhibition", text, vals });
  }

  // ── 4) Cloud base / LCL ──
  if (mlLcl != null && bestCape >= 250) {
    const vals = [{ k: "MLLCL", v: `${fmt(mlLcl)} m AGL` }];
    let text;
    if (mlLcl < 800) {
      text = "Cloud bases are very low, strongly favorable for tornadogenesis. The short sub-cloud layer allows tight low-level rotation to connect efficiently to the mesocyclone aloft.";
    } else if (mlLcl < 1200) {
      text = "Cloud bases are relatively low, supporting tornado development. The shallow sub-cloud layer maintains coherent low-level rotation.";
    } else if (mlLcl < 1800) {
      text = "Cloud bases are moderately elevated, reducing tornado potential. The deeper sub-cloud layer enhances downdraft evaporation and outflow wind damage risk.";
    } else {
      text = "Cloud bases are high, strongly limiting tornado potential but promoting vigorous downdraft development through evaporative cooling. Damaging outflow winds and microbursts become the primary risk.";
    }
    sections.push({ id: "lcl", title: "Cloud Base Height", text, vals });
  }

  // ── 5) Lapse rates ──
  if ((lr03 != null || lr36 != null) && bestCape >= 100) {
    const vals = [];
    if (lr03 != null) vals.push({ k: "LR 0-3", v: `${fmt(lr03, 1)} °C/km` });
    if (lr36 != null) vals.push({ k: "LR 3-6", v: `${fmt(lr36, 1)} °C/km` });
    const parts = [];
    if (lr03 != null) {
      if (lr03 >= 9) parts.push("exceptionally steep 0–3 km lapse rates approaching dry-adiabatic, indicating intense low-level buoyancy");
      else if (lr03 >= 8) parts.push("very steep 0–3 km lapse rates indicating a well-mixed, nearly adiabatic boundary layer");
      else if (lr03 >= 7) parts.push("moderately steep 0–3 km lapse rates");
      else if (lr03 < 5.5) parts.push("relatively weak 0–3 km lapse rates suggesting a stable boundary layer");
    }
    if (lr36 != null) {
      if (lr36 >= 8.5) parts.push("extreme mid-level lapse rates contributing to exceptionally deep CAPE");
      else if (lr36 >= 7.5) parts.push("steep mid-level lapse rates enhancing CAPE depth");
      else if (lr36 >= 6.5) parts.push("moderate mid-level lapse rates");
    }
    if (parts.length > 0) {
      sections.push({ id: "lapse", title: "Lapse Rates", text: `The temperature profile shows ${parts.join(", and ")}.`, vals });
    }
  }

  // ── 6) Moisture profile ──
  {
    const vals = [];
    if (pwat != null) vals.push({ k: "PWAT", v: `${fmt(pwat)} mm` });
    if (rh01 != null) vals.push({ k: "RH 0-1km", v: `${fmt(rh01)}%` });
    if (rh36 != null) vals.push({ k: "RH 3-6km", v: `${fmt(rh36)}%` });
    const moistParts = [];
    if (rh01 != null && bestCape >= 100) {
      if (rh01 >= 85) moistParts.push("very moist low levels");
      else if (rh01 < 60) moistParts.push("relatively dry low levels");
    }
    if (rh36 != null && bestCape >= 100) {
      if (rh36 < 30) moistParts.push("a very dry mid-level layer enhancing downdraft production");
      else if (rh36 < 50) moistParts.push("moderately dry mid-levels");
    }
    if (pwat != null) {
      if (pwat >= 50) moistParts.push("extremely high precipitable water signaling major flash flood potential");
      else if (pwat >= 40) moistParts.push("high precipitable water increasing flash flood risk");
      else if (pwat >= 30) moistParts.push("moderate precipitable water");
    }
    if (moistParts.length > 0 && vals.length > 0) {
      sections.push({ id: "moisture", title: "Moisture Profile", text: `The moisture profile features ${moistParts.join("; ")}.`, vals });
    }
  }

  // ── 7) Deep-layer shear ──
  if (bwd6 != null && bestCape >= 100) {
    const vals = [{ k: "0-6km BWD", v: `${fmt(bwd6)} kt` }];
    if (ebwd != null) vals.push({ k: "Eff BWD", v: `${fmt(ebwd)} kt` });
    let text;
    if (bwd6 >= 70) {
      text = "Deep-layer shear is exceptionally strong. This extreme kinematic environment overwhelmingly favors discrete, long-track supercells with persistent mesocyclones. Storm longevity and severity are maximized.";
    } else if (bwd6 >= 50) {
      text = "Deep-layer shear is very strong, robustly supporting supercellular convection. Discrete supercells are the expected storm mode, with strong updraft-downdraft separation.";
    } else if (bwd6 >= 35) {
      text = "Moderate-to-strong deep-layer shear supports organized convection ranging from splitting supercells to organized multicellular systems.";
    } else if (bwd6 >= 20) {
      text = "Moderate deep-layer shear provides some storm organization. Multicell clusters with embedded supercell structures are possible.";
    } else {
      text = "Weak deep-layer shear favors disorganized convection — pulse storms or loosely organized multicells without persistent updraft rotation.";
    }
    sections.push({ id: "dlshear", title: "Deep-Layer Shear", text, vals });
  }

  // ── 8) Low-level shear & SRH ──
  if (bestCape >= 250 && (srh1 != null || bwd1 != null)) {
    const vals = [];
    if (bwd1 != null) vals.push({ k: "0-1km BWD", v: `${fmt(bwd1)} kt` });
    if (srh1 != null) vals.push({ k: "0-1km SRH", v: `${fmt(srh1)} m²/s²` });
    if (srh3 != null) vals.push({ k: "0-3km SRH", v: `${fmt(srh3)} m²/s²` });
    if (esrh != null) vals.push({ k: "Eff SRH", v: `${fmt(esrh)} m²/s²` });
    const llParts = [];
    if (bwd1 != null) {
      if (bwd1 >= 30) llParts.push("extreme low-level shear");
      else if (bwd1 >= 20) llParts.push("strong low-level shear");
      else if (bwd1 >= 15) llParts.push("moderate low-level shear");
    }
    if (srh1 != null) {
      if (srh1 >= 400) llParts.push("extreme SRH strongly correlated with violent (EF4-EF5) tornadoes");
      else if (srh1 >= 200) llParts.push("very significant SRH highly supportive of strong tornadoes");
      else if (srh1 >= 100) llParts.push("notable SRH favoring mesocyclone development");
      else if (srh1 >= 50) llParts.push("modest SRH");
    }
    if (llParts.length > 0) {
      sections.push({ id: "llshear", title: "Low-Level Kinematics", text: `Low-level kinematics show ${llParts.join(", with ")}.`, vals });
    }
  }

  // ── 9) Composite indices ──
  {
    const vals = [];
    const compLines = [];
    if (stp != null) vals.push({ k: "STP", v: `${fmt(stp, 1)}` });
    if (scp != null) vals.push({ k: "SCP", v: `${fmt(scp, 1)}` });
    if (ship != null) vals.push({ k: "SHIP", v: `${fmt(ship, 1)}` });
    if (dcp != null) vals.push({ k: "DCP", v: `${fmt(dcp, 1)}` });

    if (stp != null && bestCape >= 250) {
      if (stp >= 8) compLines.push("STP is in the top percentile — high-confidence setup for violent (EF3+) tornadoes.");
      else if (stp >= 4) compLines.push("STP is a high-end value indicating a well-above-average environment for significant (EF2+) tornadoes.");
      else if (stp >= 1) compLines.push("STP exceeds the significant tornado threshold, indicating ingredients are aligned for EF2+ tornadoes.");
      else if (stp >= 0.5) compLines.push("STP suggests a non-trivial tornado risk, particularly for brief or weak (EF0-EF1) tornadoes.");
    }
    if (scp != null && bestCape >= 250) {
      if (scp >= 10) compLines.push("SCP is extreme, strongly supporting long-lived discrete supercells.");
      else if (scp >= 4) compLines.push("SCP strongly favors supercellular convection with persistent mesocyclones.");
      else if (scp >= 1) compLines.push("SCP supports supercell development with moderate mesocyclone potential.");
    }
    if (ship != null && bestCape >= 500) {
      if (ship >= 2.5) compLines.push("SHIP is very high, indicating a robust environment for significant hail (≥2 in.).");
      else if (ship >= 1.0) compLines.push("SHIP exceeds the significant-hail threshold (≥1 in. hail expected).");
      else if (ship >= 0.5) compLines.push("SHIP indicates marginal significant hail potential.");
    }
    if (dcp != null && dcp >= 2) {
      if (dcp >= 6) compLines.push("DCP is extreme, strongly signaling a derecho or widespread damaging wind event.");
      else if (dcp >= 4) compLines.push("DCP supports organized long-lived wind events with potential for widespread damage.");
      else compLines.push("DCP indicates some potential for organized damaging wind events.");
    }

    if (compLines.length > 0 && vals.length > 0) {
      sections.push({ id: "composite", title: "Composite Indices", text: compLines.join(" "), vals });
    } else if (bestCape >= 500 && bwd6 >= 20 && vals.length > 0) {
      sections.push({ id: "composite", title: "Composite Indices", text: "Composite severe parameters remain below significant-severe thresholds. General thunderstorm hazards remain possible with any convection.", vals });
    }
  }

  // ── 10) DCAPE / Downdraft threat ──
  if (dcape != null && dcape >= 600 && bestCape >= 100) {
    const vals = [{ k: "DCAPE", v: `${fmt(dcape)} J/kg` }];
    let text;
    if (dcape >= 1500) {
      text = "Extreme downdraft energy indicates very strong potential for damaging outflow winds. Microbursts and macrobursts with gusts exceeding 80 kt are likely.";
    } else if (dcape >= 1000) {
      text = "Strong downdraft energy supports damaging outflow winds. Isolated downbursts with gusts of 50–70 kt are probable.";
    } else {
      text = "Moderate downdraft energy is sufficient for gusty outflow winds (40–55 kt), particularly where mid-level dry air is entrained.";
    }
    sections.push({ id: "dcape", title: "Downdraft Potential", text, vals });
  }

  // ── 11) Hail environment ──
  if (bestCape >= 500 && frzLvl != null && wbo != null && ship >= 0.5) {
    const vals = [{ k: "FRZ Level", v: `${fmt(frzLvl)} m` }, { k: "WB Zero", v: `${fmt(wbo)} m` }];
    let text = null;
    if (frzLvl > 4000 && wbo > 3000) {
      text = "Elevated freezing level and wet-bulb zero height allow a long fall through warm air. Some melting expected but strong updrafts still support large hail aloft.";
    } else if (frzLvl < 3000 && wbo < 2500) {
      text = "Low freezing level and wet-bulb zero mean reduced melting, favoring larger hail reaching the surface relative to updraft strength.";
    }
    if (text) sections.push({ id: "hail", title: "Hail Environment", text, vals });
  }

  // ── 12) Convective mode forecast ──
  if (bestCape >= 250) {
    let mode = "";
    if (bwd6 >= 40 && (srh1 >= 100 || scp >= 1)) {
      mode = "Discrete supercells";
    } else if (bwd6 >= 35 && dcp >= 3) {
      mode = "Supercells transitioning to a bowing line segment";
    } else if (bwd6 >= 25 && dcp >= 2) {
      mode = "Organized multicells with embedded bow echoes";
    } else if (bwd6 >= 25) {
      mode = "Multicells with possible supercell structures";
    } else if (bestCape >= 1000) {
      mode = "Pulse and weakly organized multicellular storms";
    } else {
      mode = "Isolated pulse convection";
    }
    sections.push({ id: "mode", title: "Expected Convective Mode", text: `${mode}.`, vals: [] });
  }

  // ── 13) Overall threat summary ──
  const threats = [];
  if (stp != null && stp >= 1) {
    threats.push(stp >= 4 ? "Strong to violent tornadoes (EF2+)" : "Significant tornadoes");
  } else if (stp != null && stp >= 0.3 && mlLcl < 1500) {
    threats.push("Brief/weak tornadoes");
  }
  if (ship != null && ship >= 1.5) threats.push("Significant hail (≥2 in.)");
  else if (ship != null && ship >= 0.5) threats.push("Large hail");
  else if (bestCape >= 2500 && bwd6 >= 30) threats.push("Large hail");
  if (dcp != null && dcp >= 4) threats.push("Widespread damaging winds / derecho");
  else if (dcp != null && dcp >= 2) threats.push("Damaging straight-line winds");
  else if (dcape != null && dcape >= 1000 && bwd6 >= 20) threats.push("Damaging outflow gusts");
  const wmsi = n(params.wmsi);
  const mdpi = n(params.mdpi);
  if (wmsi != null && wmsi >= 3) threats.push("Wet microburst potential");
  else if (mdpi != null && mdpi >= 2) threats.push("Microburst potential");
  if (pwat != null && pwat >= 40) threats.push("Flash flooding");
  const precipType = params.precipType;
  if (precipType && precipType !== "N/A" && precipType !== "Rain") {
    threats.push(`Winter precip: ${precipType}`);
  }

  return {
    station: stName,
    date: dateStr,
    threatLevel,
    threatColor,
    threats,
    sections,
  };
}

function ParamCard({ label, value, unit, color, desc }) {
  const alertCls = getAlertClass(label, value);
  const badge = getSeverityLabel(alertCls);
  const badgeCls = getBadgeClass(alertCls);
  return (
    <div className={`param-card ${alertCls}`} title="">
      {badge && <span className={`param-severity-badge ${badgeCls}`}>{badge}</span>}
      <span className="param-label">{label}</span>
      <span className="param-value" style={!alertCls && color ? { color } : {}}>
        {value ?? "---"}
      </span>
      {unit && value != null && <span className="param-unit">{unit}</span>}
      {desc && <span className="param-tooltip">{desc}</span>}
    </div>
  );
}

function ParamSection({ title, icon, children }) {
  return (
    <div className="param-section">
      <div className="param-section-header">
        {icon}
        <h3>{title}</h3>
      </div>
      <div className="param-grid">{children}</div>
    </div>
  );
}

function CompositeCard({ label, thresholdKey, value, unit, color, highlight, desc }) {
  const alertCls = getAlertClass(thresholdKey || label, value);
  const badge = getSeverityLabel(alertCls);
  const badgeCls = getBadgeClass(alertCls);
  return (
    <div className={`composite-card${highlight ? " composite-card--highlight" : ""}`}>
      {badge && <span className={`param-severity-badge ${badgeCls}`}>{badge}</span>}
      <span className="composite-label">{label}</span>
      <span className={`composite-value ${alertCls}`} style={!alertCls && color ? { color } : {}}>
        {value ?? "---"}
      </span>
      <span className="composite-unit">{unit}</span>
      {desc && <span className="param-tooltip">{desc}</span>}
    </div>
  );
}

function RiskTable({ riskData }) {
  if (!riskData || !riskData.stations || riskData.stations.length === 0) return null;

  return (
    <div className="rv-risk-table-wrap">
      <div className="rv-risk-table-header">
        <Zap size={14} />
        <h3>Severe Weather Risk Scan</h3>
        <span className="rv-risk-table-date">{riskData.date}</span>
        <span className="rv-risk-table-count">{riskData.stations.length} stations</span>
      </div>
      <div className="rv-risk-table-scroll">
        <table className="rv-risk-table">
          <thead>
            <tr>
              <th>#</th>
              <th>Station</th>
              <th>Name</th>
              <th className="rv-rt-num">STP</th>
              <th className="rv-rt-num">SCP</th>
              <th className="rv-rt-num">SHIP</th>
              <th className="rv-rt-num">DCP</th>
              <th className="rv-rt-num">CAPE</th>
              <th className="rv-rt-num">SRH</th>
              <th className="rv-rt-num">BWD</th>
            </tr>
          </thead>
          <tbody>
            {riskData.stations.map((s, i) => (
              <tr key={s.id} className={s.stp >= 1 ? "rv-rt-high" : s.stp >= 0.3 ? "rv-rt-med" : ""}>
                <td className="rv-rt-rank">{i + 1}</td>
                <td className="rv-rt-id">{s.id}</td>
                <td className="rv-rt-name">{s.name}</td>
                <td className="rv-rt-num">
                  <span className={`rv-rt-stp ${s.stp >= 1 ? "high" : s.stp >= 0.3 ? "med" : "low"}`}>
                    {s.stp.toFixed(2)}
                  </span>
                </td>
                <td className="rv-rt-num">
                  <span className={`rv-rt-stp ${s.scp >= 4 ? "high" : s.scp >= 1 ? "med" : "low"}`}>
                    {s.scp.toFixed(2)}
                  </span>
                </td>
                <td className="rv-rt-num">
                  <span className={`rv-rt-stp ${s.ship >= 1.5 ? "high" : s.ship >= 0.5 ? "med" : "low"}`}>
                    {s.ship.toFixed(2)}
                  </span>
                </td>
                <td className="rv-rt-num">
                  <span className={`rv-rt-stp ${s.dcp >= 4 ? "high" : s.dcp >= 2 ? "med" : "low"}`}>
                    {s.dcp.toFixed(2)}
                  </span>
                </td>
                <td className="rv-rt-num">{s.cape}</td>
                <td className="rv-rt-num">{s.srh}</td>
                <td className="rv-rt-num">{s.bwd}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

export default function ResultsView({ result, loading, error, riskData, showRisk, showMap, mapProps, showTimeSeries, onCloseTimeSeries, showCompare, onCloseCompare, showVwp, onCloseVwp, compareHistoryData, onCompareHistoryConsumed, stations, selectedStation, source, lastParams }) {
  if (error) {
    return (
      <div className="results-view">
        {showMap && mapProps && <StationMap {...mapProps} />}
        {showRisk && <RiskTable riskData={riskData} />}
        {showTimeSeries && (
          <TimeSeriesChart station={selectedStation} source={source} onClose={onCloseTimeSeries} />
        )}
        {showCompare && (
          <ComparisonView stations={stations || []} onClose={onCloseCompare} historyData={compareHistoryData} onHistoryConsumed={onCompareHistoryConsumed} />
        )}
        {showVwp && (
          <VwpDisplay stations={stations || []} selectedStation={selectedStation} onClose={onCloseVwp} />
        )}
        <div className="rv-state rv-error">
          <AlertTriangle size={24} />
          <div>
            <h3>Analysis Failed</h3>
            <p>{error}</p>
          </div>
        </div>
      </div>
    );
  }

  if (loading) {
    return (
      <div className="results-view">
        {showMap && mapProps && <StationMap {...mapProps} />}
        {showRisk && <RiskTable riskData={riskData} />}
        {showTimeSeries && (
          <TimeSeriesChart station={selectedStation} source={source} onClose={onCloseTimeSeries} />
        )}
        {showCompare && (
          <ComparisonView stations={stations || []} onClose={onCloseCompare} historyData={compareHistoryData} onHistoryConsumed={onCompareHistoryConsumed} />
        )}
        {showVwp && (
          <VwpDisplay stations={stations || []} selectedStation={selectedStation} onClose={onCloseVwp} />
        )}
        <div className="rv-state rv-loading">
          <Loader2 size={24} className="spin" />
          <div>
            <h3>Fetching & Analyzing</h3>
            <p>
              Retrieving sounding data, computing thermodynamic and kinematic
              parameters, and generating the analysis plot...
            </p>
          </div>
        </div>
      </div>
    );
  }

  if (!result) {
    return (
      <div className="results-view">
        {showMap && mapProps && <StationMap {...mapProps} />}
        {showRisk && <RiskTable riskData={riskData} />}
        {showTimeSeries && (
          <TimeSeriesChart station={selectedStation} source={source} onClose={onCloseTimeSeries} />
        )}
        {showCompare && (
          <ComparisonView stations={stations || []} onClose={onCloseCompare} historyData={compareHistoryData} onHistoryConsumed={onCompareHistoryConsumed} />
        )}
        {showVwp && (
          <VwpDisplay stations={stations || []} selectedStation={selectedStation} onClose={onCloseVwp} />
        )}
        {!riskData && (
          <div className="rv-state rv-empty">
            <Wind size={32} />
            <div>
              <h3>No Sounding Loaded</h3>
              <p>
                Select a data source and station, then click Generate Sounding
                to fetch and analyze upper-air data.
              </p>
            </div>
          </div>
        )}
      </div>
    );
  }

  const { image, params, meta } = result;
  const [zoomed, setZoomed] = useState(false);
  const [linkCopied, setLinkCopied] = useState(false);
  const [exportOpen, setExportOpen] = useState(false);
  const [showSummary, setShowSummary] = useState(false);
  const exportRef = useRef(null);
  const plotRef = useRef(null);
  const dragRef = useRef({ dragging: false, startX: 0, startY: 0, scrollLeft: 0, scrollTop: 0 });

  // Close export menu on outside click
  useEffect(() => {
    if (!exportOpen) return;
    const handler = (e) => {
      if (exportRef.current && !exportRef.current.contains(e.target)) setExportOpen(false);
    };
    document.addEventListener("mousedown", handler);
    return () => document.removeEventListener("mousedown", handler);
  }, [exportOpen]);

  const handleMouseDown = useCallback((e) => {
    if (!zoomed) return;
    const el = plotRef.current;
    if (!el) return;
    dragRef.current = {
      dragging: true,
      startX: e.clientX,
      startY: e.clientY,
      scrollLeft: el.scrollLeft,
      scrollTop: el.scrollTop,
    };
    el.style.cursor = "grabbing";
    e.preventDefault();
  }, [zoomed]);

  const handleMouseMove = useCallback((e) => {
    const d = dragRef.current;
    if (!d.dragging) return;
    const el = plotRef.current;
    if (!el) return;
    el.scrollLeft = d.scrollLeft - (e.clientX - d.startX);
    el.scrollTop = d.scrollTop - (e.clientY - d.startY);
  }, []);

  const handleMouseUp = useCallback(() => {
    dragRef.current.dragging = false;
    const el = plotRef.current;
    if (el && zoomed) el.style.cursor = "grab";
  }, [zoomed]);

  const handleDownload = () => {
    const link = document.createElement("a");
    link.href = `data:image/png;base64,${image}`;
    link.download = `sounding_${meta.station || "analysis"}_${meta.date.replace(/\s/g, "_")}.png`;
    link.click();
  };

  const handleFullscreen = () => {
    const w = window.open();
    w.document.write(
      `<html><head><title>Sounding - ${meta.station}</title>
       <style>body{margin:0;background:#0a0a0a;display:flex;align-items:center;justify-content:center;min-height:100vh}
       img{max-width:100%;height:auto}</style></head>
       <body><img src="data:image/png;base64,${image}" /></body></html>`
    );
  };

  const handleCsvExport = () => {
    const rows = [["Parameter", "Value", "Unit"]];
    const entries = [
      ["SB CAPE", params.sbCape, "J/kg"],
      ["SB CIN", params.sbCin, "J/kg"],
      ["SB LCL", params.sbLclM, "m AGL"],
      ["MU CAPE", params.muCape, "J/kg"],
      ["MU CIN", params.muCin, "J/kg"],
      ["MU LCL", params.muLclM, "m AGL"],
      ["ML CAPE", params.mlCape, "J/kg"],
      ["ML CIN", params.mlCin, "J/kg"],
      ["ML LCL", params.mlLclM, "m AGL"],
      ["DCAPE", params.dcape, "J/kg"],
      ["ECAPE", params.ecape, "J/kg"],
      ["STP", params.stp, ""],
      ["SCP", params.scp, ""],
      ["SHIP", params.ship, ""],
      ["DCP", params.dcp, ""],
      ["LR 0-3 km", params.lr03, "C/km"],
      ["LR 3-6 km", params.lr36, "C/km"],
      ["PWAT", params.pwat, "mm"],
      ["FRZ Level", params.frzLevel, "m AGL"],
      ["WB Zero", params.wbo, "m AGL"],
      ["RH 0-1 km", params.rh01, "%"],
      ["RH 1-3 km", params.rh13, "%"],
      ["RH 3-6 km", params.rh36, "%"],
      ["BWD 0-500m", params.bwd500m, "kt"],
      ["BWD 0-1 km", params.bwd1km, "kt"],
      ["BWD 0-3 km", params.bwd3km, "kt"],
      ["BWD 0-6 km", params.bwd6km, "kt"],
      ["SRH 500m", params.srh500m, "m2/s2"],
      ["SRH 0-1 km", params.srh1km, "m2/s2"],
      ["SRH 0-3 km", params.srh3km, "m2/s2"],
    ];
    entries.forEach(([name, val, unit]) => {
      rows.push([name, val != null ? String(val) : "", unit]);
    });

    const csv = rows.map((r) => r.map((c) => `"${c}"`).join(",")).join("\n");
    const blob = new Blob([csv], { type: "text/csv" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = `sounding_${meta.station || "analysis"}_${meta.date.replace(/\s/g, "_")}.csv`;
    link.click();
    URL.revokeObjectURL(url);
  };

  // ── SHARPpy format export ──────────────────────────────────
  const handleSharppy = () => {
    const profile = result.profile;
    if (!profile || profile.length === 0) return;
    const lines = [];
    lines.push(`%TITLE%`);
    lines.push(` ${meta.station || "XXXX"}   ${meta.date.replace(/\s/g, "").replace("-", "").replace("Z", "")}`);
    lines.push(``);
    lines.push(`   LEVEL       HGHT       TEMP       DWPT       WDIR       WSPD`);
    lines.push(`-------------------------------------------------------------------`);
    lines.push(`%RAW%`);
    for (const lv of profile) {
      const p = lv.p != null ? lv.p.toFixed(2) : "9999.00";
      const h = lv.h != null ? lv.h.toFixed(2) : "9999.00";
      const t = lv.t != null ? lv.t.toFixed(2) : "9999.00";
      const td = lv.td != null ? lv.td.toFixed(2) : "9999.00";
      const wd = lv.wd != null ? lv.wd.toFixed(2) : "9999.00";
      const ws = lv.ws != null ? lv.ws.toFixed(2) : "9999.00";
      lines.push(`${p},${h},${t},${td},${wd},${ws}`);
    }
    lines.push(`%END%`);
    const text = lines.join("\n");
    const blob = new Blob([text], { type: "text/plain" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = `${meta.station || "sounding"}_${meta.date.replace(/\s/g, "_")}.sharppy`;
    link.click();
    URL.revokeObjectURL(url);
  };

  // ── CM1 input_sounding export ──────────────────────────────
  const handleCm1 = () => {
    const profile = result.profile;
    if (!profile || profile.length === 0) return;
    const lines = [];
    // First line: sfc pressure (mb), sfc theta (K), sfc mixing ratio (g/kg)
    const sfcP = profile[0].p;
    const sfcT = profile[0].t + 273.15; // K
    const sfcTd = profile[0].td;
    // Compute sfc potential temperature: theta = T * (1000/p)^0.286
    const sfcTheta = sfcT * Math.pow(1000.0 / sfcP, 0.286);
    // Compute sfc mixing ratio from dewpoint and pressure (Bolton 1980)
    const es = 6.112 * Math.exp((17.67 * sfcTd) / (sfcTd + 243.5));
    const sfcQv = (621.97 * es) / (sfcP - es); // g/kg
    lines.push(`${sfcP.toFixed(2)}  ${sfcTheta.toFixed(2)}  ${sfcQv.toFixed(2)}`);
    // Subsequent lines: height AGL (m), theta (K), qv (g/kg), u (m/s), v (m/s)
    const sfcH = profile[0].h;
    for (const lv of profile) {
      const hAgl = lv.h - sfcH;
      const tk = lv.t + 273.15;
      const theta = tk * Math.pow(1000.0 / lv.p, 0.286);
      const e = 6.112 * Math.exp((17.67 * lv.td) / (lv.td + 243.5));
      const qv = (621.97 * e) / (lv.p - e);
      let u = 0, v = 0;
      if (lv.wd != null && lv.ws != null) {
        const wsMs = lv.ws * 0.51444; // kt → m/s
        const wdRad = (lv.wd * Math.PI) / 180;
        u = -wsMs * Math.sin(wdRad);
        v = -wsMs * Math.cos(wdRad);
      }
      lines.push(`${hAgl.toFixed(1)}  ${theta.toFixed(2)}  ${qv.toFixed(2)}  ${u.toFixed(2)}  ${v.toFixed(2)}`);
    }
    const text = lines.join("\n");
    const blob = new Blob([text], { type: "text/plain" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = `input_sounding_${meta.station || "sounding"}_${meta.date.replace(/\s/g, "_")}`;
    link.click();
    URL.revokeObjectURL(url);
  };

  const handleCopyLink = async () => {
    try {
      await navigator.clipboard.writeText(window.location.href);
      setLinkCopied(true);
      setTimeout(() => setLinkCopied(false), 2000);
    } catch {
      // Fallback for non-HTTPS
      const ta = document.createElement("textarea");
      ta.value = window.location.href;
      document.body.appendChild(ta);
      ta.select();
      document.execCommand("copy");
      document.body.removeChild(ta);
      setLinkCopied(true);
      setTimeout(() => setLinkCopied(false), 2000);
    }
  };

  return (
    <div className="results-view">
      {/* Meta bar */}
      <div className="rv-meta-bar">
        <div className="rv-meta-items">
          <span className="rv-meta-tag rv-meta-station">{meta.station || meta.source.toUpperCase()}</span>
          <span className="rv-meta-text">{meta.stationName}</span>
          <span className="rv-meta-sep" />
          <span className="rv-meta-text">{meta.date}</span>
          <span className="rv-meta-sep" />
          <span className="rv-meta-text">{meta.levels} levels</span>
          <span className="rv-meta-text rv-meta-dim">
            {meta.sfcPressure}&ndash;{meta.topPressure} hPa
          </span>
          {meta.vadRadar && (
            <span className="rv-meta-tag rv-meta-vad" title={meta.vadTime ? `VAD valid: ${meta.vadTime}` : ""}>
              VAD: {meta.vadRadar}
            </span>
          )}
        </div>
        <div className="rv-meta-actions">
          <button
            className={`rv-btn ${zoomed ? "rv-btn-active" : ""}`}
            onClick={() => setZoomed((z) => !z)}
            title={zoomed ? "Fit to width" : "Zoom to full resolution"}
          >
            {zoomed ? <ZoomOut size={14} /> : <ZoomIn size={14} />}
          </button>
          <button className="rv-btn" onClick={handleDownload} title="Download PNG">
            <Download size={14} />
          </button>
          <div className="rv-export-wrap" ref={exportRef}>
            <button className={`rv-btn ${exportOpen ? "rv-btn-active" : ""}`} onClick={() => setExportOpen((v) => !v)} title="Export sounding data">
              <Download size={14} />
              <span className="rv-btn-label">Export</span>
              <ChevronDown size={10} />
            </button>
            {exportOpen && (
              <div className="rv-export-menu">
                <button onClick={() => { handleCsvExport(); setExportOpen(false); }}>
                  <FileSpreadsheet size={13} />
                  <div>
                    <span className="rv-export-title">CSV Parameters</span>
                    <span className="rv-export-desc">All computed parameters in spreadsheet format</span>
                  </div>
                </button>
                <button onClick={() => { handleSharppy(); setExportOpen(false); }}>
                  <FileText size={13} />
                  <div>
                    <span className="rv-export-title">SHARPpy Format</span>
                    <span className="rv-export-desc">Raw profile for SHARPpy analysis software</span>
                  </div>
                </button>
                <button onClick={() => { handleCm1(); setExportOpen(false); }}>
                  <FileText size={13} />
                  <div>
                    <span className="rv-export-title">CM1 input_sounding</span>
                    <span className="rv-export-desc">Profile for Cloud Model 1 simulations</span>
                  </div>
                </button>
              </div>
            )}
          </div>
          <button className={`rv-btn ${linkCopied ? "rv-btn-active" : ""}`} onClick={handleCopyLink} title={linkCopied ? "Link copied!" : "Copy shareable link"}>
            {linkCopied ? <Check size={14} /> : <Link2 size={14} />}
          </button>
          <button className="rv-btn" onClick={handleFullscreen} title="Fullscreen">
            <Maximize2 size={14} />
          </button>
          <button className="rv-btn" onClick={() => window.print()} title="Print multi-panel layout">
            <Printer size={14} />
          </button>
        </div>
      </div>

      {/* Map + Risk scan table */}
      {showMap && mapProps && <StationMap {...mapProps} />}
      {showRisk && <RiskTable riskData={riskData} />}

      {/* Time-Series Chart */}
      {showTimeSeries && (
        <TimeSeriesChart station={selectedStation} source={source} onClose={onCloseTimeSeries} />
      )}

      {/* Comparison View */}
      {showCompare && (
        <ComparisonView stations={stations || []} onClose={onCloseCompare} historyData={compareHistoryData} onHistoryConsumed={onCompareHistoryConsumed} />
      )}

      {/* VWP Display */}
      {showVwp && (
        <VwpDisplay stations={stations || []} selectedStation={selectedStation} onClose={onCloseVwp} />
      )}

      {/* Plot — supports drag-to-pan (desktop) and touch pan/pinch-zoom (mobile) */}
      <div
        ref={plotRef}
        className={`rv-plot-wrap ${zoomed ? "rv-plot-zoomed" : ""}`}
        onMouseDown={handleMouseDown}
        onMouseMove={handleMouseMove}
        onMouseUp={handleMouseUp}
        onMouseLeave={handleMouseUp}
        onTouchStart={(e) => {
          if (e.touches.length === 1 && zoomed) {
            const t = e.touches[0];
            dragRef.current = {
              dragging: true,
              startX: t.clientX,
              startY: t.clientY,
              scrollLeft: plotRef.current?.scrollLeft || 0,
              scrollTop: plotRef.current?.scrollTop || 0,
            };
          }
        }}
        onTouchMove={(e) => {
          const d = dragRef.current;
          if (d.dragging && e.touches.length === 1) {
            const t = e.touches[0];
            const el = plotRef.current;
            if (el) {
              el.scrollLeft = d.scrollLeft - (t.clientX - d.startX);
              el.scrollTop = d.scrollTop - (t.clientY - d.startY);
            }
          }
        }}
        onTouchEnd={() => { dragRef.current.dragging = false; }}
        style={{ touchAction: zoomed ? "none" : "auto" }}
      >
        <img
          src={`data:image/png;base64,${image}`}
          alt="Sounding analysis plot"
          className="rv-plot-img"
          draggable={false}
        />
      </div>

      {/* Parameters */}
      <div className="rv-params">
        {/* ── Row 1: Instability (SB/MU/ML parcels) + Lapse Rates & Moisture ── */}
        <ParamSection
          title="Thermodynamic"
          icon={<Thermometer size={14} />}
        >
          <ParamCard label="SB CAPE" value={params.sbCape} unit="J/kg" color="#f97316" desc="Surface-Based CAPE — total buoyant energy for a parcel lifted from the surface. Higher values indicate stronger updraft potential. >1000 notable, >3000 extreme." />
          <ParamCard label="SB CIN" value={params.sbCin} unit="J/kg" desc="Surface-Based CIN — energy needed to lift a surface parcel to its LFC. More negative = stronger cap. Values < -50 often inhibit convection initiation." />
          <ParamCard label="SB LCL" value={params.sbLclM} unit="m AGL" desc="Surface-Based Lifted Condensation Level — height where a surface parcel reaches saturation. Lower LCL (<1000m) favors tornadoes." />
          <ParamCard label="MU CAPE" value={params.muCape} unit="J/kg" desc="Most-Unstable CAPE — CAPE computed for the parcel with highest θe in the lowest 300 hPa. Represents the maximum buoyancy available." />
          <ParamCard label="MU CIN" value={params.muCin} unit="J/kg" desc="Most-Unstable CIN — inhibition for the most-unstable parcel. Useful when elevated convection is possible." />
          <ParamCard label="MU LCL" value={params.muLclM} unit="m AGL" desc="Most-Unstable LCL — condensation height for the MU parcel. May differ from SB LCL when the most unstable parcel is elevated." />
          <ParamCard label="ML CAPE" value={params.mlCape} unit="J/kg" color="#d946ef" desc="Mixed-Layer CAPE — CAPE for a parcel representing the mean of the lowest 100 hPa. Best estimate for surface-based storms in a well-mixed boundary layer." />
          <ParamCard label="ML CIN" value={params.mlCin} unit="J/kg" desc="Mixed-Layer CIN — inhibition for the ML parcel. More representative than SB CIN in the afternoon when the boundary layer is mixed." />
          <ParamCard label="ML LCL" value={params.mlLclM} unit="m AGL" desc="Mixed-Layer LCL — cloud base height for the ML parcel. Lower values increase tornado probability; <1000m is favorable." />
        </ParamSection>

        <ParamSection
          title="Lapse Rates & Moisture"
          icon={<Droplets size={14} />}
        >
          <ParamCard label="LR 0-3 km" value={params.lr03} unit="C/km" desc="0–3 km Lapse Rate — temperature decrease per km in the lowest 3 km. Values >7°C/km are steep; >8°C/km nearly absolute unstable. Key for tornado environments." />
          <ParamCard label="LR 3-6 km" value={params.lr36} unit="C/km" desc="3–6 km Lapse Rate — mid-level lapse rate. Steeper values (>7°C/km) enhance CAPE and updraft strength. >8°C/km is extreme instability." />
          <ParamCard label="PWAT" value={params.pwat} unit="mm" desc="Precipitable Water — total column water vapor. Higher values mean more moisture available for heavy rainfall. >50mm is extremely high for CONUS." />
          <ParamCard label="FRZ Level" value={params.frzLevel} unit="m AGL" desc="Freezing Level — height of the 0°C isotherm AGL. Affects hail size (higher FRZ = more melting) and snow levels." />
          <ParamCard label="WB Zero" value={params.wbo} unit="m AGL" desc="Wet-Bulb Zero Height — where the wet-bulb temperature crosses 0°C. Better predictor of hail reaching the surface than the freezing level. <2500m favors large hail." />
          <ParamCard label="WCD" value={params.wcd} unit="m" color="#22d3ee" desc="Warm Cloud Depth — distance from LCL to freezing level. Critical for hail melting: shallow WCD (<2000m) means hail survives to surface; deep WCD (>3000m) favors complete melting. Also affects precipitation efficiency." />
          <ParamCard label="ML LFC" value={params.mlLfcM} unit="m AGL" color="#fb923c" desc="Mixed-Layer LFC height — how much lifting is needed to trigger convection from the mixed-layer parcel. Lower LFC = easier initiation. Standard for CI assessment." />
          <ParamCard label="ML EL" value={params.mlElM} unit="m AGL" color="#c084fc" desc="Mixed-Layer Equilibrium Level — top of the buoyant layer for the ML parcel. Indicates maximum cloud top height and depth of the convective updraft." />
          <ParamCard label="RH 0-1 km" value={params.rh01} unit="%" desc="Relative Humidity 0–1 km — low-level moisture. Higher values (>80%) favor tornado development by reducing evaporative cooling of rain." />
          <ParamCard label="RH 1-3 km" value={params.rh13} unit="%" desc="Relative Humidity 1–3 km — mid-low moisture. Dry layers here (<50%) enhance DCAPE and outflow potential." />
          <ParamCard label="RH 3-6 km" value={params.rh36} unit="%" desc="Relative Humidity 3–6 km — mid-level moisture. Dry air here entrains into storms, increasing evaporative cooling and downdraft strength." />
        </ParamSection>

        {/* ── Row 2: Kinematic + Downburst & Microburst ── */}
        <ParamSection
          title="Kinematic"
          icon={<ArrowUpDown size={14} />}
        >
          <ParamCard label="BWD 0-1 km" value={params.bwd1km} unit="kt" color="#ef4444" desc="0–1 km Bulk Wind Difference — magnitude of the wind shear vector over the lowest 1 km. >15 kt supports organized storms; >20 kt favors tornadoes." />
          <ParamCard label="BWD 0-3 km" value={params.bwd3km} unit="kt" color="#f97316" desc="0–3 km Bulk Wind Difference — shear in the low-to-mid levels. Important for mesocyclone development. >30 kt favors supercells." />
          <ParamCard label="BWD 0-6 km" value={params.bwd6km} unit="kt" color="#eab308" desc="0–6 km Bulk Wind Difference — deep-layer shear. The primary discriminator between organized and disorganized convection. >40 kt strongly favors supercells." />
          <ParamCard label="SRH 500m" value={params.srh500m} unit="m²/s²" desc="0–500m Storm-Relative Helicity — streamwise vorticity in the lowest 500m relative to storm motion. Key for tornado potential. >150 is significant." />
          <ParamCard label="SRH 0-1 km" value={params.srh1km} unit="m²/s²" color="#ef4444" desc="0–1 km Storm-Relative Helicity — measures the rotational potential of a storm's updraft in the lowest 1 km. >100 favors mesocyclones; >300 strongly favors tornadoes." />
          <ParamCard label="SRH 0-3 km" value={params.srh3km} unit="m²/s²" color="#f97316" desc="0–3 km Storm-Relative Helicity — total low-level rotational potential. >200 favors strong mesocyclones; >400 is extreme. Used in STP and SCP composites." />
          <ParamCard label="Eff. SRH" value={params.esrh} unit="m²/s²" color="#2dd4bf" desc="Effective SRH — storm-relative helicity computed within the effective inflow layer (where CAPE ≥ 100 and CIN > -250). More physically meaningful than fixed-layer SRH." />
          <ParamCard label="Eff. BWD" value={params.ebwd} unit="kt" color="#34d399" desc="Effective Bulk Wind Difference — shear from the effective inflow base to half the MU EL height. Better discriminator for supercells than fixed 0-6 km shear." />
          <ParamCard label="EIL Base" value={params.eilBot} unit="m AGL" color="#a7f3d0" desc="Effective Inflow Layer base — lowest level where CAPE ≥ 100 J/kg and CIN > -250 J/kg. Identifies the true inflow layer for storms." />
          <ParamCard label="EIL Top" value={params.eilTop} unit="m AGL" color="#a7f3d0" desc="Effective Inflow Layer top — highest contiguous level meeting the CAPE/CIN thresholds. Deeper EIL = deeper inflow available for storms." />
          <ParamCard label="Corfidi UPW" value={params.corfidiUpSpd} unit="kt" color="orange" desc="Corfidi Upwind vector speed (Corfidi 2003) — MCS upwind propagation component. Represents back-building tendency. Slower speeds favor training echoes and flash flooding." />
          <ParamCard label="Corfidi DNW" value={params.corfidiDnSpd} unit="kt" color="#ff4444" desc="Corfidi Downwind vector speed (Corfidi 2003) — MCS forward propagation. Faster speeds = fast-moving MCS; slower = quasi-stationary. Key for flash flood risk." />
        </ParamSection>

        <ParamSection
          title="Downburst & Microburst"
          icon={<ArrowDown size={14} />}
        >
          <ParamCard label="WMSI" value={params.wmsi} unit="" color="#f97316" desc="Wet Microburst Severity Index — approximated as CAPE × Γ0-3 / 1000. Higher values indicate stronger potential for wet microbursts. >3 is significant." />
          <ParamCard label="MDPI" value={params.mdpi} unit="" color="#fb923c" desc="Microburst Day Potential Index — θe deficit (surface minus minimum in 0-6 km) divided by 20. Values >1 indicate favorable conditions for microbursts." />
          <ParamCard label="Max Gust" value={params.maxGust} unit="kt" color="#ef4444" desc="Maximum estimated downburst gust speed — derived from √(2×DCAPE). Simple theoretical maximum; actual gusts may differ. >58 kt = severe." />
        </ParamSection>

        {/* ── Row 3: Composite Indices (full‑width, prominent) ── */}
        <div className="param-section param-section--composites">
          <div className="param-section-header">
            <Gauge size={14} />
            <h3>Composite Indices & Derived Parameters</h3>
          </div>
          <div className="param-grid param-grid--composites">
            <CompositeCard label="DCAPE" thresholdKey="DCAPE" value={params.dcape} unit="J/kg" desc="Downdraft CAPE — energy available for downdrafts. Higher values (>800) indicate strong outflow winds and potential for damaging gusts." />
            <CompositeCard label="ECAPE" thresholdKey="ECAPE" value={params.ecape} unit="J/kg" color="#06b6d4" desc="Entraining CAPE (Peters et al. 2023) — CAPE adjusted for entrainment. More physically realistic than standard CAPE." />
            <CompositeCard label="3CAPE" thresholdKey="3CAPE" value={params.cape3km} unit="J/kg" color="#fb923c" desc="0–3 km CAPE — buoyant energy in the lowest 3 km. Key for tornado intensity." />
            <CompositeCard label="6CAPE" value={params.cape6km} unit="J/kg" color="#facc15" desc="0–6 km CAPE — indicates how quickly an updraft develops in the mid-levels." />
            <CompositeCard label="DCIN" value={params.dcin} unit="J/kg" color="#818cf8" desc="Downdraft CIN — inhibition of downdrafts reaching the surface. Near 0 means downdrafts easily penetrate." />
            <CompositeCard label="NCAPE" value={params.ncape} unit="J/kg/m" color="#38bdf8" desc="Normalized CAPE — buoyancy intensity per unit depth. >0.3 is very buoyant." />
            <CompositeCard label="STP" thresholdKey="STP" value={params.stp} unit="Sig Tornado" color="#60a5fa" highlight desc="Significant Tornado Parameter (fixed-layer) — values ≥1 suggest significant (EF2+) tornado environment." />
            <CompositeCard label="STP-Eff" thresholdKey="STP" value={params.stpEff} unit="Eff Tornado" color="#818cf8" highlight desc="Effective-Layer STP (Thompson et al. 2012) — SPC's operational version. Values ≥1 favor significant tornadoes." />
            <CompositeCard label="SCP" thresholdKey="SCP" value={params.scp} unit="Supercell" color="#f59e0b" highlight desc="Supercell Composite Parameter — values ≥1 support supercells; >4 strongly favors discrete supercells." />
            <CompositeCard label="SHIP" thresholdKey="SHIP" value={params.ship} unit="Sig Hail" color="#10b981" highlight desc="Significant Hail Parameter — values ≥1 indicate potential; >1.5 strongly favors significant hail." />
            <CompositeCard label="DCP" thresholdKey="DCP" value={params.dcp} unit="Derecho" color="#a78bfa" highlight desc="Derecho Composite Parameter — values ≥2 suggest potential for long-lived wind events." />
          </div>
        </div>

        {/* ── Row 4: Winter / Fire / Hazard combined ── */}
        <ParamSection
          title="Winter Wx / Precip Type"
          icon={<CloudRain size={14} />}
        >
          <ParamCard label="Precip Type" value={params.precipType} unit="" color="#22d3ee" desc="Precipitation type derived from the Bourgouin (2000) method — classifies as Rain, Snow, Ice Pellets, or Freezing Rain based on warm-nose and cold-layer energy areas." />
          <ParamCard label="Warm Area" value={params.warmLayerEnergy} unit="J/kg" color="#f97316" desc="Warm-nose energy — integrated positive area above 0°C in the melting layer. Higher values (>13.2 J/kg) indicate complete melting of ice particles." />
          <ParamCard label="Cold Area" value={params.coldLayerEnergy} unit="J/kg" color="#60a5fa" desc="Cold-layer energy — integrated negative area below 0°C beneath the warm nose. Higher values indicate refreezing of melted precipitation into ice pellets." />
        </ParamSection>

        <ParamSection
          title="Fire Weather"
          icon={<Flame size={14} />}
        >
          <ParamCard label="Fosberg FWI" value={params.fosbergFwi} unit="" color="#f97316" desc="Fosberg Fire Weather Index — combines temperature, relative humidity, and wind speed. Scale 0-100; values >50 indicate high fire weather potential, >75 extreme." />
          <ParamCard label="Haines" value={params.haines} unit="" color="#ef4444" desc="Haines Index — stability and moisture measure for wildfire growth potential. Range 2-6; 5 = high potential, 6 = very high potential for large fire growth." />
          <ParamCard label="HDW" value={params.hdw} unit="" color="#fb923c" desc="Hot-Dry-Windy Index (Srock et al. 2018) — maximum VPD × wind speed in the lowest 500m. Higher values indicate conditions favorable for rapid fire spread." />
        </ParamSection>

        {params.hazards && params.hazards.length > 0 && (
          <ParamSection
            title="Hazard Assessment"
            icon={<ShieldAlert size={14} />}
          >
            {params.hazards.map((h, i) => (
              <ParamCard
                key={i}
                label={h.type}
                value={h.level}
                unit=""
                color={h.level === "HIGH" ? "#ef4444" : h.level === "MOD" ? "#f59e0b" : "#60a5fa"}
                desc={`Sounding-derived ${h.type.toLowerCase()} threat level based on composite parameter analysis.`}
              />
            ))}
          </ParamSection>
        )}

        {params.tempAdvection && params.tempAdvection.length > 0 && (
          <ParamSection
            title="Temperature Advection"
            icon={<TrendingUp size={14} />}
          >
            {params.tempAdvection.map((layer, i) => (
              <ParamCard
                key={i}
                label={layer.layer}
                value={`${layer.type} (${layer.turn > 0 ? "+" : ""}${layer.turn}°)`}
                unit=""
                color={layer.type === "WAA" ? "#f97316" : layer.type === "CAA" ? "#60a5fa" : "#888"}
                desc={`${layer.type === "WAA" ? "Warm Air Advection — wind veering" : layer.type === "CAA" ? "Cold Air Advection — wind backing" : "Neutral — minimal directional change"} in the ${layer.layer} layer.`}
              />
            ))}
          </ParamSection>
        )}
      </div>

      {/* Text Summary */}
      <div className="rv-summary-section">
        <button className="rv-summary-toggle" onClick={() => setShowSummary((s) => !s)}>
          <FileText size={14} />
          <span>Sounding Text Summary</span>
          <ChevronDown size={12} className={showSummary ? "rv-chev-open" : ""} />
        </button>
        {showSummary && (() => {
          const summary = generateSoundingSummary(params, meta);
          return (
            <div className="rv-summary-body">
              {/* ── Header banner ── */}
              <div className="rv-sum-banner">
                <div className="rv-sum-banner-left">
                  <span className="rv-sum-station">{summary.station}</span>
                  {summary.date && <span className="rv-sum-date">{summary.date}</span>}
                </div>
                {summary.threatLevel !== "none" && (
                  <span className="rv-sum-threat-badge" style={{ background: summary.threatColor + "22", color: summary.threatColor, borderColor: summary.threatColor + "44" }}>
                    {summary.threatLevel} RISK
                  </span>
                )}
              </div>

              {/* ── Threat callout ── */}
              {summary.threats.length > 0 && (
                <div className="rv-sum-callout" style={{ borderLeftColor: summary.threatColor }}>
                  <span className="rv-sum-callout-label">Primary Hazards</span>
                  <div className="rv-sum-callout-tags">
                    {summary.threats.map((t, i) => (
                      <span key={i} className="rv-sum-hazard-tag" style={{ background: summary.threatColor + "18", color: summary.threatColor }}>{t}</span>
                    ))}
                  </div>
                </div>
              )}

              {/* ── Analysis sections ── */}
              <div className="rv-sum-sections">
                {summary.sections.map((sec) => {
                  const IconComp = SECTION_ICONS[sec.id];
                  return (
                  <div key={sec.id} className="rv-sum-card">
                    <div className="rv-sum-card-hdr">
                      {IconComp && <IconComp size={14} className="rv-sum-card-icon" />}
                      <span className="rv-sum-card-title">{sec.title}</span>
                    </div>
                    {sec.vals && sec.vals.length > 0 && (
                      <div className="rv-sum-card-vals">
                        {sec.vals.map((v) => (
                          <span key={v.k} className="rv-sum-kv">
                            <span className="rv-sum-kv-k">{v.k}</span>
                            <span className="rv-sum-kv-v">{v.v}</span>
                          </span>
                        ))}
                      </div>
                    )}
                    <p className="rv-sum-card-text">{sec.text}</p>
                  </div>
                  );
                })}
              </div>
            </div>
          );
        })()}
      </div>


    </div>
  );
}
