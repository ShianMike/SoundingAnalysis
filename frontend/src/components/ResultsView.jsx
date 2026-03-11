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
  Droplets,
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
  Play,
  RefreshCw,
} from "lucide-react";

/* â”€â”€ Icon map for summary section cards â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */
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
import { lazy, Suspense } from "react";
const StationMap       = lazy(() => import("./StationMap"));
const TimeSeriesChart  = lazy(() => import("./TimeSeriesChart"));
const ComparisonView   = lazy(() => import("./ComparisonView"));
const VwpDisplay       = lazy(() => import("./VwpDisplay"));
const InteractiveSkewT = lazy(() => import("./InteractiveSkewT"));
const InteractiveHodograph = lazy(() => import("./InteractiveHodograph"));
const SoundingAnimator = lazy(() => import("./SoundingAnimator"));
const SoundingTimeline = lazy(() => import("./SoundingTimeline"));
import "./ResultsView.css";

/* â”€â”€ Sounding text summary generator â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */
function generateSoundingSummary(params, meta) {
  const sections = [];
  const n = (v) => (v != null && typeof v === "number" ? v : null);
  const fmt = (v, d = 0) => v != null ? Number(v).toFixed(d) : "--";

  // â”€â”€ Gather all values up front â”€â”€
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
  const ehi01  = n(params.ehi01);
  const ehi03  = n(params.ehi03);
  const vgp    = n(params.vgp);
  const critAngle = n(params.criticalAngle);
  const pwat   = n(params.pwat);
  const rh01   = n(params.rh01);
  const rh36   = n(params.rh36);
  const frzLvl = n(params.frzLevel);
  const wbo    = n(params.wbo);
  const bestCape = Math.max(sbCape ?? 0, mlCape ?? 0, muCape ?? 0);

  // â”€â”€ Header â”€â”€
  const stName = meta?.stationName || meta?.station || "this location";
  const dateStr = meta?.date || "";

  // â€” Overall threat level for banner â€”
  let threatLevel = "none";
  let threatColor = "#6b7280";
  if (stp >= 4 || (bestCape >= 4000 && bwd6 >= 50)) { threatLevel = "HIGH"; threatColor = "#ef4444"; }
  else if (stp >= 1 || scp >= 4 || ship >= 2 || dcp >= 4) { threatLevel = "MODERATE"; threatColor = "#f59e0b"; }
  else if (bestCape >= 500 && bwd6 >= 25) { threatLevel = "LOW"; threatColor = "#22c55e"; }
  else if (bestCape >= 250) { threatLevel = "MARGINAL"; threatColor = "#60a5fa"; }

  // â”€â”€ 1) Thermodynamic profile overview â”€â”€
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

  // â”€â”€ 2) Low-level CAPE / updraft acceleration â”€â”€
  if (cape3 != null && bestCape >= 250 && cape3 >= 50) {
    const vals = [{ k: "0-3km CAPE", v: `${fmt(cape3)} J/kg` }];
    let text;
    if (cape3 >= 100) {
      text = "Strong low-level updraft acceleration enhances stretching of low-level vertical vorticity â€” a key ingredient for tornado intensity, dynamically stretching near-surface rotation into tighter, faster-spinning vortices.";
    } else {
      text = "Moderate low-level updraft acceleration provides meaningful support for efficient vortex stretching and tornado development.";
    }
    sections.push({ id: "llcape", title: "Low-Level Buoyancy", text, vals });
  }

  // â”€â”€ 3) Cap / CIN analysis â”€â”€
  if (mlCin != null) {
    const vals = [{ k: "MLCIN", v: `${fmt(mlCin)} J/kg` }];
    let text;
    if (mlCin > -10) {
      text = "Convective inhibition is negligible. Storms could fire readily with minimal forcing â€” this is an uncapped environment where widespread initiation is possible.";
    } else if (mlCin > -50) {
      text = "A weak cap is present. Convective initiation should occur fairly easily along boundaries, outflow, or modest terrain features.";
    } else if (mlCin > -150) {
      text = "A moderate cap is in place. Initiation requires focused mesoscale forcing. Storms that breach the cap may be explosive, and discrete supercells become more likely than cluster modes.";
    } else {
      text = "A strong cap exists, significantly inhibiting deep convection. Only intense forcing is likely to break through this inversion. If a storm initiates, the sudden release of stored energy could produce a violently explosive updraft.";
    }
    sections.push({ id: "cin", title: "Convective Inhibition", text, vals });
  }

  // â”€â”€ 4) Cloud base / LCL â”€â”€
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

  // â”€â”€ 5) Lapse rates â”€â”€
  if ((lr03 != null || lr36 != null) && bestCape >= 100) {
    const vals = [];
    if (lr03 != null) vals.push({ k: "LR 0-3", v: `${fmt(lr03, 1)} Â°C/km` });
    if (lr36 != null) vals.push({ k: "LR 3-6", v: `${fmt(lr36, 1)} Â°C/km` });
    const parts = [];
    if (lr03 != null) {
      if (lr03 >= 9) parts.push("exceptionally steep 0â€“3 km lapse rates approaching dry-adiabatic, indicating intense low-level buoyancy");
      else if (lr03 >= 8) parts.push("very steep 0â€“3 km lapse rates indicating a well-mixed, nearly adiabatic boundary layer");
      else if (lr03 >= 7) parts.push("moderately steep 0â€“3 km lapse rates");
      else if (lr03 < 5.5) parts.push("relatively weak 0â€“3 km lapse rates suggesting a stable boundary layer");
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

  // â”€â”€ 6) Moisture profile â”€â”€
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

  // â”€â”€ 7) Deep-layer shear â”€â”€
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
      text = "Weak deep-layer shear favors disorganized convection â€” pulse storms or loosely organized multicells without persistent updraft rotation.";
    }
    sections.push({ id: "dlshear", title: "Deep-Layer Shear", text, vals });
  }

  // â”€â”€ 8) Low-level shear & SRH â”€â”€
  if (bestCape >= 250 && (srh1 != null || bwd1 != null)) {
    const vals = [];
    if (bwd1 != null) vals.push({ k: "0-1km BWD", v: `${fmt(bwd1)} kt` });
    if (srh1 != null) vals.push({ k: "0-1km SRH", v: `${fmt(srh1)} mÂ²/sÂ²` });
    if (srh3 != null) vals.push({ k: "0-3km SRH", v: `${fmt(srh3)} mÂ²/sÂ²` });
    if (esrh != null) vals.push({ k: "Eff SRH", v: `${fmt(esrh)} mÂ²/sÂ²` });
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

  // â”€â”€ 9) Composite indices â”€â”€
  {
    const vals = [];
    const compLines = [];
    if (stp != null) vals.push({ k: "STP", v: `${fmt(stp, 1)}` });
    if (scp != null) vals.push({ k: "SCP", v: `${fmt(scp, 1)}` });
    if (ship != null) vals.push({ k: "SHIP", v: `${fmt(ship, 1)}` });
    if (dcp != null) vals.push({ k: "DCP", v: `${fmt(dcp, 1)}` });
    if (ehi01 != null && ehi01 > 0) vals.push({ k: "EHI 0-1", v: `${fmt(ehi01, 2)}` });
    if (ehi03 != null && ehi03 > 0) vals.push({ k: "EHI 0-3", v: `${fmt(ehi03, 2)}` });
    if (vgp != null && vgp > 0) vals.push({ k: "VGP", v: `${fmt(vgp, 3)}` });
    if (critAngle != null) vals.push({ k: "Crit Angle", v: `${fmt(critAngle, 0)}\u00b0` });

    if (stp != null && bestCape >= 250) {
      if (stp >= 8) compLines.push("STP is in the top percentile â€” high-confidence setup for violent (EF3+) tornadoes.");
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
      if (ship >= 2.5) compLines.push("SHIP is very high, indicating a robust environment for significant hail (â‰¥2 in.).");
      else if (ship >= 1.0) compLines.push("SHIP exceeds the significant-hail threshold (â‰¥1 in. hail expected).");
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

  // â”€â”€ 10) DCAPE / Downdraft threat â”€â”€
  if (dcape != null && dcape >= 600 && bestCape >= 100) {
    const vals = [{ k: "DCAPE", v: `${fmt(dcape)} J/kg` }];
    let text;
    if (dcape >= 1500) {
      text = "Extreme downdraft energy indicates very strong potential for damaging outflow winds. Microbursts and macrobursts with gusts exceeding 80 kt are likely.";
    } else if (dcape >= 1000) {
      text = "Strong downdraft energy supports damaging outflow winds. Isolated downbursts with gusts of 50â€“70 kt are probable.";
    } else {
      text = "Moderate downdraft energy is sufficient for gusty outflow winds (40â€“55 kt), particularly where mid-level dry air is entrained.";
    }
    sections.push({ id: "dcape", title: "Downdraft Potential", text, vals });
  }

  // â”€â”€ 11) Hail environment â”€â”€
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

  // â”€â”€ 12) Convective mode forecast â”€â”€
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

  // â”€â”€ 13) Overall threat summary â”€â”€
  const threats = [];
  if (stp != null && stp >= 1) {
    threats.push(stp >= 4 ? "Strong to violent tornadoes (EF2+)" : "Significant tornadoes");
  } else if (stp != null && stp >= 0.3 && mlLcl < 1500) {
    threats.push("Brief/weak tornadoes");
  }
  if (ship != null && ship >= 1.5) threats.push("Significant hail (â‰¥2 in.)");
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

function RiskTable({ riskData, onStationSelect }) {
  if (!riskData || !riskData.stations || riskData.stations.length === 0) return null;

  const handleExportPng = () => {
    const stations = riskData.stations.slice(0, 15);
    const cols = ["#", "Station", "Name", "STP", "SCP", "SHIP", "DCP", "CAPE", "SRH", "BWD"];
    const colW = [50, 80, 220, 80, 80, 80, 80, 80, 70, 70];
    const rowH = 32;
    const headerH = 70;
    const pad = 24;
    const totalW = colW.reduce((a, b) => a + b, 0) + pad * 2;
    const totalH = headerH + rowH * (stations.length + 1) + pad * 2 + 4;

    const canvas = document.createElement("canvas");
    const scale = 3;
    canvas.width = totalW * scale;
    canvas.height = totalH * scale;
    const ctx = canvas.getContext("2d");
    ctx.scale(scale, scale);

    // Background
    ctx.fillStyle = "#0f172a";
    ctx.fillRect(0, 0, totalW, totalH);

    // Header
    ctx.fillStyle = "#e2e8f0";
    ctx.font = "bold 16px system-ui, sans-serif";
    ctx.fillText("Severe Weather Risk Scan", pad, pad + 22);
    ctx.font = "13px monospace";
    ctx.fillStyle = "#94a3b8";
    let subtitle = riskData.date;
    if (riskData.model) {
      subtitle = `${riskData.model} F${riskData.fhour} Â· Init ${riskData.initTime} Â· Valid ${riskData.date}`;
    }
    ctx.fillText(subtitle, pad, pad + 46);
    ctx.fillText(`${stations.length} stations`, totalW - pad - ctx.measureText(`${stations.length} stations`).width, pad + 46);

    const tableY = headerH + pad;

    // Column header row
    ctx.fillStyle = "#1e293b";
    ctx.fillRect(pad, tableY, totalW - pad * 2, rowH);
    ctx.font = "bold 11px system-ui, sans-serif";
    ctx.fillStyle = "#94a3b8";
    let cx = pad;
    cols.forEach((c, ci) => {
      const align = ci >= 3 ? "right" : "left";
      if (align === "right") {
        ctx.textAlign = "right";
        ctx.fillText(c, cx + colW[ci] - 6, tableY + 21);
      } else {
        ctx.textAlign = "left";
        ctx.fillText(c, cx + 6, tableY + 21);
      }
      cx += colW[ci];
    });

    // Rows
    ctx.font = "13px monospace";
    stations.forEach((s, i) => {
      const y = tableY + rowH * (i + 1);
      // Alternating row bg
      if (i % 2 === 0) {
        ctx.fillStyle = "rgba(255,255,255,0.02)";
        ctx.fillRect(pad, y, totalW - pad * 2, rowH);
      }
      // Row highlight
      if (s.stp >= 1) {
        ctx.fillStyle = "rgba(239,68,68,0.08)";
        ctx.fillRect(pad, y, totalW - pad * 2, rowH);
      } else if (s.stp >= 0.3) {
        ctx.fillStyle = "rgba(234,179,8,0.06)";
        ctx.fillRect(pad, y, totalW - pad * 2, rowH);
      }

      const vals = [
        String(i + 1), s.id, s.name,
        s.stp.toFixed(2), s.scp.toFixed(2), s.ship.toFixed(2), s.dcp.toFixed(2),
        String(s.cape), String(s.srh), String(s.bwd),
      ];
      const colorMap = [
        null, null, null,
        s.stp >= 1 ? "#ef4444" : s.stp >= 0.3 ? "#eab308" : "#94a3b8",
        s.scp >= 4 ? "#ef4444" : s.scp >= 1 ? "#eab308" : "#94a3b8",
        s.ship >= 1.5 ? "#ef4444" : s.ship >= 0.5 ? "#eab308" : "#94a3b8",
        s.dcp >= 4 ? "#ef4444" : s.dcp >= 2 ? "#eab308" : "#94a3b8",
        null, null, null,
      ];

      let vx = pad;
      vals.forEach((v, ci) => {
        ctx.fillStyle = colorMap[ci] || (ci <= 2 ? "#cbd5e1" : "#94a3b8");
        ctx.font = ci >= 3 && colorMap[ci] && colorMap[ci] !== "#94a3b8"
          ? "bold 13px monospace" : "13px monospace";
        if (ci === 1) ctx.font = "bold 13px monospace";
        const align = ci >= 3 ? "right" : "left";
        if (align === "right") {
          ctx.textAlign = "right";
          ctx.fillText(v, vx + colW[ci] - 6, y + 21);
        } else {
          ctx.textAlign = "left";
          ctx.fillText(v, vx + 6, y + 21);
        }
        vx += colW[ci];
      });
    });

    // Border
    ctx.strokeStyle = "rgba(255,255,255,0.08)";
    ctx.lineWidth = 1;
    ctx.strokeRect(pad, tableY, totalW - pad * 2, rowH * (stations.length + 1));

    // Download
    const link = document.createElement("a");
    const tag = riskData.model ? `${riskData.model}_F${riskData.fhour}` : "observed";
    link.download = `risk_scan_${tag}_${riskData.date.replace(/[:\s]/g, "")}.png`;
    link.href = canvas.toDataURL("image/png");
    link.click();
  };

  return (
    <div className="rv-risk-table-wrap">
      <div className="rv-risk-table-header">
        <Zap size={14} />
        <h3>Severe Weather Risk Scan</h3>
        {riskData.model && (
          <span className="rv-risk-table-model">{riskData.model} F{riskData.fhour}</span>
        )}
        <span className="rv-risk-table-date">{riskData.date}</span>
        <span className="rv-risk-table-count">{riskData.stations.length} stations</span>
        <button
          type="button"
          className="rv-risk-export-btn"
          onClick={handleExportPng}
          title="Export as PNG"
        >
          <Download size={13} />
        </button>
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
              <tr
                key={s.id}
                className={`${s.stp >= 1 ? "rv-rt-high" : s.stp >= 0.3 ? "rv-rt-med" : ""} rv-rt-clickable`}
                onClick={() => onStationSelect?.(s.id, riskData)}
                title={`Load ${s.id} sounding${riskData.model ? ` (${riskData.model} F${riskData.fhour})` : ""}`}
              >
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

export default function ResultsView({ result, loading, error, riskData, showRisk, showMap, mapProps, showTimeSeries, onCloseTimeSeries, showCompare, onCloseCompare, showVwp, onCloseVwp, compareHistoryData, onCompareHistoryConsumed, stations, selectedStation, source, lastParams, autoRefresh, onToggleAutoRefresh, refreshInterval, onRefreshIntervalChange, theme, onTimelineSelect, onRiskStationSelect }) {
  /* â”€â”€ Hooks â€” must be called unconditionally before any return â”€â”€ */
  const [zoomed, setZoomed] = useState(false);
  const [linkCopied, setLinkCopied] = useState(false);
  const [exportOpen, setExportOpen] = useState(false);
  const [showSummary, setShowSummary] = useState(false);
  const [interactiveMode, setInteractiveMode] = useState(false);
  const [showAnimator, setShowAnimator] = useState(false);
  const [reportBusy, setReportBusy] = useState(false);
  const exportRef = useRef(null);
  const plotRef = useRef(null);
  const dragRef = useRef({ dragging: false, startX: 0, startY: 0, scrollLeft: 0, scrollTop: 0 });

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

  if (error) {
    return (
      <div className="results-view">
        <Suspense fallback={null}>
          {showMap && mapProps && <div id="section-map"><StationMap {...mapProps} /></div>}
        </Suspense>
        {showRisk && <div id="section-risk"><RiskTable riskData={riskData} onStationSelect={onRiskStationSelect} /></div>}
        <Suspense fallback={null}>
          {showTimeSeries && (
            <div id="section-timeseries"><TimeSeriesChart station={selectedStation} source={source} onClose={onCloseTimeSeries} /></div>
          )}
          {showCompare && (
            <div id="section-compare"><ComparisonView stations={stations || []} onClose={onCloseCompare} historyData={compareHistoryData} onHistoryConsumed={onCompareHistoryConsumed} /></div>
          )}
          {showVwp && (
            <div id="section-vwp"><VwpDisplay stations={stations || []} selectedStation={selectedStation} onClose={onCloseVwp} /></div>
          )}
        </Suspense>
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
        <Suspense fallback={null}>
          {showMap && mapProps && <div id="section-map"><StationMap {...mapProps} /></div>}
        </Suspense>
        {showRisk && <div id="section-risk"><RiskTable riskData={riskData} onStationSelect={onRiskStationSelect} /></div>}
        <Suspense fallback={null}>
          {showTimeSeries && (
            <div id="section-timeseries"><TimeSeriesChart station={selectedStation} source={source} onClose={onCloseTimeSeries} /></div>
          )}
          {showCompare && (
            <div id="section-compare"><ComparisonView stations={stations || []} onClose={onCloseCompare} historyData={compareHistoryData} onHistoryConsumed={onCompareHistoryConsumed} /></div>
          )}
          {showVwp && (
            <div id="section-vwp"><VwpDisplay stations={stations || []} selectedStation={selectedStation} onClose={onCloseVwp} /></div>
          )}
        </Suspense>
        <div id="section-sounding" className="rv-state rv-loading">
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
        <Suspense fallback={null}>
          {showMap && mapProps && <div id="section-map"><StationMap {...mapProps} /></div>}
        </Suspense>
        {showRisk && <div id="section-risk"><RiskTable riskData={riskData} onStationSelect={onRiskStationSelect} /></div>}
        <Suspense fallback={null}>
          {showTimeSeries && (
            <div id="section-timeseries"><TimeSeriesChart station={selectedStation} source={source} onClose={onCloseTimeSeries} /></div>
          )}
          {showCompare && (
            <div id="section-compare"><ComparisonView stations={stations || []} onClose={onCloseCompare} historyData={compareHistoryData} onHistoryConsumed={onCompareHistoryConsumed} /></div>
          )}
          {showVwp && (
            <div id="section-vwp"><VwpDisplay stations={stations || []} selectedStation={selectedStation} onClose={onCloseVwp} /></div>
          )}
        </Suspense>
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

  // â”€â”€ SHARPpy format export â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
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

  // â”€â”€ CM1 input_sounding export â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
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
        const wsMs = lv.ws * 0.51444; // kt â†’ m/s
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

  // â”€â”€ JSON export â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
  const handleJsonExport = () => {
    const payload = {
      meta: { station: meta.station, date: meta.date, source: meta.source, stationName: meta.stationName },
      params,
      profile: result.profile,
    };
    const json = JSON.stringify(payload, null, 2);
    const blob = new Blob([json], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = `sounding_${meta.station || "analysis"}_${meta.date.replace(/\s/g, "_")}.json`;
    link.click();
    URL.revokeObjectURL(url);
  };

  const handleFullReport = async () => {
    if (reportBusy) return;
    setReportBusy(true);
    try {
      // Load the sounding image
      const img = new Image();
      img.src = `data:image/png;base64,${image}`;
      await new Promise((res, rej) => { img.onload = res; img.onerror = rej; });

      // Build parameter summary rows
      const rows = [
        [`${meta.station || meta.source.toUpperCase()}  â€”  ${meta.stationName || ""}  â€”  ${meta.date}`, ""],
        ["", ""],
        ["THERMODYNAMIC", ""],
        ["SB CAPE", `${params.sbCape ?? "---"} J/kg`], ["SB CIN", `${params.sbCin ?? "---"} J/kg`],
        ["MU CAPE", `${params.muCape ?? "---"} J/kg`], ["ML CAPE", `${params.mlCape ?? "---"} J/kg`],
        ["ML CIN", `${params.mlCin ?? "---"} J/kg`], ["DCAPE", `${params.dcape ?? "---"} J/kg`],
        ["ECAPE", `${params.ecape ?? "---"} J/kg`],
        ["", ""],
        ["LAPSE RATES & MOISTURE", ""],
        ["LR 0-3 km", `${params.lr03 ?? "---"} C/km`], ["LR 3-6 km", `${params.lr36 ?? "---"} C/km`],
        ["PWAT", `${params.pwat ?? "---"} mm`], ["FRZ Level", `${params.frzLevel ?? "---"} m`],
        ["", ""],
        ["KINEMATIC", ""],
        ["BWD 0-1 km", `${params.bwd1km ?? "---"} kt`], ["BWD 0-6 km", `${params.bwd6km ?? "---"} kt`],
        ["SRH 0-1 km", `${params.srh1km ?? "---"} mÂ²/sÂ²`], ["SRH 0-3 km", `${params.srh3km ?? "---"} mÂ²/sÂ²`],
        ["", ""],
        ["COMPOSITE INDICES", ""],
        ["STP", `${params.stp ?? "---"}`], ["SCP", `${params.scp ?? "---"}`],
        ["SHIP", `${params.ship ?? "---"}`], ["DCP", `${params.dcp ?? "---"}`],
        ["BRN", `${params.brn ?? "---"}`],
        ["EHI 0-1", `${params.ehi01 ?? "---"}`], ["EHI 0-3", `${params.ehi03 ?? "---"}`],
        ["VGP", `${params.vgp ?? "---"}`], ["Crit Angle", `${params.criticalAngle != null ? params.criticalAngle + "\u00b0" : "---"}`],
      ];

      const pad = 40;
      const lineH = 20;
      const panelW = 420;
      const panelH = rows.length * lineH + pad * 2;
      const totalW = img.width + panelW;
      const totalH = Math.max(img.height, panelH);

      const canvas = document.createElement("canvas");
      canvas.width = totalW;
      canvas.height = totalH;
      const ctx = canvas.getContext("2d");

      // Background
      ctx.fillStyle = "#0a0a14";
      ctx.fillRect(0, 0, totalW, totalH);

      // Draw sounding image
      ctx.drawImage(img, 0, 0);

      // Draw parameter panel
      let y = pad;
      for (const [label, val] of rows) {
        if (val === "" && label !== "") {
          // Section header
          ctx.fillStyle = "#60a5fa";
          ctx.font = "bold 13px monospace";
          ctx.fillText(label, img.width + pad, y);
        } else if (label !== "") {
          ctx.fillStyle = "#94a3b8";
          ctx.font = "12px monospace";
          ctx.fillText(label, img.width + pad, y);
          ctx.fillStyle = "#e2e8f0";
          ctx.font = "bold 12px monospace";
          ctx.fillText(val, img.width + pad + 180, y);
        }
        y += lineH;
      }

      const link = document.createElement("a");
      link.href = canvas.toDataURL("image/png");
      link.download = `report_${meta.station || "analysis"}_${meta.date.replace(/\s/g, "_")}.png`;
      link.click();
    } catch (err) {
      console.error("Report export failed:", err);
    }
    setReportBusy(false);
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
      {/* Map + Risk scan table */}
      <Suspense fallback={null}>
        {showMap && mapProps && <div id="section-map"><StationMap {...mapProps} /></div>}
      </Suspense>
      {showRisk && <div id="section-risk"><RiskTable riskData={riskData} onStationSelect={onRiskStationSelect} /></div>}

      {/* Time-Series Chart / Comparison / VWP */}
      <Suspense fallback={null}>
        {showTimeSeries && (
          <div id="section-timeseries"><TimeSeriesChart station={selectedStation} source={source} onClose={onCloseTimeSeries} /></div>
        )}
        {showCompare && (
          <div id="section-compare"><ComparisonView stations={stations || []} onClose={onCloseCompare} historyData={compareHistoryData} onHistoryConsumed={onCompareHistoryConsumed} /></div>
        )}
        {showVwp && (
          <div id="section-vwp"><VwpDisplay stations={stations || []} selectedStation={selectedStation} onClose={onCloseVwp} /></div>
        )}
      </Suspense>

      {/* Meta bar â€” above the sounding plot */}
      <div id="section-sounding" className="rv-meta-bar">
        <div className="rv-meta-row-top">
          <span className="rv-meta-station">{meta.station || meta.source.toUpperCase()}</span>
          <span className="rv-meta-name">{meta.stationName}</span>
          <span className="rv-meta-date">{meta.date}</span>
        </div>
        <div className="rv-meta-row-bottom">
          <div className="rv-meta-details">
            <span className="rv-meta-chip">{meta.levels} <span className="rv-meta-chip-value">levels</span></span>
            <span className="rv-meta-chip">{meta.sfcPressure}&ndash;{meta.topPressure} <span className="rv-meta-chip-value">hPa</span></span>
            {meta.vadRadar && (
              <span className="rv-meta-chip rv-meta-vad" title={meta.vadTime ? `VAD valid: ${meta.vadTime}` : ""}>
                VAD <span className="rv-meta-chip-value">{meta.vadRadar}</span>
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
                <button onClick={() => { handleJsonExport(); setExportOpen(false); }}>
                  <FileText size={13} />
                  <div>
                    <span className="rv-export-title">JSON</span>
                    <span className="rv-export-desc">Full params + profile data as JSON</span>
                  </div>
                </button>
                <button onClick={() => { handleFullReport(); setExportOpen(false); }} disabled={reportBusy}>
                  <Printer size={13} />
                  <div>
                    <span className="rv-export-title">{reportBusy ? "Generatingâ€¦" : "Full Report PNG"}</span>
                    <span className="rv-export-desc">Screenshot of sounding + all parameters</span>
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
          <button
            className={`rv-btn ${interactiveMode ? "rv-btn-active" : ""}`}
            onClick={() => setInteractiveMode((v) => !v)}
            title={interactiveMode ? "Switch to static plot" : "Switch to interactive Skew-T"}
          >
            <BarChart3 size={14} />
            <span className="rv-btn-label">{interactiveMode ? "Static" : "Interactive"}</span>
          </button>
          {(source === "bufkit" || source === "psu") && (
            <button
              className={`rv-btn ${showAnimator ? "rv-btn-active" : ""}`}
              onClick={() => setShowAnimator((v) => !v)}
              title="Animate through forecast hours"
            >
              <Play size={14} />
              <span className="rv-btn-label">Animate</span>
            </button>
          )}
          {onToggleAutoRefresh && (
            <div className="rv-autorefresh-wrap">
              <button
                className={`rv-btn ${autoRefresh ? "rv-btn-active" : ""}`}
                onClick={onToggleAutoRefresh}
                title={autoRefresh ? `Auto-refreshing every ${refreshInterval / 60000} min` : "Enable auto-refresh"}
              >
                <RefreshCw size={14} className={autoRefresh ? "spin-slow" : ""} />
              </button>
              {autoRefresh && (
                <select
                  className="rv-autorefresh-sel"
                  value={refreshInterval}
                  onChange={(e) => onRefreshIntervalChange(Number(e.target.value))}
                  title="Refresh interval"
                >
                  <option value={60000}>1 min</option>
                  <option value={120000}>2 min</option>
                  <option value={300000}>5 min</option>
                  <option value={600000}>10 min</option>
                  <option value={900000}>15 min</option>
                </select>
              )}
            </div>
          )}
        </div>
        </div>
      </div>
      {/* Sounding Timeline */}
      {onTimelineSelect && (
        <Suspense fallback={null}>
          <SoundingTimeline
            station={selectedStation}
            currentDate={meta?.date}
            onSelectTime={onTimelineSelect}
            source={source}
          />
        </Suspense>
      )}
      {/* Sounding Animator */}
      {showAnimator && (
        <Suspense fallback={<div className="rv-state rv-loading"><Loader2 size={20} className="spin" /><span>Loading animatorâ€¦</span></div>}>
          <SoundingAnimator
            station={lastParams?.station || selectedStation}
            model={lastParams?.model}
            source={lastParams?.source || source}
            date={lastParams?.date}
            theme={theme || "dark"}
            onClose={() => setShowAnimator(false)}
          />
        </Suspense>
      )}
      {/* Plot â€” static PNG or interactive Skew-T */}
      {interactiveMode ? (
        <Suspense fallback={<div className="rv-state rv-loading"><Loader2 size={20} className="spin" /><span>Loading interactive Skew-Tâ€¦</span></div>}>
          <InteractiveSkewT
            profile={result.profile}
            sbParcel={result.sbParcel}
            mlParcel={result.mlParcel}
            params={params}
            theme={theme || "dark"}
          />
          <InteractiveHodograph
            profile={result.profile}
            params={params}
            theme={theme || "dark"}
          />
          {/* Piecewise CAPE strip */}
          {params.piecewiseCape && params.piecewiseCape.length > 0 && (
            <div className="rv-piecewise-cape">
              <div className="rv-pw-title">Piecewise CAPE (50 hPa layers)</div>
              <div className="rv-pw-strip">
                {params.piecewiseCape.map((layer, i) => {
                  const maxCape = Math.max(...params.piecewiseCape.map(l => l.cape), 1);
                  const barW = Math.max(2, (layer.cape / maxCape) * 100);
                  const alpha = Math.min(1, layer.cape / maxCape * 0.8 + 0.2);
                  return (
                    <div key={i} className="rv-pw-row" title={`${layer.p_bot}\u2013${layer.p_top} hPa: CAPE ${layer.cape} J/kg`}>
                      <span className="rv-pw-label">{layer.p_bot}</span>
                      <div className="rv-pw-bar-bg">
                        <div className="rv-pw-bar" style={{ width: `${barW}%`, opacity: alpha }} />
                      </div>
                      <span className="rv-pw-val">{Math.round(layer.cape)}</span>
                    </div>
                  );
                })}
              </div>
            </div>
          )}
        </Suspense>
      ) : (
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
      )}

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
              {/* â”€â”€ Header banner â”€â”€ */}
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

              {/* â”€â”€ Threat callout â”€â”€ */}
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

              {/* â”€â”€ Analysis sections â”€â”€ */}
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
