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
} from "lucide-react";
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

/* ── Sounding text summary generator ───────────────────────── */
function generateSoundingSummary(params, meta) {
  const lines = [];
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

  // ── Header / Opening ──
  const stName = meta?.stationName || meta?.station || "this location";
  const dateStr = meta?.date || "";
  lines.push(`Analysis for ${stName}${dateStr ? ` valid ${dateStr}` : ""}:`);

  // ── 1) Thermodynamic profile overview ──
  if (bestCape < 50) {
    lines.push("The thermodynamic profile is stable with virtually no buoyancy (CAPE < 50 J/kg). Thunderstorm development is not supported by this environment. Any convection would require significant mesoscale or synoptic forcing and would likely remain shallow and non-severe.");
  } else if (bestCape < 250) {
    lines.push(`Marginal instability is present with MUCAPE of ${fmt(muCape)} J/kg. Buoyancy is insufficient for robust updrafts, though isolated weak convection could occur along well-defined convergence zones or terrain-enhanced lifting. Severe weather is unlikely.`);
  } else if (bestCape < 1000) {
    lines.push(`Moderate instability exists (MUCAPE ${fmt(muCape)} J/kg, MLCAPE ${fmt(mlCape)} J/kg${ecape != null ? `, ECAPE ${fmt(ecape)}` : ""}). This level of buoyancy can support organized thunderstorms with modest updraft strength, particularly where mesoscale forcing provides reliable initiation.`);
  } else if (bestCape < 2500) {
    lines.push(`Substantial instability is present — MUCAPE ${fmt(muCape)} J/kg with MLCAPE ${fmt(mlCape)} J/kg${ecape != null ? ` and ECAPE of ${fmt(ecape)} J/kg` : ""}. This provides ample energy for vigorous updrafts, and severe weather is likely if storms develop within a supportive kinematic environment.`);
  } else if (bestCape < 4000) {
    lines.push(`Large instability is present with MUCAPE of ${fmt(muCape)} J/kg and MLCAPE of ${fmt(mlCape)} J/kg${ecape != null ? ` (ECAPE ${fmt(ecape)} J/kg)` : ""}. This is a high-end thermodynamic environment capable of producing strong to violent updrafts, very large hail, and extreme rainfall rates.`);
  } else {
    lines.push(`Extreme instability is evident — MUCAPE ${fmt(muCape)} J/kg, MLCAPE ${fmt(mlCape)} J/kg${ecape != null ? `, ECAPE ${fmt(ecape)} J/kg` : ""}. This is a rare, upper-tier thermodynamic environment where updrafts may exceed 50 m/s, supporting giant hail (≥4 in.), violent tornadoes if shear is sufficient, and flash-flood-producing rainfall.`);
  }

  // ── 2) Low-level CAPE / updraft acceleration ──
  if (cape3 != null && bestCape >= 250) {
    if (cape3 >= 100) {
      lines.push(`The 0–3 km CAPE of ${fmt(cape3)} J/kg is notably high, indicating strong low-level updraft acceleration. This enhances stretching of low-level vertical vorticity and is a key ingredient for tornado intensity — dynamically stretching near-surface rotation into tighter, faster-spinning vortices.`);
    } else if (cape3 >= 50) {
      lines.push(`The 0–3 km CAPE of ${fmt(cape3)} J/kg is moderate, providing meaningful low-level updraft acceleration that supports efficient vortex stretching for tornado development.`);
    }
  }

  // ── 3) Cap / CIN analysis ──
  if (mlCin != null) {
    if (mlCin > -10) {
      lines.push("Convective inhibition is negligible (MLCIN near zero). Storms could fire readily with minimal forcing — this is an uncapped environment where widespread initiation is possible, potentially resulting in a quick transition from clear skies to active convection.");
    } else if (mlCin > -50) {
      lines.push(`A weak cap is present (MLCIN ${fmt(mlCin)} J/kg). Convective initiation should occur fairly easily along boundaries, outflow, or modest terrain features. The weak inhibition may limit the explosive character of initial development but allows for broad storm coverage.`);
    } else if (mlCin > -150) {
      lines.push(`A moderate convective cap is in place (MLCIN ${fmt(mlCin)} J/kg). Initiation will require focused mesoscale forcing — synoptic fronts, drylines, outflow boundaries, or orographic lift. Storms that do breach the cap may be explosive due to the stored instability beneath, and discrete supercells become more likely than cluster modes.`);
    } else {
      lines.push(`A strong cap exists (MLCIN ${fmt(mlCin)} J/kg) significantly inhibiting deep convection. Only intense forcing (vigorous frontal lift, strong dryline convergence, or elevated instability mechanisms) is likely to break through this inversion. If a storm does initiate, the sudden release of pent-up energy could produce a violently explosive updraft.`);
    }
  }

  // ── 4) Cloud base height / LCL ──
  if (mlLcl != null && bestCape >= 250) {
    if (mlLcl < 800) {
      lines.push(`Cloud bases are very low (MLLCL ${fmt(mlLcl)} m AGL), which is strongly favorable for tornadogenesis. The short distance between the surface and cloud base minimizes the sub-cloud layer where rain-cooled downdraft air could disrupt the near-surface circulation, allowing tight low-level rotation to connect efficiently to the mesocyclone aloft.`);
    } else if (mlLcl < 1200) {
      lines.push(`Cloud bases are relatively low (MLLCL ${fmt(mlLcl)} m AGL), supporting tornado development. The sub-cloud layer is shallow enough to maintain coherent low-level rotation, though tornadoes may be somewhat less likely than with the very lowest LCL heights.`);
    } else if (mlLcl < 1800) {
      lines.push(`Cloud bases are moderately elevated (MLLCL ${fmt(mlLcl)} m AGL). While this reduces tornado potential, the deeper sub-cloud layer can enhance downdraft evaporation and outflow wind damage. Supercells in this regime tend to produce more hail and damaging straight-line winds than tornadoes.`);
    } else {
      lines.push(`Cloud bases are high (MLLCL ${fmt(mlLcl)} m AGL), which strongly limits tornado potential. However, the deep dry sub-cloud layer promotes vigorous downdraft development through evaporative cooling, increasing the risk of damaging outflow winds and microbursts. Hail size may also be enhanced due to reduced melting in the dry layer.`);
    }
  }

  // ── 5) Lapse rates ──
  if ((lr03 != null || lr36 != null) && bestCape >= 100) {
    const parts = [];
    if (lr03 != null) {
      if (lr03 >= 9) parts.push(`exceptionally steep 0–3 km lapse rates (${fmt(lr03, 1)} °C/km) approaching dry-adiabatic, indicating a superadiabatic or near-superadiabatic boundary layer with intense low-level buoyancy`);
      else if (lr03 >= 8) parts.push(`very steep 0–3 km lapse rates (${fmt(lr03, 1)} °C/km) indicating a well-mixed, nearly adiabatic boundary layer favorable for strong low-level stretching`);
      else if (lr03 >= 7) parts.push(`moderately steep 0–3 km lapse rates (${fmt(lr03, 1)} °C/km)`);
      else if (lr03 < 5.5) parts.push(`relatively weak 0–3 km lapse rates (${fmt(lr03, 1)} °C/km) suggesting a stable boundary layer that may resist surface-based storm development`);
    }
    if (lr36 != null) {
      if (lr36 >= 8.5) parts.push(`extreme mid-level lapse rates (${fmt(lr36, 1)} °C/km) contributing to exceptionally deep CAPE and explosive updraft growth through the troposphere`);
      else if (lr36 >= 7.5) parts.push(`steep mid-level lapse rates (${fmt(lr36, 1)} °C/km) enhancing CAPE depth and supporting vigorous updraft acceleration`);
      else if (lr36 >= 6.5) parts.push(`moderate mid-level lapse rates (${fmt(lr36, 1)} °C/km)`);
    }
    if (parts.length > 0) {
      lines.push(`The temperature profile shows ${parts.join(", and ")}.`);
    }
  }

  // ── 6) Moisture profile ──
  const moistParts = [];
  if (rh01 != null && bestCape >= 100) {
    if (rh01 >= 85) moistParts.push(`very moist low levels (0–1 km RH ${fmt(rh01)}%)`);
    else if (rh01 < 60) moistParts.push(`relatively dry low levels (0–1 km RH ${fmt(rh01)}%)`);
  }
  if (rh36 != null && bestCape >= 100) {
    if (rh36 < 30) moistParts.push(`a very dry mid-level layer (3–6 km RH ${fmt(rh36)}%) which enhances evaporative downdraft production and DCAPE`);
    else if (rh36 < 50) moistParts.push(`moderately dry mid-levels (3–6 km RH ${fmt(rh36)}%)`);
  }
  if (pwat != null) {
    if (pwat >= 50) moistParts.push(`extremely high precipitable water (${fmt(pwat)} mm) signaling a saturated column with major flash flood potential`);
    else if (pwat >= 40) moistParts.push(`high precipitable water (${fmt(pwat)} mm) increasing flash flood risk from any training or slow-moving convection`);
    else if (pwat >= 30) moistParts.push(`moderate precipitable water (${fmt(pwat)} mm)`);
  }
  if (moistParts.length > 0) {
    lines.push(`The moisture profile features ${moistParts.join("; ")}.`);
  }

  // ── 7) Deep-layer shear ──
  if (bwd6 != null && bestCape >= 100) {
    if (bwd6 >= 70) {
      lines.push(`Deep-layer shear is exceptionally strong (0–6 km BWD ${fmt(bwd6)} kt). This extreme kinematic environment overwhelmingly favors discrete, long-track supercells with persistent mesocyclones. Storm longevity and severity are maximized.`);
    } else if (bwd6 >= 50) {
      lines.push(`Deep-layer shear is very strong (0–6 km BWD ${fmt(bwd6)} kt), robustly supporting supercellular convection. Discrete supercells are the expected storm mode, with the shear vector promoting strong updraft-downdraft separation and long-lived mesocyclones.`);
    } else if (bwd6 >= 35) {
      lines.push(`Moderate-to-strong deep-layer shear exists (0–6 km BWD ${fmt(bwd6)} kt). This supports organized convection ranging from splitting supercells to organized multicellular systems, depending on storm-relative wind profiles.${ebwd != null ? ` Effective BWD of ${fmt(ebwd)} kt provides a physically-grounded shear estimate.` : ""}`);
    } else if (bwd6 >= 20) {
      lines.push(`Moderate deep-layer shear (0–6 km BWD ${fmt(bwd6)} kt) provides some storm organization. Multicell clusters with embedded supercell structures are possible, though discrete supercells are less likely. Convective mode may depend on mesoscale boundary interactions.`);
    } else {
      lines.push(`Weak deep-layer shear (0–6 km BWD ${fmt(bwd6)} kt) favors disorganized convection — pulse storms or loosely organized multicells. Without significant shear, storms will struggle to maintain persistent updraft rotation.`);
    }
  }

  // ── 8) Low-level shear & SRH — convective mode / tornado assessment ──
  if (bestCape >= 250 && (srh1 != null || bwd1 != null)) {
    const llParts = [];
    if (bwd1 != null) {
      if (bwd1 >= 30) llParts.push(`extreme 0–1 km shear (${fmt(bwd1)} kt)`);
      else if (bwd1 >= 20) llParts.push(`strong 0–1 km shear (${fmt(bwd1)} kt)`);
      else if (bwd1 >= 15) llParts.push(`moderate 0–1 km shear (${fmt(bwd1)} kt)`);
    }
    if (srh1 != null) {
      if (srh1 >= 400) llParts.push(`extreme 0–1 km SRH (${fmt(srh1)} m²/s²) — a rare value strongly correlated with violent (EF4-EF5) tornadoes`);
      else if (srh1 >= 200) llParts.push(`very significant 0–1 km SRH (${fmt(srh1)} m²/s²) highly supportive of strong to violent tornadoes`);
      else if (srh1 >= 100) llParts.push(`notable 0–1 km SRH (${fmt(srh1)} m²/s²) favoring mesocyclone development and tornado potential`);
      else if (srh1 >= 50) llParts.push(`modest 0–1 km SRH (${fmt(srh1)} m²/s²)`);
    }
    if (esrh != null && esrh > (srh1 ?? 0) * 1.2) {
      llParts.push(`effective SRH of ${fmt(esrh)} m²/s² (deeper effective inflow layer contributing additional helicity)`);
    }
    if (llParts.length > 0) {
      lines.push(`Low-level kinematics show ${llParts.join(", with ")}.`);
    }
  }

  // ── 9) Composite indices — detailed threat assessment ──
  const compLines = [];

  if (stp != null && bestCape >= 250) {
    if (stp >= 8) {
      compLines.push(`The STP of ${fmt(stp, 1)} is in the top percentile of significant tornado environments — this is a high-confidence setup for strong to violent (EF3+) tornadoes with any sustained supercell.`);
    } else if (stp >= 4) {
      compLines.push(`The STP of ${fmt(stp, 1)} is a high-end value, indicating a well-above-average environment for significant (EF2+) tornadoes. Multiple tornado-producing supercells are plausible.`);
    } else if (stp >= 1) {
      compLines.push(`The STP of ${fmt(stp, 1)} exceeds the significant tornado threshold, indicating the kinematic and thermodynamic ingredients are synergistically aligned for EF2+ tornadoes with sustained supercells.`);
    } else if (stp >= 0.5) {
      compLines.push(`The STP of ${fmt(stp, 1)} is below the canonical significant-tornado threshold of 1.0 but still suggests a non-trivial tornado risk, particularly for brief or weak (EF0-EF1) tornadoes.`);
    }
  }

  if (scp != null && bestCape >= 250) {
    if (scp >= 10) {
      compLines.push(`The SCP of ${fmt(scp, 1)} is extreme, strongly supporting long-lived discrete supercells with persistent mesocyclones and a high probability of significant severe weather.`);
    } else if (scp >= 4) {
      compLines.push(`The SCP of ${fmt(scp, 1)} strongly favors supercellular convection, with persistent mesocyclones expected in any sustained storms.`);
    } else if (scp >= 1) {
      compLines.push(`The SCP of ${fmt(scp, 1)} supports supercell development with moderate mesocyclone potential.`);
    }
  }

  if (ship != null && bestCape >= 500) {
    if (ship >= 2.5) {
      compLines.push(`The SHIP of ${fmt(ship, 1)} is very high, indicating a robust environment for significant hail (≥2 inches). Supercell updrafts will be strong enough to support large-diameter ice growth through repeated recycling above the freezing level.`);
    } else if (ship >= 1.0) {
      compLines.push(`The SHIP of ${fmt(ship, 1)} exceeds the significant-hail threshold, suggesting an environment capable of producing large hail (≥1 inch) with supercellular storms.`);
    } else if (ship >= 0.5) {
      compLines.push(`The SHIP of ${fmt(ship, 1)} indicates marginal significant hail potential — large hail is possible but may be limited in size.`);
    }
  }

  if (dcp != null && dcp >= 2) {
    if (dcp >= 6) {
      compLines.push(`The DCP of ${fmt(dcp, 1)} is extreme, strongly signaling a derecho or widespread damaging wind event. Bow echoes and book-end vortices are favored in this shear and instability regime.`);
    } else if (dcp >= 4) {
      compLines.push(`The DCP of ${fmt(dcp, 1)} is elevated, supporting organized long-lived wind events with potential for widespread wind damage.`);
    } else {
      compLines.push(`The DCP of ${fmt(dcp, 1)} indicates some potential for organized damaging wind events along squall lines.`);
    }
  }

  if (compLines.length > 0) {
    lines.push(compLines.join(" "));
  } else if (bestCape >= 500 && bwd6 >= 20) {
    lines.push("Composite severe parameters (STP, SCP, SHIP) remain below critical thresholds, suggesting that while organized convection is possible, the probability of significant severe weather (violent tornadoes, giant hail, or derecho-scale winds) is relatively low.");
  }

  // ── 10) DCAPE / Downdraft threat ──
  if (dcape != null && dcape >= 600 && bestCape >= 100) {
    if (dcape >= 1500) {
      lines.push(`DCAPE of ${fmt(dcape)} J/kg is extreme, indicating very strong potential for damaging outflow winds at the surface. Microbursts and macrobursts are likely, and wet or dry downbursts could produce wind gusts exceeding 80 kt.`);
    } else if (dcape >= 1000) {
      lines.push(`DCAPE of ${fmt(dcape)} J/kg is strong, supporting damaging outflow winds. Isolated downbursts with gusts of 50–70 kt are probable with any organized convection.`);
    } else {
      lines.push(`DCAPE of ${fmt(dcape)} J/kg provides moderate downdraft energy, sufficient for gusty outflow winds (40–55 kt) particularly where mid-level dry air is entrained into storm downdrafts.`);
    }
  }

  // ── 11) Hail environment ──
  if (bestCape >= 500 && frzLvl != null && wbo != null) {
    if (frzLvl > 4000 && wbo > 3000 && ship >= 0.5) {
      lines.push(`The freezing level (${fmt(frzLvl)} m AGL) and wet-bulb zero height (${fmt(wbo)} m AGL) are both elevated, allowing hailstones a long fall through warm air. Despite this partial melting, the strong updrafts indicated by the CAPE profile can still support large hail aloft — stones reaching the surface may be somewhat smaller but remain capable of damage.`);
    } else if (frzLvl < 3000 && wbo < 2500 && ship >= 0.5) {
      lines.push(`The low freezing level (${fmt(frzLvl)} m AGL) and wet-bulb zero (${fmt(wbo)} m AGL) mean hailstones have a shorter fall through above-freezing air, reducing melting. This favors larger hail reaching the surface relative to the updraft strength.`);
    }
  }

  // ── 12) Convective mode forecast ──
  if (bestCape >= 250) {
    let mode = "";
    if (bwd6 >= 40 && (srh1 >= 100 || scp >= 1)) {
      mode = "discrete supercells";
    } else if (bwd6 >= 35 && dcp >= 3) {
      mode = "supercells transitioning to a bowing line segment";
    } else if (bwd6 >= 25 && dcp >= 2) {
      mode = "organized multicells with embedded bow echoes";
    } else if (bwd6 >= 25) {
      mode = "multicells with possible supercell structures";
    } else if (bestCape >= 1000) {
      mode = "pulse and weakly organized multicellular storms";
    } else {
      mode = "isolated pulse convection";
    }
    lines.push(`Expected convective mode: ${mode}.`);
  }

  // ── 13) Overall threat summary ──
  const threats = [];
  if (stp != null && stp >= 1) {
    threats.push(stp >= 4 ? "strong to violent tornadoes (EF2+)" : "significant tornadoes");
  } else if (stp != null && stp >= 0.3 && mlLcl < 1500) {
    threats.push("brief/weak tornadoes");
  }
  if (ship != null && ship >= 1.5) threats.push("significant hail (≥2 inches)");
  else if (ship != null && ship >= 0.5) threats.push("large hail");
  else if (bestCape >= 2500 && bwd6 >= 30) threats.push("large hail");
  if (dcp != null && dcp >= 4) threats.push("widespread damaging winds / derecho");
  else if (dcp != null && dcp >= 2) threats.push("damaging straight-line winds");
  else if (dcape != null && dcape >= 1000 && bwd6 >= 20) threats.push("damaging outflow gusts");
  if (pwat != null && pwat >= 40) threats.push("flash flooding");

  if (threats.length > 0) {
    lines.push(`Primary hazards: ${threats.join("; ")}.`);
  } else if (bestCape >= 500) {
    lines.push("No parameters reach high-confidence thresholds for significant severe weather. General thunderstorm hazards (gusty winds, small hail, lightning, locally heavy rain) remain possible with any convection that develops.");
  } else if (bestCape < 100) {
    lines.push("No convective hazards are anticipated from this thermodynamic profile.");
  }

  return lines.join("\n\n");
}

function ParamCard({ label, value, unit, color, desc }) {
  const alertCls = getAlertClass(label, value);
  return (
    <div className={`param-card ${alertCls}`} title="">
      <span className="param-label">{label}</span>
      <span className="param-value" style={color ? { color } : {}}>
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
  const [showClimo, setShowClimo] = useState(false);
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
          <ParamCard label="DCAPE" value={params.dcape} unit="J/kg" desc="Downdraft CAPE — energy available for downdrafts. Higher values (>800) indicate strong outflow winds and potential for damaging gusts." />
          <ParamCard label="ECAPE" value={params.ecape} unit="J/kg" color="#06b6d4" desc="Entraining CAPE — CAPE adjusted for entrainment of dry environmental air (Peters et al. 2023). More physically realistic than standard CAPE." />
          <ParamCard label="3CAPE" value={params.cape3km} unit="J/kg" color="#fb923c" desc="0–3 km CAPE — buoyant energy in the lowest 3 km (MU parcel). Higher values indicate stronger low-level accelerations; key for tornado intensity." />
          <ParamCard label="6CAPE" value={params.cape6km} unit="J/kg" color="#facc15" desc="0–6 km CAPE — buoyant energy in the lowest 6 km (MU parcel). Indicates how quickly an updraft develops in the mid-levels." />
          <ParamCard label="DCIN" value={params.dcin} unit="J/kg" color="#818cf8" desc="Downdraft CIN — inhibition of downdrafts reaching the surface. More negative = stronger capping of downdrafts. Near 0 means downdrafts easily penetrate to the surface." />
          <ParamCard label="NCAPE" value={params.ncape} unit="J/kg/m" color="#38bdf8" desc="Normalized CAPE — MUCAPE divided by the LFC-to-EL depth. Measures buoyancy intensity per unit depth. >0.3 is very buoyant; indicates narrow, intense updrafts." />
          <ParamCard label="STP" value={params.stp} unit="" color="#60a5fa" desc="Significant Tornado Parameter — composite index combining CAPE, SRH, shear, and LCL. Values ≥1 suggest significant (EF2+) tornado environment. Higher = more favorable." />
          <ParamCard label="SCP" value={params.scp} unit="" color="#f59e0b" desc="Supercell Composite Parameter — combines CAPE, deep shear, and SRH. Values ≥1 support supercells; >4 strongly favors discrete supercells." />
          <ParamCard label="SHIP" value={params.ship} unit="" color="#10b981" desc="Significant Hail Parameter — composite for significant hail (≥2 in.). Values ≥1 indicate potential; >1.5 strongly favors significant hail." />
          <ParamCard label="DCP" value={params.dcp} unit="" color="#a78bfa" desc="Derecho Composite Parameter — combines DCAPE, MUCAPE, shear, and mean wind. Values ≥2 suggest potential for long-lived wind events (derechos)." />
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
          <ParamCard label="RH 0-1 km" value={params.rh01} unit="%" desc="Relative Humidity 0–1 km — low-level moisture. Higher values (>80%) favor tornado development by reducing evaporative cooling of rain." />
          <ParamCard label="RH 1-3 km" value={params.rh13} unit="%" desc="Relative Humidity 1–3 km — mid-low moisture. Dry layers here (<50%) enhance DCAPE and outflow potential." />
          <ParamCard label="RH 3-6 km" value={params.rh36} unit="%" desc="Relative Humidity 3–6 km — mid-level moisture. Dry air here entrains into storms, increasing evaporative cooling and downdraft strength." />
        </ParamSection>

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
        </ParamSection>
      </div>

      {/* Text Summary */}
      <div className="rv-summary-section">
        <button className="rv-summary-toggle" onClick={() => setShowSummary((s) => !s)}>
          <FileText size={14} />
          <span>Sounding Text Summary</span>
          <ChevronDown size={12} className={showSummary ? "rv-chev-open" : ""} />
        </button>
        {showSummary && (
          <div className="rv-summary-body">
            {generateSoundingSummary(params, meta).split("\n\n").map((p, i) => (
              <p key={i}>{p}</p>
            ))}
          </div>
        )}
      </div>

      {/* Climatology Percentiles */}
      {params.percentiles && Object.keys(params.percentiles).length > 0 && (
        <div className="rv-summary-section">
          <button className="rv-summary-toggle" onClick={() => setShowClimo((s) => !s)}>
            <BarChart3 size={14} />
            <span>Climatology Percentiles</span>
            <ChevronDown size={12} className={showClimo ? "rv-chev-open" : ""} />
          </button>
          {showClimo && (
            <div className="rv-climo-body">
              <p className="rv-climo-desc">
                Percentile rank vs. severe-weather proximity sounding climatology (SPC studies).
                Higher = more extreme environment.
              </p>
              <div className="rv-climo-bars">
                {[
                  ["SB CAPE", "sbCape"],
                  ["ML CAPE", "mlCape"],
                  ["MU CAPE", "muCape"],
                  ["ECAPE", "ecape"],
                  ["ML CIN", "mlCin"],
                  ["0-6 BWD", "bwd6km"],
                  ["0-1 BWD", "bwd1km"],
                  ["0-1 SRH", "srh1km"],
                  ["0-3 SRH", "srh3km"],
                  ["Eff SRH", "esrh"],
                  ["Eff BWD", "ebwd"],
                  ["STP", "stp"],
                  ["SCP", "scp"],
                  ["SHIP", "ship"],
                  ["DCP", "dcp"],
                  ["DCAPE", "dcape"],
                  ["LR 0-3", "lr03"],
                  ["LR 3-6", "lr36"],
                  ["PWAT", "pwat"],
                ].filter(([, k]) => params.percentiles[k] != null).map(([label, key]) => {
                  const pct = params.percentiles[key];
                  const color = pct >= 95 ? "#ef4444" : pct >= 75 ? "#f59e0b" : pct >= 50 ? "#22c55e" : "#6b7280";
                  return (
                    <div key={key} className="rv-climo-item">
                      <span className="rv-climo-label">{label}</span>
                      <div className="rv-climo-track">
                        <div className="rv-climo-fill" style={{ width: `${pct}%`, background: color }} />
                      </div>
                      <span className="rv-climo-pct">{pct}%</span>
                    </div>
                  );
                })}
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
