export function generateSoundingSummary(params, meta) {
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
  const bwd1   = n(params.bwd1km);
  const bwd6   = n(params.bwd6km);
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
    if (lr03 != null) vals.push({ k: "LR 0-3", v: `${fmt(lr03, 1)} Â°C/km` });
    if (lr36 != null) vals.push({ k: "LR 3-6", v: `${fmt(lr36, 1)} Â°C/km` });
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

  // ── 9) Composite indices ──
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
