/**
 * Canvas-based PNG export helpers.
 *
 * Pure DOM helpers — no React, no Zustand. Each function builds a canvas,
 * draws the requested artifact, and triggers a browser download.
 */

/**
 * Render a Severe-Weather Risk Scan summary table to a PNG and trigger
 * download. Mirrors the on-screen `<RiskTable>` layout.
 *
 * @param {Object} riskData - { date, model?, fhour?, initTime?, stations: [...] }
 */
export function exportRiskTablePng(riskData) {
  if (!riskData || !Array.isArray(riskData.stations) || riskData.stations.length === 0) return;

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
    subtitle = `${riskData.model} F${riskData.fhour} · Init ${riskData.initTime} · Valid ${riskData.date}`;
  }
  ctx.fillText(subtitle, pad, pad + 46);
  const countLabel = `${stations.length} stations`;
  ctx.fillText(countLabel, totalW - pad - ctx.measureText(countLabel).width, pad + 46);

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
    if (i % 2 === 0) {
      ctx.fillStyle = "rgba(255,255,255,0.02)";
      ctx.fillRect(pad, y, totalW - pad * 2, rowH);
    }
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

  ctx.strokeStyle = "rgba(255,255,255,0.08)";
  ctx.lineWidth = 1;
  ctx.strokeRect(pad, tableY, totalW - pad * 2, rowH * (stations.length + 1));

  const link = document.createElement("a");
  const tag = riskData.model ? `${riskData.model}_F${riskData.fhour}` : "observed";
  link.download = `risk_scan_${tag}_${riskData.date.replace(/[:\s]/g, "")}.png`;
  link.href = canvas.toDataURL("image/png");
  link.click();
}

/**
 * Render a "Full Report" PNG combining the sounding plot image with a panel
 * of derived parameters (CAPE / SRH / STP / etc.). Triggers download.
 *
 * @param {Object} args
 * @param {string} args.image - Base64-encoded PNG of the sounding plot (no data URI prefix).
 * @param {Object} args.params - Derived parameters (sbCape, mlCape, … as on `result.params`).
 * @param {Object} args.meta - { station, stationName, source, date }.
 */
export async function exportSoundingReportPng({ image, params, meta }) {
  // Load the sounding image
  const img = new Image();
  img.src = `data:image/png;base64,${image}`;
  await new Promise((resolve, reject) => {
    img.onload = resolve;
    img.onerror = reject;
  });

  // Build parameter summary rows
  const rows = [
    [`${meta.station || meta.source.toUpperCase()}  —  ${meta.stationName || ""}  —  ${meta.date}`, ""],
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
    ["SRH 0-1 km", `${params.srh1km ?? "---"} m²/s²`], ["SRH 0-3 km", `${params.srh3km ?? "---"} m²/s²`],
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

  // Sounding image
  ctx.drawImage(img, 0, 0);

  // Parameter panel
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
}
