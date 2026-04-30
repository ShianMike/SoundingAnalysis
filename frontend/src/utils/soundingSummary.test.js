import { describe, it, expect } from "vitest";
import { generateSoundingSummary } from "./soundingSummary";

const META = { stationName: "Norman, OK", station: "OUN", date: "2025-04-30 18:00Z" };

/** Build a `params` object stub with the given overrides. */
function p(overrides = {}) {
  return {
    sbCape: 0, mlCape: 0, muCape: 0, ecape: null, cape3km: null,
    sbCin: 0, mlCin: 0, mlLclM: null,
    dcape: null, lr03: null, lr36: null,
    bwd500m: null, bwd1km: null, bwd6km: null,
    srh500m: null, srh1km: null, srh3km: null,
    esrh: null, ebwd: null,
    stp: null, scp: null, ship: null, dcp: null, brn: null,
    ehi01: null, ehi03: null, vgp: null, criticalAngle: null,
    pwat: null, rh01: null, rh36: null,
    frzLevel: null, wbo: null, wmsi: null, mdpi: null,
    precipType: "Rain",
    ...overrides,
  };
}

describe("generateSoundingSummary — top-level shape", () => {
  it("returns the expected keys", () => {
    const r = generateSoundingSummary(p(), META);
    expect(r).toHaveProperty("station");
    expect(r).toHaveProperty("date");
    expect(r).toHaveProperty("threatLevel");
    expect(r).toHaveProperty("threatColor");
    expect(r).toHaveProperty("threats");
    expect(r).toHaveProperty("sections");
    expect(Array.isArray(r.threats)).toBe(true);
    expect(Array.isArray(r.sections)).toBe(true);
  });

  it("uses meta.stationName / station as the displayed station", () => {
    const r = generateSoundingSummary(p(), META);
    expect(r.station).toBe("Norman, OK");
    const r2 = generateSoundingSummary(p(), { station: "OUN", date: "" });
    expect(r2.station).toBe("OUN");
  });

  it("falls back to 'this location' when meta is missing identifiers", () => {
    const r = generateSoundingSummary(p(), {});
    expect(r.station).toBe("this location");
  });
});

describe("generateSoundingSummary — threat level banner", () => {
  it("classifies a stable, no-CAPE profile as 'none'", () => {
    const r = generateSoundingSummary(p(), META);
    expect(r.threatLevel).toBe("none");
    expect(r.threatColor).toBe("#6b7280");
  });

  it("classifies a marginal-CAPE profile as 'MARGINAL'", () => {
    const r = generateSoundingSummary(p({ sbCape: 350, mlCape: 300 }), META);
    expect(r.threatLevel).toBe("MARGINAL");
  });

  it("classifies a low-end severe profile (CAPE+shear) as 'LOW'", () => {
    const r = generateSoundingSummary(p({ sbCape: 800, mlCape: 700, bwd6km: 30 }), META);
    expect(r.threatLevel).toBe("LOW");
  });

  it("classifies a moderate-tornado profile (STP >= 1) as 'MODERATE'", () => {
    const r = generateSoundingSummary(p({
      sbCape: 2500, mlCape: 2200, bwd6km: 45, stp: 2.5,
    }), META);
    expect(r.threatLevel).toBe("MODERATE");
  });

  it("classifies an extreme STP profile as 'HIGH'", () => {
    const r = generateSoundingSummary(p({
      sbCape: 4500, mlCape: 4000, bwd6km: 60, stp: 6,
    }), META);
    expect(r.threatLevel).toBe("HIGH");
  });
});

describe("generateSoundingSummary — threats list", () => {
  it("flags 'Significant tornadoes' when STP exceeds 1", () => {
    const r = generateSoundingSummary(p({
      sbCape: 2500, mlCape: 2200, bwd6km: 45, stp: 2,
    }), META);
    expect(r.threats).toContain("Significant tornadoes");
  });

  it("upgrades to 'Strong to violent tornadoes (EF2+)' at STP >= 4", () => {
    const r = generateSoundingSummary(p({
      sbCape: 4000, mlCape: 3500, bwd6km: 55, stp: 5,
    }), META);
    expect(r.threats).toContain("Strong to violent tornadoes (EF2+)");
  });

  it("flags 'Significant hail (≥2 in.)' when SHIP >= 1.5", () => {
    const r = generateSoundingSummary(p({
      sbCape: 3500, mlCape: 3000, bwd6km: 40, ship: 2.0,
    }), META);
    expect(r.threats).toContain("Significant hail (≥2 in.)");
  });

  it("flags 'Widespread damaging winds / derecho' when DCP >= 4", () => {
    const r = generateSoundingSummary(p({
      sbCape: 3000, mlCape: 2500, bwd6km: 40, dcp: 5,
    }), META);
    expect(r.threats).toContain("Widespread damaging winds / derecho");
  });

  it("flags 'Flash flooding' when PWAT >= 40 mm", () => {
    const r = generateSoundingSummary(p({ sbCape: 2000, pwat: 45 }), META);
    expect(r.threats).toContain("Flash flooding");
  });

  it("flags winter precip when precipType is non-Rain", () => {
    const r = generateSoundingSummary(p({ precipType: "Snow" }), META);
    expect(r.threats).toContain("Winter precip: Snow");
  });

  it("does NOT flag winter precip for plain rain", () => {
    const r = generateSoundingSummary(p({ precipType: "Rain" }), META);
    expect(r.threats.find((t) => t.startsWith("Winter precip"))).toBeUndefined();
  });
});

describe("generateSoundingSummary — sections", () => {
  it("always emits a Thermodynamic Overview section", () => {
    const r = generateSoundingSummary(p(), META);
    const thermo = r.sections.find((s) => s.id === "thermo");
    expect(thermo).toBeDefined();
    expect(thermo.title).toBe("Thermodynamic Overview");
    expect(thermo.text).toMatch(/stable|instab|buoyancy/i);
  });

  it("emits a Convective Inhibition section when MLCIN is set", () => {
    const r = generateSoundingSummary(p({ mlCin: -100 }), META);
    expect(r.sections.find((s) => s.id === "cin")).toBeDefined();
  });

  it("emits a Cloud Base section only with sufficient CAPE", () => {
    // Below 250 J/kg CAPE the section is suppressed.
    const r1 = generateSoundingSummary(p({ mlLclM: 800 }), META);
    expect(r1.sections.find((s) => s.id === "lcl")).toBeUndefined();
    // Above 250 J/kg it should appear.
    const r2 = generateSoundingSummary(p({ sbCape: 1000, mlCape: 900, mlLclM: 800 }), META);
    expect(r2.sections.find((s) => s.id === "lcl")).toBeDefined();
  });

  it("emits a Convective Mode section above 250 J/kg CAPE", () => {
    const r = generateSoundingSummary(p({ sbCape: 1000, mlCape: 900, bwd6km: 45, srh1km: 200 }), META);
    const mode = r.sections.find((s) => s.id === "mode");
    expect(mode).toBeDefined();
    expect(mode.text).toMatch(/supercell|multicell|pulse/i);
  });

  it("emits a Composite Indices section with the right tornado tier text", () => {
    const r = generateSoundingSummary(p({
      sbCape: 4000, mlCape: 3500, bwd6km: 55, stp: 6,
    }), META);
    const comp = r.sections.find((s) => s.id === "composite");
    expect(comp).toBeDefined();
    expect(comp.text).toMatch(/significant|violent|tornado/i);
  });
});

describe("generateSoundingSummary — null safety", () => {
  it("does not crash on completely empty params", () => {
    expect(() => generateSoundingSummary({}, {})).not.toThrow();
  });

  it("does not crash on totally null params object", () => {
    expect(() => generateSoundingSummary(p({
      sbCape: null, mlCape: null, muCape: null,
    }), {})).not.toThrow();
  });
});
