import { describe, it, expect } from "vitest";
import { nearestNexrad, NEXRAD_SITES, SOURCE_META, MODEL_META, FHOUR_PRESETS } from "./constants";

describe("nearestNexrad", () => {
  it("returns null for non-finite coordinates", () => {
    expect(nearestNexrad(NaN, -100)).toBeNull();
    expect(nearestNexrad(35, Infinity)).toBeNull();
    expect(nearestNexrad(undefined, -100)).toBeNull();
    expect(nearestNexrad(35, null)).toBeNull();
  });

  it("returns an object with id, idK (K-prefixed), lat, lon", () => {
    const r = nearestNexrad(35.33, -97.28); // ~ TLX (Norman, OK)
    expect(r).not.toBeNull();
    expect(r).toHaveProperty("id");
    expect(r).toHaveProperty("idK");
    expect(r).toHaveProperty("lat");
    expect(r).toHaveProperty("lon");
    expect(r.idK).toBe(`K${r.id}`);
  });

  it("picks TLX for a point sitting on top of TLX (Norman, OK)", () => {
    // TLX: 35.33, -97.28
    const r = nearestNexrad(35.33, -97.28);
    expect(r.id).toBe("TLX");
    expect(r.idK).toBe("KTLX");
  });

  it("picks DTX for southeast Michigan", () => {
    // DTX: 42.70, -83.47 (Detroit/White Lake)
    const r = nearestNexrad(42.5, -83.0);
    expect(r.id).toBe("DTX");
  });

  it("picks one of the SoCal radars for Los Angeles", () => {
    // LA is roughly at 34.05, -118.25. Closest CONUS NEXRAD sites are
    // VTX (Vandenberg-area / Los Angeles), SOX (Santa Ana), or NKX (San Diego).
    const r = nearestNexrad(34.05, -118.25);
    expect(["VTX", "SOX", "NKX"]).toContain(r.id);
  });

  it("returns one of the NEXRAD_SITES entries (no fabricated values)", () => {
    const r = nearestNexrad(40.0, -100.0);
    const matching = NEXRAD_SITES.find(
      ([id, lat, lon]) => id === r.id && lat === r.lat && lon === r.lon,
    );
    expect(matching).toBeDefined();
  });
});

describe("NEXRAD_SITES dataset", () => {
  it("contains a non-empty list of triplets", () => {
    expect(Array.isArray(NEXRAD_SITES)).toBe(true);
    expect(NEXRAD_SITES.length).toBeGreaterThan(100);
    for (const row of NEXRAD_SITES) {
      expect(Array.isArray(row)).toBe(true);
      expect(row).toHaveLength(3);
      const [id, lat, lon] = row;
      expect(typeof id).toBe("string");
      expect(id).toHaveLength(3);
      expect(Number.isFinite(lat)).toBe(true);
      expect(Number.isFinite(lon)).toBe(true);
      expect(lat).toBeGreaterThan(20);
      expect(lat).toBeLessThan(50);
      expect(lon).toBeLessThan(-65);
      expect(lon).toBeGreaterThan(-130);
    }
  });

  it("has unique site ids", () => {
    const ids = NEXRAD_SITES.map((s) => s[0]);
    expect(new Set(ids).size).toBe(ids.length);
  });
});

describe("SOURCE_META", () => {
  it("documents each known source key", () => {
    for (const key of ["obs", "bufkit", "psu"]) {
      expect(SOURCE_META[key]).toBeDefined();
      expect(typeof SOURCE_META[key].label).toBe("string");
      expect(typeof SOURCE_META[key].desc).toBe("string");
    }
  });
});

describe("MODEL_META", () => {
  it("contains the operational mesoscale + global models with sane fhour bounds", () => {
    for (const key of ["hrrr", "rap", "nam", "namnest", "gfs"]) {
      expect(MODEL_META[key]).toBeDefined();
      const meta = MODEL_META[key];
      expect(meta.maxF).toBeGreaterThan(0);
      expect(meta.step).toBeGreaterThan(0);
      expect(meta.lag).toBeGreaterThanOrEqual(0);
      expect(meta.interval).toBeGreaterThan(0);
    }
  });

  it("HRRR maxF is at most 48 hours (operational HRRR limit)", () => {
    expect(MODEL_META.hrrr.maxF).toBeLessThanOrEqual(48);
  });

  it("GFS maxF reaches into medium range (>= 120h)", () => {
    expect(MODEL_META.gfs.maxF).toBeGreaterThanOrEqual(120);
  });
});

describe("FHOUR_PRESETS", () => {
  it("is a sorted ascending list starting at 0", () => {
    expect(FHOUR_PRESETS[0]).toBe(0);
    for (let i = 1; i < FHOUR_PRESETS.length; i++) {
      expect(FHOUR_PRESETS[i]).toBeGreaterThan(FHOUR_PRESETS[i - 1]);
    }
  });
});
