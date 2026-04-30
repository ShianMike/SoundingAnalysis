import { describe, it, expect } from "vitest";
import { bearingDeg, haversineKm, clampNum, frameToHourOffset, WIND_MAX_HOUR_OFFSET } from "./mapGeometry";

describe("clampNum", () => {
  it("clamps below the minimum", () => {
    expect(clampNum(-5, 0, 10)).toBe(0);
  });
  it("clamps above the maximum", () => {
    expect(clampNum(15, 0, 10)).toBe(10);
  });
  it("passes through values inside the range", () => {
    expect(clampNum(7, 0, 10)).toBe(7);
  });
  it("treats min and max as inclusive", () => {
    expect(clampNum(0, 0, 10)).toBe(0);
    expect(clampNum(10, 0, 10)).toBe(10);
  });
});

describe("bearingDeg", () => {
  it("returns ~0° for due north", () => {
    expect(bearingDeg(0, 0, 1, 0)).toBeCloseTo(0, 5);
  });
  it("returns ~90° for due east", () => {
    expect(bearingDeg(0, 0, 0, 1)).toBeCloseTo(90, 5);
  });
  it("returns ~180° for due south", () => {
    expect(bearingDeg(1, 0, 0, 0)).toBeCloseTo(180, 5);
  });
  it("returns ~270° for due west", () => {
    expect(bearingDeg(0, 1, 0, 0)).toBeCloseTo(270, 5);
  });
  it("normalizes the result into [0, 360)", () => {
    const b = bearingDeg(0, 0, -1, -0.1);
    expect(b).toBeGreaterThanOrEqual(0);
    expect(b).toBeLessThan(360);
  });
});

describe("haversineKm", () => {
  it("returns 0 for identical points", () => {
    expect(haversineKm(35, -97, 35, -97)).toBeCloseTo(0, 6);
  });
  it("returns ~111 km per degree of latitude", () => {
    // 1° of latitude is about 111.2 km.
    const d = haversineKm(0, 0, 1, 0);
    expect(d).toBeGreaterThan(110);
    expect(d).toBeLessThan(112);
  });
  it("Norman, OK -> Oklahoma City is roughly 30 km", () => {
    // OUN ~ (35.22, -97.44), OKC ~ (35.47, -97.52)
    const d = haversineKm(35.22, -97.44, 35.47, -97.52);
    expect(d).toBeGreaterThan(20);
    expect(d).toBeLessThan(40);
  });
});

describe("frameToHourOffset", () => {
  it("returns 0 when there is only one frame", () => {
    expect(frameToHourOffset(0, 1)).toBe(0);
  });
  it("maps the first frame to 0", () => {
    expect(frameToHourOffset(0, 13)).toBe(0);
  });
  it("maps the last frame to maxOffset", () => {
    expect(frameToHourOffset(12, 13, 12)).toBe(12);
  });
  it("scales linearly between 0 and maxOffset", () => {
    expect(frameToHourOffset(6, 13, 12)).toBe(6);
  });
  it("clamps frameIdx into the valid range", () => {
    expect(frameToHourOffset(-5, 13, 12)).toBe(0);
    expect(frameToHourOffset(99,  13, 12)).toBe(12);
  });
  it("default maxOffset is WIND_MAX_HOUR_OFFSET", () => {
    expect(frameToHourOffset(WIND_MAX_HOUR_OFFSET, WIND_MAX_HOUR_OFFSET + 1)).toBe(WIND_MAX_HOUR_OFFSET);
  });
});
