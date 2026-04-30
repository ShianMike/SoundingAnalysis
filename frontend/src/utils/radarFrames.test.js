import { describe, it, expect, beforeEach, afterEach, vi } from "vitest";
import { fmtRadarTime, buildMosaicFrames, buildSingleSiteFrames } from "./radarFrames";

describe("fmtRadarTime", () => {
  it("returns the HH:MMZ form for a known UTC instant", () => {
    // 2025-04-30T14:05:00Z → "14:05Z"
    const ts = Date.UTC(2025, 3, 30, 14, 5, 0) / 1000;
    expect(fmtRadarTime(ts)).toBe("14:05Z");
  });
  it("returns a placeholder for non-finite input", () => {
    expect(fmtRadarTime(NaN)).toBe("--:--Z");
    expect(fmtRadarTime(Infinity)).toBe("--:--Z");
    expect(fmtRadarTime(undefined)).toBe("--:--Z");
  });
});

describe("buildMosaicFrames", () => {
  beforeEach(() => {
    // Freeze "now" to 2025-04-30T14:07:42Z so the 5-minute snap is stable.
    vi.useFakeTimers();
    vi.setSystemTime(new Date("2025-04-30T14:07:42Z"));
  });
  afterEach(() => {
    vi.useRealTimers();
  });

  it("returns the requested number of frames", () => {
    expect(buildMosaicFrames(6)).toHaveLength(6);
    expect(buildMosaicFrames(24)).toHaveLength(24);
  });

  it("returns frames in ascending time order (oldest → newest)", () => {
    const frames = buildMosaicFrames(8);
    for (let i = 1; i < frames.length; i++) {
      expect(frames[i].time).toBeGreaterThan(frames[i - 1].time);
    }
  });

  it("snaps the latest frame to the most-recent 5-minute boundary", () => {
    const frames = buildMosaicFrames(1);
    // 14:07 should snap down to 14:05.
    expect(frames[0].ts).toBe("202504301405");
  });

  it("frames are spaced exactly 5 minutes apart", () => {
    const frames = buildMosaicFrames(4);
    for (let i = 1; i < frames.length; i++) {
      expect(frames[i].time - frames[i - 1].time).toBe(300);
    }
  });

  it("ts is the YYYYMMDDHHmm form (12 digits, no separators)", () => {
    const frames = buildMosaicFrames(2);
    for (const f of frames) {
      expect(f.ts).toMatch(/^\d{12}$/);
    }
  });
});

describe("buildSingleSiteFrames", () => {
  beforeEach(() => {
    vi.useFakeTimers();
    vi.setSystemTime(new Date("2025-04-30T14:07:42Z"));
  });
  afterEach(() => {
    vi.useRealTimers();
  });

  it("returns frames with a fully-qualified IEM archive URL", () => {
    const [first] = buildSingleSiteFrames("KTLX", "N0Q", 1);
    expect(first).toHaveProperty("url");
    expect(first.url).toContain("mesonet.agron.iastate.edu");
    expect(first.url).toContain("KTLX_N0Q_");
    expect(first.url).toMatch(/\.png$/);
  });

  it("encodes the 5-minute-snapped timestamp into the URL", () => {
    const [first] = buildSingleSiteFrames("KTLX", "N0Q", 1);
    expect(first.url).toContain("KTLX_N0Q_202504301405.png");
  });

  it("returns the requested number of frames in ascending order", () => {
    const frames = buildSingleSiteFrames("KOUN", "N0Q", 5);
    expect(frames).toHaveLength(5);
    for (let i = 1; i < frames.length; i++) {
      expect(frames[i].time).toBeGreaterThan(frames[i - 1].time);
    }
  });
});
