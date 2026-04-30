import { describe, it, expect } from "vitest";
import {
  WARNING_SUBCAT_STYLES,
  VALID_WARNING_EVENTS,
  WARNING_PRIORITY,
  firstParam,
  classifyWarning,
  parseEventMotion,
  escapeHtml,
  warningStyle,
} from "./warnings";

/* ─────────────────── classifyWarning ─────────────────── */

function feature(event, params = {}, extra = {}) {
  return { properties: { event, parameters: params, ...extra } };
}

describe("classifyWarning", () => {
  it("returns null for an unknown event", () => {
    expect(classifyWarning(feature("Frost Advisory"))).toBeNull();
  });

  it("returns null when the feature is malformed", () => {
    expect(classifyWarning(null)).toBeNull();
    expect(classifyWarning({})).toBeNull();
    expect(classifyWarning({ properties: null })).toBeNull();
  });

  describe("Tornado Warning sub-categories", () => {
    it("classifies a base tornado warning as TOR_RI (radar-indicated)", () => {
      expect(classifyWarning(feature("Tornado Warning"))).toBe("TOR_RI");
    });

    it("classifies an Observed tornado warning as TOR_OBS", () => {
      const f = feature("Tornado Warning", { tornadoDetection: ["Observed"] });
      expect(classifyWarning(f)).toBe("TOR_OBS");
    });

    it("upgrades to TOR_OBS via headline keywords", () => {
      const f = feature("Tornado Warning", {}, { headline: "CONFIRMED TORNADO ON THE GROUND" });
      expect(classifyWarning(f)).toBe("TOR_OBS");
    });

    it("classifies a Considerable-damage TOR as TOR_PDS", () => {
      const f = feature("Tornado Warning", { tornadoDamageThreat: ["Considerable"] });
      expect(classifyWarning(f)).toBe("TOR_PDS");
    });

    it("upgrades to TOR_PDS via PDS event prefix", () => {
      const f = feature("Particularly Dangerous Situation Tornado Warning");
      expect(classifyWarning(f)).toBe("TOR_PDS");
    });

    it("classifies a Catastrophic-damage TOR as TOR_EMG", () => {
      const f = feature("Tornado Warning", { tornadoDamageThreat: ["Catastrophic"] });
      expect(classifyWarning(f)).toBe("TOR_EMG");
    });

    it("upgrades to TOR_EMG via 'TORNADO EMERGENCY' in headline", () => {
      const f = feature("Tornado Warning", {}, { headline: "TORNADO EMERGENCY for Moore, OK" });
      expect(classifyWarning(f)).toBe("TOR_EMG");
    });
  });

  describe("Severe Thunderstorm Warning sub-categories", () => {
    it("classifies a base SVR as SVR", () => {
      expect(classifyWarning(feature("Severe Thunderstorm Warning"))).toBe("SVR");
    });
    it("classifies a Considerable SVR as SVR_CON", () => {
      const f = feature("Severe Thunderstorm Warning", { thunderstormDamageThreat: ["Considerable"] });
      expect(classifyWarning(f)).toBe("SVR_CON");
    });
    it("classifies a Destructive SVR as SVR_DST", () => {
      const f = feature("Severe Thunderstorm Warning", { thunderstormDamageThreat: ["Destructive"] });
      expect(classifyWarning(f)).toBe("SVR_DST");
    });
  });

  describe("other event types", () => {
    it("classifies Flash Flood Warning as FFW", () => {
      expect(classifyWarning(feature("Flash Flood Warning"))).toBe("FFW");
    });
    it("classifies Tornado Watch as TOA", () => {
      expect(classifyWarning(feature("Tornado Watch"))).toBe("TOA");
    });
    it("classifies Severe Thunderstorm Watch as SVA", () => {
      expect(classifyWarning(feature("Severe Thunderstorm Watch"))).toBe("SVA");
    });
    it("classifies Special Weather Statement as SPS", () => {
      expect(classifyWarning(feature("Special Weather Statement"))).toBe("SPS");
    });
    it("classifies Flood Warning as FLW", () => {
      expect(classifyWarning(feature("Flood Warning"))).toBe("FLW");
    });
  });
});

/* ─────────────────── firstParam ─────────────────── */

describe("firstParam", () => {
  it("returns the first array element", () => {
    expect(firstParam({ x: ["a", "b"] }, "x")).toBe("a");
  });
  it("returns an empty string for missing keys", () => {
    expect(firstParam({}, "missing")).toBe("");
  });
  it("returns an empty string for non-array values", () => {
    expect(firstParam({ x: "scalar" }, "x")).toBe("");
  });
  it("returns an empty string for empty arrays", () => {
    expect(firstParam({ x: [] }, "x")).toBe("");
  });
  it("trims whitespace from the result", () => {
    expect(firstParam({ x: ["  hello  "] }, "x")).toBe("hello");
  });
});

/* ─────────────────── parseEventMotion ─────────────────── */

describe("parseEventMotion", () => {
  it("parses a fully-populated motion description", () => {
    const m = parseEventMotion("2025-04-30T01:50:00-05:00...260DEG...45KT...3848 9421");
    expect(m).toEqual({ dirDeg: 260, speedKt: 45, mph: Math.round(45 * 1.15078) });
  });
  it("parses just direction when speed is missing", () => {
    const m = parseEventMotion("...260DEG...");
    expect(m.dirDeg).toBe(260);
    expect(Number.isNaN(m.speedKt)).toBe(true);
    expect(Number.isNaN(m.mph)).toBe(true);
  });
  it("parses just speed when direction is missing", () => {
    const m = parseEventMotion("...45KT...");
    expect(Number.isNaN(m.dirDeg)).toBe(true);
    expect(m.speedKt).toBe(45);
    expect(m.mph).toBe(52);
  });
  it("returns null for empty or non-string input", () => {
    expect(parseEventMotion("")).toBeNull();
    expect(parseEventMotion(null)).toBeNull();
    expect(parseEventMotion(123)).toBeNull();
  });
  it("returns null when the string contains no recognizable token", () => {
    expect(parseEventMotion("2025-04-30T01:50:00-05:00")).toBeNull();
  });
});

/* ─────────────────── escapeHtml ─────────────────── */

describe("escapeHtml", () => {
  it("escapes the five canonical HTML characters", () => {
    expect(escapeHtml(`<a href="x" class='y'>&</a>`))
      .toBe("&lt;a href=&quot;x&quot; class=&#39;y&#39;&gt;&amp;&lt;/a&gt;");
  });
  it("renders null and undefined as empty strings", () => {
    expect(escapeHtml(null)).toBe("");
    expect(escapeHtml(undefined)).toBe("");
  });
  it("coerces non-strings via String()", () => {
    expect(escapeHtml(42)).toBe("42");
  });
});

/* ─────────────────── warningStyle ─────────────────── */

describe("warningStyle", () => {
  it("returns Leaflet style values for a Tornado Warning", () => {
    const style = warningStyle(feature("Tornado Warning"));
    const subStyle = WARNING_SUBCAT_STYLES.TOR_RI;
    expect(style.color).toBe(subStyle.color);
    expect(style.weight).toBe(subStyle.weight);
    expect(style.fillColor).toBe(subStyle.color);
    expect(style.fillOpacity).toBe(subStyle.fill);
    expect(style.className).toContain("smap-warn-poly");
    expect(style.className).toContain(subStyle.cls);
  });

  it("falls back to SPS styling for unclassifiable events", () => {
    // Frost Advisory isn't in VALID_WARNING_EVENTS — classify returns null —
    // so warningStyle should still render it (faintly) using SPS.
    const style = warningStyle(feature("Frost Advisory"));
    const subStyle = WARNING_SUBCAT_STYLES.SPS;
    expect(style.color).toBe(subStyle.color);
    expect(style.fillOpacity).toBe(subStyle.fill);
  });
});

/* ─────────────────── exported tables ─────────────────── */

describe("exported metadata", () => {
  it("WARNING_PRIORITY enumerates every key in WARNING_SUBCAT_STYLES exactly once", () => {
    const styleKeys = Object.keys(WARNING_SUBCAT_STYLES).sort();
    const priorityKeys = [...WARNING_PRIORITY].sort();
    expect(priorityKeys).toEqual(styleKeys);
  });

  it("VALID_WARNING_EVENTS is a Set of event-string names", () => {
    expect(VALID_WARNING_EVENTS.has("Tornado Warning")).toBe(true);
    expect(VALID_WARNING_EVENTS.has("Severe Thunderstorm Warning")).toBe(true);
    expect(VALID_WARNING_EVENTS.has("Frost Advisory")).toBe(false);
  });
});
