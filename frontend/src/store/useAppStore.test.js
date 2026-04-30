import { describe, it, expect, beforeEach } from "vitest";
import { useAppStore } from "./useAppStore";

/**
 * Smoke tests for the sliced store. Verify that every public key the
 * pre-Phase-2 flat store exposed is still present, and that toggles
 * work end-to-end.
 */

const SNAPSHOT_KEYS = [
  // Data slice
  "stations", "sources", "models", "psuModels", "setGlobalData",
  // Result slice
  "loading", "initialLoading", "error", "result", "riskData",
  "lastParams", "compareHistoryData",
  "setLoading", "setInitialLoading", "setError", "setResult",
  "setRiskData", "setLastParams", "setCompareHistoryData",
  // UI slice
  "showHistory", "showMap", "showRisk", "showTimeSeries", "showCompare",
  "showVwp", "showFeedback",
  "toggleHistory",    "setShowHistory",
  "toggleMap",        "setShowMap",
  "toggleRisk",       "setShowRisk",
  "toggleTimeSeries", "setShowTimeSeries",
  "toggleCompare",    "setShowCompare",
  "toggleVwp",        "setShowVwp",
  "toggleFeedback",   "setShowFeedback",
  // Selection slice
  "selectedStation", "source",
  "setSelectedStation", "setSource",
  // Prefs slice
  "page", "theme", "colorblind",
  "setPage", "toggleTheme", "toggleColorblind",
];

/** Reset store between tests via `setState` since there is no provider. */
function resetStore() {
  useAppStore.setState({
    stations: [], sources: [], models: [], psuModels: [],
    loading: false, initialLoading: true, error: null,
    result: null, riskData: null,
    lastParams: null, compareHistoryData: null,
    showHistory: false, showMap: true, showRisk: false,
    showTimeSeries: false, showCompare: false, showVwp: false,
    showFeedback: false,
    selectedStation: "OUN", source: "obs",
    page: "main",
  });
}

beforeEach(() => resetStore());

describe("useAppStore — public API", () => {
  it("exposes every key the pre-slice flat store exported", () => {
    const s = useAppStore.getState();
    for (const key of SNAPSHOT_KEYS) {
      expect(s, `missing key: ${key}`).toHaveProperty(key);
    }
  });

  it("exposes a function for every setter / toggle key", () => {
    const s = useAppStore.getState();
    for (const key of SNAPSHOT_KEYS) {
      if (key.startsWith("set") || key.startsWith("toggle")) {
        expect(typeof s[key], `${key} should be a function`).toBe("function");
      }
    }
  });
});

describe("UI slice — toggles flip the matching boolean", () => {
  const pairs = [
    ["showHistory",    "toggleHistory"],
    ["showMap",        "toggleMap"],
    ["showRisk",       "toggleRisk"],
    ["showTimeSeries", "toggleTimeSeries"],
    ["showCompare",    "toggleCompare"],
    ["showVwp",        "toggleVwp"],
    ["showFeedback",   "toggleFeedback"],
  ];
  for (const [flag, toggle] of pairs) {
    it(`${toggle}() inverts ${flag}`, () => {
      const before = useAppStore.getState()[flag];
      useAppStore.getState()[toggle]();
      const after = useAppStore.getState()[flag];
      expect(after).toBe(!before);
    });
  }
});

describe("UI slice — setShow* writes the boolean", () => {
  it("setShowMap(false) sets showMap to false", () => {
    useAppStore.getState().setShowMap(false);
    expect(useAppStore.getState().showMap).toBe(false);
  });
  it("setShowHistory(true) sets showHistory to true", () => {
    useAppStore.getState().setShowHistory(true);
    expect(useAppStore.getState().showHistory).toBe(true);
  });
});

describe("Data slice — setGlobalData merges into the store", () => {
  it("merges stations / sources / models from a single payload", () => {
    useAppStore.getState().setGlobalData({
      stations: [{ id: "OUN", lat: 35.2, lon: -97.4, name: "Norman" }],
      sources: [{ id: "obs" }],
      models: [{ id: "hrrr", name: "HRRR" }],
      psuModels: [],
    });
    const s = useAppStore.getState();
    expect(s.stations).toHaveLength(1);
    expect(s.sources).toHaveLength(1);
    expect(s.models).toHaveLength(1);
  });
});

describe("Result slice — setters write the matching field", () => {
  it("setError stores the error string", () => {
    useAppStore.getState().setError("network down");
    expect(useAppStore.getState().error).toBe("network down");
  });
  it("setResult stores the payload reference", () => {
    const payload = { params: {}, image: "abc" };
    useAppStore.getState().setResult(payload);
    expect(useAppStore.getState().result).toBe(payload);
  });
  it("setLastParams stores the params object", () => {
    const params = { source: "obs", station: "OUN" };
    useAppStore.getState().setLastParams(params);
    expect(useAppStore.getState().lastParams).toBe(params);
  });
});

describe("Selection slice", () => {
  it("setSelectedStation updates the station", () => {
    useAppStore.getState().setSelectedStation("OAX");
    expect(useAppStore.getState().selectedStation).toBe("OAX");
  });
  it("setSource updates the source", () => {
    useAppStore.getState().setSource("bufkit");
    expect(useAppStore.getState().source).toBe("bufkit");
  });
});

describe("Prefs slice", () => {
  it("toggleTheme flips between dark and light", () => {
    const before = useAppStore.getState().theme;
    useAppStore.getState().toggleTheme();
    const after = useAppStore.getState().theme;
    expect(after).not.toBe(before);
    expect(["dark", "light"]).toContain(after);
  });

  it("toggleColorblind flips the boolean", () => {
    const before = useAppStore.getState().colorblind;
    useAppStore.getState().toggleColorblind();
    expect(useAppStore.getState().colorblind).toBe(!before);
  });

  it("setPage updates the page", () => {
    useAppStore.getState().setPage("ensemble");
    expect(useAppStore.getState().page).toBe("ensemble");
  });
});
