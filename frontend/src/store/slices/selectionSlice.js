/**
 * Slice creator: which station + data source the user is currently asking
 * about. Lives in its own slice so a station-list re-sort doesn't cause an
 * unrelated panel to re-render.
 */
export const createSelectionSlice = (set) => ({
  selectedStation: "OUN",
  source: "obs",
  setSelectedStation: (selectedStation) => set({ selectedStation }),
  setSource: (source) => set({ source }),
});
