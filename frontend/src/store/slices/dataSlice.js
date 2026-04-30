/**
 * Slice creator: app data fetched once at boot (stations / sources / models).
 * Pure data, no UI state.
 */
export const createDataSlice = (set) => ({
  stations: [],
  sources: [],
  models: [],
  psuModels: [],
  /** Replace any subset of the boot data in one shot. */
  setGlobalData: (data) => set(data),
});
