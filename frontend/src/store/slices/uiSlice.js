/**
 * Slice creator: panel-visibility flags and their toggles.
 *
 * Each toggle pairs a `show*` boolean with `set` and `toggle` helpers so
 * consumers can pick the narrowest API they need.
 */
export const createUiSlice = (set) => ({
  showHistory: false,
  showMap: true,
  showRisk: false,
  showTimeSeries: false,
  showCompare: false,
  showVwp: false,
  showFeedback: false,

  toggleHistory:    () => set((s) => ({ showHistory:    !s.showHistory })),
  setShowHistory:   (val) => set({ showHistory: val }),

  toggleMap:        () => set((s) => ({ showMap:        !s.showMap })),
  setShowMap:       (val) => set({ showMap: val }),

  toggleRisk:       () => set((s) => ({ showRisk:       !s.showRisk })),
  setShowRisk:      (val) => set({ showRisk: val }),

  toggleTimeSeries: () => set((s) => ({ showTimeSeries: !s.showTimeSeries })),
  setShowTimeSeries:(val) => set({ showTimeSeries: val }),

  toggleCompare:    () => set((s) => ({ showCompare:    !s.showCompare })),
  setShowCompare:   (val) => set({ showCompare: val }),

  toggleVwp:        () => set((s) => ({ showVwp:        !s.showVwp })),
  setShowVwp:       (val) => set({ showVwp: val }),

  toggleFeedback:   () => set((s) => ({ showFeedback:   !s.showFeedback })),
  setShowFeedback:  (val) => set({ showFeedback: val }),
});
