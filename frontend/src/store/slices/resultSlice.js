/**
 * Slice creator: in-flight network / result state for sounding fetches.
 *
 *   - `loading`         active sounding fetch in progress
 *   - `initialLoading`  first-ever boot fetch (stations + sources + models)
 *   - `error`           last network error (string) — null when idle
 *   - `result`          decoded sounding payload from `/api/sounding`
 *   - `riskData`        last severe-risk scan payload
 *   - `lastParams`      the params used for the *successful* sounding fetch,
 *                       used by auto-refresh + URL-share
 *   - `compareHistoryData` payload pushed from HistoryPanel into ComparisonView
 */
export const createResultSlice = (set) => ({
  loading: false,
  initialLoading: true,
  error: null,
  result: null,
  riskData: null,
  lastParams: null,
  compareHistoryData: null,

  setLoading: (loading) => set({ loading }),
  setInitialLoading: (initialLoading) => set({ initialLoading }),
  setError: (error) => set({ error }),
  setResult: (result) => set({ result }),
  setRiskData: (riskData) => set({ riskData }),
  setLastParams: (lastParams) => set({ lastParams }),
  setCompareHistoryData: (compareHistoryData) => set({ compareHistoryData }),
});
