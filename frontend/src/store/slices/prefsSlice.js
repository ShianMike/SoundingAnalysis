function loadPref(key, fallback) {
  try { return localStorage.getItem(key) || fallback; } catch { return fallback; }
}
function savePref(key, val) {
  try { localStorage.setItem(key, val); } catch { /* noop */ }
}

/**
 * Slice creator: persisted user preferences (theme, color-blind mode) plus
 * the active app page.
 *
 * Theme + colorblind are mirrored to `localStorage` and to a `data-*`
 * attribute on `<html>` so plain CSS rules can react to the change without
 * needing React to re-render the entire tree.
 */
export const createPrefsSlice = (set) => ({
  page: "main",
  theme: loadPref("sa_theme", "dark"),
  colorblind: loadPref("sa_cb", "false") === "true",

  setPage: (page) => set({ page }),

  toggleTheme: () => set((state) => {
    const newTheme = state.theme === "dark" ? "light" : "dark";
    savePref("sa_theme", newTheme);
    document.documentElement.setAttribute("data-theme", newTheme);
    return { theme: newTheme };
  }),

  toggleColorblind: () => set((state) => {
    const newCb = !state.colorblind;
    savePref("sa_cb", String(newCb));
    document.documentElement.setAttribute("data-cb", String(newCb));
    return { colorblind: newCb };
  }),
});
