import { create } from 'zustand';
import { createDataSlice }      from './slices/dataSlice';
import { createResultSlice }    from './slices/resultSlice';
import { createUiSlice }        from './slices/uiSlice';
import { createSelectionSlice } from './slices/selectionSlice';
import { createPrefsSlice }     from './slices/prefsSlice';

/**
 * Single application store, composed from feature slices.
 *
 * The flat shape and key names are preserved verbatim from the pre-slice
 * version, so every existing `useAppStore((s) => s.foo)` selector keeps
 * working. New code should prefer the per-slice creators if it only
 * needs a focused subset of the store.
 */
export const useAppStore = create((set, get, api) => ({
  ...createDataSlice(set, get, api),
  ...createResultSlice(set, get, api),
  ...createUiSlice(set, get, api),
  ...createSelectionSlice(set, get, api),
  ...createPrefsSlice(set, get, api),
}));
