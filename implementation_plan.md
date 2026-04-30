# Refactor Frontend Architecture

This implementation plan addresses critical performance bottlenecks and technical debt in the Sounding Analysis frontend. The goal is to optimize state management, break down monolithic UI components, and separate business logic from presentation. 

**Core Refactoring Rule:** If any code file has too many lines, it must be refactored and split into smaller, modular components.

**Note: Once this structural refactor is complete, we will migrate the codebase to React + TypeScript.** This will enforce type safety for the complex meteorological data structures and significantly improve maintainability going forward.

## User Review Required

> [!WARNING]
> This plan outlines a significant structural reorganization of the frontend codebase. While the visual UI and functionality will remain identical, many files will be split, moved, or updated. Please review the proposed separations to ensure they align with how you envision the project's structure.

## Proposed Changes

### State Management Optimization

#### [MODIFY] `App.jsx`, `ControlPanel.jsx`, `ResultsView.jsx`
- Replace wholesale destructuring of the Zustand `useAppStore` with individual, granular selector functions or shallow state selection. 
- This will prevent these massive components from re-rendering every time any unrelated value in the global store changes, resolving a major performance bottleneck.

---

### Component Decomposition: ControlPanel

#### [MODIFY] `ControlPanel.jsx`
- Strip out the massive inline sub-components and logic blocks. `ControlPanel` will become a simple container layout that imports its internal tabs and logic.

#### [NEW] `src/hooks/useDraggable.js`
- Extract the drag-and-drop float logic and mouse event handlers from `ControlPanel.jsx` into a reusable custom hook.

#### [NEW] `src/components/ControlPanel/DataTab.jsx`
#### [NEW] `src/components/ControlPanel/ModifyTab.jsx`
#### [NEW] `src/components/ControlPanel/ToolsTab.jsx`
#### [NEW] `src/components/ControlPanel/StationPicker.jsx`
- Isolate the specific UI and form logic for each of the Control Panel's tabs into dedicated, focused components.

---

### Component Decomposition: ResultsView

#### [MODIFY] `ResultsView.jsx`
- Remove the inline definition of the `RiskTable` component and the `generateSoundingSummary` function. Update the file to import these from their new isolated locations.

#### [NEW] `src/components/RiskTable.jsx`
- Extract the `RiskTable` UI component into its own file for better maintainability.

#### [NEW] `src/utils/soundingSummary.js`
- Extract the massive 300+ line `generateSoundingSummary` pure function from the UI layer into a dedicated utility file.

#### [NEW] `src/utils/exportImage.js`
- Extract the raw HTML5 Canvas drawing and image export logic from the `RiskTable` component into a reusable utility function.

---

### Configuration Extraction

#### [NEW] `src/config/constants.js`
- Extract all hardcoded configuration objects and arrays (e.g., `NEXRAD_SITES`, `SOURCE_META`, and `MODEL_META`) currently sitting at the top of UI files into a central configuration file.

## Verification Plan

### Manual Verification
- Launch the development server and verify that the application loads without errors.
- Test the Control Panel: Ensure tabs switch correctly, station search works, and the floating drag-and-drop functionality behaves normally.
- Test Results View: Generate a sounding, verify the text summary loads correctly, and ensure the risk table rendering and image export functions still work.
- Use React Developer Tools (Profiler) to verify that interacting with inputs (like typing in the search box) no longer triggers re-renders of the entire `ResultsView` or `App`.
