import { useState, useRef, useCallback, useEffect } from "react";

/**
 * `useDraggable` — encapsulates the floating-panel drag/collapse behavior
 * that ControlPanel uses while StationMap is in fullscreen mode.
 *
 * Behavior preserved verbatim from the original ControlPanel implementation:
 *   - Position is x/y in viewport pixels.
 *   - Panel snaps to (0, 0) and auto-collapses when fullscreen activates,
 *     and resets to (16, 16) and auto-expands when fullscreen deactivates.
 *   - While fullscreen is active, the panel keeps its top edge below the
 *     map toolbar (which is dynamically resized via a ResizeObserver).
 *   - Drag is initiated only when the user grabs the drag handle (not a
 *     button inside it) and only while fullscreen is active.
 *
 * Returns:
 *   {
 *     dragPos,                   // { x, y } current position
 *     floatCollapsed,            // bool — collapsed (header-only) state
 *     toggleFloatCollapsed,      // () => void
 *     dragRef,                   // ref to attach to the panel element
 *     dragHandlers: {            // spread onto the drag handle element
 *       onPointerDown,
 *       onPointerMove,
 *       onPointerUp,
 *     }
 *   }
 *
 * @param {Object} [opts]
 * @param {string} [opts.fullscreenClass="smap-fullscreen-active"] - body class
 *   that signals fullscreen mode is active.
 * @param {string} [opts.toolbarSelector=".smap-fullscreen .smap-toolbar"] -
 *   CSS selector for the toolbar whose bottom edge defines the panel's
 *   minimum Y position.
 * @param {{x:number,y:number}} [opts.initialPos={x:16,y:16}] - starting position.
 * @param {number} [opts.fallbackPanelWidth=270] - width to use when measuring
 *   the panel for clamping before its real width is known.
 */
export function useDraggable({
  fullscreenClass = "smap-fullscreen-active",
  toolbarSelector = ".smap-fullscreen .smap-toolbar",
  initialPos = { x: 16, y: 16 },
  fallbackPanelWidth = 270,
} = {}) {
  const [floatCollapsed, setFloatCollapsed] = useState(false);
  const [dragPos, setDragPos] = useState(initialPos);

  const dragRef = useRef(null);
  const isDragging = useRef(false);
  const dragStart = useRef({ x: 0, y: 0, posX: 0, posY: 0 });

  /** Minimum Y = bottom of the map toolbar + gap so the panel never overlaps it */
  const getMinY = useCallback(() => {
    const tb = document.querySelector(toolbarSelector);
    if (tb) {
      const r = tb.getBoundingClientRect();
      return Math.max(r.bottom + 8, 0);
    }
    return 90; // safe fallback
  }, [toolbarSelector]);

  const onPointerDown = useCallback((e) => {
    if (!document.body.classList.contains(fullscreenClass)) return;
    // Don't start drag if clicking a button (let the click event fire)
    if (e.target.closest("button")) return;
    isDragging.current = true;
    dragStart.current = {
      x: e.clientX,
      y: e.clientY,
      posX: dragPos.x,
      posY: dragPos.y,
    };
    e.currentTarget.setPointerCapture(e.pointerId);
  }, [dragPos, fullscreenClass]);

  const onPointerMove = useCallback((e) => {
    if (!isDragging.current) return;
    const dx = e.clientX - dragStart.current.x;
    const dy = e.clientY - dragStart.current.y;
    const minY = getMinY();
    const pad = 12;
    const panelW = dragRef.current ? dragRef.current.offsetWidth : fallbackPanelWidth;
    setDragPos({
      x: Math.max(pad, Math.min(window.innerWidth - panelW - pad, dragStart.current.posX + dx)),
      y: Math.max(minY, Math.min(window.innerHeight - 40, dragStart.current.posY + dy)),
    });
  }, [getMinY, fallbackPanelWidth]);

  const onPointerUp = useCallback(() => {
    isDragging.current = false;
  }, []);

  // Auto-position when fullscreen activates / deactivates. Only reacts to
  // *changes* in the fullscreen class (not every class mutation).
  const wasFullscreen = useRef(false);
  useEffect(() => {
    const obs = new MutationObserver(() => {
      const active = document.body.classList.contains(fullscreenClass);
      if (active === wasFullscreen.current) return;
      wasFullscreen.current = active;
      if (active) {
        setDragPos({ x: 0, y: 0 });
        setFloatCollapsed(true);
      } else {
        setDragPos(initialPos);
        setFloatCollapsed(false);
      }
    });
    obs.observe(document.body, { attributes: true, attributeFilter: ["class"] });
    return () => obs.disconnect();
    // initialPos & fullscreenClass intentionally captured once at mount.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // Keep the panel below the toolbar when the toolbar height changes.
  useEffect(() => {
    const tb = document.querySelector(toolbarSelector);
    if (!tb) return;
    const ro = new ResizeObserver(() => {
      if (!document.body.classList.contains(fullscreenClass)) return;
      const bottom = tb.getBoundingClientRect().bottom;
      setDragPos((prev) => ({ x: prev.x, y: (bottom > 0 ? bottom : 0) + 8 }));
    });
    ro.observe(tb);
    return () => ro.disconnect();
  });

  const toggleFloatCollapsed = useCallback(
    () => setFloatCollapsed((v) => !v),
    []
  );

  return {
    dragPos,
    floatCollapsed,
    setFloatCollapsed,
    toggleFloatCollapsed,
    dragRef,
    dragHandlers: { onPointerDown, onPointerMove, onPointerUp },
  };
}
