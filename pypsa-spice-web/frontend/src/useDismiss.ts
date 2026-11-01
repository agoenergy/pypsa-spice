import { useEffect } from "react";
import type { RefObject } from "react";

/** Calls onDismiss when Escape is pressed, while active. */
export function useDismissOnEscape(active: boolean, onDismiss: () => void) {
  useEffect(() => {
    if (!active) return;
    const dismiss = (event: KeyboardEvent) => {
      if (event.key === "Escape") onDismiss();
    };
    document.addEventListener("keydown", dismiss);
    return () => document.removeEventListener("keydown", dismiss);
  }, [active, onDismiss]);
}

/** Calls onDismiss on Escape, or on a pointer press outside container, while active. */
export function useDismissOnEscapeOrOutside(
  active: boolean,
  container: RefObject<HTMLElement | null>,
  onDismiss: () => void,
) {
  useDismissOnEscape(active, onDismiss);
  useEffect(() => {
    if (!active) return;
    const dismiss = (event: PointerEvent) => {
      if (!container.current?.contains(event.target as Node)) onDismiss();
    };
    document.addEventListener("pointerdown", dismiss);
    return () => document.removeEventListener("pointerdown", dismiss);
  }, [active, container, onDismiss]);
}
