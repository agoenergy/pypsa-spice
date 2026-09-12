import { useEffect, useState } from "react";
import type { RefObject } from "react";

/**
 * True once the element has scrolled to within `rootMargin` of the viewport, and
 * stays true afterwards. Results pages mount up to 16 chart cards at once, so
 * gating each card's request on this keeps a page load from firing all of them.
 */
export default function useNearViewport(ref: RefObject<Element | null>, rootMargin = "600px"): boolean {
  const [near, setNear] = useState(false);
  useEffect(() => {
    if (near) return;
    const element = ref.current;
    if (!element) return;
    if (typeof IntersectionObserver === "undefined") {
      setNear(true);
      return;
    }
    const observer = new IntersectionObserver(
      (entries) => {
        if (entries.some((entry) => entry.isIntersecting)) setNear(true);
      },
      { rootMargin },
    );
    observer.observe(element);
    return () => observer.disconnect();
  }, [ref, rootMargin, near]);
  return near;
}
