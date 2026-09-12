import { useEffect, useState } from "react";

import chartSurfaceStyles from "./ChartSurface.module.scss";
import workspaceFeedbackStyles from "./WorkspaceFeedback.module.scss";

const LOADING_INDICATOR_DELAY_MS = 500;

/**
 * True once `active` has held for `delayMs`, false as soon as it clears. Keeps a
 * spinner from flashing on requests that return quickly.
 */
export function useDelayedFlag(active: boolean, delayMs = LOADING_INDICATOR_DELAY_MS): boolean {
  const [raised, setRaised] = useState(false);
  useEffect(() => {
    if (!active) {
      setRaised(false);
      return;
    }
    const timer = window.setTimeout(() => setRaised(true), delayMs);
    return () => window.clearTimeout(timer);
  }, [active, delayMs]);
  return raised;
}

/** Placeholder shown while a chart's result tables are still being read. */
export function ChartLoadingState({ message }: { message: string }) {
  return (
    <div className={workspaceFeedbackStyles["state"]}>
      <span className={workspaceFeedbackStyles["spinner"]} />
      {message}
    </div>
  );
}

/** Placeholder shown when a chart request failed. */
export function ChartErrorState({ message }: { message: string }) {
  return (
    <div className={[workspaceFeedbackStyles["state"], workspaceFeedbackStyles["empty"]].join(" ")}>
      <b>No chart data</b>
      <span>{message}</span>
    </div>
  );
}

/** Overlay shown while an already-drawn chart refetches, for example on a new time range. */
export function ChartRefreshOverlay() {
  return (
    <div className={chartSurfaceStyles["hourly-loading-overlay"]} role="status" aria-live="polite">
      <span className={workspaceFeedbackStyles["spinner"]} aria-hidden="true" />
      <span>Updating chart…</span>
    </div>
  );
}

/** Eyebrow and title naming which scenario a plot belongs to. */
export function ChartContextHeading({ label, title }: { label: string; title: string }) {
  return (
    <div className={chartSurfaceStyles["scenario-label"]}>
      <small>{label}</small>
      <h4>{title}</h4>
    </div>
  );
}
