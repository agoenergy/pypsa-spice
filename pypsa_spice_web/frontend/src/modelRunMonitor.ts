import { getLatestModelRun, getModelRun } from "./api";
import type { ModelRun, ModelRunMonitorState } from "./types";

const POLL_DELAY_MS = 1000;
const MAX_RETRY_DELAY_MS = 30000;
const REQUEST_TIMEOUT_MS = 15000;

export function isActiveModelRun(run: ModelRun | null): boolean {
  return Boolean(run && ["queued", "running", "canceling"].includes(run.status));
}

// One monitor belongs to App and survives navigation between workspace pages.
export function createModelRunMonitor(
  onChange: (state: ModelRunMonitorState) => void,
  readLatest = getLatestModelRun,
  readRun = getModelRun,
) {
  let state: ModelRunMonitorState = { run: null, error: "", loading: true };
  let stopped = true;
  let generation = 0;
  let failures = 0;
  let timer: ReturnType<typeof setTimeout> | undefined;
  let deadline: ReturnType<typeof setTimeout> | undefined;
  let controller: AbortController | undefined;

  const publish = (next: ModelRunMonitorState) => {
    state = next;
    onChange(state);
  };

  const invalidate = () => {
    generation += 1;
    clearTimeout(timer);
    clearTimeout(deadline);
    controller?.abort();
    controller = undefined;
  };

  const schedule = (delay: number, latest = false) => {
    if (!stopped) timer = setTimeout(() => void poll(latest), delay);
  };

  const poll = async (latest = false) => {
    const requestGeneration = generation;
    const requestController = new AbortController();
    controller = requestController;
    const isCurrent = () => !stopped && generation === requestGeneration;
    let timedOut = false;
    deadline = setTimeout(() => {
      timedOut = true;
      requestController.abort();
    }, REQUEST_TIMEOUT_MS);

    try {
      const run =
        latest || !state.run
          ? await readLatest(requestController.signal)
          : await readRun(state.run.id, requestController.signal);
      if (!isCurrent()) return;
      if (timedOut) throw new Error("The status request timed out.");
      failures = 0;
      publish({ run, error: "", loading: false });
      if (isActiveModelRun(run)) schedule(POLL_DELAY_MS);
    } catch {
      if (!isCurrent()) return;
      failures += 1;
      publish({ ...state, error: "Connection lost; retrying. Showing the last known run status.", loading: false });
      schedule(Math.min(POLL_DELAY_MS * 2 ** Math.min(failures, 5), MAX_RETRY_DELAY_MS), latest);
    } finally {
      if (isCurrent()) {
        clearTimeout(deadline);
        controller = undefined;
      }
    }
  };

  return {
    start() {
      invalidate();
      stopped = false;
      failures = 0;
      void poll(true);
    },
    stop() {
      stopped = true;
      invalidate();
    },
    updateRun(run: ModelRun | null) {
      if (stopped) return;
      // A start/cancel response supersedes any pending status request.
      invalidate();
      failures = 0;
      publish({ run, error: "", loading: false });
      if (isActiveModelRun(run)) schedule(POLL_DELAY_MS);
    },
  };
}
