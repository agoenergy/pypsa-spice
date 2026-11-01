import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { createModelRunMonitor } from "./modelRunMonitor";
import type { ModelRun, ModelRunMonitorState } from "./types";

function run(status: ModelRun["status"] = "running", id = "run-1"): ModelRun {
  return {
    id,
    status,
    progress: 20,
    message: "Solving",
    current_rule: null,
    created_at: "2026-09-11T10:00:00Z",
    started_at: null,
    ended_at: null,
    exit_code: null,
    pid: null,
    dataset: "example",
    project: "project_01",
    input_scenario: "scenario_01",
    output_scenario: "scenario_01_run",
    cores: 1,
    target: "solve_all_networks",
    config_file: "config.yaml",
    log_file: "run.log",
    manifest_file: "manifest.json",
    log: "Solving",
  };
}

function deferred<T>() {
  let resolve!: (value: T) => void;
  let reject!: (reason: Error) => void;
  const promise = new Promise<T>((onResolve, onReject) => {
    resolve = onResolve;
    reject = onReject;
  });
  return { promise, resolve, reject };
}

function setup(initial: ModelRun | null = run()) {
  const onChange = vi.fn<(state: ModelRunMonitorState) => void>();
  const readLatest = vi.fn<(signal?: AbortSignal) => Promise<ModelRun | null>>().mockResolvedValue(initial);
  const readRun = vi.fn<(id: string, signal?: AbortSignal) => Promise<ModelRun>>().mockResolvedValue(run());
  const monitor = createModelRunMonitor(onChange, readLatest, readRun);
  const state = () => onChange.mock.lastCall![0];
  return { monitor, onChange, readLatest, readRun, state };
}

beforeEach(() => vi.useFakeTimers());
afterEach(() => {
  vi.clearAllTimers();
  vi.useRealTimers();
});

describe("model run monitoring", () => {
  it("waits for the response and then one second before requesting again", async () => {
    const { monitor, readRun } = setup();
    const pending = deferred<ModelRun>();
    readRun.mockReturnValueOnce(pending.promise);
    monitor.start();
    await vi.advanceTimersByTimeAsync(1000);
    expect(readRun).toHaveBeenCalledTimes(1);
    await vi.advanceTimersByTimeAsync(5000);
    expect(readRun).toHaveBeenCalledTimes(1);
    pending.resolve(run());
    await vi.advanceTimersByTimeAsync(999);
    expect(readRun).toHaveBeenCalledTimes(1);
    await vi.advanceTimersByTimeAsync(1);
    expect(readRun).toHaveBeenCalledTimes(2);
  });

  it.each(["succeeded", "failed", "canceled"] as const)("stops after the run is %s", async (status) => {
    const { monitor, readRun, state } = setup();
    readRun.mockResolvedValue(run(status));
    monitor.start();
    await vi.advanceTimersByTimeAsync(60000);
    expect(readRun).toHaveBeenCalledTimes(1);
    expect(state().run?.status).toBe(status);
    expect(vi.getTimerCount()).toBe(0);
  });

  it.each([null, run("succeeded")])("does not poll an idle workspace", async (initial) => {
    const { monitor, readLatest, readRun } = setup(initial);
    monitor.start();
    await vi.advanceTimersByTimeAsync(60000);
    expect(readLatest).toHaveBeenCalledTimes(1);
    expect(readRun).not.toHaveBeenCalled();
    expect(vi.getTimerCount()).toBe(0);
  });

  it("continues across queued, running and canceling statuses", async () => {
    const { monitor, readRun, state } = setup(run("queued"));
    readRun.mockResolvedValueOnce(run()).mockResolvedValueOnce(run("canceling")).mockResolvedValueOnce(run("canceled"));
    monitor.start();
    await vi.advanceTimersByTimeAsync(3000);
    expect(readRun).toHaveBeenCalledTimes(3);
    expect(state().run?.status).toBe("canceled");
  });

  it("aborts on cleanup and ignores late responses and mutations", async () => {
    const { monitor, readRun, onChange } = setup();
    const pending = deferred<ModelRun>();
    readRun.mockReturnValue(pending.promise);
    monitor.start();
    await vi.advanceTimersByTimeAsync(1000);
    const signal = readRun.mock.lastCall![1]!;
    monitor.stop();
    expect(signal.aborted).toBe(true);
    onChange.mockClear();
    pending.resolve(run("failed"));
    monitor.updateRun(run("queued", "run-2"));
    await vi.advanceTimersByTimeAsync(60000);
    expect(onChange).not.toHaveBeenCalled();
    expect(vi.getTimerCount()).toBe(0);
  });

  it("does not let a pending poll overwrite a cancellation response", async () => {
    const { monitor, readRun, state } = setup();
    const pending = deferred<ModelRun>();
    readRun.mockReturnValueOnce(pending.promise).mockResolvedValue(run("canceled"));
    monitor.start();
    await vi.advanceTimersByTimeAsync(1000);
    monitor.updateRun(run("canceling"));
    expect(readRun.mock.lastCall![1]!.aborted).toBe(true);
    pending.resolve(run());
    await vi.advanceTimersByTimeAsync(0);
    expect(state().run?.status).toBe("canceling");
    await vi.advanceTimersByTimeAsync(1000);
    expect(state().run?.status).toBe("canceled");
  });

  it("does not let the initial lookup overwrite a newly started run", async () => {
    const { monitor, readLatest, readRun, state } = setup();
    const pending = deferred<ModelRun | null>();
    readLatest.mockReturnValue(pending.promise);
    monitor.start();
    monitor.updateRun(run("queued", "new-run"));
    pending.resolve(run("succeeded", "old-run"));
    await vi.advanceTimersByTimeAsync(0);
    expect(state().run?.id).toBe("new-run");
    await vi.advanceTimersByTimeAsync(1000);
    expect(readRun.mock.lastCall![0]).toBe("new-run");
  });

  it("retains status on failure, caps retry delays, and resets the delay after recovery", async () => {
    const { monitor, readRun, state } = setup();
    readRun.mockRejectedValue(new Error("offline"));
    monitor.start();
    await vi.advanceTimersByTimeAsync(1000);
    expect(state().run?.status).toBe("running");
    expect(state().error).toContain("Connection lost; retrying");
    let calls = 1;
    for (const delay of [2000, 4000, 8000, 16000, 30000, 30000]) {
      await vi.advanceTimersByTimeAsync(delay - 1);
      expect(readRun).toHaveBeenCalledTimes(calls);
      await vi.advanceTimersByTimeAsync(1);
      expect(readRun).toHaveBeenCalledTimes(++calls);
    }
    readRun.mockResolvedValueOnce(run());
    await vi.advanceTimersByTimeAsync(30000);
    expect(state().error).toBe("");
    await vi.advanceTimersByTimeAsync(1000);
    const callsAfterRecovery = readRun.mock.calls.length;
    await vi.advanceTimersByTimeAsync(1999);
    expect(readRun).toHaveBeenCalledTimes(callsAfterRecovery);
    await vi.advanceTimersByTimeAsync(1);
    expect(readRun).toHaveBeenCalledTimes(callsAfterRecovery + 1);
  });

  it("retries a failed initial lookup", async () => {
    const { monitor, readLatest, state } = setup(null);
    readLatest.mockRejectedValueOnce(new Error("offline"));
    monitor.start();
    await vi.advanceTimersByTimeAsync(0);
    expect(state().loading).toBe(false);
    expect(state().error).not.toBe("");
    await vi.advanceTimersByTimeAsync(2000);
    expect(readLatest).toHaveBeenCalledTimes(2);
    expect(state().error).toBe("");
    expect(vi.getTimerCount()).toBe(0);
  });

  it("aborts a timed-out request and retries", async () => {
    const { monitor, readRun, state } = setup();
    readRun.mockImplementationOnce(
      (_id, signal) =>
        new Promise((_resolve, reject) => {
          signal!.addEventListener("abort", () => reject(new Error("aborted")), { once: true });
        }),
    );
    monitor.start();
    await vi.advanceTimersByTimeAsync(16000);
    expect(readRun.mock.lastCall![1]!.aborted).toBe(true);
    expect(state().error).not.toBe("");
    await vi.advanceTimersByTimeAsync(2000);
    expect(readRun).toHaveBeenCalledTimes(2);
    expect(state().error).toBe("");
  });

  it("ignores an old request after a stop/start cycle", async () => {
    const { monitor, readLatest, state } = setup(run("queued", "new-run"));
    const pending = deferred<ModelRun | null>();
    readLatest.mockReturnValueOnce(pending.promise);
    monitor.start();
    monitor.stop();
    monitor.start();
    await vi.advanceTimersByTimeAsync(0);
    pending.reject(new Error("late failure"));
    await vi.advanceTimersByTimeAsync(0);
    expect(state().run?.id).toBe("new-run");
    expect(state().error).toBe("");
  });
});
