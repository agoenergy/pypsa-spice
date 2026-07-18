import { useEffect, useMemo, useRef, useState } from "react";
import { AlertTriangle, CheckCircle2, CircleStop, Clock3, ExternalLink, FileCode2, GitBranch, Pencil, Play, RotateCcw, SquareTerminal } from "lucide-react";
import { cancelModelRun, getLatestModelRun, getModelRun, getModelRunOptions, getScenarioConfig, getScenarioWorkspaceStatus, startModelRun } from "./api";
import PageHeader from "./PageHeader";
import type { InputSelection, ModelRun, ModelRunOptions, ScenarioConfigResponse, ScenarioWorkspaceStatus } from "./types";

const activeStatuses = new Set(["queued", "running", "canceling"]);
const CHILLI_COUNT = 12;

function readableTime(value: string | null): string {
  if (!value) return "—";
  return new Intl.DateTimeFormat(undefined, {
    dateStyle: "medium",
    timeStyle: "medium",
  }).format(new Date(value));
}

function statusLabel(status: ModelRun["status"]): string {
  return status.charAt(0).toUpperCase() + status.slice(1);
}

export default function RunModel({ selection, onEditConfiguration, onOpenResults }: { selection: InputSelection; onEditConfiguration: () => void; onOpenResults: (runName: string, dataset: string, project: string) => void }) {
  const [options, setOptions] = useState<ModelRunOptions | null>(null);
  const [scenarioConfig, setScenarioConfig] = useState<ScenarioConfigResponse | null>(null);
  const [workspace, setWorkspace] = useState<ScenarioWorkspaceStatus | null>(null);
  const [run, setRun] = useState<ModelRun | null>(null);
  const [outputScenario, setOutputScenario] = useState("");
  const [cores, setCores] = useState(1);
  const [loading, setLoading] = useState(true);
  const [submitting, setSubmitting] = useState(false);
  const [error, setError] = useState("");
  const logRef = useRef<HTMLPreElement>(null);

  useEffect(() => {
    let current = true;
    setLoading(true);
    Promise.all([getModelRunOptions(), getLatestModelRun(), getScenarioConfig(selection), getScenarioWorkspaceStatus(selection.dataset)])
      .then(([nextOptions, latest, nextConfig, nextWorkspace]) => {
        if (!current) return;
        setOptions(nextOptions);
        setScenarioConfig(nextConfig);
        setWorkspace(nextWorkspace);
        setRun(latest);
        const defaultsMatch =
          nextOptions.defaults.dataset === selection.dataset
          && nextOptions.defaults.project === selection.project
          && nextOptions.defaults.input_scenario === selection.scenario;
        setOutputScenario(defaultsMatch ? nextOptions.defaults.output_scenario : `${selection.scenario}_run`);
        setError("");
      })
      .catch((reason) => {
        if (current) setError(reason instanceof Error ? reason.message : "Could not load the run workspace.");
      })
      .finally(() => { if (current) setLoading(false); });
    return () => { current = false; };
  }, [selection.dataset, selection.project, selection.scenario]);

  useEffect(() => {
    if (!run || !activeStatuses.has(run.status)) return;
    const timer = window.setInterval(() => {
      getModelRun(run.id)
        .then((nextRun) => { setRun(nextRun); setError(""); })
        .catch((reason) => setError(reason instanceof Error ? reason.message : "Could not refresh the model run."));
    }, 1000);
    return () => window.clearInterval(timer);
  }, [run?.id, run?.status]);

  useEffect(() => {
    if (logRef.current) logRef.current.scrollTop = logRef.current.scrollHeight;
  }, [run?.log]);

  const active = Boolean(run && activeStatuses.has(run.status));
  const runMatchesSelection = Boolean(
    run
    && run.dataset === selection.dataset
    && run.project === selection.project
    && run.input_scenario === selection.scenario,
  );
  const dimensions = useMemo(() => {
    if (!options) return [];
    return [
      { label: "Years", value: options.years.join(", ") || "Not configured" },
      { label: "Sector", value: options.sectors.join(", ") || "Not configured" },
      { label: "Regions", value: options.regions.join(", ") || "Not configured" },
      { label: "Currency", value: options.currency || "Not configured" },
    ];
  }, [options]);
  const scenarioSummary = useMemo(() => {
    const settings = scenarioConfig?.sections.scenario_configs || {};
    const snapshots = (settings.snapshots || {}) as Record<string, unknown>;
    const resolution = (settings.resolution || {}) as Record<string, unknown>;
    const resolutionValue = resolution.method === "clustered"
      ? `${resolution.number_of_days || "—"} representative days`
      : `Every ${resolution.stepsize || "—"} hour${Number(resolution.stepsize) === 1 ? "" : "s"}`;
    return {
      modelYear: String(snapshots.start || "").slice(0, 4) || "Not configured",
      resolution: resolutionValue,
    };
  }, [scenarioConfig]);

  const start = async (event: React.FormEvent) => {
    event.preventDefault();
    setSubmitting(true);
    try {
      const nextRun = await startModelRun({
        dataset: selection.dataset,
        project: selection.project,
        input_scenario: selection.scenario,
        output_scenario: outputScenario.trim(),
        cores,
      });
      setRun(nextRun);
      setError("");
    } catch (reason) {
      setError(reason instanceof Error ? reason.message : "Could not start the model run.");
    } finally {
      setSubmitting(false);
    }
  };

  const cancel = async () => {
    if (!run) return;
    try {
      setRun(await cancelModelRun(run.id));
      setError("");
    } catch (reason) {
      setError(reason instanceof Error ? reason.message : "Could not stop the model run.");
    }
  };

  if (loading) return <div className="editor-loading"><span className="spinner" />Reading base_config.yaml…</div>;

  return <>
    <PageHeader title={`Review & run ${selection.scenario}`} className="run-page-title" />

    {error && <div className="notice error run-notice"><AlertTriangle aria-hidden="true" />{error}</div>}

    <section className="run-review editor-panel">
      <header className="editor-panel-head"><div><h2>Run summary</h2><code>{scenarioConfig?.path || "scenario_config.yaml"}</code></div><CheckCircle2 aria-hidden="true" /></header>
      <div className="run-review-grid">
        <ReviewItem label="Data folder" value={selection.dataset} />
        <ReviewItem label="Project" value={selection.project} />
        <ReviewItem label="Input scenario" value={selection.scenario} />
        <ReviewItem label="Model year" value={scenarioSummary.modelYear} />
        <ReviewItem label="Resolution" value={scenarioSummary.resolution} />
        <ReviewItem label="Regions" value={options?.regions.join(", ") || "Not configured"} />
        <ReviewItem label="Sectors" value={options?.sectors.join(", ") || "Not configured"} />
        <ReviewItem label="Currency" value={options?.currency || "Not configured"} />
      </div>
      <div className="repository-review">
        <GitBranch aria-hidden="true" />
        {workspace?.repository.available ? <><div><span>Data repository</span><b>{workspace.repository.branch || "detached"} · {workspace.repository.commit || "unknown commit"}</b><small>{workspace.repository.remote || workspace.repository.root}</small></div><strong className={workspace.repository.dirty ? "dirty" : "clean"}>{workspace.repository.dirty ? `${workspace.repository.changes.length} local change${workspace.repository.changes.length === 1 ? "" : "s"}` : "Clean worktree"}</strong></> : <><div><span>Data repository</span><b>Local files only</b><small>No dataset-level Git repository was found.</small></div><strong>Not connected</strong></>}
      </div>
      <div className="run-requirements"><span><CheckCircle2 aria-hidden="true" /><b>{options?.environment || "hotpot"}</b> Conda environment required</span><span><FileCode2 aria-hidden="true" /><b>{options?.config_file || "base_config.yaml"}</b> copied per run</span><span><SquareTerminal aria-hidden="true" /><b>{options?.target || "solve_all_networks"}</b> target</span></div>
    </section>

    <section className="run-layout">
      <form className="run-setup editor-panel" onSubmit={start}>
        <header className="editor-panel-head">
          <div><h2>Run configuration</h2><code>{options?.config_file || "base_config.yaml"}</code></div>
          <FileCode2 aria-hidden="true" />
        </header>
        <div className="run-form">
          <p className="field-help">The values below replace the four <code>path_configs</code> entries in a run-local copy of <code>base_config.yaml</code>. The repository file remains the source template.</p>
          <div className="run-path-grid">
            <ReadOnlyField label="Data folder" value={selection.dataset} />
            <ReadOnlyField label="Project" value={selection.project} />
            <ReadOnlyField label="Input scenario" value={selection.scenario} />
            <label className="field"><span>Output scenario</span><input value={outputScenario} onChange={(event) => setOutputScenario(event.target.value)} required maxLength={255} disabled={active} /></label>
          </div>
          <div className="run-dimensions">
            {dimensions.map((item) => <div key={item.label}><span>{item.label}</span><b>{item.value}</b></div>)}
          </div>
          <div className="run-advanced">
            <label className="field"><span>CPU cores</span><input type="number" min={1} max={128} value={cores} onChange={(event) => setCores(Number(event.target.value))} disabled={active} /></label>
            <div><span>Workflow target</span><code>{options?.target || "solve_all_networks"}</code></div>
          </div>
        </div>
        <footer className="run-actions">
          <div><b>Snakemake writes results into</b><code>data/{selection.dataset}/{selection.project}/results/{outputScenario || "…"}</code></div>
          <button className="button primary" type="submit" disabled={active || submitting || !outputScenario.trim()}><Play aria-hidden="true" />{submitting ? "Starting…" : active ? "Run in progress" : "Run model"}</button>
        </footer>
      </form>

      <section className="run-monitor editor-panel" aria-live="polite">
        <header className="editor-panel-head">
          <div><p className="eyebrow">Workflow monitor</p><h2>{run ? statusLabel(run.status) : "Ready to run"}</h2>{run && <code>{run.id}</code>}</div>
          <StatusIcon run={run} />
        </header>
        {run ? <>
          <div className="run-progress-block">
            {!runMatchesSelection && <div className="run-selection-warning"><AlertTriangle aria-hidden="true" />This status belongs to {run.project} / {run.input_scenario}.</div>}
            <div className="run-progress-label"><span>{run.message}</span><b>{Math.round(run.progress)}%</b></div>
            <ChilliProgress value={run.progress} />
            <dl className="run-meta">
              <div><dt>Current rule</dt><dd>{run.current_rule || "Waiting for Snakemake"}</dd></div>
              <div><dt>Started</dt><dd>{readableTime(run.started_at || run.created_at)}</dd></div>
              <div><dt>Output</dt><dd>{run.output_scenario}</dd></div>
              <div><dt>Log file</dt><dd>{run.log_file}</dd></div>
              <div><dt>Manifest</dt><dd>{run.manifest_file || "Not recorded"}</dd></div>
            </dl>
          </div>
          <div className="run-log-head"><span><SquareTerminal aria-hidden="true" />Snakemake log</span>{active && <button className="button secondary" type="button" onClick={cancel} disabled={run.status === "canceling"}><CircleStop aria-hidden="true" />{run.status === "canceling" ? "Stopping…" : "Stop run"}</button>}</div>
          <pre className="run-log" ref={logRef}>{run.log || "Waiting for Snakemake output…"}</pre>
          {!active && <div className="run-complete-actions"><button className="button secondary" type="button" onClick={onEditConfiguration}><Pencil aria-hidden="true" />Edit configuration</button><button className="button secondary" type="button" onClick={() => setRun(null)}><RotateCcw aria-hidden="true" />Run again</button>{run.status === "succeeded" && <button className="button primary" type="button" onClick={() => onOpenResults(run.output_scenario, run.dataset, run.project)}><ExternalLink aria-hidden="true" />Open results</button>}</div>}
        </> : <div className="run-empty"><SquareTerminal aria-hidden="true" /><b>No web run yet</b><span>Confirm the output scenario and start the model to see progress and logs here.</span></div>}
      </section>
    </section>
  </>;
}

function ReviewItem({ label, value }: { label: string; value: string }) {
  return <div><span>{label}</span><b>{value}</b></div>;
}

function ReadOnlyField({ label, value }: { label: string; value: string }) {
  return <label className="field"><span>{label}</span><input value={value} readOnly aria-readonly="true" /></label>;
}

function ChilliProgress({ value }: { value: number }) {
  const progress = Math.max(0, Math.min(100, value));
  return <div className="chilli-progress" role="progressbar" aria-label="Model run progress" aria-valuemin={0} aria-valuemax={100} aria-valuenow={Math.round(progress)} aria-valuetext={`${Math.round(progress)} percent complete`}>
    {Array.from({ length: CHILLI_COUNT }, (_, index) => {
      const fill = Math.max(0, Math.min(1, progress / 100 * CHILLI_COUNT - index));
      return <span className="chilli-step" aria-hidden="true" key={index}>
        <img className="chilli-outline" src="/ui/favicon.svg" alt="" />
        <span className="chilli-fill" style={{ clipPath: `inset(0 ${100 - fill * 100}% 0 0)` }}><img src="/ui/favicon.svg" alt="" /></span>
      </span>;
    })}
  </div>;
}

function StatusIcon({ run }: { run: ModelRun | null }) {
  if (!run || run.status === "queued" || run.status === "running" || run.status === "canceling") return <Clock3 className={run ? "run-status-icon running" : "run-status-icon"} aria-hidden="true" />;
  if (run.status === "succeeded") return <CheckCircle2 className="run-status-icon succeeded" aria-hidden="true" />;
  return <AlertTriangle className="run-status-icon failed" aria-hidden="true" />;
}
