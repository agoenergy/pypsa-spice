import { useEffect, useMemo, useState, type ReactNode } from "react";
import {
  ArrowRight,
  BarChart3,
  FilePlus2,
  FileInput,
  FolderOpen,
  GitCompareArrows,
  LayoutDashboard,
  Settings2,
} from "lucide-react";
import { LocalDashboardStore, type DashboardDefinition } from "./dashboard";
import type { Catalog, InputCatalog } from "./types";

type HomeDestination = "inputs" | "configure" | "compare" | "outputs" | "dashboard";

interface WorkspaceInventoryItem {
  key: string;
  dataset: string;
  project: string;
  inputScenarios: string[];
  resultRuns: string[];
  countries: string[];
}

export function workspaceInventory(catalog: Catalog | null, inputCatalog: InputCatalog | null): WorkspaceInventoryItem[] {
  const workspaces = new Map<string, WorkspaceInventoryItem>();
  const ensureWorkspace = (dataset: string, project: string) => {
    const key = `${dataset}::${project}`;
    const existing = workspaces.get(key);
    if (existing) return existing;
    const created = { key, dataset, project, inputScenarios: [], resultRuns: [], countries: [] };
    workspaces.set(key, created);
    return created;
  };

  inputCatalog?.datasets.forEach((dataset) => dataset.projects.forEach((project) => {
    const workspace = ensureWorkspace(dataset.name, project.name);
    workspace.inputScenarios = project.scenarios;
    workspace.countries = project.countries;
  }));
  catalog?.datasets.forEach((dataset) => dataset.projects.forEach((project) => {
    ensureWorkspace(dataset.name, project.name).resultRuns = project.scenarios.map((scenario) => scenario.name);
  }));

  return [...workspaces.values()].sort((left, right) => left.project.localeCompare(right.project) || left.dataset.localeCompare(right.dataset));
}

export default function HomePage({
  catalog,
  inputCatalog,
  resultError,
  inputError,
  onNavigate,
}: {
  catalog: Catalog | null;
  inputCatalog: InputCatalog | null;
  resultError: string;
  inputError: string;
  onNavigate: (destination: HomeDestination) => void;
}) {
  const inventory = useMemo(() => workspaceInventory(catalog, inputCatalog), [catalog, inputCatalog]);
  const dashboardStore = useMemo(() => new LocalDashboardStore(), []);
  const [dashboards, setDashboards] = useState<DashboardDefinition[] | null>(null);
  const [dashboardError, setDashboardError] = useState("");
  const loading = !catalog && !inputCatalog && !resultError && !inputError;

  useEffect(() => {
    let current = true;
    dashboardStore.list()
      .then((summaries) => Promise.all(summaries.map((summary) => dashboardStore.get(summary.id))))
      .then((stored) => {
        if (!current) return;
        setDashboards(stored.filter((dashboard): dashboard is DashboardDefinition => Boolean(dashboard)));
        setDashboardError("");
      })
      .catch((reason) => {
        if (!current) return;
        setDashboards([]);
        setDashboardError(reason instanceof Error ? reason.message : "Local dashboards could not be opened.");
      });
    return () => { current = false; };
  }, [dashboardStore]);

  const openDashboard = (id: string) => {
    dashboardStore.setLastOpenedId(id);
    onNavigate("dashboard");
  };

  return <div className="home-page">
    <section className="home-guide-panel">
      <header className="home-panel-heading"><div><span className="eyebrow">Workflow</span><h2>How to use the workspace</h2></div></header>
      <ul className="home-workflow">
        <WorkflowStep icon={<FilePlus2 />} title="Build input data skeleton" description="Create the initial project and input-file structure." disabled />
        <WorkflowStep icon={<FileInput />} title="Inputs" description="Inspect technologies and edit the selected scenario's CSV inputs." onClick={() => onNavigate("inputs")} />
        <WorkflowStep icon={<GitCompareArrows />} title="Scenario differences" description="Review configuration and input differences before running the model." onClick={() => onNavigate("compare")} />
        <WorkflowStep icon={<Settings2 />} title="Configure & run" description="Set scenario parameters, review the model scope, and start Snakemake." onClick={() => onNavigate("configure")} />
        <WorkflowStep icon={<BarChart3 />} title="Results" description="Inspect Power, Industry, Transport, Emissions, and Costs." onClick={() => onNavigate("outputs")} />
        <WorkflowStep icon={<LayoutDashboard />} title="Dashboards" description="Open, create, and arrange saved result dashboards." onClick={() => onNavigate("dashboard")} />
      </ul>
    </section>

    <div className="home-primary">
      <section className="home-workspaces">
        <header className="home-panel-heading">
          <div><span className="eyebrow">Local data</span><h2>Available projects</h2></div>
          <span>{inventory.length} {inventory.length === 1 ? "project" : "projects"}</span>
        </header>

        {loading && <div className="home-loading"><span className="spinner" />Discovering projects and scenarios…</div>}
        {!loading && inventory.length === 0 && <div className="home-empty">
          <FolderOpen aria-hidden="true" />
          <b>No projects found</b>
          <span>{inputError || resultError || "Add a project under the local data directory, then refresh the workspace."}</span>
        </div>}
        {inventory.map((workspace) => <WorkspaceCard key={workspace.key} workspace={workspace} />)}
      </section>

      <DashboardList dashboards={dashboards} error={dashboardError} onOpen={openDashboard} onOpenWorkspace={() => onNavigate("dashboard")} />

      <p className="home-source-note"><b>Local files are the source of truth.</b> Input changes are written only when saved. Result files remain read-only.</p>
    </div>
  </div>;
}

function WorkspaceCard({ workspace }: { workspace: WorkspaceInventoryItem }) {
  return <article className="home-project-card">
    <header>
      <div className="home-project-title"><FolderOpen aria-hidden="true" /><div><h3>{workspace.project}</h3><span>{workspace.dataset}</span></div></div>
      {workspace.countries.length > 0 && <span className="home-project-country-count">{workspace.countries.length} {workspace.countries.length === 1 ? "country" : "countries"}</span>}
    </header>
    <div className="home-project-columns">
      <ScenarioList
        label="Input scenarios"
        emptyLabel="No editable scenarios found"
        items={workspace.inputScenarios}
        icon={<FileInput />}
        href={(scenario) => `?view=inputs&dataset=${encodeURIComponent(workspace.dataset)}&project=${encodeURIComponent(workspace.project)}&scenario=${encodeURIComponent(scenario)}`}
      />
      <ScenarioList
        label="Result runs"
        emptyLabel="No result runs found"
        items={workspace.resultRuns}
        icon={<BarChart3 />}
        href={(run) => `?section=power&dataset=${encodeURIComponent(workspace.dataset)}&project=${encodeURIComponent(workspace.project)}&run=${encodeURIComponent(run)}`}
      />
    </div>
  </article>;
}

function ScenarioList({ label, emptyLabel, items, icon, href }: { label: string; emptyLabel: string; items: string[]; icon: ReactNode; href: (item: string) => string }) {
  return <section className="home-scenario-list">
    <header><span>{icon}{label}</span><small>{items.length}</small></header>
    {items.length ? <div>{items.map((item) => <a href={href(item)} key={item}>{item}<ArrowRight aria-hidden="true" /></a>)}</div> : <p>{emptyLabel}</p>}
  </section>;
}

function WorkflowStep({ icon, title, description, onClick, disabled = false }: { icon: ReactNode; title: string; description: string; onClick?: () => void; disabled?: boolean }) {
  return <li>
    <button type="button" className="home-workflow-main" onClick={onClick} disabled={disabled}><span className="home-step-icon">{icon}</span><span><b>{title}{disabled && <em>Not implemented</em>}</b><small>{description}</small></span>{!disabled && <ArrowRight aria-hidden="true" />}</button>
  </li>;
}

function DashboardList({ dashboards, error, onOpen, onOpenWorkspace }: { dashboards: DashboardDefinition[] | null; error: string; onOpen: (id: string) => void; onOpenWorkspace: () => void }) {
  return <section className="home-dashboards">
    <header className="home-panel-heading">
      <div><span className="eyebrow">Saved locally in this browser</span><h2>Available dashboards</h2></div>
      <div className="home-dashboard-heading-actions">
        {dashboards && <span>{dashboards.length} {dashboards.length === 1 ? "dashboard" : "dashboards"}</span>}
        <button type="button" className="button secondary" onClick={onOpenWorkspace}><LayoutDashboard aria-hidden="true" />Open dashboard workspace</button>
      </div>
    </header>
    {dashboards === null && <div className="home-dashboard-state"><span className="spinner" />Opening saved dashboards…</div>}
    {dashboards?.length === 0 && <div className="home-dashboard-state"><LayoutDashboard aria-hidden="true" /><b>No saved dashboards</b><span>{error || "Open the dashboard workspace to create one from existing result figures."}</span><button type="button" className="button primary" onClick={onOpenWorkspace}>Create a dashboard</button></div>}
    {dashboards && dashboards.length > 0 && <div className="home-dashboard-list">
      {dashboards.map((dashboard) => {
        const chartCount = dashboard.rows.filter((row) => row.type === "chart").length;
        const headingCount = dashboard.rows.filter((row) => row.type === "heading").length;
        return <button type="button" key={dashboard.id} onClick={() => onOpen(dashboard.id)}>
          <span className="home-dashboard-icon"><LayoutDashboard aria-hidden="true" /></span>
          <span className="home-dashboard-copy"><b>{dashboard.title}</b><small>{dashboard.description || "No description"}</small></span>
          <span className="home-dashboard-project"><small>Result project</small><b>{dashboard.project}</b><em>{dashboard.dataset}</em></span>
          <span className="home-dashboard-contents"><small>Contents</small><b>{chartCount} {chartCount === 1 ? "chart" : "charts"}</b><em>{headingCount} {headingCount === 1 ? "heading" : "headings"}</em></span>
          <span className="home-dashboard-updated"><small>Updated</small><b>{formatDashboardDate(dashboard.updatedAt)}</b></span>
          <ArrowRight aria-hidden="true" />
        </button>;
      })}
    </div>}
  </section>;
}

function formatDashboardDate(value: string): string {
  const date = new Date(value);
  if (!Number.isFinite(date.getTime())) return "Unknown";
  return new Intl.DateTimeFormat(undefined, { day: "numeric", month: "short", year: "numeric" }).format(date);
}
