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
import "./HomePage.css";
import { LocalDashboardStore, type DashboardDefinition } from "../utility";
import PageHeader from "../components/PageHeader";
import type { Catalog, InputCatalog } from "../types";

type HomeDestination = "inputs" | "configure" | "compare" | "outputs" | "dashboard";

interface WorkspaceInventoryItem {
  key: string;
  dataset: string;
  project: string;
  inputScenarios: string[];
  resultRuns: string[];
  countries: string[];
  dashboards: Pick<DashboardDefinition, "id" | "title">[];
}

type DashboardWorkspaceSource = Pick<DashboardDefinition, "dataset" | "project" | "id" | "title">;

export function workspaceInventory(catalog: Catalog | null, inputCatalog: InputCatalog | null, dashboards: DashboardWorkspaceSource[] = []): WorkspaceInventoryItem[] {
  const workspaces = new Map<string, WorkspaceInventoryItem>();
  const ensureWorkspace = (dataset: string, project: string) => {
    const key = `${dataset}::${project}`;
    const existing = workspaces.get(key);
    if (existing) return existing;
    const created: WorkspaceInventoryItem = { key, dataset, project, inputScenarios: [], resultRuns: [], countries: [], dashboards: [] };
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
  dashboards.forEach((dashboard) => {
    ensureWorkspace(dashboard.dataset, dashboard.project).dashboards.push({ id: dashboard.id, title: dashboard.title });
  });

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
  const dashboardStore = useMemo(() => new LocalDashboardStore(), []);
  const [dashboards, setDashboards] = useState<DashboardDefinition[] | null>(null);
  const [dashboardError, setDashboardError] = useState("");
  const inventory = useMemo(() => workspaceInventory(catalog, inputCatalog, dashboards ?? []), [catalog, inputCatalog, dashboards]);
  const [selectedWorkspaceKey, setSelectedWorkspaceKey] = useState("");
  const loading = !catalog && !inputCatalog && !resultError && !inputError;
  const selectedWorkspace = inventory.find((workspace) => workspace.key === selectedWorkspaceKey) || inventory[0];

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
    <PageHeader title="Workspace overview" />
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
          {inventory.length > 0 && <label className="home-project-picker">
            <span>Project · {inventory.length} available</span>
            <select value={selectedWorkspace?.key ?? ""} onChange={(event) => setSelectedWorkspaceKey(event.target.value)}>
              {inventory.map((workspace) => <option value={workspace.key} key={workspace.key}>{workspace.project} · {workspace.dataset}</option>)}
            </select>
          </label>}
        </header>

        {loading && <div className="home-loading"><span className="spinner" />Discovering projects and scenarios…</div>}
        {!loading && inventory.length === 0 && <div className="home-empty">
          <FolderOpen aria-hidden="true" />
          <b>No projects found</b>
          <span>{inputError || resultError || "Add a project under the local data directory, then refresh the workspace."}</span>
        </div>}
        {!loading && selectedWorkspace && <WorkspaceCard
          key={selectedWorkspace.key}
          workspace={selectedWorkspace}
          dashboardsLoading={dashboards === null}
          dashboardError={dashboardError}
          onOpenDashboard={openDashboard}
          onOpenDashboardWorkspace={() => onNavigate("dashboard")}
        />}
      </section>

      <p className="home-source-note"><b>Local files are the source of truth.</b> Input changes are written only when saved. Result files remain read-only.</p>
    </div>
  </div>;
}

function WorkspaceCard({ workspace, dashboardsLoading, dashboardError, onOpenDashboard, onOpenDashboardWorkspace }: {
  workspace: WorkspaceInventoryItem;
  dashboardsLoading: boolean;
  dashboardError: string;
  onOpenDashboard: (id: string) => void;
  onOpenDashboardWorkspace: () => void;
}) {
  return <article className="home-project-card">
    <header>
      <div className="home-project-title"><FolderOpen aria-hidden="true" /><div><h3>{workspace.project}</h3><span>{workspace.dataset}</span></div></div>
      {workspace.countries.length > 0 && <span className="home-project-country-count">{workspace.countries.length} {workspace.countries.length === 1 ? "country" : "countries"}</span>}
    </header>
    <div className="home-project-columns">
      <ScenarioList
        className="home-input-scenarios"
        label="Input scenarios"
        emptyLabel="No editable scenarios found"
        items={workspace.inputScenarios}
        icon={<FileInput />}
        href={(scenario) => `?view=inputs&dataset=${encodeURIComponent(workspace.dataset)}&project=${encodeURIComponent(workspace.project)}&scenario=${encodeURIComponent(scenario)}`}
      />
      <ScenarioList
        className="home-result-runs"
        label="Result runs"
        emptyLabel="No result runs found"
        items={workspace.resultRuns}
        icon={<BarChart3 />}
        href={(run) => `?section=power&dataset=${encodeURIComponent(workspace.dataset)}&project=${encodeURIComponent(workspace.project)}&run=${encodeURIComponent(run)}`}
      />
      <ProjectDashboardList
        dashboards={workspace.dashboards}
        loading={dashboardsLoading}
        error={dashboardError}
        onOpen={onOpenDashboard}
        onOpenWorkspace={onOpenDashboardWorkspace}
      />
    </div>
  </article>;
}

function ScenarioList({ className, label, emptyLabel, items, icon, href }: { className: "home-input-scenarios" | "home-result-runs"; label: string; emptyLabel: string; items: string[]; icon: ReactNode; href: (item: string) => string }) {
  return <section className={`home-scenario-list ${className}`}>
    <header><span>{icon}{label}</span><small>{items.length}</small></header>
    {items.length ? <div>{items.map((item) => <a href={href(item)} key={item}>{item}<ArrowRight aria-hidden="true" /></a>)}</div> : <p>{emptyLabel}</p>}
  </section>;
}

function ProjectDashboardList({ dashboards, loading, error, onOpen, onOpenWorkspace }: {
  dashboards: Pick<DashboardDefinition, "id" | "title">[];
  loading: boolean;
  error: string;
  onOpen: (id: string) => void;
  onOpenWorkspace: () => void;
}) {
  return <section className="home-scenario-list home-project-dashboards">
    <header><span><LayoutDashboard aria-hidden="true" />Saved dashboards</span><small>{dashboards.length}</small></header>
    {loading ? <p>Opening saved dashboards…</p> : dashboards.length ? <div>{dashboards.map((dashboard) => <button type="button" onClick={() => onOpen(dashboard.id)} key={dashboard.id}>{dashboard.title}<ArrowRight aria-hidden="true" /></button>)}</div> : <p>{error || "No saved dashboards found"}<button type="button" onClick={onOpenWorkspace}>Create one</button></p>}
  </section>;
}

function WorkflowStep({ icon, title, description, onClick, disabled = false }: { icon: ReactNode; title: string; description: string; onClick?: () => void; disabled?: boolean }) {
  return <li>
    <button type="button" className="home-workflow-main" onClick={onClick} disabled={disabled}><span className="home-step-icon">{icon}</span><span><b>{title}{disabled && <em>Not implemented</em>}</b><small>{description}</small></span>{!disabled && <ArrowRight aria-hidden="true" />}</button>
  </li>;
}
