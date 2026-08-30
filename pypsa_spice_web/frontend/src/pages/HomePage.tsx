import { useEffect, useMemo, useState } from "react";
import {
  BarChart3,
  FilePlus2,
  FileInput,
  FolderOpen,
  GitCompareArrows,
  LayoutDashboard,
  Settings2,
} from "lucide-react";
import "./HomePage.css";
import { LocalDashboardStore } from "../utility";
import PageHeader from "../components/PageHeader";
import WorkspaceActionCard from "../components/WorkspaceActionCard";
import WorkspaceCard from "../components/WorkspaceCard";
import type { Catalog, DashboardDefinition, InputCatalog, WorkspaceInventoryItem } from "../types";

type HomeDestination = "inputs" | "configure" | "compare" | "outputs" | "dashboard";

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
        <WorkspaceActionCard icon={<FilePlus2 />} title="Build input data skeleton" description="Create the initial project and input-file structure." disabled />
        <WorkspaceActionCard icon={<FileInput />} title="Inputs" description="Inspect technologies and edit the selected scenario's CSV inputs." onClick={() => onNavigate("inputs")} />
        <WorkspaceActionCard icon={<GitCompareArrows />} title="Scenario differences" description="Review configuration and input differences before running the model." onClick={() => onNavigate("compare")} />
        <WorkspaceActionCard icon={<Settings2 />} title="Configure & run" description="Set scenario parameters, review the model scope, and start Snakemake." onClick={() => onNavigate("configure")} />
        <WorkspaceActionCard icon={<BarChart3 />} title="Results" description="Inspect Power, Industry, Transport, Emissions, and Costs." onClick={() => onNavigate("outputs")} />
        <WorkspaceActionCard icon={<LayoutDashboard />} title="Dashboards" description="Open, create, and arrange saved result dashboards." onClick={() => onNavigate("dashboard")} />
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
