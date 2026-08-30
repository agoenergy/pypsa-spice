import type { ReactNode } from "react";
import { ArrowRight, BarChart3, FileInput, FolderOpen, LayoutDashboard } from "lucide-react";
import type { DashboardDefinition, WorkspaceInventoryItem } from "../types";
import styles from "./WorkspaceCard.module.scss";

interface WorkspaceCardProps {
  workspace: WorkspaceInventoryItem;
  dashboardsLoading: boolean;
  dashboardError: string;
  onOpenDashboard: (id: string) => void;
  onOpenDashboardWorkspace: () => void;
}

export default function WorkspaceCard({
  workspace,
  dashboardsLoading,
  dashboardError,
  onOpenDashboard,
  onOpenDashboardWorkspace,
}: WorkspaceCardProps) {
  return <article className={styles["project-card"]}>
    <header>
      <div className={styles["project-title"]}>
        <FolderOpen aria-hidden="true" />
        <div><h3>{workspace.project}</h3><span>{workspace.dataset}</span></div>
      </div>
      {workspace.countries.length > 0 && <span className={styles["country-count"]}>
        {workspace.countries.length} {workspace.countries.length === 1 ? "country" : "countries"}
      </span>}
    </header>
    <div className={styles["project-columns"]}>
      <ScenarioList
        variant="inputs"
        label="Input scenarios"
        emptyLabel="No editable scenarios found"
        items={workspace.inputScenarios}
        icon={<FileInput />}
        href={(scenario) => `?view=inputs&dataset=${encodeURIComponent(workspace.dataset)}&project=${encodeURIComponent(workspace.project)}&scenario=${encodeURIComponent(scenario)}`}
      />
      <ScenarioList
        variant="results"
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

function ScenarioList({ variant, label, emptyLabel, items, icon, href }: {
  variant: "inputs" | "results";
  label: string;
  emptyLabel: string;
  items: string[];
  icon: ReactNode;
  href: (item: string) => string;
}) {
  return <section className={`${styles["scenario-list"]} ${styles[variant]}`}>
    <ListHeader icon={icon} label={label} count={items.length} />
    {items.length
      ? <div>{items.map((item) => <a href={href(item)} key={item}>{item}<ArrowRight aria-hidden="true" /></a>)}</div>
      : <p>{emptyLabel}</p>}
  </section>;
}

function ProjectDashboardList({ dashboards, loading, error, onOpen, onOpenWorkspace }: {
  dashboards: Pick<DashboardDefinition, "id" | "title">[];
  loading: boolean;
  error: string;
  onOpen: (id: string) => void;
  onOpenWorkspace: () => void;
}) {
  return <section className={`${styles["scenario-list"]} ${styles["dashboards"]}`}>
    <ListHeader icon={<LayoutDashboard />} label="Saved dashboards" count={dashboards.length} />
    {loading
      ? <p>Opening saved dashboards…</p>
      : dashboards.length
        ? <div>{dashboards.map((dashboard) => <button type="button" onClick={() => onOpen(dashboard.id)} key={dashboard.id}>{dashboard.title}<ArrowRight aria-hidden="true" /></button>)}</div>
        : <p>{error || "No saved dashboards found"}<button type="button" onClick={onOpenWorkspace}>Create one</button></p>}
  </section>;
}

function ListHeader({ icon, label, count }: { icon: ReactNode; label: string; count: number }) {
  return <header><span>{icon}{label}</span><small>{count}</small></header>;
}
