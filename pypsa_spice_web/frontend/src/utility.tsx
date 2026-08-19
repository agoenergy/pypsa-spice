import type { Catalog, InputCatalog, InputSelection, Selection } from "./types";

export const DASHBOARD_FORMAT = "pypsa-spice-dashboard";
export const DASHBOARD_SCHEMA_VERSION = 2;
const DASHBOARD_STORAGE_KEY = "pypsa-spice-dashboards-v2";
const LEGACY_DASHBOARD_STORAGE_KEY = "pypsa-spice-dashboards-v1";
const LAST_DASHBOARD_KEY = "pypsa-spice-last-dashboard";

export type DashboardMode = "scenario" | "difference";

export interface DashboardChartConfig {
  sectionId: string;
  chartId: string;
  sector: string;
  scenarios: string[];
  mode: DashboardMode;
  country: string;
  year?: string;
  filterValue?: string;
  startTime?: string;
  endTime?: string;
  customTitle?: string;
}

export interface DashboardHeadingRow {
  id: string;
  type: "heading";
  title: string;
}

export interface DashboardChartRow {
  id: string;
  type: "chart";
  chart: DashboardChartConfig;
}

export type DashboardRow = DashboardHeadingRow | DashboardChartRow;

export interface DashboardDefinition {
  schemaVersion: 2;
  id: string;
  title: string;
  description: string;
  dataset: string;
  project: string;
  rows: DashboardRow[];
  createdAt: string;
  updatedAt: string;
}

export interface DashboardSummary {
  id: string;
  title: string;
  chartCount: number;
  updatedAt: string;
}

export interface DashboardExport {
  format: typeof DASHBOARD_FORMAT;
  schemaVersion: 2;
  exportedAt: string;
  dashboard: DashboardDefinition;
}

export interface DashboardStore {
  list(): Promise<DashboardSummary[]>;
  get(id: string): Promise<DashboardDefinition | null>;
  save(dashboard: DashboardDefinition): Promise<void>;
  delete(id: string): Promise<void>;
  getLastOpenedId(): string | null;
  setLastOpenedId(id: string): void;
}

type StoredDashboards = {
  schemaVersion: 2;
  dashboards: Record<string, DashboardDefinition>;
};

function identifier(): string {
  if (typeof crypto !== "undefined" && typeof crypto.randomUUID === "function") return crypto.randomUUID();
  return `${Date.now().toString(36)}-${Math.random().toString(36).slice(2)}`;
}

export function dashboardIdentifier(): string {
  return identifier();
}

function isRecord(value: unknown): value is Record<string, unknown> {
  return typeof value === "object" && value !== null && !Array.isArray(value);
}

function stringValue(record: Record<string, unknown>, key: string, fallback = ""): string {
  return typeof record[key] === "string" ? record[key] as string : fallback;
}

function validIsoDate(value: string): boolean {
  return Boolean(value && Number.isFinite(new Date(value).getTime()));
}

function parseHeadingRow(value: Record<string, unknown>): DashboardHeadingRow {
  return {
    id: stringValue(value, "id") || identifier(),
    type: "heading",
    title: stringValue(value, "title").trim() || "Untitled heading",
  };
}

function parseChartConfig(value: unknown): DashboardChartConfig {
  if (!isRecord(value)) throw new Error("A dashboard chart entry is invalid.");
  const sectionId = stringValue(value, "sectionId");
  const chartId = stringValue(value, "chartId");
  const sector = stringValue(value, "sector");
  const scenarios = [...new Set(Array.isArray(value.scenarios)
    ? value.scenarios.filter((item): item is string => typeof item === "string" && Boolean(item))
    : [])].slice(0, 2);
  const mode = value.mode === "difference" ? "difference" : "scenario";
  if (!sectionId || !chartId || !sector) throw new Error("A dashboard chart is missing its source.");
  if (scenarios.length === 0) throw new Error("Every dashboard chart needs at least one scenario.");
  if (mode === "difference" && scenarios.length !== 2) {
    throw new Error("A difference chart must contain exactly two scenarios.");
  }
  return {
    sectionId,
    chartId,
    sector,
    scenarios,
    mode,
    country: stringValue(value, "country", "ALL"),
    ...(stringValue(value, "year") ? { year: stringValue(value, "year") } : {}),
    ...(stringValue(value, "filterValue") ? { filterValue: stringValue(value, "filterValue") } : {}),
    ...(stringValue(value, "startTime") ? { startTime: stringValue(value, "startTime") } : {}),
    ...(stringValue(value, "endTime") ? { endTime: stringValue(value, "endTime") } : {}),
    ...(stringValue(value, "customTitle") ? { customTitle: stringValue(value, "customTitle") } : {}),
  };
}

function parseSchema2Rows(value: unknown): DashboardRow[] {
  if (!Array.isArray(value)) throw new Error("The dashboard row list is invalid.");
  return value.map((row) => {
    if (!isRecord(row)) throw new Error("A dashboard row is invalid.");
    if (row.type === "heading") return parseHeadingRow(row);
    if (row.type === "chart") {
      return {
        id: stringValue(row, "id") || identifier(),
        type: "chart",
        chart: parseChartConfig(row.chart),
      };
    }
    throw new Error("A dashboard row has an unsupported type.");
  });
}

function migrateSchema1Rows(value: Record<string, unknown>): DashboardRow[] {
  if (!Array.isArray(value.items)) throw new Error("The dashboard chart list is invalid.");
  const charts = value.items.map((item) => ({
    config: parseChartConfig(item),
    rowId: isRecord(item) ? stringValue(item, "rowId") : "",
    headingId: isRecord(item) ? stringValue(item, "headingId") : "",
  }));
  const rows: DashboardRow[] = [];
  const assignedCharts = new Set<number>();

  if (Array.isArray(value.rows)) {
    for (const rawRow of value.rows) {
      if (!isRecord(rawRow)) throw new Error("A dashboard row is invalid.");
      if (rawRow.type === "heading") {
        rows.push(parseHeadingRow(rawRow));
        continue;
      }
      if (rawRow.type !== "chart" && rawRow.type !== "charts") continue;
      const rowId = stringValue(rawRow, "id");
      const matches = charts
        .map((chart, index) => ({ ...chart, index }))
        .filter((chart) => chart.rowId === rowId && !assignedCharts.has(chart.index));
      matches.forEach((chart, index) => {
        rows.push({
          id: index === 0 && rowId ? rowId : identifier(),
          type: "chart",
          chart: chart.config,
        });
        assignedCharts.add(chart.index);
      });
    }
  } else {
    const headings = Array.isArray(value.headingSections)
      ? value.headingSections.filter(isRecord).map(parseHeadingRow)
      : [];
    const headingIds = new Set(headings.map((heading) => heading.id));
    const appendCharts = (headingId: string) => {
      charts.forEach((chart, index) => {
        if (assignedCharts.has(index)) return;
        const validHeadingId = headingIds.has(chart.headingId) ? chart.headingId : "";
        if (validHeadingId !== headingId) return;
        rows.push({ id: identifier(), type: "chart", chart: chart.config });
        assignedCharts.add(index);
      });
    };
    appendCharts("");
    headings.forEach((heading) => {
      rows.push(heading);
      appendCharts(heading.id);
    });
  }

  charts.forEach((chart, index) => {
    if (!assignedCharts.has(index)) rows.push({ id: identifier(), type: "chart", chart: chart.config });
  });
  return rows;
}

export function parseDashboardDefinition(value: unknown): DashboardDefinition {
  if (!isRecord(value)) throw new Error("The dashboard configuration is not an object.");
  if (value.schemaVersion !== 1 && value.schemaVersion !== DASHBOARD_SCHEMA_VERSION) {
    throw new Error(`Dashboard schema version ${String(value.schemaVersion)} is not supported.`);
  }
  const title = stringValue(value, "title").trim();
  const dataset = stringValue(value, "dataset");
  const project = stringValue(value, "project");
  if (!title || !dataset || !project) throw new Error("The dashboard title, dataset, or project is missing.");
  const now = new Date().toISOString();
  const createdAt = stringValue(value, "createdAt");
  const updatedAt = stringValue(value, "updatedAt");
  return {
    schemaVersion: DASHBOARD_SCHEMA_VERSION,
    id: stringValue(value, "id") || identifier(),
    title,
    description: stringValue(value, "description"),
    dataset,
    project,
    rows: value.schemaVersion === DASHBOARD_SCHEMA_VERSION
      ? parseSchema2Rows(value.rows)
      : migrateSchema1Rows(value),
    createdAt: validIsoDate(createdAt) ? createdAt : now,
    updatedAt: validIsoDate(updatedAt) ? updatedAt : now,
  };
}

export function parseDashboardExport(text: string): DashboardDefinition {
  let value: unknown;
  try {
    value = JSON.parse(text);
  } catch {
    throw new Error("This file does not contain valid JSON.");
  }
  if (!isRecord(value) || value.format !== DASHBOARD_FORMAT) {
    throw new Error("This is not a PyPSA-SPICE dashboard configuration.");
  }
  if (value.schemaVersion !== 1 && value.schemaVersion !== DASHBOARD_SCHEMA_VERSION) {
    throw new Error(`Dashboard export version ${String(value.schemaVersion)} is not supported.`);
  }
  return parseDashboardDefinition(value.dashboard);
}

export function serializeDashboard(dashboard: DashboardDefinition): string {
  const exported: DashboardExport = {
    format: DASHBOARD_FORMAT,
    schemaVersion: DASHBOARD_SCHEMA_VERSION,
    exportedAt: new Date().toISOString(),
    dashboard,
  };
  return `${JSON.stringify(exported, null, 2)}\n`;
}

export function importedDashboardCopy(dashboard: DashboardDefinition): DashboardDefinition {
  const now = new Date().toISOString();
  return {
    ...dashboard,
    id: identifier(),
    rows: dashboard.rows.map((row) => row.type === "heading"
      ? { ...row, id: identifier() }
      : { ...row, id: identifier(), chart: { ...row.chart, scenarios: [...row.chart.scenarios] } }),
    createdAt: now,
    updatedAt: now,
  };
}

export function duplicateDashboard(dashboard: DashboardDefinition): DashboardDefinition {
  return {
    ...importedDashboardCopy(dashboard),
    title: `${dashboard.title} copy`,
  };
}

export function createDashboard(catalog: Catalog, title = "Untitled dashboard"): DashboardDefinition {
  const dataset = catalog.datasets[0];
  const project = dataset?.projects[0];
  if (!dataset || !project) throw new Error("No result project is available for a dashboard.");
  const now = new Date().toISOString();
  return {
    schemaVersion: DASHBOARD_SCHEMA_VERSION,
    id: identifier(),
    title,
    description: "",
    dataset: dataset.name,
    project: project.name,
    rows: [],
    createdAt: now,
    updatedAt: now,
  };
}

function emptyStorage(): StoredDashboards {
  return { schemaVersion: DASHBOARD_SCHEMA_VERSION, dashboards: {} };
}

export class LocalDashboardStore implements DashboardStore {
  constructor(private readonly storage: Storage = window.localStorage) {}

  private read(): StoredDashboards {
    const raw = this.storage.getItem(DASHBOARD_STORAGE_KEY)
      || this.storage.getItem(LEGACY_DASHBOARD_STORAGE_KEY);
    if (!raw) return emptyStorage();
    try {
      const parsed = JSON.parse(raw) as unknown;
      if (!isRecord(parsed) || !isRecord(parsed.dashboards)) return emptyStorage();
      const dashboards: Record<string, DashboardDefinition> = {};
      for (const value of Object.values(parsed.dashboards)) {
        try {
          const dashboard = parseDashboardDefinition(value);
          dashboards[dashboard.id] = dashboard;
        } catch {
          // Keep valid dashboards available if one stored draft is malformed.
        }
      }
      return { schemaVersion: DASHBOARD_SCHEMA_VERSION, dashboards };
    } catch {
      return emptyStorage();
    }
  }

  private write(value: StoredDashboards): void {
    this.storage.setItem(DASHBOARD_STORAGE_KEY, JSON.stringify(value));
  }

  async list(): Promise<DashboardSummary[]> {
    return Object.values(this.read().dashboards)
      .map((dashboard) => ({
        id: dashboard.id,
        title: dashboard.title,
        chartCount: dashboard.rows.filter((row) => row.type === "chart").length,
        updatedAt: dashboard.updatedAt,
      }))
      .sort((left, right) => right.updatedAt.localeCompare(left.updatedAt));
  }

  async get(id: string): Promise<DashboardDefinition | null> {
    return this.read().dashboards[id] || null;
  }

  async save(dashboard: DashboardDefinition): Promise<void> {
    const stored = this.read();
    stored.dashboards[dashboard.id] = parseDashboardDefinition(dashboard);
    this.write(stored);
  }

  async delete(id: string): Promise<void> {
    const stored = this.read();
    delete stored.dashboards[id];
    this.write(stored);
    if (this.getLastOpenedId() === id) this.storage.removeItem(LAST_DASHBOARD_KEY);
  }

  getLastOpenedId(): string | null {
    return this.storage.getItem(LAST_DASHBOARD_KEY);
  }

  setLastOpenedId(id: string): void {
    this.storage.setItem(LAST_DASHBOARD_KEY, id);
  }
}

export function dashboardFilename(title: string): string {
  const stem = title
    .normalize("NFKD")
    .replace(/[^\w\s-]/g, "")
    .trim()
    .replace(/[\s_]+/g, "-")
    .replace(/-+/g, "-")
    .toLowerCase();
  return `${stem || "dashboard"}.json`;
}

const dirtyEditors = new Set<string>();

export function setEditorDirty(id: string, dirty: boolean): void {
  if (dirty) dirtyEditors.add(id);
  else dirtyEditors.delete(id);
}

export function confirmDiscardChanges(): boolean {
  return dirtyEditors.size === 0 || window.confirm("You have unsaved changes. Discard them?");
}

export type ViewMode = "home" | "outputs" | "inputs" | "configure" | "compare" | "dashboard";
export type WorkspaceOption = { value: string; label: string };

const WORKSPACE_SEPARATOR = "::";

export function locationParams(): URLSearchParams {
  return new URLSearchParams(window.location.search);
}

export function sectionFromLocation(): string {
  return locationParams().get("section")?.toLowerCase() || "power";
}

export function viewFromLocation(): ViewMode {
  const params = locationParams();
  const value = params.get("view");
  if (value === "home" || value === "inputs" || value === "configure" || value === "compare" || value === "dashboard") return value;
  if (value === "outputs" || params.has("section") || params.has("run") || params.has("sector")) return "outputs";
  return "home";
}

export function countryFromLocation(): string {
  return locationParams().get("country") || "ALL";
}

export function workspaceValue(dataset: string, project: string): string {
  return `${dataset}${WORKSPACE_SEPARATOR}${project}`;
}

export function splitWorkspace(value: string): { dataset: string; project: string } {
  const [dataset, ...project] = value.split(WORKSPACE_SEPARATOR);
  return { dataset, project: project.join(WORKSPACE_SEPARATOR) };
}

export function resolveOutputSelection(data: Catalog, current: Selection, params = locationParams()): Selection {
  const dataset = data.datasets.find((item) => item.name === (params.get("dataset") || current.dataset)) || data.datasets[0];
  const project = dataset.projects.find((item) => item.name === (params.get("project") || current.project)) || dataset.projects[0];
  const scenario = project.scenarios.find((item) => item.name === (params.get("run") || current.scenario)) || project.scenarios[0];
  const sector = scenario.sectors.find((item) => item.name === (params.get("sector") || current.sector)) || scenario.sectors[0];
  const comparisonName = params.get("compare") || current.comparison;
  return {
    dataset: dataset.name,
    project: project.name,
    scenario: scenario.name,
    comparison: project.scenarios.some((item) => item.name === comparisonName && item.name !== scenario.name) ? comparisonName : "",
    sector: sector.name,
    year: sector.years.includes(current.year) ? current.year : sector.years[0] || "",
  };
}

export function resolveInputSelection(data: InputCatalog, current: InputSelection, params = locationParams()): InputSelection {
  const dataset = data.datasets.find((item) => item.name === (params.get("dataset") || current.dataset)) || data.datasets[0];
  const project = dataset.projects.find((item) => item.name === (params.get("project") || current.project)) || dataset.projects[0];
  const requestedScenario = params.get("scenario") || current.scenario;
  return { dataset: dataset.name, project: project.name, scenario: project.scenarios.includes(requestedScenario) ? requestedScenario : project.scenarios[0] || "" };
}

export function workspaceOptions(datasets: { name: string; projects: { name: string }[] }[]): WorkspaceOption[] {
  const duplicateNames = new Set<string>();
  const seen = new Set<string>();
  datasets.flatMap((dataset) => dataset.projects).forEach((project) => {
    if (seen.has(project.name)) duplicateNames.add(project.name);
    seen.add(project.name);
  });
  return datasets.flatMap((dataset) => dataset.projects.map((project) => ({
    value: workspaceValue(dataset.name, project.name),
    label: duplicateNames.has(project.name) ? `${project.name} · ${dataset.name}` : project.name,
  })));
}
