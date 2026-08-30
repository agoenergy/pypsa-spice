export type ResultRow = Record<string, string | number | null>;

export interface ChartDefinition {
  id: string;
  name: string;
  units?: string;
  summary: "sum" | "none";
  table_name: string;
  leg_col: string;
  fil_col?: string;
  type: "bar" | "grouped_bar" | "filtered_bar" | "area_share" | "hourly_bar" | "filtered_hourly_bar" | "hourly_line" | "hourly_dual";
  hourly: boolean;
  secondary_y_lab?: string[];
}

export interface Sector {
  name: string;
  years: string[];
}
export interface Scenario { name: string; sectors: Sector[] }
export interface Project { name: string; scenarios: Scenario[] }
export interface Dataset { name: string; projects: Project[] }
export interface Section {
  id: string;
  label: string;
  title: string;
  charts: ChartDefinition[];
}
export interface Mapping { label: string; color: string }
export interface Catalog { datasets: Dataset[]; sections: Section[]; mappings: Record<string, Mapping> }

export interface ChartResponse {
  rows: ResultRow[];
  dimensions: Record<string, string[]>;
  meta: {
    source_rows: number;
    returned_rows: number;
    sampled: boolean;
    files: number;
    available_start: string | null;
    available_end: string | null;
  };
}

export interface Selection {
  dataset: string;
  project: string;
  scenario: string;
  comparison: string;
  sector: string;
  year: string;
}

export interface InputTableDefinition {
  id: string;
  label: string;
  scope: "global" | "scenario";
  sector?: "power" | "industry" | "transport";
  timeseries?: boolean;
}

export interface InputProject {
  name: string;
  scenarios: string[];
  countries: string[];
  technologies: InputTechnology[];
}

export interface InputTechnology {
  id: string;
  label: string;
  sector: "power" | "industry" | "transport" | "other";
  classes: string[];
  carriers: string[];
}

export interface InputDataset { name: string; projects: InputProject[] }

export interface InputCatalog {
  table_query_version: number;
  datasets: InputDataset[];
  global_tables: InputTableDefinition[];
  sector_tables: Record<string, InputTableDefinition[]>;
}

export interface InputSelection {
  dataset: string;
  project: string;
  scenario: string;
}

export interface InputColumn {
  name: string;
  label: string;
  kind: "boolean" | "number" | "string";
  editable: boolean;
}

export type InputCell = string | number | boolean | null;
export type InputRow = Record<string, InputCell> & { __row_id: number };

export interface InputTableResponse {
  path: string;
  revision: string;
  columns: InputColumn[];
  rows: InputRow[];
  total_rows: number;
  total_filtered_rows: number;
  filter_column?: string;
  country_options: string[];
  filter_options: string[];
}

export interface ScenarioConfigResponse {
  path: string;
  revision: string;
  sections: Record<string, Record<string, unknown>>;
}

export type ScenarioDifferenceStatus = "changed" | "added" | "removed" | "ambiguous";

export interface ScenarioDifference {
  status: ScenarioDifferenceStatus;
  item: string;
  country: string;
  parameter: string;
  reference: unknown;
  comparison: unknown;
  delta: number | null;
}

export interface ScenarioDifferenceSection {
  id: string;
  label: string;
  category: "configuration" | "power" | "industry" | "transport";
  kind: "config" | "constraint" | "input";
  changes: ScenarioDifference[];
}

export interface ScenarioComparisonResponse {
  dataset: string;
  project: string;
  reference: string;
  comparison: string;
  global_inputs_shared: boolean;
  summary: {
    changes: number;
    groups: number;
    input_tables: number;
    configuration_groups: number;
  };
  countries: string[];
  sections: ScenarioDifferenceSection[];
}

export interface ModelRunOptions {
  config_file: string;
  target: string;
  defaults: {
    dataset: string;
    project: string;
    input_scenario: string;
    output_scenario: string;
  };
  years: string[];
  sectors: string[];
  regions: string[];
  currency: string;
  environment: string;
}

export interface DataRepositoryStatus {
  available: boolean;
  root: string;
  remote: string;
  branch: string;
  commit: string;
  dirty: boolean;
  changes: string[];
}

export interface ScenarioWorkspaceStatus {
  mutation_locked: boolean;
  active_run: { id: string; status: ModelRunStatus } | null;
  repository: DataRepositoryStatus;
}

export interface CreateScenarioRequest {
  dataset: string;
  project: string;
  source_scenario: string;
  new_scenario: string;
}

export interface CreatedScenario {
  dataset: string;
  project: string;
  scenario: string;
  source_scenario: string;
  path: string;
  repository: DataRepositoryStatus;
}

export type ModelRunStatus = "queued" | "running" | "canceling" | "succeeded" | "failed" | "canceled";

export interface ModelRun {
  id: string;
  status: ModelRunStatus;
  progress: number;
  message: string;
  current_rule: string | null;
  created_at: string;
  started_at: string | null;
  ended_at: string | null;
  exit_code: number | null;
  pid: number | null;
  dataset: string;
  project: string;
  input_scenario: string;
  output_scenario: string;
  cores: number;
  target: string;
  config_file: string;
  log_file: string;
  manifest_file: string;
  log: string;
}

export interface StartModelRunRequest {
  dataset: string;
  project: string;
  input_scenario: string;
  output_scenario: string;
  cores: number;
}

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
  format: "pypsa-spice-dashboard";
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

export type ViewMode = "home" | "outputs" | "inputs" | "configure" | "compare" | "dashboard";

export interface WorkspaceOption {
  value: string;
  label: string;
}

export interface WorkspaceInventoryItem {
  key: string;
  dataset: string;
  project: string;
  inputScenarios: string[];
  resultRuns: string[];
  countries: string[];
  dashboards: Pick<DashboardDefinition, "id" | "title">[];
}
