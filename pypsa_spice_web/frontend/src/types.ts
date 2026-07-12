export type ResultRow = Record<string, string | number | null>;

export interface ChartDefinition {
  id: string;
  name: string;
  units?: string;
  table_name: string;
  leg_col: string;
  fil_col?: string;
  type: "bar" | "filtered_bar" | "area_share" | "hourly_bar" | "filtered_hourly_bar" | "hourly_line" | "hourly_dual";
  hourly: boolean;
  primary_y_lab?: string[];
  secondary_y_lab?: string[];
}

export interface Sector {
  name: string;
  years: string[];
  chart_count: number;
  section_counts: Record<string, number>;
}
export interface Scenario { name: string; sectors: Sector[] }
export interface Project { name: string; scenarios: Scenario[] }
export interface Dataset { name: string; projects: Project[] }
export interface Section {
  id: string;
  label: string;
  icon: string;
  title: string;
  eyebrow: string;
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
  csv_name: string;
  identifier: string;
  filter_col: string;
  with_charts: boolean;
  timeseries?: boolean;
  editable: boolean;
  uniform_unit?: string;
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
  countries: string[];
}

export interface InputDataset { name: string; projects: InputProject[] }

export interface InputCatalog {
  datasets: InputDataset[];
  global_tables: InputTableDefinition[];
  sector_tables: Record<string, InputTableDefinition[]>;
  config_sections: string[];
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
  truncated: boolean;
  filter_column?: string;
  with_charts: boolean;
  timeseries: boolean;
}

export interface ScenarioConfigResponse {
  path: string;
  revision: string;
  sections: Record<string, Record<string, unknown>>;
}
