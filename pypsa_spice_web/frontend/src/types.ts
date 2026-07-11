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
