import type { Catalog, ChartDefinition, ChartResponse, CreatedScenario, CreateScenarioRequest, InputCatalog, InputCell, InputSelection, InputTableDefinition, InputTableResponse, InputTechnology, ModelRun, ModelRunOptions, ScenarioComparisonResponse, ScenarioConfigResponse, ScenarioWorkspaceStatus, Selection, StartModelRunRequest } from "./types";

async function apiJson<T>(response: Response, fallback: string): Promise<T> {
  if (!response.ok) {
    const body = await response.json().catch(() => ({}));
    throw new Error(body.detail || fallback);
  }
  return response.json();
}

export async function getCatalog(refresh = false): Promise<Catalog> {
  const response = await fetch(`/api/catalog${refresh ? `?t=${Date.now()}` : ""}`);
  if (!response.ok) throw new Error("The local result catalog could not be read.");
  return response.json();
}

export function chartParams(
  chart: ChartDefinition,
  selection: Selection,
  scenario: string,
  country = "ALL",
  filterValue = "ALL",
  startTime = "",
  endTime = "",
): URLSearchParams {
  const params = new URLSearchParams({
    dataset: selection.dataset,
    project: selection.project,
    scenario,
    sector: selection.sector,
    table: chart.table_name,
    legend: chart.leg_col,
    hourly: String(chart.hourly),
    country,
  });
  if (chart.hourly) params.set("year", selection.year);
  if (chart.hourly && startTime) params.set("start_time", startTime);
  if (chart.hourly && endTime) params.set("end_time", endTime);
  if (chart.fil_col) {
    params.set("filter_column", chart.fil_col);
    params.set("filter_value", filterValue);
  }
  return params;
}

export async function getChart(
  chart: ChartDefinition,
  selection: Selection,
  scenario: string,
  country: string,
  filterValue: string,
  startTime: string,
  endTime: string,
  signal: AbortSignal,
): Promise<ChartResponse> {
  const response = await fetch(
    `/api/chart?${chartParams(chart, selection, scenario, country, filterValue, startTime, endTime)}`,
    { signal },
  );
  if (!response.ok) {
    const body = await response.json().catch(() => ({}));
    throw new Error(body.detail || `No data found for ${chart.name}.`);
  }
  return response.json();
}

export function downloadUrl(chart: ChartDefinition, selection: Selection): string {
  const params = chartParams(chart, selection, selection.scenario);
  ["legend", "country", "filter_column", "filter_value"].forEach((key) => params.delete(key));
  return `/api/download?${params}`;
}

export async function getInputCatalog(): Promise<InputCatalog> {
  const catalog = await apiJson<InputCatalog>(
    await fetch(`/api/input/catalog?t=${Date.now()}`),
    "The input catalog could not be read.",
  );
  if (catalog.table_query_version !== 2) {
    throw new Error("The local backend is out of date. Restart the app with ./run_web.sh.");
  }
  return catalog;
}

export interface InputTableQuery {
  technology?: InputTechnology;
  country?: string;
  filterValue?: string;
  query?: string;
  offset?: number;
  limit?: number;
}

function inputParams(selection: InputSelection, definition: InputTableDefinition, query: InputTableQuery = {}): URLSearchParams {
  const params = new URLSearchParams({
    dataset: selection.dataset,
    project: selection.project,
    scenario: selection.scenario,
    scope: definition.scope,
    sector: definition.sector || "power",
    table: definition.id,
  });
  if (query.technology) {
    params.set("technology", query.technology.id);
    query.technology.carriers.forEach((carrier) => params.append("technology_carrier", carrier));
  }
  if (query.country && query.country !== "ALL") params.set("country", query.country);
  if (query.filterValue && query.filterValue !== "ALL") params.set("filter_value", query.filterValue);
  if (query.query) params.set("query", query.query);
  if (query.offset) params.set("offset", String(query.offset));
  if (query.limit) params.set("limit", String(query.limit));
  return params;
}

export async function getInputTable(selection: InputSelection, definition: InputTableDefinition, query: InputTableQuery = {}, signal?: AbortSignal): Promise<InputTableResponse> {
  return apiJson(
    await fetch(`/api/input/table?${inputParams(selection, definition, query)}`, { signal }),
    `Could not read ${definition.label}.`,
  );
}

export async function saveInputTable(
  selection: InputSelection,
  definition: InputTableDefinition,
  revision: string,
  changes: { row: number; column: string; value: InputCell }[],
  query: InputTableQuery = {},
): Promise<InputTableResponse> {
  return apiJson(
    await fetch(`/api/input/table?${inputParams(selection, definition, query)}`, {
      method: "PUT",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ revision, changes }),
    }),
    `Could not save ${definition.label}.`,
  );
}

const FUEL_SUPPLY_TABLE: InputTableDefinition = {
  id: "Fuel_costs",
  label: "Fuel costs",
  scope: "scenario",
  sector: "power",
};

export async function getFuelSupplyTable(selection: InputSelection, signal?: AbortSignal): Promise<InputTableResponse> {
  const first = await getInputTable(selection, FUEL_SUPPLY_TABLE, { limit: 500 }, signal);
  const rows = [...first.rows];
  let offset = rows.length;
  while (offset < first.total_filtered_rows) {
    const page = await getInputTable(selection, FUEL_SUPPLY_TABLE, { offset, limit: 500 }, signal);
    rows.push(...page.rows);
    if (!page.rows.length) break;
    offset += page.rows.length;
  }
  return { ...first, rows };
}

export async function saveFuelSupplyLimits(
  selection: InputSelection,
  revision: string,
  changes: { row: number; column: string; value: InputCell }[],
): Promise<InputTableResponse> {
  return saveInputTable(selection, FUEL_SUPPLY_TABLE, revision, changes, { limit: 500 });
}

function configParams(selection: InputSelection): URLSearchParams {
  return new URLSearchParams({
    dataset: selection.dataset,
    project: selection.project,
    scenario: selection.scenario,
  });
}

export async function getScenarioConfig(selection: InputSelection, signal?: AbortSignal): Promise<ScenarioConfigResponse> {
  return apiJson(
    await fetch(`/api/input/scenario-config?${configParams(selection)}`, { signal }),
    "Could not read the scenario configuration.",
  );
}

export async function getScenarioComparison(
  selection: InputSelection,
  comparison: string,
  signal?: AbortSignal,
): Promise<ScenarioComparisonResponse> {
  const params = new URLSearchParams({
    dataset: selection.dataset,
    project: selection.project,
    reference: selection.scenario,
    comparison,
  });
  return apiJson(
    await fetch(`/api/input/compare?${params}`, { signal }),
    "Could not compare the selected scenarios.",
  );
}

export async function getScenarioWorkspaceStatus(dataset: string): Promise<ScenarioWorkspaceStatus> {
  return apiJson(
    await fetch(`/api/input/workspace?dataset=${encodeURIComponent(dataset)}&t=${Date.now()}`),
    "Could not read the local scenario workspace status.",
  );
}

export async function createScenario(request: CreateScenarioRequest): Promise<CreatedScenario> {
  return apiJson(
    await fetch("/api/input/scenarios", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(request),
    }),
    "Could not create the local scenario.",
  );
}

export async function saveScenarioConfigSection(
  selection: InputSelection,
  section: string,
  revision: string,
  value: Record<string, unknown>,
): Promise<ScenarioConfigResponse> {
  return apiJson(
    await fetch(`/api/input/scenario-config/${encodeURIComponent(section)}?${configParams(selection)}`, {
      method: "PUT",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ revision, value }),
    }),
    "Could not save the scenario configuration.",
  );
}

export async function saveScenarioConfigSections(
  selection: InputSelection,
  revision: string,
  sections: Record<string, Record<string, unknown>>,
): Promise<ScenarioConfigResponse> {
  const response = await fetch(`/api/input/scenario-config?${configParams(selection)}`, {
    method: "PUT",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ revision, sections }),
  });
  if (response.status !== 405) {
    return apiJson(response, "Could not save the scenario configuration.");
  }

  // Compatibility with a web server process started before the combined-save route
  // was added. Its established per-section endpoints still provide revision checks.
  let currentRevision = revision;
  let saved: ScenarioConfigResponse | null = null;
  for (const [section, value] of Object.entries(sections)) {
    saved = await saveScenarioConfigSection(selection, section, currentRevision, value);
    currentRevision = saved.revision;
  }
  if (!saved) throw new Error("No scenario configuration changes were provided.");
  return saved;
}

export async function getModelRunOptions(): Promise<ModelRunOptions> {
  return apiJson(
    await fetch(`/api/runs/options?t=${Date.now()}`),
    "Could not read base_config.yaml.",
  );
}

export async function getLatestModelRun(): Promise<ModelRun | null> {
  return apiJson(
    await fetch(`/api/runs/latest?t=${Date.now()}`),
    "Could not read the latest model run.",
  );
}

export async function getModelRun(runId: string): Promise<ModelRun> {
  return apiJson(
    await fetch(`/api/runs/${encodeURIComponent(runId)}?t=${Date.now()}`),
    "Could not refresh the model run.",
  );
}

export async function startModelRun(request: StartModelRunRequest): Promise<ModelRun> {
  return apiJson(
    await fetch("/api/runs", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(request),
    }),
    "Could not start the model run.",
  );
}

export async function cancelModelRun(runId: string): Promise<ModelRun> {
  return apiJson(
    await fetch(`/api/runs/${encodeURIComponent(runId)}`, { method: "DELETE" }),
    "Could not stop the model run.",
  );
}
