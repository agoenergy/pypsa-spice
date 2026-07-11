import type { Catalog, ChartDefinition, ChartResponse, Selection } from "./types";

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
