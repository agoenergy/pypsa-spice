import { useEffect, useMemo, useState, type Dispatch, type ReactNode, type SetStateAction } from "react";
import {
  ArrowLeftRight,
  Expand,
  Minimize2,
  Settings2,
  Table2,
} from "lucide-react";
import "./DashboardChartCard.css";
import { getChart } from "../api";
import type { DashboardChartConfig } from "../types";
import Plot, { buildDifferenceRows, ChartLegend, getLegendValues } from "./Plot";
import type { Catalog, ChartDefinition, ChartResponse, Project, ResultRow, Selection } from "../types";

interface Props {
  config: DashboardChartConfig;
  chart: ChartDefinition | null;
  project: Project;
  datasetName: string;
  projectName: string;
  mappings: Catalog["mappings"];
  sections: Catalog["sections"];
  darkMode: boolean;
  onChange: (config: DashboardChartConfig) => void;
  onInspect: (chart: ChartDefinition, rows: ResultRow[], sourceCount: number) => void;
  rowActions: ReactNode;
}

type LoadedSeries = { scenario: string; response: ChartResponse };
const LOADING_INDICATOR_DELAY_MS = 500;

function selectionFor(config: DashboardChartConfig, dataset: string, project: string, scenario: string, year: string): Selection {
  return {
    dataset,
    project,
    scenario,
    comparison: "",
    sector: config.sector,
    year,
  };
}

function scenarioYears(project: Project, scenarios: string[], sector: string): string[] {
  const values = new Set<string>();
  for (const scenarioName of scenarios) {
    const scenario = project.scenarios.find((candidate) => candidate.name === scenarioName);
    scenario?.sectors.find((candidate) => candidate.name === sector)?.years.forEach((year) => values.add(year));
  }
  return [...values].sort((left, right) => left.localeCompare(right, undefined, { numeric: true }));
}

export default function DashboardChartCard({
  config,
  chart,
  project,
  datasetName,
  projectName,
  mappings,
  sections,
  darkMode,
  onChange,
  onInspect,
  rowActions,
}: Props) {
  const [series, setSeries] = useState<LoadedSeries[]>([]);
  const [error, setError] = useState("");
  const [loading, setLoading] = useState(false);
  const [showLoadingIndicator, setShowLoadingIndicator] = useState(false);
  const [expanded, setExpanded] = useState(false);
  const [editing, setEditing] = useState(false);
  const [hiddenLegendValues, setHiddenLegendValues] = useState<Set<string>>(() => new Set());
  const years = useMemo(
    () => scenarioYears(project, config.scenarios, config.sector),
    [project, config.scenarios.join("\u0000"), config.sector],
  );
  const selectedYear = config.year && years.includes(config.year) ? config.year : years[0] || "";
  const eligibleScenarios = project.scenarios.filter((scenario) => (
    scenario.sectors.some((sector) => sector.name === config.sector)
  ));
  const eligibleScenarioNames = new Set(eligibleScenarios.map((scenario) => scenario.name));
  const missingScenarios = config.scenarios.filter((name) => !project.scenarios.some((scenario) => scenario.name === name));
  const incompatibleScenarios = config.scenarios.filter((name) => (
    !missingScenarios.includes(name) && !eligibleScenarioNames.has(name)
  ));
  const scenarioOptions = [
    ...missingScenarios.map((name) => ({ name, missing: true })),
    ...eligibleScenarios.map((scenario) => ({ name: scenario.name, missing: false })),
  ];
  const sectors = [...new Set(project.scenarios.flatMap((scenario) => scenario.sectors.map((sector) => sector.name)))];

  useEffect(() => {
    setHiddenLegendValues(new Set());
  }, [config.chartId, config.scenarios.join("\u0000"), config.sector]);

  useEffect(() => {
    if (!chart || missingScenarios.length || incompatibleScenarios.length || !config.scenarios.length) {
      setSeries([]);
      return;
    }
    if (chart.hourly && !selectedYear) {
      setSeries([]);
      setError("No model year is available for this hourly chart.");
      return;
    }
    const controller = new AbortController();
    setLoading(true);
    setError("");
    Promise.all(config.scenarios.map(async (scenario) => ({
      scenario,
      response: await getChart(
        chart,
        selectionFor(config, datasetName, projectName, scenario, selectedYear),
        scenario,
        config.country || "ALL",
        config.filterValue || "ALL",
        config.startTime || "",
        config.endTime || "",
        controller.signal,
      ),
    })))
      .then(setSeries)
      .catch((reason) => {
        if (reason.name !== "AbortError") setError(reason instanceof Error ? reason.message : "Could not load this chart.");
      })
      .finally(() => {
        if (!controller.signal.aborted) setLoading(false);
      });
    return () => controller.abort();
  }, [
    chart,
    project,
    datasetName,
    projectName,
    config.chartId,
    config.scenarios.join("\u0000"),
    config.sector,
    config.country,
    config.filterValue,
    config.startTime,
    config.endTime,
    selectedYear,
    missingScenarios.join("\u0000"),
    incompatibleScenarios.join("\u0000"),
  ]);

  useEffect(() => {
    if (!loading) {
      setShowLoadingIndicator(false);
      return;
    }
    const timer = window.setTimeout(() => setShowLoadingIndicator(true), LOADING_INDICATOR_DELAY_MS);
    return () => window.clearTimeout(timer);
  }, [loading]);

  if (!chart) {
    return <article className="dashboard-chart-card chart-card">
      <header className="chart-head">
        <div className="chart-title"><h3>{config.customTitle || config.chartId}</h3></div>
        <div className="dashboard-chart-toolbar">{rowActions}</div>
      </header>
      <div className="dashboard-source-missing"><b>Chart unavailable</b><span>The chart definition “{config.chartId}” is no longer in the Results catalog. Choose a replacement below.</span></div>
      <div className="dashboard-card-config dashboard-source-repair"><ChartSourceSelect config={config} sections={sections} onChange={onChange} /></div>
    </article>;
  }

  const countries = [...new Set(series.flatMap((entry) => entry.response.dimensions.country || []))];
  const filters = chart.fil_col
    ? [...new Set(series.flatMap((entry) => entry.response.dimensions[chart.fil_col!] || []))]
    : [];
  const legendValues = getLegendValues(chart, ...series.map((entry) => entry.response));
  const sourceCount = series.reduce((sum, entry) => sum + entry.response.meta.source_rows, 0);
  const displayedRows = config.mode === "difference" && series.length === 2
    ? buildDifferenceRows(series[0].response, series[1].response, chart, series[0].scenario, series[1].scenario)
    : series.flatMap((entry) => entry.response.rows.map((row) => ({ ...row, scenario: entry.scenario })));
  const hasRows = series.some((entry) => entry.response.rows.length > 0);
  const title = config.customTitle || chart.name;
  const update = (changes: Partial<DashboardChartConfig>) => onChange({ ...config, ...changes });
  const toggleScenario = (scenarioName: string) => {
    const selected = config.scenarios.includes(scenarioName);
    if (selected && config.scenarios.length === 1) return;
    if (!selected && config.scenarios.length >= 2) return;
    update({ scenarios: selected ? config.scenarios.filter((name) => name !== scenarioName) : [...config.scenarios, scenarioName] });
  };
  const changeMode = (mode: DashboardChartConfig["mode"]) => {
    if (mode === "scenario") {
      update({ mode });
      return;
    }
    const scenarios = [...config.scenarios];
    for (const candidate of eligibleScenarios) {
      if (scenarios.length >= 2) break;
      if (!scenarios.includes(candidate.name)) scenarios.push(candidate.name);
    }
    update({ mode, scenarios: scenarios.slice(0, 2) });
  };

  return <article className={`dashboard-chart-card chart-card ${expanded ? "expanded" : ""} ${series.length > 1 ? "comparing" : ""}`}>
    <header className="chart-head">
      <div className="chart-title">
        <small>{config.mode === "difference" && config.scenarios.length === 2 ? `${config.scenarios[1]} − ${config.scenarios[0]}` : config.scenarios.join(" · ")}</small>
        <h3>{title}</h3>
      </div>
      <div className="dashboard-chart-toolbar">
        <div className="dashboard-card-actions">
          <button className={editing ? "active" : ""} onClick={() => setEditing((current) => !current)} title="Configure chart" aria-label={`Configure ${title}`}><Settings2 aria-hidden="true" /></button>
          <button onClick={() => onInspect({ ...chart, name: title }, displayedRows, sourceCount)} disabled={!displayedRows.length} title="View source data" aria-label={`View ${title} source data`}><Table2 aria-hidden="true" /></button>
          <button onClick={() => setExpanded((current) => !current)} title={expanded ? "Close expanded chart" : "Expand chart"} aria-label={`${expanded ? "Close" : "Expand"} ${title}`}>{expanded ? <Minimize2 aria-hidden="true" /> : <Expand aria-hidden="true" />}</button>
        </div>
        {rowActions}
      </div>
    </header>

    {editing && <div className="dashboard-card-config">
      <label className="field"><span>Chart title</span><input value={config.customTitle || ""} placeholder={chart.name} onChange={(event) => update({ customTitle: event.target.value })} /></label>
      <label className="field"><span>Display</span><select value={config.mode} onChange={(event) => changeMode(event.target.value as DashboardChartConfig["mode"])}><option value="scenario">Scenarios</option><option value="difference" disabled={eligibleScenarios.length < 2}>Difference</option></select></label>
      <ChartSourceSelect config={config} sections={sections} onChange={onChange} />
      <label className="field"><span>Sector run</span><select value={config.sector} onChange={(event) => {
        const sector = event.target.value;
        const available = project.scenarios.filter((scenario) => scenario.sectors.some((candidate) => candidate.name === sector)).map((scenario) => scenario.name);
        const selected = config.scenarios.filter((scenario) => available.includes(scenario));
        const mode = config.mode === "difference" && available.length < 2 ? "scenario" : config.mode;
        update({ sector, mode, scenarios: selected.length ? selected.slice(0, 2) : available.slice(0, mode === "difference" ? 2 : 1), year: undefined, startTime: undefined, endTime: undefined });
      }}>{sectors.map((sector) => <option key={sector}>{sector}</option>)}</select></label>
      <fieldset className="dashboard-scenario-picker"><legend>{config.mode === "difference" ? "Reference and comparison" : "Scenarios (maximum two)"}</legend>{scenarioOptions.map((scenario) => {
        const checked = config.scenarios.includes(scenario.name);
        const disabled = !checked && config.scenarios.length >= 2;
        return <label key={scenario.name}><input type="checkbox" checked={checked} disabled={disabled} onChange={() => toggleScenario(scenario.name)} /><span>{scenario.name}</span>{scenario.missing && <small>Missing · deselect to replace</small>}</label>;
      })}</fieldset>
      {config.mode === "difference" && config.scenarios.length === 2 && <button className="button secondary dashboard-swap" onClick={() => update({ scenarios: [config.scenarios[1], config.scenarios[0]] })}><ArrowLeftRight aria-hidden="true" />Swap direction</button>}
      <label className="field"><span>Country</span><select value={config.country} onChange={(event) => update({ country: event.target.value })}><option value="ALL">All countries</option>{countries.map((country) => <option key={country}>{country}</option>)}</select></label>
      {chart.fil_col && <label className="field"><span>{chart.fil_col}</span><select value={config.filterValue || "ALL"} onChange={(event) => update({ filterValue: event.target.value })}><option value="ALL">All</option>{filters.map((filter) => <option key={filter}>{filter}</option>)}</select></label>}
      {chart.hourly && <label className="field"><span>Year</span><select value={selectedYear} onChange={(event) => update({ year: event.target.value, startTime: "", endTime: "" })}>{years.map((year) => <option key={year}>{year}</option>)}</select></label>}
      {chart.hourly && <><label className="field"><span>Start time</span><input type="datetime-local" value={config.startTime || ""} onChange={(event) => update({ startTime: event.target.value })} /></label><label className="field"><span>End time</span><input type="datetime-local" value={config.endTime || ""} onChange={(event) => update({ endTime: event.target.value })} /></label></>}
    </div>}

    {missingScenarios.length > 0 && <div className="dashboard-source-missing"><b>Scenario unavailable</b><span>{missingScenarios.join(", ")} could not be found in this project. Open chart settings to choose a replacement.</span></div>}
    {!missingScenarios.length && incompatibleScenarios.length > 0 && <div className="dashboard-source-missing"><b>Sector unavailable</b><span>{incompatibleScenarios.join(", ")} does not contain the selected sector “{config.sector}”.</span></div>}
    {!missingScenarios.length && !incompatibleScenarios.length && <div className={`chart-body dashboard-chart-body ${hasRows ? "with-plot" : ""} ${config.mode === "scenario" && series.length > 1 ? "dashboard-multi-layout" : ""}`} aria-busy={loading}>
      {showLoadingIndicator && !series.length && <div className="state"><span className="spinner" />Reading result tables…</div>}
      {!loading && error && <div className="state empty"><b>No chart data</b><span>{error}</span></div>}
      {!loading && !error && series.length > 0 && !hasRows && <div className="state empty"><b>No values in this result table</b><span>The selected scenarios contain no values for this chart.</span></div>}
      {!error && hasRows && config.mode === "scenario" && series.length === 1 && <Plot chart={chart} primary={series[0].response} comparison={null} primaryName={series[0].scenario} comparisonName="" mappings={mappings} darkMode={darkMode} expanded={expanded} hiddenLegendValues={hiddenLegendValues} onLegendToggle={(value) => toggleLegend(value, setHiddenLegendValues)} />}
      {!error && hasRows && config.mode === "scenario" && series.length > 1 && <>
        {series.map((entry) => <div className="scenario-plot" key={entry.scenario}><ScenarioHeading title={entry.scenario} /><Plot chart={chart} primary={entry.response} comparison={null} primaryName={entry.scenario} comparisonName="" mappings={mappings} darkMode={darkMode} expanded={expanded} legendValues={legendValues} showLegend={false} hiddenLegendValues={hiddenLegendValues} onLegendToggle={(value) => toggleLegend(value, setHiddenLegendValues)} /></div>)}
        <ChartLegend values={legendValues} mappings={mappings} hiddenValues={hiddenLegendValues} onToggle={(value) => toggleLegend(value, setHiddenLegendValues)} />
      </>}
      {!error && hasRows && config.mode === "difference" && series.length === 2 && <div className="difference-plot"><ScenarioHeading label="Difference" title={`${series[1].scenario} − ${series[0].scenario}`} /><Plot chart={chart} primary={series[0].response} comparison={series[1].response} primaryName={series[0].scenario} comparisonName={series[1].scenario} mappings={mappings} darkMode={darkMode} expanded={expanded} difference hiddenLegendValues={hiddenLegendValues} onLegendToggle={(value) => toggleLegend(value, setHiddenLegendValues)} /></div>}
      {chart.hourly && showLoadingIndicator && series.length > 0 && <div className="hourly-loading-overlay" role="status" aria-live="polite"><span className="spinner" aria-hidden="true" /><span>Updating chart…</span></div>}
    </div>}
  </article>;
}

function toggleLegend(value: string, setter: Dispatch<SetStateAction<Set<string>>>) {
  setter((current) => {
    const next = new Set(current);
    if (next.has(value)) next.delete(value);
    else next.add(value);
    return next;
  });
}

function ScenarioHeading({ label = "Scenario", title }: { label?: string; title: string }) {
  return <div className="scenario-label"><small>{label}</small><h4>{title}</h4></div>;
}

function ChartSourceSelect({ config, sections, onChange }: {
  config: DashboardChartConfig;
  sections: Catalog["sections"];
  onChange: (config: DashboardChartConfig) => void;
}) {
  const value = `${config.sectionId}\u0000${config.chartId}`;
  return <label className="field"><span>Results chart</span><select value={value} onChange={(event) => {
    const [sectionId, chartId] = event.target.value.split("\u0000");
    onChange({
      ...config,
      sectionId,
      chartId,
      country: "ALL",
      filterValue: "ALL",
      year: undefined,
      startTime: undefined,
      endTime: undefined,
    });
  }}>{sections.map((section) => <optgroup label={section.label} key={section.id}>{section.charts.map((chart) => <option value={`${section.id}\u0000${chart.id}`} key={chart.id}>{chart.name}</option>)}</optgroup>)}</select></label>;
}
