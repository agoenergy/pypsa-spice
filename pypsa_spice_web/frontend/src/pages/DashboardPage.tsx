import { useEffect, useMemo, useRef, useState } from "react";
import { createPortal } from "react-dom";
import {
  Copy,
  ChevronDown,
  ChevronUp,
  Download,
  FileUp,
  Heading2,
  LayoutDashboard,
  Plus,
  Trash2,
  X,
} from "lucide-react";
import "./DashboardPage.css";
import DashboardChartCard from "../components/DashboardChartCard";
import PageHeader from "../components/PageHeader";
import {
  LocalDashboardStore,
  createDashboard,
  dashboardFilename,
  dashboardIdentifier,
  duplicateDashboard,
  importedDashboardCopy,
  parseDashboardExport,
  serializeDashboard,
} from "../utility";
import type {
  Catalog,
  ChartDefinition,
  DashboardChartConfig,
  DashboardChartRow,
  DashboardDefinition,
  DashboardHeadingRow,
  DashboardRow,
  DashboardSummary,
  Project,
  ResultRow,
} from "../types";

interface Props {
  catalog: Catalog;
  darkMode: boolean;
  onInspect: (chart: ChartDefinition, rows: ResultRow[], sourceCount: number) => void;
}

type SaveStatus = "saved" | "saving" | "error";

function workspaceKey(dataset: string, project: string): string {
  return `${dataset}\u0000${project}`;
}

function splitWorkspaceKey(value: string): { dataset: string; project: string } {
  const [dataset, project] = value.split("\u0000");
  return { dataset, project };
}

function projectFor(catalog: Catalog, dashboard: DashboardDefinition): Project | null {
  return catalog.datasets
    .find((dataset) => dataset.name === dashboard.dataset)
    ?.projects.find((project) => project.name === dashboard.project) || null;
}

export default function DashboardPage({ catalog, darkMode, onInspect }: Props) {
  const store = useMemo(() => new LocalDashboardStore(), []);
  const [dashboard, setDashboard] = useState<DashboardDefinition | null>(null);
  const [summaries, setSummaries] = useState<DashboardSummary[]>([]);
  const [saveStatus, setSaveStatus] = useState<SaveStatus>("saved");
  const [storageError, setStorageError] = useState("");
  const [chartPickerOpen, setChartPickerOpen] = useState(false);
  const [topbarTarget, setTopbarTarget] = useState<HTMLElement | null>(null);
  const [imported, setImported] = useState<DashboardDefinition | null>(null);
  const [importError, setImportError] = useState("");
  const fileInput = useRef<HTMLInputElement>(null);
  const hydrated = useRef(false);
  const newHeadingId = useRef("");

  const refreshSummaries = async () => setSummaries(await store.list());

  useEffect(() => {
    setTopbarTarget(document.getElementById("dashboard-topbar-controls"));
  }, []);

  useEffect(() => {
    let current = true;
    const load = async () => {
      try {
        let available = await store.list();
        let selected = store.getLastOpenedId();
        let next = selected ? await store.get(selected) : null;
        if (!next && available.length) next = await store.get(available[0].id);
        if (!next) {
          next = createDashboard(catalog);
          await store.save(next);
          available = await store.list();
        }
        if (!current) return;
        setDashboard(next);
        setSummaries(available);
        store.setLastOpenedId(next.id);
        setSaveStatus("saved");
        hydrated.current = true;
      } catch (reason) {
        if (!current) return;
        setStorageError(reason instanceof Error ? reason.message : "Local dashboards could not be opened.");
        setSaveStatus("error");
      }
    };
    void load();
    return () => { current = false; };
  }, [catalog, store]);

  useEffect(() => {
    if (!dashboard || !hydrated.current) return;
    setSaveStatus("saving");
    const timer = window.setTimeout(() => {
      store.save(dashboard)
        .then(async () => {
          store.setLastOpenedId(dashboard.id);
          await refreshSummaries();
          setSaveStatus("saved");
          setStorageError("");
        })
        .catch((reason) => {
          setStorageError(reason instanceof Error ? reason.message : "The dashboard could not be saved locally.");
          setSaveStatus("error");
        });
    }, 500);
    return () => window.clearTimeout(timer);
  }, [dashboard, store]);

  if (!dashboard) {
    return <div className="dashboard-boot"><span className="spinner" />{storageError || "Opening local dashboards…"}</div>;
  }

  const project = projectFor(catalog, dashboard);
  const workspaceOptions = catalog.datasets.flatMap((dataset) => dataset.projects.map((candidate) => ({
    value: workspaceKey(dataset.name, candidate.name),
    label: `${candidate.name} · ${dataset.name}`,
  })));
  const updateDashboard = (changes: Partial<DashboardDefinition>) => {
    setDashboard((current) => current ? { ...current, ...changes, updatedAt: new Date().toISOString() } : current);
  };
  const chooseDashboard = async (id: string) => {
    const next = await store.get(id);
    if (!next) return;
    setChartPickerOpen(false);
    store.setLastOpenedId(id);
    setDashboard(next);
    setSaveStatus("saved");
  };
  const newDashboard = async () => {
    const next = createDashboard(catalog);
    await store.save(next);
    setChartPickerOpen(false);
    store.setLastOpenedId(next.id);
    setDashboard(next);
    await refreshSummaries();
  };
  const copyDashboard = async () => {
    const next = duplicateDashboard(dashboard);
    await store.save(next);
    store.setLastOpenedId(next.id);
    setDashboard(next);
    await refreshSummaries();
  };
  const removeDashboard = async () => {
    if (!window.confirm(`Delete “${dashboard.title}” from this browser?`)) return;
    await store.delete(dashboard.id);
    const available = await store.list();
    if (available.length) {
      const next = await store.get(available[0].id);
      setDashboard(next);
      if (next) store.setLastOpenedId(next.id);
    } else {
      const next = createDashboard(catalog);
      await store.save(next);
      store.setLastOpenedId(next.id);
      setDashboard(next);
    }
    await refreshSummaries();
  };
  const chooseProject = (value: string) => {
    const next = splitWorkspaceKey(value);
    if ((next.dataset === dashboard.dataset && next.project === dashboard.project)
      || (dashboard.rows.some((row) => row.type === "chart") && !window.confirm("Changing project removes all chart rows from this dashboard. Continue?"))) return;
    setChartPickerOpen(false);
    updateDashboard({
      ...next,
      rows: dashboard.rows.filter((row) => row.type === "heading"),
    });
  };
  const addHeadingRow = () => {
    const row: DashboardHeadingRow = {
      id: dashboardIdentifier(),
      type: "heading",
      title: "",
    };
    newHeadingId.current = row.id;
    updateDashboard({ rows: [...dashboard.rows, row] });
  };
  const addChartRow = () => {
    setChartPickerOpen(true);
  };
  const updateRow = (index: number, row: DashboardRow) => {
    updateDashboard({ rows: dashboard.rows.map((current, rowIndex) => rowIndex === index ? row : current) });
  };
  const moveRow = (index: number, direction: -1 | 1) => {
    const target = index + direction;
    if (target < 0 || target >= dashboard.rows.length) return;
    const rows = [...dashboard.rows];
    [rows[index], rows[target]] = [rows[target], rows[index]];
    updateDashboard({ rows });
  };
  const removeRow = (index: number) => {
    const row = dashboard.rows[index];
    const label = row.type === "heading" ? `heading “${row.title}”` : "chart row";
    if (!window.confirm(`Remove this ${label}?`)) return;
    updateDashboard({ rows: dashboard.rows.filter((_, rowIndex) => rowIndex !== index) });
  };
  const exportConfiguration = () => {
    const blob = new Blob([serializeDashboard(dashboard)], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = dashboardFilename(dashboard.title);
    link.click();
    URL.revokeObjectURL(url);
  };
  const readImport = async (file: File | undefined) => {
    if (!file) return;
    setImportError("");
    try {
      setImported(parseDashboardExport(await file.text()));
    } catch (reason) {
      setImported(null);
      setImportError(reason instanceof Error ? reason.message : "The dashboard configuration could not be read.");
    } finally {
      if (fileInput.current) fileInput.current.value = "";
    }
  };
  const confirmImport = async () => {
    if (!imported) return;
    const next = importedDashboardCopy(imported);
    await store.save(next);
    store.setLastOpenedId(next.id);
    setDashboard(next);
    setImported(null);
    await refreshSummaries();
  };
  const renderChartRow = (row: DashboardChartRow, rowIndex: number, resultProject: Project) => {
    const chart = catalog.sections.find((section) => section.id === row.chart.sectionId)?.charts.find((candidate) => candidate.id === row.chart.chartId) || null;
    return <DashboardChartCard
      config={row.chart}
      chart={chart}
      project={resultProject}
      datasetName={dashboard.dataset}
      projectName={dashboard.project}
      mappings={catalog.mappings}
      sections={catalog.sections}
      darkMode={darkMode}
      onChange={(next) => updateRow(rowIndex, { ...row, chart: next })}
      onInspect={onInspect}
      rowActions={<RowActions
        index={rowIndex}
        rowCount={dashboard.rows.length}
        label={`chart row ${rowIndex + 1}`}
        onMove={(direction) => moveRow(rowIndex, direction)}
        onRemove={() => removeRow(rowIndex)}
      />}
    />;
  };

  return <>
    {topbarTarget && createPortal(<div className="dashboard-topbar-content" aria-label="Dashboard controls">
      <label className="context-control dashboard-selector"><span>Dashboard</span><select value={dashboard.id} onChange={(event) => void chooseDashboard(event.target.value)}>{summaries.map((summary) => <option value={summary.id} key={summary.id}>{summary.title} · {summary.chartCount} charts</option>)}</select></label>
      <label className="context-control dashboard-project"><span>Result project</span><select value={workspaceKey(dashboard.dataset, dashboard.project)} onChange={(event) => chooseProject(event.target.value)}>{workspaceOptions.map((option) => <option value={option.value} key={option.value}>{option.label}</option>)}</select></label>
      <div className="dashboard-topbar-actions">
        <button className="button secondary" onClick={() => void newDashboard()} aria-label="New dashboard" title="New dashboard"><Plus aria-hidden="true" /></button>
        <button className="button secondary" onClick={() => void copyDashboard()} aria-label="Duplicate dashboard" title="Duplicate dashboard"><Copy aria-hidden="true" /></button>
        <button className="button secondary" onClick={() => fileInput.current?.click()} aria-label="Import dashboard" title="Import dashboard"><FileUp aria-hidden="true" /></button>
        <input className="sr-only" ref={fileInput} type="file" accept=".json,application/json" onChange={(event) => void readImport(event.target.files?.[0])} />
        <button className="button secondary" onClick={exportConfiguration} aria-label="Export dashboard" title="Export dashboard"><Download aria-hidden="true" /></button>
        <button className="button secondary danger" onClick={() => void removeDashboard()} aria-label="Delete dashboard" title="Delete dashboard"><Trash2 aria-hidden="true" /></button>
      </div>
    </div>, topbarTarget)}

    <PageHeader title="Custom dashboards">
      <div className="dashboard-save-state" role="status" data-status={saveStatus}>
        <i aria-hidden="true" />
        {saveStatus === "saving" ? "Saving…" : saveStatus === "error" ? "Save failed" : "Saved locally"}
      </div>
    </PageHeader>

    {storageError && <div className="notice error">{storageError}</div>}
    {importError && <div className="notice error">{importError}</div>}

    <section className="dashboard-heading">
      <label><span className="sr-only">Dashboard title</span><input className="dashboard-title-input" value={dashboard.title} maxLength={120} onChange={(event) => updateDashboard({ title: event.target.value || "Untitled dashboard" })} /></label>
      <label><span className="sr-only">Dashboard description</span><textarea value={dashboard.description} maxLength={600} rows={1} placeholder="Add a short description for this dashboard…" onChange={(event) => updateDashboard({ description: event.target.value })} /></label>
    </section>

    {!project && <div className="dashboard-project-missing"><LayoutDashboard aria-hidden="true" /><b>Result project unavailable</b><span>The saved project “{dashboard.project}” is not present in the current Results catalog. Select another project above to rebuild this dashboard.</span></div>}

    {project && dashboard.rows.length === 0 && <div className="dashboard-empty">
      <LayoutDashboard aria-hidden="true" />
      <b>Build the first row</b>
      <span>Add a full-width chart or start with a heading row.</span>
      <div className="dashboard-empty-actions"><button className="button primary" onClick={addChartRow}><Plus aria-hidden="true" />Add chart</button><button className="button secondary" onClick={addHeadingRow}><Heading2 aria-hidden="true" />Add heading row</button></div>
    </div>}

    {project && dashboard.rows.map((row, index) => row.type === "heading"
      ? <section className="dashboard-builder-row dashboard-heading-row" key={row.id}>
        <label className="dashboard-heading-row-field">
          <span className="sr-only">Section heading</span>
          <textarea
            value={row.title}
            maxLength={120}
            rows={1}
            aria-label="Section heading"
            placeholder="Type a section heading…"
            autoFocus={newHeadingId.current === row.id}
            onFocus={() => { newHeadingId.current = ""; }}
            onBlur={() => { if (!row.title.trim()) updateRow(index, { ...row, title: "Untitled heading" }); }}
            onChange={(event) => updateRow(index, { ...row, title: event.target.value })}
          />
        </label>
        <RowActions index={index} rowCount={dashboard.rows.length} label={`${row.title || "Untitled"} heading row`} onMove={(direction) => moveRow(index, direction)} onRemove={() => removeRow(index)} />
      </section>
      : <section className="dashboard-builder-row dashboard-chart-row" aria-label={`Chart row ${index + 1}`} key={row.id}>
        <div className="dashboard-row-content">{renderChartRow(row, index, project)}</div>
      </section>)}

    {project && dashboard.rows.length > 0 && <div className="dashboard-builder-footer"><span>Add the next row</span><button className="button primary" onClick={addChartRow}><Plus aria-hidden="true" />Chart</button><button className="button secondary" onClick={addHeadingRow}><Heading2 aria-hidden="true" />Heading row</button></div>}

    {chartPickerOpen && project && <AddChartDialog catalog={catalog} project={project} onClose={() => setChartPickerOpen(false)} onAdd={(chart) => {
      updateDashboard({
        rows: [...dashboard.rows, { id: dashboardIdentifier(), type: "chart", chart }],
      });
      setChartPickerOpen(false);
    }} />}
    {imported && <ImportDashboardDialog catalog={catalog} dashboard={imported} onClose={() => setImported(null)} onImport={() => void confirmImport()} />}
  </>;
}

function RowActions({ index, rowCount, label, onMove, onRemove }: {
  index: number;
  rowCount: number;
  label: string;
  onMove: (direction: -1 | 1) => void;
  onRemove: () => void;
}) {
  return <div className="dashboard-row-actions">
    <button className="icon-button" disabled={index === 0} onClick={() => onMove(-1)} title="Move row up" aria-label={`Move ${label} up`}><ChevronUp aria-hidden="true" /></button>
    <button className="icon-button" disabled={index === rowCount - 1} onClick={() => onMove(1)} title="Move row down" aria-label={`Move ${label} down`}><ChevronDown aria-hidden="true" /></button>
    <button className="icon-button danger" onClick={onRemove} title="Remove row" aria-label={`Remove ${label}`}><Trash2 aria-hidden="true" /></button>
  </div>;
}

function AddChartDialog({ catalog, project, onClose, onAdd }: {
  catalog: Catalog;
  project: Project;
  onClose: () => void;
  onAdd: (config: DashboardChartConfig) => void;
}) {
  const [sectionId, setSectionId] = useState(catalog.sections[0]?.id || "");
  const section = catalog.sections.find((candidate) => candidate.id === sectionId) || catalog.sections[0];
  const [chartId, setChartId] = useState(section?.charts[0]?.id || "");
  const chart = section?.charts.find((candidate) => candidate.id === chartId) || section?.charts[0];
  const sectors = [...new Set(project.scenarios.flatMap((scenario) => scenario.sectors.map((sector) => sector.name)))];
  const [sector, setSector] = useState(sectors[0] || "");
  const eligibleScenarios = project.scenarios.filter((scenario) => scenario.sectors.some((candidate) => candidate.name === sector));
  const [scenarios, setScenarios] = useState<string[]>(eligibleScenarios[0] ? [eligibleScenarios[0].name] : []);
  const [mode, setMode] = useState<DashboardChartConfig["mode"]>("scenario");
  const years = scenarioYearsForPicker(project, scenarios, sector);

  useEffect(() => {
    setChartId(section?.charts[0]?.id || "");
  }, [sectionId, section]);
  useEffect(() => {
    const valid = scenarios.filter((name) => eligibleScenarios.some((candidate) => candidate.name === name));
    setScenarios(valid.length ? valid : eligibleScenarios[0] ? [eligibleScenarios[0].name] : []);
  }, [sector]);

  const toggleScenario = (name: string) => {
    if (scenarios.includes(name)) {
      if (scenarios.length > 1) setScenarios(scenarios.filter((scenario) => scenario !== name));
      return;
    }
    if (scenarios.length >= 2) return;
    setScenarios([...scenarios, name]);
  };
  const changeMode = (next: DashboardChartConfig["mode"]) => {
    setMode(next);
    if (next === "difference") {
      const selected = [...scenarios];
      for (const candidate of eligibleScenarios) {
        if (selected.length >= 2) break;
        if (!selected.includes(candidate.name)) selected.push(candidate.name);
      }
      setScenarios(selected.slice(0, 2));
    }
  };
  const valid = Boolean(chart && sector && scenarios.length > 0 && (mode !== "difference" || scenarios.length === 2));
  const add = () => {
    if (!chart || !valid) return;
    onAdd({
      sectionId: section.id,
      chartId: chart.id,
      sector,
      scenarios,
      mode,
      country: "ALL",
      filterValue: "ALL",
      ...(chart.hourly && years[0] ? { year: years[0] } : {}),
    });
  };

  return <div className="dialog-backdrop" role="presentation" onMouseDown={(event) => { if (event.target === event.currentTarget) onClose(); }}>
    <section className="dialog dashboard-dialog" role="dialog" aria-modal="true" aria-labelledby="add-dashboard-chart-title">
      <header><h2 id="add-dashboard-chart-title">Choose a Results chart</h2><button className="icon-button" onClick={onClose} aria-label="Close"><X aria-hidden="true" /></button></header>
      <div className="dashboard-dialog-body">
        <div className="dashboard-dialog-grid">
          <label className="field"><span>Results section</span><select value={sectionId} onChange={(event) => setSectionId(event.target.value)}>{catalog.sections.map((candidate) => <option value={candidate.id} key={candidate.id}>{candidate.label}</option>)}</select></label>
          <label className="field"><span>Chart</span><select value={chart?.id || ""} onChange={(event) => setChartId(event.target.value)}>{section?.charts.map((candidate) => <option value={candidate.id} key={candidate.id}>{candidate.name}</option>)}</select></label>
          <label className="field"><span>Sector run</span><select value={sector} onChange={(event) => setSector(event.target.value)}>{sectors.map((candidate) => <option key={candidate}>{candidate}</option>)}</select></label>
          <label className="field"><span>Display</span><select value={mode} onChange={(event) => changeMode(event.target.value as DashboardChartConfig["mode"])}><option value="scenario">Selected scenarios</option><option value="difference" disabled={eligibleScenarios.length < 2}>Difference</option></select></label>
        </div>
        <fieldset className="dashboard-scenario-picker large"><legend>{mode === "difference" ? "Choose reference, then comparison" : "Choose one or two scenarios"}</legend>{eligibleScenarios.map((scenario) => {
          const checked = scenarios.includes(scenario.name);
          return <label key={scenario.name}><input type="checkbox" checked={checked} disabled={!checked && scenarios.length >= 2} onChange={() => toggleScenario(scenario.name)} /><span>{scenario.name}</span>{checked && mode === "difference" && <small>{scenarios.indexOf(scenario.name) === 0 ? "Reference" : "Comparison"}</small>}</label>;
        })}</fieldset>
        {mode === "difference" && scenarios.length === 2 && <div className="dashboard-difference-preview"><span>Difference calculation</span><b>{scenarios[1]} − {scenarios[0]}</b><button className="button secondary" onClick={() => setScenarios([scenarios[1], scenarios[0]])}>Swap</button></div>}
      </div>
      <footer><span>{chart?.hourly ? "Hourly chart" : "Yearly chart"} · {project.name}</span><button className="button secondary" onClick={onClose}>Cancel</button><button className="button primary" onClick={add} disabled={!valid}><Plus aria-hidden="true" />Use chart</button></footer>
    </section>
  </div>;
}

function ImportDashboardDialog({ catalog, dashboard, onClose, onImport }: {
  catalog: Catalog;
  dashboard: DashboardDefinition;
  onClose: () => void;
  onImport: () => void;
}) {
  const project = projectFor(catalog, dashboard);
  const charts = dashboard.rows.filter((row): row is DashboardChartRow => row.type === "chart");
  const availableScenarios = new Set(project?.scenarios.map((scenario) => scenario.name) || []);
  const missingScenarios = [...new Set(charts.flatMap((row) => row.chart.scenarios).filter((scenario) => !availableScenarios.has(scenario)))];
  const missingCharts = charts.filter((row) => !catalog.sections
    .find((section) => section.id === row.chart.sectionId)
    ?.charts.some((chart) => chart.id === row.chart.chartId));
  return <div className="dialog-backdrop" role="presentation" onMouseDown={(event) => { if (event.target === event.currentTarget) onClose(); }}>
    <section className="dialog dashboard-import-dialog" role="dialog" aria-modal="true" aria-labelledby="import-dashboard-title">
      <header><div><p className="eyebrow">Configuration import</p><h2 id="import-dashboard-title">{dashboard.title}</h2></div><button className="icon-button" onClick={onClose} aria-label="Close"><X aria-hidden="true" /></button></header>
      <div className="dashboard-import-summary">
        <dl><div><dt>Project</dt><dd>{dashboard.project} · {dashboard.dataset}</dd></div><div><dt>Charts</dt><dd>{charts.length}</dd></div><div><dt>Scenarios</dt><dd>{new Set(charts.flatMap((row) => row.chart.scenarios)).size}</dd></div></dl>
        {!project && <div className="notice">This result project is not available locally. The configuration can still be imported and repaired later.</div>}
        {project && (missingScenarios.length > 0 || missingCharts.length > 0) && <div className="notice">{missingScenarios.length} scenario references and {missingCharts.length} chart references are unavailable. They will remain in the imported configuration.</div>}
        <p>Importing creates a new local dashboard and does not overwrite an existing dashboard.</p>
      </div>
      <footer><button className="button secondary" onClick={onClose}>Cancel</button><button className="button primary" onClick={onImport}><FileUp aria-hidden="true" />Import dashboard</button></footer>
    </section>
  </div>;
}

function scenarioYearsForPicker(project: Project, scenarios: string[], sector: string): string[] {
  const years = new Set<string>();
  project.scenarios
    .filter((scenario) => scenarios.includes(scenario.name))
    .forEach((scenario) => scenario.sectors.find((candidate) => candidate.name === sector)?.years.forEach((year) => years.add(year)));
  return [...years].sort((left, right) => left.localeCompare(right, undefined, { numeric: true }));
}
