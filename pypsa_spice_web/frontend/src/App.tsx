import { useEffect, useMemo, useRef, useState } from "react";
import { ArrowLeftRight, BarChart3, CarFront, CircleDollarSign, Cloud, Factory, FileInput, GitCompareArrows, List, Menu, Moon, Plus, RefreshCw, Settings2, Sun, X, Zap } from "lucide-react";
import { getCatalog, getInputCatalog, getLatestModelRun } from "./api";
import ChartCard from "./ChartCard";
import DataDialog from "./DataDialog";
import InputEditor from "./InputEditor";
import NewScenarioDialog from "./NewScenarioDialog";
import PageHeader from "./PageHeader";
import ScenarioConfigEditor from "./ScenarioConfigEditor";
import ScenarioComparison from "./ScenarioComparison";
import { confirmDiscardChanges } from "./dirtyState";
import type { Catalog, ChartDefinition, InputCatalog, InputSelection, ModelRunStatus, ResultRow, Selection } from "./types";

type ViewMode = "outputs" | "inputs" | "configure" | "compare";
type WorkspaceOption = { value: string; label: string };

const WORKSPACE_SEPARATOR = "::";
const emptySelection: Selection = { dataset: "", project: "", scenario: "", comparison: "", sector: "", year: "" };
const emptyInputSelection: InputSelection = { dataset: "", project: "", scenario: "" };
const sectionIcons = { power: Zap, industry: Factory, transport: CarFront, emissions: Cloud, costs: CircleDollarSign };
const activeModelRunStatuses = new Set<ModelRunStatus>(["queued", "running", "canceling"]);

function locationParams() { return new URLSearchParams(window.location.search); }
function sectionFromLocation() { return locationParams().get("section")?.toLowerCase() || "power"; }
function viewFromLocation(): ViewMode { const value = locationParams().get("view"); return value === "inputs" || value === "configure" || value === "compare" ? value : "outputs"; }
function countryFromLocation() { return locationParams().get("country") || "ALL"; }
function workspaceValue(dataset: string, project: string) { return `${dataset}${WORKSPACE_SEPARATOR}${project}`; }
function splitWorkspace(value: string) { const [dataset, ...project] = value.split(WORKSPACE_SEPARATOR); return { dataset, project: project.join(WORKSPACE_SEPARATOR) }; }

function resolveOutputSelection(data: Catalog, current: Selection, params = locationParams()): Selection {
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

function resolveInputSelection(data: InputCatalog, current: InputSelection, params = locationParams()): InputSelection {
  const dataset = data.datasets.find((item) => item.name === (params.get("dataset") || current.dataset)) || data.datasets[0];
  const project = dataset.projects.find((item) => item.name === (params.get("project") || current.project)) || dataset.projects[0];
  const requestedScenario = params.get("scenario") || current.scenario;
  return { dataset: dataset.name, project: project.name, scenario: project.scenarios.includes(requestedScenario) ? requestedScenario : project.scenarios[0] || "" };
}

function workspaceOptions(datasets: { name: string; projects: { name: string }[] }[]): WorkspaceOption[] {
  const duplicateNames = new Set<string>();
  const seen = new Set<string>();
  datasets.flatMap((dataset) => dataset.projects).forEach((project) => { if (seen.has(project.name)) duplicateNames.add(project.name); seen.add(project.name); });
  return datasets.flatMap((dataset) => dataset.projects.map((project) => ({
    value: workspaceValue(dataset.name, project.name),
    label: duplicateNames.has(project.name) ? `${project.name} · ${dataset.name}` : project.name,
  })));
}

export default function App() {
  const [catalog, setCatalog] = useState<Catalog | null>(null);
  const [inputCatalog, setInputCatalog] = useState<InputCatalog | null>(null);
  const [selection, setSelection] = useState<Selection>(emptySelection);
  const [inputSelection, setInputSelection] = useState<InputSelection>(emptyInputSelection);
  const [inputComparison, setInputComparison] = useState(() => locationParams().get("compare") || "");
  const [view, setView] = useState<ViewMode>(viewFromLocation);
  const [sectionId, setSectionId] = useState(sectionFromLocation);
  const [country, setCountry] = useState(countryFromLocation);
  const [sidebarOpen, setSidebarOpen] = useState(false);
  const [error, setError] = useState("");
  const [inputError, setInputError] = useState("");
  const [inputRefresh, setInputRefresh] = useState(0);
  const [inspector, setInspector] = useState<{ chart: ChartDefinition; rows: ResultRow[]; sourceCount: number } | null>(null);
  const [dark, setDark] = useState(() => localStorage.getItem("spice-theme") === "dark");
  const [newScenarioOpen, setNewScenarioOpen] = useState(false);
  const [modelRunStatus, setModelRunStatus] = useState<ModelRunStatus | null>(null);

  const loadCatalog = async () => {
    try {
      const data = await getCatalog(true);
      if (!data.datasets.length) throw new Error("No result CSV folders were found under data/.");
      setCatalog(data); setError(""); setSelection((current) => resolveOutputSelection(data, current));
    } catch (reason) { setError(reason instanceof Error ? reason.message : "Could not load results."); }
  };
  const loadInputs = async (preferred?: InputSelection) => {
    try {
      const data = await getInputCatalog();
      if (!data.datasets.length) throw new Error("No input projects were found under data/.");
      setInputCatalog(data); setInputError("");
      setInputSelection((current) => preferred ? resolveInputSelection(data, preferred, new URLSearchParams()) : resolveInputSelection(data, current));
    } catch (reason) { setInputError(reason instanceof Error ? reason.message : "Could not load model inputs."); }
  };

  useEffect(() => { void loadCatalog(); void loadInputs(); }, []);
  useEffect(() => {
    let current = true;
    const refreshRunStatus = () => {
      getLatestModelRun()
        .then((run) => { if (current) setModelRunStatus(run && activeModelRunStatuses.has(run.status) ? run.status : null); })
        .catch(() => { /* Keep the last known status during a transient refresh failure. */ });
    };
    refreshRunStatus();
    const timer = window.setInterval(refreshRunStatus, 2000);
    return () => { current = false; window.clearInterval(timer); };
  }, []);
  useEffect(() => { document.documentElement.dataset.theme = dark ? "dark" : "light"; localStorage.setItem("spice-theme", dark ? "dark" : "light"); }, [dark]);
  useEffect(() => {
    const syncLocation = () => {
      const params = locationParams();
      setView(viewFromLocation()); setSectionId(sectionFromLocation()); setCountry(countryFromLocation());
      if (catalog) setSelection((current) => resolveOutputSelection(catalog, current, params));
      if (inputCatalog) {
        setInputSelection((current) => resolveInputSelection(inputCatalog, current, params));
        setInputComparison(params.get("compare") || "");
      }
    };
    window.addEventListener("popstate", syncLocation); return () => window.removeEventListener("popstate", syncLocation);
  }, [catalog, inputCatalog]);

  const dataset = catalog?.datasets.find((item) => item.name === selection.dataset);
  const project = dataset?.projects.find((item) => item.name === selection.project);
  const scenario = project?.scenarios.find((item) => item.name === selection.scenario);
  const sector = scenario?.sectors.find((item) => item.name === selection.sector);
  const sections = catalog?.sections || [];
  const section = sections.find((item) => item.id === sectionId) || sections[0];
  const charts = section?.charts || [];
  const inputDataset = inputCatalog?.datasets.find((item) => item.name === inputSelection.dataset);
  const inputProject = inputDataset?.projects.find((item) => item.name === inputSelection.project);
  const activeDatasetName = view === "outputs" ? selection.dataset : inputSelection.dataset;
  const activeProjectName = view === "outputs" ? selection.project : inputSelection.project;
  const configurationCountries = inputProject?.countries || [];
  const projectChoices = useMemo(() => workspaceOptions(view === "outputs" ? catalog?.datasets || [] : inputCatalog?.datasets || []), [view, catalog, inputCatalog]);

  useEffect(() => { if (!catalog || !section || section.id === sectionId) return; setSectionId(section.id); }, [catalog, section, sectionId]);
  useEffect(() => { if (country !== "ALL" && inputProject && !configurationCountries.includes(country)) setCountry("ALL"); }, [country, inputProject, configurationCountries]);
  useEffect(() => {
    if (!inputProject) return;
    const alternatives = inputProject.scenarios.filter((name) => name !== inputSelection.scenario);
    if (!alternatives.includes(inputComparison)) setInputComparison(alternatives[0] || "");
  }, [inputProject, inputSelection.scenario, inputComparison]);
  useEffect(() => {
    if (!(selection.project || inputSelection.project)) return;
    const params = new URLSearchParams();
    if (view !== "outputs") params.set("view", view); else params.set("section", section?.id || sectionId);
    if (activeDatasetName) params.set("dataset", activeDatasetName);
    if (activeProjectName) params.set("project", activeProjectName);
    if (inputSelection.scenario) params.set("scenario", inputSelection.scenario);
    if (selection.scenario) params.set("run", selection.scenario);
    if (view === "outputs" && selection.comparison) params.set("compare", selection.comparison);
    if (view === "compare" && inputComparison) params.set("compare", inputComparison);
    if (selection.sector) params.set("sector", selection.sector);
    if (view === "configure" && country !== "ALL") params.set("country", country);
    const next = `${window.location.pathname}?${params.toString()}`;
    if (`${window.location.pathname}${window.location.search}` !== next) window.history.replaceState(null, "", next);
  }, [view, section?.id, sectionId, activeDatasetName, activeProjectName, inputSelection.scenario, inputComparison, selection.scenario, selection.comparison, selection.sector, country]);

  const setOutputWorkspace = (datasetName: string, projectName: string) => {
    if (!catalog) return;
    const nextDataset = catalog.datasets.find((item) => item.name === datasetName);
    const nextProject = nextDataset?.projects.find((item) => item.name === projectName);
    if (!nextDataset || !nextProject) return;
    const nextScenario = nextProject.scenarios[0]; const nextSector = nextScenario.sectors[0];
    setSelection({ dataset: datasetName, project: projectName, scenario: nextScenario.name, comparison: "", sector: nextSector.name, year: nextSector.years[0] || "" });
  };
  const setInputWorkspace = (datasetName: string, projectName: string) => {
    if (!inputCatalog) return;
    const nextDataset = inputCatalog.datasets.find((item) => item.name === datasetName);
    const nextProject = nextDataset?.projects.find((item) => item.name === projectName);
    if (!nextDataset || !nextProject) return;
    setInputSelection({ dataset: datasetName, project: projectName, scenario: nextProject.scenarios[0] || "" });
    setInputComparison(nextProject.scenarios[1] || "");
  };
  const chooseWorkspace = (value: string) => {
    if (!confirmDiscardChanges()) return;
    const next = splitWorkspace(value);
    setOutputWorkspace(next.dataset, next.project);
    setInputWorkspace(next.dataset, next.project);
  };
  const chooseScenario = (name: string) => { if (!confirmDiscardChanges()) return; const nextScenario = project!.scenarios.find((item) => item.name === name)!; const nextSector = nextScenario.sectors[0]; setSelection((current) => ({ ...current, scenario: name, comparison: current.comparison === name ? "" : current.comparison, sector: nextSector.name, year: nextSector.years[0] || "" })); };
  const chooseInputScenario = (name: string) => {
    if (!confirmDiscardChanges()) return;
    setInputSelection((current) => ({ ...current, scenario: name }));
    if (name === inputComparison) {
      setInputComparison(inputProject?.scenarios.find((scenarioName) => scenarioName !== name) || "");
    }
  };
  const chooseCountry = (name: string) => { if (confirmDiscardChanges()) setCountry(name); };
  const chooseSector = (name: string) => { if (!confirmDiscardChanges()) return; const next = scenario!.sectors.find((item) => item.name === name)!; setSelection((current) => ({ ...current, sector: name, year: next.years[0] || "" })); };
  const chooseSection = (name: string) => { if (!confirmDiscardChanges()) return; setView("outputs"); setSectionId(name); setSidebarOpen(false); window.scrollTo({ top: 0, behavior: "smooth" }); };
  const chooseView = (next: ViewMode) => { if (!confirmDiscardChanges()) return; const params = locationParams(); if (next === "outputs") { params.delete("view"); params.set("section", section?.id || "power"); } else { params.set("view", next); params.delete("section"); } if (next !== "configure") params.delete("country"); if (next !== "outputs" && next !== "compare") params.delete("compare"); window.history.pushState(null, "", `?${params.toString()}`); setView(next); setSidebarOpen(false); window.scrollTo({ top: 0, behavior: "smooth" }); };
  const refreshCurrent = () => { if (!confirmDiscardChanges()) return; if (view === "outputs") void loadCatalog(); else { void loadInputs(); setInputRefresh((current) => current + 1); } };

  const outputReady = Boolean(catalog && selection.scenario && section && sector);
  const inputReady = Boolean(inputCatalog && inputSelection.dataset && inputSelection.project && inputSelection.scenario);
  const comparisonReady = Boolean(inputReady && inputComparison);
  const contextReady = view === "outputs" ? outputReady : inputReady;
  const activeWorkspace = workspaceValue(activeDatasetName, activeProjectName);

  return <div className="app-shell">
    <a className="skip-link" href="#workspace">Skip to workspace</a>
    <aside className={`sidebar ${sidebarOpen ? "open" : ""}`}>
      <div className="brand"><img src="/brand/pypsa-logo.svg" alt="PyPSA-SPICE" /></div>
      <div className="workspace-nav" role="navigation" aria-label="Workflow">
        <div className="sidebar-section">
          <a className={`sidebar-primary ${view === "inputs" ? "active" : ""}`} href="?view=inputs" onClick={(event) => { event.preventDefault(); chooseView("inputs"); }}><FileInput aria-hidden="true" />Inputs</a>
          {view === "inputs" && <div className="sidebar-submenu-slot" id="input-table-menu" />}
        </div>
        <div className="sidebar-section">
          <a className={`sidebar-primary ${view === "configure" ? "active" : ""}`} href="?view=configure" onClick={(event) => { event.preventDefault(); chooseView("configure"); }}><Settings2 aria-hidden="true" />Configure &amp; run</a>
          {view === "configure" && <div className="sidebar-submenu-slot" id="config-section-tabs" />}
        </div>
        <div className="sidebar-section">
          <a className={`sidebar-primary ${view === "compare" ? "active" : ""}`} href="?view=compare" onClick={(event) => { event.preventDefault(); chooseView("compare"); }}><GitCompareArrows aria-hidden="true" />Scenario differences</a>
        </div>
        <div className="sidebar-section">
          <a className={`sidebar-primary ${view === "outputs" ? "active" : ""}`} href={`?section=${section?.id || "power"}`} onClick={(event) => { event.preventDefault(); chooseView("outputs"); }}><BarChart3 aria-hidden="true" />Results</a>
          {view === "outputs" && section && <nav className="sidebar-submenu-list" aria-label="Result pages">{sections.map((item) => { const SectionIcon = sectionIcons[item.id as keyof typeof sectionIcons] || Zap; return <a className={`sidebar-submenu-item ${item.id === section.id ? "active" : ""}`} href={`?section=${item.id}`} key={item.id} aria-current={item.id === section.id ? "page" : undefined} onClick={(event) => { event.preventDefault(); chooseSection(item.id); }}><SectionIcon aria-hidden="true" /><b>{item.label}</b><small>{item.charts.length}</small></a>; })}</nav>}
        </div>
      </div>
      <div className="sidebar-foot"><div className="sidebar-statuses"><span className="sidebar-status-indicator" data-label="Local files connected" aria-label="Local files connected" tabIndex={0}><i className="status-dot" aria-hidden="true" /></span>{modelRunStatus && <span className="sidebar-status-indicator model-run" data-label={modelRunStatus === "queued" ? "Model run queued" : modelRunStatus === "canceling" ? "Stopping model run" : "Model running"} aria-label={modelRunStatus === "queued" ? "Model run queued" : modelRunStatus === "canceling" ? "Stopping model run" : "Model running"} role="status" aria-live="polite" tabIndex={0}><i className="status-dot" aria-hidden="true" /></span>}</div><a href="/docs" target="_blank">API</a><button onClick={() => setDark(!dark)} aria-label="Toggle dark mode">{dark ? <Sun aria-hidden="true" /> : <Moon aria-hidden="true" />}</button></div>
    </aside>
    <div className="scrim" onClick={() => setSidebarOpen(false)} />
    <div className="main-column">
      <header className="workspace-bar">
        <button className="menu" onClick={() => setSidebarOpen(!sidebarOpen)} aria-label="Open workspace navigation"><Menu aria-hidden="true" /></button>
        {contextReady && <div className="workspace-context">
          <ContextControl className="project-control" label="Project" value={activeWorkspace} onChange={chooseWorkspace} options={projectChoices} />
          {view === "outputs" ? <>
            <ContextControl className="scenario-control" label="Result run" value={selection.scenario} onChange={chooseScenario} options={project!.scenarios.map((item) => ({ value: item.name, label: item.name }))} />
            <ContextControl className={`compare-control ${selection.comparison ? "is-active" : ""}`} label="Compare with" value={selection.comparison} onChange={(comparison) => setSelection((current) => ({ ...current, comparison }))} options={[{ value: "", label: "No comparison" }, ...project!.scenarios.filter((item) => item.name !== selection.scenario).map((item) => ({ value: item.name, label: item.name }))]} />
          </> : view === "compare" ? <>
            <ContextControl className="scenario-control" label="Reference scenario" value={inputSelection.scenario} onChange={chooseInputScenario} options={inputProject!.scenarios.map((item) => ({ value: item, label: item }))} />
            <button className="comparison-swap" onClick={() => { const previous = inputSelection.scenario; setInputSelection((current) => ({ ...current, scenario: inputComparison })); setInputComparison(previous); }} aria-label="Swap reference and comparison scenarios" title="Swap scenarios" disabled={!inputComparison}><ArrowLeftRight aria-hidden="true" /></button>
            <ContextControl className="compare-control is-active" label="Comparison scenario" value={inputComparison} onChange={setInputComparison} options={inputProject!.scenarios.filter((item) => item !== inputSelection.scenario).map((item) => ({ value: item, label: item }))} />
          </> : <>
            <div className="scenario-context-group"><ContextControl className="scenario-control" label="Scenario" value={inputSelection.scenario} onChange={chooseInputScenario} options={inputProject!.scenarios.map((item) => ({ value: item, label: item }))} />{view === "configure" && <button className="context-add" onClick={() => setNewScenarioOpen(true)} aria-label="Create new scenario" title="Create new scenario"><Plus aria-hidden="true" /></button>}</div>
            {view === "inputs" && <div className="input-topbar-controls" id="input-topbar-controls" />}
            {view === "configure" && <ContextControl className="country-control" label="Country" value={country} onChange={chooseCountry} options={[{ value: "ALL", label: "All countries" }, ...configurationCountries.map((item) => ({ value: item, label: item }))]} />}
          </>}
        </div>}
        <div className="top-actions"><button className="button secondary" onClick={refreshCurrent} aria-label="Refresh data" title="Refresh data"><RefreshCw aria-hidden="true" /></button></div>
      </header>
      <main id="workspace">
        {view === "outputs" ? outputReady ? <><PageHeader title={section!.title}>{scenario!.sectors.length > 1 && <div className="page-controls"><ContextControl label="Sector run" value={selection.sector} onChange={chooseSector} options={scenario!.sectors.map((item) => ({ value: item.name, label: item.name }))} /></div>}</PageHeader><section className="analysis">{error && <div className="notice">{error}</div>}<div className={`chart-grid ${selection.comparison ? "comparison-active" : ""}`}>{charts.map((chart) => <ChartCard key={chart.id} chart={chart} selection={selection} years={sector!.years} mappings={catalog!.mappings} darkMode={dark} onInspect={(inspectedChart, rows, sourceCount) => setInspector({ chart: inspectedChart, rows, sourceCount })} />)}</div>{charts.length === 0 && <div className="no-results">No visualisations are configured for this section.</div>}</section><ResultsToc key={section!.id} charts={charts} /></> : <Boot message={error || "Discovering local results…"} /> : view === "compare" ? comparisonReady ? <ScenarioComparison key={`${inputRefresh}:${inputSelection.scenario}:${inputComparison}`} selection={inputSelection} comparison={inputComparison} /> : inputReady ? <div className="comparison-empty"><GitCompareArrows aria-hidden="true" /><b>Two scenarios are required</b><span>Create another scenario before opening the comparison workspace.</span></div> : <Boot message={inputError || "Discovering model inputs…"} /> : inputReady ? view === "inputs" ? <InputEditor key={inputRefresh} catalog={inputCatalog!} selection={inputSelection} onNavigate={() => setSidebarOpen(false)} /> : <ScenarioConfigEditor key={inputRefresh} selection={inputSelection} country={country} onNavigate={() => setSidebarOpen(false)} onOpenResults={(runName, datasetName, projectName) => { window.location.href = `/?section=power&dataset=${encodeURIComponent(datasetName)}&project=${encodeURIComponent(projectName)}&run=${encodeURIComponent(runName)}`; }} /> : <Boot message={inputError || "Discovering model inputs…"} />}
      </main>
    </div>
    {inspector && <DataDialog {...inspector} onClose={() => setInspector(null)} />}
    {newScenarioOpen && inputProject && <NewScenarioDialog selection={inputSelection} scenarios={inputProject.scenarios} onClose={() => setNewScenarioOpen(false)} onCreated={(created) => { setNewScenarioOpen(false); void loadInputs({ dataset: created.dataset, project: created.project, scenario: created.scenario }).then(() => setInputRefresh((current) => current + 1)); }} />}
  </div>;
}

function ContextControl({ className = "", label, value, onChange, options }: { className?: string; label: string; value: string; onChange: (value: string) => void; options: WorkspaceOption[] }) {
  return <label className={`context-control ${className}`.trim()}><span>{label}</span><select value={value} onChange={(event) => onChange(event.target.value)}>{options.map((option) => <option value={option.value} key={option.value || "empty"}>{option.label}</option>)}</select></label>;
}
function ResultsToc({ charts }: { charts: ChartDefinition[] }) {
  const [open, setOpen] = useState(false);
  const root = useRef<HTMLDivElement>(null);
  useEffect(() => {
    if (!open) return;
    const closeOutside = (event: PointerEvent) => { if (!root.current?.contains(event.target as Node)) setOpen(false); };
    const closeWithEscape = (event: KeyboardEvent) => { if (event.key === "Escape") setOpen(false); };
    document.addEventListener("pointerdown", closeOutside);
    document.addEventListener("keydown", closeWithEscape);
    return () => { document.removeEventListener("pointerdown", closeOutside); document.removeEventListener("keydown", closeWithEscape); };
  }, [open]);
  if (!charts.length) return null;
  return <div className={`results-toc ${open ? "open" : ""}`} ref={root}>
    {open && <nav className="results-toc-panel" id="results-figure-list" aria-label="Figures on this page">
      <header><h2>Figures</h2><button className="icon-button" onClick={() => setOpen(false)} aria-label="Close figure list"><X aria-hidden="true" /></button></header>
      <ol>{charts.map((chart, index) => <li key={chart.id}><a href={`#figure-${chart.id}`} onClick={() => setOpen(false)}><span>{String(index + 1).padStart(2, "0")}</span><b>{chart.name}</b></a></li>)}</ol>
    </nav>}
    <button className="results-toc-trigger" onClick={() => setOpen((current) => !current)} aria-label="Open figure list" aria-expanded={open} aria-controls="results-figure-list"><List aria-hidden="true" /></button>
  </div>;
}
function Boot({ message }: { message: string }) { return <div className="boot inline"><img src="/brand/pypsa-logo.svg" alt="PyPSA" /><span className="spinner" />{message}</div>; }
