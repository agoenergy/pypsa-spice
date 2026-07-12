import { useEffect, useMemo, useState } from "react";
import { BarChart3, CarFront, CircleDollarSign, Cloud, Factory, FileInput, Menu, Moon, Play, RefreshCw, Search, Settings2, Sun, Zap } from "lucide-react";
import { getCatalog, getInputCatalog } from "./api";
import ChartCard from "./ChartCard";
import DataDialog from "./DataDialog";
import InputEditor from "./InputEditor";
import ScenarioConfigEditor from "./ScenarioConfigEditor";
import type { Catalog, InputCatalog, InputSelection, ResultRow, Selection } from "./types";

type ViewMode = "outputs" | "inputs" | "configure";
const emptySelection: Selection = { dataset: "", project: "", scenario: "", comparison: "", sector: "", year: "" };
const emptyInputSelection: InputSelection = { dataset: "", project: "", scenario: "" };
const sectionIcons = { power: Zap, industry: Factory, transport: CarFront, emissions: Cloud, costs: CircleDollarSign };

function sectionFromLocation() { return new URLSearchParams(window.location.search).get("section")?.toLowerCase() || "power"; }
function viewFromLocation(): ViewMode { const value = new URLSearchParams(window.location.search).get("view"); return value === "inputs" || value === "configure" ? value : "outputs"; }

export default function App() {
  const [catalog, setCatalog] = useState<Catalog | null>(null);
  const [inputCatalog, setInputCatalog] = useState<InputCatalog | null>(null);
  const [selection, setSelection] = useState<Selection>(emptySelection);
  const [inputSelection, setInputSelection] = useState<InputSelection>(emptyInputSelection);
  const [view, setView] = useState<ViewMode>(viewFromLocation);
  const [sectionId, setSectionId] = useState(sectionFromLocation);
  const [search, setSearch] = useState("");
  const [sidebarOpen, setSidebarOpen] = useState(false);
  const [error, setError] = useState("");
  const [inputError, setInputError] = useState("");
  const [inputRefresh, setInputRefresh] = useState(0);
  const [inspector, setInspector] = useState<{ title: string; rows: ResultRow[]; sourceCount: number } | null>(null);
  const [dark, setDark] = useState(() => localStorage.getItem("spice-theme") === "dark");

  const loadCatalog = async () => {
    try {
      const data = await getCatalog(true);
      if (!data.datasets.length) throw new Error("No result CSV folders were found under data/.");
      setCatalog(data); setError("");
      const dataset = data.datasets.find((item) => item.name === selection.dataset) || data.datasets[0];
      const project = dataset.projects.find((item) => item.name === selection.project) || dataset.projects[0];
      const scenario = project.scenarios.find((item) => item.name === selection.scenario) || project.scenarios[0];
      const sector = scenario.sectors.find((item) => item.name === selection.sector) || scenario.sectors[0];
      setSelection((current) => ({ ...current, dataset: dataset.name, project: project.name, scenario: scenario.name, sector: sector.name, year: sector.years.includes(current.year) ? current.year : sector.years[0] || "", comparison: project.scenarios.some((item) => item.name === current.comparison) ? current.comparison : "" }));
    } catch (reason) { setError(reason instanceof Error ? reason.message : "Could not load results."); }
  };
  const loadInputs = async () => {
    try {
      const data = await getInputCatalog();
      if (!data.datasets.length) throw new Error("No input projects were found under data/.");
      setInputCatalog(data); setInputError("");
      const dataset = data.datasets.find((item) => item.name === inputSelection.dataset) || data.datasets[0];
      const project = dataset.projects.find((item) => item.name === inputSelection.project) || dataset.projects[0];
      const scenario = project.scenarios.includes(inputSelection.scenario) ? inputSelection.scenario : project.scenarios[0] || "";
      setInputSelection({ dataset: dataset.name, project: project.name, scenario });
    } catch (reason) { setInputError(reason instanceof Error ? reason.message : "Could not load model inputs."); }
  };

  useEffect(() => { void loadCatalog(); void loadInputs(); }, []);
  useEffect(() => { document.documentElement.dataset.theme = dark ? "dark" : "light"; localStorage.setItem("spice-theme", dark ? "dark" : "light"); }, [dark]);
  useEffect(() => {
    const syncLocation = () => { setView(viewFromLocation()); setSectionId(sectionFromLocation()); };
    window.addEventListener("popstate", syncLocation); return () => window.removeEventListener("popstate", syncLocation);
  }, []);

  const dataset = catalog?.datasets.find((item) => item.name === selection.dataset);
  const project = dataset?.projects.find((item) => item.name === selection.project);
  const scenario = project?.scenarios.find((item) => item.name === selection.scenario);
  const sector = scenario?.sectors.find((item) => item.name === selection.sector);
  const sections = catalog?.sections || [];
  const section = sections.find((item) => item.id === sectionId) || sections[0];
  const charts = useMemo(() => section?.charts.filter((chart) => chart.name.toLowerCase().includes(search.toLowerCase().trim())) || [], [section, search]);
  const inputDataset = inputCatalog?.datasets.find((item) => item.name === inputSelection.dataset);
  const inputProject = inputDataset?.projects.find((item) => item.name === inputSelection.project);

  useEffect(() => { if (!catalog || !section || section.id === sectionId) return; window.history.replaceState(null, "", `?section=${section.id}`); setSectionId(section.id); }, [catalog, section, sectionId]);
  useEffect(() => { setSearch(""); }, [sectionId]);

  const chooseDataset = (name: string) => { const nextDataset = catalog!.datasets.find((item) => item.name === name)!; const nextProject = nextDataset.projects[0]; const nextScenario = nextProject.scenarios[0]; const nextSector = nextScenario.sectors[0]; setSelection({ dataset: name, project: nextProject.name, scenario: nextScenario.name, comparison: "", sector: nextSector.name, year: nextSector.years[0] || "" }); };
  const chooseProject = (name: string) => { const nextProject = dataset!.projects.find((item) => item.name === name)!; const nextScenario = nextProject.scenarios[0]; const nextSector = nextScenario.sectors[0]; setSelection((current) => ({ ...current, project: name, scenario: nextScenario.name, comparison: "", sector: nextSector.name, year: nextSector.years[0] || "" })); };
  const chooseScenario = (name: string) => { const nextScenario = project!.scenarios.find((item) => item.name === name)!; const nextSector = nextScenario.sectors[0]; setSelection((current) => ({ ...current, scenario: name, comparison: current.comparison === name ? "" : current.comparison, sector: nextSector.name, year: nextSector.years[0] || "" })); };
  const chooseSector = (name: string) => { const next = scenario!.sectors.find((item) => item.name === name)!; setSelection((current) => ({ ...current, sector: name, year: next.years[0] || "" })); };
  const chooseSection = (name: string) => { window.history.pushState(null, "", `?section=${name}`); setView("outputs"); setSectionId(name); setSidebarOpen(false); };
  const chooseView = (next: ViewMode) => { window.history.pushState(null, "", next === "outputs" ? `?section=${section?.id || "power"}` : `?view=${next}`); setView(next); setSidebarOpen(false); };
  const chooseInputDataset = (name: string) => { const nextDataset = inputCatalog!.datasets.find((item) => item.name === name)!; const nextProject = nextDataset.projects[0]; setInputSelection({ dataset: name, project: nextProject.name, scenario: nextProject.scenarios[0] || "" }); };
  const chooseInputProject = (name: string) => { const next = inputDataset!.projects.find((item) => item.name === name)!; setInputSelection((current) => ({ ...current, project: name, scenario: next.scenarios[0] || "" })); };
  const refreshCurrent = () => { if (view === "outputs") void loadCatalog(); else { void loadInputs(); setInputRefresh((current) => current + 1); } };

  const outputReady = Boolean(catalog && selection.scenario && section && sector);
  const inputReady = Boolean(inputCatalog && inputSelection.dataset && inputSelection.project && inputSelection.scenario);

  return <div className="app-shell">
    <a className="skip-link" href="#workspace">Skip to workspace</a>
    <aside className={`sidebar ${sidebarOpen ? "open" : ""}`}>
      <div className="brand"><img src="/brand/pypsa-logo.svg" alt="PyPSA-SPICE" /></div>
      <nav><p className="eyebrow">Workspace</p><a className={view === "outputs" ? "active" : ""} href={`?section=${section?.id || "power"}`} onClick={(event) => { event.preventDefault(); chooseView("outputs"); }}><BarChart3 aria-hidden="true" />Outputs</a><a className={view === "inputs" ? "active" : ""} href="?view=inputs" onClick={(event) => { event.preventDefault(); chooseView("inputs"); }}><FileInput aria-hidden="true" />Inputs</a><a className={view === "configure" ? "active" : ""} href="?view=configure" onClick={(event) => { event.preventDefault(); chooseView("configure"); }}><Settings2 aria-hidden="true" />Scenario config</a><span className="disabled"><Play aria-hidden="true" />Run model<small>Next</small></span></nav>
      {view === "outputs" && outputReady ? <div className="source-controls"><p className="eyebrow">Result source</p><Control label="Dataset" value={selection.dataset} onChange={chooseDataset} options={catalog!.datasets.map((item) => item.name)} /><Control label="Project" value={selection.project} onChange={chooseProject} options={dataset!.projects.map((item) => item.name)} /><Control label="Primary scenario" value={selection.scenario} onChange={chooseScenario} options={project!.scenarios.map((item) => item.name)} /><Control label="Compare with" value={selection.comparison} onChange={(comparison) => setSelection((current) => ({ ...current, comparison }))} options={["", ...project!.scenarios.map((item) => item.name).filter((name) => name !== selection.scenario)]} emptyLabel="No comparison" /><Control label="Sector run" value={selection.sector} onChange={chooseSector} options={scenario!.sectors.map((item) => item.name)} /></div> : inputReady ? <div className="source-controls"><p className="eyebrow">Input source</p><Control label="Dataset" value={inputSelection.dataset} onChange={chooseInputDataset} options={inputCatalog!.datasets.map((item) => item.name)} /><Control label="Project" value={inputSelection.project} onChange={chooseInputProject} options={inputDataset!.projects.map((item) => item.name)} /><Control label="Input scenario" value={inputSelection.scenario} onChange={(scenarioName) => setInputSelection((current) => ({ ...current, scenario: scenarioName }))} options={inputProject!.scenarios} /></div> : null}
      <div className="sidebar-foot"><span className="status-dot" />Local data connected<button onClick={() => setDark(!dark)} aria-label="Toggle dark mode">{dark ? <Sun aria-hidden="true" /> : <Moon aria-hidden="true" />}</button></div>
    </aside>
    <div className="scrim" onClick={() => setSidebarOpen(false)} />
    <div className="main-column">
      <header className="topbar"><button className="menu" onClick={() => setSidebarOpen(!sidebarOpen)} aria-label="Open source controls"><Menu aria-hidden="true" /></button>{view === "outputs" && section ? <div className="section-tabs top-section-tabs" role="tablist" aria-label="Output sections">{sections.map((item) => { const SectionIcon = sectionIcons[item.id as keyof typeof sectionIcons] || Zap; return <a className={`section-tab ${item.id === section.id ? "active" : ""}`} href={`?section=${item.id}`} key={item.id} role="tab" aria-selected={item.id === section.id} onClick={(event) => { event.preventDefault(); chooseSection(item.id); }}><SectionIcon aria-hidden="true" /><b>{item.label}</b><small>{item.charts.length}</small></a>; })}</div> : view === "configure" ? <div className="config-topbar-slot" id="config-section-tabs" /> : <div className="topbar-context"><b>Input data editor</b><span>{inputSelection.project}{inputSelection.scenario ? ` · ${inputSelection.scenario}` : ""}</span></div>}<div className="top-actions"><button className="button secondary" onClick={refreshCurrent}><RefreshCw aria-hidden="true" />Refresh data</button><a className="button primary" href="/docs" target="_blank">API</a></div></header>
      <main id="workspace">
        {view === "outputs" ? outputReady ? <><section className="page-title"><h1>{section!.title}</h1><label className="search"><Search aria-hidden="true" /><input aria-label={`Search ${section!.label.toLowerCase()} visualisations`} value={search} onChange={(event) => setSearch(event.target.value)} type="search" placeholder={`Find a ${section!.label.toLowerCase()} visualisation`} /></label></section><section className="analysis">{error && <div className="notice">{error}</div>}<div className={`chart-grid ${selection.comparison ? "comparison-active" : ""}`}>{charts.map((chart) => <ChartCard key={chart.id} chart={chart} selection={selection} years={sector!.years} onYearChange={(year) => setSelection((current) => ({ ...current, year }))} mappings={catalog!.mappings} darkMode={dark} onInspect={(title, rows, sourceCount) => setInspector({ title, rows, sourceCount })} />)}</div>{charts.length === 0 && <div className="no-results">No visualisation matches “{search}”.</div>}</section></> : <Boot message={error || "Discovering local results…"} /> : inputReady ? view === "inputs" ? <InputEditor key={inputRefresh} catalog={inputCatalog!} selection={inputSelection} /> : <ScenarioConfigEditor key={inputRefresh} selection={inputSelection} /> : <Boot message={inputError || "Discovering model inputs…"} />}
      </main>
      <footer><span>PyPSA-SPICE web workspace</span><span>React · FastAPI · Local CSV and YAML source of truth</span></footer>
    </div>
    {inspector && <DataDialog {...inspector} onClose={() => setInspector(null)} />}
  </div>;
}

function Control({ label, value, onChange, options, emptyLabel }: { label: string; value: string; onChange: (value: string) => void; options: string[]; emptyLabel?: string }) { return <label>{label}<select value={value} onChange={(event) => onChange(event.target.value)}>{options.map((option) => <option value={option} key={option}>{option || emptyLabel}</option>)}</select></label>; }
function Boot({ message }: { message: string }) { return <div className="boot inline"><img src="/brand/pypsa-logo.svg" alt="PyPSA" /><span className="spinner" />{message}</div>; }
