import { useEffect, useMemo, useState } from "react";
import { getCatalog } from "./api";
import ChartCard from "./ChartCard";
import DataDialog from "./DataDialog";
import type { Catalog, ResultRow, Selection } from "./types";

const emptySelection: Selection = { dataset: "", project: "", scenario: "", comparison: "", sector: "", year: "" };

export default function App() {
  const [catalog, setCatalog] = useState<Catalog | null>(null);
  const [selection, setSelection] = useState<Selection>(emptySelection);
  const [search, setSearch] = useState("");
  const [sidebarOpen, setSidebarOpen] = useState(false);
  const [error, setError] = useState("");
  const [inspector, setInspector] = useState<{ title: string; rows: ResultRow[]; sourceCount: number } | null>(null);
  const [dark, setDark] = useState(() => localStorage.getItem("spice-theme") === "dark");

  const loadCatalog = async (refresh = false) => {
    try {
      const data = await getCatalog(refresh);
      if (!data.datasets.length) throw new Error("No result CSV folders were found under data/.");
      setCatalog(data); setError("");
      const dataset = data.datasets.find((item) => item.name === selection.dataset) || data.datasets[0];
      const project = dataset.projects.find((item) => item.name === selection.project) || dataset.projects[0];
      const scenario = project.scenarios.find((item) => item.name === selection.scenario) || project.scenarios[0];
      const sector = scenario.sectors.find((item) => item.name === selection.sector) || scenario.sectors[0];
      setSelection((current) => ({ ...current, dataset: dataset.name, project: project.name, scenario: scenario.name, sector: sector.name, year: sector.years.includes(current.year) ? current.year : sector.years[0] || "", comparison: project.scenarios.some((item) => item.name === current.comparison) ? current.comparison : "" }));
    } catch (reason) { setError(reason instanceof Error ? reason.message : "Could not load results."); }
  };

  useEffect(() => { loadCatalog(); }, []);
  useEffect(() => { document.documentElement.dataset.theme = dark ? "dark" : "light"; localStorage.setItem("spice-theme", dark ? "dark" : "light"); }, [dark]);

  const dataset = catalog?.datasets.find((item) => item.name === selection.dataset);
  const project = dataset?.projects.find((item) => item.name === selection.project);
  const scenario = project?.scenarios.find((item) => item.name === selection.scenario);
  const sector = scenario?.sectors.find((item) => item.name === selection.sector);
  const power = catalog?.sections.find((item) => item.id === "power");
  const charts = useMemo(() => power?.charts.filter((chart) => chart.name.toLowerCase().includes(search.toLowerCase().trim())) || [], [power, search]);

  const chooseDataset = (name: string) => {
    const nextDataset = catalog!.datasets.find((item) => item.name === name)!; const nextProject = nextDataset.projects[0]; const nextScenario = nextProject.scenarios[0]; const nextSector = nextScenario.sectors[0];
    setSelection({ dataset: name, project: nextProject.name, scenario: nextScenario.name, comparison: "", sector: nextSector.name, year: nextSector.years[0] || "" });
  };
  const chooseProject = (name: string) => {
    const nextProject = dataset!.projects.find((item) => item.name === name)!; const nextScenario = nextProject.scenarios[0]; const nextSector = nextScenario.sectors[0];
    setSelection((current) => ({ ...current, project: name, scenario: nextScenario.name, comparison: "", sector: nextSector.name, year: nextSector.years[0] || "" }));
  };
  const chooseScenario = (name: string) => {
    const nextScenario = project!.scenarios.find((item) => item.name === name)!; const nextSector = nextScenario.sectors[0];
    setSelection((current) => ({ ...current, scenario: name, comparison: current.comparison === name ? "" : current.comparison, sector: nextSector.name, year: nextSector.years[0] || "" }));
  };
  const chooseSector = (name: string) => { const next = scenario!.sectors.find((item) => item.name === name)!; setSelection((current) => ({ ...current, sector: name, year: next.years[0] || "" })); };

  if (!catalog || !selection.scenario) return <div className="boot"><img src="/brand/pypsa-logo.svg" alt="PyPSA" /><span className="spinner" />{error || "Discovering local results…"}</div>;

  return <div className="app-shell">
    <a className="skip-link" href="#workspace">Skip to visualisations</a>
    <aside className={`sidebar ${sidebarOpen ? "open" : ""}`}>
      <div className="brand"><img src="/brand/pypsa-logo.svg" alt="PyPSA" /><div><b>SPICE</b><span>Model workspace</span></div></div>
      <nav><p className="eyebrow">Workspace</p><a className="active" href="/"><i>⌁</i>Outputs</a><span className="disabled"><i>≡</i>Configure<small>Next</small></span><span className="disabled"><i>▷</i>Run model<small>Next</small></span></nav>
      <div className="source-controls"><p className="eyebrow">Result source</p>
        <Control label="Dataset" value={selection.dataset} onChange={chooseDataset} options={catalog.datasets.map((item) => item.name)} />
        <Control label="Project" value={selection.project} onChange={chooseProject} options={dataset!.projects.map((item) => item.name)} />
        <Control label="Primary scenario" value={selection.scenario} onChange={chooseScenario} options={project!.scenarios.map((item) => item.name)} />
        <Control label="Compare with" value={selection.comparison} onChange={(comparison) => setSelection((current) => ({ ...current, comparison }))} options={["", ...project!.scenarios.map((item) => item.name).filter((name) => name !== selection.scenario)]} emptyLabel="No comparison" />
        <Control label="Sector run" value={selection.sector} onChange={chooseSector} options={scenario!.sectors.map((item) => item.name)} />
        <Control label="Hourly year" value={selection.year} onChange={(year) => setSelection((current) => ({ ...current, year }))} options={sector!.years} />
      </div>
      <div className="sidebar-foot"><span className="status-dot" />Local data connected<button onClick={() => setDark(!dark)} aria-label="Toggle dark mode">{dark ? "☀" : "◐"}</button></div>
    </aside>
    <div className="scrim" onClick={() => setSidebarOpen(false)} />
    <div className="main-column">
      <header className="topbar"><button className="menu" onClick={() => setSidebarOpen(!sidebarOpen)} aria-label="Open result filters">☰</button><div className="crumbs"><span>PyPSA-SPICE</span><i>/</i><b>Power outputs</b></div><div className="top-actions"><button className="button secondary" onClick={() => loadCatalog(true)}>↻ Refresh data</button><a className="button primary" href="/docs" target="_blank">API</a></div></header>
      <main id="workspace">
        <section className="page-title"><h1>Power Sector</h1></section>
        <section className="analysis"><div className="section-tab"><span>ϟ</span><b>Power</b><small>{power?.charts.length}</small></div><div className="analysis-head"><div><p className="eyebrow">Power system</p><h2>Output visualisations</h2></div><div className="analysis-tools"><label className="search"><span>⌕</span><input value={search} onChange={(event) => setSearch(event.target.value)} type="search" placeholder="Find a visualisation" /></label></div></div>{error && <div className="notice">{error}</div>}<div className={`chart-grid ${selection.comparison ? "comparison-active" : ""}`}>{charts.map((chart) => <ChartCard key={chart.id} chart={chart} selection={selection} mappings={catalog.mappings} onInspect={(title, rows, sourceCount) => setInspector({ title, rows, sourceCount })} />)}</div>{charts.length === 0 && <div className="no-results">No visualisation matches “{search}”.</div>}</section>
      </main>
      <footer><span>PyPSA-SPICE power output explorer</span><span>React · FastAPI · Local CSV source of truth</span></footer>
    </div>
    {inspector && <DataDialog {...inspector} onClose={() => setInspector(null)} />}
  </div>;
}

function Control({ label, value, onChange, options, emptyLabel }: { label: string; value: string; onChange: (value: string) => void; options: string[]; emptyLabel?: string }) { return <label>{label}<select value={value} onChange={(event) => onChange(event.target.value)}>{options.map((option) => <option value={option} key={option}>{option || emptyLabel}</option>)}</select></label>; }
