import pageHeaderStyles from "./components/PageHeader.module.scss";
import styles from "./App.module.scss";
import workspaceFeedbackStyles from "./components/WorkspaceFeedback.module.scss";
import scenarioComparisonStyles from "./pages/ScenarioComparison.module.scss";
import workspaceBarStyles from "./components/WorkspaceBar.module.scss";
import dashboardPageStyles from "./pages/DashboardPage.module.scss";
import IconButton from "./components/IconButton";
import { useEffect, useMemo, useState } from "react";
import { ArrowLeftRight, GitCompareArrows, Plus, RefreshCw } from "lucide-react";
import { getCatalog, getInputCatalog } from "./api";
import useModelRun from "./useModelRun";
import { isActiveModelRun } from "./modelRunMonitor";
import ChartCard from "./components/ChartCard";
import DataDialog from "./components/DataDialog";
import { SelectField } from "./components/FormControls";
import NewScenarioDialog from "./components/NewScenarioDialog";
import PageHeader from "./components/PageHeader";
import ResultsToc from "./components/ResultsToc";
import Sidebar from "./components/Sidebar";
import DashboardPage from "./pages/DashboardPage";
import HomePage from "./pages/HomePage";
import InputEditor from "./pages/InputEditor";
import ScenarioConfigEditor from "./pages/ScenarioConfigEditor";
import ScenarioComparison from "./pages/ScenarioComparison";
import {
  confirmDiscardChanges,
  countryFromLocation,
  locationParams,
  resolveInputSelection,
  resolveOutputSelection,
  sectionFromLocation,
  splitWorkspace,
  viewFromLocation,
  workspaceOptions,
  workspaceValue,
} from "./utility";
import type {
  Catalog,
  ChartDefinition,
  InputCatalog,
  InputSelection,
  ModelRunStatus,
  ResultRow,
  Selection,
  ViewMode,
  WorkspaceOption,
} from "./types";

const emptySelection: Selection = { dataset: "", project: "", scenario: "", comparison: "", sector: "", year: "" };
const emptyInputSelection: InputSelection = { dataset: "", project: "", scenario: "" };
const modelRunLabels: Partial<Record<ModelRunStatus, string>> = {
  queued: "Model run queued",
  running: "Model running",
  canceling: "Stopping model run",
};

export default function App() {
  const [catalog, setCatalog] = useState<Catalog | null>(null);
  const [inputCatalog, setInputCatalog] = useState<InputCatalog | null>(null);
  const [selection, setSelection] = useState<Selection>(emptySelection);
  const [inputSelection, setInputSelection] = useState<InputSelection>(emptyInputSelection);
  const [inputComparison, setInputComparison] = useState(() => locationParams().get("compare") || "");
  const [view, setView] = useState<ViewMode>(viewFromLocation);
  const [sectionId, setSectionId] = useState(sectionFromLocation);
  const [country, setCountry] = useState(countryFromLocation);
  const [error, setError] = useState("");
  const [inputError, setInputError] = useState("");
  const [inputRefresh, setInputRefresh] = useState(0);
  const [inspector, setInspector] = useState<{ chart: ChartDefinition; rows: ResultRow[]; sourceCount: number } | null>(
    null,
  );
  const [dark, setDark] = useState(() => localStorage.getItem("spice-theme") === "dark");
  const [newScenarioOpen, setNewScenarioOpen] = useState(false);
  const runMonitor = useModelRun();
  const modelRunStatus = isActiveModelRun(runMonitor.run) ? runMonitor.run?.status : null;

  const loadCatalog = async () => {
    try {
      const data = await getCatalog(true);
      if (!data.datasets.length) throw new Error("No result CSV folders were found under data/.");
      setCatalog(data);
      setError("");
      setSelection((current) => resolveOutputSelection(data, current));
    } catch (reason) {
      setError(reason instanceof Error ? reason.message : "Could not load results.");
    }
  };
  const loadInputs = async (preferred?: InputSelection) => {
    try {
      const data = await getInputCatalog();
      if (!data.datasets.length) throw new Error("No input projects were found under data/.");
      setInputCatalog(data);
      setInputError("");
      setInputSelection((current) =>
        preferred
          ? resolveInputSelection(data, preferred, new URLSearchParams())
          : resolveInputSelection(data, current),
      );
    } catch (reason) {
      setInputError(reason instanceof Error ? reason.message : "Could not load model inputs.");
    }
  };

  useEffect(() => {
    void loadCatalog();
    void loadInputs();
  }, []);
  useEffect(() => {
    document.documentElement.dataset.theme = dark ? "dark" : "light";
    localStorage.setItem("spice-theme", dark ? "dark" : "light");
  }, [dark]);
  useEffect(() => {
    const syncLocation = () => {
      const params = locationParams();
      setView(viewFromLocation());
      setSectionId(sectionFromLocation());
      setCountry(countryFromLocation());
      if (catalog) setSelection((current) => resolveOutputSelection(catalog, current, params));
      if (inputCatalog) {
        setInputSelection((current) => resolveInputSelection(inputCatalog, current, params));
        setInputComparison(params.get("compare") || "");
      }
    };
    window.addEventListener("popstate", syncLocation);
    return () => window.removeEventListener("popstate", syncLocation);
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
  const projectChoices = useMemo(
    () => workspaceOptions(view === "outputs" ? catalog?.datasets || [] : inputCatalog?.datasets || []),
    [view, catalog, inputCatalog],
  );

  useEffect(() => {
    if (!catalog || !section || section.id === sectionId) return;
    setSectionId(section.id);
  }, [catalog, section, sectionId]);
  useEffect(() => {
    if (country !== "ALL" && inputProject && !configurationCountries.includes(country)) setCountry("ALL");
  }, [country, inputProject, configurationCountries]);
  useEffect(() => {
    if (!inputProject) return;
    const alternatives = inputProject.scenarios.filter((name) => name !== inputSelection.scenario);
    if (!alternatives.includes(inputComparison)) setInputComparison(alternatives[0] || "");
  }, [inputProject, inputSelection.scenario, inputComparison]);
  useEffect(() => {
    if (!(selection.project || inputSelection.project)) return;
    const params = new URLSearchParams();
    if (view === "home") {
      params.set("view", "home");
    } else if (view === "dashboard") {
      params.set("view", "dashboard");
    } else {
      if (view !== "outputs") params.set("view", view);
      else params.set("section", section?.id || sectionId);
      if (activeDatasetName) params.set("dataset", activeDatasetName);
      if (activeProjectName) params.set("project", activeProjectName);
      if (inputSelection.scenario) params.set("scenario", inputSelection.scenario);
      if (selection.scenario) params.set("run", selection.scenario);
      if (view === "outputs" && selection.comparison) params.set("compare", selection.comparison);
      if (view === "compare" && inputComparison) params.set("compare", inputComparison);
      if (selection.sector) params.set("sector", selection.sector);
      if (view === "configure" && country !== "ALL") params.set("country", country);
    }
    const next = `${window.location.pathname}?${params.toString()}`;
    if (`${window.location.pathname}${window.location.search}` !== next) window.history.replaceState(null, "", next);
  }, [
    view,
    section?.id,
    sectionId,
    activeDatasetName,
    activeProjectName,
    inputSelection.scenario,
    inputComparison,
    selection.scenario,
    selection.comparison,
    selection.sector,
    country,
  ]);

  const setOutputWorkspace = (datasetName: string, projectName: string) => {
    if (!catalog) return;
    const nextDataset = catalog.datasets.find((item) => item.name === datasetName);
    const nextProject = nextDataset?.projects.find((item) => item.name === projectName);
    if (!nextDataset || !nextProject) return;
    const nextScenario = nextProject.scenarios[0];
    const nextSector = nextScenario.sectors[0];
    setSelection({
      dataset: datasetName,
      project: projectName,
      scenario: nextScenario.name,
      comparison: "",
      sector: nextSector.name,
      year: nextSector.years[0] || "",
    });
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
  const chooseScenario = (name: string) => {
    if (!confirmDiscardChanges()) return;
    const nextScenario = project!.scenarios.find((item) => item.name === name)!;
    const nextSector = nextScenario.sectors[0];
    setSelection((current) => ({
      ...current,
      scenario: name,
      comparison: current.comparison === name ? "" : current.comparison,
      sector: nextSector.name,
      year: nextSector.years[0] || "",
    }));
  };
  const chooseInputScenario = (name: string) => {
    if (!confirmDiscardChanges()) return;
    setInputSelection((current) => ({ ...current, scenario: name }));
    if (name === inputComparison) {
      setInputComparison(inputProject?.scenarios.find((scenarioName) => scenarioName !== name) || "");
    }
  };
  const chooseCountry = (name: string) => {
    if (confirmDiscardChanges()) setCountry(name);
  };
  const chooseSector = (name: string) => {
    if (!confirmDiscardChanges()) return;
    const next = scenario!.sectors.find((item) => item.name === name)!;
    setSelection((current) => ({ ...current, sector: name, year: next.years[0] || "" }));
  };
  const chooseSection = (name: string) => {
    if (!confirmDiscardChanges()) return;
    setView("outputs");
    setSectionId(name);
    window.scrollTo({ top: 0, behavior: "smooth" });
  };
  const chooseView = (next: ViewMode) => {
    if (!confirmDiscardChanges()) return;
    const params = locationParams();
    if (next === "outputs") {
      params.delete("view");
      params.set("section", section?.id || "power");
    } else {
      params.set("view", next);
      params.delete("section");
    }
    if (next !== "configure") params.delete("country");
    if (next !== "outputs" && next !== "compare") params.delete("compare");
    window.history.pushState(null, "", `?${params.toString()}`);
    setView(next);
    window.scrollTo({ top: 0, behavior: "smooth" });
  };
  const refreshCurrent = () => {
    if (!confirmDiscardChanges()) return;
    if (view === "home") {
      void loadCatalog();
      void loadInputs();
    } else if (view === "outputs" || view === "dashboard") void loadCatalog();
    else {
      void loadInputs();
      setInputRefresh((current) => current + 1);
    }
  };

  const outputReady = Boolean(catalog && selection.scenario && section && sector);
  const inputReady = Boolean(
    inputCatalog && inputSelection.dataset && inputSelection.project && inputSelection.scenario,
  );
  const comparisonReady = Boolean(inputReady && inputComparison);
  const contextReady = view === "home" || view === "dashboard" ? false : view === "outputs" ? outputReady : inputReady;
  const activeWorkspace = workspaceValue(activeDatasetName, activeProjectName);
  const modelRunLabel = runMonitor.error
    ? "Run status unavailable; retrying"
    : modelRunStatus
      ? modelRunLabels[modelRunStatus]
      : undefined;

  const renderWorkspace = () => {
    if (view === "home") {
      return (
        <HomePage
          catalog={catalog}
          inputCatalog={inputCatalog}
          resultError={error}
          inputError={inputError}
          onNavigate={chooseView}
        />
      );
    }
    if (view === "dashboard") {
      return catalog ? (
        <DashboardPage
          catalog={catalog}
          darkMode={dark}
          onInspect={(inspectedChart, rows, sourceCount) => setInspector({ chart: inspectedChart, rows, sourceCount })}
        />
      ) : (
        <Boot message={error || "Discovering local results…"} />
      );
    }
    if (view === "outputs") {
      if (!outputReady) return <Boot message={error || "Discovering local results…"} />;
      return (
        <>
          <PageHeader title={section!.title}>
            {scenario!.sectors.length > 1 && (
              <div className={pageHeaderStyles["page-controls"]}>
                <ContextControl
                  label="Sector run"
                  value={selection.sector}
                  onChange={chooseSector}
                  options={scenario!.sectors.map((item) => ({ value: item.name, label: item.name }))}
                />
              </div>
            )}
          </PageHeader>
          <section className={styles["analysis"]}>
            {error && <div className={workspaceFeedbackStyles["notice"]}>{error}</div>}
            <div
              className={[styles["chart-grid"], selection.comparison ? styles["comparison-active"] : ""]
                .filter(Boolean)
                .join(" ")}
            >
              {charts.map((chart) => (
                <ChartCard
                  key={chart.id}
                  chart={chart}
                  selection={selection}
                  years={sector!.years}
                  mappings={catalog!.mappings}
                  darkMode={dark}
                  onInspect={(inspectedChart, rows, sourceCount) =>
                    setInspector({ chart: inspectedChart, rows, sourceCount })
                  }
                />
              ))}
            </div>
            {charts.length === 0 && (
              <div className={workspaceFeedbackStyles["no-results"]}>
                No visualisations are configured for this section.
              </div>
            )}
          </section>
          <ResultsToc
            key={section!.id}
            id="results-figure-list"
            heading="Figures"
            panelLabel="Figures on this page"
            listName="figure list"
            entries={charts.map((chart) => ({ key: chart.id, href: `#figure-${chart.id}`, label: chart.name }))}
          />
        </>
      );
    }
    if (view === "compare") {
      if (comparisonReady) {
        return (
          <ScenarioComparison
            key={`${inputRefresh}:${inputSelection.scenario}:${inputComparison}`}
            selection={inputSelection}
            comparison={inputComparison}
          />
        );
      }
      if (inputReady) {
        return (
          <div className={scenarioComparisonStyles["comparison-empty"]}>
            <GitCompareArrows aria-hidden="true" />
            <b>Two scenarios are required</b>
            <span>Create another scenario before opening the comparison workspace.</span>
          </div>
        );
      }
      return <Boot message={inputError || "Discovering model inputs…"} />;
    }
    if (!inputReady) return <Boot message={inputError || "Discovering model inputs…"} />;
    if (view === "inputs") {
      return <InputEditor key={inputRefresh} catalog={inputCatalog!} selection={inputSelection} />;
    }
    return (
      <ScenarioConfigEditor
        key={inputRefresh}
        selection={inputSelection}
        country={country}
        runMonitor={runMonitor}
        onOpenResults={(runName, datasetName, projectName) => {
          window.location.href = `/?section=power&dataset=${encodeURIComponent(datasetName)}&project=${encodeURIComponent(projectName)}&run=${encodeURIComponent(runName)}`;
        }}
      />
    );
  };

  return (
    <div className={styles["app-shell"]}>
      <a className={styles["skip-link"]} href="#workspace">
        Skip to workspace
      </a>
      <Sidebar
        view={view}
        sections={sections}
        activeSectionId={section?.id}
        darkMode={dark}
        modelRunLabel={modelRunLabel}
        onSelectView={chooseView}
        onSelectSection={chooseSection}
        onToggleDarkMode={() => setDark((current) => !current)}
      />
      <div className={styles["main-column"]}>
        <header className={workspaceBarStyles["workspace-bar"]}>
          {view === "home" && (
            <div className={workspaceBarStyles["home-bar-copy"]}>
              <b>Workspace overview</b>
              <span>Open a local project, select a scenario, or continue with a saved dashboard.</span>
            </div>
          )}
          {contextReady && (
            <div className={workspaceBarStyles["workspace-context"]}>
              <ContextControl
                className={workspaceBarStyles["project-control"]}
                label="Project"
                value={activeWorkspace}
                onChange={chooseWorkspace}
                options={projectChoices}
              />
              {view === "outputs" ? (
                <>
                  <ContextControl
                    className={workspaceBarStyles["scenario-control"]}
                    label="Result run"
                    value={selection.scenario}
                    onChange={chooseScenario}
                    options={project!.scenarios.map((item) => ({ value: item.name, label: item.name }))}
                  />
                  <ContextControl
                    className={[
                      workspaceBarStyles["compare-control"],
                      selection.comparison ? workspaceBarStyles["is-active"] : "",
                    ]
                      .filter(Boolean)
                      .join(" ")}
                    label="Compare with"
                    value={selection.comparison}
                    onChange={(comparison) => setSelection((current) => ({ ...current, comparison }))}
                    options={[
                      { value: "", label: "No comparison" },
                      ...project!.scenarios
                        .filter((item) => item.name !== selection.scenario)
                        .map((item) => ({ value: item.name, label: item.name })),
                    ]}
                  />
                </>
              ) : view === "compare" ? (
                <>
                  <ContextControl
                    className={workspaceBarStyles["scenario-control"]}
                    label="Reference scenario"
                    value={inputSelection.scenario}
                    onChange={chooseInputScenario}
                    options={inputProject!.scenarios.map((item) => ({ value: item, label: item }))}
                  />
                  <IconButton
                    variant="surface"
                    className={workspaceBarStyles["comparison-swap"]}
                    onClick={() => {
                      const previous = inputSelection.scenario;
                      setInputSelection((current) => ({ ...current, scenario: inputComparison }));
                      setInputComparison(previous);
                    }}
                    aria-label="Swap reference and comparison scenarios"
                    title="Swap scenarios"
                    disabled={!inputComparison}
                  >
                    <ArrowLeftRight aria-hidden="true" />
                  </IconButton>
                  <ContextControl
                    className={[workspaceBarStyles["compare-control"], workspaceBarStyles["is-active"]].join(" ")}
                    label="Comparison scenario"
                    value={inputComparison}
                    onChange={setInputComparison}
                    options={inputProject!.scenarios
                      .filter((item) => item !== inputSelection.scenario)
                      .map((item) => ({ value: item, label: item }))}
                  />
                </>
              ) : (
                <>
                  <div className={workspaceBarStyles["scenario-context-group"]}>
                    <ContextControl
                      className={workspaceBarStyles["scenario-control"]}
                      label="Scenario"
                      value={inputSelection.scenario}
                      onChange={chooseInputScenario}
                      options={inputProject!.scenarios.map((item) => ({ value: item, label: item }))}
                    />
                    {view === "configure" && (
                      <IconButton
                        variant="surface"
                        className={workspaceBarStyles["context-add"]}
                        onClick={() => setNewScenarioOpen(true)}
                        aria-label="Create new scenario"
                        title="Create new scenario"
                      >
                        <Plus aria-hidden="true" />
                      </IconButton>
                    )}
                  </div>
                  {view === "inputs" && (
                    <div className={workspaceBarStyles["input-topbar-controls"]} id="input-topbar-controls" />
                  )}
                  {view === "configure" && (
                    <ContextControl
                      className={workspaceBarStyles["country-control"]}
                      label="Country"
                      value={country}
                      onChange={chooseCountry}
                      options={[
                        { value: "ALL", label: "All countries" },
                        ...configurationCountries.map((item) => ({ value: item, label: item })),
                      ]}
                    />
                  )}
                </>
              )}
            </div>
          )}
          {view === "dashboard" && (
            <div className={workspaceBarStyles["workspace-context"]}>
              <div id="dashboard-topbar-controls" className={dashboardPageStyles["dashboard-topbar-controls"]} />
            </div>
          )}
          {view === "home" && modelRunLabel && (
            <span className={workspaceBarStyles["home-run-status"]} role="status" aria-live="polite">
              <i />
              {modelRunLabel}
            </span>
          )}
          <div className={workspaceBarStyles["top-actions"]}>
            <IconButton variant="surface" onClick={refreshCurrent} aria-label="Refresh data" title="Refresh data">
              <RefreshCw aria-hidden="true" />
            </IconButton>
          </div>
        </header>
        <main className={styles["main-content"]} id="workspace">
          {renderWorkspace()}
        </main>
      </div>
      {inspector && <DataDialog {...inspector} onClose={() => setInspector(null)} />}
      {newScenarioOpen && inputProject && (
        <NewScenarioDialog
          selection={inputSelection}
          scenarios={inputProject.scenarios}
          onClose={() => setNewScenarioOpen(false)}
          onCreated={(created) => {
            setNewScenarioOpen(false);
            void loadInputs({ dataset: created.dataset, project: created.project, scenario: created.scenario }).then(
              () => setInputRefresh((current) => current + 1),
            );
          }}
        />
      )}
    </div>
  );
}

function ContextControl({
  className = "",
  label,
  value,
  onChange,
  options,
}: {
  className?: string;
  label: string;
  value: string;
  onChange: (value: string) => void;
  options: WorkspaceOption[];
}) {
  return (
    <SelectField
      variant="context"
      className={className}
      label={label}
      value={value}
      onChange={onChange}
      options={options}
    />
  );
}
function Boot({ message }: { message: string }) {
  return (
    <div className={[styles["boot"], styles["inline"]].join(" ")}>
      <img src="/ui/pypsa-logo.svg" alt="PyPSA" />
      <span className={workspaceFeedbackStyles["spinner"]} />
      {message}
    </div>
  );
}
