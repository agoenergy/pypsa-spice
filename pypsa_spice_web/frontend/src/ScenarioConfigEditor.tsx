import { useEffect, useId, useMemo, useRef, useState } from "react";
import { createPortal } from "react-dom";
import { Check, ClipboardCheck, Code2, List, Plus, RotateCcw, Save, Settings2, Trash2, X } from "lucide-react";
import { getFuelSupplyTable, getScenarioConfig, saveFuelSupplyLimits, saveScenarioConfigSection, saveScenarioConfigSections } from "./api";
import type { InputRow, InputSelection, InputTableResponse, ScenarioConfigResponse } from "./types";
import { confirmDiscardChanges, setEditorDirty } from "./dirtyState";
import RunModel from "./RunModel";

const COMBINED_SECTION = "co2_constraints";
const labels: Record<string, string> = { scenario_configs: "Scenario settings", [COMBINED_SECTION]: "CO₂ & constraints", review_run: "Review & run" };
const icons = { scenario_configs: Settings2, [COMBINED_SECTION]: Code2, review_run: ClipboardCheck };

function normaliseSection(requested: string) {
  if (requested === "co2_management" || requested === "custom_constraints") return COMBINED_SECTION;
  return Object.hasOwn(labels, requested) ? requested : "scenario_configs";
}

function draftForSection(config: ScenarioConfigResponse, section: string): Record<string, unknown> {
  if (section === COMBINED_SECTION) {
    return {
      co2_management: structuredClone(config.sections.co2_management || {}),
      custom_constraints: structuredClone(config.sections.custom_constraints || {}),
    };
  }
  return structuredClone(config.sections[section] || {});
}

function fuelConstraintSnapshot(value: Record<string, unknown>): Record<string, unknown> {
  return Object.fromEntries(Object.entries(objectValue(value.custom_constraints)).map(([country, raw]) => [country, objectValue(objectValue(raw).production_constraint_fuels)]));
}

export default function ScenarioConfigEditor({ selection, country, initialSection, onNavigate, onOpenResults }: { selection: InputSelection; country: string; initialSection?: string; onNavigate: () => void; onOpenResults: (runName: string, dataset: string, project: string) => void }) {
  const editorId = useId();
  const [config, setConfig] = useState<ScenarioConfigResponse | null>(null);
  const [fuelTable, setFuelTable] = useState<InputTableResponse | null>(null);
  const [fuelRows, setFuelRows] = useState<InputRow[]>([]);
  const [fuelError, setFuelError] = useState("");
  const [section, setSection] = useState(() => {
    const requested = new URLSearchParams(window.location.search).get("step") || initialSection || "scenario_configs";
    return normaliseSection(requested);
  });
  const [draft, setDraft] = useState<Record<string, unknown>>({});
  const [loading, setLoading] = useState(true);
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState("");
  const [success, setSuccess] = useState("");
  const [navigationTarget, setNavigationTarget] = useState<HTMLElement | null>(null);

  const load = async (signal?: AbortSignal) => {
    setLoading(true); setError(""); setFuelError(""); setSuccess("");
    try {
      const [data, fuelResult] = await Promise.all([
        getScenarioConfig(selection, signal),
        getFuelSupplyTable(selection, signal).then((value) => ({ value, error: "" })).catch((reason) => ({ value: null, error: reason instanceof Error ? reason.message : "Could not read fuel supply limits." })),
      ]);
      setConfig(data); setDraft(draftForSection(data, section));
      setFuelTable(fuelResult.value); setFuelRows(structuredClone(fuelResult.value?.rows || [])); setFuelError(fuelResult.error);
    }
    catch (reason) { if (!(reason instanceof DOMException && reason.name === "AbortError")) setError(reason instanceof Error ? reason.message : "Could not load the scenario config."); }
    finally { if (!signal?.aborted) setLoading(false); }
  };
  useEffect(() => {
    const controller = new AbortController();
    void load(controller.signal);
    return () => controller.abort();
  }, [selection.dataset, selection.project, selection.scenario]);
  useEffect(() => { if (config) { setDraft(draftForSection(config, section)); setFuelRows(structuredClone(fuelTable?.rows || [])); setSuccess(""); setError(""); } }, [section]);
  useEffect(() => { setNavigationTarget(document.getElementById("config-section-tabs")); }, []);

  const original = useMemo(() => config ? draftForSection(config, section) : {}, [config, section]);
  const configDirty = useMemo(() => JSON.stringify(draft) !== JSON.stringify(original), [draft, original]);
  const fuelChanges = useMemo(() => {
    if (section !== COMBINED_SECTION || !fuelTable) return [];
    const originals = new Map(fuelTable.rows.map((row) => [row.__row_id, row.max_supply__mwh_year]));
    return fuelRows.flatMap((row) => JSON.stringify(row.max_supply__mwh_year) === JSON.stringify(originals.get(row.__row_id)) ? [] : [{ row: row.__row_id, column: "max_supply__mwh_year", value: row.max_supply__mwh_year }]);
  }, [fuelRows, fuelTable, section]);
  const fuelDirty = fuelChanges.length > 0;
  const dirty = configDirty || fuelDirty;
  const fuelConstraintDirty = section === COMBINED_SECTION && JSON.stringify(fuelConstraintSnapshot(draft)) !== JSON.stringify(fuelConstraintSnapshot(original));
  const validationError = useMemo(() => validateSection(section, draft) || (fuelError ? fuelConstraintDirty ? fuelError : "" : validateFuelLimits(section, draft, fuelRows)), [section, draft, fuelRows, fuelError, fuelConstraintDirty]);
  const invalid = Boolean(validationError);
  useEffect(() => { const warn = (event: BeforeUnloadEvent) => { if (dirty) event.preventDefault(); }; window.addEventListener("beforeunload", warn); return () => window.removeEventListener("beforeunload", warn); }, [dirty]);
  useEffect(() => { setEditorDirty(editorId, dirty); return () => setEditorDirty(editorId, false); }, [editorId, dirty]);
  const save = async () => {
    if (!config || !dirty) return; setSaving(true); setError(""); setSuccess("");
    try {
      let savedFuel: InputTableResponse | null = null;
      if (fuelDirty && fuelTable) {
        savedFuel = await saveFuelSupplyLimits(selection, fuelTable.revision, fuelChanges);
        setFuelTable({ ...fuelTable, revision: savedFuel.revision, rows: structuredClone(fuelRows) });
      }
      const changedSections = section === COMBINED_SECTION ? Object.fromEntries(
        ["co2_management", "custom_constraints"]
          .filter((name) => JSON.stringify(draft[name]) !== JSON.stringify(original[name]))
          .map((name) => [name, draft[name] as Record<string, unknown>]),
      ) : {};
      const data = configDirty
        ? section === COMBINED_SECTION
          ? await saveScenarioConfigSections(selection, config.revision, changedSections)
          : await saveScenarioConfigSection(selection, section, config.revision, draft)
        : config;
      setConfig(data); setDraft(draftForSection(data, section));
      setSuccess("Changes saved.");
    }
    catch (reason) { setError(reason instanceof Error ? reason.message : "Could not save this section."); }
    finally { setSaving(false); }
  };

  const chooseSection = (name: string) => {
    if (!confirmDiscardChanges()) return;
    setSection(name); onNavigate();
    const params = new URLSearchParams(window.location.search);
    params.set("view", "configure"); params.set("step", name); params.delete("section");
    window.history.pushState(null, "", `?${params.toString()}`);
    window.scrollTo({ top: 0, behavior: "smooth" });
  };

  const navigation = navigationTarget && createPortal(<nav className="sidebar-submenu-list" aria-label="Configuration and run pages">{Object.keys(labels).map((name) => { const SectionIcon = icons[name as keyof typeof icons]; return <button key={name} className={`sidebar-submenu-item ${section === name ? "active" : ""}`} onClick={() => chooseSection(name)} aria-current={section === name ? "page" : undefined}><SectionIcon aria-hidden="true" /><b>{labels[name]}</b></button>; })}</nav>, navigationTarget);

  if (section === "review_run") return <>{navigation}<RunModel selection={selection} onEditConfiguration={() => chooseSection("scenario_configs")} onOpenResults={onOpenResults} /></>;

  return <>
    {navigation}
    <section className={`page-title editor-title ${section === COMBINED_SECTION ? "compact-config-title" : ""}`}><div><p className="eyebrow pink">Configure &amp; run</p><h1>Configure {selection.scenario}</h1><p>Edit the model settings and constraints stored in this scenario’s input files.</p></div></section>
    <div className="config-layout">
      <section className={`editor-panel config-panel ${section === COMBINED_SECTION ? "combined-config-panel" : ""}`}>
        <header className="editor-panel-head"><div><p className="eyebrow">{selection.project} · {selection.scenario}</p><h2>{labels[section]}</h2>{config && <code>{config.path}</code>}</div></header>
        {error && <div className="notice error">{error}<button onClick={() => void load()}>Reload</button></div>}
        {validationError && <div className="notice error">{validationError}</div>}
        {loading ? <div className="editor-loading"><span className="spinner" />Reading configuration…</div> : section === "scenario_configs" ? <ScenarioSettings value={draft} country={country} onChange={setDraft} /> : <CombinedConstraintsEditor value={draft} country={country} fuelRows={fuelRows} fuelError={fuelError} onFuelRowsChange={setFuelRows} onChange={setDraft} />}
        <footer className="config-actions floating-config-actions"><button className="button secondary" disabled={!dirty || saving} onClick={() => { setDraft(structuredClone(original)); setFuelRows(structuredClone(fuelTable?.rows || [])); }}><RotateCcw aria-hidden="true" />Discard</button><button className="button primary" disabled={!dirty || saving || invalid} onClick={save}><Save aria-hidden="true" />{saving ? "Saving…" : "Save changes"}</button>{success && <span className="save-message" role="status"><Check aria-hidden="true" />{success}</span>}</footer>
      </section>
    </div>
  </>;
}

function CombinedConstraintsEditor({ value, country, fuelRows, fuelError, onFuelRowsChange, onChange }: { value: Record<string, unknown>; country: string; fuelRows: InputRow[]; fuelError: string; onFuelRowsChange: (rows: InputRow[]) => void; onChange: (value: Record<string, unknown>) => void }) {
  const co2 = (value.co2_management || {}) as Record<string, unknown>;
  const constraints = (value.custom_constraints || {}) as Record<string, unknown>;
  const visibleConstraints = country === "ALL" ? Object.values(constraints) : [constraints[country]];
  const constraintNames = [...new Set([...Object.keys(CONSTRAINT_DEFAULTS), ...visibleConstraints.flatMap((raw) => Object.keys(objectValue(raw)))])];
  const tocItems = [{ id: "config-co2-management", label: "CO₂ management" }, ...constraintNames.map((name) => ({ id: constraintAnchor(name), label: prettyConfigLabel(name) }))];
  return <><div className="combined-config">
    <section className="combined-config-section" id="config-co2-management" aria-labelledby="co2-management-heading"><header><p className="eyebrow pink">Carbon policy</p><h3 id="co2-management-heading">CO₂ management</h3><p>Country carbon cap or price by year.</p></header><Co2Editor value={co2} country={country} onChange={(next) => onChange({ ...value, co2_management: next })} /></section>
    <section className="combined-config-section" aria-labelledby="custom-constraints-heading"><header><p className="eyebrow pink">Model rules</p><h3 id="custom-constraints-heading">Custom constraints</h3><p>Activate constraints and edit their parameters directly.</p></header><CustomConstraintsEditor value={constraints} country={country} countries={Object.keys(co2)} fuelRows={fuelRows} fuelError={fuelError} onFuelRowsChange={onFuelRowsChange} onChange={(next) => onChange({ ...value, custom_constraints: next })} /></section>
  </div><ConfigToc items={tocItems} /></>;
}

function ScenarioSettings({ value, country, onChange }: { value: Record<string, unknown>; country: string; onChange: (value: Record<string, unknown>) => void }) {
  const snapshots = (value.snapshots || {}) as Record<string, unknown>;
  const resolution = (value.resolution || {}) as Record<string, unknown>;
  const interest = (value.interest || {}) as Record<string, unknown>;
  const start = String(snapshots.start || ""); const end = String(snapshots.end || "");
  const modelYear = Number(start.slice(0, 4)) || new Date().getFullYear();
  const [manual, setManual] = useState(!start.endsWith("-01-01") || end !== `${modelYear + 1}-01-01` || snapshots.inclusive !== "left");
  const patch = (key: string, child: Record<string, unknown>) => onChange({ ...value, [key]: child });
  return <div className="config-form">
    <div className="form-grid"><label className="field"><span>Model year</span><input type="number" min="2010" max="3000" value={modelYear} onChange={(event) => { const year = Number(event.target.value); patch("snapshots", { ...snapshots, start: `${year}-01-01`, end: `${year + 1}-01-01`, inclusive: "left" }); }} /></label><label className="field"><span>Remove assets below (MW)</span><input type="number" min="0" step="0.1" value={String(value.remove_threshold ?? 0)} onChange={(event) => onChange({ ...value, remove_threshold: Number(event.target.value) })} /></label></div>
    <label className="toggle-row"><input type="checkbox" checked={manual} onChange={(event) => { setManual(event.target.checked); if (!event.target.checked) patch("snapshots", { start: `${modelYear}-01-01`, end: `${modelYear + 1}-01-01`, inclusive: "left" }); }} /><span>Edit snapshot range manually</span></label>
    {manual ? <div className="form-grid three"><label className="field"><span>Snapshot start</span><input type="date" value={start} onChange={(event) => patch("snapshots", { ...snapshots, start: event.target.value })} /></label><label className="field"><span>Snapshot end</span><input type="date" value={end} onChange={(event) => patch("snapshots", { ...snapshots, end: event.target.value })} /></label><label className="field"><span>Inclusive</span><select value={String(snapshots.inclusive || "left")} onChange={(event) => patch("snapshots", { ...snapshots, inclusive: event.target.value })}>{["both", "neither", "left", "right"].map((item) => <option key={item}>{item}</option>)}</select></label></div> : <p className="field-help">Full-year hourly range: {modelYear}-01-01 to {modelYear + 1}-01-01, inclusive left.</p>}
    <div className="form-section"><h3>Temporal resolution</h3><div className="form-grid"><label className="field"><span>Method</span><select value={String(resolution.method || "nth_hour")} onChange={(event) => patch("resolution", { ...resolution, method: event.target.value })}><option value="nth_hour">Every nth hour</option><option value="clustered">Clustered representative days</option></select></label>{resolution.method === "clustered" ? <label className="field"><span>Number of days</span><input type="number" min="1" value={String(resolution.number_of_days ?? 3)} onChange={(event) => patch("resolution", { ...resolution, number_of_days: Number(event.target.value) })} /></label> : <label className="field"><span>Step size</span><input type="number" min="1" value={String(resolution.stepsize ?? 25)} onChange={(event) => patch("resolution", { ...resolution, stepsize: Number(event.target.value) })} /></label>}</div></div>
    <div className="form-section"><h3>Interest rates</h3><p className="field-help">Country-specific decimal rates; 0.05 means 5%.</p><CountryValueEditor value={interest} country={country} onChange={(next) => onChange({ ...value, interest: next })} /></div>
  </div>;
}

function Co2Editor({ value, country, onChange }: { value: Record<string, unknown>; country: string; onChange: (value: Record<string, unknown>) => void }) {
  const updateCountry = (country: string, next: Record<string, unknown>) => onChange({ ...value, [country]: next });
  const entries = Object.entries(value).filter(([name]) => country === "ALL" || name === country);
  return <div className="config-form">{entries.map(([name, raw]) => { const config = (raw || {}) as Record<string, unknown>; return <article className="country-config" key={name}><header><h3>{name}</h3><label className="field compact"><span>Instrument</span><select value={String(config.option || "co2_cap")} onChange={(event) => updateCountry(name, { ...config, option: event.target.value })}><option value="co2_cap">CO₂ cap</option><option value="co2_price">CO₂ price</option></select></label></header><MappingTable label="Year" value={(config.value || {}) as Record<string, unknown>} yearKeys allowAdd={false} onChange={(next) => updateCountry(name, { ...config, value: next })} /></article>; })}{entries.length === 0 && <div className="editor-empty compact">No CO₂ configuration is defined for {country}.</div>}</div>;
}

function CountryValueEditor({ value, country, onChange }: { value: Record<string, unknown>; country: string; onChange: (value: Record<string, unknown>) => void }) {
  if (country === "ALL") return <KeyValueEditor value={value} onChange={onChange} />;
  return <div className="key-value-grid"><label className="field"><span>{country}</span><input type="number" step="any" value={String(value[country] ?? "")} onChange={(event) => onChange({ ...value, [country]: Number(event.target.value) })} /></label></div>;
}

function KeyValueEditor({ value, onChange }: { value: Record<string, unknown>; onChange: (value: Record<string, unknown>) => void }) {
  const remove = (key: string) => { const next = { ...value }; delete next[key]; onChange(next); };
  return <div className="mapping-editor"><div className="key-value-grid">{Object.entries(value).map(([key, raw]) => <div className="mapping-field" key={key}><label className="field"><span>{key}</span><input type="number" step="any" value={String(raw ?? "")} onChange={(event) => onChange({ ...value, [key]: Number(event.target.value) })} /></label><button className="icon-button" aria-label={`Remove ${key}`} onClick={() => remove(key)}><Trash2 aria-hidden="true" /></button></div>)}</div></div>;
}

const YEAR_FIELDS = new Set(["ei_fraction", "res_generation_share"]);
const CONFIG_LABELS: Record<string, string> = {
  energy_independence: "Energy independence",
  production_constraint_fuels: "Fuel production limits",
  reserve_margin: "Reserve margin",
  res_generation: "Renewable generation share",
  thermal_must_run: "Thermal must-run",
  capacity_factor_constraint: "Capacity factor",
  maximum_power_generation_constraint: "Maximum power generation",
  pe_conv_fraction: "Primary energy conversion fractions",
  ei_fraction: "Energy independence fraction",
  fuels: "Restricted fuels",
  epsilon_load: "Load reserve fraction",
  epsilon_vre: "VRE reserve contribution",
  contingency: "Contingency (MW)",
  res_generation_share: "Renewable generation share",
  min_must_run_ratio: "Minimum must-run ratio",
};
const CONSTRAINT_DEFAULTS: Record<string, Record<string, unknown>> = {
  energy_independence: { activate: false, pe_conv_fraction: { Solar: null, Wind: null, Geothermal: null, Water: null }, ei_fraction: {} },
  production_constraint_fuels: { activate: false, fuels: [] },
  reserve_margin: { activate: false, epsilon_load: null, epsilon_vre: null, contingency: null, method: "static" },
  res_generation: { activate: false, math_symbol: "<=", res_generation_share: {} },
  thermal_must_run: { activate: false, min_must_run_ratio: null },
  capacity_factor_constraint: { activate: false, value: {} },
  maximum_power_generation_constraint: { activate: false, value: {} },
};

function objectValue(value: unknown): Record<string, unknown> {
  return value && typeof value === "object" && !Array.isArray(value) ? value as Record<string, unknown> : {};
}

function prettyConfigLabel(value: string): string {
  return CONFIG_LABELS[value] || value.replaceAll("_", " ").replace(/\b\w/g, (letter) => letter.toUpperCase());
}

function constraintAnchor(name: string): string {
  return `config-constraint-${name.replace(/[^a-zA-Z0-9_-]/g, "-")}`;
}

function isYearKey(value: string): boolean {
  return /^\d{4}$/.test(value);
}

function isScalar(value: unknown): boolean {
  return value === null || ["string", "number", "boolean"].includes(typeof value);
}

function isYearMatrix(value: Record<string, unknown>): boolean {
  const entries = Object.values(value);
  return entries.length > 0 && entries.every((raw) => {
    if (!raw || typeof raw !== "object" || Array.isArray(raw)) return false;
    const row = objectValue(raw);
    return Object.keys(row).every(isYearKey);
  });
}

function inputNumber(value: string): number | null {
  return value === "" ? null : Number(value);
}

function CustomConstraintsEditor({ value, country, countries, fuelRows, fuelError, onFuelRowsChange, onChange }: { value: Record<string, unknown>; country: string; countries: string[]; fuelRows: InputRow[]; fuelError: string; onFuelRowsChange: (rows: InputRow[]) => void; onChange: (value: Record<string, unknown>) => void }) {
  const countryNames = country === "ALL" ? [...new Set([...countries, ...Object.keys(value)])] : [country];
  const constraintNames = [...new Set([...Object.keys(CONSTRAINT_DEFAULTS), ...countryNames.flatMap((name) => Object.keys(objectValue(value[name])))])];
  const updateConstraint = (countryName: string, constraintName: string, next: Record<string, unknown>) => {
    const countryConstraints = objectValue(value[countryName]);
    onChange({ ...value, [countryName]: { ...countryConstraints, [constraintName]: next } });
  };
  return <div className="config-form constraints-form">
    {constraintNames.map((constraintName) => {
      const activeCountries = countryNames.filter((countryName) => objectValue(objectValue(value[countryName])[constraintName]).activate === true).length;
      return <section className="constraint-group" id={constraintAnchor(constraintName)} key={constraintName}>
        <header><div><p className="eyebrow">Constraint</p><h4>{prettyConfigLabel(constraintName)}</h4></div><span>{activeCountries} of {countryNames.length} active</span></header>
        <div className="constraint-country-list">{countryNames.map((countryName) => <ConstraintCard key={countryName} name={constraintName} country={countryName} value={objectValue(objectValue(value[countryName])[constraintName])} fuelRows={fuelRows} fuelError={fuelError} onFuelRowsChange={onFuelRowsChange} onChange={(next) => updateConstraint(countryName, constraintName, next)} />)}</div>
      </section>;
    })}
    {countryNames.length === 0 && <div className="editor-empty compact">No countries are available for custom constraints.</div>}
  </div>;
}

function ConfigToc({ items }: { items: { id: string; label: string }[] }) {
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
  return <div className={`results-toc ${open ? "open" : ""}`} ref={root}>
    {open && <nav className="results-toc-panel" id="config-section-list" aria-label="Configuration sections on this page">
      <header><div><p className="eyebrow pink">On this page</p><h2>Sections</h2></div><button className="icon-button" onClick={() => setOpen(false)} aria-label="Close section list"><X aria-hidden="true" /></button></header>
      <ol>{items.map((item, index) => <li key={item.id}><a href={`#${item.id}`} onClick={() => setOpen(false)}><span>{String(index + 1).padStart(2, "0")}</span><b>{item.label}</b></a></li>)}</ol>
    </nav>}
    <button className="results-toc-trigger" onClick={() => setOpen((current) => !current)} aria-label="Open section list" aria-expanded={open} aria-controls="config-section-list"><List aria-hidden="true" /></button>
  </div>;
}

function ConstraintCard({ name, country, value, fuelRows, fuelError, onFuelRowsChange, onChange }: { name: string; country: string; value: Record<string, unknown>; fuelRows: InputRow[]; fuelError: string; onFuelRowsChange: (rows: InputRow[]) => void; onChange: (value: Record<string, unknown>) => void }) {
  const editableValue = { ...(CONSTRAINT_DEFAULTS[name] || { activate: false }), ...value };
  const active = editableValue.activate === true;
  const fields = Object.keys(editableValue).filter((key) => key !== "activate");
  return <section className={`constraint-card ${active ? "active" : ""}`}>
    <header>
      <h5>{country}</h5>
      <label className="constraint-toggle"><input type="checkbox" aria-label={`Activate ${prettyConfigLabel(name)} for ${country}`} checked={active} onChange={(event) => onChange({ ...editableValue, activate: event.target.checked })} /><i aria-hidden="true" /><span>{active ? "Active" : "Inactive"}</span></label>
    </header>
    {active && (name === "production_constraint_fuels" ? <FuelProductionLimitsEditor country={country} value={editableValue} rows={fuelRows} error={fuelError} onRowsChange={onFuelRowsChange} onChange={onChange} /> : fields.length ? <ConstraintFields constraintName={name} value={editableValue} onChange={onChange} /> : <p className="constraint-empty">No parameters required.</p>)}
  </section>;
}

function FuelProductionLimitsEditor({ country, value, rows, error, onRowsChange, onChange }: { country: string; value: Record<string, unknown>; rows: InputRow[]; error: string; onRowsChange: (rows: InputRow[]) => void; onChange: (value: Record<string, unknown>) => void }) {
  const countryRows = rows.filter((row) => String(row.country) === country);
  const selected = new Set(Array.isArray(value.fuels) ? value.fuels.map(String) : []);
  const carriers = [...new Set([...countryRows.map((row) => String(row.carrier)), ...selected])].filter(Boolean).sort();
  const years = [...new Set(countryRows.map((row) => String(row.year)).filter(isYearKey))].sort((left, right) => Number(left) - Number(right));
  const rowsByCell = new Map(countryRows.map((row) => [`${String(row.carrier)}\u0000${String(row.year)}`, row]));
  const toggleCarrier = (carrier: string, checked: boolean) => {
    const next = new Set(selected);
    if (checked) next.add(carrier); else next.delete(carrier);
    onChange({ ...value, fuels: carriers.filter((name) => next.has(name)) });
  };
  const updateCell = (rowId: number, raw: string) => onRowsChange(rows.map((row) => row.__row_id === rowId ? { ...row, max_supply__mwh_year: raw } : row));
  if (error) return <div className="notice error fuel-limit-notice">{error}</div>;
  if (!countryRows.length && !selected.size) return <p className="constraint-empty">No fuel supply rows are available for {country}.</p>;
  return <div className="constraint-fields fuel-limit-fields"><div className="config-entry-block constraint-wide"><p className="field-label">Maximum annual fuel supply (MWh/year)</p><p className="field-help">Select the carriers to constrain, then enter their maximum supply for each year. Use <code>inf</code> for no numerical limit. Values are saved to <code>power/fuel_supplies.csv</code>.</p>
    <div className="config-table-wrap"><table className="config-entry-table matrix-table fuel-limit-table"><thead><tr><th>Apply</th><th>Fuel</th>{years.map((year) => <th key={year}>{year}</th>)}</tr></thead><tbody>{carriers.map((carrier) => <tr key={carrier}><td><input className="fuel-limit-check" type="checkbox" aria-label={`Apply fuel production limit to ${carrier} in ${country}`} checked={selected.has(carrier)} onChange={(event) => toggleCarrier(carrier, event.target.checked)} /></td><th scope="row">{carrier}</th>{years.map((year) => { const row = rowsByCell.get(`${carrier}\u0000${year}`); return <td key={year}>{row ? <input aria-label={`${carrier} ${year} maximum supply in MWh per year`} inputMode="decimal" value={String(row.max_supply__mwh_year ?? "")} disabled={!selected.has(carrier)} onChange={(event) => updateCell(row.__row_id, event.target.value)} /> : <span className="fuel-limit-missing" title="No fuel supply row for this year">—</span>}</td>; })}</tr>)}</tbody></table></div>
  </div></div>;
}

function ConstraintFields({ constraintName, value, onChange }: { constraintName: string; value: Record<string, unknown>; onChange: (value: Record<string, unknown>) => void }) {
  return <div className="constraint-fields">{Object.entries(value).filter(([key]) => key !== "activate").map(([key, raw]) => <ConstraintField key={key} constraintName={constraintName} name={key} value={raw} onChange={(next) => onChange({ ...value, [key]: next })} />)}</div>;
}

function ConstraintField({ constraintName, name, value, onChange }: { constraintName: string; name: string; value: unknown; onChange: (value: unknown) => void }) {
  if (Array.isArray(value)) {
    return <label className="field constraint-wide"><span>{prettyConfigLabel(name)}</span><input value={value.join(", ")} placeholder="Comma-separated values" onChange={(event) => onChange(event.target.value.split(",").map((item) => item.trim()).filter(Boolean))} /></label>;
  }
  if (value && typeof value === "object") {
    const mapping = objectValue(value);
    const matrix = (constraintName === "maximum_power_generation_constraint" && name === "value") || isYearMatrix(mapping);
    const years = YEAR_FIELDS.has(name) || (Object.keys(mapping).length > 0 && Object.keys(mapping).every(isYearKey));
    const label = name === "value" && constraintName === "maximum_power_generation_constraint" ? "Generation limits (TWh)" : name === "value" && constraintName === "capacity_factor_constraint" ? "Technology capacity factors" : prettyConfigLabel(name);
    if (matrix) return <YearMatrixEditor label={label} value={mapping} onChange={onChange} />;
    if (years || Object.values(mapping).every(isScalar)) return <MappingTable label={years ? "Year" : label} value={mapping} yearKeys={years} onChange={onChange} />;
    return <section className="constraint-object"><h6>{prettyConfigLabel(name)}</h6><ConstraintFields constraintName={constraintName} value={mapping} onChange={onChange} /></section>;
  }
  if (typeof value === "boolean") {
    return <label className="toggle-row constraint-boolean"><input type="checkbox" checked={value} onChange={(event) => onChange(event.target.checked)} /><span>{prettyConfigLabel(name)}</span></label>;
  }
  if (name === "method") {
    return <label className="field"><span>Method</span><select value={String(value || "static")} onChange={(event) => onChange(event.target.value)}><option value="static">Static</option><option value="dynamic">Dynamic</option></select></label>;
  }
  if (name === "math_symbol") {
    return <label className="field"><span>Comparison</span><select value={String(value || "<=")} onChange={(event) => onChange(event.target.value)}><option value="<=">At most (≤)</option><option value=">=">At least (≥)</option><option value="==">Exactly (=)</option></select></label>;
  }
  if (typeof value === "string") {
    return <label className="field"><span>{prettyConfigLabel(name)}</span><input value={value} onChange={(event) => onChange(event.target.value)} /></label>;
  }
  return <label className="field"><span>{prettyConfigLabel(name)}</span><input type="number" step="any" value={String(value ?? "")} onChange={(event) => onChange(inputNumber(event.target.value))} /></label>;
}

function MappingTable({ label, value, yearKeys = false, allowAdd = true, onChange }: { label: string; value: Record<string, unknown>; yearKeys?: boolean; allowAdd?: boolean; onChange: (value: Record<string, unknown>) => void }) {
  const [newKey, setNewKey] = useState("");
  const entries = Object.entries(value).sort(([left], [right]) => yearKeys ? Number(left) - Number(right) : left.localeCompare(right));
  const remove = (key: string) => { const next = { ...value }; delete next[key]; onChange(next); };
  const add = () => { const key = newKey.trim(); if (!key || Object.hasOwn(value, key)) return; onChange({ ...value, [key]: 0 }); setNewKey(""); };
  return <div className="config-entry-block">
    <div className="config-table-wrap"><table className="config-entry-table"><thead><tr><th>{label}</th><th>Value</th><th><span className="sr-only">Actions</span></th></tr></thead><tbody>{entries.map(([key, raw]) => <tr key={key}><th scope="row">{key}</th><td><input aria-label={`${label} ${key} value`} type={typeof raw === "string" ? "text" : "number"} step="any" value={String(raw ?? "")} onChange={(event) => onChange({ ...value, [key]: typeof raw === "string" ? event.target.value : inputNumber(event.target.value) })} /></td><td><button className="icon-button" aria-label={`Remove ${key}`} onClick={() => remove(key)}><Trash2 aria-hidden="true" /></button></td></tr>)}</tbody></table></div>
    {allowAdd && <div className="config-table-add"><label className="field"><span>New {label.toLowerCase()}</span><input type={yearKeys ? "number" : "text"} value={newKey} placeholder={yearKeys ? "2035" : "Name"} onChange={(event) => setNewKey(event.target.value)} onKeyDown={(event) => { if (event.key === "Enter") { event.preventDefault(); add(); } }} /></label><button className="button secondary" disabled={!newKey.trim() || Object.hasOwn(value, newKey.trim())} onClick={add}><Plus aria-hidden="true" />Add</button></div>}
  </div>;
}

function YearMatrixEditor({ label, value, onChange }: { label: string; value: Record<string, unknown>; onChange: (value: Record<string, unknown>) => void }) {
  const [newTechnology, setNewTechnology] = useState("");
  const [newYear, setNewYear] = useState("");
  const technologies = Object.keys(value).sort();
  const years = [...new Set(Object.values(value).flatMap((raw) => Object.keys(objectValue(raw))).filter(isYearKey))].sort((left, right) => Number(left) - Number(right));
  const setCell = (technology: string, year: string, raw: string) => onChange({ ...value, [technology]: { ...objectValue(value[technology]), [year]: inputNumber(raw) } });
  const removeTechnology = (technology: string) => { const next = { ...value }; delete next[technology]; onChange(next); };
  const addTechnology = () => { const name = newTechnology.trim(); if (!name || Object.hasOwn(value, name)) return; onChange({ ...value, [name]: Object.fromEntries(years.map((year) => [year, null])) }); setNewTechnology(""); };
  const addYear = () => { const year = newYear.trim(); if (!isYearKey(year) || years.includes(year) || !technologies.length) return; onChange(Object.fromEntries(Object.entries(value).map(([technology, raw]) => [technology, { ...objectValue(raw), [year]: null }]))); setNewYear(""); };
  return <div className="config-entry-block matrix-entry-block"><p className="field-label">{label}</p>
    {technologies.length ? <div className="config-table-wrap"><table className="config-entry-table matrix-table"><thead><tr><th>Technology</th>{years.map((year) => <th key={year}>{year}</th>)}<th><span className="sr-only">Actions</span></th></tr></thead><tbody>{technologies.map((technology) => { const row = objectValue(value[technology]); return <tr key={technology}><th scope="row">{technology}</th>{years.map((year) => <td key={year}><input aria-label={`${technology} ${year} value`} type="number" step="any" value={String(row[year] ?? "")} onChange={(event) => setCell(technology, year, event.target.value)} /></td>)}<td><button className="icon-button" aria-label={`Remove ${technology}`} onClick={() => removeTechnology(technology)}><Trash2 aria-hidden="true" /></button></td></tr>; })}</tbody></table></div> : <p className="constraint-empty">Add a technology to start this table.</p>}
    <div className="matrix-add-row"><div className="config-table-add"><label className="field"><span>New technology</span><input value={newTechnology} placeholder="Technology" onChange={(event) => setNewTechnology(event.target.value)} /></label><button className="button secondary" disabled={!newTechnology.trim() || Object.hasOwn(value, newTechnology.trim())} onClick={addTechnology}><Plus aria-hidden="true" />Add</button></div><div className="config-table-add"><label className="field"><span>New year</span><input type="number" value={newYear} placeholder="2035" onChange={(event) => setNewYear(event.target.value)} /></label><button className="button secondary" disabled={!technologies.length || !isYearKey(newYear.trim()) || years.includes(newYear.trim())} onClick={addYear}><Plus aria-hidden="true" />Add</button></div></div>
  </div>;
}

function validateFuelLimits(section: string, value: Record<string, unknown>, rows: InputRow[]): string {
  if (section !== COMBINED_SECTION) return "";
  const constraints = objectValue(value.custom_constraints);
  for (const [country, rawConstraints] of Object.entries(constraints)) {
    const production = objectValue(objectValue(rawConstraints).production_constraint_fuels);
    if (production.activate !== true) continue;
    const fuels = Array.isArray(production.fuels) ? production.fuels.map(String) : [];
    if (!fuels.length) return `Select at least one fuel production limit for ${country}.`;
    for (const fuel of fuels) {
      const fuelRows = rows.filter((row) => String(row.country) === country && String(row.carrier) === fuel);
      if (!fuelRows.length) return `No fuel supply data is available for ${fuel} in ${country}.`;
      for (const row of fuelRows) {
        const raw = row.max_supply__mwh_year;
        if (typeof raw === "string" && raw.trim().toLowerCase() === "inf") continue;
        if (raw === "" || raw === null || !Number.isFinite(Number(raw)) || Number(raw) < 0) return `${fuel} limits for ${country} must be zero or greater, or inf.`;
      }
    }
  }
  return "";
}

function validateSection(section: string, value: Record<string, unknown>): string {
  if (section === "scenario_configs") {
    const snapshots = (value.snapshots || {}) as Record<string, unknown>;
    if (!snapshots.start || !snapshots.end) return "Snapshot start and end dates are required.";
    if (String(snapshots.end) < String(snapshots.start)) return "Snapshot end must be on or after snapshot start.";
    const resolution = (value.resolution || {}) as Record<string, unknown>;
    if (!["nth_hour", "clustered"].includes(String(resolution.method))) return "Choose a supported temporal resolution method.";
    if (resolution.method === "clustered" && Number(resolution.number_of_days) < 1) return "Number of days must be at least 1.";
    if (resolution.method !== "clustered" && Number(resolution.stepsize) < 1) return "Step size must be at least 1.";
    if (!Number.isFinite(Number(value.remove_threshold)) || Number(value.remove_threshold) < 0) return "Remove threshold must be zero or greater.";
  }
  if (section === COMBINED_SECTION || section === "co2_management") {
    const co2 = section === COMBINED_SECTION ? (value.co2_management || {}) as Record<string, unknown> : value;
    for (const raw of Object.values(co2)) {
      const config = (raw || {}) as Record<string, unknown>;
      if (!["co2_cap", "co2_price"].includes(String(config.option))) return "Choose a supported CO₂ instrument for every country.";
      if (Object.values((config.value || {}) as Record<string, unknown>).some((item) => !Number.isFinite(Number(item)))) return "CO₂ year values must be numbers.";
    }
  }
  return "";
}
