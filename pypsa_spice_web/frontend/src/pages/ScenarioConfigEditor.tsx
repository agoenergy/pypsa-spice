import { useEffect, useRef, useState } from "react";
import { createPortal } from "react-dom";
import { ClipboardCheck, Code2, List, Plus, Settings2, Trash2, X } from "lucide-react";
import "./ScenarioConfigEditor.css";
import Co2Editor from "./Co2Editor";
import PageHeader from "../components/PageHeader";
import RunModel from "../components/RunModel";
import SaveDiscardActions from "../components/SaveDiscardActions";
import sidebarStyles from "../components/Sidebar.module.scss";
import type { InputRow, InputSelection, ModelRunStatus } from "../types";
import { confirmDiscardChanges } from "../utility";
import { MappingTable } from "./ScenarioConfigControls";
import {
  COMBINED_SECTION,
  CONSTRAINT_DEFAULTS,
  SECTION_LABELS,
  YEAR_FIELDS,
  constraintAnchor,
  inputNumber,
  isScalar,
  isYearKey,
  isYearMatrix,
  normaliseSection,
  objectValue,
  prettyConfigLabel,
} from "./ScenarioConfigEditorUtils";
import ScenarioSettings from "./ScenarioSettings";
import useScenarioConfigEditor from "./useScenarioConfigEditor";

const icons = { scenario_configs: Settings2, [COMBINED_SECTION]: Code2, review_run: ClipboardCheck };

export default function ScenarioConfigEditor({ selection, country, onNavigate, onRunStatusChange, onOpenResults }: { selection: InputSelection; country: string; onNavigate: () => void; onRunStatusChange?: (status: ModelRunStatus | null) => void; onOpenResults: (runName: string, dataset: string, project: string) => void }) {
  const [section, setSection] = useState(() => {
    const requested = new URLSearchParams(window.location.search).get("step") || "scenario_configs";
    return normaliseSection(requested);
  });
  const [navigationTarget, setNavigationTarget] = useState<HTMLElement | null>(null);
  const editor = useScenarioConfigEditor(selection, section);
  useEffect(() => { setNavigationTarget(document.getElementById("config-section-tabs")); }, []);

  const chooseSection = (name: string) => {
    if (!confirmDiscardChanges()) return;
    setSection(name); onNavigate();
    const params = new URLSearchParams(window.location.search);
    params.set("view", "configure"); params.set("step", name); params.delete("section");
    window.history.pushState(null, "", `?${params.toString()}`);
    window.scrollTo({ top: 0, behavior: "smooth" });
  };

  const navigation = navigationTarget && createPortal(<nav className={sidebarStyles["submenu-list"]} aria-label="Configuration and run pages">{Object.keys(SECTION_LABELS).map((name) => { const SectionIcon = icons[name as keyof typeof icons]; return <button key={name} className={`${sidebarStyles["submenu-item"]} ${section === name ? sidebarStyles["active"] : ""}`} onClick={() => chooseSection(name)} aria-current={section === name ? "page" : undefined}><SectionIcon aria-hidden="true" /><b>{SECTION_LABELS[name]}</b></button>; })}</nav>, navigationTarget);

  if (section === "review_run") return <>{navigation}<RunModel selection={selection} onEditConfiguration={() => chooseSection("scenario_configs")} onRunStatusChange={onRunStatusChange} onOpenResults={onOpenResults} /></>;

  return <>
    {navigation}
    <PageHeader title={`Configure ${selection.scenario}`} />
    <div className="config-layout">
      <section className={`editor-panel config-panel ${section === COMBINED_SECTION ? "combined-config-panel" : ""}`}>
        <header className="editor-panel-head"><div><h2>{SECTION_LABELS[section]}</h2>{editor.config && <code>{editor.config.path}</code>}</div></header>
        {editor.error && <div className="notice error">{editor.error}<button onClick={() => void editor.reload()}>Reload</button></div>}
        {!editor.loading && editor.validationError && <div className="notice error">{editor.validationError}</div>}
        {editor.loading ? <div className="editor-loading"><span className="spinner" />Reading configuration…</div> : section === "scenario_configs" ? <ScenarioSettings value={editor.draft} country={country} onChange={editor.setDraft} /> : <CombinedConstraintsEditor value={editor.draft} country={country} fuelRows={editor.fuelRows} fuelError={editor.fuelError} onFuelRowsChange={editor.setFuelRows} onChange={editor.setDraft} />}
        <SaveDiscardActions
          hasChanges={editor.dirty}
          saving={editor.saving}
          saveDisabled={Boolean(editor.validationError)}
          status={editor.success}
          floating
          avoidSideControl={section === COMBINED_SECTION}
          onDiscard={editor.discard}
          onSave={() => void editor.save()}
        />
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
    <section className="combined-config-section" id="config-co2-management" aria-labelledby="co2-management-heading"><header><h3 id="co2-management-heading">CO₂ management</h3><p>Country carbon cap or price by year.</p></header><Co2Editor value={co2} country={country} onChange={(next) => onChange({ ...value, co2_management: next })} /></section>
    <section className="combined-config-section" aria-labelledby="custom-constraints-heading"><header><h3 id="custom-constraints-heading">Custom constraints</h3><p>Activate constraints and edit their parameters directly.</p></header><CustomConstraintsEditor value={constraints} country={country} countries={Object.keys(co2)} fuelRows={fuelRows} fuelError={fuelError} onFuelRowsChange={onFuelRowsChange} onChange={(next) => onChange({ ...value, custom_constraints: next })} /></section>
  </div><ConfigToc items={tocItems} /></>;
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
        <header><h4>{prettyConfigLabel(constraintName)}</h4><span>{activeCountries} of {countryNames.length} active</span></header>
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
      <header><h2>Sections</h2><button className="icon-button" onClick={() => setOpen(false)} aria-label="Close section list"><X aria-hidden="true" /></button></header>
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
