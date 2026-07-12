import { useEffect, useMemo, useState } from "react";
import { Check, Plus, RotateCcw, Save, Trash2 } from "lucide-react";
import { getScenarioConfig, saveScenarioConfigSection } from "./api";
import type { InputSelection, ScenarioConfigResponse } from "./types";

const labels: Record<string, string> = { scenario_configs: "Scenario settings", co2_management: "CO₂ management", custom_constraints: "Custom constraints" };

export default function ScenarioConfigEditor({ selection }: { selection: InputSelection }) {
  const [config, setConfig] = useState<ScenarioConfigResponse | null>(null);
  const [section, setSection] = useState("scenario_configs");
  const [draft, setDraft] = useState<Record<string, unknown>>({});
  const [loading, setLoading] = useState(true);
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState("");
  const [success, setSuccess] = useState("");
  const [editorInvalid, setEditorInvalid] = useState(false);

  const load = async () => {
    const controller = new AbortController(); setLoading(true); setError(""); setSuccess("");
    try { const data = await getScenarioConfig(selection, controller.signal); setConfig(data); setDraft(structuredClone(data.sections[section] || {})); }
    catch (reason) { if (!(reason instanceof DOMException && reason.name === "AbortError")) setError(reason instanceof Error ? reason.message : "Could not load the scenario config."); }
    finally { setLoading(false); }
    return () => controller.abort();
  };
  useEffect(() => { void load(); }, [selection.dataset, selection.project, selection.scenario]);
  useEffect(() => { if (config) { setDraft(structuredClone(config.sections[section] || {})); setSuccess(""); setError(""); setEditorInvalid(false); } }, [section]);

  const original = config?.sections[section] || {};
  const dirty = useMemo(() => JSON.stringify(draft) !== JSON.stringify(original), [draft, original]);
  const validationError = useMemo(() => validateSection(section, draft), [section, draft]);
  const invalid = editorInvalid || Boolean(validationError);
  useEffect(() => { const warn = (event: BeforeUnloadEvent) => { if (dirty) event.preventDefault(); }; window.addEventListener("beforeunload", warn); return () => window.removeEventListener("beforeunload", warn); }, [dirty]);
  const save = async () => {
    if (!config || !dirty) return; setSaving(true); setError(""); setSuccess("");
    try { const data = await saveScenarioConfigSection(selection, section, config.revision, draft); setConfig(data); setDraft(structuredClone(data.sections[section] || {})); setSuccess(`${labels[section]} saved directly to scenario_config.yaml.`); }
    catch (reason) { setError(reason instanceof Error ? reason.message : "Could not save this section."); }
    finally { setSaving(false); }
  };

  return <>
    <section className="page-title editor-title"><div><p className="eyebrow pink">Scenario configuration</p><h1>Configure {selection.scenario}</h1><p>Edit the model settings and constraints stored in this scenario’s YAML file.</p></div></section>
    <div className="config-layout">
      <nav className="config-nav" aria-label="Configuration sections">{Object.keys(labels).map((name) => <button key={name} className={section === name ? "active" : ""} onClick={() => { if (!dirty || window.confirm("Discard unsaved changes in this section?")) setSection(name); }}>{labels[name]}</button>)}</nav>
      <section className="editor-panel config-panel">
        <header className="editor-panel-head"><div><p className="eyebrow">{selection.project} · {selection.scenario}</p><h2>{labels[section]}</h2>{config && <code>{config.path}</code>}</div><div className="editor-actions"><button className="button secondary" disabled={!dirty || saving} onClick={() => { setDraft(structuredClone(original)); setEditorInvalid(false); }}><RotateCcw aria-hidden="true" />Discard</button><button className="button primary" disabled={!dirty || saving || invalid} onClick={save}><Save aria-hidden="true" />{saving ? "Saving…" : "Save changes"}</button></div></header>
        {error && <div className="notice error">{error}<button onClick={() => void load()}>Reload</button></div>}
        {success && <div className="notice success"><Check aria-hidden="true" />{success}</div>}
        {validationError && <div className="notice error">{validationError}</div>}
        {loading ? <div className="editor-loading"><span className="spinner" />Reading YAML…</div> : section === "scenario_configs" ? <ScenarioSettings value={draft} onChange={setDraft} /> : section === "co2_management" ? <Co2Editor value={draft} onChange={setDraft} /> : <JsonSectionEditor value={draft} onChange={setDraft} onInvalid={setEditorInvalid} />}
      </section>
    </div>
  </>;
}

function ScenarioSettings({ value, onChange }: { value: Record<string, unknown>; onChange: (value: Record<string, unknown>) => void }) {
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
    <div className="form-section"><h3>Interest rates</h3><p className="field-help">Country-specific decimal rates; 0.05 means 5%.</p><KeyValueEditor value={interest} onChange={(next) => onChange({ ...value, interest: next })} /></div>
  </div>;
}

function Co2Editor({ value, onChange }: { value: Record<string, unknown>; onChange: (value: Record<string, unknown>) => void }) {
  const updateCountry = (country: string, next: Record<string, unknown>) => onChange({ ...value, [country]: next });
  return <div className="config-form"><p className="field-help">Choose a carbon cap or price for each country and edit its year-specific values.</p>{Object.entries(value).map(([country, raw]) => { const config = (raw || {}) as Record<string, unknown>; return <article className="country-config" key={country}><header><h3>{country}</h3><label className="field compact"><span>Instrument</span><select value={String(config.option || "co2_cap")} onChange={(event) => updateCountry(country, { ...config, option: event.target.value })}><option value="co2_cap">CO₂ cap</option><option value="co2_price">CO₂ price</option></select></label></header><KeyValueEditor value={(config.value || {}) as Record<string, unknown>} onChange={(next) => updateCountry(country, { ...config, value: next })} /></article>; })}</div>;
}

function KeyValueEditor({ value, onChange }: { value: Record<string, unknown>; onChange: (value: Record<string, unknown>) => void }) {
  const [newKey, setNewKey] = useState("");
  const add = () => { const key = newKey.trim(); if (!key || Object.hasOwn(value, key)) return; onChange({ ...value, [key]: 0 }); setNewKey(""); };
  const remove = (key: string) => { const next = { ...value }; delete next[key]; onChange(next); };
  return <div className="mapping-editor"><div className="key-value-grid">{Object.entries(value).map(([key, raw]) => <div className="mapping-field" key={key}><label className="field"><span>{key}</span><input type="number" step="any" value={String(raw ?? "")} onChange={(event) => onChange({ ...value, [key]: Number(event.target.value) })} /></label><button className="icon-button" aria-label={`Remove ${key}`} onClick={() => remove(key)}><Trash2 aria-hidden="true" /></button></div>)}</div><div className="mapping-add"><label className="field"><span>New key</span><input value={newKey} onChange={(event) => setNewKey(event.target.value)} onKeyDown={(event) => { if (event.key === "Enter") { event.preventDefault(); add(); } }} /></label><button className="button secondary" disabled={!newKey.trim() || Object.hasOwn(value, newKey.trim())} onClick={add}><Plus aria-hidden="true" />Add</button></div></div>;
}

function JsonSectionEditor({ value, onChange, onInvalid }: { value: Record<string, unknown>; onChange: (value: Record<string, unknown>) => void; onInvalid: (invalid: boolean) => void }) {
  const [text, setText] = useState(() => JSON.stringify(value, null, 2)); const [error, setError] = useState("");
  useEffect(() => setText(JSON.stringify(value, null, 2)), [value]);
  const edit = (next: string) => { setText(next); try { const parsed = JSON.parse(next); if (!parsed || Array.isArray(parsed) || typeof parsed !== "object") throw new Error("The root must be an object."); setError(""); onInvalid(false); onChange(parsed); } catch (reason) { setError(reason instanceof Error ? reason.message : "Invalid JSON"); onInvalid(true); } };
  return <div className="config-form"><p className="field-help">Edit the complete nested constraint mapping. JSON year keys are restored as numeric YAML keys when saved.</p><label className="field json-field"><span>Constraint configuration</span><textarea spellCheck="false" value={text} onChange={(event) => edit(event.target.value)} aria-invalid={Boolean(error)} /></label>{error && <p className="validation-error">{error}</p>}</div>;
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
  if (section === "co2_management") {
    for (const raw of Object.values(value)) {
      const config = (raw || {}) as Record<string, unknown>;
      if (!["co2_cap", "co2_price"].includes(String(config.option))) return "Choose a supported CO₂ instrument for every country.";
      if (Object.values((config.value || {}) as Record<string, unknown>).some((item) => !Number.isFinite(Number(item)))) return "CO₂ year values must be numbers.";
    }
  }
  return "";
}
