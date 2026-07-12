import { useEffect, useId, useMemo, useState } from "react";
import { createPortal } from "react-dom";
import { AlertTriangle, Check, ChevronLeft, ChevronRight, Cpu, RotateCcw, Save, Search, Table2 } from "lucide-react";
import { getInputTable, saveInputTable } from "./api";
import type { InputCatalog, InputCell, InputRow, InputSelection, InputTableDefinition, InputTableResponse, InputTechnology } from "./types";
import { confirmDiscardChanges, setEditorDirty } from "./dirtyState";

const PAGE_SIZE = 100;

export default function InputEditor({ catalog, selection, country }: { catalog: InputCatalog; selection: InputSelection; country: string }) {
  const [view, setView] = useState<"table" | "technology">("technology");
  const [scope, setScope] = useState<"global" | "scenario">("global");
  const [sector, setSector] = useState("power");
  const project = catalog.datasets.find((item) => item.name === selection.dataset)?.projects.find((item) => item.name === selection.project);
  const technologies = (project?.technologies || []).filter((item) => item.sector === sector && (country === "ALL" || item.countries.includes(country)));
  const [technologyId, setTechnologyId] = useState("");
  const [menuTarget, setMenuTarget] = useState<HTMLElement | null>(null);
  const definitions = scope === "global" ? catalog.global_tables : (catalog.sector_tables[sector] || []);
  const [tableId, setTableId] = useState(catalog.global_tables.find((item) => item.id === "Technologies")?.id || definitions[0]?.id || "");

  useEffect(() => {
    const available = scope === "global" ? catalog.global_tables : (catalog.sector_tables[sector] || []);
    const preferred = scope === "global" ? "Technologies" : ({ power: "Power_generators", industry: "Heat_generators", transport: "Transport_loads" } as Record<string, string>)[sector];
    if (!available.some((item) => item.id === tableId)) setTableId(available.find((item) => item.id === preferred)?.id || available[0]?.id || "");
  }, [catalog, scope, sector, tableId]);

  useEffect(() => {
    if (!technologies.some((item) => item.id === technologyId)) setTechnologyId(technologies[0]?.id || "");
  }, [technologies, technologyId]);
  useEffect(() => { setMenuTarget(document.getElementById("input-table-menu")); }, []);

  const definition = definitions.find((item) => item.id === tableId);
  const technology = technologies.find((item) => item.id === technologyId);
  const guarded = (action: () => void) => { if (confirmDiscardChanges()) action(); };
  return <>
    {menuTarget && createPortal(<nav className="section-tabs top-section-tabs input-section-tabs" aria-label="Input sections" role="tablist">
      <button className={`section-tab ${view === "technology" ? "active" : ""}`} onClick={() => guarded(() => setView("technology"))} role="tab" aria-selected={view === "technology"}><Cpu aria-hidden="true" /><b>By technology</b></button>
      <button className={`section-tab ${view === "table" ? "active" : ""}`} onClick={() => guarded(() => setView("table"))} role="tab" aria-selected={view === "table"}><Table2 aria-hidden="true" /><b>By table</b></button>
    </nav>, menuTarget)}
    <section className="page-title editor-title"><div><p className="eyebrow pink">Model inputs</p><h1>Input data</h1><p>Explore and edit the model’s source CSV files by table or by technology. Changes are written only when you select Save changes.</p></div></section>
    <section className="editor-primary-select" aria-label={view === "table" ? "Table selection" : "Technology selection"}>
      <div><p className="eyebrow pink">Current selection</p><h2>{view === "table" ? "Choose a table" : "Choose a technology"}</h2><p>{view === "table" ? "Select the source CSV you want to inspect and edit." : "Select a technology to review its shared and scenario-specific inputs."}</p></div>
      <div className="editor-selection-fields">
        {view === "table" && <div className="segmented" role="tablist" aria-label="Input scope"><button className={scope === "global" ? "active" : ""} onClick={() => guarded(() => setScope("global"))} role="tab" aria-selected={scope === "global"}>Global inputs</button><button className={scope === "scenario" ? "active" : ""} onClick={() => guarded(() => setScope("scenario"))} role="tab" aria-selected={scope === "scenario"}>Scenario inputs</button></div>}
        {(view === "technology" || scope === "scenario") && <label className="field sector-select"><span>Sector</span><select value={sector} onChange={(event) => guarded(() => setSector(event.target.value))}>{["power", "industry", "transport"].map((item) => <option value={item} key={item}>{item}</option>)}</select></label>}
        {view === "table" ? <label className="field primary-select"><span>Table</span><select value={tableId} onChange={(event) => guarded(() => setTableId(event.target.value))}>{definitions.map((item) => <option value={item.id} key={item.id}>{item.label}</option>)}</select></label> : <label className="field primary-select"><span>Technology</span><select value={technologyId} onChange={(event) => guarded(() => setTechnologyId(event.target.value))}>{technologies.map((item) => <option value={item.id} key={item.id}>{item.label} ({item.id})</option>)}</select></label>}
      </div>
    </section>
    {view === "table" ? definition ? <TableEditor key={`${selection.dataset}:${selection.project}:${selection.scenario}:${scope}:${sector}:${definition.id}`} definition={definition} selection={selection} country={country} /> : <div className="editor-empty">No configured table is available for this selection.</div> : technology ? <TechnologyEditor catalog={catalog} selection={selection} sector={sector} technology={technology} country={country} /> : <div className="editor-empty">No mapped technology is available for this sector and country.</div>}
  </>;
}

function TechnologyEditor({ catalog, selection, sector, technology, country }: { catalog: InputCatalog; selection: InputSelection; sector: string; technology: InputTechnology; country: string }) {
  const globalDefinitions = catalog.global_tables.filter((item) => item.id !== "Demand_Profiles");
  const scenarioDefinitions = catalog.sector_tables[sector] || [];
  return <div className="technology-view">
    <section className="technology-summary"><div><p className="eyebrow pink">Selected technology</p><h2>{technology.label}</h2><code>{technology.id}</code></div><dl><div><dt>PyPSA class</dt><dd>{technology.classes.join(", ") || "—"}</dd></div><div><dt>Carrier</dt><dd>{technology.carriers.join(", ") || "—"}</dd></div></dl></section>
    <section className="technology-group"><header><p className="eyebrow">Shared assumptions</p><h2>Global input</h2><span>Changes here apply to every scenario in this project{country === "ALL" ? "" : ` for ${country}`}.</span></header><div className="technology-panels">{globalDefinitions.map((definition) => <TableEditor key={`global:${definition.id}:${technology.id}`} definition={definition} selection={selection} technology={technology} country={country} hideWhenEmpty />)}</div></section>
    <section className="technology-group"><header><p className="eyebrow">{selection.scenario}</p><h2>Scenario input</h2><span>Regional assets and constraints that reference this technology{country === "ALL" ? "" : ` in ${country}`}.</span></header><div className="technology-panels">{scenarioDefinitions.map((definition) => <TableEditor key={`scenario:${definition.id}:${technology.id}`} definition={definition} selection={selection} technology={technology} country={country} hideWhenEmpty />)}</div></section>
  </div>;
}

function technologyMatches(row: InputRow, definition: InputTableDefinition, technology: InputTechnology): boolean {
  if (definition.id === "Direct_air_capture") return technology.id === "DAC";
  const raw = String(row[definition.filter_col] ?? "").trim();
  if (!raw) return false;
  if (technology.id === "PEVCH" && (raw.startsWith("EVCH") || raw.startsWith("EVST"))) return true;
  if (definition.id.toLowerCase().includes("decommission")) {
    const parts = raw.split("_");
    return parts[parts.length - 1] === technology.id;
  }
  if (definition.filter_col === "carrier") return technology.carriers.includes(raw);
  return raw === technology.id || raw.split(/[;,|]/).map((value) => value.trim()).includes(technology.id);
}

function TableEditor({ definition, selection, country, technology, hideWhenEmpty = false }: { definition: InputTableDefinition; selection: InputSelection; country: string; technology?: InputTechnology; hideWhenEmpty?: boolean }) {
  const editorId = useId();
  const [table, setTable] = useState<InputTableResponse | null>(null);
  const [rows, setRows] = useState<InputRow[]>([]);
  const [changes, setChanges] = useState<Map<string, InputCell>>(new Map());
  const [query, setQuery] = useState("");
  const [filter, setFilter] = useState("ALL");
  const [page, setPage] = useState(0);
  const [loading, setLoading] = useState(true);
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState("");
  const [success, setSuccess] = useState("");

  const load = async () => {
    const controller = new AbortController();
    setLoading(true); setError(""); setSuccess("");
    try { const data = await getInputTable(selection, definition, controller.signal); setTable(data); setRows(data.rows); setChanges(new Map()); }
    catch (reason) { if (!(reason instanceof DOMException && reason.name === "AbortError")) setError(reason instanceof Error ? reason.message : "Could not load this table."); }
    finally { setLoading(false); }
    return () => controller.abort();
  };

  useEffect(() => { void load(); }, [definition.id, selection.dataset, selection.project, selection.scenario]);
  useEffect(() => { setPage(0); }, [query, filter]);
  useEffect(() => {
    const warn = (event: BeforeUnloadEvent) => { if (changes.size) event.preventDefault(); };
    window.addEventListener("beforeunload", warn); return () => window.removeEventListener("beforeunload", warn);
  }, [changes.size]);
  useEffect(() => { setEditorDirty(editorId, changes.size > 0); return () => setEditorDirty(editorId, false); }, [editorId, changes.size]);

  const filterColumn = table?.filter_column;
  const showFilter = Boolean(filterColumn && !technology);
  const countryRows = useMemo(() => rows.filter((row) => country === "ALL" || !("country" in row) || String(row.country) === country), [rows, country]);
  const filterOptions = useMemo(() => {
    if (!showFilter || !filterColumn) return [];
    return [...new Set(countryRows.map((row) => String(row[filterColumn] ?? "")).filter(Boolean))].sort();
  }, [countryRows, filterColumn, showFilter]);
  const technologyRows = useMemo(() => technology ? countryRows.filter((row) => technologyMatches(row, definition, technology)) : countryRows, [countryRows, definition, technology]);
  const filtered = useMemo(() => technologyRows.filter((row) => {
    if (showFilter && filterColumn && filter !== "ALL" && String(row[filterColumn]) !== filter) return false;
    const needle = query.trim().toLowerCase();
    return !needle || Object.values(row).some((value) => String(value ?? "").toLowerCase().includes(needle));
  }), [technologyRows, showFilter, filterColumn, filter, query]);
  const pageCount = Math.max(1, Math.ceil(filtered.length / PAGE_SIZE));
  const visibleRows = filtered.slice(page * PAGE_SIZE, (page + 1) * PAGE_SIZE);

  const edit = (rowId: number, column: string, value: InputCell) => {
    const key = `${rowId}:${column}`;
    setRows((current) => current.map((row) => row.__row_id === rowId ? { ...row, [column]: value } : row));
    setChanges((current) => { const next = new Map(current); next.set(key, value); return next; });
    setSuccess("");
  };
  const discard = () => { if (table) setRows(table.rows); setChanges(new Map()); setSuccess(""); };
  const save = async () => {
    if (!table || !changes.size) return;
    setSaving(true); setError(""); setSuccess("");
    try {
      const payload = [...changes].map(([key, value]) => { const split = key.indexOf(":"); return { row: Number(key.slice(0, split)), column: key.slice(split + 1), value }; });
      const data = await saveInputTable(selection, definition, table.revision, payload);
      setTable(data); setRows(data.rows); setChanges(new Map()); setSuccess(`Saved ${payload.length} ${payload.length === 1 ? "cell" : "cells"} directly to the CSV.`);
    } catch (reason) { setError(reason instanceof Error ? reason.message : "Could not save changes."); }
    finally { setSaving(false); }
  };

  if (hideWhenEmpty && (loading || error || technologyRows.length === 0)) return null;
  return <section className={`editor-panel${technology ? " technology-panel" : ""}`}>
    <header className="editor-panel-head"><div><p className="eyebrow">{definition.scope === "global" ? "Global source" : `${definition.sector} · ${selection.scenario}`}</p><h2>{definition.label}</h2>{table && <code>{table.path}</code>}</div><div className="editor-actions"><button className="button secondary" disabled={!changes.size || saving} onClick={discard}><RotateCcw aria-hidden="true" />Discard</button><button className="button primary" disabled={!changes.size || saving || definition.timeseries} onClick={save}><Save aria-hidden="true" />{saving ? "Saving…" : `Save changes${changes.size ? ` (${changes.size})` : ""}`}</button></div></header>
    {definition.timeseries && <div className="editor-warning"><AlertTriangle aria-hidden="true" /><span><b>Read-only timeseries.</b> Large hourly inputs are shown for inspection and edited locally outside the browser.</span></div>}
    {error && <div className="notice error">{error}<button onClick={() => void load()}>Reload</button></div>}
    {success && <div className="notice success"><Check aria-hidden="true" />{success}</div>}
    {loading ? <div className="editor-loading"><span className="spinner" />Reading CSV…</div> : table && <>
      <div className="table-tools"><label className="search"><Search aria-hidden="true" /><input value={query} onChange={(event) => setQuery(event.target.value)} type="search" placeholder="Find a row or value" /></label>{showFilter && filterColumn && <label className="field compact"><span>{filterColumn.replaceAll("_", " ")}</span><select value={filter} onChange={(event) => setFilter(event.target.value)}><option value="ALL">All values</option>{filterOptions.map((value) => <option key={value}>{value}</option>)}</select></label>}<span className="row-count">{filtered.length.toLocaleString()} of {table.total_rows.toLocaleString()} rows</span></div>
      <div className="editable-table-wrap"><table className="editable-table"><thead><tr>{table.columns.map((column) => <th key={column.name}><span>{column.label}</span>{column.editable && <small>Editable</small>}</th>)}</tr></thead><tbody>{visibleRows.map((row) => <tr key={row.__row_id}>{table.columns.map((column) => <td key={column.name} className={column.editable ? "editable-cell" : "locked-cell"}>{column.editable ? <CellEditor value={row[column.name]} kind={column.kind} onChange={(value) => edit(row.__row_id, column.name, value)} /> : <span>{String(row[column.name] ?? "")}</span>}</td>)}</tr>)}</tbody></table></div>
      <footer className="table-pagination"><span>Page {Math.min(page + 1, pageCount)} of {pageCount}{table.truncated ? " · first 10,000 rows loaded" : ""}</span><div><button className="icon-button" aria-label="Previous page" disabled={page === 0} onClick={() => setPage((current) => Math.max(0, current - 1))}><ChevronLeft /></button><button className="icon-button" aria-label="Next page" disabled={page >= pageCount - 1} onClick={() => setPage((current) => Math.min(pageCount - 1, current + 1))}><ChevronRight /></button></div></footer>
    </>}
  </section>;
}

function CellEditor({ value, kind, onChange }: { value: InputCell; kind: string; onChange: (value: InputCell) => void }) {
  if (kind === "boolean") return <label className="cell-check"><input type="checkbox" checked={Boolean(value)} onChange={(event) => onChange(event.target.checked)} /><span aria-hidden="true" /></label>;
  return <input className="cell-input" value={String(value ?? "")} inputMode={kind === "number" ? "decimal" : undefined} onChange={(event) => onChange(event.target.value)} aria-label="Editable cell" />;
}
