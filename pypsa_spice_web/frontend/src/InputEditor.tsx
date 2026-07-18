import { useDeferredValue, useEffect, useId, useMemo, useRef, useState } from "react";
import { createPortal } from "react-dom";
import { AlertTriangle, Check, ChevronLeft, ChevronRight, Cpu, RotateCcw, Save, Search, Table2 } from "lucide-react";
import { getInputTable, saveInputTable } from "./api";
import type { InputCatalog, InputCell, InputRow, InputSelection, InputTableDefinition, InputTableResponse, InputTechnology } from "./types";
import { confirmDiscardChanges, setEditorDirty } from "./dirtyState";
import PageHeader from "./PageHeader";

const PAGE_SIZE = 100;

export default function InputEditor({ catalog, selection, onNavigate }: { catalog: InputCatalog; selection: InputSelection; onNavigate: () => void }) {
  const [view, setView] = useState<"table" | "technology">("technology");
  const [sector, setSector] = useState("power");
  const project = catalog.datasets.find((item) => item.name === selection.dataset)?.projects.find((item) => item.name === selection.project);
  const technologies = (project?.technologies || []).filter((item) => item.sector === sector);
  const [technologyId, setTechnologyId] = useState("");
  const [menuTarget, setMenuTarget] = useState<HTMLElement | null>(null);
  const [topbarTarget, setTopbarTarget] = useState<HTMLElement | null>(null);
  const definitions = useMemo(() => [...catalog.global_tables, ...(catalog.sector_tables[sector] || [])], [catalog, sector]);
  const [tableId, setTableId] = useState(catalog.global_tables.find((item) => item.id === "Technologies")?.id || definitions[0]?.id || "");

  useEffect(() => {
    if (definitions.some((item) => item.id === tableId)) return;
    const preferred = ({ power: "Power_generators", industry: "Heat_generators", transport: "Transport_loads" } as Record<string, string>)[sector];
    setTableId(definitions.find((item) => item.id === preferred)?.id || catalog.global_tables.find((item) => item.id === "Technologies")?.id || definitions[0]?.id || "");
  }, [catalog, definitions, sector, tableId]);

  useEffect(() => {
    if (!technologies.some((item) => item.id === technologyId)) setTechnologyId(technologies[0]?.id || "");
  }, [technologies, technologyId]);
  useEffect(() => {
    setMenuTarget(document.getElementById("input-table-menu"));
    setTopbarTarget(document.getElementById("input-topbar-controls"));
  }, []);

  const definition = definitions.find((item) => item.id === tableId);
  const technology = technologies.find((item) => item.id === technologyId) || technologies[0];
  const guarded = (action: () => void) => { if (confirmDiscardChanges()) action(); };
  return <>
    {menuTarget && createPortal(<nav className="sidebar-submenu-list" aria-label="Input pages">
      <button className={`sidebar-submenu-item ${view === "technology" ? "active" : ""}`} onClick={() => guarded(() => { setView("technology"); onNavigate(); })} aria-current={view === "technology" ? "page" : undefined}><Cpu aria-hidden="true" /><b>By technology</b></button>
      <button className={`sidebar-submenu-item ${view === "table" ? "active" : ""}`} onClick={() => guarded(() => { setView("table"); onNavigate(); })} aria-current={view === "table" ? "page" : undefined}><Table2 aria-hidden="true" /><b>By table</b></button>
    </nav>, menuTarget)}
    {topbarTarget && createPortal(<>
      <label className="context-control input-sector-control"><span>Sector</span><select value={sector} onChange={(event) => guarded(() => setSector(event.target.value))}>{["power", "industry", "transport"].map((item) => <option value={item} key={item}>{item}</option>)}</select></label>
      {view === "technology" ? <label className="context-control input-technology-control"><span>Technology</span><select value={technologyId} onChange={(event) => guarded(() => setTechnologyId(event.target.value))}>{technologies.map((item) => <option value={item.id} key={item.id}>{item.label} ({item.id})</option>)}</select></label> : <label className="context-control input-table-control"><span>Table</span><select value={tableId} onChange={(event) => guarded(() => setTableId(event.target.value))}>{definitions.map((item) => <option value={item.id} key={item.id}>{item.label}</option>)}</select></label>}
    </>, topbarTarget)}
    {view === "technology" ? technology && <TechnologyTitle technology={technology} /> : definition && <TableTitle definition={definition} />}
    {view === "table" ? definition ? <TableView key={`${selection.dataset}:${selection.project}:${selection.scenario}:${sector}:${definition.id}`} definition={definition} selection={selection} /> : <div className="editor-empty">No configured table is available for this selection.</div> : technology ? <TechnologyEditor catalog={catalog} selection={selection} sector={sector} technology={technology} /> : <div className="editor-empty">No mapped technology is available for this sector.</div>}
  </>;
}

function TechnologyTitle({ technology }: { technology: InputTechnology }) {
  return <PageHeader title={technology.label} className="selection-title"><dl><div><dt>PyPSA class</dt><dd>{technology.classes.join(", ") || "—"}</dd></div><div><dt>Carrier</dt><dd>{technology.carriers.join(", ") || "—"}</dd></div></dl></PageHeader>;
}

function TableTitle({ definition }: { definition: InputTableDefinition }) {
  return <PageHeader title={definition.label} className="selection-title"><dl><div><dt>Scope</dt><dd>{definition.scope === "global" ? "Global input" : "Scenario input"}</dd></div><div><dt>Sector</dt><dd>{definition.sector || "All sectors"}</dd></div></dl></PageHeader>;
}

function TableView({ definition, selection }: { definition: InputTableDefinition; selection: InputSelection }) {
  const global = definition.scope === "global";
  return <div className="table-view">
    <section className="technology-group"><header><h2>{global ? "Global input" : "Scenario input"}</h2><span>{global ? "Changes here apply to every country and every scenario in this project." : "Assets and constraints for this scenario. Country filters appear only on tables with country-specific rows."}</span></header><TableEditor definition={definition} selection={selection} /></section>
  </div>;
}

function TechnologyEditor({ catalog, selection, sector, technology }: { catalog: InputCatalog; selection: InputSelection; sector: string; technology: InputTechnology }) {
  const globalDefinitions = catalog.global_tables.filter((item) => item.id !== "Demand_Profiles");
  const scenarioDefinitions = catalog.sector_tables[sector] || [];
  return <div className="technology-view">
    <section className="technology-group"><header><h2>Global input</h2><span>Changes here apply to every country and every scenario in this project.</span></header><div className="technology-panels">{globalDefinitions.map((definition) => <TableEditor key={`${selection.dataset}:${selection.project}:global:${definition.id}:${technology.id}`} definition={definition} selection={selection} technology={technology} hideWhenEmpty />)}</div></section>
    <section className="technology-group"><header><h2>Scenario input</h2><span>Assets and constraints for this scenario. Country filters appear only on tables with country-specific rows.</span></header><div className="technology-panels">{scenarioDefinitions.map((definition) => <TableEditor key={`${selection.dataset}:${selection.project}:${selection.scenario}:${definition.id}:${technology.id}`} definition={definition} selection={selection} technology={technology} hideWhenEmpty />)}</div></section>
  </div>;
}

function TableEditor({ definition, selection, technology, hideWhenEmpty = false }: { definition: InputTableDefinition; selection: InputSelection; technology?: InputTechnology; hideWhenEmpty?: boolean }) {
  const editorId = useId();
  const [table, setTable] = useState<InputTableResponse | null>(null);
  const [rows, setRows] = useState<InputRow[]>([]);
  const [changes, setChanges] = useState<Map<string, InputCell>>(new Map());
  const changesRef = useRef(changes);
  const [query, setQuery] = useState("");
  const deferredQuery = useDeferredValue(query);
  const [filter, setFilter] = useState("ALL");
  const [country, setCountry] = useState("ALL");
  const [page, setPage] = useState(0);
  const [loading, setLoading] = useState(true);
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState("");
  const [success, setSuccess] = useState("");

  const load = async (signal?: AbortSignal) => {
    setLoading(true); setError(""); setSuccess("");
    try {
      const data = await getInputTable(selection, definition, {
        technology,
        country,
        filterValue: filter,
        query: deferredQuery,
        offset: page * PAGE_SIZE,
        limit: PAGE_SIZE,
      }, signal);
      const pending = changesRef.current;
      const loadedRows = data.rows.map((row) => {
        const next = { ...row };
        for (const [key, value] of pending) {
          const separator = key.indexOf(":");
          if (Number(key.slice(0, separator)) === row.__row_id) next[key.slice(separator + 1)] = value;
        }
        return next;
      });
      setTable(data); setRows(loadedRows);
    }
    catch (reason) { if (!(reason instanceof DOMException && reason.name === "AbortError")) setError(reason instanceof Error ? reason.message : "Could not load this table."); }
    finally { if (!signal?.aborted) setLoading(false); }
  };

  useEffect(() => {
    const controller = new AbortController();
    void load(controller.signal);
    return () => controller.abort();
  }, [definition.id, selection.dataset, selection.project, selection.scenario, technology?.id, deferredQuery, filter, country, page]);
  useEffect(() => { setPage(0); }, [query, filter, country, technology?.id]);
  useEffect(() => {
    const warn = (event: BeforeUnloadEvent) => { if (changes.size) event.preventDefault(); };
    window.addEventListener("beforeunload", warn); return () => window.removeEventListener("beforeunload", warn);
  }, [changes.size]);
  useEffect(() => { setEditorDirty(editorId, changes.size > 0); return () => setEditorDirty(editorId, false); }, [editorId, changes.size]);

  const filterColumn = table?.filter_column;
  const countryOptions = table?.country_options || [];
  const showCountryFilter = countryOptions.length > 1;
  const showFilter = Boolean(filterColumn && filterColumn.toLowerCase() !== "country" && !technology);
  useEffect(() => { if (country !== "ALL" && !countryOptions.includes(country)) setCountry("ALL"); }, [country, countryOptions]);
  const filterOptions = table?.filter_options || [];
  useEffect(() => { if (filter !== "ALL" && !filterOptions.includes(filter)) setFilter("ALL"); }, [filter, filterOptions]);
  const pageCount = Math.max(1, Math.ceil((table?.total_filtered_rows || 0) / PAGE_SIZE));

  const edit = (rowId: number, column: string, value: InputCell) => {
    const key = `${rowId}:${column}`;
    setRows((current) => current.map((row) => row.__row_id === rowId ? { ...row, [column]: value } : row));
    setChanges((current) => { const next = new Map(current); next.set(key, value); changesRef.current = next; return next; });
    setSuccess("");
  };
  const discard = () => { if (table) setRows(table.rows); const next = new Map<string, InputCell>(); changesRef.current = next; setChanges(next); setSuccess(""); };
  const save = async () => {
    if (!table || !changes.size) return;
    setSaving(true); setError(""); setSuccess("");
    try {
      const payload = [...changes].map(([key, value]) => { const split = key.indexOf(":"); return { row: Number(key.slice(0, split)), column: key.slice(split + 1), value }; });
      const data = await saveInputTable(selection, definition, table.revision, payload, {
        technology,
        country,
        filterValue: filter,
        query: deferredQuery,
        offset: page * PAGE_SIZE,
        limit: PAGE_SIZE,
      });
      const next = new Map<string, InputCell>();
      changesRef.current = next;
      setTable(data); setRows(data.rows); setChanges(next); setSuccess(`Saved ${payload.length} ${payload.length === 1 ? "cell" : "cells"} directly to the CSV.`);
    } catch (reason) { setError(reason instanceof Error ? reason.message : "Could not save changes."); }
    finally { setSaving(false); }
  };

  if (hideWhenEmpty && (loading || (!error && table?.total_filtered_rows === 0))) return null;
  return <section className={`editor-panel${technology ? " technology-panel" : ""}`}>
    <header className="editor-panel-head"><div><p className="eyebrow">{definition.scope === "global" ? "Global source" : `${definition.sector} · ${selection.scenario}`}</p><h2>{definition.label}</h2>{table && <code>{table.path}</code>}</div><div className="editor-actions"><button className="button secondary" disabled={!changes.size || saving} onClick={discard}><RotateCcw aria-hidden="true" />Discard</button><button className="button primary" disabled={!changes.size || saving || definition.timeseries} onClick={save}><Save aria-hidden="true" />{saving ? "Saving…" : `Save changes${changes.size ? ` (${changes.size})` : ""}`}</button></div></header>
    {definition.timeseries && <div className="editor-warning"><AlertTriangle aria-hidden="true" /><span><b>Read-only timeseries.</b> Large hourly inputs are shown for inspection and edited locally outside the browser.</span></div>}
    {error && <div className="notice error">{error}<button onClick={() => void load()}>Reload</button></div>}
    {success && <div className="notice success"><Check aria-hidden="true" />{success}</div>}
    {loading ? <div className="editor-loading"><span className="spinner" />Reading CSV…</div> : table && <>
      {(!technology || showCountryFilter) && <div className="table-tools">{!technology && <label className="search"><Search aria-hidden="true" /><input value={query} onChange={(event) => setQuery(event.target.value)} type="search" placeholder="Find a row or value" /></label>}{showCountryFilter && <label className="field compact"><span>Country</span><select value={country} onChange={(event) => setCountry(event.target.value)}><option value="ALL">All countries</option>{countryOptions.map((value) => <option key={value}>{value}</option>)}</select></label>}{showFilter && filterColumn && <label className="field compact"><span>{filterColumn.replaceAll("_", " ")}</span><select value={filter} onChange={(event) => setFilter(event.target.value)}><option value="ALL">All values</option>{filterOptions.map((value) => <option key={value}>{value}</option>)}</select></label>}{!technology && <span className="row-count">{table.total_filtered_rows.toLocaleString()} of {table.total_rows.toLocaleString()} rows</span>}</div>}
      <div className="editable-table-wrap"><table className="editable-table"><thead><tr>{table.columns.map((column) => <th key={column.name}><span>{column.label}</span>{column.editable && <small>Editable</small>}</th>)}</tr></thead><tbody>{rows.map((row) => <tr key={row.__row_id}>{table.columns.map((column) => <td key={column.name} className={column.editable ? "editable-cell" : "locked-cell"}>{column.editable ? <CellEditor value={row[column.name]} kind={column.kind} onChange={(value) => edit(row.__row_id, column.name, value)} /> : <span>{String(row[column.name] ?? "")}</span>}</td>)}</tr>)}</tbody></table></div>
      <footer className="table-pagination"><span>Page {Math.min(page + 1, pageCount)} of {pageCount} · {table.total_filtered_rows.toLocaleString()} matching rows</span><div><button className="icon-button" aria-label="Previous page" disabled={page === 0} onClick={() => setPage((current) => Math.max(0, current - 1))}><ChevronLeft /></button><button className="icon-button" aria-label="Next page" disabled={page >= pageCount - 1} onClick={() => setPage((current) => Math.min(pageCount - 1, current + 1))}><ChevronRight /></button></div></footer>
    </>}
  </section>;
}

function CellEditor({ value, kind, onChange }: { value: InputCell; kind: string; onChange: (value: InputCell) => void }) {
  if (kind === "boolean") return <label className="cell-check"><input type="checkbox" checked={Boolean(value)} onChange={(event) => onChange(event.target.checked)} /><span aria-hidden="true" /></label>;
  return <input className="cell-input" value={String(value ?? "")} inputMode={kind === "number" ? "decimal" : undefined} onChange={(event) => onChange(event.target.value)} aria-label="Editable cell" />;
}
