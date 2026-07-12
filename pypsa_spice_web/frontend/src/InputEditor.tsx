import { useEffect, useMemo, useState } from "react";
import { AlertTriangle, Check, ChevronLeft, ChevronRight, RotateCcw, Save, Search } from "lucide-react";
import { getInputTable, saveInputTable } from "./api";
import type { InputCatalog, InputCell, InputRow, InputSelection, InputTableDefinition, InputTableResponse } from "./types";

const PAGE_SIZE = 100;

export default function InputEditor({ catalog, selection }: { catalog: InputCatalog; selection: InputSelection }) {
  const [scope, setScope] = useState<"global" | "scenario">("global");
  const [sector, setSector] = useState("power");
  const definitions = scope === "global" ? catalog.global_tables : (catalog.sector_tables[sector] || []);
  const [tableId, setTableId] = useState(catalog.global_tables.find((item) => item.id === "Technologies")?.id || definitions[0]?.id || "");

  useEffect(() => {
    const available = scope === "global" ? catalog.global_tables : (catalog.sector_tables[sector] || []);
    const preferred = scope === "global" ? "Technologies" : ({ power: "Power_generators", industry: "Heat_generators", transport: "Transport_loads" } as Record<string, string>)[sector];
    if (!available.some((item) => item.id === tableId)) setTableId(available.find((item) => item.id === preferred)?.id || available[0]?.id || "");
  }, [catalog, scope, sector, tableId]);

  const definition = definitions.find((item) => item.id === tableId);
  return <>
    <section className="page-title editor-title"><div><p className="eyebrow pink">Model inputs</p><h1>Input data</h1><p>Edit the model’s source CSV files. Changes are written only when you select Save changes.</p></div></section>
    <section className="editor-toolbar" aria-label="Input table selection">
      <div className="segmented" role="tablist" aria-label="Input scope">
        <button className={scope === "global" ? "active" : ""} onClick={() => setScope("global")} role="tab" aria-selected={scope === "global"}>Global inputs</button>
        <button className={scope === "scenario" ? "active" : ""} onClick={() => setScope("scenario")} role="tab" aria-selected={scope === "scenario"}>Scenario inputs</button>
      </div>
      {scope === "scenario" && <div className="segmented" role="tablist" aria-label="Input sector">{["power", "industry", "transport"].map((item) => <button key={item} className={sector === item ? "active" : ""} onClick={() => setSector(item)} role="tab" aria-selected={sector === item}>{item}</button>)}</div>}
      <label className="field compact"><span>Table</span><select value={tableId} onChange={(event) => setTableId(event.target.value)}>{definitions.map((item) => <option value={item.id} key={item.id}>{item.label}</option>)}</select></label>
    </section>
    {definition ? <TableEditor key={`${selection.dataset}:${selection.project}:${selection.scenario}:${scope}:${sector}:${definition.id}`} definition={definition} selection={selection} /> : <div className="editor-empty">No configured table is available for this selection.</div>}
  </>;
}

function TableEditor({ definition, selection }: { definition: InputTableDefinition; selection: InputSelection }) {
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

  const filterColumn = table?.filter_column;
  const filterOptions = useMemo(() => {
    if (!filterColumn) return [];
    return [...new Set(rows.map((row) => String(row[filterColumn] ?? "")).filter(Boolean))].sort();
  }, [rows, filterColumn]);
  const filtered = useMemo(() => rows.filter((row) => {
    if (filterColumn && filter !== "ALL" && String(row[filterColumn]) !== filter) return false;
    const needle = query.trim().toLowerCase();
    return !needle || Object.values(row).some((value) => String(value ?? "").toLowerCase().includes(needle));
  }), [rows, filterColumn, filter, query]);
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

  return <section className="editor-panel">
    <header className="editor-panel-head"><div><p className="eyebrow">{definition.scope === "global" ? "Global source" : `${definition.sector} · ${selection.scenario}`}</p><h2>{definition.label}</h2>{table && <code>{table.path}</code>}</div><div className="editor-actions"><button className="button secondary" disabled={!changes.size || saving} onClick={discard}><RotateCcw aria-hidden="true" />Discard</button><button className="button primary" disabled={!changes.size || saving || definition.timeseries} onClick={save}><Save aria-hidden="true" />{saving ? "Saving…" : `Save changes${changes.size ? ` (${changes.size})` : ""}`}</button></div></header>
    {definition.timeseries && <div className="editor-warning"><AlertTriangle aria-hidden="true" /><span><b>Read-only timeseries.</b> As in the existing app, large hourly inputs are shown for inspection and edited locally outside the browser.</span></div>}
    {error && <div className="notice error">{error}<button onClick={() => void load()}>Reload</button></div>}
    {success && <div className="notice success"><Check aria-hidden="true" />{success}</div>}
    {loading ? <div className="editor-loading"><span className="spinner" />Reading CSV…</div> : table && <>
      <div className="table-tools"><label className="search"><Search aria-hidden="true" /><input value={query} onChange={(event) => setQuery(event.target.value)} type="search" placeholder="Find a row or value" /></label>{filterColumn && <label className="field compact"><span>{filterColumn.replaceAll("_", " ")}</span><select value={filter} onChange={(event) => setFilter(event.target.value)}><option value="ALL">All values</option>{filterOptions.map((value) => <option key={value}>{value}</option>)}</select></label>}<span className="row-count">{filtered.length.toLocaleString()} of {table.total_rows.toLocaleString()} rows</span></div>
      <div className="editable-table-wrap"><table className="editable-table"><thead><tr>{table.columns.map((column) => <th key={column.name}><span>{column.label}</span>{column.editable && <small>Editable</small>}</th>)}</tr></thead><tbody>{visibleRows.map((row) => <tr key={row.__row_id}>{table.columns.map((column) => <td key={column.name} className={column.editable ? "editable-cell" : "locked-cell"}>{column.editable ? <CellEditor value={row[column.name]} kind={column.kind} onChange={(value) => edit(row.__row_id, column.name, value)} /> : <span>{String(row[column.name] ?? "")}</span>}</td>)}</tr>)}</tbody></table></div>
      <footer className="table-pagination"><span>Page {Math.min(page + 1, pageCount)} of {pageCount}{table.truncated ? " · first 10,000 rows loaded" : ""}</span><div><button className="icon-button" aria-label="Previous page" disabled={page === 0} onClick={() => setPage((current) => Math.max(0, current - 1))}><ChevronLeft /></button><button className="icon-button" aria-label="Next page" disabled={page >= pageCount - 1} onClick={() => setPage((current) => Math.min(pageCount - 1, current + 1))}><ChevronRight /></button></div></footer>
    </>}
  </section>;
}

function CellEditor({ value, kind, onChange }: { value: InputCell; kind: string; onChange: (value: InputCell) => void }) {
  if (kind === "boolean") return <label className="cell-check"><input type="checkbox" checked={Boolean(value)} onChange={(event) => onChange(event.target.checked)} /><span aria-hidden="true" /></label>;
  return <input className="cell-input" value={String(value ?? "")} inputMode={kind === "number" ? "decimal" : undefined} onChange={(event) => onChange(event.target.value)} aria-label="Editable cell" />;
}
