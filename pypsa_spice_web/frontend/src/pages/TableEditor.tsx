import { useDeferredValue, useEffect, useId, useRef, useState } from "react";
import { AlertTriangle, Check, ChevronLeft, ChevronRight } from "lucide-react";
import { getInputTable, saveInputTable } from "../api";
import { SearchField, SelectField } from "../components/FormControls";
import PageHeader from "../components/PageHeader";
import SaveDiscardActions from "../components/SaveDiscardActions";
import type {
  InputCell,
  InputRow,
  InputSelection,
  InputTableDefinition,
  InputTableResponse,
  InputTechnology,
} from "../types";
import { setEditorDirty } from "../utility";

const PAGE_SIZE = 100;

export function TableTitle({ definition }: { definition: InputTableDefinition }) {
  return <PageHeader title={definition.label} className="selection-title">
    <dl>
      <div><dt>Scope</dt><dd>{definition.scope === "global" ? "Global input" : "Scenario input"}</dd></div>
      <div><dt>Sector</dt><dd>{definition.sector || "All sectors"}</dd></div>
    </dl>
  </PageHeader>;
}

export function TableView({ definition, selection }: {
  definition: InputTableDefinition;
  selection: InputSelection;
}) {
  const global = definition.scope === "global";
  return <div className="table-view">
    <section className="technology-group">
      <header>
        <h2>{global ? "Global input" : "Scenario input"}</h2>
        <span>{global
          ? "Changes here apply to every country and every scenario in this project."
          : "Assets and constraints for this scenario. Country filters appear only on tables with country-specific rows."}</span>
      </header>
      <TableEditor definition={definition} selection={selection} />
    </section>
  </div>;
}

export default function TableEditor({ definition, selection, technology, hideWhenEmpty = false }: {
  definition: InputTableDefinition;
  selection: InputSelection;
  technology?: InputTechnology;
  hideWhenEmpty?: boolean;
}) {
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
    setLoading(true);
    setError("");
    setSuccess("");
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
      setTable(data);
      setRows(loadedRows);
    } catch (reason) {
      if (!(reason instanceof DOMException && reason.name === "AbortError")) {
        setError(reason instanceof Error ? reason.message : "Could not load this table.");
      }
    } finally {
      if (!signal?.aborted) setLoading(false);
    }
  };

  useEffect(() => {
    const controller = new AbortController();
    void load(controller.signal);
    return () => controller.abort();
  }, [definition.id, selection.dataset, selection.project, selection.scenario, technology?.id, deferredQuery, filter, country, page]);
  useEffect(() => { setPage(0); }, [query, filter, country, technology?.id]);
  useEffect(() => {
    const warn = (event: BeforeUnloadEvent) => { if (changes.size) event.preventDefault(); };
    window.addEventListener("beforeunload", warn);
    return () => window.removeEventListener("beforeunload", warn);
  }, [changes.size]);
  useEffect(() => {
    setEditorDirty(editorId, changes.size > 0);
    return () => setEditorDirty(editorId, false);
  }, [editorId, changes.size]);

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
    setChanges((current) => {
      const next = new Map(current);
      next.set(key, value);
      changesRef.current = next;
      return next;
    });
    setSuccess("");
  };
  const discard = () => {
    if (table) setRows(table.rows);
    const next = new Map<string, InputCell>();
    changesRef.current = next;
    setChanges(next);
    setSuccess("");
  };
  const save = async () => {
    if (!table || !changes.size) return;
    setSaving(true);
    setError("");
    setSuccess("");
    try {
      const payload = [...changes].map(([key, value]) => {
        const split = key.indexOf(":");
        return { row: Number(key.slice(0, split)), column: key.slice(split + 1), value };
      });
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
      setTable(data);
      setRows(data.rows);
      setChanges(next);
      setSuccess(`Saved ${payload.length} ${payload.length === 1 ? "cell" : "cells"} directly to the CSV.`);
    } catch (reason) {
      setError(reason instanceof Error ? reason.message : "Could not save changes.");
    } finally {
      setSaving(false);
    }
  };

  if (hideWhenEmpty && (loading || (!error && table?.total_filtered_rows === 0))) return null;
  return <section className={`editor-panel${technology ? " technology-panel" : ""}`}>
    <header className="editor-panel-head">
      <div>
        <p className="eyebrow">{definition.scope === "global" ? "Global source" : `${definition.sector} · ${selection.scenario}`}</p>
        <h2>{definition.label}</h2>
        {table && <code>{table.path}</code>}
      </div>
      <SaveDiscardActions
        hasChanges={changes.size > 0}
        saving={saving}
        saveDisabled={definition.timeseries}
        saveLabel={`Save changes${changes.size ? ` (${changes.size})` : ""}`}
        onDiscard={discard}
        onSave={() => void save()}
      />
    </header>
    {definition.timeseries && <div className="editor-warning">
      <AlertTriangle aria-hidden="true" />
      <span><b>Read-only timeseries.</b> Large hourly inputs are shown for inspection and edited locally outside the browser.</span>
    </div>}
    {error && <div className="notice error">{error}<button onClick={() => void load()}>Reload</button></div>}
    {success && <div className="notice success"><Check aria-hidden="true" />{success}</div>}
    {loading ? <div className="editor-loading"><span className="spinner" />Reading CSV…</div> : table && <>
      {(!technology || showCountryFilter) && <div className="table-tools">
        {!technology && <SearchField value={query} onChange={setQuery} placeholder="Find a row or value" />}
        {showCountryFilter && <SelectField
          label="Country"
          value={country}
          onChange={setCountry}
          compact
          options={[{ value: "ALL", label: "All countries" }, ...countryOptions.map((value) => ({ value, label: value }))]}
        />}
        {showFilter && filterColumn && <SelectField
          label={filterColumn.replaceAll("_", " ")}
          value={filter}
          onChange={setFilter}
          compact
          options={[{ value: "ALL", label: "All values" }, ...filterOptions.map((value) => ({ value, label: value }))]}
        />}
        {!technology && <span className="row-count">{table.total_filtered_rows.toLocaleString()} of {table.total_rows.toLocaleString()} rows</span>}
      </div>}
      <div className="editable-table-wrap"><table className="editable-table"><thead><tr>
        {table.columns.map((column) => <th key={column.name}><span>{column.label}</span>{column.editable && <small>Editable</small>}</th>)}
      </tr></thead><tbody>
        {rows.map((row) => <tr key={row.__row_id}>{table.columns.map((column) => <td key={column.name} className={column.editable ? "editable-cell" : "locked-cell"}>
          {column.editable
            ? <CellEditor value={row[column.name]} kind={column.kind} onChange={(value) => edit(row.__row_id, column.name, value)} />
            : <span>{String(row[column.name] ?? "")}</span>}
        </td>)}</tr>)}
      </tbody></table></div>
      <footer className="table-pagination">
        <span>Page {Math.min(page + 1, pageCount)} of {pageCount} · {table.total_filtered_rows.toLocaleString()} matching rows</span>
        <div>
          <button className="icon-button" aria-label="Previous page" disabled={page === 0} onClick={() => setPage((current) => Math.max(0, current - 1))}><ChevronLeft /></button>
          <button className="icon-button" aria-label="Next page" disabled={page >= pageCount - 1} onClick={() => setPage((current) => Math.min(pageCount - 1, current + 1))}><ChevronRight /></button>
        </div>
      </footer>
    </>}
  </section>;
}

function CellEditor({ value, kind, onChange }: {
  value: InputCell;
  kind: string;
  onChange: (value: InputCell) => void;
}) {
  if (kind === "boolean") {
    return <label className="cell-check">
      <input type="checkbox" checked={Boolean(value)} onChange={(event) => onChange(event.target.checked)} />
      <span aria-hidden="true" />
    </label>;
  }
  return <input
    className="cell-input"
    value={String(value ?? "")}
    inputMode={kind === "number" ? "decimal" : undefined}
    onChange={(event) => onChange(event.target.value)}
    aria-label="Editable cell"
  />;
}
