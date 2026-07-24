import { useEffect, useMemo, useState } from "react";
import {
  getCoreRowModel,
  getFilteredRowModel,
  getPaginationRowModel,
  getSortedRowModel,
  useReactTable,
  type ColumnDef,
  type ColumnFiltersState,
  type ColumnSizingState,
  type FilterFn,
  type PaginationState,
  type SortingState,
  type VisibilityState,
} from "@tanstack/react-table";
import { ChevronLeft, ChevronRight, Columns, Download, ListFilter, RotateCcw, Search, X } from "lucide-react";
import type { ChartDefinition, ResultRow } from "./types";

interface Props {
  chart: ChartDefinition;
  rows: ResultRow[];
  sourceCount: number;
  onClose: () => void;
}

export interface DisplayTable {
  columns: string[];
  rows: ResultRow[];
  yearColumns: Set<string>;
  pivoted: boolean;
}

const PAGE_SIZES = [50, 100, 250] as const;

const categoricalFilter: FilterFn<ResultRow> = (row, columnId, selected: string[]) => (
  !selected?.length || selected.includes(String(row.getValue(columnId) ?? ""))
);
categoricalFilter.autoRemove = (selected: string[]) => !selected?.length;

function compareYears(left: string, right: string): number {
  const leftNumber = Number(left);
  const rightNumber = Number(right);
  if (Number.isFinite(leftNumber) && Number.isFinite(rightNumber)) return leftNumber - rightNumber;
  return left.localeCompare(right, undefined, { numeric: true });
}

function columnsWithScenarioFirst(rows: ResultRow[]): string[] {
  const columns = [...new Set(rows.flatMap((row) => Object.keys(row)))];
  const scenarioIndex = columns.indexOf("scenario");
  if (scenarioIndex > 0) columns.unshift(...columns.splice(scenarioIndex, 1));
  return columns;
}

function flatTable(rows: ResultRow[]): DisplayTable {
  return {
    columns: [...new Set(rows.flatMap((row) => Object.keys(row)))],
    rows,
    yearColumns: new Set(),
    pivoted: false,
  };
}

function yearlyFlatTable(rows: ResultRow[]): DisplayTable {
  return {
    columns: columnsWithScenarioFirst(rows),
    rows,
    yearColumns: new Set(),
    pivoted: false,
  };
}

export function pivotYearlyRows(rows: ResultRow[]): DisplayTable {
  if (!rows.length || rows.some((row) => row.year === null || row.year === undefined)) return yearlyFlatTable(rows);

  const allColumns = [...new Set(rows.flatMap((row) => Object.keys(row)))];
  const dimensionColumns = allColumns.filter((column) => column !== "year" && column !== "value");
  const scenarioIndex = dimensionColumns.indexOf("scenario");
  if (scenarioIndex > 0) dimensionColumns.unshift(...dimensionColumns.splice(scenarioIndex, 1));

  const years = [...new Set(rows.map((row) => String(row.year)))].sort(compareYears);
  const pivoted = new Map<string, ResultRow>();

  for (const sourceRow of rows) {
    const dimensions = dimensionColumns.map((column) => sourceRow[column] ?? null);
    const key = JSON.stringify(dimensions);
    let targetRow = pivoted.get(key);
    if (!targetRow) {
      targetRow = Object.fromEntries(dimensionColumns.map((column, index) => [column, dimensions[index]]));
      pivoted.set(key, targetRow);
    }

    const year = String(sourceRow.year);
    const value = sourceRow.value;
    if (typeof value === "number" && typeof targetRow[year] === "number") targetRow[year] += value;
    else targetRow[year] = value;
  }

  return {
    columns: [...dimensionColumns, ...years],
    rows: [...pivoted.values()],
    yearColumns: new Set(years),
    pivoted: true,
  };
}

function formatCell(value: ResultRow[string]): string | number | null {
  return typeof value === "number"
    ? value.toLocaleString(undefined, { maximumFractionDigits: 3 })
    : value;
}

function displayLabel(column: string): string {
  return column.replaceAll("_", " ");
}

function isNumericColumn(table: DisplayTable, column: string): boolean {
  return table.rows.some((row) => typeof row[column] === "number");
}

export function summarizeRows(rows: ResultRow[], columns: string[]): Record<string, number> {
  return Object.fromEntries(columns.map((column) => [
    column,
    rows.reduce((sum, row) => sum + (typeof row[column] === "number" ? row[column] : 0), 0),
  ]));
}

function csvCell(value: ResultRow[string]): string {
  const text = value === null || value === undefined ? "" : String(value);
  return /[",\r\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

export function tableToCsv(rows: ResultRow[], columns: string[]): string {
  return [
    columns.map(csvCell).join(","),
    ...rows.map((row) => columns.map((column) => csvCell(row[column])).join(",")),
  ].join("\r\n");
}

export function limitHourlyRows(rows: ResultRow[]): ResultRow[] {
  return rows.slice(0, 500);
}

function HourlyDataDialog({ title, rows, sourceCount, onClose }: Omit<Props, "chart"> & { title: string }) {
  const table = flatTable(rows);
  const visibleRows = limitHourlyRows(table.rows);
  return <div className="dialog-backdrop" role="presentation" onMouseDown={(event) => { if (event.target === event.currentTarget) onClose(); }}>
    <section className="dialog" role="dialog" aria-modal="true" aria-labelledby="dialog-title">
      <header><h2 id="dialog-title">{title}</h2><button className="icon-button" onClick={onClose} aria-label="Close data table"><X aria-hidden="true" /></button></header>
      <div className="table-wrap"><table><thead><tr>{table.columns.map((column) => <th key={column}>{displayLabel(column)}</th>)}</tr></thead><tbody>{visibleRows.map((row, index) => <tr key={index}>{table.columns.map((column) => <td className={typeof row[column] === "number" ? "number" : ""} key={column}>{formatCell(row[column])}</td>)}</tr>)}</tbody></table></div>
      <footer><span>Showing {visibleRows.length.toLocaleString()} of {table.rows.length.toLocaleString()} displayed rows ({sourceCount.toLocaleString()} source rows).</span><button className="button secondary" onClick={onClose}>Close</button></footer>
    </section>
  </div>;
}

export default function DataDialog({ chart, rows, sourceCount, onClose }: Props) {
  if (chart.hourly) {
    return <HourlyDataDialog title={chart.name} rows={rows} sourceCount={sourceCount} onClose={onClose} />;
  }
  return <YearlyDataDialog chart={chart} rows={rows} sourceCount={sourceCount} onClose={onClose} />;
}

function YearlyDataDialog({ chart, rows, sourceCount, onClose }: Props) {
  const canPivot = rows.length > 0 && rows.every((row) => row.year !== null && row.year !== undefined);
  const isDifferenceView = rows.some(
    (row) => typeof row.scenario === "string" && row.scenario.includes(" − "),
  );
  const [view, setView] = useState<"pivot" | "raw">(canPivot ? "pivot" : "raw");
  const [globalFilter, setGlobalFilter] = useState("");
  const [columnFilters, setColumnFilters] = useState<ColumnFiltersState>([]);
  const [sorting, setSorting] = useState<SortingState>([]);
  const [columnVisibility, setColumnVisibility] = useState<VisibilityState>({});
  const [columnSizing, setColumnSizing] = useState<ColumnSizingState>({});
  const [pagination, setPagination] = useState<PaginationState>({ pageIndex: 0, pageSize: 50 });

  const pivotedTable = useMemo(() => pivotYearlyRows(rows), [rows]);
  const rawTable = useMemo(() => yearlyFlatTable(rows), [rows]);
  const displayTable = view === "pivot" ? pivotedTable : rawTable;
  const firstColumn = displayTable.columns[0];
  const numericColumns = useMemo(
    () => displayTable.columns.filter((column) => isNumericColumn(displayTable, column)),
    [displayTable],
  );
  const columnDefinitions = useMemo<ColumnDef<ResultRow>[]>(() => (
    displayTable.columns.map((column) => {
      const numeric = numericColumns.includes(column);
      return {
        id: column,
        accessorFn: (row) => row[column],
        enableHiding: column !== firstColumn,
        filterFn: numeric ? "inNumberRange" : categoricalFilter,
        sortingFn: numeric ? "basic" : "alphanumeric",
        sortDescFirst: false,
        sortUndefined: "last",
        minSize: 84,
        maxSize: 420,
        size: numeric ? 116 : column === "scenario" ? (isDifferenceView ? 280 : 170) : 140,
      };
    })
  ), [displayTable.columns, firstColumn, isDifferenceView, numericColumns]);

  const table = useReactTable({
    data: displayTable.rows,
    columns: columnDefinitions,
    state: {
      globalFilter,
      columnFilters,
      sorting,
      columnVisibility,
      columnSizing,
      pagination,
      columnPinning: { left: firstColumn ? [firstColumn] : [], right: [] },
    },
    onGlobalFilterChange: setGlobalFilter,
    onColumnFiltersChange: setColumnFilters,
    onSortingChange: setSorting,
    onColumnVisibilityChange: setColumnVisibility,
    onColumnSizingChange: setColumnSizing,
    onPaginationChange: setPagination,
    globalFilterFn: "includesString",
    enableMultiSort: false,
    columnResizeMode: "onChange",
    getCoreRowModel: getCoreRowModel(),
    getFilteredRowModel: getFilteredRowModel(),
    getSortedRowModel: getSortedRowModel(),
    getPaginationRowModel: getPaginationRowModel(),
  });

  const filteredRows = table.getFilteredRowModel().rows.map((row) => row.original);
  const sortedRows = table.getSortedRowModel().rows.map((row) => row.original);
  const pageRows = table.getRowModel().rows;
  const visibleColumns = table.getVisibleLeafColumns();
  const summaryColumns = displayTable.pivoted
    ? [...displayTable.yearColumns]
    : numericColumns.filter((column) => column === "value");
  const totals = useMemo(
    () => chart.summary === "sum" ? summarizeRows(filteredRows, summaryColumns) : {},
    [chart.summary, filteredRows, summaryColumns.join("\u0000")],
  );
  const activeFilterCount = columnFilters.length + (globalFilter.trim() ? 1 : 0);

  useEffect(() => {
    setColumnFilters([]);
    setSorting([]);
    setPagination((current) => ({ ...current, pageIndex: 0 }));
  }, [view]);
  useEffect(() => {
    const openTableMenus = () => (
      [...document.querySelectorAll<HTMLDetailsElement>(".yearly-table-dialog details[open]")]
    );
    const closeMenusOutside = (event: PointerEvent) => {
      if (!(event.target instanceof Node)) return;
      for (const menu of openTableMenus()) {
        if (!menu.contains(event.target)) menu.removeAttribute("open");
      }
    };
    const closeWithEscape = (event: KeyboardEvent) => {
      if (event.key !== "Escape") return;
      const menus = openTableMenus();
      if (menus.length > 0) {
        event.preventDefault();
        for (const menu of menus) menu.removeAttribute("open");
        return;
      }
      onClose();
    };
    document.addEventListener("pointerdown", closeMenusOutside);
    document.addEventListener("keydown", closeWithEscape);
    return () => {
      document.removeEventListener("pointerdown", closeMenusOutside);
      document.removeEventListener("keydown", closeWithEscape);
    };
  }, [onClose]);

  const categoricalOptions = (column: string) => [...new Set(
    displayTable.rows.map((row) => String(row[column] ?? "")).filter(Boolean),
  )].sort((left, right) => left.localeCompare(right, undefined, { numeric: true }));

  const updateCategoricalFilter = (column: string, value: string, checked: boolean) => {
    const tableColumn = table.getColumn(column);
    const selected = new Set((tableColumn?.getFilterValue() as string[] | undefined) || []);
    if (checked) selected.add(value);
    else selected.delete(value);
    tableColumn?.setFilterValue([...selected]);
  };

  const updateNumericFilter = (column: string, index: 0 | 1, value: string) => {
    const tableColumn = table.getColumn(column);
    const range = [...((tableColumn?.getFilterValue() as [number | undefined, number | undefined] | undefined) || [undefined, undefined])] as [number | undefined, number | undefined];
    range[index] = value === "" ? undefined : Number(value);
    tableColumn?.setFilterValue(range);
  };

  const clearAllFilters = () => {
    setGlobalFilter("");
    setColumnFilters([]);
  };

  const exportFiltered = () => {
    const columns = visibleColumns.map((column) => column.id);
    const csv = tableToCsv(sortedRows, columns);
    const url = URL.createObjectURL(new Blob([csv], { type: "text/csv;charset=utf-8" }));
    const anchor = document.createElement("a");
    anchor.href = url;
    anchor.download = `${chart.name.toLocaleLowerCase().replace(/[^a-z0-9]+/g, "_").replace(/^_|_$/g, "") || "results"}_${view}_filtered.csv`;
    anchor.click();
    URL.revokeObjectURL(url);
  };

  return <div className="dialog-backdrop yearly-table-backdrop" role="presentation" onMouseDown={(event) => { if (event.target === event.currentTarget) onClose(); }}>
    <section className="dialog yearly-table-dialog" role="dialog" aria-modal="true" aria-labelledby="dialog-title">
      <header className="yearly-table-header">
        <div>
          <p className="eyebrow">Results table · {chart.units || "Values"}</p>
          <h2 id="dialog-title">{chart.name}</h2>
        </div>
        <button className="icon-button" onClick={onClose} aria-label="Close data table"><X aria-hidden="true" /></button>
      </header>

      <div className="yearly-table-toolbar">
        <label className="yearly-table-search">
          <Search aria-hidden="true" />
          <span className="sr-only">Search table</span>
          <input value={globalFilter} onChange={(event) => setGlobalFilter(event.target.value)} placeholder="Search all columns…" />
          {globalFilter && <button type="button" onClick={() => setGlobalFilter("")} aria-label="Clear search"><X aria-hidden="true" /></button>}
        </label>
        <div className="yearly-view-toggle" aria-label="Table view">
          <button className={view === "pivot" ? "active" : ""} type="button" onClick={() => setView("pivot")} disabled={!pivotedTable.pivoted}>Pivoted</button>
          <button className={view === "raw" ? "active" : ""} type="button" onClick={() => setView("raw")}>Raw</button>
        </div>
        <details className="yearly-columns-menu">
          <summary><Columns aria-hidden="true" />Columns</summary>
          <div>
            <b>Visible columns</b>
            {table.getAllLeafColumns().map((column) => <label key={column.id}>
              <input
                type="checkbox"
                checked={column.getIsVisible()}
                disabled={!column.getCanHide()}
                onChange={(event) => column.toggleVisibility(event.target.checked)}
              />
              <span>{displayLabel(column.id)}</span>
            </label>)}
          </div>
        </details>
        <button className="button secondary yearly-export" type="button" onClick={exportFiltered}><Download aria-hidden="true" />Export filtered view</button>
        {activeFilterCount > 0 && <button className="yearly-reset" type="button" onClick={clearAllFilters}><RotateCcw aria-hidden="true" />Clear filters</button>}
      </div>

      {activeFilterCount > 0 && <div className="yearly-filter-chips" aria-label="Active filters">
        {globalFilter.trim() && <button type="button" onClick={() => setGlobalFilter("")}>Search: “{globalFilter.trim()}” <X aria-hidden="true" /></button>}
        {columnFilters.map((filter) => {
          const numeric = numericColumns.includes(filter.id);
          const value = filter.value as string[] | [number | undefined, number | undefined];
          const label = numeric
            ? `${value[0] !== undefined ? `≥ ${value[0]}` : ""}${value[0] !== undefined && value[1] !== undefined ? " · " : ""}${value[1] !== undefined ? `≤ ${value[1]}` : ""}`
            : `${value.length} selected`;
          return <button type="button" key={filter.id} onClick={() => table.getColumn(filter.id)?.setFilterValue(undefined)}>{displayLabel(filter.id)}: {label} <X aria-hidden="true" /></button>;
        })}
      </div>}

      <div className="table-wrap yearly-table-wrap">
        <table
          className={`yearly-results-table ${displayTable.pivoted ? "pivot-table" : ""}`}
          style={{ width: table.getTotalSize() }}
        >
          <thead>{table.getHeaderGroups().map((headerGroup) => <tr key={headerGroup.id}>{headerGroup.headers.map((header) => {
            const column = header.column;
            const numeric = numericColumns.includes(column.id);
            const selectedValues = (column.getFilterValue() as string[] | undefined) || [];
            const range = (column.getFilterValue() as [number | undefined, number | undefined] | undefined) || [undefined, undefined];
            const filtered = column.getIsFiltered();
            const pinned = column.getIsPinned();
            const sorted = column.getIsSorted();
            return <th
              className={`${numeric ? "number" : ""} ${pinned ? "pinned-column" : ""}`}
              key={header.id}
              style={{
                width: header.getSize(),
                minWidth: header.getSize(),
                maxWidth: header.getSize(),
                left: pinned === "left" ? column.getStart("left") : undefined,
              }}
            >
              <div className="yearly-column-heading">
                <button className="yearly-sort" type="button" onClick={column.getToggleSortingHandler()} aria-label={`Sort by ${displayLabel(column.id)}`}>
                  <span>{displayLabel(column.id)}{(displayTable.yearColumns.has(column.id) || column.id === "value") && chart.units ? <small>{chart.units}</small> : null}</span>
                  <i aria-hidden="true">{sorted === "asc" ? "↑" : sorted === "desc" ? "↓" : "↕"}</i>
                </button>
                <details className={`yearly-column-filter ${filtered ? "active" : ""}`}>
                  <summary aria-label={`Filter ${displayLabel(column.id)}`} title={`Filter ${displayLabel(column.id)}`}>
                    <ListFilter aria-hidden="true" />
                  </summary>
                  <div onClick={(event) => event.stopPropagation()}>
                    <header><b>{displayLabel(column.id)}</b>{filtered && <button type="button" onClick={() => column.setFilterValue(undefined)}>Clear</button>}</header>
                    {numeric ? <div className="yearly-number-filter">
                      <label><span>Minimum</span><input type="number" value={range[0] ?? ""} onChange={(event) => updateNumericFilter(column.id, 0, event.target.value)} /></label>
                      <label><span>Maximum</span><input type="number" value={range[1] ?? ""} onChange={(event) => updateNumericFilter(column.id, 1, event.target.value)} /></label>
                    </div> : <div className="yearly-value-filter">
                      {categoricalOptions(column.id).map((value) => <label key={value}>
                        <input type="checkbox" checked={selectedValues.includes(value)} onChange={(event) => updateCategoricalFilter(column.id, value, event.target.checked)} />
                        <span>{value}</span>
                      </label>)}
                    </div>}
                  </div>
                </details>
              </div>
              <span
                className={`yearly-column-resizer ${column.getIsResizing() ? "is-resizing" : ""}`}
                role="separator"
                aria-orientation="vertical"
                onMouseDown={header.getResizeHandler()}
                onTouchStart={header.getResizeHandler()}
              />
            </th>;
          })}</tr>)}</thead>
          <tbody>{pageRows.map((row) => <tr key={row.id}>{row.getVisibleCells().map((cell) => {
            const pinned = cell.column.getIsPinned();
            const value = cell.getValue<ResultRow[string]>();
            return <td
              className={`${typeof value === "number" ? "number" : ""} ${pinned ? "pinned-column" : ""}`}
              key={cell.id}
              style={{
                width: cell.column.getSize(),
                minWidth: cell.column.getSize(),
                maxWidth: cell.column.getSize(),
                left: pinned === "left" ? cell.column.getStart("left") : undefined,
              }}
              title={value === null || value === undefined ? undefined : String(value)}
            >{formatCell(value)}</td>;
          })}</tr>)}</tbody>
          {chart.summary === "sum" && filteredRows.length > 0 && <tfoot><tr>{visibleColumns.map((column, index) => {
            const pinned = column.getIsPinned();
            return <td
              className={`${summaryColumns.includes(column.id) ? "number" : ""} ${pinned ? "pinned-column" : ""}`}
              key={column.id}
              style={{
                width: column.getSize(),
                minWidth: column.getSize(),
                maxWidth: column.getSize(),
                left: pinned === "left" ? column.getStart("left") : undefined,
              }}
            >{index === 0 ? <b>Total</b> : summaryColumns.includes(column.id) ? <b>{formatCell(totals[column.id])}</b> : null}</td>;
          })}</tr></tfoot>}
        </table>
        {pageRows.length === 0 && <div className="yearly-table-empty"><b>No matching rows</b><span>Change or clear the active filters.</span></div>}
      </div>

      <footer className="yearly-table-footer">
        <span>Showing {pageRows.length.toLocaleString()} of {filteredRows.length.toLocaleString()} filtered rows · {displayTable.rows.length.toLocaleString()} total{displayTable.pivoted ? ` from ${rows.length.toLocaleString()} values` : ""} · {sourceCount.toLocaleString()} source rows</span>
        <label>Rows per page <select value={pagination.pageSize} onChange={(event) => table.setPageSize(Number(event.target.value))}>{PAGE_SIZES.map((size) => <option value={size} key={size}>{size}</option>)}</select></label>
        <div className="yearly-pagination">
          <button type="button" onClick={() => table.previousPage()} disabled={!table.getCanPreviousPage()} aria-label="Previous page"><ChevronLeft aria-hidden="true" /></button>
          <span>Page {pagination.pageIndex + 1} of {Math.max(1, table.getPageCount())}</span>
          <button type="button" onClick={() => table.nextPage()} disabled={!table.getCanNextPage()} aria-label="Next page"><ChevronRight aria-hidden="true" /></button>
        </div>
        <button className="button secondary" onClick={onClose}>Close</button>
      </footer>
    </section>
  </div>;
}
