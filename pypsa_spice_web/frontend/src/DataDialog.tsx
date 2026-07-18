import type { ResultRow } from "./types";
import { X } from "lucide-react";

interface Props { title: string; rows: ResultRow[]; sourceCount: number; hourly: boolean; onClose: () => void }

interface DisplayTable {
  columns: string[];
  rows: ResultRow[];
  yearColumns: Set<string>;
  pivoted: boolean;
}

function compareYears(left: string, right: string): number {
  const leftNumber = Number(left);
  const rightNumber = Number(right);
  if (Number.isFinite(leftNumber) && Number.isFinite(rightNumber)) return leftNumber - rightNumber;
  return left.localeCompare(right, undefined, { numeric: true });
}

function flatTable(rows: ResultRow[]): DisplayTable {
  return {
    columns: [...new Set(rows.flatMap((row) => Object.keys(row)))],
    rows,
    yearColumns: new Set(),
    pivoted: false,
  };
}

export function pivotYearlyRows(rows: ResultRow[]): DisplayTable {
  if (!rows.length || rows.some((row) => row.year === null || row.year === undefined)) return flatTable(rows);

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
  return typeof value === "number" ? value.toLocaleString(undefined, { maximumFractionDigits: 3 }) : value;
}

export default function DataDialog({ title, rows, sourceCount, hourly, onClose }: Props) {
  const table = hourly ? flatTable(rows) : pivotYearlyRows(rows);
  const visibleRows = table.rows.slice(0, 500);
  return <div className="dialog-backdrop" role="presentation" onMouseDown={(event) => { if (event.target === event.currentTarget) onClose(); }}>
    <section className="dialog" role="dialog" aria-modal="true" aria-labelledby="dialog-title">
      <header><h2 id="dialog-title">{title}</h2><button className="icon-button" onClick={onClose} aria-label="Close data table"><X aria-hidden="true" /></button></header>
      <div className="table-wrap"><table className={table.pivoted ? "pivot-table" : undefined}><thead><tr>{table.columns.map((column) => <th className={table.yearColumns.has(column) ? "number" : undefined} key={column}>{column.replaceAll("_", " ")}</th>)}</tr></thead><tbody>{visibleRows.map((row, index) => <tr key={index}>{table.columns.map((column) => <td className={typeof row[column] === "number" ? "number" : ""} key={column}>{formatCell(row[column])}</td>)}</tr>)}</tbody></table></div>
      <footer><span>{table.pivoted ? <>Showing {visibleRows.length.toLocaleString()} of {table.rows.length.toLocaleString()} pivoted rows from {rows.length.toLocaleString()} displayed values ({sourceCount.toLocaleString()} source rows).</> : <>Showing {visibleRows.length.toLocaleString()} of {table.rows.length.toLocaleString()} displayed rows ({sourceCount.toLocaleString()} source rows).</>}</span><button className="button secondary" onClick={onClose}>Close</button></footer>
    </section>
  </div>;
}
