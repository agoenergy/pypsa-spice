import type { ResultRow } from "./types";
import { X } from "lucide-react";

interface Props { title: string; rows: ResultRow[]; sourceCount: number; onClose: () => void }

export default function DataDialog({ title, rows, sourceCount, onClose }: Props) {
  const columns = [...new Set(rows.flatMap((row) => Object.keys(row)))];
  return <div className="dialog-backdrop" role="presentation" onMouseDown={(event) => { if (event.target === event.currentTarget) onClose(); }}>
    <section className="dialog" role="dialog" aria-modal="true" aria-labelledby="dialog-title">
      <header><div><p className="eyebrow">Displayed data</p><h2 id="dialog-title">{title}</h2></div><button className="icon-button" onClick={onClose} aria-label="Close data table"><X aria-hidden="true" /></button></header>
      <div className="table-wrap"><table><thead><tr>{columns.map((column) => <th key={column}>{column.replaceAll("_", " ")}</th>)}</tr></thead><tbody>{rows.slice(0, 500).map((row, index) => <tr key={index}>{columns.map((column) => <td className={typeof row[column] === "number" ? "number" : ""} key={column}>{typeof row[column] === "number" ? Number(row[column]).toLocaleString(undefined, { maximumFractionDigits: 3 }) : row[column]}</td>)}</tr>)}</tbody></table></div>
      <footer><span>Showing {Math.min(500, rows.length).toLocaleString()} of {rows.length.toLocaleString()} displayed rows ({sourceCount.toLocaleString()} source rows).</span><button className="button secondary" onClick={onClose}>Close</button></footer>
    </section>
  </div>;
}
