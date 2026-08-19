import { useEffect, useMemo, useRef, useState } from "react";
import { ArrowRight, CheckCircle2, Info, List, Search, X } from "lucide-react";
import { getScenarioComparison } from "../api";
import PageHeader from "../components/PageHeader";
import "./ScenarioComparison.css";
import type {
  InputSelection,
  ScenarioComparisonResponse,
  ScenarioDifference,
  ScenarioDifferenceSection,
  ScenarioDifferenceStatus,
} from "../types";

export interface ScenarioComparisonFilters {
  query: string;
}

const categoryLabels: Record<string, string> = {
  configuration: "Configuration",
  power: "Power inputs",
  industry: "Industry inputs",
  transport: "Transport inputs",
};
const categoryOrder = ["configuration", "power", "industry", "transport"];

function searchableValue(value: unknown): string {
  if (value === null || value === undefined) return "";
  if (Array.isArray(value)) return value.join(" ");
  if (typeof value === "object") return JSON.stringify(value);
  return String(value);
}

export function filterComparisonSections(
  sections: ScenarioDifferenceSection[],
  filters: ScenarioComparisonFilters,
): ScenarioDifferenceSection[] {
  const needle = filters.query.trim().toLowerCase();
  return sections.flatMap((section) => {
    const changes = section.changes.filter((change) => {
      if (!needle) return true;
      return [
        section.label,
        change.item,
        change.country,
        change.parameter,
        searchableValue(change.reference),
        searchableValue(change.comparison),
      ].some((value) => value.toLowerCase().includes(needle));
    });
    return changes.length ? [{ ...section, changes }] : [];
  });
}

function formatValue(value: unknown): string {
  if (value === null || value === undefined) return "—";
  if (value === "") return "Empty";
  if (value === true) return "True";
  if (value === false) return "False";
  if (typeof value === "number") return value.toLocaleString();
  if (Array.isArray(value)) return value.length ? value.map(formatValue).join(", ") : "None";
  if (typeof value === "object") return JSON.stringify(value);
  if (String(value).toLowerCase() === "inf") return "∞";
  if (String(value).toLowerCase() === "-inf") return "−∞";
  return String(value);
}

function formatDelta(change: ScenarioDifference): string {
  if (change.delta === null) return "—";
  if (change.delta === 0) return "0";
  return `${change.delta > 0 ? "+" : "−"}${Math.abs(change.delta).toLocaleString()}`;
}

function statusLabel(status: ScenarioDifferenceStatus): string {
  if (status === "ambiguous") return "Review";
  return `${status[0].toUpperCase()}${status.slice(1)}`;
}

export default function ScenarioComparison({
  selection,
  comparison,
}: {
  selection: InputSelection;
  comparison: string;
}) {
  const [data, setData] = useState<ScenarioComparisonResponse | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState("");
  const [filters, setFilters] = useState<ScenarioComparisonFilters>({
    query: "",
  });

  const load = async (signal?: AbortSignal) => {
    setLoading(true); setError("");
    try {
      setData(await getScenarioComparison(selection, comparison, signal));
    } catch (reason) {
      if (!(reason instanceof DOMException && reason.name === "AbortError")) {
        setError(reason instanceof Error ? reason.message : "Could not compare the selected scenarios.");
      }
    } finally {
      if (!signal?.aborted) setLoading(false);
    }
  };

  useEffect(() => {
    const controller = new AbortController();
    void load(controller.signal);
    return () => controller.abort();
  }, [selection.dataset, selection.project, selection.scenario, comparison]);

  const filteredSections = useMemo(
    () => filterComparisonSections(data?.sections || [], filters),
    [data, filters],
  );
  const filteredChanges = filteredSections.reduce((total, section) => total + section.changes.length, 0);
  const groupedSections = categoryOrder.flatMap((category) => {
    const sections = filteredSections.filter((section) => section.category === category);
    return sections.length ? [{ category, sections }] : [];
  });

  return <>
    <PageHeader title="Scenario differences" className="comparison-page-title">
      <div className="comparison-direction" aria-label={`${selection.scenario} compared with ${comparison}`}>
        <span><small>Reference</small><b>{selection.scenario}</b></span>
        <ArrowRight aria-hidden="true" />
        <span><small>Comparison</small><b>{comparison}</b></span>
      </div>
    </PageHeader>

    {error && <div className="notice error comparison-notice">{error}<button onClick={() => void load()}>Reload</button></div>}
    {loading ? <div className="editor-loading comparison-loading"><span className="spinner" />Comparing scenario inputs…</div> : data && <>
      <p className="comparison-scope-note"><Info aria-hidden="true" />Only scenario-specific differences are shown. Global inputs are shared by both scenarios and excluded.</p>

      {data.summary.changes > 0 && <section className="comparison-tools" aria-label="Filter scenario differences">
        <label className="comparison-search"><span className="sr-only">Search differences</span><Search aria-hidden="true" /><input type="search" value={filters.query} placeholder="Find an item, parameter, or value" onChange={(event) => setFilters((current) => ({ ...current, query: event.target.value }))} /></label>
        <span className="comparison-filter-count">Showing {filteredChanges.toLocaleString()} of {data.summary.changes.toLocaleString()}</span>
      </section>}

      {data.summary.changes === 0 ? <div className="comparison-empty"><CheckCircle2 aria-hidden="true" /><b>No differences found</b><span>The scenario configuration and scenario-scoped input tables are equivalent.</span></div> : groupedSections.length === 0 ? <div className="comparison-empty"><Search aria-hidden="true" /><b>No matching differences</b><span>Try clearing or changing the filters.</span></div> : <div className="comparison-results">
        {groupedSections.map(({ category, sections }) => <section className="comparison-category" key={category}>
          <header><h2>{categoryLabels[category]}</h2></header>
          <div className="comparison-section-list">{sections.map((section) => <DifferenceTable key={section.id} section={section} reference={data.reference} comparison={data.comparison} />)}</div>
        </section>)}
      </div>}
      {groupedSections.length > 0 && <ComparisonToc sections={filteredSections} />}
    </>}
  </>;
}

function ComparisonToc({ sections }: { sections: ScenarioDifferenceSection[] }) {
  const [open, setOpen] = useState(false);
  const root = useRef<HTMLDivElement>(null);
  useEffect(() => {
    if (!open) return;
    const closeOutside = (event: PointerEvent) => {
      if (!root.current?.contains(event.target as Node)) setOpen(false);
    };
    const closeWithEscape = (event: KeyboardEvent) => {
      if (event.key === "Escape") setOpen(false);
    };
    document.addEventListener("pointerdown", closeOutside);
    document.addEventListener("keydown", closeWithEscape);
    return () => {
      document.removeEventListener("pointerdown", closeOutside);
      document.removeEventListener("keydown", closeWithEscape);
    };
  }, [open]);

  return <div className={`results-toc ${open ? "open" : ""}`} ref={root}>
    {open && <nav className="results-toc-panel" id="comparison-section-list" aria-label="Changed groups on this page">
      <header><h2>Changed groups</h2><button className="icon-button" onClick={() => setOpen(false)} aria-label="Close changed-group list"><X aria-hidden="true" /></button></header>
      <ol>{sections.map((section, index) => <li key={section.id}><a href={`#${section.id.replaceAll(":", "-")}`} onClick={() => setOpen(false)}><span>{String(index + 1).padStart(2, "0")}</span><b>{section.label}</b></a></li>)}</ol>
    </nav>}
    <button className="results-toc-trigger" onClick={() => setOpen((current) => !current)} aria-label="Open changed-group list" aria-expanded={open} aria-controls="comparison-section-list"><List aria-hidden="true" /></button>
  </div>;
}

function DifferenceTable({
  section,
  reference,
  comparison,
}: {
  section: ScenarioDifferenceSection;
  reference: string;
  comparison: string;
}) {
  return <section className="editor-panel comparison-panel" id={section.id.replaceAll(":", "-")}>
    <header className="editor-panel-head"><div><p className="eyebrow">{section.kind === "input" ? "Scenario input" : section.kind === "constraint" ? "Constraint" : "Configuration"}</p><h3>{section.label}</h3></div></header>
    <div className="comparison-table-wrap"><table className="comparison-table">
      <thead><tr><th>Change</th><th>Item</th><th>Parameter</th><th><span>Reference</span><small>{reference}</small></th><th><span>Comparison</span><small>{comparison}</small></th><th>Δ</th></tr></thead>
      <tbody>{section.changes.map((change, index) => <tr key={`${change.item}:${change.parameter}:${index}`}>
        <td><span className={`difference-status ${change.status}`}>{statusLabel(change.status)}</span></td>
        <th scope="row"><span>{change.item}</span>{change.country && change.country !== change.item && <small>{change.country}</small>}</th>
        <td>{change.parameter}</td>
        <td className={change.status === "removed" ? "difference-emphasis removed" : ""}>{formatValue(change.reference)}</td>
        <td className={change.status === "added" ? "difference-emphasis added" : ""}>{formatValue(change.comparison)}</td>
        <td className="difference-delta">{formatDelta(change)}</td>
      </tr>)}</tbody>
    </table></div>
  </section>;
}
