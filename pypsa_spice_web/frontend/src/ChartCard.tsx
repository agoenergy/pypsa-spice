import { useEffect, useMemo, useState } from "react";
import { Download, Expand, Minimize2, RotateCcw, Table2 } from "lucide-react";
import { downloadUrl, getChart } from "./api";
import Plot, { ChartLegend, getLegendValues } from "./Plot";
import type { Catalog, ChartDefinition, ChartResponse, ResultRow, Selection } from "./types";

interface Props {
  chart: ChartDefinition;
  selection: Selection;
  years: string[];
  onYearChange: (year: string) => void;
  mappings: Catalog["mappings"];
  darkMode: boolean;
  onInspect: (title: string, rows: ResultRow[], sourceCount: number) => void;
}

const HOUR_MS = 60 * 60 * 1000;

function timestampToMs(value: string | null | undefined): number | null {
  if (!value) return null;
  const timestamp = new Date(value.replace(" ", "T")).getTime();
  return Number.isFinite(timestamp) ? timestamp : null;
}

function msToTimestamp(value: number): string {
  const date = new Date(value);
  const pad = (part: number) => String(part).padStart(2, "0");
  return `${date.getFullYear()}-${pad(date.getMonth() + 1)}-${pad(date.getDate())}T${pad(date.getHours())}:${pad(date.getMinutes())}`;
}

function readableTimestamp(value: number): string {
  const date = new Date(value);
  const day = date.toLocaleDateString(undefined, { day: "numeric", month: "short", year: "numeric" });
  const time = date.toLocaleTimeString(undefined, { hour: "2-digit", minute: "2-digit", hour12: false });
  return `${day} · ${time}`;
}

export default function ChartCard({ chart, selection, years, onYearChange, mappings, darkMode, onInspect }: Props) {
  const [primary, setPrimary] = useState<ChartResponse | null>(null);
  const [comparison, setComparison] = useState<ChartResponse | null>(null);
  const [country, setCountry] = useState("ALL");
  const [filterValue, setFilterValue] = useState("ALL");
  const [error, setError] = useState("");
  const [loading, setLoading] = useState(true);
  const [expanded, setExpanded] = useState(false);
  const [showDifference, setShowDifference] = useState(false);
  const [startTime, setStartTime] = useState("");
  const [endTime, setEndTime] = useState("");

  useEffect(() => {
    setPrimary(null); setComparison(null);
    setCountry("ALL"); setFilterValue("ALL"); setStartTime(""); setEndTime("");
  }, [selection.dataset, selection.project, selection.scenario, selection.sector]);
  useEffect(() => { setComparison(null); setShowDifference(false); }, [selection.comparison]);

  useEffect(() => {
    const controller = new AbortController();
    setLoading(true); setError("");
    Promise.all([
      getChart(chart, selection, selection.scenario, country, filterValue, startTime, endTime, controller.signal),
      selection.comparison ? getChart(chart, selection, selection.comparison, country, filterValue, startTime, endTime, controller.signal) : Promise.resolve(null),
    ]).then(([first, second]) => { setPrimary(first); setComparison(second); })
      .catch((reason) => { if (reason.name !== "AbortError") setError(reason.message); })
      .finally(() => { if (!controller.signal.aborted) setLoading(false); });
    return () => controller.abort();
  }, [chart, selection.dataset, selection.project, selection.scenario, selection.comparison, selection.sector, chart.hourly ? selection.year : "", country, filterValue, startTime, endTime]);

  const rows = useMemo(() => [
    ...(primary?.rows.map((row) => ({ ...row, scenario: selection.scenario })) || []),
    ...(comparison?.rows.map((row) => ({ ...row, scenario: selection.comparison })) || []),
  ], [primary, comparison, selection]);
  const sourceCount = (primary?.meta.source_rows || 0) + (comparison?.meta.source_rows || 0);
  const countries = primary?.dimensions.country || [];
  const filters = chart.fil_col ? primary?.dimensions[chart.fil_col] || [] : [];
  const filterLabel = chart.fil_col === "to" ? "destinations" : `${chart.fil_col}s`;

  const comparing = Boolean(selection.comparison && comparison);
  const comparisonLegendValues = useMemo(
    () => getLegendValues(chart, primary, comparison),
    [chart, primary, comparison],
  );
  const availableStart = timestampToMs(primary?.meta.available_start);
  const availableEnd = timestampToMs(primary?.meta.available_end);
  const selectedStart = timestampToMs(startTime) ?? availableStart;
  const selectedEnd = timestampToMs(endTime) ?? availableEnd;
  const hasTimeRange = Boolean(startTime || endTime);
  const rangeSpan = availableStart !== null && availableEnd !== null ? Math.max(1, availableEnd - availableStart) : 1;
  const minimumRange = Math.min(HOUR_MS, rangeSpan);
  const startPercent = availableStart !== null && selectedStart !== null ? ((selectedStart - availableStart) / rangeSpan) * 100 : 0;
  const endPercent = availableStart !== null && selectedEnd !== null ? ((selectedEnd - availableStart) / rangeSpan) * 100 : 100;

  const changeStart = (value: number) => {
    if (selectedEnd === null) return;
    setStartTime(msToTimestamp(Math.min(value, selectedEnd - minimumRange)));
  };
  const changeEnd = (value: number) => {
    if (selectedStart === null) return;
    setEndTime(msToTimestamp(Math.max(value, selectedStart + minimumRange)));
  };
  const changeYear = (year: string) => {
    setStartTime("");
    setEndTime("");
    onYearChange(year);
  };

  return <article className={`chart-card ${expanded ? "expanded" : ""} ${comparing ? "comparing" : ""} ${showDifference ? "showing-difference" : ""}`}>
    <header className="chart-head">
      <div className="chart-title"><h3>{chart.name}</h3></div>
      <div className="chart-toolbar">
        <div className="chart-controls">
          {countries.length > 1 && <select aria-label={`Country for ${chart.name}`} value={country} onChange={(event) => setCountry(event.target.value)}><option value="ALL">All countries</option>{countries.map((item) => <option key={item}>{item}</option>)}</select>}
          {chart.fil_col && filters.length > 0 && <select aria-label={`${chart.fil_col} for ${chart.name}`} value={filterValue} onChange={(event) => setFilterValue(event.target.value)}><option value="ALL">All {filterLabel}</option>{filters.map((item) => <option key={item}>{item}</option>)}</select>}
          {chart.hourly && years.length > 0 && <select aria-label={`Hourly year for ${chart.name}`} value={selection.year} onChange={(event) => changeYear(event.target.value)}>{years.map((year) => <option key={year}>{year}</option>)}</select>}
          {selection.comparison && <label className="chart-difference-toggle" title={`${selection.comparison} − ${selection.scenario}`}><input type="checkbox" checked={showDifference} onChange={(event) => setShowDifference(event.target.checked)} /><i aria-hidden="true" /><span>Difference</span></label>}
        </div>
        <div className="chart-actions">
          <button title="View source data" aria-label={`View ${chart.name} source data`} onClick={() => onInspect(chart.name, rows, sourceCount)}><Table2 aria-hidden="true" /></button>
          <a title="Download complete CSV" aria-label={`Download ${chart.name} CSV`} href={downloadUrl(chart, selection)} download><Download aria-hidden="true" /></a>
          <button title={expanded ? "Close expanded chart" : "Expand chart"} aria-label={`${expanded ? "Close" : "Expand"} ${chart.name}`} onClick={() => setExpanded(!expanded)}>{expanded ? <Minimize2 aria-hidden="true" /> : <Expand aria-hidden="true" />}</button>
        </div>
      </div>
    </header>
    {chart.hourly && availableStart !== null && availableEnd !== null && selectedStart !== null && selectedEnd !== null && <div className="time-range-controls" aria-label={`Time range for ${chart.name}`}>
      <div className="time-range-head"><span>Time range</span>{hasTimeRange && <button type="button" onClick={() => { setStartTime(""); setEndTime(""); }}><RotateCcw aria-hidden="true" />Reset</button>}</div>
      <div className="time-range-values" aria-live="polite"><output>{readableTimestamp(selectedStart)}</output><output>{readableTimestamp(selectedEnd)}</output></div>
      <div className="dual-range">
        <div className="range-track" aria-hidden="true"><i style={{ left: `${startPercent}%`, right: `${100 - endPercent}%` }} /></div>
        <input type="range" aria-label={`Start time for ${chart.name}`} min={availableStart} max={availableEnd} step={HOUR_MS} value={selectedStart} onInput={(event) => changeStart(Number(event.currentTarget.value))} style={{ zIndex: startPercent > 50 ? 4 : 2 }} />
        <input type="range" aria-label={`End time for ${chart.name}`} min={availableStart} max={availableEnd} step={HOUR_MS} value={selectedEnd} onInput={(event) => changeEnd(Number(event.currentTarget.value))} />
      </div>
    </div>}
    <div className={`chart-body ${primary && primary.rows.length > 0 ? "with-plot" : ""} ${comparing && !showDifference ? "comparison-layout" : ""}`} aria-busy={loading}>
      {loading && !primary && <div className="state"><span className="spinner" />Reading result table…</div>}
      {!loading && error && <div className="state empty"><b>No chart data</b><span>{error}</span></div>}
      {!loading && !error && primary && primary.rows.length === 0 && <div className="state empty"><b>No values in this result table</b><span>The chart will appear when the selected run contains data.</span></div>}
      {!error && primary && primary.rows.length > 0 && !comparing && <Plot chart={chart} primary={primary} comparison={null} primaryName={selection.scenario} comparisonName="" mappings={mappings} darkMode={darkMode} expanded={expanded} />}
      {!error && primary && primary.rows.length > 0 && comparing && !showDifference && <>
        <div className="scenario-plot"><ChartContextHeading label="Primary scenario" title={selection.scenario} /><Plot chart={chart} primary={primary} comparison={null} primaryName={selection.scenario} comparisonName="" mappings={mappings} darkMode={darkMode} expanded={expanded} legendValues={comparisonLegendValues} showLegend={false} /></div>
        <div className="scenario-plot"><ChartContextHeading label="Comparison scenario" title={selection.comparison} /><Plot chart={chart} primary={comparison!} comparison={null} primaryName={selection.comparison} comparisonName="" mappings={mappings} darkMode={darkMode} expanded={expanded} legendValues={comparisonLegendValues} showLegend={false} /></div>
        <ChartLegend values={comparisonLegendValues} mappings={mappings} />
      </>}
      {!error && primary && primary.rows.length > 0 && comparing && showDifference && <div className="difference-plot"><ChartContextHeading label="Difference" title={`${selection.comparison} − ${selection.scenario}`} /><Plot chart={chart} primary={primary} comparison={comparison} primaryName={selection.scenario} comparisonName={selection.comparison} mappings={mappings} darkMode={darkMode} expanded={expanded} difference /></div>}
    </div>
  </article>;
}

function ChartContextHeading({ label, title }: { label: string; title: string }) {
  return <div className="scenario-label"><small>{label}</small><h4>{title}</h4></div>;
}
