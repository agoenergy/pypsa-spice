import { useEffect, useMemo, useRef } from "react";
import type { Catalog, ChartDefinition, ChartResponse, ResultRow } from "./types";

declare global {
  interface Window { Plotly: any }
}

const fallbackColors = ["#e6007e", "#005ca9", "#60a917", "#ec6608", "#7553a6", "#009e8e", "#c33c54", "#79848d", "#d5a400", "#3f7c85", "#9b4b96", "#86a6c2"];

interface Props {
  chart: ChartDefinition;
  primary: ChartResponse;
  comparison: ChartResponse | null;
  primaryName: string;
  comparisonName: string;
  mappings: Catalog["mappings"];
  darkMode: boolean;
  expanded: boolean;
  difference?: boolean;
  legendValues?: string[];
  showLegend?: boolean;
  hiddenLegendValues: ReadonlySet<string>;
  onLegendToggle: (value: string) => void;
}

function pretty(value: string, mappings: Catalog["mappings"]): string {
  return mappings[value]?.label || value.replaceAll("_", " ").replace(/\b\w/g, (letter) => letter.toUpperCase());
}

function aggregate(rows: ResultRow[], chart: ChartDefinition): Map<string, { x: string | number; y: number }[]> {
  const xKey = chart.hourly ? "snapshot" : "year";
  const values = new Map<string, number>();
  for (const row of rows) {
    const x = row[xKey];
    const legend = String(row[chart.leg_col] ?? "Series");
    if (x === null || x === undefined) continue;
    const key = `${x}\u0000${legend}`;
    values.set(key, (values.get(key) || 0) + Number(row.value || 0));
  }
  const groups = new Map<string, { x: string | number; y: number }[]>();
  for (const [key, y] of values) {
    const [rawX, legend] = key.split("\u0000");
    if (!groups.has(legend)) groups.set(legend, []);
    groups.get(legend)!.push({ x: chart.hourly ? rawX : Number(rawX), y });
  }
  for (const points of groups.values()) points.sort((a, b) => String(a.x).localeCompare(String(b.x)));
  return groups;
}

export function getLegendValues(chart: ChartDefinition, ...responses: (ChartResponse | null)[]): string[] {
  const values = new Set<string>();
  for (const response of responses) {
    if (!response) continue;
    for (const value of aggregate(response.rows, chart).keys()) values.add(value);
  }
  return [...values];
}

function legendColor(value: string, index: number, mappings: Catalog["mappings"]): string {
  return mappings[value]?.color || fallbackColors[Math.max(0, index) % fallbackColors.length];
}

export function ChartLegend({ values, mappings, hiddenValues, onToggle }: {
  values: string[];
  mappings: Catalog["mappings"];
  hiddenValues: ReadonlySet<string>;
  onToggle: (value: string) => void;
}) {
  return <div className="html-legend" aria-label="Chart legend">
    {values.map((value, index) => <button
      type="button"
      className={`html-legend-item ${hiddenValues.has(value) ? "is-hidden" : ""}`}
      key={value}
      aria-pressed={!hiddenValues.has(value)}
      title={`${hiddenValues.has(value) ? "Show" : "Hide"} ${pretty(value, mappings)}`}
      onClick={() => onToggle(value)}
    >
      <i style={{ backgroundColor: legendColor(value, index, mappings) }} aria-hidden="true" />
      <span>{pretty(value, mappings)}</span>
    </button>)}
  </div>;
}

function isSecondarySeries(chart: ChartDefinition, legend: string): boolean {
  return Boolean(chart.secondary_y_lab?.includes(legend));
}

function traces(
  response: ChartResponse,
  chart: ChartDefinition,
  mappings: Catalog["mappings"],
  comparison: boolean,
  legendValues: string[],
  hiddenLegendValues: ReadonlySet<string>,
) {
  return [...aggregate(response.rows, chart).entries()].map(([legend, points]) => {
    const color = legendColor(legend, legendValues.indexOf(legend), mappings);
    const isArea = chart.type === "area_share";
    const isBar = chart.type.includes("bar");
    const trace: Record<string, unknown> = {
      name: pretty(legend, mappings),
      x: points.map((point) => point.x),
      y: points.map((point) => point.y),
      marker: { color },
      line: { color, width: comparison ? 1.5 : 2, dash: comparison ? "dot" : "solid" },
      opacity: comparison ? 0.5 : 0.94,
      hovertemplate: `<b>${pretty(legend, mappings)}</b>: %{y:,.2f} ${chart.units || ""}<extra></extra>`,
      yaxis: chart.secondary_y_lab?.includes(legend) ? "y2" : "y",
      visible: hiddenLegendValues.has(legend) ? "legendonly" : true,
    };
    if (isArea) Object.assign(trace, { type: "scatter", mode: "lines", fill: comparison ? "none" : "tonexty", stackgroup: comparison ? undefined : "one" });
    else if (chart.hourly && !isBar) Object.assign(trace, { type: "scatter", mode: "lines" });
    else trace.type = "bar";
    return trace;
  });
}

function differenceTraces(
  primary: ChartResponse,
  comparison: ChartResponse,
  chart: ChartDefinition,
  mappings: Catalog["mappings"],
  legendValues: string[],
  hiddenLegendValues: ReadonlySet<string>,
) {
  const first = aggregate(primary.rows, chart);
  const second = aggregate(comparison.rows, chart);
  const legends = [...new Set([...first.keys(), ...second.keys()])];
  return legends.map((legend) => {
    const firstPoints = new Map((first.get(legend) || []).map((point) => [String(point.x), point.y]));
    const secondPoints = new Map((second.get(legend) || []).map((point) => [String(point.x), point.y]));
    const xValues = [...new Set([...firstPoints.keys(), ...secondPoints.keys()])].sort();
    const color = legendColor(legend, legendValues.indexOf(legend), mappings);
    const isBar = chart.type.includes("bar") || !chart.hourly;
    return {
      type: isBar ? "bar" : "scatter",
      mode: isBar ? undefined : "lines",
      name: pretty(legend, mappings),
      x: xValues.map((value) => chart.hourly ? value : Number(value)),
      y: xValues.map((value) => (secondPoints.get(value) || 0) - (firstPoints.get(value) || 0)),
      marker: { color }, line: { color, width: 2 },
      hovertemplate: `<b>${pretty(legend, mappings)}</b>: %{y:+,.2f} ${chart.units || ""}<extra></extra>`,
      visible: hiddenLegendValues.has(legend) ? "legendonly" : true,
    };
  });
}

function formatTotal(value: number): string {
  return new Intl.NumberFormat(undefined, {
    notation: Math.abs(value) >= 10_000 ? "compact" : "standard",
    maximumFractionDigits: 2,
  }).format(value);
}

function stackedBarTotalTrace(
  response: ChartResponse,
  chart: ChartDefinition,
  hiddenLegendValues: ReadonlySet<string>,
) {
  if (chart.hourly || chart.type === "grouped_bar" || !chart.type.includes("bar")) return null;
  const totals = new Map<string, { x: string | number; total: number; positive: number; negative: number }>();
  for (const [legend, points] of aggregate(response.rows, chart)) {
    if (hiddenLegendValues.has(legend) || isSecondarySeries(chart, legend)) continue;
    for (const point of points) {
      const key = String(point.x);
      const current = totals.get(key) || { x: point.x, total: 0, positive: 0, negative: 0 };
      current.total += point.y;
      if (point.y >= 0) current.positive += point.y;
      else current.negative += point.y;
      totals.set(key, current);
    }
  }
  const values = [...totals.values()].sort((a, b) => String(a.x).localeCompare(String(b.x)));
  if (values.length === 0) return null;
  return {
    type: "scatter",
    mode: "text",
    name: "Column total",
    x: values.map((value) => value.x),
    y: values.map((value) => value.positive > 0 ? value.positive : value.negative),
    text: values.map((value) => formatTotal(value.total)),
    textposition: values.map((value) => value.positive > 0 ? "top center" : "bottom center"),
    textfont: { size: 10 },
    cliponaxis: false,
    hoverinfo: "skip",
    showlegend: false,
  };
}

export default function Plot({ chart, primary, comparison, primaryName, comparisonName, mappings, darkMode, expanded, difference = false, legendValues: sharedLegendValues, showLegend = true, hiddenLegendValues, onLegendToggle }: Props) {
  const ref = useRef<HTMLDivElement>(null);
  const derivedLegendValues = useMemo(() => getLegendValues(chart, primary, comparison), [primary, comparison, chart]);
  const legendValues = sharedLegendValues || derivedLegendValues;
  useEffect(() => {
    if (!ref.current || !window.Plotly) return;
    const grid = "#e2e6e4";
    const text = darkMode ? "#a9b5b1" : "#65717d";
    const chartTraces = difference && comparison
      ? differenceTraces(primary, comparison, chart, mappings, legendValues, hiddenLegendValues)
      : [...traces(primary, chart, mappings, false, legendValues, hiddenLegendValues), ...(comparison ? traces(comparison, chart, mappings, true, legendValues, hiddenLegendValues) : [])];
    const totalTrace = difference ? null : stackedBarTotalTrace(primary, chart, hiddenLegendValues);
    const allTraces = totalTrace ? [...chartTraces, totalTrace] : chartTraces;
    window.Plotly.react(ref.current, allTraces, {
      autosize: true,
      margin: { l: 58, r: chart.secondary_y_lab ? 58 : 18, t: totalTrace ? 28 : 14, b: 38 },
      paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
      font: { family: "Flexo, sans-serif", size: 10, color: text },
      showlegend: false,
      hovermode: "x unified", barmode: chart.type === "grouped_bar" || comparison ? "group" : "relative",
      xaxis: { showgrid: !darkMode, gridcolor: grid, zeroline: false, tickfont: { size: 9 }, unifiedhovertitle: { text: chart.hourly ? "%{x|%d %b · %H:%M}" : "%{x}" } },
      yaxis: { title: { text: difference ? `Difference (${chart.units || "value"})` : chart.units || "", font: { size: 10 } }, showgrid: !darkMode, gridcolor: grid, zeroline: difference, zerolinecolor: darkMode ? "rgba(255,255,255,.18)" : text, zerolinewidth: 1, rangemode: "tozero" },
      yaxis2: { overlaying: "y", side: "right", showgrid: false, title: "State of charge" },
      hoverlabel: { bgcolor: darkMode ? "#222b28" : "#fff", bordercolor: darkMode ? "#34403c" : grid, font: { size: 10 }, align: "left" },
      uirevision: `${chart.id}-${primaryName}-${difference}`,
    }, {
      responsive: true, displaylogo: false,
      modeBarButtonsToRemove: ["lasso2d", "select2d"],
      toImageButtonOptions: { format: "png", filename: `${primaryName}_${chart.table_name}`, scale: 2 },
    });
    return () => { if (ref.current) window.Plotly.purge(ref.current); };
  }, [chart, primary, comparison, primaryName, comparisonName, mappings, darkMode, difference, legendValues, hiddenLegendValues]);

  useEffect(() => { if (ref.current && window.Plotly) window.Plotly.Plots.resize(ref.current); }, [expanded]);
  return <div className="plot-with-legend">
    <div ref={ref} className="plot" />
    {showLegend && <ChartLegend values={legendValues} mappings={mappings} hiddenValues={hiddenLegendValues} onToggle={onLegendToggle} />}
  </div>;
}
