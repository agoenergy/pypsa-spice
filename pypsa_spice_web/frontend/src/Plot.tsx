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
  expanded: boolean;
  difference?: boolean;
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

function traces(
  response: ChartResponse,
  chart: ChartDefinition,
  scenario: string,
  mappings: Catalog["mappings"],
  comparison: boolean,
) {
  return [...aggregate(response.rows, chart).entries()].map(([legend, points], index) => {
    const color = mappings[legend]?.color || fallbackColors[index % fallbackColors.length];
    const isArea = chart.type === "area_share";
    const isBar = chart.type.includes("bar");
    const trace: Record<string, unknown> = {
      name: `${pretty(legend, mappings)}${comparison ? ` · ${scenario}` : ""}`,
      x: points.map((point) => point.x),
      y: points.map((point) => point.y),
      marker: { color },
      line: { color, width: comparison ? 1.5 : 2, dash: comparison ? "dot" : "solid" },
      opacity: comparison ? 0.5 : 0.94,
      hovertemplate: `<b>${pretty(legend, mappings)}</b><br>%{x}<br>%{y:,.2f} ${chart.units || ""}<extra>${scenario}</extra>`,
      yaxis: chart.secondary_y_lab?.includes(legend) ? "y2" : "y",
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
  primaryName: string,
  comparisonName: string,
  mappings: Catalog["mappings"],
) {
  const first = aggregate(primary.rows, chart);
  const second = aggregate(comparison.rows, chart);
  const legends = [...new Set([...first.keys(), ...second.keys()])];
  return legends.map((legend, index) => {
    const firstPoints = new Map((first.get(legend) || []).map((point) => [String(point.x), point.y]));
    const secondPoints = new Map((second.get(legend) || []).map((point) => [String(point.x), point.y]));
    const xValues = [...new Set([...firstPoints.keys(), ...secondPoints.keys()])].sort();
    const color = mappings[legend]?.color || fallbackColors[index % fallbackColors.length];
    const isBar = chart.type.includes("bar") || !chart.hourly;
    return {
      type: isBar ? "bar" : "scatter",
      mode: isBar ? undefined : "lines",
      name: pretty(legend, mappings),
      x: xValues.map((value) => chart.hourly ? value : Number(value)),
      y: xValues.map((value) => (secondPoints.get(value) || 0) - (firstPoints.get(value) || 0)),
      marker: { color }, line: { color, width: 2 },
      hovertemplate: `<b>${pretty(legend, mappings)}</b><br>%{x}<br>Difference: %{y:+,.2f} ${chart.units || ""}<extra>${comparisonName} − ${primaryName}</extra>`,
    };
  });
}

export default function Plot({ chart, primary, comparison, primaryName, comparisonName, mappings, expanded, difference = false }: Props) {
  const ref = useRef<HTMLDivElement>(null);
  const legendValues = useMemo(() => {
    const values = new Set(aggregate(primary.rows, chart).keys());
    if (comparison) {
      for (const value of aggregate(comparison.rows, chart).keys()) values.add(value);
    }
    return [...values];
  }, [primary, comparison, chart]);
  useEffect(() => {
    if (!ref.current || !window.Plotly) return;
    const dark = document.documentElement.dataset.theme === "dark";
    const grid = dark ? "#34403c" : "#e2e6e4";
    const text = dark ? "#a9b5b1" : "#65717d";
    const allTraces = difference && comparison
      ? differenceTraces(primary, comparison, chart, primaryName, comparisonName, mappings)
      : [...traces(primary, chart, primaryName, mappings, false), ...(comparison ? traces(comparison, chart, comparisonName, mappings, true) : [])];
    window.Plotly.react(ref.current, allTraces, {
      autosize: true,
      margin: { l: 58, r: chart.secondary_y_lab ? 58 : 18, t: 14, b: 38 },
      paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
      font: { family: "Flexo, sans-serif", size: 10, color: text },
      showlegend: false,
      hovermode: "x unified", barmode: difference ? "relative" : comparison ? "group" : "relative",
      xaxis: { gridcolor: grid, zeroline: false, tickfont: { size: 9 } },
      yaxis: { title: { text: difference ? `Difference (${chart.units || "value"})` : chart.units || "", font: { size: 10 } }, gridcolor: grid, zeroline: difference, zerolinecolor: text, zerolinewidth: 1, rangemode: "tozero" },
      yaxis2: { overlaying: "y", side: "right", showgrid: false, title: "State of charge" },
      hoverlabel: { bgcolor: dark ? "#222b28" : "#fff", bordercolor: grid },
      uirevision: `${chart.id}-${primaryName}-${difference}`,
    }, {
      responsive: true, displaylogo: false,
      modeBarButtonsToRemove: ["lasso2d", "select2d"],
      toImageButtonOptions: { format: "png", filename: `${primaryName}_${chart.table_name}`, scale: 2 },
    });
    return () => { if (ref.current) window.Plotly.purge(ref.current); };
  }, [chart, primary, comparison, primaryName, comparisonName, mappings, difference]);

  useEffect(() => { if (ref.current && window.Plotly) window.Plotly.Plots.resize(ref.current); }, [expanded]);
  return <div className="plot-with-legend">
    <div ref={ref} className="plot" />
    <div className="html-legend" aria-label="Chart legend">
      {legendValues.map((value, index) => <span className="html-legend-item" key={value}>
        <i style={{ backgroundColor: mappings[value]?.color || fallbackColors[index % fallbackColors.length] }} aria-hidden="true" />
        <span>{pretty(value, mappings)}</span>
      </span>)}
    </div>
  </div>;
}
