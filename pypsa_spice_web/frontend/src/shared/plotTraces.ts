import type { PlotData } from "plotly.js";
import type { PlotTrace, ScatterTextTrace } from "./types";
import type { Catalog, ChartDefinition, ChartResponse } from "../types";
import { aggregate, differenceAggregates } from "./chartData";
import { chartFont, formatLegendLabel, getLegendColour } from "./chartPresentation";

interface PlotTraceOptions {
  primary: ChartResponse;
  comparison: ChartResponse | null;
  chart: ChartDefinition;
  mappings: Catalog["mappings"];
  legendValues: string[];
  hiddenLegendValues: ReadonlySet<string>;
  difference: boolean;
}

export function buildPlotTraces({
  primary,
  comparison,
  chart,
  mappings,
  legendValues,
  hiddenLegendValues,
  difference,
}: PlotTraceOptions): { data: PlotTrace[]; hasColumnTotals: boolean } {
  const chartTraces =
    difference && comparison
      ? differenceTraces(primary, comparison, chart, mappings, legendValues, hiddenLegendValues)
      : [
          ...traces(primary, chart, mappings, false, legendValues, hiddenLegendValues),
          ...(comparison ? traces(comparison, chart, mappings, true, legendValues, hiddenLegendValues) : []),
        ];
  const totalTrace = difference ? null : stackedBarTotalTrace(primary, chart, hiddenLegendValues);
  return {
    data: totalTrace ? [...chartTraces, totalTrace] : chartTraces,
    hasColumnTotals: totalTrace !== null,
  };
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
): Partial<PlotData>[] {
  return [...aggregate(response.rows, chart).entries()].map(([legend, points]): Partial<PlotData> => {
    const color = getLegendColour(legend, legendValues.indexOf(legend), mappings);
    const isArea = chart.type === "area_share";
    const isBar = chart.type.includes("bar");
    const trace: Partial<PlotData> = {
      name: formatLegendLabel(legend, mappings),
      x: points.map((point) => point.x),
      y: points.map((point) => point.y),
      marker: { color },
      line: { color, width: comparison ? 1.5 : 2, dash: comparison ? "dot" : "solid" },
      opacity: comparison ? 0.5 : 0.94,
      hovertemplate: `<b>${formatLegendLabel(legend, mappings)}</b>: %{y:,.2f} ${chart.units || ""}<extra></extra>`,
      yaxis: chart.secondary_y_lab?.includes(legend) ? "y2" : "y",
      visible: hiddenLegendValues.has(legend) ? "legendonly" : true,
    };
    if (isArea) {
      trace.type = "scatter";
      trace.mode = "lines";
      trace.fill = comparison ? "none" : "tonexty";
      trace.stackgroup = comparison ? undefined : "one";
    } else if (chart.hourly && !isBar) {
      trace.type = "scatter";
      trace.mode = "lines";
    } else trace.type = "bar";
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
): Partial<PlotData>[] {
  return [...differenceAggregates(primary, comparison, chart)].map(([legend, points]): Partial<PlotData> => {
    const color = getLegendColour(legend, legendValues.indexOf(legend), mappings);
    const isBar = chart.type.includes("bar") || !chart.hourly;
    return {
      type: isBar ? "bar" : "scatter",
      mode: isBar ? undefined : "lines",
      name: formatLegendLabel(legend, mappings),
      x: points.map((point) => point.x),
      y: points.map((point) => point.y),
      marker: { color },
      line: { color, width: 2 },
      hovertemplate: `<b>${formatLegendLabel(legend, mappings)}</b>: %{y:+,.2f} ${chart.units || ""}<extra></extra>`,
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
): ScatterTextTrace | null {
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
    y: values.map((value) => (value.positive > 0 ? value.positive : value.negative)),
    text: values.map((value) => formatTotal(value.total)),
    textposition: values.map((value) => (value.positive > 0 ? "top center" : "bottom center")),
    textfont: { size: chartFont.body },
    cliponaxis: false,
    hoverinfo: "skip",
    showlegend: false,
  };
}
