import styles from "./Plot.module.scss";
import chartSurfaceStyles from "./ChartSurface.module.scss";
import { useEffect, useMemo, useRef, useState } from "react";

import { ChartLegend } from "./ChartLegend";
import { aggregate, differenceAggregates, getLegendValues } from "../shared/chartData";
import { formatLegendLabel, getLegendColour } from "../shared/chartPresentation";
import { loadPlotly } from "../plotly";
import type { Catalog, ChartDefinition, ChartResponse } from "../types";

// Plotly draws its own text, so it cannot read the CSS type scale in global.scss.
// These mirror --text-xs and --text-sm so chart type matches the surrounding interface.
const chartFont = { body: 12, hover: 13 };

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
    const color = getLegendColour(legend, legendValues.indexOf(legend), mappings);
    const isArea = chart.type === "area_share";
    const isBar = chart.type.includes("bar");
    const trace: Record<string, unknown> = {
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
    if (isArea)
      Object.assign(trace, {
        type: "scatter",
        mode: "lines",
        fill: comparison ? "none" : "tonexty",
        stackgroup: comparison ? undefined : "one",
      });
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
  return [...differenceAggregates(primary, comparison, chart)].map(([legend, points]) => {
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
    y: values.map((value) => (value.positive > 0 ? value.positive : value.negative)),
    text: values.map((value) => formatTotal(value.total)),
    textposition: values.map((value) => (value.positive > 0 ? "top center" : "bottom center")),
    textfont: { size: chartFont.body },
    cliponaxis: false,
    hoverinfo: "skip",
    showlegend: false,
  };
}

export default function Plot({
  chart,
  primary,
  comparison,
  primaryName,
  comparisonName,
  mappings,
  darkMode,
  expanded,
  difference = false,
  legendValues: sharedLegendValues,
  showLegend = true,
  hiddenLegendValues,
  onLegendToggle,
}: Props) {
  const ref = useRef<HTMLDivElement>(null);
  const [plotlyReady, setPlotlyReady] = useState(() => Boolean(window.Plotly));
  const [plotlyError, setPlotlyError] = useState("");
  const derivedLegendValues = useMemo(() => getLegendValues(chart, primary, comparison), [primary, comparison, chart]);
  const legendValues = sharedLegendValues || derivedLegendValues;
  useEffect(() => {
    let current = true;
    loadPlotly()
      .then(() => {
        if (current) setPlotlyReady(true);
      })
      .catch((reason) => {
        if (!current) return;
        setPlotlyError(reason instanceof Error ? reason.message : "The chart library could not be loaded.");
      });
    return () => {
      current = false;
    };
  }, []);
  useEffect(() => {
    if (!plotlyReady || !ref.current || !window.Plotly) return;
    const plotly = window.Plotly;
    const plotElement = ref.current;
    const grid = "#e2e6e4";
    const text = darkMode ? "#a9b5b1" : "#65717d";
    const chartTraces =
      difference && comparison
        ? differenceTraces(primary, comparison, chart, mappings, legendValues, hiddenLegendValues)
        : [
            ...traces(primary, chart, mappings, false, legendValues, hiddenLegendValues),
            ...(comparison ? traces(comparison, chart, mappings, true, legendValues, hiddenLegendValues) : []),
          ];
    const totalTrace = difference ? null : stackedBarTotalTrace(primary, chart, hiddenLegendValues);
    const allTraces = totalTrace ? [...chartTraces, totalTrace] : chartTraces;
    plotly.react(
      plotElement,
      allTraces,
      {
        autosize: true,
        margin: { l: 66, r: chart.secondary_y_lab ? 66 : 20, t: totalTrace ? 34 : 16, b: 46 },
        paper_bgcolor: "rgba(0,0,0,0)",
        plot_bgcolor: "rgba(0,0,0,0)",
        font: { family: "Flexo, sans-serif", size: chartFont.body, color: text },
        showlegend: false,
        hovermode: "x unified",
        barmode: chart.type === "grouped_bar" || comparison ? "group" : "relative",
        xaxis: {
          showgrid: !darkMode,
          gridcolor: grid,
          zeroline: false,
          tickfont: { size: chartFont.body },
          unifiedhovertitle: { text: chart.hourly ? "%{x|%d %b · %H:%M}" : "%{x}" },
        },
        yaxis: {
          title: {
            text: difference ? `Difference (${chart.units || "value"})` : chart.units || "",
            font: { size: chartFont.body },
          },
          showgrid: !darkMode,
          gridcolor: grid,
          zeroline: difference,
          zerolinecolor: darkMode ? "rgba(255,255,255,.18)" : text,
          zerolinewidth: 1,
          rangemode: "tozero",
        },
        yaxis2: { overlaying: "y", side: "right", showgrid: false, title: "State of charge" },
        hoverlabel: {
          bgcolor: darkMode ? "#222b28" : "#fff",
          bordercolor: darkMode ? "#34403c" : grid,
          font: { size: chartFont.hover },
          align: "left",
        },
        uirevision: `${chart.id}-${primaryName}-${difference}`,
      },
      {
        responsive: true,
        displaylogo: false,
        modeBarButtonsToRemove: ["lasso2d", "select2d"],
        toImageButtonOptions: { format: "png", filename: `${primaryName}_${chart.table_name}`, scale: 2 },
      },
    );
    return () => {
      plotly.purge(plotElement);
    };
  }, [
    chart,
    primary,
    comparison,
    primaryName,
    comparisonName,
    mappings,
    darkMode,
    difference,
    legendValues,
    hiddenLegendValues,
    plotlyReady,
  ]);

  useEffect(() => {
    if (ref.current && window.Plotly) window.Plotly.Plots.resize(ref.current);
  }, [expanded]);
  return (
    <div className={styles["plot-with-legend"]}>
      <div ref={ref} className={chartSurfaceStyles["plot"]}>
        {!plotlyReady && !plotlyError && <span className={styles["plot-message"]}>Loading chart…</span>}
        {plotlyError && <span className={[styles["plot-message"], styles["error"]].join(" ")}>{plotlyError}</span>}
      </div>
      {showLegend && (
        <ChartLegend
          values={legendValues}
          mappings={mappings}
          hiddenValues={hiddenLegendValues}
          onToggle={onLegendToggle}
        />
      )}
    </div>
  );
}
