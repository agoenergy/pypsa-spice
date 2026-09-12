import styles from "./Plot.module.scss";
import chartSurfaceStyles from "./ChartSurface.module.scss";
import { useEffect, useMemo, useRef, useState } from "react";

import { ChartLegend } from "./ChartLegend";
import { getLegendValues } from "../shared/chartData";
import { chartFont } from "../shared/chartPresentation";
import { buildPlotTraces } from "../shared/plotTraces";
import { loadPlotly } from "../plotly";
import type { Catalog, ChartDefinition, ChartResponse } from "../types";

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
    const { data, hasColumnTotals } = buildPlotTraces({
      primary,
      comparison,
      chart,
      mappings,
      legendValues,
      hiddenLegendValues,
      difference,
    });
    plotly.react(
      plotElement,
      data,
      {
        autosize: true,
        margin: { l: 66, r: chart.secondary_y_lab ? 66 : 20, t: hasColumnTotals ? 34 : 16, b: 46 },
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
