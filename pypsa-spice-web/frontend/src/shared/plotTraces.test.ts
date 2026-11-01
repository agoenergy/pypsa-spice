import { describe, expect, it } from "vitest";
import { buildPlotTraces } from "./plotTraces";
import type { ChartDefinition, ChartResponse, ResultRow } from "../types";

const chart: ChartDefinition = {
  id: "capacity",
  name: "Capacity",
  units: "GW",
  summary: "sum",
  table_name: "capacity_yearly",
  leg_col: "technology",
  type: "bar",
  hourly: false,
};

function response(rows: ResultRow[]): ChartResponse {
  return {
    rows,
    dimensions: {},
    meta: {
      source_rows: rows.length,
      returned_rows: rows.length,
      sampled: false,
      files: 1,
      available_start: null,
      available_end: null,
    },
  };
}

const defaults = {
  chart,
  primary: response([{ technology: "wind", year: 2030, value: 2 }]),
  comparison: null,
  mappings: { wind: { label: "Wind power", color: "#123456" } },
  legendValues: ["wind", "solar", "storage"],
  hiddenLegendValues: new Set<string>(),
  difference: false,
};

describe("buildPlotTraces", () => {
  it("preserves mapped labels, colours, and comparison styling", () => {
    const { data, hasColumnTotals } = buildPlotTraces({
      ...defaults,
      comparison: response([{ technology: "wind", year: 2030, value: 5 }]),
    });
    expect(hasColumnTotals).toBe(true);
    expect(data).toHaveLength(3);
    expect(data[0]).toMatchObject({
      type: "bar",
      name: "Wind power",
      x: [2030],
      y: [2],
      marker: { color: "#123456" },
      line: { color: "#123456", width: 2, dash: "solid" },
      opacity: 0.94,
      visible: true,
    });
    expect(data[1]).toMatchObject({
      type: "bar",
      name: "Wind power",
      y: [5],
      marker: { color: "#123456" },
      line: { width: 1.5, dash: "dot" },
      opacity: 0.5,
    });
    expect(data[2]).toMatchObject({ name: "Column total", text: ["2"] });
  });

  it("excludes hidden and secondary series from totals while retaining their traces", () => {
    const { data, hasColumnTotals } = buildPlotTraces({
      ...defaults,
      chart: { ...chart, secondary_y_lab: ["storage"] },
      primary: response([
        { technology: "wind", year: 2030, value: 2 },
        { technology: "wind", year: 2030, value: 3 },
        { technology: "solar", year: 2030, value: 7 },
        { technology: "storage", year: 2030, value: 100 },
        { technology: "wind", year: 2040, value: -4 },
      ]),
      hiddenLegendValues: new Set(["solar"]),
    });
    expect(hasColumnTotals).toBe(true);
    expect(data[1]).toMatchObject({ visible: "legendonly", marker: { color: "#005ca9" } });
    expect(data[2]).toMatchObject({ yaxis: "y2" });
    expect(data[3]).toMatchObject({
      type: "scatter",
      mode: "text",
      x: [2030, 2040],
      y: [5, -4],
      text: ["5", "-4"],
      textposition: ["top center", "bottom center"],
      hoverinfo: "skip",
      showlegend: false,
    });
  });

  it("uses the net total as text and the positive stack as its position", () => {
    const { data } = buildPlotTraces({
      ...defaults,
      primary: response([
        { technology: "wind", year: 2030, value: 8 },
        { technology: "solar", year: 2030, value: -3 },
      ]),
    });
    expect(data[2]).toMatchObject({ y: [8], text: ["5"], textposition: ["top center"] });
  });

  it("calculates comparison minus primary including series missing from either scenario", () => {
    const { data, hasColumnTotals } = buildPlotTraces({
      ...defaults,
      comparison: response([{ technology: "solar", year: 2040, value: 5 }]),
      difference: true,
      hiddenLegendValues: new Set(["wind"]),
    });
    expect(hasColumnTotals).toBe(false);
    expect(data).toHaveLength(2);
    expect(data[0]).toMatchObject({ type: "bar", x: [2030], y: [-2], visible: "legendonly" });
    expect(data[1]).toMatchObject({ type: "bar", x: [2040], y: [5], visible: true });
  });

  it("stacks primary hourly areas and keeps comparison and difference lines unfilled", () => {
    const options = {
      ...defaults,
      chart: { ...chart, hourly: true, type: "area_share" } satisfies ChartDefinition,
      primary: response([{ technology: "wind", snapshot: "2030-01-01 00:00", value: 2 }]),
      comparison: response([{ technology: "wind", snapshot: "2030-01-01 00:00", value: 5 }]),
    };
    const areas = buildPlotTraces(options);
    expect(areas.hasColumnTotals).toBe(false);
    expect(areas.data[0]).toMatchObject({ type: "scatter", mode: "lines", fill: "tonexty", stackgroup: "one" });
    expect(areas.data[1]).toMatchObject({ type: "scatter", fill: "none", stackgroup: undefined });
    const difference = buildPlotTraces({ ...options, difference: true });
    expect(difference.data).toHaveLength(1);
    expect(difference.data[0]).toMatchObject({ type: "scatter", mode: "lines", x: ["2030-01-01 00:00"], y: [3] });
    expect(difference.data[0]).not.toHaveProperty("stackgroup");
    expect(difference.data[0]).not.toHaveProperty("fill");
    const lines = buildPlotTraces({ ...options, chart: { ...options.chart, type: "hourly_line" } });
    expect(lines.data[0]).toMatchObject({ type: "scatter", mode: "lines" });
    expect(lines.data[0]).not.toHaveProperty("fill");
  });

  it("omits total labels for grouped bars, empty data, and fully hidden data", () => {
    expect(buildPlotTraces({ ...defaults, chart: { ...chart, type: "grouped_bar" } }).hasColumnTotals).toBe(false);
    expect(buildPlotTraces({ ...defaults, primary: response([]) })).toEqual({ data: [], hasColumnTotals: false });
    expect(buildPlotTraces({ ...defaults, hiddenLegendValues: new Set(["wind"]) }).hasColumnTotals).toBe(false);
  });
});
