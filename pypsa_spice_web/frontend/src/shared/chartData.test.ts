import { describe, expect, it } from "vitest";
import { buildDifferenceRows, getSharedXAxisRange, getSharedYAxisRanges } from "./chartData";
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

function response(rows: ResultRow[], availability: { start?: string; end?: string } = {}): ChartResponse {
  return {
    rows,
    dimensions: {},
    meta: {
      source_rows: rows.length,
      returned_rows: rows.length,
      sampled: false,
      files: 1,
      available_start: availability.start || null,
      available_end: availability.end || null,
    },
  };
}

describe("buildDifferenceRows", () => {
  it("matches the chart aggregation and calculates comparison minus primary", () => {
    const primary = response([
      { country: "AA", technology: "wind", unit: "GW", year: 2030, value: 2 },
      { country: "BB", technology: "wind", unit: "GW", year: 2030, value: 3 },
      { country: "AA", technology: "solar", unit: "GW", year: 2030, value: 4 },
    ]);
    const comparison = response([
      { country: "AA", technology: "wind", unit: "GW", year: 2030, value: 8 },
      { country: "AA", technology: "gas", unit: "GW", year: 2030, value: 1 },
    ]);

    expect(buildDifferenceRows(primary, comparison, chart, "baseline", "policy")).toEqual([
      { scenario: "policy − baseline", technology: "wind", unit: "GW", year: 2030, value: 3 },
      { scenario: "policy − baseline", technology: "solar", unit: "GW", year: 2030, value: -4 },
      { scenario: "policy − baseline", technology: "gas", unit: "GW", year: 2030, value: 1 },
    ]);
  });

  it("uses snapshots for hourly difference rows", () => {
    const hourlyChart: ChartDefinition = {
      ...chart,
      id: "generation-hourly",
      table_name: "generation_hourly",
      type: "hourly_line",
      hourly: true,
    };
    const primary = response([{ snapshot: "2030-01-01 00:00", technology: "wind", value: 2 }]);
    const comparison = response([{ snapshot: "2030-01-01 00:00", technology: "wind", value: 5 }]);

    expect(buildDifferenceRows(primary, comparison, hourlyChart, "baseline", "policy")).toEqual([
      {
        scenario: "policy − baseline",
        technology: "wind",
        snapshot: "2030-01-01 00:00",
        value: 3,
      },
    ]);
  });

  it("omits series whose difference is zero everywhere", () => {
    const primary = response([
      { technology: "unchanged", year: 2030, value: 5 },
      { technology: "changing", year: 2030, value: 2 },
      { technology: "changing", year: 2040, value: 4 },
    ]);
    const comparison = response([
      { technology: "unchanged", year: 2030, value: 5 },
      { technology: "changing", year: 2030, value: 2 },
      { technology: "changing", year: 2040, value: 7 },
    ]);

    expect(buildDifferenceRows(primary, comparison, chart, "baseline", "policy")).toEqual([
      { scenario: "policy − baseline", technology: "changing", year: 2030, value: 0 },
      { scenario: "policy − baseline", technology: "changing", year: 2040, value: 3 },
    ]);
  });
});

describe("getSharedYAxisRanges", () => {
  it("uses exact source extents when hourly rows were sampled", () => {
    const hourlyChart: ChartDefinition = { ...chart, type: "hourly_line", hourly: true };
    const sampled = response([{ technology: "wind", snapshot: "2030-01-01 00:00", value: 4 }]);
    sampled.meta.sampled = true;
    sampled.meta.source_rows = 8760;
    sampled.meta.axis_extents = { primary: [-20, 100], secondary: null };

    expect(getSharedYAxisRanges([sampled], hourlyChart, new Set())).toEqual({
      primary: [-26, 106],
      secondary: undefined,
    });
  });

  it("uses the largest stacked extent across both scenarios", () => {
    const primary = response([
      { technology: "wind", year: 2030, value: 40 },
      { technology: "solar", year: 2030, value: 20 },
    ]);
    const comparison = response([
      { technology: "wind", year: 2030, value: 80 },
      { technology: "solar", year: 2030, value: 20 },
      { technology: "gas", year: 2040, value: -10 },
    ]);

    expect(getSharedYAxisRanges([primary, comparison], chart, new Set())).toEqual({
      primary: [-15.5, 105.5],
      secondary: undefined,
    });
  });

  it("excludes hidden series and calculates a shared secondary-axis range", () => {
    const dualAxisChart: ChartDefinition = {
      ...chart,
      type: "hourly_dual",
      hourly: true,
      secondary_y_lab: ["stateOfCharge"],
    };
    const primary = response([
      { technology: "power", snapshot: "2030-01-01 00:00", value: 4 },
      { technology: "hidden", snapshot: "2030-01-01 00:00", value: 400 },
      { technology: "stateOfCharge", snapshot: "2030-01-01 00:00", value: 40 },
    ]);
    const comparison = response([
      { technology: "power", snapshot: "2030-01-01 00:00", value: 10 },
      { technology: "stateOfCharge", snapshot: "2030-01-01 00:00", value: 80 },
    ]);

    expect(getSharedYAxisRanges([primary, comparison], dualAxisChart, new Set(["hidden"]))).toEqual({
      primary: [0, 10.5],
      secondary: [38, 82],
    });
  });
});

describe("getSharedXAxisRange", () => {
  const hourlyChart: ChartDefinition = { ...chart, type: "hourly_line", hourly: true };
  const primary = response([{ technology: "wind", snapshot: "2030-01-01 00:00:00", value: 2 }], {
    start: "2030-01-01 00:00:00",
    end: "2030-12-31 23:00:00",
  });
  const comparison = response([{ technology: "wind", snapshot: "2030-02-01 00:00:00", value: 5 }], {
    start: "2030-02-01 00:00:00",
    end: "2031-01-31 23:00:00",
  });

  it("uses the full hourly extent across both scenarios", () => {
    expect(getSharedXAxisRange([primary, comparison], hourlyChart)).toEqual([
      "2030-01-01 00:00:00",
      "2031-01-31 23:00:00",
    ]);
  });

  it("clamps the shared extent to the selected time range", () => {
    expect(getSharedXAxisRange([primary, comparison], hourlyChart, "2030-03-01T00:00", "2030-10-01T00:00")).toEqual([
      "2030-03-01T00:00",
      "2030-10-01T00:00",
    ]);
  });

  it("leaves yearly charts on their existing automatic x-axis", () => {
    expect(getSharedXAxisRange([primary, comparison], chart)).toBeUndefined();
  });
});
