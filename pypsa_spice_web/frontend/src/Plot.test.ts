import { describe, expect, it } from "vitest";
import { buildDifferenceRows } from "./Plot";
import type { ChartDefinition, ChartResponse, ResultRow } from "./types";

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
