import { describe, expect, it } from "vitest";
import {
  limitHourlyRows,
  pivotYearlyRows,
  summarizeRows,
  tableToCsv,
} from "./DataDialog";
import type { ResultRow } from "./types";

const yearlyRows: ResultRow[] = [
  { scenario: "base", country: "DE", technology: "wind", year: 2030, value: 4, unit: "GW" },
  { scenario: "base", country: "DE", technology: "wind", year: 2030, value: 1, unit: "GW" },
  { scenario: "base", country: "FR", technology: "solar", year: 2030, value: 3, unit: "GW" },
  { scenario: "base", country: "DE", technology: "wind", year: 2040, value: 7, unit: "GW" },
];

describe("yearly result-table transformations", () => {
  it("pivots years, keeps scenario first, and combines duplicate values", () => {
    const table = pivotYearlyRows(yearlyRows);

    expect(table.pivoted).toBe(true);
    expect(table.columns[0]).toBe("scenario");
    expect(table.columns.slice(-2)).toEqual(["2030", "2040"]);
    expect(table.rows.find((row) => row.country === "DE")?.["2030"]).toBe(5);
  });

  it("falls back to a scenario-first flat table when rows have no year", () => {
    const table = pivotYearlyRows([{ country: "DE", value: 4, scenario: "base", bus: "DE1" }]);

    expect(table.pivoted).toBe(false);
    expect(table.columns[0]).toBe("scenario");
  });

  it("sums all filtered rows rather than a paginated subset", () => {
    const table = pivotYearlyRows(yearlyRows);
    const totals = summarizeRows(table.rows, ["2030", "2040"]);

    expect(totals).toEqual({ "2030": 8, "2040": 7 });
  });

  it("exports the requested columns in their displayed order with CSV escaping", () => {
    const csv = tableToCsv([{ country: "DE", technology: 'wind, "onshore"', value: 5 }], ["technology", "country"]);

    expect(csv).toBe('technology,country\r\n"wind, ""onshore""",DE');
  });
});

describe("hourly result-table regression", () => {
  it("keeps the existing 500-row display limit", () => {
    const rows = Array.from({ length: 620 }, (_, index) => ({ snapshot: String(index), value: index }));

    expect(limitHourlyRows(rows)).toHaveLength(500);
    expect(limitHourlyRows(rows)[499].value).toBe(499);
  });
});
