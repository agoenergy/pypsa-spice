import { describe, expect, it } from "vitest";
import {
  DASHBOARD_FORMAT,
  LocalDashboardStore,
  dashboardFilename,
  importedDashboardCopy,
  parseDashboardDefinition,
  parseDashboardExport,
  serializeDashboard,
} from "./utility";
import type { DashboardChartConfig, DashboardDefinition } from "./types";

const chart: DashboardChartConfig = {
  sectionId: "power",
  chartId: "p1",
  sector: "p-i-t",
  scenarios: ["baseline", "policy"],
  mode: "difference",
  country: "ALL",
  year: "2030",
};

const dashboard: DashboardDefinition = {
  schemaVersion: 2,
  id: "dashboard-1",
  title: "Energy transition",
  description: "Selected results",
  dataset: "dataset",
  project: "project",
  rows: [{
    id: "heading-1",
    type: "heading",
    title: "Power system",
  }, {
    id: "row-1",
    type: "chart",
    chart,
  }],
  createdAt: "2026-07-25T10:00:00.000Z",
  updatedAt: "2026-07-25T10:00:00.000Z",
};

function schema1Dashboard() {
  return {
    ...dashboard,
    schemaVersion: 1,
    rows: [
      dashboard.rows[0],
      { id: "old-chart-row", type: "charts", chartCount: 2 },
    ],
    items: [
      { id: "left-chart", rowId: "old-chart-row", column: "left", ...chart },
      { id: "right-chart", rowId: "old-chart-row", column: "right", ...chart, chartId: "p2" },
    ],
  };
}

class MemoryStorage implements Storage {
  private values = new Map<string, string>();
  get length() { return this.values.size; }
  clear() { this.values.clear(); }
  getItem(key: string) { return this.values.get(key) ?? null; }
  key(index: number) { return [...this.values.keys()][index] ?? null; }
  removeItem(key: string) { this.values.delete(key); }
  setItem(key: string, value: string) { this.values.set(key, value); }
}

describe("dashboard configuration", () => {
  it("round-trips the canonical schema-2 export format", () => {
    const text = serializeDashboard(dashboard);
    expect(JSON.parse(text)).toMatchObject({ format: DASHBOARD_FORMAT, schemaVersion: 2 });
    expect(parseDashboardExport(text)).toEqual(dashboard);
  });

  it("rejects unrelated JSON without changing the message into a syntax error", () => {
    expect(() => parseDashboardExport('{"title":"not a dashboard"}')).toThrow(
      "This is not a PyPSA-SPICE dashboard configuration.",
    );
  });

  it("creates new dashboard and row identifiers when importing", () => {
    const copy = importedDashboardCopy(dashboard);
    expect(copy.id).not.toBe(dashboard.id);
    expect(copy.rows[0].id).not.toBe(dashboard.rows[0].id);
    expect(copy.rows[1].id).not.toBe(dashboard.rows[1].id);
    expect(copy.title).toBe(dashboard.title);
  });

  it("migrates a schema-1 two-chart row into full-width chart rows", () => {
    const parsed = parseDashboardDefinition(schema1Dashboard());

    expect(parsed.schemaVersion).toBe(2);
    expect(parsed.rows.map((row) => row.type)).toEqual(["heading", "chart", "chart"]);
    expect(parsed.rows.filter((row) => row.type === "chart").map((row) => row.chart.chartId)).toEqual(["p1", "p2"]);
  });

  it("migrates the pre-row heading-section format without losing charts", () => {
    const legacy: Record<string, unknown> = schema1Dashboard();
    const legacyItems = legacy.items as Record<string, unknown>[];
    delete legacy.rows;
    legacy.headingSections = [{
      id: "old-heading",
      title: "Old heading",
    }];
    legacy.items = [
      { ...legacyItems[0], rowId: undefined, headingId: "old-heading" },
    ];

    const parsed = parseDashboardDefinition(legacy);
    expect(parsed.rows.map((row) => row.type)).toEqual(["heading", "chart"]);
    expect(parsed.rows[0]).toMatchObject({ title: "Old heading" });
  });

  it("limits imported chart selections to two scenarios", () => {
    const overLimit = JSON.parse(JSON.stringify(dashboard));
    overLimit.rows[1].chart.mode = "scenario";
    overLimit.rows[1].chart.scenarios = ["baseline", "policy", "high-demand"];

    const parsed = parseDashboardDefinition(overLimit);
    expect(parsed.rows[1]).toMatchObject({
      type: "chart",
      chart: { scenarios: ["baseline", "policy"] },
    });
  });

  it("uses a portable filename", () => {
    expect(dashboardFilename("Germany: Power & Heat")).toBe("germany-power-heat.json");
  });
});

describe("LocalDashboardStore", () => {
  it("reads schema-1 browser storage and writes subsequent saves as schema 2", async () => {
    const storage = new MemoryStorage();
    storage.setItem("pypsa-spice-dashboards-v1", JSON.stringify({
      schemaVersion: 1,
      dashboards: { [dashboard.id]: schema1Dashboard() },
    }));
    const store = new LocalDashboardStore(storage);

    const migrated = await store.get(dashboard.id);
    expect(migrated?.schemaVersion).toBe(2);
    expect(migrated?.rows.filter((row) => row.type === "chart")).toHaveLength(2);

    await store.save(migrated!);
    expect(JSON.parse(storage.getItem("pypsa-spice-dashboards-v2") || "{}")).toMatchObject({
      schemaVersion: 2,
    });
  });

  it("saves, lists, loads, and deletes dashboards", async () => {
    const storage = new MemoryStorage();
    const store = new LocalDashboardStore(storage);
    await store.save(dashboard);
    store.setLastOpenedId(dashboard.id);

    expect(await store.list()).toEqual([{
      id: dashboard.id,
      title: dashboard.title,
      chartCount: 1,
      updatedAt: dashboard.updatedAt,
    }]);
    expect(await store.get(dashboard.id)).toEqual(dashboard);

    await store.delete(dashboard.id);
    expect(await store.list()).toEqual([]);
    expect(store.getLastOpenedId()).toBeNull();
  });
});
